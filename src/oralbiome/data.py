"""Data layer: download-or-synthesize, caching, parsing and the quality report.

The real dataset is a BIOM 1.0 (JSON) ASV table from Guo et al. (2026),
deposited on figshare under CC BY 4.0. It is downloaded once and cached as
the original zip in ``data/raw``; every later run reads the cached file, so
the project runs offline. If the download fails and no cache exists, a
Dirichlet-multinomial synthetic table is generated instead and labelled
SYNTHETIC everywhere.
"""
from __future__ import annotations

import json
import re
import zipfile
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import config

RANKS: list[str] = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
_UNINFORMATIVE = re.compile(r"(uncultured|unclassified|unidentified|unassigned|metagenome|^$|^unknown)", re.I)


@dataclass
class Dataset:
    """Counts (taxa x samples), taxonomy (taxa x ranks) and sample metadata."""

    counts: pd.DataFrame
    taxonomy: pd.DataFrame
    metadata: pd.DataFrame
    is_synthetic: bool

    def subset(self, groups: tuple[str, ...]) -> "Dataset":
        """Return a copy restricted to the given group codes, dropping empty taxa."""
        meta = self.metadata[self.metadata["group_code"].isin(groups)].copy()
        counts = self.counts[meta.index]
        counts = counts.loc[counts.sum(axis=1) > 0]
        return Dataset(counts, self.taxonomy.loc[counts.index], meta, self.is_synthetic)


def download_raw(timeout: int = 60) -> bool:
    """Download the raw zip into ``data/raw`` unless it is already cached.

    Returns True if the cached file is available afterwards.
    """
    config.ensure_dirs()
    if config.RAW_ZIP.exists() and config.RAW_ZIP.stat().st_size > 0:
        return True
    try:
        import requests

        response = requests.get(config.DATA_URL, timeout=timeout)
        response.raise_for_status()
        config.RAW_ZIP.write_bytes(response.content)
        zipfile.ZipFile(config.RAW_ZIP).testzip()
        return True
    except Exception as exc:  # network errors, bad zip, HTTP errors
        print(f"[data] download failed ({exc!r}); falling back to synthetic data")
        if config.RAW_ZIP.exists():
            config.RAW_ZIP.unlink()
        return False


def _clean_name(value: str) -> str:
    """Strip the SILVA rank prefix (e.g. 'g__') from a taxonomy string."""
    return re.sub(r"^[a-z]__", "", value.strip())


def _informative(value: str) -> bool:
    return not _UNINFORMATIVE.search(value)


def clean_taxonomy(raw: list[list[str]], index: list[str]) -> pd.DataFrame:
    """Turn SILVA lineages into a tidy table with readable labels.

    Uninformative names ('uncultured', 'unclassified_X', ...) are replaced by
    'Unclassified <nearest informative higher rank>' so that unrelated
    unclassified ASVs are not merged into a single fake taxon.
    """
    rows = []
    for lineage in raw:
        names = [_clean_name(x) for x in lineage] + [""] * (len(RANKS) - len(lineage))
        cleaned, last = [], "Bacteria"
        for name in names[: len(RANKS)]:
            if _informative(name):
                cleaned.append(name)
                last = name
            else:
                cleaned.append(f"Unclassified {last}")
        rows.append(cleaned)
    return pd.DataFrame(rows, index=index, columns=RANKS)


def _read_group_table(text: str) -> dict[str, str]:
    """Extract the 'SamplesID<TAB>treat1' block from the vendor metadata file."""
    start = text.find("SamplesID\\ttreat1")
    if start < 0:
        return {}
    block = text[start : text.find('"', start)]
    pairs = [p.split("\\t") for p in block.split("\\n")[1:]]
    return {p[0].strip(): p[1].strip() for p in pairs if len(p) == 2 and p[0].strip()}


def parse_raw_zip() -> Dataset:
    """Parse the cached figshare zip into a :class:`Dataset`."""
    with zipfile.ZipFile(config.RAW_ZIP) as zf:
        names = {n.rsplit("/", 1)[-1]: n for n in zf.namelist()}
        biom = json.loads(zf.read(names["OTU_Feature_Table.tsv"]))
        meta_text = zf.read(names["Sample_Metadata.csv"]).decode("utf-8", "replace")

    n_rows, n_cols = biom["shape"]
    dense = np.zeros((n_rows, n_cols), dtype=np.int64)
    for r, c, v in biom["data"]:
        dense[int(r), int(c)] = int(round(v))
    taxa = [row["id"] for row in biom["rows"]]
    samples = [col["id"] for col in biom["columns"]]
    counts = pd.DataFrame(dense, index=taxa, columns=samples)
    taxonomy = clean_taxonomy([row["metadata"]["taxonomy"] for row in biom["rows"]], taxa)

    groups = _read_group_table(meta_text)
    if not groups:  # fall back to the sample-ID prefix, which encodes the group
        groups = {s: re.match(r"[A-Z]+", s).group(0) for s in samples}
    meta = pd.DataFrame({"sample_id": samples, "group_code": [groups[s] for s in samples]})
    meta["group"] = meta["group_code"].map(config.GROUP_LABELS)
    meta["read_depth"] = counts.sum(axis=0).to_numpy()
    meta = meta.set_index("sample_id")
    return Dataset(counts, taxonomy, meta, is_synthetic=False)


def load_dataset(force_synthetic: bool = False) -> Dataset:
    """Load the real dataset if possible, otherwise a labelled synthetic one."""
    config.ensure_dirs()
    marker = config.DATA_PROCESSED / "DATA_IS_SYNTHETIC"
    if not force_synthetic and download_raw():
        if marker.exists():
            marker.unlink()
        return parse_raw_zip()
    from .synthetic import generate_synthetic_dataset

    marker.write_text("The pipeline is running on SYNTHETIC data.\n")
    return generate_synthetic_dataset()


def save_processed(ds: Dataset) -> None:
    """Write the parsed tables to ``data/processed`` (gzipped CSV)."""
    tag = f"{config.SYNTHETIC_TAG}_" if ds.is_synthetic else ""
    ds.counts.to_csv(config.DATA_PROCESSED / f"{tag}asv_counts.csv.gz")
    ds.taxonomy.to_csv(config.DATA_PROCESSED / f"{tag}asv_taxonomy.csv.gz")
    ds.metadata.to_csv(config.DATA_PROCESSED / f"{tag}sample_metadata.csv")


def quality_report(ds: Dataset) -> dict:
    """Summarise read depth, sparsity and group sizes."""
    depth = ds.counts.sum(axis=0)
    per_group = {}
    for code, sub in ds.metadata.groupby("group_code"):
        d = depth[sub.index]
        per_group[code] = {
            "label": config.GROUP_LABELS.get(code, code),
            "n": int(len(sub)),
            "depth_min": int(d.min()),
            "depth_median": float(d.median()),
            "depth_max": int(d.max()),
            "observed_asvs_median": float((ds.counts[sub.index] > 0).sum(axis=0).median()),
        }
    return {
        "is_synthetic": ds.is_synthetic,
        "n_samples": int(ds.counts.shape[1]),
        "n_taxa": int(ds.counts.shape[0]),
        "total_reads": int(depth.sum()),
        "depth_min": int(depth.min()),
        "depth_median": float(depth.median()),
        "depth_max": int(depth.max()),
        "sparsity_fraction_zero": float((ds.counts.to_numpy() == 0).mean()),
        "singleton_taxa": int((ds.counts.sum(axis=1) == 1).sum()),
        "taxa_in_one_sample_only": int(((ds.counts > 0).sum(axis=1) == 1).sum()),
        "groups": per_group,
    }
