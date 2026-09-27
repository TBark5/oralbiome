"""Provenance checks: citations are byte-identical everywhere, licenses are present."""
import hashlib
import re
import subprocess

import pytest

from oralbiome import config

PAPER = ('Guo Z., Yu X., Liu Y., Hu Q., Zhang Z., Zhang C., Li J. (2026). "A comparative analysis '
         'of oral microbial communities in hypertensive patients with and without chronic '
         'periodontitis." BMC Oral Health 26(1):836. doi:10.1186/s12903-026-08144-6')
DATASET = ('Guo, Ziyin (2025). OTU Abundance Data: Oral Microbiome in Hypertensive Patients '
           'with/without Periodontitis. figshare. Dataset. '
           'https://doi.org/10.6084/m9.figshare.29897750.v1')
TEXT_SUFFIXES = {".md", ".py", ".txt", ".toml", ""}


def _tracked_text_files():
    out = subprocess.run(["git", "ls-files"], cwd=config.ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        pytest.skip("not a git checkout")
    files = []
    for name in out.stdout.splitlines():
        path = config.ROOT / name
        if name.startswith(("data/", "results/")) or path.suffix not in TEXT_SUFFIXES:
            continue
        if path.is_file() and path.name != "test_provenance.py":
            files.append(path)
    return files


def test_every_paper_citation_is_canonical():
    found = 0
    for path in _tracked_text_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        # any full citation starts with the first author list; it must be the canonical string
        for match in re.finditer(r"Guo Z\., Yu X\.", text):
            found += 1
            assert text[match.start():match.start() + len(PAPER)] == PAPER, path.name
        assert "Guo et al" not in text, f"partial citation in {path.name}"
        assert "Guo Z. et al" not in text, f"partial citation in {path.name}"
    assert found >= 4


def test_every_dataset_citation_is_canonical():
    found = 0
    for path in _tracked_text_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        for match in re.finditer(r"Guo, Ziyin \(", text):
            found += 1
            assert text[match.start():match.start() + len(DATASET)] == DATASET, path.name
    assert found >= 3


def test_licenses_present_and_distinct():
    code = (config.ROOT / "LICENSE").read_text(encoding="utf-8")
    data = (config.ROOT / "LICENSE-DATA").read_text(encoding="utf-8")
    assert code.startswith("MIT License")
    assert "CC BY 4.0" in data and PAPER in data and DATASET in data


@pytest.mark.skipif(not config.RAW_ZIP.exists(), reason="raw data not cached")
def test_cached_zip_matches_published_checksum():
    digest = hashlib.md5(config.RAW_ZIP.read_bytes()).hexdigest()
    assert digest == "6a49330fadc5e773059c52a84f6af67c"
