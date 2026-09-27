"""Regenerate every result table and figure from the raw data.

Usage:
    python run_all.py              # real data (downloads once, then cached)
    python run_all.py --synthetic  # force the synthetic fallback dataset
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from oralbiome import composition, config, data, preprocessing, qc, style  # noqa: E402


def step(name: str):
    """Print a header and return a timer callback for one pipeline step."""
    print(f"\n=== {name} ===", flush=True)
    start = time.perf_counter()
    return lambda: print(f"    done in {time.perf_counter() - start:.1f} s", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic", action="store_true", help="use synthetic data")
    args = parser.parse_args()

    config.ensure_dirs()
    style.apply_style()

    done = step("Data layer")
    full = data.load_dataset(force_synthetic=args.synthetic)
    data.save_processed(full)
    report = qc.run(full)
    print(f"    {report['n_samples']} samples, {report['n_taxa']} taxa, "
          f"synthetic={report['is_synthetic']}")
    flags = preprocessing.flag_atypical_samples(full)
    full.metadata["flagged_atypical"] = flags.loc[full.metadata.index]
    full.metadata.to_csv(config.RESULTS / "sample_metadata_with_flags.csv")
    print(f"    flagged atypical samples: {list(flags[flags].index)}")
    done()

    summary: dict = {"is_synthetic": full.is_synthetic,
                     "flagged_samples": list(flags[flags].index)}
    primary = full.subset(config.PRIMARY_GROUPS)

    done = step("M1 preprocessing")
    pre = preprocessing.preprocess(primary, "primary")
    summary["m1"] = pre.summary
    done()

    done = step("M2 composition")
    summary["m2"] = composition.run(pre)
    done()

    (config.RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\nAll steps finished. Key numbers: results/summary.json")


if __name__ == "__main__":
    main()
