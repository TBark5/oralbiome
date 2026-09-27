"""Tests for the data layer (parsing, taxonomy cleanup, flags) and the dashboard."""
import pytest

from oralbiome import config, data
from oralbiome.preprocessing import flag_atypical_samples

needs_raw = pytest.mark.skipif(not config.RAW_ZIP.exists(), reason="raw data not cached")
needs_results = pytest.mark.skipif(not (config.RESULTS / "summary.json").exists(),
                                   reason="run_all.py has not been run")


def test_clean_taxonomy_replaces_uninformative_names():
    raw = [["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__uncultured", "f__unclassified_x",
            "g__uncultured", "s__uncultured_bacterium"]]
    tax = data.clean_taxonomy(raw, ["ASV1"])
    assert tax.loc["ASV1", "class"] == "Bacilli"
    assert tax.loc["ASV1", "genus"] == "Unclassified Bacilli"
    assert tax.loc["ASV1", "species"] == "Unclassified Bacilli"


def test_group_table_parser():
    text = 'x "SamplesID\\ttreat1\\nH1\\tH\\nP2\\tP\\nTP3\\tTP" y'
    assert data._read_group_table(text) == {"H1": "H", "P2": "P", "TP3": "TP"}


@needs_raw
def test_real_dataset_structure():
    ds = data.parse_raw_zip()
    assert ds.counts.shape == (25540, 67)
    assert ds.metadata["group_code"].value_counts().to_dict() == {"P": 18, "T": 17, "H": 16, "TP": 16}
    assert (ds.counts.to_numpy() >= 0).all()
    for sp in config.RED_COMPLEX_SPECIES:
        assert (ds.taxonomy["species"] == sp).any()


@needs_raw
def test_flagged_samples_are_stable():
    flags = flag_atypical_samples(data.parse_raw_zip())
    assert sorted(flags[flags].index) == sorted(["H7", "P7", "P9", "P18", "T9", "T10", "TP7", "TP15"])


def test_synthetic_fallback_is_labelled(monkeypatch, tmp_path):
    monkeypatch.setattr(config, "DATA_PROCESSED", tmp_path)
    monkeypatch.setattr(config, "RESULTS", tmp_path)
    ds = data.load_dataset(force_synthetic=True)
    assert ds.is_synthetic
    assert (tmp_path / "DATA_IS_SYNTHETIC").exists()
    assert (tmp_path / f"{config.SYNTHETIC_TAG}_ground_truth.csv").exists()


@needs_results
def test_dashboard_runs_without_errors():
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(config.ROOT / "app.py"), default_timeout=90).run()
    assert not at.exception
    assert len(at.tabs) == 9
    at.toggle(key="flagged").set_value(False).run()
    at.selectbox(key="level").select("phylum").run()
    at.segmented_control(key="ordination").set_value("NMDS").run()
    assert not at.exception
