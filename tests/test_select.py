"""Tests for repseq.select -- output-writing helpers with synthetic data."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from repseq.select import _write_coverage_summary, _write_report


class TestWriteReport:
    """Tests for _write_report TSV generation."""

    def _make_matrix(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "AMR:Bla:TEM-1": [1, 0, 1],
                "AMR:AGly:aac": [0, 1, 0],
                "REP:IncFII": [1, 0, 0],
            },
            index=["s1", "s2", "s3"],
        )

    def test_report_columns(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        _write_report(m, phylo_selected=["s1"], amr_selected=["s2"], output_dir=str(tmp_path))
        report = pd.read_csv(tmp_path / "report.tsv", sep="\t")
        expected_cols = {
            "sample_id", "selected", "slot_type",
            "amr_genes", "replicons", "n_amr_genes", "n_replicons",
        }
        assert set(report.columns) == expected_cols

    def test_report_slot_types(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        _write_report(m, phylo_selected=["s1"], amr_selected=["s2"], output_dir=str(tmp_path))
        report = pd.read_csv(tmp_path / "report.tsv", sep="\t")
        report = report.set_index("sample_id")
        assert report.loc["s1", "slot_type"] == "phylo"
        assert report.loc["s2", "slot_type"] == "amr"
        assert report.loc["s3", "slot_type"] == "not_selected"

    def test_report_selected_flag(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        _write_report(m, phylo_selected=["s1"], amr_selected=["s2"], output_dir=str(tmp_path))
        report = pd.read_csv(tmp_path / "report.tsv", sep="\t")
        report = report.set_index("sample_id")
        assert report.loc["s1", "selected"] == "yes"
        assert report.loc["s3", "selected"] == "no"

    def test_report_amr_gene_list(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        _write_report(m, phylo_selected=["s1"], amr_selected=[], output_dir=str(tmp_path))
        report = pd.read_csv(tmp_path / "report.tsv", sep="\t")
        report = report.set_index("sample_id")
        # s1 has Bla:TEM-1 -> amr_genes should contain "Bla:TEM-1"
        assert "Bla:TEM-1" in report.loc["s1", "amr_genes"]
        # s2 has no replicons
        assert report.loc["s2", "replicons"] == "-"


class TestWriteCoverageSummary:
    """Tests for _write_coverage_summary text generation."""

    def _make_matrix(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "AMR:Bla:TEM-1": [1, 0, 1],
                "AMR:AGly:aac": [0, 1, 0],
                "REP:IncFII": [1, 0, 0],
            },
            index=["s1", "s2", "s3"],
        )

    def test_summary_file_created(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        features = list(m.columns)
        _write_coverage_summary(m, ["s1", "s2"], features, str(tmp_path))
        summary = (tmp_path / "coverage_summary.txt").read_text()
        assert "repseq Coverage Summary" in summary

    def test_summary_counts(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        features = list(m.columns)
        _write_coverage_summary(m, ["s1"], features, str(tmp_path))
        summary = (tmp_path / "coverage_summary.txt").read_text()
        # s1 covers AMR:Bla:TEM-1 and REP:IncFII, but not AMR:AGly:aac
        assert "Total unique: 2" in summary  # 2 AMR features
        assert "Covered by selection: 1" in summary  # 1 AMR feature covered

    def test_full_coverage(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        features = list(m.columns)
        _write_coverage_summary(m, ["s1", "s2"], features, str(tmp_path))
        summary = (tmp_path / "coverage_summary.txt").read_text()
        assert "Coverage: 100.0%" in summary

    def test_selected_list_in_summary(self, tmp_path: Path) -> None:
        m = self._make_matrix()
        features = list(m.columns)
        _write_coverage_summary(m, ["s1", "s3"], features, str(tmp_path))
        summary = (tmp_path / "coverage_summary.txt").read_text()
        assert "s1" in summary
        assert "s3" in summary

    def test_empty_features(self, tmp_path: Path) -> None:
        m = pd.DataFrame(index=["s1", "s2"])
        _write_coverage_summary(m, ["s1"], [], str(tmp_path))
        summary = (tmp_path / "coverage_summary.txt").read_text()
        assert "Coverage: 100.0%" in summary
