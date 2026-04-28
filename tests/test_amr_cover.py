"""Tests for repseq.amr_cover -- pure-logic functions only (no shelling out)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import click
import pandas as pd
import pytest

from repseq.amr_cover import (
    _find_plasmidfinder_db,
    greedy_set_cover,
    parse_kleborate,
    parse_plasmidfinder,
)
from repseq.plots import _col_drug_class

# ---------------------------------------------------------------------------
# parse_kleborate
# ---------------------------------------------------------------------------

class TestParseKleborate:
    """Tests for parse_kleborate with synthetic Kleborate TSVs."""

    def test_v2_single_sample(self, tmp_path: Path) -> None:
        tsv = tmp_path / "kleb.tsv"
        tsv.write_text(textwrap.dedent("""\
            strain\tAGly_genes\tBla_genes
            sample_A\taac(3)-IIa\tSHV-11;TEM-1
            sample_B\t-\tSHV-11
        """))
        matrix, features = parse_kleborate(str(tsv))
        assert "sample_A" in matrix.index
        assert "sample_B" in matrix.index
        assert len(features) > 0
        # sample_A should have AGly and both Bla genes
        agly_cols = [f for f in features if "AGly" in f]
        assert any(matrix.loc["sample_A", c] == 1 for c in agly_cols)
        # sample_B has no AGly
        assert all(matrix.loc["sample_B", c] == 0 for c in agly_cols)

    def test_empty_features(self, tmp_path: Path) -> None:
        tsv = tmp_path / "kleb.tsv"
        tsv.write_text(textwrap.dedent("""\
            strain\tsome_col
            sample_X\tfoo
        """))
        matrix, features = parse_kleborate(str(tsv))
        assert "sample_X" in matrix.index
        assert features == []

    def test_strips_path_from_strain(self, tmp_path: Path) -> None:
        tsv = tmp_path / "kleb.tsv"
        tsv.write_text(textwrap.dedent("""\
            strain\tAGly_genes
            /data/assemblies/sample_Z.fasta\taac(3)-IIa
        """))
        matrix, _ = parse_kleborate(str(tsv))
        assert "sample_Z" in matrix.index

    def test_v3_prefixed_columns(self, tmp_path: Path) -> None:
        tsv = tmp_path / "kleb.tsv"
        tsv.write_text(textwrap.dedent("""\
            strain\tklebsiella_pneumo_complex__amr__AGly_acquired
            sampleV3\taac(3)-IIa
        """))
        matrix, features = parse_kleborate(str(tsv))
        assert len(features) == 1
        assert matrix.loc["sampleV3", features[0]] == 1

    def test_missing_id_column_raises(self, tmp_path: Path) -> None:
        tsv = tmp_path / "kleb.tsv"
        tsv.write_text("col_a\tcol_b\n1\t2\n")
        with pytest.raises(click.ClickException):
            parse_kleborate(str(tsv))


# ---------------------------------------------------------------------------
# greedy_set_cover
# ---------------------------------------------------------------------------

class TestGreedySetCover:
    """Tests for the greedy set-cover algorithm."""

    def _make_matrix(self) -> pd.DataFrame:
        """3 samples, 4 features."""
        return pd.DataFrame(
            {
                "AMR:Bla:TEM-1": [1, 0, 0],
                "AMR:Bla:SHV-11": [1, 1, 0],
                "AMR:AGly:aac": [0, 1, 0],
                "REP:IncF": [0, 0, 1],
            },
            index=["s1", "s2", "s3"],
        )

    def test_selects_best_cover(self) -> None:
        m = self._make_matrix()
        selected = greedy_set_cover(m, exclude_samples=[], n_amr=2)
        # s1 covers 2 features, s2 covers 2 (1 new), s3 covers 1 (new)
        assert len(selected) == 2
        # First pick should be s1 or s2 (both cover 2); second should add new features
        assert set(selected).issubset({"s1", "s2", "s3"})

    def test_respects_exclude(self) -> None:
        m = self._make_matrix()
        selected = greedy_set_cover(m, exclude_samples=["s1", "s2"], n_amr=1)
        assert selected == ["s3"]

    def test_zero_budget(self) -> None:
        m = self._make_matrix()
        assert greedy_set_cover(m, [], 0) == []

    def test_empty_matrix(self) -> None:
        m = pd.DataFrame(index=["s1", "s2"])
        selected = greedy_set_cover(m, [], 1)
        assert len(selected) == 1

    def test_all_excluded(self) -> None:
        m = self._make_matrix()
        selected = greedy_set_cover(m, ["s1", "s2", "s3"], 2)
        assert selected == []


# ---------------------------------------------------------------------------
# _col_drug_class (from plots)
# ---------------------------------------------------------------------------

class TestColDrugClass:
    """Tests for the column-name to drug-class mapper."""

    def test_amr_column(self) -> None:
        assert _col_drug_class("AMR:Bla_acquired:TEM-1") == "Bla"

    def test_rep_column(self) -> None:
        assert _col_drug_class("REP:IncFII") == "REP"

    def test_unknown_column(self) -> None:
        assert _col_drug_class("something_else") == "Other"

    def test_amr_no_acquired_suffix(self) -> None:
        assert _col_drug_class("AMR:AGly:aac(3)") == "AGly"


# ---------------------------------------------------------------------------
# _find_plasmidfinder_db
# ---------------------------------------------------------------------------

class TestFindPlasminderDb:
    """Tests for database path resolution."""

    def test_returns_string(self) -> None:
        result = _find_plasmidfinder_db()
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# parse_plasmidfinder
# ---------------------------------------------------------------------------

class TestParsePlasmidFinder:
    """Tests for parse_plasmidfinder with synthetic data."""

    def test_adds_replicon_columns(self, tmp_path: Path) -> None:
        tsv = tmp_path / "pf.tsv"
        tsv.write_text(textwrap.dedent("""\
            sample_id\tPlasmid
            s1\tIncFII
            s1\tIncX4
            s2\tIncFII
        """))
        matrix = pd.DataFrame(
            {"AMR:Bla:TEM-1": [1, 0]},
            index=["s1", "s2"],
        )
        result = parse_plasmidfinder(str(tsv), matrix)
        assert "REP:IncFII" in result.columns
        assert "REP:IncX4" in result.columns
        assert result.loc["s1", "REP:IncFII"] == 1
        assert result.loc["s1", "REP:IncX4"] == 1
        assert result.loc["s2", "REP:IncFII"] == 1
        assert result.loc["s2", "REP:IncX4"] == 0

    def test_ignores_unknown_samples(self, tmp_path: Path) -> None:
        tsv = tmp_path / "pf.tsv"
        tsv.write_text(textwrap.dedent("""\
            sample_id\tPlasmid
            unknown_sample\tIncFII
        """))
        matrix = pd.DataFrame(
            {"AMR:Bla:TEM-1": [1]},
            index=["s1"],
        )
        result = parse_plasmidfinder(str(tsv), matrix)
        # Should not crash, and should not add unknown_sample
        assert "unknown_sample" not in result.index

    def test_empty_tsv(self, tmp_path: Path) -> None:
        tsv = tmp_path / "pf.tsv"
        tsv.write_text("sample_id\tPlasmid\n")
        matrix = pd.DataFrame(
            {"AMR:Bla:TEM-1": [1]},
            index=["s1"],
        )
        result = parse_plasmidfinder(str(tsv), matrix)
        assert result.shape == matrix.shape
