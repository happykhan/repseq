"""Tests for repseq.ghru -- GHRU sampling prioritisation method."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from repseq.ghru import (
    TIER_ORDER,
    assign_pp_codes,
    assign_rp_codes,
    assign_tier,
    compute_priority_rank,
    parse_plasmid_profile,
    parse_resist_pattern,
    run_ghru_from_csv,
    select_primary_batch,
    select_site_guarantees,
)


# ---------------------------------------------------------------------------
# assign_tier
# ---------------------------------------------------------------------------


class TestAssignTier:
    """Tests for resistance tier assignment."""

    def test_carbapenem_nonsus_is_tier1(self) -> None:
        assert assign_tier(["AMP", "ETP", "CIP"], carb_nonsus=True, colistin_r=False) == "Tier1_CR"

    def test_colistin_r_is_tier2(self) -> None:
        assert assign_tier(["AMP", "COL"], carb_nonsus=False, colistin_r=True) == "Tier2_MDR"

    def test_mdr_three_classes_is_tier2(self) -> None:
        # CIP=FQN, GEN=AMG, SXT=FOL => 3 classes (excluding PEN)
        assert assign_tier(["AMP", "CIP", "GEN", "SXT"], carb_nonsus=False, colistin_r=False) == "Tier2_MDR"

    def test_two_classes_is_tier3(self) -> None:
        # CIP=FQN, TCY=TET => 2 classes
        assert assign_tier(["AMP", "CIP", "TCY"], carb_nonsus=False, colistin_r=False) == "Tier3_nonMDR_resistant"

    def test_one_non_pen_class_is_tier3(self) -> None:
        # CIP=FQN => 1 class (not counting PEN)
        assert assign_tier(["AMP", "CIP"], carb_nonsus=False, colistin_r=False) == "Tier3_nonMDR_resistant"

    def test_amp_only_is_pansusceptible(self) -> None:
        # AMP is PEN class, excluded from MDR counting
        assert assign_tier(["AMP"], carb_nonsus=False, colistin_r=False) == "Tier4_pansusceptible"

    def test_no_drugs_is_pansusceptible(self) -> None:
        assert assign_tier([], carb_nonsus=False, colistin_r=False) == "Tier4_pansusceptible"

    def test_carbapenem_overrides_mdr(self) -> None:
        # Even with MDR drugs, carbapenem non-sus puts it in Tier 1
        assert assign_tier(["AMP", "CIP", "GEN", "SXT", "IPM"], carb_nonsus=True, colistin_r=False) == "Tier1_CR"


# ---------------------------------------------------------------------------
# parse_resist_pattern and parse_plasmid_profile
# ---------------------------------------------------------------------------


class TestParseResistPattern:
    """Tests for resist_pattern parsing."""

    def test_normal_drugs(self) -> None:
        assert parse_resist_pattern("AMP, CIP, GEN") == "AMP;CIP;GEN"

    def test_exclude_drugs(self) -> None:
        assert parse_resist_pattern("AMP, CTT, CIP", exclude_drugs={"CTT"}) == "AMP;CIP"

    def test_pansusceptible(self) -> None:
        assert parse_resist_pattern("pan-susceptible") == "pan-susceptible"
        assert parse_resist_pattern("") == "pan-susceptible"
        assert parse_resist_pattern(None) == "pan-susceptible"  # type: ignore[arg-type]

    def test_all_excluded(self) -> None:
        assert parse_resist_pattern("CTT, TGC", exclude_drugs={"CTT", "TGC"}) == "pan-susceptible"


class TestParsePlasmidProfile:
    """Tests for plasmid profile parsing."""

    def test_multiple_plasmids_sorted(self) -> None:
        assert parse_plasmid_profile("IncFII, IncFIB(K)") == "IncFIB(K), IncFII"

    def test_single_plasmid(self) -> None:
        assert parse_plasmid_profile("IncFII") == "IncFII"

    def test_empty_is_none(self) -> None:
        assert parse_plasmid_profile("") == "none"
        assert parse_plasmid_profile(None) == "none"  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# assign_rp_codes / assign_pp_codes
# ---------------------------------------------------------------------------


class TestAssignCodes:
    """Tests for RP and PP code assignment."""

    def test_rp_codes_by_frequency(self) -> None:
        df = pd.DataFrame({
            "resist_pattern": ["A;B", "A;B", "A;B", "C;D", "C;D", "E"],
        })
        codes = assign_rp_codes(df, tier_col=None)
        # A;B appears 3x -> RP1, C;D appears 2x -> RP2, E appears 1x -> RP3
        assert codes.iloc[0] == "RP1"
        assert codes.iloc[3] == "RP2"
        assert codes.iloc[5] == "RP3"

    def test_rp_codes_per_tier(self) -> None:
        df = pd.DataFrame({
            "resist_pattern": ["A;B", "A;B", "A;B", "A;B", "C;D"],
            "tier": ["Tier1_CR", "Tier1_CR", "Tier2_MDR", "Tier2_MDR", "Tier1_CR"],
        })
        codes = assign_rp_codes(df, tier_col="tier")
        # Tier1: A;B(2x) -> RP1, C;D(1x) -> RP2
        # Tier2: A;B(2x) -> RP3 (continues numbering)
        assert codes.iloc[0] == "RP1"  # Tier1 A;B
        assert codes.iloc[4] == "RP2"  # Tier1 C;D
        assert codes.iloc[2] == "RP3"  # Tier2 A;B

    def test_pp_codes_by_frequency(self) -> None:
        df = pd.DataFrame({
            "plasmid_profile": ["IncFII", "IncFII", "IncX4", "none"],
        })
        codes = assign_pp_codes(df)
        assert codes.iloc[0] == "PP1"  # IncFII most frequent
        assert codes.iloc[2] == "PP2"  # IncX4
        assert codes.iloc[3] == "PP3"  # none


# ---------------------------------------------------------------------------
# select_primary_batch
# ---------------------------------------------------------------------------


class TestSelectPrimaryBatch:
    """Tests for Site_Guarantee/Primary/Secondary batch assignment."""

    def test_one_per_site_per_combo_no_guarantee(self) -> None:
        """With guarantee_sites=False, behaves as original Primary/Secondary."""
        df = pd.DataFrame({
            "isolate_id": ["a", "b", "c", "d"],
            "laboratory": ["LAB1", "LAB1", "LAB2", "LAB1"],
            "spec_date": ["2024-01-01", "2024-06-01", "2024-03-01", "2024-09-01"],
            "tier": ["Tier1_CR", "Tier1_CR", "Tier1_CR", "Tier2_MDR"],
            "rp_pp_combo": ["RP1_PP1", "RP1_PP1", "RP1_PP1", "RP2_PP2"],
        })
        batch = select_primary_batch(df, guarantee_sites=False)
        # For RP1_PP1 at LAB1: b (2024-06-01) is most recent -> Primary
        # a is Secondary, c is Primary (only LAB2 entry for RP1_PP1)
        # d is Primary (only LAB1 entry for RP2_PP2)
        assert batch.loc[df[df["isolate_id"] == "b"].index[0]] == "Primary"
        assert batch.loc[df[df["isolate_id"] == "a"].index[0]] == "Secondary"
        assert batch.loc[df[df["isolate_id"] == "c"].index[0]] == "Primary"
        assert batch.loc[df[df["isolate_id"] == "d"].index[0]] == "Primary"

    def test_site_guarantee_one_per_site(self) -> None:
        """With guarantee_sites=True, every site gets at least one guarantee."""
        df = pd.DataFrame({
            "isolate_id": ["a", "b", "c", "d"],
            "laboratory": ["LAB1", "LAB1", "LAB2", "LAB3"],
            "spec_date": ["2024-01-01", "2024-06-01", "2024-03-01", "2024-09-01"],
            "tier": ["Tier1_CR", "Tier2_MDR", "Tier4_pansusceptible", "Tier4_pansusceptible"],
            "rp_pp_combo": ["RP1_PP1", "RP2_PP2", "RP3_PP3", "RP3_PP3"],
        })
        batch = select_primary_batch(df, guarantee_sites=True)
        # All 3 sites should have exactly one Site_Guarantee
        guarantee_sites_set = set()
        for idx, label in batch.items():
            if label == "Site_Guarantee":
                guarantee_sites_set.add(df.loc[idx, "laboratory"])
        assert guarantee_sites_set == {"LAB1", "LAB2", "LAB3"}

    def test_site_guarantee_best_tier_selected(self) -> None:
        """Site guarantee should pick the best-tier isolate from each site."""
        df = pd.DataFrame({
            "isolate_id": ["a", "b", "c"],
            "laboratory": ["LAB1", "LAB1", "LAB1"],
            "spec_date": ["2024-01-01", "2024-06-01", "2024-03-01"],
            "tier": ["Tier4_pansusceptible", "Tier1_CR", "Tier2_MDR"],
            "rp_pp_combo": ["RP3_PP3", "RP1_PP1", "RP2_PP2"],
        })
        batch = select_primary_batch(df, guarantee_sites=True)
        # b is Tier1_CR (best tier) -> Site_Guarantee
        assert batch.loc[df[df["isolate_id"] == "b"].index[0]] == "Site_Guarantee"

    def test_guarantee_not_duplicated_as_primary(self) -> None:
        """An isolate tagged Site_Guarantee should not also be tagged Primary."""
        df = pd.DataFrame({
            "isolate_id": ["a", "b"],
            "laboratory": ["LAB1", "LAB2"],
            "spec_date": ["2024-06-01", "2024-03-01"],
            "tier": ["Tier1_CR", "Tier4_pansusceptible"],
            "rp_pp_combo": ["RP1_PP1", "RP2_PP2"],
        })
        batch = select_primary_batch(df, guarantee_sites=True)
        # Both should be Site_Guarantee (one per site, and they're the only ones)
        assert (batch == "Site_Guarantee").sum() == 2
        # No Primary since all are guaranteed
        assert (batch == "Primary").sum() == 0


# ---------------------------------------------------------------------------
# compute_priority_rank (integration on synthetic data)
# ---------------------------------------------------------------------------


class TestComputePriorityRank:
    """Tests for the full hierarchical sort."""

    def _make_df(self, guarantee_sites: bool = False) -> pd.DataFrame:
        """Build a small synthetic dataset with known expected ranking."""
        rows = [
            # Tier1_CR, RP1_PP1 combo (total=3, sites=2)
            {"isolate_id": "T1A", "laboratory": "SiteA", "spec_date": "2024-06-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "T1B", "laboratory": "SiteB", "spec_date": "2024-05-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "T1C", "laboratory": "SiteA", "spec_date": "2024-03-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            # Tier1_CR, RP2_PP2 combo (total=1, sites=1)
            {"isolate_id": "T1D", "laboratory": "SiteC", "spec_date": "2024-07-01",
             "tier": "Tier1_CR", "rp_code": "RP2", "pp_code": "PP2"},
            # Tier2_MDR, RP3_PP3 combo (total=1, sites=1)
            {"isolate_id": "T2A", "laboratory": "SiteA", "spec_date": "2024-04-01",
             "tier": "Tier2_MDR", "rp_code": "RP3", "pp_code": "PP3"},
            # Tier4_pansusceptible, RP4_PP4 combo (total=2, sites=2)
            {"isolate_id": "T4A", "laboratory": "SiteA", "spec_date": "2024-02-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP4", "pp_code": "PP4"},
            {"isolate_id": "T4B", "laboratory": "SiteB", "spec_date": "2024-01-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP4", "pp_code": "PP4"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        df["selection_batch"] = select_primary_batch(df, guarantee_sites=guarantee_sites)
        return df

    def test_tier_order_no_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        # All Tier1 before Tier2, all Tier2 before Tier4
        tier1_ranks = result[result["tier"] == "Tier1_CR"]["priority_rank"]
        tier2_ranks = result[result["tier"] == "Tier2_MDR"]["priority_rank"]
        tier4_ranks = result[result["tier"] == "Tier4_pansusceptible"]["priority_rank"]
        assert tier1_ranks.max() < tier2_ranks.min()
        assert tier2_ranks.max() < tier4_ranks.min()

    def test_guarantee_before_primary_before_secondary(self) -> None:
        df = self._make_df(guarantee_sites=True)
        result = compute_priority_rank(df)
        guarantee_ranks = result[result["selection_batch"] == "Site_Guarantee"]["priority_rank"]
        primary_ranks = result[result["selection_batch"] == "Primary"]["priority_rank"]
        secondary_ranks = result[result["selection_batch"] == "Secondary"]["priority_rank"]
        if len(guarantee_ranks) > 0 and len(primary_ranks) > 0:
            assert guarantee_ranks.max() < primary_ranks.min()
        if len(primary_ranks) > 0 and len(secondary_ranks) > 0:
            assert primary_ranks.max() < secondary_ranks.min()

    def test_primary_before_secondary_no_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        # Within Tier1: all Primary before all Secondary
        t1 = result[result["tier"] == "Tier1_CR"]
        primary_ranks = t1[t1["selection_batch"] == "Primary"]["priority_rank"]
        secondary_ranks = t1[t1["selection_batch"] == "Secondary"]["priority_rank"]
        if len(secondary_ranks) > 0:
            assert primary_ranks.max() < secondary_ranks.min()

    def test_higher_total_first_no_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        # RP1_PP1 (total=3) primary entries should rank before RP2_PP2 (total=1)
        rp1_primary = result[
            (result["rp_pp_combo"] == "RP1_PP1") & (result["selection_batch"] == "Primary")
        ]
        rp2_primary = result[
            (result["rp_pp_combo"] == "RP2_PP2") & (result["selection_batch"] == "Primary")
        ]
        assert rp1_primary["priority_rank"].max() < rp2_primary["priority_rank"].min()

    def test_site_alphabetical_within_combo_no_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        # Within RP1_PP1 Primary: SiteA (T1A) before SiteB (T1B)
        rp1_primary = result[
            (result["rp_pp_combo"] == "RP1_PP1") & (result["selection_batch"] == "Primary")
        ].sort_values("priority_rank")
        labs = list(rp1_primary["laboratory"])
        assert labs == sorted(labs)

    def test_most_recent_date_selected_no_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        # For RP1_PP1 at SiteA: T1A (2024-06-01) should be Primary, T1C (2024-03-01) Secondary
        t1a = result[result["isolate_id"] == "T1A"].iloc[0]
        t1c = result[result["isolate_id"] == "T1C"].iloc[0]
        assert t1a["selection_batch"] == "Primary"
        assert t1c["selection_batch"] == "Secondary"

    def test_all_isolates_ranked(self) -> None:
        df = self._make_df(guarantee_sites=False)
        result = compute_priority_rank(df)
        assert list(result["priority_rank"]) == list(range(1, len(df) + 1))

    def test_all_isolates_ranked_with_guarantee(self) -> None:
        df = self._make_df(guarantee_sites=True)
        result = compute_priority_rank(df)
        assert list(result["priority_rank"]) == list(range(1, len(df) + 1))


# ---------------------------------------------------------------------------
# Site guarantee logic
# ---------------------------------------------------------------------------


class TestSiteGuarantee:
    """Tests for the site guarantee feature."""

    def test_every_site_represented(self) -> None:
        """Every unique site must have at least one Site_Guarantee isolate."""
        rows = [
            {"isolate_id": "a", "laboratory": "Manila", "spec_date": "2024-01-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "b", "laboratory": "Manila", "spec_date": "2024-06-01",
             "tier": "Tier2_MDR", "rp_code": "RP2", "pp_code": "PP2"},
            {"isolate_id": "c", "laboratory": "Cebu", "spec_date": "2024-03-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP3", "pp_code": "PP3"},
            {"isolate_id": "d", "laboratory": "Davao", "spec_date": "2024-09-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP3", "pp_code": "PP3"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        batch = select_primary_batch(df, guarantee_sites=True)
        guaranteed_sites = {df.loc[i, "laboratory"] for i, v in batch.items() if v == "Site_Guarantee"}
        assert guaranteed_sites == {"Manila", "Cebu", "Davao"}

    def test_tier4_only_site_still_guaranteed(self) -> None:
        """A site with only Tier4 isolates must still contribute one sequence."""
        rows = [
            {"isolate_id": "a", "laboratory": "SiteA", "spec_date": "2024-01-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "b", "laboratory": "SiteB", "spec_date": "2024-06-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP4", "pp_code": "PP4"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        batch = select_primary_batch(df, guarantee_sites=True)
        assert batch.loc[df[df["isolate_id"] == "b"].index[0]] == "Site_Guarantee"

    def test_guarantee_picks_best_tier_per_site(self) -> None:
        """Within a site, the guarantee should go to the best-tier isolate."""
        rows = [
            {"isolate_id": "a", "laboratory": "SiteA", "spec_date": "2024-01-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP4", "pp_code": "PP4"},
            {"isolate_id": "b", "laboratory": "SiteA", "spec_date": "2024-06-01",
             "tier": "Tier2_MDR", "rp_code": "RP2", "pp_code": "PP2"},
            {"isolate_id": "c", "laboratory": "SiteA", "spec_date": "2024-03-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        batch = select_primary_batch(df, guarantee_sites=True)
        # c is Tier1_CR -> should be the guarantee
        assert batch.loc[df[df["isolate_id"] == "c"].index[0]] == "Site_Guarantee"
        # a and b should not be guarantees
        assert batch.loc[df[df["isolate_id"] == "a"].index[0]] != "Site_Guarantee"
        assert batch.loc[df[df["isolate_id"] == "b"].index[0]] != "Site_Guarantee"

    def test_guarantee_ranks_before_primary(self) -> None:
        """Site_Guarantee isolates should rank before Primary in the output."""
        rows = [
            {"isolate_id": "a", "laboratory": "SiteA", "spec_date": "2024-01-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "b", "laboratory": "SiteB", "spec_date": "2024-06-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "c", "laboratory": "SiteA", "spec_date": "2024-03-01",
             "tier": "Tier1_CR", "rp_code": "RP2", "pp_code": "PP2"},
            {"isolate_id": "d", "laboratory": "SiteC", "spec_date": "2024-09-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP3", "pp_code": "PP3"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        df["selection_batch"] = select_primary_batch(df, guarantee_sites=True)
        result = compute_priority_rank(df)
        guarantee_max = result[result["selection_batch"] == "Site_Guarantee"]["priority_rank"].max()
        primary_rows = result[result["selection_batch"] == "Primary"]
        if len(primary_rows) > 0:
            primary_min = primary_rows["priority_rank"].min()
            assert guarantee_max < primary_min

    def test_no_guarantee_flag_disables(self) -> None:
        """With guarantee_sites=False, no Site_Guarantee labels appear."""
        rows = [
            {"isolate_id": "a", "laboratory": "SiteA", "spec_date": "2024-01-01",
             "tier": "Tier1_CR", "rp_code": "RP1", "pp_code": "PP1"},
            {"isolate_id": "b", "laboratory": "SiteB", "spec_date": "2024-06-01",
             "tier": "Tier4_pansusceptible", "rp_code": "RP4", "pp_code": "PP4"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        batch = select_primary_batch(df, guarantee_sites=False)
        assert "Site_Guarantee" not in batch.values

    def test_single_isolate_per_site_is_guarantee_and_not_secondary(self) -> None:
        """If a site has exactly one isolate, it should be Site_Guarantee, not Secondary."""
        rows = [
            {"isolate_id": "a", "laboratory": "SiteA", "spec_date": "2024-01-01",
             "tier": "Tier3_nonMDR_resistant", "rp_code": "RP1", "pp_code": "PP1"},
        ]
        df = pd.DataFrame(rows)
        df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]
        batch = select_primary_batch(df, guarantee_sites=True)
        assert batch.iloc[0] == "Site_Guarantee"


# ---------------------------------------------------------------------------
# Flat CSV input
# ---------------------------------------------------------------------------


class TestRunGhruFromCsv:
    """Tests for run_ghru_from_csv -- flat CSV input format."""

    def _write_csv(self, tmp_path: Path, content: str) -> Path:
        csv_path = tmp_path / "input.csv"
        csv_path.write_text(content)
        return csv_path

    def test_basic_csv_with_tier(self, tmp_path: Path) -> None:
        """Minimal CSV with tier column works."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons\n"
            "ISO001,Manila,Tier1_CR,\"IncFIB(K), IncFII(K)\"\n"
            "ISO002,Cebu,Tier2_MDR,IncX4\n"
            "ISO003,Davao,Tier4_pansusceptible,\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        assert len(result) == 3
        assert (out / "priority_full.tsv").exists()
        assert (out / "primary_batch.tsv").exists()
        assert (out / "selected.txt").exists()

    def test_csv_all_sites_guaranteed(self, tmp_path: Path) -> None:
        """Every site in the CSV gets a Site_Guarantee."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons\n"
            "ISO001,Manila,Tier1_CR,IncFIB(K)\n"
            "ISO002,Cebu,Tier4_pansusceptible,\n"
            "ISO003,Davao,Tier4_pansusceptible,\n"
            "ISO004,Manila,Tier2_MDR,IncX4\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        guarantee = result[result["selection_batch"] == "Site_Guarantee"]
        sites = set(guarantee["laboratory"])
        assert sites == {"Manila", "Cebu", "Davao"}

    def test_csv_with_precomputed_rp_pp(self, tmp_path: Path) -> None:
        """Pre-computed rp_code and pp_code columns skip derivation."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons,rp_code,pp_code\n"
            "ISO001,Manila,Tier1_CR,IncFIB(K),RP1,PP1\n"
            "ISO002,Cebu,Tier2_MDR,IncX4,RP2,PP2\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        assert result[result["isolate_id"] == "ISO001"].iloc[0]["rp_code"] == "RP1"
        assert result[result["isolate_id"] == "ISO002"].iloc[0]["pp_code"] == "PP2"

    def test_csv_with_spec_date(self, tmp_path: Path) -> None:
        """spec_date column is used as tiebreaker."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons,spec_date\n"
            "ISO001,Manila,Tier1_CR,IncFIB(K),2024-06-01\n"
            "ISO002,Manila,Tier1_CR,IncFIB(K),2024-01-01\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        # Both have the same profile, same site, same tier.
        # ISO001 (more recent) should rank higher
        iso1_rank = result[result["isolate_id"] == "ISO001"].iloc[0]["priority_rank"]
        iso2_rank = result[result["isolate_id"] == "ISO002"].iloc[0]["priority_rank"]
        assert iso1_rank < iso2_rank

    def test_csv_without_spec_date(self, tmp_path: Path) -> None:
        """Missing spec_date defaults to 2000-01-01."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons\n"
            "ISO001,Manila,Tier1_CR,IncFIB(K)\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        assert len(result) == 1

    def test_csv_missing_required_column_raises(self, tmp_path: Path) -> None:
        """Missing required column raises ValueError."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,tier,replicons\n"
            "ISO001,Tier1_CR,IncFIB(K)\n"
        ))
        out = tmp_path / "out"
        with pytest.raises(ValueError, match="site"):
            run_ghru_from_csv(str(csv), str(out))

    def test_csv_no_guarantee(self, tmp_path: Path) -> None:
        """With guarantee_sites=False, no Site_Guarantee labels."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,tier,replicons\n"
            "ISO001,Manila,Tier1_CR,IncFIB(K)\n"
            "ISO002,Cebu,Tier4_pansusceptible,\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out), guarantee_sites=False)
        assert "Site_Guarantee" not in result["selection_batch"].values

    def test_csv_tier_derivation(self, tmp_path: Path) -> None:
        """When tier column is absent but derivation columns present, tier is derived."""
        csv = self._write_csv(tmp_path, (
            "isolate_id,site,resist_pattern,drug_classes,carb_nonsus,replicons\n"
            "ISO001,Manila,AMP;CIP;GEN;SXT,3,False,IncFIB(K)\n"
            "ISO002,Cebu,pan-susceptible,0,False,\n"
            "ISO003,Davao,AMP;ETP,1,True,IncX4\n"
        ))
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(csv), str(out))
        # ISO003 has carb_nonsus=True -> Tier1_CR
        assert result[result["isolate_id"] == "ISO003"].iloc[0]["tier"] == "Tier1_CR"
        # ISO001 has 3 non-PEN classes -> Tier2_MDR
        assert result[result["isolate_id"] == "ISO001"].iloc[0]["tier"] == "Tier2_MDR"
        # ISO002 is pan-susceptible -> Tier4
        assert result[result["isolate_id"] == "ISO002"].iloc[0]["tier"] == "Tier4_pansusceptible"

    def test_csv_tsv_format(self, tmp_path: Path) -> None:
        """TSV format (tab-separated) works when file extension is .tsv."""
        tsv_path = tmp_path / "input.tsv"
        tsv_path.write_text(
            "isolate_id\tsite\ttier\treplicons\n"
            "ISO001\tManila\tTier1_CR\tIncFIB(K)\n"
        )
        out = tmp_path / "out"
        result = run_ghru_from_csv(str(tsv_path), str(out))
        assert len(result) == 1
