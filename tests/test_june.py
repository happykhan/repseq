"""Tests for repseq.june -- June Gayeta's sampling prioritisation method."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from repseq.june import (
    TIER_ORDER,
    assign_pp_codes,
    assign_rp_codes,
    assign_tier,
    compute_priority_rank,
    parse_plasmid_profile,
    parse_resist_pattern,
    select_primary_batch,
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
    """Tests for Primary/Secondary batch assignment."""

    def test_one_per_site_per_combo(self) -> None:
        df = pd.DataFrame({
            "isolate_id": ["a", "b", "c", "d"],
            "laboratory": ["LAB1", "LAB1", "LAB2", "LAB1"],
            "spec_date": ["2024-01-01", "2024-06-01", "2024-03-01", "2024-09-01"],
            "rp_pp_combo": ["RP1_PP1", "RP1_PP1", "RP1_PP1", "RP2_PP2"],
        })
        batch = select_primary_batch(df)
        # For RP1_PP1 at LAB1: b (2024-06-01) is most recent -> Primary
        # a is Secondary, c is Primary (only LAB2 entry for RP1_PP1)
        # d is Primary (only LAB1 entry for RP2_PP2)
        assert batch.loc[df[df["isolate_id"] == "b"].index[0]] == "Primary"
        assert batch.loc[df[df["isolate_id"] == "a"].index[0]] == "Secondary"
        assert batch.loc[df[df["isolate_id"] == "c"].index[0]] == "Primary"
        assert batch.loc[df[df["isolate_id"] == "d"].index[0]] == "Primary"


# ---------------------------------------------------------------------------
# compute_priority_rank (integration on synthetic data)
# ---------------------------------------------------------------------------


class TestComputePriorityRank:
    """Tests for the full hierarchical sort."""

    def _make_df(self) -> pd.DataFrame:
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
        df["selection_batch"] = select_primary_batch(df)
        return df

    def test_tier_order(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        # All Tier1 before Tier2, all Tier2 before Tier4
        tier1_ranks = result[result["tier"] == "Tier1_CR"]["priority_rank"]
        tier2_ranks = result[result["tier"] == "Tier2_MDR"]["priority_rank"]
        tier4_ranks = result[result["tier"] == "Tier4_pansusceptible"]["priority_rank"]
        assert tier1_ranks.max() < tier2_ranks.min()
        assert tier2_ranks.max() < tier4_ranks.min()

    def test_primary_before_secondary(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        # Within Tier1: all Primary before all Secondary
        t1 = result[result["tier"] == "Tier1_CR"]
        primary_ranks = t1[t1["selection_batch"] == "Primary"]["priority_rank"]
        secondary_ranks = t1[t1["selection_batch"] == "Secondary"]["priority_rank"]
        if len(secondary_ranks) > 0:
            assert primary_ranks.max() < secondary_ranks.min()

    def test_higher_total_first(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        # RP1_PP1 (total=3) primary entries should rank before RP2_PP2 (total=1)
        rp1_primary = result[
            (result["rp_pp_combo"] == "RP1_PP1") & (result["selection_batch"] == "Primary")
        ]
        rp2_primary = result[
            (result["rp_pp_combo"] == "RP2_PP2") & (result["selection_batch"] == "Primary")
        ]
        assert rp1_primary["priority_rank"].max() < rp2_primary["priority_rank"].min()

    def test_site_alphabetical_within_combo(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        # Within RP1_PP1 Primary: SiteA (T1A) before SiteB (T1B)
        rp1_primary = result[
            (result["rp_pp_combo"] == "RP1_PP1") & (result["selection_batch"] == "Primary")
        ].sort_values("priority_rank")
        labs = list(rp1_primary["laboratory"])
        assert labs == sorted(labs)

    def test_most_recent_date_selected_as_primary(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        # For RP1_PP1 at SiteA: T1A (2024-06-01) should be Primary, T1C (2024-03-01) Secondary
        t1a = result[result["isolate_id"] == "T1A"].iloc[0]
        t1c = result[result["isolate_id"] == "T1C"].iloc[0]
        assert t1a["selection_batch"] == "Primary"
        assert t1c["selection_batch"] == "Secondary"

    def test_all_isolates_ranked(self) -> None:
        df = self._make_df()
        result = compute_priority_rank(df)
        assert list(result["priority_rank"]) == list(range(1, len(df) + 1))
