"""Fail-closed two-parameter *reference* HiggsSignals screening tests.

Archived numerical inputs are real pinned HiggsTools runs, NOT mock chi2.
Synthetic points exercise the logic but imply no physics exclusions.
"""
import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

from ZeeHiggs.higgs_subset_screen import (
    REFERENCE, LOOP, HIGGS_BOUNDS, PRESCRIPTIONS, PRIMARY,
    TWO_PARAMETER_95, NOBS, screen_point, evaluate_archived,
    selected_higgsbounds_diagnostic, validate_receipts,
)


class HiggsReferenceScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=evaluate_archived()

    def test_published_two_parameter_statistical_constant(self):
        self.assertAlmostEqual(TWO_PARAMETER_95,5.991464547107982,places=13)
        self.assertEqual(NOBS,159)
        self.assertEqual(self.r["reference_data"]["HS_observable_count_from_original_executed_output"],159)
        self.assertFalse(self.r["criteria_have_true_95pct_coverage_for_current_Zee_scan"])
        self.assertFalse(self.r["scanned_parameter_dimension_known"])
        self.assertFalse(self.r["model_nuisances_profiled"])

    def test_actual_native_NO_and_IO_are_robust_reference_screen_fails(self):
        self.assertEqual(self.r["status"],
                         "HIGGS_SCOPED_REFERENCE_CHI2_SCREEN_VALIDATED_NO_FULL_95CL")
        self.assertEqual([x["point_id"] for x in self.r["points"]],
                         ["s00017_NO","s00017_IO"])
        for x, expected_primary, expected_least in zip(
                self.r["points"], (27.685262515290475,31.17256523233496),
                (15.894871451690278,19.663821287799607)):
            with self.subTest(point=x["point_id"]):
                self.assertEqual(x["status"],"HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL")
                self.assertAlmostEqual(x["primary_delta_chi2"],expected_primary,places=10)
                self.assertAlmostEqual(x["min_delta_chi2"],expected_least,places=10)
                self.assertTrue(x["all_prescriptions_above_threshold"])
                self.assertFalse(x["all_prescriptions_below_or_equal_threshold"])
                self.assertGreater(x["reference_chi2_gap_to_clear_screen"],0)
                self.assertEqual(len(x["chi2_by_prescription"]),4)
                self.assertAlmostEqual(x["chi2_SM_div_observable_count"],.953,delta=.01)
                self.assertEqual(x["physical_higgs"],"UNKNOWN")
                self.assertIsNone(x["higgsSignals_physical_95CL"])
                self.assertIsNone(x["globally_allowed"])
                self.assertFalse(x["is_formal_95CL_model_exclusion"])
                self.assertTrue(x["selection_cannot_be_used_as_global_Higgs_cut"])
                self.assertFalse(x["profile_likelihood_for_full_Zee_performed"])
                self.assertFalse(x["HiggsBounds_historic_diagnostic"][
                    "selected_limit_input_matches_preferred_HXSWG_ggF"])
                self.assertEqual(x["HiggsBounds_historic_diagnostic"]["selected_particle_ids"],
                                 ["25","35","36"])
                self.assertEqual(x["HiggsBounds_historic_diagnostic"]["physical_HiggsBounds_constraint"],
                                 "UNKNOWN")

    def test_robust_synthetic_reference_pass_fail_borderline(self):
        sm=150.0
        p={k:150.2 for k in PRESCRIPTIONS}
        r=screen_point("toy_PASS",sm,p)
        self.assertEqual(r["status"],"HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_PASS")
        self.assertEqual(r["physical_higgs"],"UNKNOWN")
        self.assertIsNone(r["higgsSignals_physical_95CL"])
        q={k:160.0 for k in PRESCRIPTIONS}
        r=screen_point("toy_FAIL",sm,q)
        self.assertEqual(r["status"],"HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL")
        q[PRIMARY]=152.0
        r=screen_point("toy_BORDERLINE",sm,q)
        self.assertEqual(r["status"],"HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_BORDERLINE")
        self.assertFalse(r["all_prescriptions_above_threshold"])

    def test_malformed_values_or_inconsistent_method_set_fails_closed(self):
        valid={k:160. for k in PRESCRIPTIONS}
        for bad in (float("nan"),float("inf"),-1.,"0",True):
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                screen_point("toy",150.,{**valid,PRIMARY:bad})
        with self.assertRaises(ValueError):
            screen_point("",150.,valid)
        with self.assertRaises(ValueError):
            screen_point("test",150.,{k:v for k,v in valid.items() if k!=PRIMARY})
        with self.assertRaises(ValueError):
            screen_point("test",150.,valid,threshold=3.84)
        with self.assertRaises(ValueError):
            screen_point("test",float("nan"),valid)

    def test_provenance_mutations_cannot_promote_physical_exclusion(self):
        r=json.loads(REFERENCE.read_text())
        loop=json.loads(LOOP.read_text())
        hb=json.loads(HIGGS_BOUNDS.read_text())
        changes=[
            lambda x: x.__setitem__("physical_higgs","PASS"),
            lambda x: x.__setitem__("theory_sigma",.1),
            lambda x: x["source"].__setitem__("higgsSignals_dataset","v999"),
            lambda x: x["source"].__setitem__("original_reference_JSON_sha256","missing"),
            lambda x: x.__setitem__("all_chi2_same_HS_dataset_and_spectrum",False),
            lambda x: x["results"]["NO"]["HiggsSignals_chi2"].__setitem__(
                "spheno_xs_effective_BRs",99.),
            lambda x: x["results"]["NO"]["HiggsSignals_chi2_delta_vs_SM"].__setitem__(
                "hxswg_xs_effective_BRs",0.),
        ]
        for mutate in changes:
            changed=copy.deepcopy(r)
            mutate(changed)
            with self.subTest(mutate=str(mutate)),self.assertRaises(ValueError):
                validate_receipts(changed,loop,hb)
        loop_bad=copy.deepcopy(loop)
        loop_bad["source"]["flags"]["calcHgamgam"]=True
        with self.assertRaisesRegex(ValueError,"native"):
            validate_receipts(r,loop_bad,hb)

    def test_historic_higgsbounds_neutral_only_not_full_Zee_pass(self):
        hb=json.loads(HIGGS_BOUNDS.read_text())
        for p in hb["points"]:
            d=selected_higgsbounds_diagnostic(p)
            self.assertEqual(d["status"],"HB_HISTORIC_NEUTRAL_SELECTED_NOT_EXCLUDED")
            self.assertEqual(d["selected_particle_ids"],["25","35","36"])
            self.assertFalse(d["selected_limit_input_matches_preferred_HXSWG_ggF"])
            self.assertFalse(d["is_full_Zee_HiggsBounds_pass"])
        bad=copy.deepcopy(hb["points"][0])
        bad["branches"]["effective_coupling_h125_BR"]["selectedLimits"][0]["obs_ratio"]=1.5
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            selected_higgsbounds_diagnostic(bad)

    def test_cli_creates_non_overwritten_source_bound_output(self):
        from ZeeHiggs.higgs_subset_screen import main
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"screen.json"
            with patch("sys.argv",["higgs_subset_screen","--output",str(p)]):
                main()
                d=json.loads(p.read_text())
                self.assertEqual(d["status"],
                         "HIGGS_SCOPED_REFERENCE_CHI2_SCREEN_VALIDATED_NO_FULL_95CL")
                self.assertEqual(d["points"][0]["status"],
                         "HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL")
                self.assertEqual(d["physical_higgs"],"UNKNOWN")
                with self.assertRaises(FileExistsError):
                    main()


if __name__=="__main__":
    unittest.main()
