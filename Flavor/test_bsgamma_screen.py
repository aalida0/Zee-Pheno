"""Screening tests: archived real NLO evidence plus adversarial falsifiers."""
import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

from Flavor.bsgamma_screen import (
    S00017, REFERENCE, evaluate_archive, screen_prediction, validate_archived_data
)


class BsgammaScreenTests(unittest.TestCase):
    def test_2026_reference_contract(self):
        self.assertEqual(REFERENCE["sm_central"],3.54)
        self.assertEqual(REFERENCE["sm_sigma"],.14)
        self.assertEqual(REFERENCE["exp_central"],3.49)
        self.assertEqual(REFERENCE["exp_sigma"],.19)
        self.assertFalse(REFERENCE["published_2HDM_bound_reproduced_with_this_prescription"])

    def test_real_native_benchmarks_borderline(self):
        result = evaluate_archive()
        self.assertEqual(result["physical_bsgamma"],"UNKNOWN")
        self.assertIsNone(result["theory_sigma"])
        self.assertFalse(result["nominal_confidence_interval_is_validated_95CL"])
        self.assertEqual(len(result["points"]),2)
        for row in result["points"]:
            with self.subTest(point=row["point_id"]):
                self.assertEqual(row["screen_status"],"B_SCREEN_BORDERLINE")
                self.assertEqual(row["physical_constraint"],"B_PHYSICAL_UNKNOWN")
                self.assertIsNone(row["physically_allowed"])
                self.assertIsNone(row["experimental_exclusion_at_95CL"])
                self.assertFalse(row["method_spread_is_theory_uncertainty"])
                self.assertAlmostEqual(row["modern_SM_plus_HEPfit_NLO_delta_BR"],
                                       3.969996,delta=.000002)
                self.assertAlmostEqual(row["published_2026_linear_LO_total_in_1e-4"],
                                       4.046545,delta=.000002)
                self.assertAlmostEqual(row["method_spread_in_1e-4"],
                                       .076549,delta=.000002)
                self.assertAlmostEqual(row["fraction_of_NP_shift_needed_to_enter_band"],
                                       .04051,delta=.00001)
                interval=row["illustrative_gaussian_95_interval_in_1e-4"]
                self.assertAlmostEqual(interval[0],3.0274233901287273)
                self.assertAlmostEqual(interval[1],3.952576609871273)

    def test_null_np_returns_screen_pass_not_physical_pass(self):
        r = screen_prediction(0.,0.,0.,label="SM_like")
        self.assertEqual(r["screen_status"],"B_SCREEN_PASS")
        self.assertEqual(r["physical_constraint"],"B_PHYSICAL_UNKNOWN")
        self.assertIsNone(r["fraction_of_NP_shift_needed_to_enter_band"])

    def test_large_np_shift_screen_fail_not_physical_exclusion(self):
        # Identical LO and NLO shifts: no method-spread guard by construction.
        shift4=2.0
        c7=-shift4/8.51
        r=screen_prediction(shift4*1e-4,c7,0.,label="large_positive_shift")
        self.assertEqual(r["screen_status"],"B_SCREEN_FAIL")
        self.assertIsNone(r["experimental_exclusion_at_95CL"])

    def test_guard_stops_false_exclusion_from_small_method_difference(self):
        # Both methods outside naive band, but shift to cross it is smaller
        # than their disagreement; therefore screening is borderline.
        r=screen_prediction(4.299957011767789e-5,-.0464087353537987,
                            -.051669556843575014,label="NO")
        self.assertGreater(min(r["method_range_in_1e-4"]),
                           r["illustrative_gaussian_95_interval_in_1e-4"][1])
        self.assertLess(r["guarded_range_in_1e-4"][0],
                        r["illustrative_gaussian_95_interval_in_1e-4"][1])
        self.assertEqual(r["screen_status"],"B_SCREEN_BORDERLINE")

    def test_method_spread_is_not_gaussian_error(self):
        p=screen_prediction(2e-5,-.03,-.02,label="w")
        self.assertIsNone(p["NP_theory_sigma"])
        self.assertEqual(p["formal_statistical_interpretation"],"NONE_FOR_SCREEN")
        self.assertFalse(p["method_spread_is_theory_uncertainty"])

    def test_synthetic_fail_closed(self):
        for values in ((float("nan"),0,0),(0,float("inf"),0),("x",0,0)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                screen_prediction(*values,label="invalid")
        with self.assertRaises(ValueError):
            screen_prediction(0.,0.,0.,label=None)

    def test_archived_physical_gate_falsification(self):
        src=json.loads(S00017.read_text())
        mutations=(
            lambda o:o.__setitem__("physical_bsgamma","PASS"),
            lambda o:o.__setitem__("theory_sigma",.14e-4),
            lambda o:o.__setitem__("matching_scale_GeV",120.),
            lambda o:o["source"]["restored_backend"].__setitem__("validated_runs",1),
            lambda o:o["points"][0].__setitem__("BR_NLO_Egamma_gt_1p6GeV",.0005),
            lambda o:o["points"][0]["scalar"].__setitem__("card_sha256","a"*64),
            lambda o:o["points"][0]["matching"].__setitem__("C7_LO",float("nan")),
        )
        for mutate in mutations:
            with self.subTest(mutate=repr(mutate)):
                altered=copy.deepcopy(src)
                mutate(altered)
                with self.assertRaises(ValueError):
                    validate_archived_data(altered)

    def test_cli_output_complete_and_rejects_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"screen.json"
            from Flavor.bsgamma_screen import main
            import sys
            prior=sys.argv
            try:
                sys.argv=["bsgamma_screen","--output",str(out)]
                main()
                result=json.loads(out.read_text())
                self.assertEqual(result["status"],
                                 "TWO_METHOD_CONSERVATIVE_EXPLORATORY_SCREEN_ONLY")
                with self.assertRaises(FileExistsError):
                    main()
            finally:
                sys.argv=prior


if __name__=="__main__":
    unittest.main()
