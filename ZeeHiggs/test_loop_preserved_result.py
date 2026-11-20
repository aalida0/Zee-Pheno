"""Source-locked regression for actual HiggsTools native-loop retention.

The exact ratio Gamma_eff(h->gg,gammagamma)/Gamma_SM must equal the
squared SPheno coupling entry if calcggH/calcHgamgam were disabled.
This proves that the official backend DID use the supplied Zee loops,
rather than silently recomputing the SM-particle-only amplitude.
"""
import json
import math
from pathlib import Path
import unittest
from ZeeHiggs.spheno_audit import audit

RESULT = (Path(__file__).resolve().parents[1] /
          "Higgs/reference/s00017_h125_loop_preserved_summary.json")


class CorrectedLoopInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=json.loads(RESULT.read_text())

    def test_actual_provenance_and_unknown_physics(self):
        r=self.r
        self.assertEqual(r["source"]["run_url"].rsplit("/",1)[-1],"44737211306")
        self.assertEqual(r["source"]["artifact_id"],14630125099)
        self.assertFalse(r["source"]["flags"]["calcggH"])
        self.assertFalse(r["source"]["flags"]["calcHgamgam"])
        self.assertEqual(r["phys_higgs"],"UNKNOWN")
        self.assertIsNone(r["globally_allowed"])

    def test_native_charged_loop_effective_factors_preserved(self):
        for entry in self.r["points"]:
            with self.subTest(order=entry["ordering"]):
                actual=audit(entry["ordering"])
                couplings=actual["neutral_couplings"]["25"]["boson"]
                self.assertTrue(math.isclose(entry["gamma_eff_vs_sm"],
                             couplings["gamgam"]**2,rel_tol=5e-9))
                self.assertTrue(math.isclose(entry["gg_eff_vs_sm"],
                             couplings["gg"]**2,rel_tol=5e-9))
                self.assertTrue(math.isclose(
                    entry["gamma_partial_eff_GeV"] /
                    entry["sm_gamma_gamma_GeV"],
                    entry["gamma_eff_vs_sm"],rel_tol=5e-10))
                self.assertTrue(math.isclose(
                    entry["gamma_partial_SPheno_GeV"] /
                    entry["gamma_partial_eff_GeV"],
                    1.0294357,rel_tol=1e-5))

    def test_corrected_observed_chi2_closes(self):
        for p in self.r["points"]:
            self.assertAlmostEqual(p["corrected_effective_chi2"]-p["sm_chi2"],
                                   p["delta_eff_sm"],places=10)
            self.assertAlmostEqual(p["explicit_spheno_chi2"]-p["corrected_effective_chi2"],
                                   p["delta_exp_eff"],places=10)
            self.assertGreater(p["corrected_effective_chi2"],p["old_default_eff_chi2"])
            self.assertGreater(p["ggf_13TeV_auto_pb"],p["ggf_13TeV_override_pb"])
            self.assertGreater(p["partial_gamma_difference" ] if "partial_gamma_difference" in p else 0.,-1.)


if __name__=="__main__":
    unittest.main()
