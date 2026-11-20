"""Regression contract for actual official HiggsTools v1.3 experimental output.

Reuses the saved output from the successful official workflow; this is NOT
a replacement for a fresh external-backend run, nor a 95% CL Zee verdict.
"""
import json
from pathlib import Path
import unittest

from ZeeHiggs.spheno_audit import audit

RESULT = Path(__file__).resolve().parents[1] / "Higgs/reference/s00017_higgstools_v1_3_executed.json"


class OfficialNumericalReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = json.loads(RESULT.read_text())

    def test_external_source_and_datasets_are_pinned(self):
        r=self.r
        self.assertEqual(r["status"],"HIGGSTOOLS_EXECUTED_BUT_PHYSICS_INCOMPLETE")
        self.assertEqual(r["_archival_provenance"]["original_workflow_artifact_id"],14627701307)
        self.assertEqual(len(r["points"]),2)
        self.assertEqual(r["backend"]["HBdataset"]["git_head"],
                         "3539f9adf3724f8752701208bf059d1658824cd2")
        self.assertEqual(r["backend"]["HSdataset"]["git_head"],
                         "3c6c2538d9ba0549ea0559cf9b38c18eb80c3af1")
        self.assertEqual(r["backend"]["HBdataset"]["n_json_files"],310)
        self.assertEqual(r["backend"]["HSdataset"]["n_json_files"],22)
        self.assertEqual(r["physical_higgs"],"UNKNOWN")
        self.assertIsNone(r["globally_allowed"])

    def test_source_bound_inputs_and_charged_coverage(self):
        for r in self.r["points"]:
            source=audit(r["point_id"].split("_")[-1])
            self.assertEqual(r["source"],source["source"])
            self.assertEqual(r["PDG_to_HiggsTools_ID"],{
                "25":"25","35":"35","36":"36","37":"37","9000037":"9000037"})
            self.assertEqual(r["physical_higgs"],"UNKNOWN")
            self.assertEqual(r["higgssignals_confidence"],"UNKNOWN")
            self.assertEqual(r["higgsbounds_full_model_constraint"],"UNKNOWN")
            self.assertIsNone(r["globally_allowed"])
            selected=r["branches"]["effective_coupling_h125_BR"]["selectedLimits"]
            self.assertEqual([r["particle"] for r in selected],["25","35","36"])
            self.assertTrue(r["branches"]["effective_coupling_h125_BR"][
                "HiggsBounds_allowed_with_incomplete_inputs"])
            self.assertTrue(all(0 <= z["obs_ratio"] < 1 for z in selected))
            self.assertTrue(any("Charged-Higgs production rates absent" in x
                                for x in r["missing"]))

    def test_observed_higgs_signals_chi2_reference(self):
        for p, reference, expected_eff, expected_exp in zip(
                self.r["points"],
                (151.55556002183712,151.55556065245997),
                (11.454614038732046,15.195553807522913),
                (26.577376039475894,30.004387013862072)):
            sm=p["higgsSignals_SM_chi2_same_h125_mass"]
            eff=p["branches"]["effective_coupling_h125_BR"]
            explicit=p["branches"]["all_SPheno_explicit_BR"]
            self.assertAlmostEqual(sm,reference,places=10)
            self.assertAlmostEqual(eff["chi2_HiggsSignals"]-sm,
                                   eff["delta_chi2_vs_same_mass_SM"],places=10)
            self.assertAlmostEqual(explicit["chi2_HiggsSignals"]-sm,
                                   explicit["delta_chi2_vs_same_mass_SM"],places=10)
            self.assertAlmostEqual(eff["delta_chi2_vs_same_mass_SM"],expected_eff,places=10)
            self.assertAlmostEqual(explicit["delta_chi2_vs_same_mass_SM"],expected_exp,places=10)
            self.assertGreater(explicit["chi2_HiggsSignals"],eff["chi2_HiggsSignals"])


if __name__=="__main__":
    unittest.main()
