"""Regression witnesses for archived *actually executed* official references.

The full numeric source artifacts can be reproduced independently from
the two successful original GitHub Actions runs. These tests never turn
a SM-only HDECAY control or HXSWG approximation into a model exclusion.
"""
import json
import math
from pathlib import Path
import unittest

BASE=Path(__file__).resolve().parents[1]/"Higgs/reference"


class ExecutedH125ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hx=json.loads((BASE/"s00017_HXSWG_ggF_width_validation_summary.json").read_text())
        cls.hd=json.loads((BASE/"s00017_independent_HDECAY_SM_control.json").read_text())

    def test_exact_external_provenance(self):
        h=self.hd
        self.assertEqual(h["source"]["HDECAY_git_commit"],"4884213286e1edbdf3b2bb9fe30095b29e336dff")
        self.assertEqual(h["source"]["artifact_id"],14703608040)
        self.assertEqual(h["source"]["original_JSON_sha256"],
                         "718567d56ac1d0ab171c34d41ca7c3cb8633d5cfba578d1494339db922da52e7")
        self.assertEqual(h["physical_higgs"],"UNKNOWN")
        self.assertIsNone(h["theory_sigma"])
        self.assertIsNone(h["full_model_allowed"])
        self.assertFalse(h["physical_accuracy"]["HDECAY_calculates_Zee"])
        self.assertIsNone(self.hx["theory_sigma"])
        self.assertEqual(self.hx["physical_higgs"],"UNKNOWN")
        self.assertIsNone(self.hx["higgsSignals_95CL_exclusion"])

    def test_official_hdecay_sm_widths(self):
        h=self.hd["results"]["NO"]
        self.assertAlmostEqual(h["SM_total_width_GeV"],.004101,places=10)
        self.assertAlmostEqual(h["SM_partial_width_GeV"]["ZZ"],.00011060397,places=10)
        self.assertAlmostEqual(h["SM_partial_width_GeV"]["bb"],.0023605356,places=10)
        self.assertAlmostEqual(h["SM_partial_width_GeV"]["gamgam"],.000010014642,places=10)
        self.assertAlmostEqual(self.hd["comparison"]["HDECAY_SM_div_HXSWG_SM_width"],
                               .004101/.00412288504122885,places=10)
        self.assertAlmostEqual(self.hd["comparison"][
            "native_Zee_SPheno_div_HDECAY_SM_kappa_proxy_NO"]["ZZ"],.8498484194913585,places=10)

    def test_hxswg_source_loop_coupling_and_four_chi2_corner(self):
        h=self.hx
        for order in ("NO","IO"):
            row=h["results"][order]
            p=row["HiggsSignals_chi2"]
            self.assertAlmostEqual(row["HXSWG_normalized_Zee_ggF_pb"],
                                   row["HXSWG_SM_ggF_pb"]*row["kappa_g2"],places=7)
            self.assertGreater(row["HXSWG_normalized_Zee_ggF_pb"],
                               row["native_spheno_Zee_ggF_pb"])
            for name,value in p.items():
                self.assertAlmostEqual(row["HiggsSignals_chi2_delta_vs_SM"][name],
                                       value-row["higgsSignals_SM_chi2"],places=9)
        self.assertAlmostEqual(h["results"]["NO"]["HiggsSignals_chi2"]["hxswg_xs_effective_BRs"],
                               179.2408225371276,places=8)
        self.assertAlmostEqual(h["results"]["IO"]["HiggsSignals_chi2"]["hxswg_xs_effective_BRs"],
                               182.72812588479493,places=8)


if __name__=="__main__":
    unittest.main()
