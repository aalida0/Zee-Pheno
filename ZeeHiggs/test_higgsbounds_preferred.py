"""Official HiggsBounds re-evaluation input-prescription mocks.

This test does not replace the pinned official experimental runtime;
it verifies that h125 alone switches 13 TeV ggF schemes and that
H/A and both charged mass eigenstates remain in HiggsTools inputs.
"""
import sys
import types
import unittest
from unittest.mock import patch

from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.higgsbounds_preferred import run_one


class HiggsBoundsPreferredInterfaceTest(unittest.TestCase):
    def test_same_spectrum_two_production_policies(self):
        data={
            "cxn_25_ggH_LHC13":48.3793449,
            "cxn_35_ggH_LHC13":.017,
            "cxn_36_ggH_LHC13":.015
        }
        ids={25:"25",35:"35",36:"36",37:"37",9000037:"9000037"}
        received=[]
        class FakeP:
            def __init__(self,inputs):
                self.old="cxn_25_ggH_LHC13" in inputs
                self.inputs=inputs
            def particle(self,identifier):
                if identifier!="25":
                    raise AssertionError(identifier)
                return self
            def cxn(self,collider,process):
                if collider!="LHC13" or process!="ggH":
                    raise ValueError()
                return 48.3793449 if self.old else 51.494807038
        def fake_make_inputs(report,hinput):
            return dict(data), dict(ids), [], {"cxn_25_ggH_LHC13":48.3793449}
        def fake_prediction(d,identifiers,hinput,*,explicit):
            self.assertFalse(explicit)
            self.assertEqual(identifiers[9000037],"9000037")
            self.assertEqual(d["cxn_35_ggH_LHC13"],.017)
            self.assertEqual(d["cxn_36_ggH_LHC13"],.015)
            received.append(dict(d))
            return FakeP(dict(d))
        class BoundOutput:
            allowed=True
        def bounds(pred):
            return BoundOutput()
        def fake_limits(x):
            return [
                {"particle":"25","obs_ratio":.73,"exp_ratio":1.0},
                {"particle":"35","obs_ratio":.19,"exp_ratio":.21},
                {"particle":"36","obs_ratio":.05,"exp_ratio":.05}
            ]
        m1=types.ModuleType("Higgs");m2=types.ModuleType("Higgs.tools")
        m3=types.ModuleType("Higgs.tools.Input")
        m2.Input=m3;m1.tools=m2
        with patch.dict(sys.modules,{"Higgs":m1,"Higgs.tools":m2,
                                     "Higgs.tools.Input":m3}),\
             patch("ZeeHiggs.higgsbounds_preferred.make_inputs",side_effect=fake_make_inputs),\
             patch("ZeeHiggs.higgsbounds_preferred.predictions_for",side_effect=fake_prediction),\
             patch("ZeeHiggs.higgsbounds_preferred.result_limits",side_effect=fake_limits):
            r=run_one("NO",bounds)
        self.assertEqual(len(received),2)
        self.assertIn("cxn_25_ggH_LHC13",received[0])
        self.assertNotIn("cxn_25_ggH_LHC13",received[1])
        self.assertEqual(r["old_spheno_ggF_13TeV_pb"],48.3793449)
        self.assertAlmostEqual(r["preferred_HXSWG_ggF_13TeV_pb"],51.494807038)
        self.assertEqual(r["preferred_HXSWG_ggF"]["selected_charged_particle_ids"],[])
        self.assertFalse(r["charged_production_complete"])
        self.assertEqual(r["neutral_selected_limits_physical_acceptance"],"UNKNOWN")
        self.assertEqual(r["physical_higgs"],"UNKNOWN")
        self.assertIsNone(r["globally_allowed"])


if __name__=="__main__":
    unittest.main()
