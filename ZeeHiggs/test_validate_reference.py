"""Tests for HXSWG-rescaled Zee h125 production and BR crosscheck."""
import math
import unittest
from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.validate_reference import expected_width_scaling, COUPLINGS


class ValidationTests(unittest.TestCase):
    def test_source_squared_couplings(self):
        for order in ("NO","IO"):
            source=audit(order)
            with self.subTest(order=order):
                h=source["neutral_couplings"]["25"]
                for ch, (typ,key) in COUPLINGS.items():
                    k2=expected_width_scaling(source,ch)
                    self.assertTrue(math.isfinite(k2))
                    self.assertGreaterEqual(k2,0)
                    if typ=="boson":
                        self.assertAlmostEqual(k2,h["boson"][key]**2,delta=1e-12)
                    else:
                        s,p=h["fermion_re_im"][str(key)]
                        self.assertAlmostEqual(k2,s*s+p*p,delta=1e-12)
                self.assertAlmostEqual(expected_width_scaling(source,"gg"),
                                       h["boson"]["gg"]**2,delta=1e-12)
                self.assertAlmostEqual(expected_width_scaling(source,"gamgam"),
                                       h["boson"]["gamgam"]**2,delta=1e-12)

    def test_known_input_sensitivity(self):
        expected=1.03142251**2
        obs=expected_width_scaling(audit("NO"),"gg")
        self.assertAlmostEqual(obs,expected,delta=5e-8)

    def test_no_95CL_from_reference_validation(self):
        # SM-like kappa rescaling is a prediction method, not an experimental
        # likelihood.  The numerical runner must keep its physical gate UNKNOWN.
        from inspect import getsource
        from ZeeHiggs.validate_reference import validate_reference
        src=getsource(validate_reference)
        self.assertIn('"physical_higgs":"UNKNOWN"',src)
        self.assertIn('"physical_HiggsSignals_95CL":None',src)


if __name__=="__main__":
    unittest.main()
