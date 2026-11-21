"""Independent general scalar mixing S,T,U and Feynman integral checks."""
import math
import unittest

import numpy as np

from Scan.zee_oblique_T_independent import source_point, independent_t
from Scan.zee_oblique_STU_independent import (
    G,G_hat,G_tilde,independent_S_U,evaluate,MAX_DELTA_ABSOLUTE
)


class GeneralZeeObliqueSTUTests(unittest.TestCase):
    def test_G_is_symmetric_and_numerically_stable(self):
        for x,y,q in ((150**2,200**2,91**2),
                      (666.305**2,1184.55**2,91**2),
                      (125.25**2,125.25**2,91**2),
                      (160**2,180**2,80**2)):
            with self.subTest(x=x,y=y):
                self.assertAlmostEqual(G(x,y,q),G(y,x,q),delta=2e-12)
                self.assertAlmostEqual(G(x,y,q,64),G(x,y,q,128),delta=1e-10)
                self.assertTrue(math.isfinite(G(x,y,q)))
        with self.assertRaises(ValueError):
            G(0.,100.,90.)
        with self.assertRaises(ValueError):
            G(25.,36.,2e4)

    def test_Ghat_asymptotic_expansion_from_paper_C11(self):
        # For I/Q=1000, Eq C11 up to 1/I has O(1/I^2) residual.
        i,q=1000.,1.
        asym=math.log(i)-5/3-17/(2*i)
        self.assertAlmostEqual(G_hat(i,q),asym,delta=1.5e-4)
        self.assertLess(G(1000,1000,1),0)
        self.assertLess(G_tilde(1000,1,1),0)

    def test_NO_IO_full_two_charged_Grimus_STU(self):
        for order in ("NO","IO"):
            with self.subTest(order=order):
                res=evaluate(order)
                self.assertEqual(res["status"],"FULL_ONE_LOOP_SCALAR_STU_INDEPENDENT_CROSSCHECK")
                self.assertTrue(res["two_charged_states"])
                self.assertTrue(res["independent_one_loop_scalar_STU_pragmatic_consistency"])
                for scheme in ("running_gauge","on_shell"):
                    r=res["results"][scheme]
                    self.assertTrue(r["source_and_independent_agree_at_pragmatic_tolerance"])
                    for key in ("S","T","U"):
                        self.assertLess(abs(r["absolute_differences"][key]),MAX_DELTA_ABSOLUTE[key])
                    self.assertFalse(r["same_exact_renormalization_scheme_validated"])
                grimus=res["results"]["running_gauge"]["independent_Grimus_scalar_STU"]
                self.assertAlmostEqual(grimus["S"],.00303747136979,delta=1e-6)
                self.assertAlmostEqual(grimus["T"],.012179048717659,delta=1e-6)
                self.assertAlmostEqual(grimus["U"],.000027648552952,delta=2e-6)

    def test_rephasing_invariant_complete_STU(self):
        p=source_point("NO")
        q=independent_S_U(p)
        rng=np.random.default_rng(13)
        for _ in range(4):
            altered=dict(p)
            altered["U"]=p["U"]@np.diag(np.exp(1j*rng.uniform(-math.pi,math.pi,3)))
            altered["V"]=p["V"]@np.diag(rng.choice([-1,1],4))
            a=independent_S_U(altered)
            self.assertAlmostEqual(a["S"],q["S"],delta=1e-10)
            self.assertAlmostEqual(a["U"],q["U"],delta=1e-10)


if __name__=="__main__":unittest.main()
