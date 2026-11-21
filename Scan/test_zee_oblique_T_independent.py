"""Independent scalar T with both charged eigenstates / all neutral states."""
import math
import unittest

import numpy as np

from Scan.zee_oblique_T_independent import (
    f_kernel,source_point,independent_t,evaluate
)


class TOneLoopTests(unittest.TestCase):
    def test_kernel_positive_symmetric_and_nondegenerate(self):
        self.assertEqual(f_kernel(17.,17.),0.)
        for x,y in ((1.,2.),(1e4,1e6),(1234.,1234.000002)):
            self.assertAlmostEqual(f_kernel(x,y),f_kernel(y,x),
                                   delta=1e-8*max(1,abs(f_kernel(x,y))))
            self.assertGreaterEqual(f_kernel(x,y),0)
        for bad in (-1.,0.,float("nan")):
            with self.assertRaises(ValueError):
                f_kernel(bad,2.)

    def test_real_source_bound_NO_IO_two_charged_T(self):
        for ordering in ("NO","IO"):
            with self.subTest(ordering=ordering):
                p=source_point(ordering)
                out=evaluate(ordering)
                for setting in ("running_gauge","on_shell"):
                    row=out["results"][setting]
                    self.assertEqual(row["charged_states_included"],[37,9000037])
                    self.assertEqual(row["neutral_states_included"],[25,35,36])
                    self.assertAlmostEqual(row["SPheno_T"],.0120201035,delta=2e-5)
                    self.assertLess(abs(row["difference"]),.002)
                    self.assertGreater(row["independent_scalar_T"],0)
                    self.assertFalse(row["independently_validated_STU_complete"])
                # This is a narrow numerical crosscheck, NOT a full
                # official self-energy counterterm/renorm scheme certification.
                self.assertEqual(out["physical_full_EWPO"],"UNKNOWN")
                self.assertFalse(out["full_STU_independently_checked"])

    def test_rephasing_neutral_and_charged_mass_eigenstates(self):
        p=source_point("NO")
        orig=independent_t(p)["independent_scalar_T"]
        for seed in (1,4,22):
            rng=np.random.default_rng(seed)
            alt=dict(p)
            alt["U"]=p["U"]@np.diag(np.exp(1j*rng.uniform(-math.pi,math.pi,3)))
            # Real eigenstate rephasings allowed for CP even; signs
            # sufficient to check invariance of the quadratic kernels.
            alt["V"]=p["V"]@np.diag(rng.choice([-1,1],4))
            self.assertAlmostEqual(
                independent_t(alt)["independent_scalar_T"],
                orig,places=10)

    def test_running_vs_onshell_scheme_diff_is_nonzero(self):
        p=source_point("NO")
        r=independent_t(p,"running_gauge")
        os=independent_t(p,"on_shell")
        self.assertNotAlmostEqual(r["weak_angle_sin2"],
                                  os["weak_angle_sin2"],places=4)
        self.assertGreater(abs(r["independent_scalar_T"]-
                               os["independent_scalar_T"]),1e-5)
        with self.assertRaisesRegex(ValueError,"convention"):
            independent_t(p,"hybrid")


if __name__=="__main__":
    unittest.main()
