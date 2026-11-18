"""Step 1 regression: NO/IO targets and exact inverse at native running point."""
import unittest
from pathlib import Path

import numpy as np

from Neutrino.fit import fit, target, complex_matrix, NUFIT


REFERENCE = Path(__file__).resolve().parents[1]/"Model/reference/native_neutrino_snapshot.json"


class FitTests(unittest.TestCase):
    def test_nufit_targets_and_both_inverse_orders(self):
        for ordering, mlight in (("NO", .0093), ("IO", .001)):
            with self.subTest(ordering=ordering):
                mnu, mass, u = target(ordering, mlight)
                np.testing.assert_allclose(
                    u.T@mnu@u, np.diag(mass), rtol=1e-12, atol=1e-14)
                self.assertAlmostEqual(mass[1]**2-mass[0]**2, NUFIT[ordering]["dm21"], places=14)
                self.assertAlmostEqual(
                    mass[2]**2-mass[0 if ordering=="NO" else 1]**2,
                    NUFIT[ordering]["dm3l"], places=14)
                report = fit(REFERENCE,ordering,mlight)
                self.assertEqual(report["status"],"ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN")
                self.assertTrue(report["native_revalidation_required"])
                self.assertLess(report["forward_residual"],2e-9)
                self.assertLessEqual(report["max_abs_D"],1.)
                f = complex_matrix(report["F_gauge"])
                np.testing.assert_allclose(f+f.T,np.zeros((3,3)),atol=1e-14)

    def test_numerical_and_phase_input_rejection(self):
        for value in (-1e-3,float("nan"),float("inf"),.3):
            with self.assertRaises(ValueError):
                target("NO",value)
        with self.assertRaises(ValueError):
            target("NO",.001,float("nan"))


if __name__=="__main__":
    unittest.main()
