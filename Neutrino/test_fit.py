"""Step 1 regression: NO/IO targets and exact inverse at native running point."""
import unittest
import tempfile
import json
from pathlib import Path

import numpy as np

from Neutrino.fit import fit, target, complex_matrix, NUFIT, inspect_native


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

    def test_known_pmns_is_recovered_from_native_takagi(self):
        # An independently constructed one-loop Takagi output must recover
        # the NO and IO NuFIT input angles and Dirac phase, not just masses.
        for ordering,mlight in (("NO",.0093),("IO",.001)):
            with self.subTest(ordering=ordering):
                mnu,m,u = target(ordering,mlight)
                native_order = [0,1,2] if ordering=="NO" else [2,0,1]
                u_sorted=u[:,native_order]
                mock={"snapshots":[{
                    "mat1":{"real":(mnu.real/1e9).tolist(),
                            "imag":(mnu.imag/1e9).tolist()},
                    "ZEL":{"real":np.eye(3).tolist(),
                           "imag":np.zeros((3,3)).tolist()},
                    "UV_1L":{"real":u_sorted.conj().T.real.tolist(),
                             "imag":u_sorted.conj().T.imag.tolist()},
                    "MFv_1L":(m[native_order]/1e9).tolist(),
                }],"spectrum_snapshot_index":0}
                with tempfile.TemporaryDirectory() as td:
                    path=Path(td)/"native.json"
                    path.write_text(json.dumps(mock))
                    screened=inspect_native(path,ordering)
                self.assertEqual(screened["status"],"NATIVE_3SIGMA_SCREEN_PASS")
                for key in ("s12","s13","s23","dm21","dm3l","delta_deg"):
                    self.assertAlmostEqual(
                        screened["parameters"][key], NUFIT[ordering][key], places=10)

    def test_known_unmatched_native_snapshot_is_rejected(self):
        # A matching matrix calculation is not an oscillation acceptance.
        snapshot = json.loads(REFERENCE.read_text())["snapshot"]
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"native.json"
            p.write_text(json.dumps({"snapshots":[snapshot],"spectrum_snapshot_index":0}))
            result = inspect_native(p,"NO")
        self.assertEqual(result["status"],"NATIVE_3SIGMA_SCREEN_FAIL")
        self.assertFalse(result["checks"]["dm21"])


    def test_numerical_and_phase_input_rejection(self):
        for value in (-1e-3,float("nan"),float("inf"),.3):
            with self.assertRaises(ValueError):
                target("NO",value)
        with self.assertRaises(ValueError):
            target("NO",.001,float("nan"))


if __name__=="__main__":
    unittest.main()
