"""Step 1 regression: NO/IO targets and exact inverse at native running point."""
import unittest
import tempfile
import json
import io
from contextlib import redirect_stdout
from unittest.mock import patch
from pathlib import Path

import numpy as np

from Neutrino.fit import fit, inverse, zee_matrix, target, complex_matrix, NUFIT, inspect_native


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
                self.assertAlmostEqual(report["sum_mnu_eV"],sum(mass),places=13)
                self.assertAlmostEqual(report["m_bb_eV"],abs(mnu[0,0]),places=13)
                self.assertAlmostEqual(report["m_beta_eV"]**2,
                                       float(np.sum(abs(u[0])**2 * mass**2)),places=13)
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


    def test_mass_and_majorana_phase_freedom(self):
        # Oscillations cannot fix the lightest mass or Majorana phases.
        # No full scan; vary each once and require exact forward closure.
        for ordering in ("NO", "IO"):
            for light, alpha21, alpha31 in ((.001,0.,0.),
                                              (.0093,np.pi/2,0.),
                                              (.03,np.pi/3,-np.pi/2)):
                with self.subTest(ordering=ordering, light=light):
                    target_m, mass, u = target(ordering,light,alpha21,alpha31)
                    np.testing.assert_allclose(u.T@target_m@u,
                                               np.diag(mass),atol=1e-14)
                    report = fit(REFERENCE,ordering,light,alpha21,alpha31)
                    self.assertLess(report["forward_residual"],1e-8)
                    self.assertTrue(report["native_revalidation_required"])

    def test_both_null_roots_and_four_complex_yukawa_directions(self):
        snap=json.loads(REFERENCE.read_text())["snapshot"]
        zp=np.asarray(snap["ZP"],float)
        beta=np.arctan2(snap["vu"],snap["vd"])
        sine=2*zp[1,2]*(-np.sin(beta)*zp[1,0]+np.cos(beta)*zp[1,1])
        leptons=np.asarray(snap["MFe"],float)
        charged=np.asarray(snap["MHm"],float)[1:]
        mnu,_,_=target("NO",.0093,.4,.9)
        solutions=[]
        for root in (0,1):
            for z in (np.zeros(4,complex),
                      np.array([1e-4+3e-5j,2e-4j,-3e-5,.0001])):
                F,D,res=inverse(mnu,leptons,charged,sine,1e-4,t=.25+.1j,
                                root=root,homogeneous=z,f_phase=.37)
                self.assertLess(res,1e-8)
                self.assertLess(np.linalg.norm(zee_matrix(F,D,leptons,charged,sine)-mnu)
                                /np.linalg.norm(mnu),1e-8)
                solutions.append((F,D))
        self.assertGreater(np.linalg.norm(solutions[0][1]-solutions[1][1]),1e-5)
        self.assertGreater(np.linalg.norm(solutions[0][0]-solutions[2][0]),1e-7)
        for chart in ((1,2,0),(2,0,1)):
            F,D,res=inverse(mnu,leptons,charged,sine,1e-4,t=.25+.1j,
                            root=0,chart=chart)
            self.assertLess(res,1e-8)
            self.assertLess(np.linalg.norm(zee_matrix(F,D,leptons,charged,sine)-mnu)
                            /np.linalg.norm(mnu),1e-8)
        with self.assertRaises(ValueError):
            inverse(mnu,leptons,charged,sine,1e-4,chart=(0,0,1))
        with self.assertRaises(ValueError):
            inverse(mnu,leptons,charged,sine,1e-4,root=2)
        with self.assertRaises(ValueError):
            inverse(mnu,leptons,charged,sine,1e-4,homogeneous=[1,2,3])

    def test_native_card_hash_is_required_for_bound_pass(self):
        from Neutrino.fit import main
        from Model.audit_neutrino_native import sha
        snapshot=json.loads(REFERENCE.read_text())["snapshot"]
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            card=root/"candidate.in"
            card.write_text("fixed candidate input\\n")
            native=root/"result.json"
            native.write_text(json.dumps({"input_sha256":sha(card),
                            "spectrum_snapshot_index":0,"snapshots":[snapshot]}))
            cmd=["python","--ordering","NO","--mlight-ev",".0093",
                 "--native-result",str(native),"--checked-card",str(card)]
            with patch("sys.argv",cmd),redirect_stdout(io.StringIO()) as stream:
                main()
            result=json.loads(stream.getvalue())
            self.assertEqual(result["status"],"NATIVE_3SIGMA_FAIL")
            self.assertTrue(result["native_revalidation_required"])
            card.write_text("tampered candidate input\\n")
            with patch("sys.argv",cmd),self.assertRaises(ValueError):
                main()

    def test_numerical_and_phase_input_rejection(self):
        for value in (-1e-3,float("nan"),float("inf"),.3):
            with self.assertRaises(ValueError):
                target("NO",value)
        with self.assertRaises(ValueError):
            target("NO",.001,float("nan"))


if __name__=="__main__":
    unittest.main()
