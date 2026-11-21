"""Regression tests of source-preserving global-vacuum witness and Yukawa audits."""
import math
import unittest

import numpy as np

from Scan.zee_chargebreaking_vacuum import (
    VDIM, normal_tree_mass_squared, source_inputs, potential,check_tadpole,
    ew_normal,search_witness,evaluate,
)
from Scan.zee_yukawa_perturbativity import (
    read_native,analyze_matrices,evaluate_trajectory,MAX_2BODY_YUKAWA
)
from ZeeHiggs.spheno_audit import audit,read_slha
from pathlib import Path


class VacuumTests(unittest.TestCase):
    def test_exact_neutral_tree_tadpole_and_potential(self):
        p=source_inputs("NO")
        t=check_tadpole(p)
        self.assertFalse(t["native_matches_tree_stationarity"])
        self.assertAlmostEqual(t["native_minus_tree_GeV2"][0],-38576.3624,delta=1)
        self.assertAlmostEqual(t["native_minus_tree_GeV2"][1],-47246.6747,delta=1)
        x,v=ew_normal(p)
        self.assertAlmostEqual(v,potential(x,p))
        self.assertFalse(evaluate("NO",nstarts=8)["original_native_tree_global_EW_vacuum_fully_certified"])
        self.assertEqual(evaluate("IO",nstarts=8)["original_native_tree_global_EW_vacuum"],
                         "UNKNOWN_NEEDS_MATCHED_TREE_TADPOLES")

    def test_cubic_phase_reduction_is_exact(self):
        p=source_inputs("NO")
        r1,r2,rc=(100.,170.,290.)
        theta=.6
        phi=.8
        r12=r1*r2*math.cos(theta)
        y1,y2=(0.+0j,r1+0j)
        # Construct a generic doublet with given |Phi2|, overlap, det.
        phi1=np.array([0.+0j,r1+0j])
        phi2=np.array([r2*math.sin(theta)+0j,r2*math.cos(theta)*np.exp(1j*phi)])
        det=phi1[0]*phi2[1]-phi1[1]*phi2[0]
        chi=rc*det/abs(det) if det!=0 else rc
        # For real positive mu, choose conjugate chi with opposite sign
        # so cubic term 2mu Re(det*conj(chi)) = -2 |mu det| rc.
        chi=-rc*det/abs(det)
        u=np.vdot(phi1,phi2)
        q=p["quartics"]
        x,y,z=np.vdot(phi1,phi1).real,np.vdot(phi2,phi2).real,abs(chi)**2
        reconstructed=(
            p["m1sq"]*x+p["m2sq"]*y-2*p["m12sq"]*u.real+p["mCsq"]*z
            +.5*q["Lambda1"]*x*x+.5*q["Lambda2"]*y*y
            +q["Lambda3"]*x*y+q["Lambda4"]*abs(u)**2
            +q["Lambda5"]*(u*u).real+q["LambdaC"]*z*z
            +(q["K1"]*x+q["K2"]*y-2*q["K12"]*u.real)*z
            +2*p["mu4"]*np.real(det*np.conj(chi))
        )/VDIM**4
        reduced=potential([r1/VDIM,r2/VDIM,rc/VDIM,theta,phi],p)
        self.assertAlmostEqual(reconstructed,reduced,delta=2e-12*max(1,abs(reduced)))

    def test_constructive_deeper_chargebreaking_witness(self):
        p=source_inputs("NO")
        # Make the pure singlet classically tachyonic on purpose.
        p["mCsq"]=-3.e5
        expected=normal_tree_mass_squared(p["quartics"],p["vevs"],p["m12sq"])
        p["m1sq"],p["m2sq"]=expected
        out=search_witness(p,nstarts=8,seed=21)
        self.assertEqual(out["status"],"DEEPER_WITNESS_FOUND")
        self.assertTrue(out["charge_breaking_witness"])
        self.assertLess(out["delta_potential_div_v4"],0)
        self.assertEqual(out["global_EW_vacuum"],"FAIL_DEEPER_TREE_WITNESS")

    def test_non_exclusion_from_finite_seed_search(self):
        p=source_inputs("NO")
        p["m1sq"],p["m2sq"]=normal_tree_mass_squared(
            p["quartics"],p["vevs"],p["m12sq"])
        out=search_witness(p,nstarts=8,seed=13)
        self.assertFalse(out["certified_exhaustive_global_search"])
        self.assertIn(out["global_EW_vacuum"],("FAIL_DEEPER_TREE_WITNESS","UNKNOWN"))
        with self.assertRaises(ValueError):
            search_witness(p,nstarts=1)


class YukawaTests(unittest.TestCase):
    def test_native_all_complex_lepton_and_quark_Yukawas(self):
        for ordering in ("NO","IO"):
            with self.subTest(ordering=ordering):
                r=read_native(ordering)
                self.assertEqual(set(r["yukawas"]["singular_values"]),
                                 {"Yu","Yd","Ye","YEE2","FZEE"})
                self.assertEqual(r["yukawas"]["local_Yukawa_perturbativity"],"PASS")
                self.assertLess(r["yukawas"]["max_Yukawa_singular_value"],MAX_2BODY_YUKAWA)
                self.assertLess(r["yukawas"]["singular_values"]["FZEE"]["antisymmetry_residual"],1e-8)
                self.assertEqual(r["full_UV_perturbativity_to_10TeV"],
                                 "UNKNOWN_REQUIRES_SOURCE_RGE_TRAJECTORY")
                self.assertIsNone(r["globally_allowed"])

    def test_fake_oversized_Yukawa_is_detected(self):
        r=audit("NO")
        blocks,_=read_slha(Path(r["spectrum_path"]).read_text())
        modified=dict(blocks)
        modified["YEE2"]=[(f.copy(),raw) for f,raw in blocks["YEE2"]]
        for f,raw in modified["YEE2"]:
            if f[:2]==["1","1"]:
                f[2]="5.0"
        result=analyze_matrices(modified,native_Q_GeV=160.,cutoff_Q_GeV=160.)
        self.assertEqual(result["local_Yukawa_perturbativity"],"FAIL")
        self.assertTrue(result["RGE_to_cutoff_proven"])

    def test_high_scale_cannot_be_inferred_from_low_scale(self):
        r=read_native("NO")
        scale0=r["yukawas"]["native_Q_GeV"]
        row={"native_Q_GeV":scale0,"local_Yukawa_perturbativity":"PASS"}
        with self.assertRaises(ValueError):
            evaluate_trajectory([row],1e4)
        with self.assertRaises(ValueError):
            evaluate_trajectory([row,dict(row,native_Q_GeV=1000.)],1e4)
        d=evaluate_trajectory([row,dict(row,native_Q_GeV=1e4)],1e4)
        self.assertEqual(d["RGE_final_status"],"UNKNOWN_CONTINUOUS_FLOW_NOT_CERTIFIED")


if __name__=="__main__":
    unittest.main()
