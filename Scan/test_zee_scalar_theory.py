"""Analytic BFB proof, full Zee high-energy unitarity, and PDG2025 STU tests."""
import copy
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from Scan.zee_scalar_theory import (NAMES,bfb,scattering_matrices,
    independently_known_eigenvalues,unitarity,evaluate as scalar_evaluate)
from Scan.zee_ewpo import (PDG2025,correlated_pdg2025,read_attested_theory,
    evaluate as native_evaluate)


def sample(**kw):
    d={"Lambda1":.6,"Lambda2":.4,"Lambda3":.25,"Lambda4":.1,
       "Lambda5":-.05,"LambdaC":.3,"K1":.25,"K2":.4,"K12":.11}
    d.update(kw)
    return d


class TheoryTests(unittest.TestCase):
    def test_analytic_portal_matrix_PSD_is_weaker_than_old_proxy(self):
        d=sample(K1=.04,K2=.64,K12=.12)
        # Old k1>=abs(k12), k2>=abs(k12) fails, but determinant
        # k1*k2-k12**2>0 proves non-negative for every complex Phi.
        self.assertLess(d["K1"],abs(d["K12"]))
        r=bfb(d)
        self.assertEqual(r["full_Zee_quartic_BFB"],"CERTIFIED_SUFFICIENT_PASS")
        self.assertGreater(r["portal_matrix_determinant"],0)
        self.assertFalse(r["BFB_is_global_EW_vacuum_proof"])

    def test_necessary_fail_vs_unproven_distinct(self):
        for d in [sample(Lambda1=-.4),sample(LambdaC=-.1),
                  sample(K1=-10)]:
            with self.subTest(p=d):
                self.assertEqual(bfb(d)["full_Zee_quartic_BFB"],"FAIL_NECESSARY")
        d=sample(K1=.01,K2=.01,K12=.3)
        r=bfb(d)
        self.assertEqual(r["full_Zee_quartic_BFB"],"UNKNOWN")
        self.assertEqual(r["portal_hermitian_PSD_sufficient"],"NOT_PROVEN")
        self.assertEqual(r["doublet_BFB_necessary"],"PASS")
        self.assertEqual(r["pure_singlet_BFB_necessary"],"PASS")

    def test_flat_quartic_boundary_is_not_falsely_excluded(self):
        # Eq.41 allows equality; BFB can be semipositive on flat rays.
        # Source V4 alone cannot rule out a full potential with LC=0.
        for p in (sample(LambdaC=0.), sample(Lambda1=0.)):
            b=bfb(p)
            self.assertEqual(b["full_Zee_quartic_BFB"],"UNKNOWN")
            self.assertNotEqual(b["doublet_BFB_necessary"],"FAIL")
            self.assertNotEqual(b["pure_singlet_BFB_necessary"],"FAIL")

    def test_random_scalar_fields_show_sufficient_quartic_nonnegative(self):
        # Test the exact Eq17 V4 on random complex doublet/singlet
        # directions to independently challenge the symbolic proof.
        rng=np.random.default_rng(3443)
        for p in (sample(),sample(K1=.04,K2=.64,K12=.12)):
            self.assertEqual(bfb(p)["full_Zee_quartic_BFB"],"CERTIFIED_SUFFICIENT_PASS")
            for _ in range(500):
                phi1=rng.normal(size=2)+1j*rng.normal(size=2)
                phi2=rng.normal(size=2)+1j*rng.normal(size=2)
                chi=rng.normal()+1j*rng.normal()
                x=float(np.vdot(phi1,phi1).real)
                y=float(np.vdot(phi2,phi2).real)
                u=np.vdot(phi1,phi2)
                z=abs(chi)**2
                v4=(.5*p["Lambda1"]*x*x+.5*p["Lambda2"]*y*y+
                    p["Lambda3"]*x*y+p["Lambda4"]*abs(u)**2+
                    p["Lambda5"]*(u*u).real+p["LambdaC"]*z*z+
                    z*(p["K1"]*x+p["K2"]*y-2*p["K12"]*u.real))
                self.assertGreaterEqual(v4,-1e-10)

    def test_published_E22_E29_all_five_scattering_sectors(self):
        p=sample()
        matrices=scattering_matrices(p)
        self.assertEqual({k:v.shape[0] for k,v in matrices.items()},
          {"Q2_Y1_E22":6,"Q1_Y1_E23":6,"Q1_Y0_E24":6,
           "Q0_Y1_E25":3,"Q0_Y0_E26_to_E29":13})
        raw=np.concatenate([np.linalg.eigvalsh(m) for m in matrices.values()])
        for lam in independently_known_eigenvalues(p):
            self.assertTrue(np.any(np.isclose(raw,lam,rtol=1e-10,atol=1e-8)),lam)
        u=unitarity(p)
        self.assertEqual(u["tree_high_energy_full_Zee_unitarity"],"PASS")
        self.assertEqual(u["analytic_15_eigenvalues_crosschecked"],"PASS")
        self.assertEqual(u["limit_abs_eigenvalue"],8*math.pi)
        self.assertFalse(u["unitarity_is_absolute_quartic_4pi_test"])

    def test_large_quartic_under_4pi_but_tree_unitarity_fails(self):
        # A conventional 4pi individual coupling bound cannot substitute
        # for the high-energy two-body scattering eigenvalue test.
        d=sample(Lambda1=10.)
        self.assertLess(max(abs(x) for x in d.values()),4*math.pi)
        self.assertEqual(unitarity(d)["tree_high_energy_full_Zee_unitarity"],"FAIL")
        self.assertEqual(scalar_evaluate(d)["max_abs_quartic_lt_4pi_conventional_prior"],"PASS")

    def test_doublet_relabeling_invariance_of_scattering(self):
        p=sample()
        q=dict(p, Lambda1=p["Lambda2"],Lambda2=p["Lambda1"],
               K1=p["K2"],K2=p["K1"])
        a=unitarity(p)
        b=unitarity(q)
        self.assertAlmostEqual(a["max_abs_scattering_eigenvalue"],
                               b["max_abs_scattering_eigenvalue"],places=10)

    def test_missing_nonfinite_or_wrong_scalar_input_rejected(self):
        for p in ({"Lambda1":.1},sample(K1=float("inf")),
                  sample(K2=float("nan")),sample(K12="0.1")):
            with self.assertRaises(ValueError):
                scalar_evaluate(p)


class STUTests(unittest.TestCase):
    def test_pdg2025_three_parameter_correlated_fit(self):
        mean={k:v for k,v in zip(("S","T","U"),PDG2025["mean_STU"])}
        fit=correlated_pdg2025(mean)
        self.assertAlmostEqual(fit["PDG2025_3D_correlated_chi2"],0,places=14)
        self.assertEqual(fit["status"],"PDG2025_STU_CORRELATED_3D_SCREEN_PASS")
        self.assertTrue(fit["PDG2025_U_floats"])
        self.assertEqual(fit["physical_EWPO"],"UNKNOWN")
        self.assertIsNone(fit["globally_allowed"])
        self.assertEqual(
            correlated_pdg2025({"S":.3,"T":-.5,"U":.8})["status"],
            "PDG2025_STU_CORRELATED_3D_SCREEN_FAIL")
        for bad in ({"S":0,"T":0},{"S":0,"T":0,"U":float("nan")}):
            with self.assertRaises(ValueError):
                correlated_pdg2025(bad)

    def test_real_source_bound_spheno_quartics_and_STU_NO_IO(self):
        for order in ("NO","IO"):
            with self.subTest(order=order):
                native=read_attested_theory(order)
                self.assertEqual(native["MINPAR_match_running_ZEEQUARTIC"],"PASS")
                self.assertEqual(native["MINPAR_match_ZEEPOT"],"PASS")
                self.assertEqual(set(native["quartics"]),set(NAMES))
                self.assertTrue(native["spheno_low_energy_enabled"])
                self.assertFalse(native["STU_independently_recomputed_with_two_mixed_charged_scalars"])
                self.assertAlmostEqual(native["quartics"]["Lambda1"],.991396137,places=7)
                self.assertAlmostEqual(native["quartics"]["LambdaC"],.180298685,places=7)
                self.assertAlmostEqual(native["STU"]["T"],.0120201035,delta=.001)
                self.assertAlmostEqual(native["STU"]["S"],.00296941567,delta=.001)
                row=native_evaluate(order)
                self.assertEqual(row["BFB_quartic_asymptotic"],"CERTIFIED_SUFFICIENT_PASS")
                self.assertEqual(row["tree_high_energy_Zee_unitarity"],"PASS")
                self.assertEqual(row["EWPO_correlated_3D_screen"],
                                 "PDG2025_STU_CORRELATED_3D_SCREEN_PASS")
                self.assertEqual(row["EWPO"]["physical_EWPO"],"UNKNOWN")
                self.assertEqual(row["global_vacuum"],"UNKNOWN")
                self.assertIsNone(row["globally_allowed"])

    def test_source_label_and_value_corruption_fail_closed(self):
        from Scan.zee_ewpo import block_1d
        with self.assertRaisesRegex(ValueError,"duplicate"):
            block_1d([(["1","0.5"],""),(["1","0.6"],"")],"STU")
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            block_1d([(["1","nan"],"")],"STU")


if __name__=="__main__":
    unittest.main()
