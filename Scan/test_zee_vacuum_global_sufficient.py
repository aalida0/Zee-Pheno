"""Mathematical and randomized witnesses for sufficient global Zee EW tree."""
import math
import unittest
import numpy as np

from Scan.zee_chargebreaking_vacuum import (
    source_inputs, normal_tree_mass_squared, potential, ew_normal, VDIM
)
from Scan.zee_vacuum_global_sufficient import auxiliary_tree_vacuum, from_native_point


class GlobalSufficientTests(unittest.TestCase):
    def test_native_tadpole_mismatch_prevents_false_vacuum_certification(self):
        for order in ("NO","IO"):
            p=source_inputs(order)
            a=from_native_point(p)
            self.assertEqual(a["matched_native_potential"]["status"],
                             "UNKNOWN_SOURCE_TADPOLE_SCHEME_MISMATCH")
            self.assertFalse(a["matched_native_potential"]["certified_global_tree_EW"])
            self.assertEqual(a["physical_global_vacuum"],"UNKNOWN")

    def test_tree_consistent_surrogate_proved_globally_minimal(self):
        for order in ("NO","IO"):
            p=source_inputs(order)
            r=from_native_point(p)["defined_tree_tadpole_potential"]
            self.assertEqual(r["status"],"CERTIFIED_SUFFICIENT_GLOBAL_TREE_EW")
            self.assertTrue(r["certified_global_tree_EW"])
            self.assertTrue(r["portal_K_PSD"])
            self.assertGreater(r["auxiliary_charged_mass_squared_GeV2"],0)
            self.assertGreater(r["auxiliary_CP_odd_mass_squared_GeV2"],0)
            self.assertGreater(min(r["auxiliary_CP_even_masses_squared_GeV2"]),0)
            self.assertGreater(r["2HDM_global_discriminant_GeV4"],0)
            self.assertFalse(r["native_potential_certified"])
            self.assertFalse(r["one_loop_or_thermal_global_stability_certified"])

    def test_exact_pointwise_lower_bound_over_random_complex_orbits(self):
        p=source_inputs("NO")
        t=normal_tree_mass_squared(p["quartics"],p["vevs"],p["m12sq"])
        p=dict(p,m1sq=t[0],m2sq=t[1])
        d=p["mu4"]**2/p["mCsq"]
        q=p["quartics"]
        rng=np.random.default_rng(2021)
        for _ in range(1800):
            radii=rng.uniform(0,2300,3)
            theta=float(rng.uniform(0,math.pi/2))
            phase=float(rng.uniform(0,math.pi))
            r1,r2,rc=radii
            z1,z2=r1*r1,r2*r2
            overlap=r1*r2*math.cos(theta)
            det2=z1*z2-overlap*overlap
            low=(
                p["m1sq"]*z1+p["m2sq"]*z2-2*p["m12sq"]*overlap*math.cos(phase)
                +.5*q["Lambda1"]*z1*z1+.5*q["Lambda2"]*z2*z2
                +q["Lambda3"]*z1*z2
                +(q["Lambda4"]+q["Lambda5"]*math.cos(2*phase))*overlap**2
                -d*det2
            )/VDIM**4
            actual=potential([r1/VDIM,r2/VDIM,rc/VDIM,theta,phase],p)
            self.assertGreaterEqual(actual-low,-2e-9*max(1,abs(actual),abs(low)))
        x,v_ew=ew_normal(p)
        r1,r2,rc=x[:3]*VDIM
        self.assertAlmostEqual(v_ew,potential(x,p),places=10)

    def test_counterexamples_do_not_false_pass(self):
        p=source_inputs("NO")
        args=(p["quartics"],*p["vevs"],p["m12sq"],p["mCsq"],p["mu4"])
        q=dict(p["quartics"],K12=3.)
        fail=auxiliary_tree_vacuum(q,*args[1:])
        self.assertEqual(fail["status"],"UNKNOWN_SUFFICIENT_CONDITIONS_NOT_MET")
        self.assertFalse(fail["certified_global_tree_EW"])
        q=dict(p["quartics"],Lambda3=-5.)
        fail=auxiliary_tree_vacuum(q,*args[1:])
        self.assertFalse(fail["certified_global_tree_EW"])
        self.assertEqual(
            auxiliary_tree_vacuum(p["quartics"],*p["vevs"],p["m12sq"],-1,p["mu4"])["status"],
            "UNKNOWN_SUFFICIENT_PRECONDITIONS_NOT_MET")
        fail=auxiliary_tree_vacuum(p["quartics"],*p["vevs"],p["m12sq"],p["mCsq"],p["mu4"],
                                   quad_masses=[p["m1sq"],p["m2sq"]])
        self.assertEqual(fail["status"],"UNKNOWN_SOURCE_TADPOLE_SCHEME_MISMATCH")


if __name__=="__main__":unittest.main()
