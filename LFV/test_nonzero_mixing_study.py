"""Nonzero Zee mixing and nonzero F: deterministic LFV truncation tests."""
import copy
import math
import unittest

import numpy as np

from LFV.nonzero_mixing_study import (
    DEFAULT_GRID, grid_points, evaluate_grid, _validate_scheme,
)
from LFV.test_constraints import sample


def inversion():
    return dict(soft_mass_sq_GeV2=640.**2,alignment_delta_rad=0.,
                lambda_c=.25,K1=.2,K2=.4,K12=.01)


def grid():
    return dict(tan_beta=[2.],gamma_rad=[.4],
                mH2_GeV=[900.],F_multiplier=[1.,2.])


class ZeeNonzeroMixingTests(unittest.TestCase):
    def test_fixed_grid_potential_reconstruction_without_tuning(self):
        p=sample()
        original=copy.deepcopy(p)
        points=list(grid_points(p,inversion(),grid()))
        self.assertEqual(len(points),2)
        self.assertEqual(p,original)
        self.assertEqual(points[0][0]["scalar"]["gamma"],.4)
        for point,flags in points:
            self.assertAlmostEqual(point["scalar"]["alpha"],point["scalar"]["beta"])
            self.assertAlmostEqual(flags["m12sq_GeV2"],
                                   640.**2*math.sin(point["scalar"]["beta"])*
                                   math.cos(point["scalar"]["beta"]),places=9)
            self.assertTrue(flags["simple_quartic_magnitude_below_4pi"])
            self.assertTrue(flags["doublet_only_BFB_necessary"])
        self.assertEqual(points[1][0]["scalar"]["gamma"],.4)
        self.assertAlmostEqual(
            points[0][0]["potential"]["Mu4"],
            points[1][0]["potential"]["Mu4"],
        )
        self.assertAlmostEqual(
            points[0][0]["scalar"]["mH2"],900.
        )
        for point,axis in points:
            self.assertGreater(point["scalar"]["mH2"]-point["scalar"]["mH1"],91.1876)
            self.assertNotEqual(point["potential"]["Mu4"],0.)
            self.assertTrue(axis["charged_cascade_open"])

    def test_nonzero_zee_and_F_squared_scaling(self):
        output=evaluate_grid(sample(),inversion(),grid(),.23129)
        self.assertEqual(output["status"],
                         "NONZERO_MIXING_ZEE_TRUNCATION_DIAGNOSTIC_ONLY")
        self.assertEqual(output["n_points"],2)
        self.assertEqual(output["n_simple_quartic_magnitude_below_4pi"],2)
        self.assertEqual(output["n_doublet_BFB_necessary_pass"],2)
        self.assertEqual(output["n_channel_rows"],6)
        for channel in ("mu_to_e_gamma","tau_to_e_gamma","tau_to_mu_gamma"):
            self.assertEqual(output["summaries"][channel]["n_points"],2)
            self.assertEqual(output["summaries"][channel]["full_Zee_LFV_gate"],"UNKNOWN")
        self.assertIsNone(output["physical_full_Zee_LFV_allowed_points"])
        for r in output["rows"]:
            self.assertGreater(r["Z12_over_e_abs"],0.)
            self.assertEqual(r["full_Zee_LFV_gate"],"UNKNOWN")
            self.assertIsNone(
                r["full_Zee_gauge_amplitude_upper_bound_GeV_minus2"])
            self.assertTrue(math.isfinite(r["retained_partial_BR"]))
            self.assertGreaterEqual(r["retained_partial_BR"],0.)
        for channel in ("mu_to_e_gamma","tau_to_e_gamma","tau_to_mu_gamma"):
            x=[r for r in output["rows"] if r["channel"]==channel]
            self.assertEqual(len(x),2)
            # F enters the computed dipole quadratically for m_nu=0.
            self.assertAlmostEqual(x[1]["known_F_one_loop_norm_GeV_minus2"]/
                                   x[0]["known_F_one_loop_norm_GeV_minus2"],
                                   4.,places=12)
            self.assertAlmostEqual(x[0]["Z12_over_e_abs"],
                                   x[1]["Z12_over_e_abs"],places=13)

    def test_refuse_closed_cascade_or_zero_singlet(self):
        p=sample()
        g=grid()
        g["mH2_GeV"]=[620.]
        with self.assertRaises(ValueError):
            list(grid_points(p,inversion(),g))
        g=grid()
        g["gamma_rad"]=[0.]
        with self.assertRaises(ValueError):
            list(grid_points(p,inversion(),g))
        p=sample()
        p["F"]={"real":np.zeros((3,3)).tolist(),
                "imag":np.zeros((3,3)).tolist()}
        with self.assertRaises(ValueError):
            list(grid_points(p,inversion(),grid()))

    def test_do_not_reuse_stale_potential_or_bad_inversion(self):
        p=sample()
        p["potential"]={"Mu4":0.}
        with self.assertRaises(ValueError):
            list(grid_points(p,inversion(),grid()))
        p=sample()
        bad=inversion()
        bad["K1"]=float("nan")
        with self.assertRaises(ValueError):
            list(grid_points(p,bad,grid()))
        bad=inversion()
        bad.pop("lambda_c")
        with self.assertRaises(ValueError):
            list(grid_points(p,bad,grid()))


if __name__=="__main__":
    unittest.main()
