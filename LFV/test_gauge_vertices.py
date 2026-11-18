"""Tree-level Zee charged electroweak current constraints."""
import math
import unittest

import numpy as np

from LFV.gauge_vertices import charged_currents


class GaugeCurrentTests(unittest.TestCase):
    def test_2hdm_decoupling_and_singlet_charge(self):
        sw2=.23129
        q=charged_currents(0.,sw2)
        g=q["Z_over_e"]
        fac=math.sqrt(sw2*(1-sw2))
        np.testing.assert_allclose(g,np.diag([(.5-sw2)/fac,-sw2/fac]))
        np.testing.assert_array_equal(q["photon_over_e"],np.eye(2))
        self.assertAlmostEqual(g[0,1],0.,places=15)
        np.testing.assert_array_equal(q["W_doublet_projector"],np.diag([1.,0.]))

    def test_mixed_current_off_diagonal_and_sum_rules(self):
        sw2=.23129
        for gamma in (.15,.7,-.6):
            with self.subTest(gamma=gamma):
                q=charged_currents(gamma,sw2)
                r=q["charged_rotation"]
                s=q["singlet_fractions"]
                d=q["doublet_fractions"]
                self.assertTrue(np.allclose(r@r.T,np.eye(2),atol=1e-15))
                self.assertAlmostEqual(d@s,0.,places=15)
                np.testing.assert_allclose(np.outer(d,d)+np.outer(s,s),np.eye(2),
                                           atol=1e-15)
                np.testing.assert_allclose(
                    q["Z_over_e"],
                    r@np.diag([(.5-sw2)/math.sqrt(sw2*(1-sw2)),
                               -sw2/math.sqrt(sw2*(1-sw2))])@r.T,atol=1e-15)
                offdiag=-math.sin(2*gamma)/(4*math.sqrt(sw2*(1-sw2)))
                self.assertAlmostEqual(q["Z_over_e"][0,1],offdiag,places=14)
                self.assertEqual(q["status"],"TREE_EW_VERTEX_IDENTITIES_ONLY")
                self.assertAlmostEqual(np.trace(q["T3_mass_basis"]),.5)

    def test_bad_weak_angle_and_nans_rejected(self):
        for gamma,sw2 in ((float("nan"),.23),(0.,0.),(0.,1.),
                          (0.,-1.),(0.,float("inf"))):
            with self.assertRaises(ValueError):
                charged_currents(gamma,sw2)


if __name__=="__main__":
    unittest.main()
