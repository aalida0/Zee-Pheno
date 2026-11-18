"""Tree-level Zee charged electroweak current constraints."""
import math
import unittest

import numpy as np

from LFV.gauge_vertices import charged_currents, charged_lepton_vertices


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

    def test_native_charged_lepton_F_D_vertices_and_completeness(self):
        F=np.array([[0,2e-5,3e-5],[-2e-5,0,4e-5],
                    [-3e-5,-4e-5,0]],complex)
        D=np.array([[.002,.0001j,0],[.0002,.004,.0003j],
                    [0,0,.01]],complex)
        for angle in (0.,.37,-.6):
            v=charged_lepton_vertices(F,D,angle)
            L,R=v["L"],v["R"]
            self.assertEqual(L.shape,(2,3,3))
            self.assertEqual(R.shape,(2,3,3))
            np.testing.assert_allclose(
                np.sum(L.conj()*L,axis=0),abs(F)**2,atol=1e-22)
            np.testing.assert_allclose(
                np.sum(R.conj()*R,axis=0),abs(D.conj().T)**2,atol=1e-22)
            np.testing.assert_allclose(
                np.sum(L.conj()*R,axis=0),np.zeros((3,3)),atol=1e-22)
            if angle==0.:
                np.testing.assert_array_equal(L[0],np.zeros((3,3)))
                np.testing.assert_array_equal(R[1],np.zeros((3,3)))
                np.testing.assert_array_equal(L[1],-F)
                np.testing.assert_array_equal(R[0],-D.conj().T)
        with self.assertRaises(ValueError):
            charged_lepton_vertices(F+np.eye(3),D,0.)

    def test_bad_weak_angle_and_nans_rejected(self):
        for gamma,sw2 in ((float("nan"),.23),(0.,0.),(0.,1.),
                          (0.,-1.),(0.,float("inf"))):
            with self.assertRaises(ValueError):
                charged_currents(gamma,sw2)


if __name__=="__main__":
    unittest.main()
