"""Scope-limited identities; NOT native-SPheno or LFV tests."""
import unittest
import numpy as np
from Higgs.neutral_lfv_tree import (neutral_yukawas,inclusive_pair_width,
                                     fermion_width,H_hh_trilinear,included_widths)

class TreeNeutralTests(unittest.TestCase):
    def test_alignment_h_lfv_vanishes(self):
        D=np.ones((3,3),complex)*.02
        y=neutral_yukawas(D,.4,.4)
        self.assertEqual(inclusive_pair_width(125,y['h'],1,2),0.)
        self.assertGreater(inclusive_pair_width(600,y['H'],1,2),0.)
    def test_yukawa_scaling_quadratic(self):
        D=np.zeros((3,3),complex);D[1,2]=.01+0.02j
        y=neutral_yukawas(D,.4,.4)['H']
        self.assertAlmostEqual(inclusive_pair_width(600,2*y,1,2)/inclusive_pair_width(600,y,1,2),4)
    def test_thresholds(self):
        self.assertEqual(fermion_width(300,1.,172.4,3,'scalar'),0.)
        self.assertGreater(fermion_width(600,1.,172.4,3,'pseudo'),fermion_width(600,1.,172.4,3,'scalar'))
    def test_partial_widths_nonnegative_and_unknown(self):
        p=dict(mh1=125.,mh2=603.6039389818,ma=622.40011586555,alpha=.6088373713,
               beta=.5963066963,v=246.21965)
        quartics=dict(Lambda1=.6483415613,Lambda2=1.4216039111,Lambda3=1.0961328959,
                      Lambda4=-.5137425074,Lambda5=-1.0145961353)
        D=np.zeros((3,3),complex);D[1,2]=.01
        r=included_widths(p,quartics,D)
        for boson in r:
            self.assertEqual(r[boson]['acceptance'],'UNKNOWN')
            self.assertTrue(all(x>=0 for x in r[boson]['partial_widths_GeV'].values()))
            self.assertLess(r[boson]['BR_mu_tau_optimistic_upper'],1.)
if __name__=='__main__':unittest.main()
