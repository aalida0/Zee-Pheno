"""Mathematical doublet-charge projection from mixed Zee H1/H2."""
import copy
import math
import unittest

import numpy as np

from LFV.doublet_proxy_study import doublet_proxy
from LFV.nonzero_mixing_study import grid_points
from LFV.test_constraints import sample
from LFV.test_nonzero_mixing_study import inversion,grid


class DoubletProxyTests(unittest.TestCase):
    def test_exact_unmixed_doublet_mass_and_quartics(self):
        g=grid()
        g["gamma_rad"]=[.95]
        p,_=next(grid_points(sample(),inversion(),g))
        initial=copy.deepcopy(p)
        ref,info=doublet_proxy(p)
        self.assertEqual(p,initial)
        self.assertEqual(ref["scalar"]["gamma"],0.)
        self.assertEqual(ref["potential"]["Mu4"],0.)
        self.assertEqual(ref["potential"]["K1"],0.)
        self.assertEqual(ref["potential"]["K2"],0.)
        self.assertEqual(ref["potential"]["K12"],0.)
        c,s=math.cos(.95),math.sin(.95)
        self.assertAlmostEqual(ref["scalar"]["mH1"]**2,
            c*c*p["scalar"]["mH1"]**2+s*s*p["scalar"]["mH2"]**2,places=7)
        for k in ("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"):
            self.assertAlmostEqual(ref["potential"][k],p["potential"][k])
        self.assertLess(info["doublet_mass_tree_consistency_residual"],1e-9)
        self.assertEqual(ref["mass_yukawa_scheme"],
                         "AUTHOR_2HDM_CONTROL_WITH_POLE_LEPTONS_NOT_ZEE")
        self.assertEqual(ref["internal_lepton_masses_GeV"],
                         p["pole_lepton_masses_GeV"])
        self.assertFalse(any(any(row) for row in ref["F"]["real"]))

    def test_fails_on_inconsistent_or_zero_mixed_input(self):
        p,_=next(grid_points(sample(),inversion(),grid()))
        bad=copy.deepcopy(p)
        bad["potential"]["Lambda5"]+=.1
        with self.assertRaisesRegex(ValueError,"consistent 2HDM doublet"):
            doublet_proxy(bad)
        bad=copy.deepcopy(p)
        bad["scalar"]["gamma"]=0.
        with self.assertRaises(ValueError):
            doublet_proxy(bad)
        bad=copy.deepcopy(p)
        bad["F"]={"real":np.zeros((3,3)).tolist(),
                  "imag":np.zeros((3,3)).tolist()}
        with self.assertRaises(ValueError):
            doublet_proxy(bad)


if __name__=="__main__":
    unittest.main()
