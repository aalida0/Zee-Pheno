"""Lepton basis matching check with a complex, off-diagonal second doublet."""
import math
import unittest
from Model.check_basis import inspect


def output(broken=False):
    beta = 1.24
    c, s = math.cos(beta), math.sin(beta)
    masses = (.000511, .105658, 1.77686)
    v = 246.21965
    ye2 = [[complex(.002*(i+1),.0005*j) if i != j else complex(.015*(i+1), .001)
            for j in range(3)] for i in range(3)]
    ye = [[(complex(math.sqrt(2)*masses[i]/v) if i==j else 0j)-s*ye2[i][j]
            for j in range(3)] for i in range(3)]
    ye = [[z/c for z in row] for row in ye]
    if broken: ye[0][1] += .01
    f = [[0j, 1e-4+2e-4j, -3e-5+4e-5j],
         [-1e-4-2e-4j,0j,5e-5-2e-5j],
         [3e-5-4e-5j,-5e-5+2e-5j,0j]]
    text = 'BLOCK MINPAR\n 10 %.17e # tanbeta\n' % math.tan(beta)
    for name, arr in [('YE',ye),('YEE2',ye2),('FZEE',f)]:
        for imag in (False,True):
            text+='Block '+('IM' if imag else '')+name+'\n'
            for i,row in enumerate(arr):
                for j,z in enumerate(row):
                    text+=f' {i+1} {j+1} {(z.imag if imag else z.real):.17e} # entry\n'
    return text


class BasisTests(unittest.TestCase):
    def test_mass_basis_with_general_ye2(self):
        r=inspect(output())
        self.assertEqual(r['status'],'TREE_LEPTON_BASIS_PASS')
        self.assertLess(r['mass_yukawa_offdiag_relative'],1e-13)
        self.assertEqual(r['radiative_neutrino_mass'],'NOT_EVALUATED')
        self.assertEqual(r['pmns_from_tree_UV'],'NOT_ESTABLISHED')

    def test_detects_charged_lepton_misalignment(self):
        self.assertEqual(inspect(output(broken=True))['status'],'TREE_LEPTON_BASIS_FAIL')


if __name__=='__main__':
    unittest.main()
