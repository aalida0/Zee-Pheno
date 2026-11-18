"""Source-equation comparison: Zee Higgs basis and f/2 conventions."""
import math
import unittest

from Model.model import zee_f


def mul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def transpose(a):
    return tuple(zip(*a))


def combine(a, b, coefficient=1):
    return tuple(tuple(a[i][j]+coefficient*b[i][j] for j in range(3)) for i in range(3))


def scale(a, z):
    return tuple(tuple(z*v for v in row) for row in a)


class PublishedFormulaConvention(unittest.TestCase):
    def test_2017_eq230_2023_eq4_and_sarah_half_f_agree(self):
        # Nonzero imaginary and off-diagonal Yukawas prevent accidental
        # agreement from a real/diagonal-only control point.
        beta, gamma, v = 1.24, -0.71, 246.21965
        m1, m2 = 300., 700.
        masses = (.000511, .105658, 1.77686)
        mass = tuple(tuple(masses[i] if i == j else 0. for j in range(3))
                     for i in range(3))
        mass2 = mul(mass, mass)
        f_sarah = zee_f(1e-4+3e-5j, -4e-5+8e-6j, 2e-5-9e-6j)
        f_paper = scale(f_sarah, .5)  # 2017 vertex is 2 f; SPheno is Fzee.
        ye2 = ((.01+.001j, .002-.004j, .0009+.0004j),
               (-.002+.003j, .025-.002j, .003-.001j),
               (.0011-.0008j, -.004+.002j, .015+.001j))
        c, s = math.cos(beta), math.sin(beta)
        sm = scale(mass, math.sqrt(2)/v)
        ye = scale(combine(sm, ye2, -s), 1/c)
        d = combine(scale(ye, -s), scale(ye2, c))
        self.assertLess(max(abs(a-b) for x,y in zip(combine(scale(ye,c),scale(ye2,s)),sm)
                            for a,b in zip(x,y)),1e-15)
        kern = math.sin(2*gamma)/(16*math.pi**2)*math.log((m1/m2)**2)
        paper = mul(mul(f_paper,mass),d)
        paper = scale(combine(paper,transpose(paper)),kern)
        native = mul(mul(f_sarah,mass),d)
        native = scale(combine(native,transpose(native)),kern/2)
        # 2017 eqs. (2.30, 2.31), with the generic-basis Ye2=Y2.
        a2017 = (math.sin(2*gamma)*math.tan(beta)/(8*math.sqrt(2)*math.pi**2*v)
                 *math.log((m2/m1)**2))
        fm2 = mul(f_paper,mass2)
        fmye2 = mul(mul(f_paper,mass),ye2)
        bridge = combine(combine(fm2,transpose(fm2)),
                         combine(fmye2,transpose(fmye2)),
                         -v/(math.sqrt(2)*s))
        old = scale(bridge,a2017)
        for target in (native, old):
            self.assertLess(max(abs(a-b) for x,y in zip(paper,target)
                                for a,b in zip(x,y)),1e-17)
        self.assertLess(max(abs(paper[i][j]-paper[j][i])
                            for i in range(3) for j in range(3)),1e-17)

    def test_expected_zero_limits(self):
        # The dominant mass vanishes if the charged scalar mixing or mass
        # splitting vanishes. This does not test the finite-lepton loop terms.
        m1, m2 = 300., 700.
        def coeff(gamma, a, b):
            return math.sin(2*gamma)/(32*math.pi**2)*math.log((a/b)**2)
        self.assertEqual(coeff(0.,m1,m2),0.)
        self.assertEqual(coeff(0.52,m1,m1),0.)


if __name__ == '__main__':
    unittest.main()
