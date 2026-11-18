"""Independent tree-level spectrum closure tests; no SPheno claims."""
import math
import unittest
from Model.model import scalar_parameters, zee_f


def eigenvalues2(a, b, c):
    delta = math.sqrt((a-c)**2 + 4*b*b)
    return sorted(((a+c-delta)/2, (a+c+delta)/2))


class ModelTests(unittest.TestCase):
    def test_complex_f_is_antisymmetric(self):
        f = zee_f(1e-4+2e-4j, -3e-5+4e-5j, 5e-6-1e-6j)
        for i in range(3):
            for j in range(3):
                self.assertEqual(f[i][j], -f[j][i])
        with self.assertRaises(ValueError):
            zee_f(float('nan'), 0, 0)

    def test_scalar_inverse_closes_neutral_charged_and_pseudoscalar_masses(self):
        x = dict(mh1=125., mh2=580.7, ma=633.7, mH1=300., mH2=700.,
                 alpha=1.398, beta=1.394, gamma=1.089, m12sq=57700.,
                 lc=4.473, k1=1.082, k2=.00398, k12=-.001266, v=246.21965)
        p = scalar_parameters(**x)
        sb, cb = math.sin(x['beta']), math.cos(x['beta'])
        v2 = x['v']**2
        l345 = p['Lambda3']+p['Lambda4']+p['Lambda5']
        neutral = eigenvalues2(
            p['M12Sq']*sb/cb + v2*cb**2*p['Lambda1'],
            -p['M12Sq'] + l345*v2*sb*cb,
            p['M12Sq']*cb/sb + v2*sb**2*p['Lambda2'])
        charged = eigenvalues2(
            p['M12Sq']/(sb*cb) - v2/2*(p['Lambda4']+p['Lambda5']),
            -p['Mu4']*x['v']/math.sqrt(2),
            p['MCS2'] + v2/2*(p['K1']*cb**2 + p['K2']*sb**2-2*p['K12']*sb*cb))
        for found, expected in zip(neutral, [x['mh1']**2,x['mh2']**2]):
            self.assertAlmostEqual(found/expected, 1., places=9)
        for found, expected in zip(charged, [x['mH1']**2,x['mH2']**2]):
            self.assertAlmostEqual(found/expected, 1., places=9)
        ma2 = p['M12Sq']/(sb*cb) - p['Lambda5']*v2
        self.assertAlmostEqual(ma2 / x['ma']**2, 1., places=9)

    def test_tadpoles_and_charged_goldstone(self):
        x = dict(mh1=125., mh2=580.7, ma=633.7, mH1=300., mH2=700.,
                 alpha=1.398, beta=1.394, gamma=1.089, m12sq=57700.,
                 lc=4.473, k1=1.082, k2=.00398, k12=-.001266)
        p = scalar_parameters(**x)
        v = 246.21965
        cb, sb = math.cos(x["beta"]), math.sin(x["beta"])
        v1, v2 = v*cb, v*sb
        l345 = p["Lambda3"]+p["Lambda4"]+p["Lambda5"]
        t1 = p["M112"]*v1-p["M12Sq"]*v2 + (
            p["Lambda1"]*v1**3+l345*v1*v2**2)/2
        t2 = p["M222"]*v2-p["M12Sq"]*v1 + (
            p["Lambda2"]*v2**3+l345*v2*v1**2)/2
        self.assertLess(max(abs(t1), abs(t2)), 1e-6)
        # Full charged 3x3 matrix, independent of the projected 2x2 test.
        d = p["M12Sq"]/(sb*cb)-v*v*(p["Lambda4"]+p["Lambda5"])/2
        mixing = p["Mu4"]*v/math.sqrt(2)
        charged = ((d*sb*sb, -d*sb*cb, mixing*sb),
                   (-d*sb*cb, d*cb*cb, -mixing*cb),
                   (mixing*sb, -mixing*cb, p["MCS2"]))
        goldstone = (cb, sb, 0.)
        for row in charged:
            self.assertLess(abs(sum(a*b for a, b in zip(row, goldstone))), 1e-7)

    def test_signed_gamma_flips_only_charged_offdiagonal(self):
        args = dict(mh1=125.,mh2=580.7,ma=633.7,mH1=300.,mH2=700.,
                    alpha=1.398,beta=1.394,m12sq=57700.,lc=4.473,
                    k1=1.082,k2=.00398,k12=-.001266)
        p = scalar_parameters(gamma=1.089,**args)
        n = scalar_parameters(gamma=-1.089,**args)
        self.assertAlmostEqual(p['Mu4'], -n['Mu4'])
        for k in p:
            if k != 'Mu4':
                self.assertAlmostEqual(p[k], n[k], places=9)
        with self.assertRaises(ValueError):
            scalar_parameters(gamma=math.pi/2,**args)


if __name__ == '__main__':
    unittest.main()
