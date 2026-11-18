"""Regression against the actual generated native complex matrix, not mock masses."""
import json
import unittest
from pathlib import Path

import mpmath as mp
import numpy as np

from Model.compare_neutrino_masses import finite_kernel
from Model.diagnose_neutrino_native import analyze, instrument


class NativeMatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads((Path(__file__).parent/'reference/native_neutrino_snapshot.json').read_text())
        cls.snapshot = {k: np.asarray(v['real'])+1j*np.asarray(v['imag'])
                        if isinstance(v, dict) else np.asarray(v)
                        for k, v in cls.reference['snapshot'].items()}

    def test_actual_native_matrix_and_takagi(self):
        result = analyze(self.snapshot)
        for key in ('matrix_relative_residual', 'vertex_relative_residual',
                    'majorana_symmetry_relative_residual', 'takagi_relative_residual',
                    'charged_lepton_mass_matrix_relative_residual'):
            self.assertLess(result[key], 1e-11, key)
        np.testing.assert_allclose(result['matched_analytic_masses_eV'],
                                   self.reference['spheno_pole_masses_eV'], rtol=1e-8)
        self.assertAlmostEqual(result['sin2gamma_ZP'], result['sin2gamma_Mu4'], places=13)

    def test_original_unmatched_matrix_really_fails(self):
        old = self.reference['original_comparison']
        a = old['analytic_Mnu_eV']
        original = np.array(a['real'])+1j*np.array(a['imag'])
        native = self.snapshot['mat1']*1e9
        self.assertGreater(np.max(abs(original-native))/np.max(abs(native)), .5)
        self.assertGreater(old['max_abs_mass_difference_over_analytic_max'], .56)
        self.assertFalse(old['matched_native_running_parameter_scheme'])

    def test_complex_congruence_preserves_singular_values(self):
        rng = np.random.default_rng(20261118)
        u, _ = np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
        m = self.snapshot['mat1']
        np.testing.assert_allclose(np.linalg.svd(u@m@u.T, compute_uv=False),
                                   np.linalg.svd(m, compute_uv=False), rtol=1e-13)

    def test_finite_kernel_against_direct_quadrature(self):
        m1, m2 = self.snapshot['MHm'][1:]
        for ml in self.snapshot['MFe']:
            with mp.workdps(45):
                a, b, z = (mp.mpf(float(x))**2 for x in (m1, m2, ml))
                integral = mp.quad(lambda x: mp.log((x*z+(1-x)*a)/(x*z+(1-x)*b)), [0, 1])
            self.assertAlmostEqual(finite_kernel(m1, m2, ml), float(integral), places=13)
        self.assertAlmostEqual(finite_kernel(m1, m2, 0), np.log((m1/m2)**2), places=14)

    def test_instrumentation_rejects_unrecognized_source(self):
        with self.assertRaises(ValueError):
            instrument('not the generated source')


if __name__ == '__main__':
    unittest.main()
