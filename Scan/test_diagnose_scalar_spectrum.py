"""Regressions against witnessed native matrices and the independent potential."""
import json
import math
import unittest
from pathlib import Path

import numpy as np

from Scan.diagnose_scalar_spectrum import potential_hessian, tree_comparison, instrument_tree
from Scan.physical_proposals import scalar_updates


class ScalarNativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = json.loads((Path(__file__).parent/
            'reference/native_scalar_investigation.json').read_text())

    def test_failure_is_in_firstguess_with_reference_quartics_except_lambda3(self):
        cases = self.evidence['cases']
        ref = cases['reference']['first_tree'][0]
        failures = dict(s00001=-74437.055789404374,s00010=-99521.371844103094,
                        s00017=-35290.085776078282,s00020=-130954.31698470516)
        for sid,expected in failures.items():
            actual = cases[sid+'_original']['first_tree'][0]
            self.assertEqual(actual['stack'],'FirstGuess|TreeMassesZeeTypeII|CalculateMhh|')
            self.assertAlmostEqual(np.linalg.eigvalsh(actual['matrix_GeV2'])[0],expected,places=7)
            changed = {k for k in ref['inputs'] if actual['inputs'][k]!=ref['inputs'][k]}
            self.assertEqual(changed,{'Lam3','M112','M222'})
        self.assertEqual(cases['s00001_original']['first_tree'],
                         cases['s00001_reference_yukawa']['first_tree'])

    def test_requested_spectrum_from_direct_potential_hessian(self):
        for row in self.evidence['scalar_proposals']:
            x,p = row['scalar'],row['lagrangian']
            beta,alpha,gamma = (x[k] for k in ('beta','alpha','gamma'))
            cb,sb = math.cos(beta),math.sin(beta)
            names = dict(M12S='M12Sq',MC2='MCS2',Mu4='Mu4',M112='M112',M222='M222',
                         **{f'Lam{i}':f'Lambda{i}' for i in range(1,6)},
                         LamC='LambdaC',K1='K1',K2='K2',K12='K12')
            native = {k:p[v] for k,v in names.items()}
            native.update(vd=x['v']*cb,vu=x['v']*sb)
            gradient,h = potential_hessian(native)
            self.assertLess(np.max(abs(gradient)),2e-6)
            ca,sa = math.cos(alpha),math.sin(alpha)
            rotation = np.array([[ca,sa],[-sa,ca]])
            np.testing.assert_allclose(rotation @ h[:2,:2] @ rotation.T,
                np.diag([x['mh1']**2,x['mh2']**2]),rtol=1e-12,atol=1e-6)
            odd = np.array([-sb,cb])
            self.assertAlmostEqual(float(odd @ h[2:4,2:4] @ odd)/x['ma']**2,1,places=11)
            u = np.array([[-sb,cb,0],[0,0,1]])
            cg,sg = math.cos(gamma),math.sin(gamma)
            rotation = np.array([[cg,sg],[-sg,cg]])
            np.testing.assert_allclose(rotation @ u @ h[4:,4:] @ u.T @ rotation.T,
                np.diag([x['mH1']**2,x['mH2']**2]),rtol=1e-12,atol=1e-6)

    def test_native_tree_matrices_and_input_propagation(self):
        for row in self.evidence['scalar_proposals']:
            case = self.evidence['cases'][row['scalar_id']+'_NO_f00_corrected']
            self.assertEqual(tree_comparison(case['first_tree'])['status'],'PASS')
            native = case['first_tree'][0]['inputs']
            updates = scalar_updates(row['lagrangian'],math.tan(row['scalar']['beta']))
            for i in range(1,6):
                self.assertAlmostEqual(native[f'Lam{i}'],updates['MINPAR',(i,)],places=12)
            self.assertAlmostEqual(native['vu']/native['vd'],updates['MINPAR',(10,)],places=12)
            self.assertEqual(case['scalar_input_propagation']['status'],'PASS')
            final = case['scalar_input_propagation']['HMIX']
            self.assertAlmostEqual(final['103']/final['102'],updates['MINPAR',(10,)],places=7)

    def test_all_corrected_spectra_are_not_physical_passes(self):
        cases = {k:v for k,v in self.evidence['cases'].items()
                 if k.endswith('_corrected')}
        self.assertEqual(len(cases),16)
        for case in cases.values():
            self.assertEqual(case['spectrum_status'],'PASS')
            self.assertEqual(case['neutrino_status'],'PASS')
            self.assertLess(case['matched_comparison']['matrix_relative_residual'],1e-11)
            self.assertEqual(case['cascade_status'],'PASS')
            self.assertEqual(case['h125_mass_screen'],'FAIL')
            self.assertEqual(case['oscillation_screen']['status'],'NATIVE_3SIGMA_SCREEN_FAIL')

    def test_instrumentation_rejects_unrecognized_source(self):
        with self.assertRaises(ValueError):
            instrument_tree('not the generated source')


if __name__ == '__main__':
    unittest.main()
