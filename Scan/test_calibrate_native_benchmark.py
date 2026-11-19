"""Checks against actual calibrated native spectra, not fabricated fitted masses."""
import json
import math
import unittest
from pathlib import Path

import numpy as np

from Model.audit_neutrino_native import parse_entries, sha
from Model.check_native import blocks, fix_card
from Model.diagnose_neutrino_native import analyze
from Neutrino.fit import complex_matrix, inspect_native, target, zee_matrix
from Scan.calibrate_native_benchmark import (
    THEORY, failure_causes, fit_pass, potential, refit_update, running,
    snapshot, verify_running_inputs)
from Scan.diagnose_scalar_spectrum import potential_hessian
from Scan.physical_proposals import cheap_scalar
from Scan.native_gates import evaluate_entry

ROOT = Path(__file__).parent/'reference/native_benchmark'


class NativeBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=json.loads((ROOT/'ledger.json').read_text())
        cls.row=next(r for r in cls.report['original_scalar_proposals']
                     if r['scalar_id']=='s00017')

    def test_native_final_cards_are_source_bound_and_both_orderings_pass(self):
        for b in self.report['benchmarks']:
            d=ROOT/b['directory']
            for name,expected in b['artifacts_sha256'].items():
                self.assertEqual(sha(d/name),expected)
            r=json.loads((d/'native.json').read_text())
            self.assertEqual(sha(d/'input.in'),r['input_sha256'])
            self.assertTrue(fit_pass(b['observables']))
            self.assertEqual(inspect_native(d/'native.json',b['ordering'])['status'],
                             'NATIVE_3SIGMA_SCREEN_PASS')
            s={k:complex_matrix(v) if isinstance(v,dict) else np.asarray(v)
               for k,v in snapshot(r).items()}
            self.assertLess(analyze(s)['matrix_relative_residual'],1e-10)
            verify_running_inputs((d/'input.in').read_text(),r)
            self.assertEqual(b['production_control'],'PASS')
            self.assertEqual(b['production_max_mass_difference_GeV'],0.)
            self.assertTrue(b['production_spectrum_identical'])
            self.assertIsNone(b['globally_allowed'])

    def test_scalar_calibration_preserves_other_coordinates_and_theory(self):
        mh=self.report['calibration']['tree_mass_GeV']
        p=potential(self.row,mh)
        self.assertEqual(p,self.report['calibration']['potential'])
        for gate in THEORY:
            self.assertEqual(cheap_scalar(p)[gate],'PASS')
        for name in ('Lambda4','Lambda5','LambdaC','K1','K2','K12','Mu4','MCS2','M12Sq'):
            self.assertEqual(p[name],self.row['lagrangian'][name])
        beta=self.row['scalar']['beta'];v=self.row['scalar']['v']
        names=dict(M12S='M12Sq',MC2='MCS2',Mu4='Mu4',M112='M112',M222='M222',
                   **{f'Lam{i}':f'Lambda{i}' for i in range(1,6)},
                   LamC='LambdaC',K1='K1',K2='K2',K12='K12')
        native={k:p[v] for k,v in names.items()}
        native.update(vd=v*math.cos(beta),vu=v*math.sin(beta))
        _,h=potential_hessian(native)
        np.testing.assert_allclose(np.linalg.eigvalsh(h[:2,:2]),
            [mh**2,self.row['scalar']['mh2']**2],rtol=1e-11,atol=1e-6)
        for b in self.report['benchmarks']:
            text=(ROOT/b['directory']/'input.in').read_text()
            flags=dict(parse_entries(text,'SPHENOINPUT'))
            self.assertEqual(flags[19],-2.)
            self.assertEqual([x for k,x in parse_entries(text,'SPHENOINPUT') if k==55],[1.,1.])
            self.assertEqual(flags[67],1.)

    def test_refit_uses_measured_running_direction_and_full_complex_inverse(self):
        for b in self.report['benchmarks']:
            start='higgs_03' if b['ordering']=='NO' else 'IO_start'
            r=json.loads((ROOT/(start+'.json')).read_text())
            family=b['selected_family']
            update,proposal=refit_update(r,family)
            self.assertLess(proposal['forward_residual'],1e-10)
            cb,sb,el,er,_,_=running(snapshot(r))
            f,ye2=complex_matrix(proposal['F_gauge']),complex_matrix(proposal['Ye2_gauge'])
            s=snapshot(r)
            sy=cb*complex_matrix(s['Ye'])+sb*complex_matrix(s['Ye2'])
            d=(ye2-sb*sy)/cb
            matrix=zee_matrix(el.conj()@f@el.conj().T,er.conj()@d@el.conj().T,
                              s['MFe'],s['MHm'][1:],r['matched_comparison']['sin2gamma_ZP'])
            expected,_,_=target(b['ordering'],family['mlight_eV'],
                               family['alpha21_rad'],family['alpha31_rad'])
            self.assertLess(np.linalg.norm(matrix-expected)/np.linalg.norm(expected),1e-10)
            self.assertGreater(np.max(abs(f.imag)),0.)
            self.assertGreater(np.max(abs(ye2.imag)),0.)
            self.assertEqual(len(update),36)

    def test_failure_statistics_and_quantitative_matching_cause_are_retained(self):
        index=json.loads((ROOT/'execution_statistics.json').read_text())
        self.assertEqual(index['total_native_executions'],len(index['executions']))
        self.assertEqual(index['total_native_executions'],104)
        self.assertEqual(index['total_native_spectrum_failures'],0)
        self.assertEqual(self.report['statistics']['baseline'],{'TARGET_FAIL':16})
        self.assertEqual(len(self.report['failure_causes']),16)
        self.assertEqual(self.report['statistics']['neutrino_refit'],
                         {'TARGET_FAIL':2,'TARGET_PASS':2})
        r=json.loads((ROOT/'s00017_NO_f00_baseline.json').read_text())
        ref=json.loads((ROOT.parents[2]/'Model/reference/native_neutrino_snapshot.json').read_text())['snapshot']
        family=self.report['benchmarks'][0]['selected_family']
        # The source family F,D values are reconstructed by the original inverse.
        from Neutrino.fit import inverse, encode
        target_matrix,_,_=target('NO',family['mlight_eV'],family['alpha21_rad'],family['alpha31_rad'])
        fm,dm,_=inverse(target_matrix,ref['MFe'],
            [self.row['scalar']['mH1'],self.row['scalar']['mH2']],
            math.sin(2*self.row['scalar']['gamma']),family['f_scale'],complex(*family['t']),
            root=family['root'],chart=family['chart'],f_phase=family['f_phase_rad'])
        family=dict(family,F_mass_basis=encode(fm),D_mass_basis=encode(dm))
        cause=failure_causes(self.row,family,r,ref)
        self.assertLess(cause['kernel_only_target_residual'],.004)
        self.assertGreater(cause['fully_matched_target_residual'],4.)
        self.assertLess(cause['D_matching_identity_abs_residual'],1e-15)

    def test_boundary_and_fit_failures_are_not_promoted(self):
        b=self.report['benchmarks'][0];d=ROOT/b['directory']
        r=json.loads((d/'native.json').read_text());text=(d/'input.in').read_text()
        value=blocks(text)['YEE2IN'][2,2]
        altered=fix_card(text,{('YEE2IN',(2,2)):value+1e-6})
        with self.assertRaisesRegex(ValueError,'Ye2 propagation'):
            verify_running_inputs(altered,r)
        failed=dict(b['observables'],target_matrix_relative_residual=.01)
        self.assertFalse(fit_pass(failed))
        self.assertFalse(fit_pass(dict(b['observables'],Higgs_mass_GeV=126.)))
        self.assertFalse(fit_pass(dict(b['observables'],cascade_status='FAIL')))
        self.assertTrue(all(status=='UNKNOWN' for status in self.report['outstanding_gates'].values()))

    def test_downstream_pipeline_accepts_native_row_phases_without_global_promotion(self):
        for b in self.report['benchmarks']:
            d=ROOT/b['directory']
            result=evaluate_entry(dict(point_id=b['directory'],ordering=b['ordering'],
                card=str((d/'input.in').resolve()),native=str((d/'native.json').resolve()),
                spectrum=str((d/'spectrum.slha').resolve())),ROOT)
            self.assertEqual(result['verdict'],'UNKNOWN_PENDING_GATES')
            self.assertIsNone(result['globally_allowed'])
            self.assertEqual(result['LFV_full_gate'],'UNKNOWN')
            self.assertTrue(all(result['gates'][g]['status']=='PASS' for g in
                ('neutrino_1D_3sigma','h125_mass_screen','cascade_on_shell')))


if __name__=='__main__':
    unittest.main()
