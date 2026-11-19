"""Rebuilt-backend provenance falsifiers and regressions on actual HEPfit logs.

CI checks archived numerical evidence; it does not compile or execute HEPfit.
"""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Flavor.bsgamma import (HEPFIT_COMMIT, backend, injected_coeff_check,
                            native_scalar, parse_driver_output, sha)
from Flavor.validate_hepfit import compare_runs, verify_attestation

DATA=Path(__file__).parent/'reference/hepfit_restored'


class RestoredBackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.replay=json.loads((DATA/'historical-replay/HEPFIT_EXECUTION.json').read_text())
        cls.attestation=json.loads((DATA/'RESTORED_BACKEND.json').read_text())
        cls.points=json.loads((DATA/'s00017/bsgamma_nlo.json').read_text())

    def test_historical_sm_d7_f13_and_new_build_binding(self):
        r=self.attestation
        self.assertEqual(r['validated_runs'],37)
        self.assertEqual(r['maximum_output_absolute_residual'],0.)
        self.assertEqual(r['HEPfit_commit'],HEPFIT_COMMIT)
        build=json.loads((DATA/'BUILD_RECEIPT.json').read_text())
        self.assertEqual(r['binary_sha256'],build['artifacts']['/home/amine/zee-hepfit-build/driver'])
        self.assertEqual(r['source_lock_sha256'],sha(DATA/'SOURCE_LOCK.json'))
        self.assertEqual(r['execution_receipt_sha256'],sha(DATA/'historical-replay/HEPFIT_EXECUTION.json'))
        for name,expected in [('SM_reference',3.491423243740102e-4),
                ('d7_NO_neutral_near_heavy_07_central',4.17955384241e-4),
                ('F13_heavy_0001330_central',3.83711422745e-4)]:
            self.assertAlmostEqual(r['original_central_BR'][name],expected,places=15)
        self.assertIsNone(r['theory_sigma'])
        self.assertIsNone(r['experimental_exclusion'])
        self.assertEqual(r['physical_bsgamma'],'UNKNOWN')

    def test_actual_37_run_logs_inputs_initialization_and_operator_map(self):
        directory=DATA/'historical-replay'
        for name,r in self.replay['runs'].items():
            self.assertEqual(sha(directory/r['stdout']),r['stdout_sha256'])
            self.assertEqual(sha(directory/(name+'.txt')),r['input_sha256'])
            raw=(directory/r['stdout']).read_text()
            v=parse_driver_output(raw,matching_scale=r['input']['muw'])
            self.assertEqual(v,r['values'])
            self.assertEqual(v['cut_GeV'],1.6)
            if r['mode']!='sm':
                self.assertEqual(v['initialized_unused_LO_C1'],0.)
                self.assertEqual(v['initialized_unused_NNLO_C1'],0.)
            if r['mode']=='bridge':
                np=[[r['input'].get(f'NP{o}_{i+1}',0.) for i in range(8)] for o in range(2)]
                expected=injected_coeff_check(v,np)
                for key,value in expected.items():
                    self.assertLessEqual(abs(v[key]-value),2e-14)
            if r['input']['muw']!=160.:
                with self.assertRaisesRegex(ValueError,'matching scale'):
                    parse_driver_output(raw)

    def test_native_predictions_sm_zero_np_and_no_uncertainty_promotion(self):
        directory=DATA/'s00017'
        expected={'NO':3.9214198285658576e-4,'IO':3.921419354070257e-4}
        for p in self.points['points']:
            name=p['ordering']
            self.assertEqual(p['scalar'],native_scalar(name))
            values={label:parse_driver_output((directory/(name+'_'+label+'.log')).read_text())
                    for label in ('SM','zero_NP','Zee_NLO')}
            self.assertEqual(values['SM']['BR_NLO'],values['zero_NP']['BR_NLO'])
            self.assertEqual(values['zero_NP']['mt_MSbar_mut'],164.29247164143644)
            self.assertAlmostEqual(values['Zee_NLO']['BR_NLO'],expected[name],places=16)
            self.assertEqual(values['Zee_NLO']['BR_NLO'],p['BR_NLO_Egamma_gt_1p6GeV'])
            self.assertEqual(p['delta_BR_same_scheme'],p['BR_NLO_Egamma_gt_1p6GeV']-p['BR_SM_same_scheme'])
            self.assertEqual(p['historical_FlavorKit_BR'],3.63471789e-4)
            self.assertGreater(p['BR_NLO_Egamma_gt_1p6GeV'],p['BR_SM_same_scheme'])
            injected_coeff_check(values['Zee_NLO'],p['ordinary_NP_vector'])
            self.assertEqual(p['physical_bsgamma_gate'],'UNKNOWN')
        self.assertEqual(self.points['physical_bsgamma'],'UNKNOWN')
        self.assertIsNone(self.points['theory_sigma'])
        self.assertIsNone(self.points['experimental_exclusion'])

    def test_numeric_equivalence_rejects_changed_inputs_missing_runs_nan_and_coefficients(self):
        self.assertTrue(all(v==0 for v in compare_runs(self.replay,self.replay).values()))
        name='d7_NO_neutral_near_heavy_07_central'
        for kind in ('input','missing_run','nan','coefficient','rate'):
            changed=copy.deepcopy(self.replay)
            if kind=='input':changed['runs'][name]['input']['muw']=161.
            if kind=='missing_run':changed['runs'].pop(name)
            if kind=='nan':changed['runs'][name]['values']['BR_NLO']=float('nan')
            if kind=='coefficient':changed['runs'][name]['values']['high_1_1_7']+=1e-12
            if kind=='rate':changed['runs'][name]['values']['BR_NLO']+=1e-10
            with self.subTest(kind=kind),self.assertRaises(ValueError):
                compare_runs(self.replay,changed)

    def test_changed_binary_cannot_use_the_legacy_provenance_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'research/full_constraints/bsgamma_hepfit'
            (base/'execution').mkdir(parents=True);(root/'scripts').mkdir()
            match=root/'scripts/zee_module5_bmu_nlo.py';match.write_text('matching source')
            (base/'driver.cpp').write_text('original driver')
            (base/'SOURCE_LOCK.json').write_text(json.dumps({'HEPfit_commit':HEPFIT_COMMIT}))
            (base/'execution/HEPFIT_EXECUTION.json').write_text(json.dumps(dict(
                HEPfit_commit=HEPFIT_COMMIT,driver_sha256=sha(base/'driver.cpp'),binary_sha256='0'*64)))
            binary=root/'driver';binary.write_text('changed executable')
            with patch('Flavor.bsgamma.MATCHER_SHA',sha(match)):
                with self.assertRaisesRegex(ValueError,'requires a numerically validated'):
                    backend(root,binary,root/'source')

    def test_restored_attestation_rechecks_evidence_instead_of_trusting_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);build=root/'build.json';execution=root/'execution.json'
            build.write_text('build bytes');execution.write_text('execution bytes')
            record=dict(schema_version=1,status='RESTORED_HEPFIT_NUMERICALLY_EQUIVALENT',
                build_receipt=str(build),execution_receipt=str(execution),
                build_receipt_sha256=sha(build),execution_receipt_sha256=sha(execution))
            path=root/'attestation.json';path.write_text(json.dumps(record))
            with patch('Flavor.validate_hepfit.validate',return_value=record) as verify:
                self.assertEqual(verify_attestation(path,root,root/'driver',root/'source'),record)
                verify.assert_called_once()
            with patch('Flavor.validate_hepfit.validate',return_value=dict(record,binary_sha256='invalid')):
                with self.assertRaisesRegex(ValueError,'actual evidence'):
                    verify_attestation(path,root,root/'driver',root/'source')
            build.write_text('tampered build bytes')
            with self.assertRaisesRegex(ValueError,'receipt changed'):
                verify_attestation(path,root,root/'driver',root/'source')


if __name__=='__main__':
    unittest.main()
