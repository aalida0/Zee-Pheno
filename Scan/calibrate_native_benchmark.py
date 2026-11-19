"""Calibrate one scalar coordinate and refit F/Ye2 in the frozen native scheme.

This is an explicitly fitted benchmark, not a scan or a globally allowed point.
Every solver evaluation is retained. No SPheno input flag is changed in the fit.
"""
import argparse
import hashlib
import json
import math
import platform
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq

from Model.audit_neutrino_native import sha
from Model.check_native import blocks, fix_card
from Model.diagnose_neutrino_native import BINARY_SHA
from Model.model import scalar_parameters
from Neutrino.fit import (NUFIT, SOURCE, complex_matrix, encode, inspect_native,
                          inverse, target, zee_matrix)
from Scan.diagnose_scalar_spectrum import build, run_case
from Scan.physical_proposals import cheap_scalar, scalar_updates

HIGGS_TARGET = 125.25
TREE_BRACKET = (70., 125.25)  # GeV; fixed before evaluating the selected benchmark
HIGGS_ACCURACY = .02  # numerical calibration accuracy in GeV, not the Higgs screen
MATRIX_ACCURACY = 1e-8  # relative fit closure; not a claimed perturbative precision
MAX_REFITS = 6
THEORY = ('doublet_BFB_necessary', 'zee_quartic_BFB_sufficient',
          'quartic_lt_4pi_convention')
UNKNOWN = ('full_Zee_LFV', 'B_flavor', 'Higgs_rates', 'Higgs_direct_searches',
           'full_tree_unitarity', 'global_vacuum', 'EWPO', 'cosmology')


def potential(row, tree_mass):
    x = dict(row['scalar'], mh1=tree_mass)
    p = row['lagrangian']
    return scalar_parameters(**x, m12sq=p['M12Sq'], lc=p['LambdaC'],
                             k1=p['K1'], k2=p['K2'], k12=p['K12'])


def snapshot(result):
    i = result['spectrum_snapshot_index']
    if i is None:
        raise ValueError('No unique spectrum-associated native snapshot')
    return result['snapshots'][i]


def running(s):
    v = math.hypot(s['vd'], s['vu'])
    cb, sb = s['vd']/v, s['vu']/v
    el, er = complex_matrix(s['ZEL']), complex_matrix(s['ZER'])
    ye, ye2 = complex_matrix(s['Ye']), complex_matrix(s['Ye2'])
    fm = el.conj() @ complex_matrix(s['fzee']) @ el.conj().T
    dm = er.conj() @ (-sb*ye+cb*ye2) @ el.conj().T
    return cb, sb, el, er, fm, dm


def matrix_updates(f, ye2):
    return {((('IM' if imag else '')+name), (i+1,j+1)):
            float(matrix[i,j].imag if imag else matrix[i,j].real)
            for name,matrix in (('FZEEIN',f),('YEE2IN',ye2))
            for imag in (False,True) for i in range(3) for j in range(3)}


def refit_update(result, family):
    """Invert at this card's running masses, scalar mixing and lepton basis."""
    s = snapshot(result)
    cb,sb,el,er,_,_ = running(s)
    mnu,_,_ = target(family['ordering'],family['mlight_eV'],
                    family['alpha21_rad'],family['alpha31_rad'])
    fm,dm,residual = inverse(
        mnu,s['MFe'],s['MHm'][1:],result['matched_comparison']['sin2gamma_ZP'],
        family['f_scale'],complex(*family['t']),root=family['root'],
        chart=family['chart'],f_phase=family['f_phase_rad'],
        homogeneous=complex_matrix(family['homogeneous_coefficients']).reshape(4))
    # LOW boundary assignments reset fzee and Ye2 to their input matrices.
    # Ye is solved by native matching: preserve its measured SM mass direction.
    sy = cb*complex_matrix(s['Ye'])+sb*complex_matrix(s['Ye2'])
    fg, dg = el.T @ fm @ el, er.T @ dm @ el
    ye2 = sb*sy+cb*dg
    ye = cb*sy-sb*dg
    largest = max(float(np.max(abs(a))) for a in (fg,ye2,ye))
    return matrix_updates(fg,ye2), dict(
        forward_residual=residual, max_abs_yukawa=largest,
        f_scale=family['f_scale'], F_gauge=encode(fg), Ye2_gauge=encode(ye2))


def observables(result, directory, family):
    s = snapshot(result)
    _,_,el,_,_,_ = running(s)
    mnu,_,_ = target(family['ordering'],family['mlight_eV'],
                    family['alpha21_rad'],family['alpha31_rad'])
    native = el.conj() @ complex_matrix(s['mat1']) @ el.conj().T*1e9
    screen = inspect_native(directory/'result.json',family['ordering'])
    c = result['matched_comparison']
    return dict(
        Higgs_mass_GeV=result['masses_GeV']['25'],
        h125_mass_screen=result['h125_mass_screen'],
        cascade_status=result['cascade_status'],
        cascade_margin_GeV=result['cascade_margin_GeV'],
        oscillation_screen=screen,
        target_matrix_relative_residual=float(np.linalg.norm(native-mnu)/np.linalg.norm(mnu)),
        matched_calculator_relative_residual=c['matrix_relative_residual'],
        charged_lepton_matrix_relative_residual=c['charged_lepton_mass_matrix_relative_residual'],
        scalar_propagation=result['scalar_input_propagation']['status'],
        max_abs_running_yukawa=max(float(np.max(abs(complex_matrix(s[n]))))
                                    for n in ('fzee','Ye','Ye2')),
        running_v_GeV=math.hypot(s['vd'],s['vu']))


def verify_running_inputs(text, result):
    """Check the generated LOW-boundary contract before treating inputs as running values."""
    inputs = blocks(text)
    s = snapshot(result)
    for stem,name in (('FZEEIN','fzee'),('YEE2IN','Ye2')):
        matrix=np.array([[complex(inputs[stem][i+1,j+1],inputs['IM'+stem][i+1,j+1])
                          for j in range(3)] for i in range(3)])
        if not np.allclose(matrix,complex_matrix(s[name]),rtol=1e-12,atol=1e-18):
            raise ValueError(f'LOW input-to-running {name} propagation changed')


def fit_pass(o):
    return (o['h125_mass_screen']==o['cascade_status']==o['scalar_propagation']=='PASS'
            and abs(o['Higgs_mass_GeV']-HIGGS_TARGET)<=HIGGS_ACCURACY
            and o['oscillation_screen']['status']=='NATIVE_3SIGMA_SCREEN_PASS'
            and o['target_matrix_relative_residual']<MATRIX_ACCURACY
            and o['matched_calculator_relative_residual']<1e-10
            and o['charged_lepton_matrix_relative_residual']<1e-10
            and o['max_abs_running_yukawa']<1.)


def failure_causes(row, family, result, reference):
    """Freeze F,D and replace only kernel masses/mixing, then use actual F,D."""
    s = snapshot(result)
    cb,sb,el,er,fm,dm = running(s)
    f0,d0 = complex_matrix(family['F_mass_basis']),complex_matrix(family['D_mass_basis'])
    mnu,_,_ = target(family['ordering'],family['mlight_eV'],
                    family['alpha21_rad'],family['alpha31_rad'])
    sine = result['matched_comparison']['sin2gamma_ZP']
    kernel = zee_matrix(f0,d0,s['MFe'],s['MHm'][1:],sine)
    actual = zee_matrix(fm,dm,s['MFe'],s['MHm'][1:],sine)
    cr,sr,elr,err,_,_ = running(reference)
    sy_ref = cr*complex_matrix(reference['Ye'])+sr*complex_matrix(reference['Ye2'])
    sy_native = cb*complex_matrix(s['Ye'])+sb*complex_matrix(s['Ye2'])
    transported = er.conj() @ err.T @ d0 @ elr @ el.conj().T
    matching_shift = er.conj() @ ((sb/cb)*(sy_ref-sy_native)) @ el.conj().T
    trees = [t for t in result['scalar_trace'] if t['name']=='Mhh']
    matched_tree = [t for t in trees if t['stack'].startswith('Match_and_Run|OneLoopMasses|')][-1]
    return dict(
        scalar_id=row['scalar_id'], ordering=family['ordering'], card=family['card'],
        requested_tree_higgs_GeV=row['scalar']['mh1'],
        first_native_tree_higgs_GeV=math.sqrt(trees[0]['eigenvalues_GeV2'][0]),
        last_recorded_tree_higgs_GeV=math.sqrt(trees[-1]['eigenvalues_GeV2'][0]),
        final_matching_tree_higgs_GeV=math.sqrt(matched_tree['eigenvalues_GeV2'][0]),
        native_higgs_pole_GeV=result['masses_GeV']['25'],
        kernel_only_target_residual=float(np.linalg.norm(kernel-mnu)/np.linalg.norm(mnu)),
        fully_matched_target_residual=float(np.linalg.norm(actual-mnu)/np.linalg.norm(mnu)),
        D_proposal_diagonal=encode(np.diag(d0)), D_native_diagonal=encode(np.diag(dm)),
        D_shift_frobenius=float(np.linalg.norm(dm-d0)),
        predicted_matching_D_shift=encode(np.diag(matching_shift)),
        D_matching_identity_abs_residual=float(np.max(abs(dm-transported-matching_shift))),
        F_shift_frobenius=float(np.linalg.norm(fm-f0)))


def flag67_control(text):
    """Diagnostic only; this card can never be exported as a fitted benchmark."""
    current = None
    out = []
    count = 0
    for line in text.splitlines(keepends=True):
        header = re.match(r'\s*Block\s+(\w+)',line,re.I)
        if header:
            current = header[1].upper()
        if current=='SPHENOINPUT' and re.match(r'\s*67\s',line):
            line = ' 67 0 # diagnostic full one-loop Higgs, not a fit configuration\n'
            count += 1
        out.append(line)
    if count!=1:
        raise ValueError('Expected exactly one flag 67 for this controlled comparison')
    return ''.join(out)


def run(spheno, pilot, reference, output):
    spheno,pilot,reference,output = (p.resolve() for p in (spheno,pilot,reference,output))
    if sha(spheno/'bin/SPhenoZeeTypeII')!=BINARY_SHA:
        raise ValueError('Unrecognized original SPheno executable')
    catalog = json.loads((pilot/'proposals.json').read_text())
    if (catalog.get('schema_version')!=2 or
        catalog.get('status')!='SCALAR_NEUTRINO_PROPOSALS_NATIVE_UNKNOWN' or
        catalog.get('snapshot_sha256')!=sha(reference)):
        raise ValueError('Pilot schema or original running-reference provenance mismatch')
    ref_json = json.loads(reference.read_text())
    ref = (snapshot(ref_json) if 'snapshots' in ref_json else ref_json['snapshot'])
    rows = [r for r in catalog['outputs'] if r['families']]
    if len(rows)!=4:
        raise ValueError('This experiment requires the frozen four-configuration pilot')
    output.mkdir(parents=True,exist_ok=False)
    exe = build(spheno,output/'build')
    ledger = dict(schema_version=1, status='IN_PROGRESS',
                  selection='Smallest original max_abs_quartic, then scalar_id; no LFV selection',
                  original_binary_sha256=BINARY_SHA, diagnostic_binary_sha256=sha(exe),
                  pilot_sha256=sha(pilot/'proposals.json'),
                  reference_snapshot_sha256=sha(reference),
                  python_packages=dict(python=platform.python_version(),
                      numpy=np.__version__,scipy=scipy.__version__),
                  sources_sha256={n:sha(spheno/'ZeeTypeII'/n) for n in
                      ('SPhenoZeeTypeII.f90','TreeLevelMasses_ZeeTypeII.f90',
                       'LoopMasses_ZeeTypeII.f90','Boundaries_ZeeTypeII.f90',
                       'Model_Data_ZeeTypeII.f90','InputOutput_ZeeTypeII.f90',
                       'TadpoleEquations_ZeeTypeII.f90')},
                  libraries_sha256={n:sha(spheno/'lib'/n) for n in
                      ('libSPhenoZeeTypeII.a','libSPheno.a')},
                  policy=dict(higgs_target_GeV=HIGGS_TARGET,tree_bracket_GeV=TREE_BRACKET,
                      higgs_numerical_accuracy_GeV=HIGGS_ACCURACY,
                      target_matrix_accuracy=MATRIX_ACCURACY,max_refits=MAX_REFITS,
                      nufit_source=SOURCE,nufit_targets=NUFIT,
                      theory_requirements=THEORY, max_abs_yukawa=1.),
                  original_scalar_proposals=[{k:r[k] for k in
                      ('scalar_id','scalar','lagrangian','gates')} for r in rows],
                  trials=[],failure_causes=[],benchmarks=[],
                  outstanding_gates={name:'UNKNOWN' for name in UNKNOWN},
                  globally_allowed=None)

    def save():
        counts = {}
        for t in ledger['trials']:
            bucket = counts.setdefault(t['phase'],{})
            bucket[t['status']] = bucket.get(t['status'],0)+1
        ledger['statistics'] = counts
        (output/'ledger.json').write_text(json.dumps(ledger,indent=2,allow_nan=False)+'\n')

    def trial(phase,name,text,family=None,parameters=None,proposal=None,binary=None):
        d = output/name
        entry = dict(phase=phase,name=name,status='STARTED',
                     proposed_input_sha256=hashlib.sha256(text.encode()).hexdigest(),
                     executed_binary_sha256=sha(binary or exe))
        if proposal is not None:
            entry['inverse_update'] = proposal
        ledger['trials'].append(entry)
        if parameters is not None:
            entry['potential']=parameters
            entry['theory']=cheap_scalar(parameters)
            if not all(entry['theory'][key]=='PASS' for key in THEORY):
                entry['status']='REJECTED_THEORY'
                save()
                return None,d
        save()
        try:
            result = run_case(binary or exe,text,d)
            entry.update(input_sha256=result['input_sha256'],
                         spectrum_status=result['spectrum_status'])
            if (result['spectrum_status']!='PASS' or
                (binary is None and result['neutrino_status']!='PASS')):
                entry['status']='FAILED_NATIVE'
            else:
                entry['status']='NATIVE_COMPLETE'
                if family:
                    entry['observables']=observables(result,d,family)
                    entry['status']=('TARGET_PASS' if fit_pass(entry['observables'])
                                     else 'TARGET_FAIL')
                entry['masses_GeV']=result['masses_GeV']
                entry['cascade_status']=result['cascade_status']
                entry['cascade_margin_GeV']=result['cascade_margin_GeV']
                entry['scalar_propagation']=result['scalar_input_propagation']['status']
        except (subprocess.SubprocessError,ValueError,ArithmeticError) as exc:
            entry.update(status='FAILED_NATIVE_OR_ANALYSIS',error=str(exc))
            result = None
        entry['artifacts_sha256']={n:sha(d/n) for n in
            ('input.in','output.spc','result.json','run.log','snapshots.dat','scalar_trace.dat')
            if (d/n).is_file()}
        save()
        print(name,entry['status'],entry.get('masses_GeV',{}).get('25'),flush=True)
        return result,d

    baseline = {}
    for row in rows:
        for family in row['families']:
            card = pilot/'cards'/family['card']
            if sha(card)!=family['card_sha256']:
                raise ValueError('Frozen pilot card changed')
            text = fix_card(card.read_text(),scalar_updates(
                row['lagrangian'],math.tan(row['scalar']['beta'])))
            result,d = trial('baseline',Path(family['card']).stem+'_baseline',text,family,
                             row['lagrangian'])
            baseline[family['card']] = (text,result,d)
            if result and result['spectrum_status']=='PASS':
                ledger['failure_causes'].append(failure_causes(row,family,result,ref))
        family = row['families'][0]
        text,result,_ = baseline[family['card']]
        if result and result['spectrum_status']=='PASS':
            control,_ = trial('higgs_control',row['scalar_id']+'_full_one_loop',
                              flag67_control(text))
            if control and control['spectrum_status']=='PASS':
                ledger['failure_causes'][-len(row['families'])][
                    'native_higgs_without_effective_SM_GeV']=control['masses_GeV']['25']

    row = min(rows,key=lambda r:(r['gates']['max_abs_quartic'],r['scalar_id']))
    ledger['selected_scalar_id']=row['scalar_id']
    family = row['families'][0]
    source,initial,initial_dir = baseline[family['card']]
    evaluated = {row['scalar']['mh1']**2:(source,initial,initial_dir)}

    def objective(x):
        if x not in evaluated:
            par = potential(row,math.sqrt(x))
            text = fix_card(source,scalar_updates(par,math.tan(row['scalar']['beta'])))
            result,d = trial('higgs_calibration',f"higgs_{len(evaluated):02d}",text,
                             family,par)
            evaluated[x]=(text,result,d)
        _,result,_ = evaluated[x]
        if not result or result['spectrum_status']!='PASS' or result['cascade_status']!='PASS':
            raise ValueError('Higgs calibration blocked by theory, native spectrum or cascade')
        return result['masses_GeV']['25']**2-HIGGS_TARGET**2

    try:
        root = brentq(objective,*(mass**2 for mass in TREE_BRACKET),xtol=.1,maxiter=12)
        objective(root)
        _,calibrated,calibrated_dir = evaluated[root]
        if abs(calibrated['masses_GeV']['25']-HIGGS_TARGET)>HIGGS_ACCURACY:
            raise ValueError('Higgs calibration missed its declared numerical accuracy')
    except (ValueError,RuntimeError) as exc:
        ledger.update(status='BLOCKED_SCALAR_CALIBRATION',blocker=str(exc))
        save()
        return ledger
    par = potential(row,math.sqrt(root))
    ledger['calibration']=dict(tree_mass_GeV=math.sqrt(root),potential=par,
        unchanged_coordinates={k:v for k,v in row['scalar'].items() if k!='mh1'},
        calibrated_Higgs_mass_GeV=calibrated['masses_GeV']['25'],
        method='Bracketed root in tree mh1 squared, with all other scalar coordinates fixed')

    for family in row['families'][:2]:
        ordering = family['ordering']
        text = fix_card(baseline[family['card']][0],scalar_updates(
            par,math.tan(row['scalar']['beta'])))
        if family is row['families'][0]:
            result,d = calibrated,calibrated_dir
        else:
            result,d = trial('neutrino_start',f'{ordering}_start',text,family,par)
        converged = False
        for iteration in range(MAX_REFITS):
            if not result or result['neutrino_status']!='PASS':
                break
            try:
                verify_running_inputs(text,result)
            except ValueError as exc:
                ledger['trials'].append(dict(phase='neutrino_refit',name=ordering,
                    status='BLOCKED_YUKAWA_PROPAGATION',error=str(exc)))
                save()
                break
            update,proposal = refit_update(result,family)
            if proposal['max_abs_yukawa']>=1:
                ledger['trials'].append(dict(phase='neutrino_refit',
                    name=f'{ordering}_{iteration:02d}',status='REJECTED_YUKAWA_MAGNITUDE',
                    inverse_update=proposal))
                save()
                break
            text=fix_card(text,update)
            result,d=trial('neutrino_refit',f'{row["scalar_id"]}_{ordering}_refit_{iteration:02d}',
                           text,family,par,proposal)
            if result and result['neutrino_status']=='PASS':
                o=observables(result,d,family)
                if fit_pass(o):
                    converged=True
                    break
        if not converged:
            ledger['benchmarks'].append(dict(ordering=ordering,status='FAIL_REFIT',
                                             last_trial=str(d.relative_to(output))))
            save()
            continue
        # Independently rerun the exact final card on the uninstrumented binary.
        production_result,production=trial('production_control',f'{ordering}_production',
            text,parameters=par,binary=spheno/'bin/SPhenoZeeTypeII')
        mass_delta=(max(abs(production_result['masses_GeV'][k]-result['masses_GeV'][k])
                       for k in result['masses_GeV'])
                    if production_result and production_result['spectrum_status']=='PASS' else None)
        identical=(production/'output.spc').is_file() and sha(production/'output.spc')==sha(d/'output.spc')
        production_ok=mass_delta==0. and identical
        if not production_ok:
            ledger['benchmarks'].append(dict(ordering=ordering,
                status='FAIL_PRODUCTION_CONTROL',last_trial=production.name))
            save()
            continue
        destination=output/f'{row["scalar_id"]}_{ordering}_benchmark'
        destination.mkdir()
        shutil.copyfile(d/'input.in',destination/'input.in')
        shutil.copyfile(production/'output.spc',destination/'spectrum.slha')
        native={k:v for k,v in result.items() if k not in ('scalar_trace','snapshots','comparisons')}
        native.update(snapshots=[snapshot(result)],spectrum_snapshot_index=0,
            original_snapshot_index=result['spectrum_snapshot_index'],
            production_binary_sha256=BINARY_SHA,diagnostic_binary_sha256=sha(exe),
            sources_sha256=ledger['sources_sha256'],libraries_sha256=ledger['libraries_sha256'])
        (destination/'native.json').write_text(json.dumps(native,indent=2,allow_nan=False)+'\n')
        from Scan.native_gates import evaluate_entry
        partial=evaluate_entry(dict(point_id=f'{row["scalar_id"]}_{ordering}',
            ordering=ordering,card=str(destination/'input.in'),
            native=str(destination/'native.json'),spectrum=str(destination/'spectrum.slha')),output)
        (destination/'partial_gates.json').write_text(json.dumps(partial,indent=2,allow_nan=False)+'\n')
        ledger['benchmarks'].append(dict(
            scalar_id=row['scalar_id'],ordering=ordering,
            status='PASS_NECESSARY_SCREENS_ONLY' if production_ok else 'FAIL_PRODUCTION_CONTROL',
            production_control='PASS' if production_ok else 'FAIL',
            production_max_mass_difference_GeV=mass_delta,
            production_spectrum_identical=identical,
            selected_family={k:v for k,v in family.items() if k not in
                ('F_mass_basis','D_mass_basis','Ye2_gauge_proposal')},
            observables=o,potential=par,theory=cheap_scalar(par),
            partial_pipeline_verdict=partial['verdict'],
            directory=destination.name,
            artifacts_sha256={n:sha(destination/n) for n in
                ('input.in','spectrum.slha','native.json','partial_gates.json')},
            globally_allowed=None))
        save()
    ledger['status']=('BENCHMARK_NECESSARY_SCREENS_PASS' if any(
        b['status']=='PASS_NECESSARY_SCREENS_ONLY' for b in ledger['benchmarks']) else
        'BLOCKED_NATIVE_REFIT')
    unchanged=(sha(spheno/'bin/SPhenoZeeTypeII')==BINARY_SHA and
        all(sha(spheno/'ZeeTypeII'/n)==h for n,h in ledger['sources_sha256'].items()) and
        all(sha(spheno/'lib'/n)==h for n,h in ledger['libraries_sha256'].items()))
    ledger['original_installation_unchanged']='PASS' if unchanged else 'FAIL'
    if not unchanged:
        ledger['status']='BLOCKED_INSTALLATION_CHANGED'
    save()
    return ledger


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('spheno','pilot','output-dir'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--reference-snapshot',type=Path,default=Path(__file__).resolve().parents[1]/
                   'Model/reference/native_neutrino_snapshot.json')
    a=p.parse_args()
    report=run(a.spheno,a.pilot,a.reference_snapshot,a.output_dir)
    print(json.dumps({k:report[k] for k in ('status','statistics','benchmarks')},indent=2))
    return 0 if report['status']=='BENCHMARK_NECESSARY_SCREENS_PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
