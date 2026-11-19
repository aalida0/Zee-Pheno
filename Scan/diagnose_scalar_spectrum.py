"""Isolated, output-only scalar/neutrino instrumentation of the frozen SPheno build.

No production sources, libraries, cards or binaries are changed. This is a
diagnostic experiment, not a physical acceptance gate.
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from Model.audit_neutrino_native import parse_entries, sha
from Model.check_native import fix_card
from Model.diagnose_neutrino_native import (
    BINARY_SHA, analyze, encode, instrument, read_snapshots)
from Neutrino.fit import REFERENCE_CARD_SHA

TREE_INPUTS = {
    'Mhh': 'M12S Lam5 M112 M222 Lam1 Lam3 Lam4 Lam2 vd vu',
    'MAh': 'M12S Lam5 M112 M222 Lam1 Lam3 Lam4 Lam2 vd vu',
    'MHm': 'M12S Mu4 Lam5 K12 M112 M222 MC2 Lam1 Lam3 Lam4 Lam2 K1 K2 vd vu',
}


def potential_hessian(p):
    """Differentiate the polynomial Lagrangian, independently of mass formulas.

    Symmetric differences with Richardson cancellation are exact for a quartic
    polynomial in exact arithmetic. The step is v/8, avoiding tiny-step loss.
    Coordinates are (phi1,phi2,sigma1,sigma2,Re H1+,Re H2+,Re S+),
    all canonically normalized. No gauge-fixing potential is included.
    """
    def potential(x):
        h1 = np.array([x[4], p['vd']+x[0]+1j*x[2]])/np.sqrt(2)
        h2 = np.array([x[5], p['vu']+x[1]+1j*x[3]])/np.sqrt(2)
        singlet = x[6]/np.sqrt(2)
        a,b = np.vdot(h1,h1).real,np.vdot(h2,h2).real
        z = np.vdot(h1,h2)
        c = singlet**2
        return (p['M112']*a+p['M222']*b-2*p['M12S']*z.real+p['MC2']*c
                +p['Lam1']/2*a*a+p['Lam2']/2*b*b+p['Lam3']*a*b
                +p['Lam4']*abs(z)**2+p['Lam5']*(z*z).real
                +p.get('LamC',0)*c*c+p['K1']*c*a+p['K2']*c*b
                -2*p['K12']*c*z.real
                +2*p['Mu4']*singlet*(h1[0]*h2[1]-h1[1]*h2[0]).real)
    origin = np.zeros(7)
    f0 = potential(origin)
    def differentiate(step):
        basis = np.eye(7)*step
        hessian = np.zeros((7,7))
        gradient = np.zeros(7)
        for i,a in enumerate(basis):
            gradient[i] = (potential(a)-potential(-a))/(2*step)
            hessian[i,i] = (potential(a)+potential(-a)-2*f0)/step**2
            for j,b in enumerate(basis[:i]):
                hessian[i,j] = hessian[j,i] = (
                    potential(a+b)-potential(a-b)-potential(-a+b)+potential(-a-b)
                )/(4*step**2)
        return gradient,hessian
    step = np.hypot(p['vd'],p['vu'])/8
    g,h = differentiate(step)
    g2,h2 = differentiate(step/2)
    return (4*g2-g)/3,(4*h2-h)/3


def tree_comparison(trace):
    first = {name:next(t for t in trace if t['name']==name) for name in TREE_INPUTS}
    p = first['MHm']['inputs']
    v = np.hypot(p['vd'],p['vu'])
    cb,sb = p['vd']/v,p['vu']/v
    gradient,hessian = potential_hessian(p)
    u = np.array([[-sb,cb,0],[0,0,1]])
    odd = np.array([-sb,cb])
    native_even = np.array(first['Mhh']['matrix_GeV2'])
    native_odd = odd @ np.array(first['MAh']['matrix_GeV2']) @ odd
    native_charged = u @ np.array(first['MHm']['matrix_GeV2']) @ u.T
    residual = max(np.max(abs(hessian[:2,:2]-native_even)),
                   abs(odd @ hessian[2:4,2:4] @ odd-native_odd),
                   np.max(abs(u @ hessian[4:,4:] @ u.T-native_charged)))
    return dict(status='PASS' if residual<1e-6 else 'FAIL',
                method='Direct potential Hessian; project out gauge-fixing Goldstones',
                max_matrix_abs_residual_GeV2=float(residual),
                max_tadpole_abs_GeV3=float(np.max(abs(gradient))),
                CP_even_eigenvalues_GeV2=np.linalg.eigvalsh(hessian[:2,:2]).tolist(),
                CP_odd_physical_GeV2=float(odd @ hessian[2:4,2:4] @ odd),
                charged_physical_matrix_GeV2=(u @ hessian[4:,4:] @ u.T).tolist(),
                charged_physical_eigenvalues_GeV2=np.linalg.eigvalsh(
                    u @ hessian[4:,4:] @ u.T).tolist())


def instrument_tree(source):
    for name, inputs in TREE_INPUTS.items():
        start = source.index('Subroutine Calculate' + name + '(')
        end = source.index('End Subroutine Calculate' + name, start)
        routine = source[start:end]
        marker = f'Call EigenSystem(mat,{name}2,'
        if routine.count(marker) != 1:
            raise ValueError('Unexpected tree diagonalization site')
        before = ("open(93,file='scalar_trace.dat',status='unknown',position='append')\n"
                  + f"write(93,*) '{name}'\n"
                  + "write(93,*) (Trim(NameOfUnit(i1))//'|',i1=1,Iname)\n"
                  + f"write(93,'(*(ES26.17E3,1X))') {','.join(inputs.split())}\n"
                  + "write(93,'(*(ES26.17E3,1X))') mat\n")
        routine = routine.replace(marker, before + marker)
        match = re.search(r'Call EigenSystem\(mat,' + name + r'2,[^\n]+\n', routine)
        routine = routine[:match.end()] + (
            f"write(93,'(*(ES26.17E3,1X))') {name}2\nclose(93)\n"
        ) + routine[match.end():]
        source = source[:start] + routine + source[end:]
    return source


def read_trace(path):
    if not path.exists():
        return []
    lines = iter(path.read_text().splitlines())
    result = []
    for line in lines:
        name = line.strip()
        stack = next(lines).strip()
        values = np.fromstring(next(lines), sep=' ')
        n = 3 if name == 'MHm' else 2
        matrix = np.fromstring(next(lines), sep=' ').reshape((n, n), order='F')
        eigenvalues = np.fromstring(next(lines), sep=' ')
        result.append(dict(name=name, stack=stack,
                           inputs=dict(zip(TREE_INPUTS[name].split(), values.tolist())),
                           matrix_GeV2=matrix.tolist(), eigenvalues_GeV2=eigenvalues.tolist(),
                           eigenvalue_residual_GeV2=float(np.max(abs(
                               np.linalg.eigvalsh(matrix)-np.sort(eigenvalues))))))
    return result


def command(args, directory, log):
    run = subprocess.run([str(a) for a in args], cwd=directory, capture_output=True,
                         text=True, timeout=300)
    (directory/log).write_text(run.stdout+run.stderr)
    run.check_returncode()
    return run.stdout+run.stderr


def scalar_propagation(text, spc):
    minimum = dict(parse_entries(text, 'MINPAR'))
    quartic = dict(parse_entries(spc, 'ZEEQUARTIC'))
    pot = dict(parse_entries(spc, 'ZEEPOT'))
    hmix = dict(parse_entries(spc,'HMIX'))
    residual = max(
        [abs(quartic[i]-minimum[i])/max(1,abs(minimum[i])) for i in range(1,10)]
        + [abs(pot[k]-minimum[i])/max(1,abs(minimum[i]))
           for k,i in ((3,11),(4,12),(5,13))])
    tanbeta_residual = abs(hmix[103]/hmix[102]-minimum[10])/max(1,abs(minimum[10]))
    return dict(status='PASS' if max(residual,tanbeta_residual)<5e-9 else 'FAIL',
                max_scaled_residual=residual,
                tanbeta_scaled_residual=tanbeta_residual,
                method='12 independent scalar inputs against 9-significant-digit SLHA output',
                MINPAR=minimum, ZEEQUARTIC=quartic, ZEEPOT=pot,
                HMIX=hmix)


def build(spheno, directory):
    directory.mkdir()
    sources = ('TreeLevelMasses_ZeeTypeII.f90', 'LoopMasses_ZeeTypeII.f90',
               'SPhenoZeeTypeII.f90')
    for name in sources:
        source = (spheno/'ZeeTypeII'/name).read_text()
        if name.startswith('Tree'):
            source = instrument_tree(source)
        elif name.startswith('Loop'):
            source = instrument(source)
        (directory/name).write_text(source)
        command(['gfortran', '-c', '-g', '-ffree-line-length-none',
                 '-I'+str(directory), '-I'+str(spheno/'include'),
                 '-J'+str(directory), name, '-o', name.replace('.f90', '.o')],
                directory, name+'.build.log')
    command(['gfortran', '-o', 'SPhenoScalarDiagnostic', 'SPhenoZeeTypeII.o',
             'LoopMasses_ZeeTypeII.o', 'TreeLevelMasses_ZeeTypeII.o',
             spheno/'lib/libSPhenoZeeTypeII.a', spheno/'lib/libSPheno.a'],
            directory, 'link.log')
    return directory/'SPhenoScalarDiagnostic'


def run_case(exe, text, directory):
    directory.mkdir()
    card = directory/'input.in'
    card.write_text(text)
    log = command([exe, card, directory/'output.spc'], directory, 'run.log')
    masses = (dict(parse_entries((directory/'output.spc').read_text(), 'MASS'))
              if (directory/'output.spc').exists() else {})
    trace = read_trace(directory/'scalar_trace.dat')
    result = dict(input_sha256=sha(card), scalar_trace=trace,
                  masses_GeV={str(k): v for k, v in masses.items()},
                  spectrum_status='PASS' if all(masses.get(k, 0)>0 for k in
                      (25,35,36,37,9000037)) and 'TerminateProgram' not in log else 'FAIL',
                  neutrino_status='UNKNOWN')
    if set(TREE_INPUTS) <= {t['name'] for t in trace}:
        result['independent_tree'] = tree_comparison(trace)
    if (directory/'snapshots.dat').exists():
        snapshots = read_snapshots(directory/'snapshots.dat')
        pole = np.sort([masses.get(k, 0)*1e9 for k in (12,14,16)])
        selected = [i for i,s in enumerate(snapshots) if
                    np.allclose(s['MFv_1L']*1e9, pole, rtol=1e-7, atol=1e-15)]
        result.update(snapshots=[encode(s) for s in snapshots],
                      comparisons=[analyze(s) for s in snapshots],
                      spheno_pole_masses_eV=pole.tolist(),
                      spectrum_snapshot_index=selected[0] if len(selected)==1 else None,
                      neutrino_status='PASS' if len(selected)==1 else 'UNKNOWN')
        if len(selected)==1:
            result['matched_comparison'] = result['comparisons'][selected[0]]
    if result['spectrum_status']=='PASS':
        spc = (directory/'output.spc').read_text()
        result['scalar_input_propagation'] = scalar_propagation(text,spc)
        result['h125_mass_screen'] = 'PASS' if abs(masses[25]-125.25)<=3 else 'FAIL'
        gap = masses[9000037]-masses[37]-masses[23]
        result.update(cascade_status='PASS' if gap>0 else 'FAIL',
                      cascade_margin_GeV=gap)
    ordering = re.search(r'_(NO|IO)_',directory.name)
    if ordering and result['neutrino_status']=='PASS':
        # The existing source-matched screen uses this persisted snapshot.
        from Neutrino.fit import inspect_native
        path = directory/'result.json'
        path.write_text(json.dumps(result,indent=2)+'\n')
        result['oscillation_screen'] = inspect_native(path,ordering[1])
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def matching(text, order):
    current = None
    out = []
    for line in text.splitlines(keepends=True):
        header = re.match(r'\s*Block\s+(\w+)', line, re.I)
        if header:
            current = header[1].upper()
        if current == 'SPHENOINPUT' and re.match(r'\s*19\s', line):
            line = f' 19 {order} # controlled matching diagnostic\n'
        out.append(line)
    return ''.join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('spheno', 'reference', 'pilot', 'output-dir'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    spheno, reference, pilot, output = (x.resolve() for x in
                                       (a.spheno,a.reference,a.pilot,a.output_dir))
    if sha(spheno/'bin/SPhenoZeeTypeII') != BINARY_SHA or sha(reference) != REFERENCE_CARD_SHA:
        raise ValueError('Frozen reference/binary hash mismatch')
    output.mkdir(parents=True, exist_ok=False)
    catalog = json.loads((pilot/'proposals.json').read_text())
    rows = [r for r in catalog['outputs'] if r['families']]
    ref = reference.read_text()
    cases = [('reference', ref)]
    first = rows[0]
    card = pilot/'cards'/first['families'][0]['card']
    if sha(card) != first['families'][0]['card_sha256']:
        raise ValueError('Pilot card hash mismatch')
    old = card.read_text()
    cases.append(('s00001_original', old))
    # Preserve the failed scalar input while replacing only the lepton inputs.
    from Model.check_native import MATRICES, blocks
    data = blocks(ref)
    cases.append(('s00001_reference_yukawa', fix_card(old, {
        (name,key):value for name in MATRICES for key,value in data[name].items()})))
    production_controls = {}
    for name,text in cases[:2]:
        production_controls[name] = run_case(
            spheno/'bin/SPhenoZeeTypeII',text,output/('production_'+name))
    exe = build(spheno, output/'build')
    from Scan.physical_proposals import scalar_updates
    for row in rows:
        if row is not first:
            cases.append((row['scalar_id']+'_original',
                          (pilot/'cards'/row['families'][0]['card']).read_text()))
        for family in row['families']:
            path = pilot/'cards'/family['card']
            if sha(path) != family['card_sha256']:
                raise ValueError('Pilot card hash mismatch')
            corrected = fix_card(path.read_text(),scalar_updates(
                row['lagrangian'], np.tan(row['scalar']['beta'])))
            name = path.stem+'_corrected'
            cases.append((name,corrected))
            if row is first and family is row['families'][0]:
                for order in (0,1,2,-1):
                    cases.append((name+f'_matching_{order}', matching(corrected,order)))
    report = dict(production_binary_sha256=BINARY_SHA,
                  production_controls=production_controls,
                  diagnostic_binary_sha256=sha(exe), reference_sha256=sha(reference),
                  pilot_sha256=sha(pilot/'proposals.json'),
                  sources_sha256={n:sha(spheno/'ZeeTypeII'/n) for n in
                      ('InputOutput_ZeeTypeII.f90','Boundaries_ZeeTypeII.f90',
                       'TreeLevelMasses_ZeeTypeII.f90','TadpoleEquations_ZeeTypeII.f90',
                       'LoopMasses_ZeeTypeII.f90','SPhenoZeeTypeII.f90')},
                  scalar_proposals=[{k:r[k] for k in ('scalar_id','scalar','lagrangian')}
                                    for r in rows],
                  libraries_sha256={n:sha(spheno/'lib'/n) for n in
                      ('libSPhenoZeeTypeII.a','libSPheno.a')}, cases={})
    for name,text in cases:
        result = run_case(exe, text, output/name)
        report['cases'][name] = {k:v for k,v in result.items()
                                  if k not in ('snapshots','scalar_trace','comparisons')}
        report['cases'][name]['first_tree'] = result['scalar_trace'][:3]
        report['cases'][name]['artifacts_sha256'] = {
            n:sha(output/name/n) for n in
            ('input.in','run.log','output.spc','result.json','scalar_trace.dat','snapshots.dat')
            if (output/name/n).exists()}
        (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        print(name, result['spectrum_status'], result['masses_GeV'].get('25'), flush=True)


if __name__ == '__main__':
    main()
