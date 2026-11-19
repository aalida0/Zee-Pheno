"""Build an isolated, instrumented SPheno copy; compare actual one-loop inputs.

Never writes to the supplied installation. Requires its compiled libraries and
modules. Each snapshot is a separate OneLoopFv invocation (later Higgs coupling
calculations may call it again with different parameters).
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

from Model.audit_neutrino_native import parse_entries, sha
from Neutrino.fit import zee_matrix
BINARY_SHA = "2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52"

# Ordered Fortran output schema: shape, complex.
PARAMETERS = {'vd': ((), False), 'vu': ((), False), 'Mu4': ((), False),
              'ZP': ((3, 3), False),
              **{n: ((3, 3), True) for n in ('fzee', 'Ye', 'Ye2', 'ZEL', 'ZER', 'UV')}}
SELFENERGY = {**{n: ((3,), False) for n in ('MHm', 'MFe', 'MFv')},
              **{n: ((3, 3), True) for n in ('sigSL', 'sigSR', 'mat1')},
              **{n: ((3, 3, 3), True) for n in
                 ('cplUFvFecHmL', 'cplUFvFecHmR', 'cplcFeUFvHmL', 'cplcFeUFvHmR')}}
RETURNED = {'MFv_1L': ((3,), False), 'UV_1L': ((3, 3), True)}
SCHEMA = {**PARAMETERS, **SELFENERGY, **RETURNED}


def dump_code(schema):
    result = ''
    for name, (_, complex_) in schema.items():
        result += f"write(91,*) '{name}'\n"
        result += f"write(91,'(*(ES26.17E3,1X))') Real({name},dp)\n"
        if complex_:
            result += f"write(91,'(*(ES26.17E3,1X))') Aimag({name})\n"
    return result


def instrument(source, rotate_tree=False):
    before = "open(91,file='snapshots.dat',status='unknown',position='append')\n"
    close = 'close(91)\n'
    marker = 'Call OneLoopFv(MHm'
    if source.count(marker) != 1:
        raise ValueError('Unexpected generated OneLoopFv call site')
    source = source.replace(marker, before + "write(91,*) 'BEGIN'\n"
                            + dump_code(PARAMETERS) + close + marker, 1)
    if rotate_tree:
        source = source.replace('Call CouplingsForLoopMasses(',
                                'Do i1=1,3\nDo i2=1,3\n'
                                'UV(i1,i2)=Exp(Cmplx(0._dp,2*pi*(i1-1)*(i2-1)/3,dp))/Sqrt(3._dp)\n'
                                'End Do\nEnd Do\nCall CouplingsForLoopMasses(', 1)
    start = source.index('Subroutine OneLoopFv(')
    end = source.index('End Subroutine OneLoopFv', start)
    routine = source[start:end]
    marker = 'If (ForceRealMatrices) mat1'
    if routine.count(marker) != 3:
        raise ValueError('Unexpected generated mass assembly')
    routine = routine.replace(marker, before + dump_code(SELFENERGY) + close + marker, 1)
    # Retain every momentum evaluation to test the pole iteration, too.
    routine = routine.replace(marker, "open(92,file='momenta.dat',status='unknown',position='append')\n"
                              "write(92,'(*(ES26.17E3,1X))') p2, Real(mat1,dp), Aimag(mat1)\n"
                              'close(92)\n' + marker)
    routine = routine.replace('Iname = Iname -1', before + dump_code(RETURNED)
                              + "write(91,*) 'END'\n" + close + 'Iname = Iname -1')
    return source[:start] + routine + source[end:]


def read_snapshots(path):
    lines = iter(path.read_text().splitlines())
    snapshots = []
    current = None
    for line in lines:
        name = line.strip()
        if name == 'BEGIN':
            if current is not None:
                raise ValueError('Incomplete native snapshot')
            current = {}
        elif name == 'END':
            if current is None or set(current) != set(SCHEMA):
                raise ValueError('Incomplete native snapshot')
            snapshots.append(current)
            current = None
        else:
            shape, complex_ = SCHEMA[name]
            a = np.fromstring(next(lines), sep=' ')
            if complex_:
                a = a + 1j*np.fromstring(next(lines), sep=' ')
            current[name] = a.reshape(shape, order='F')
    if current is not None or not snapshots:
        raise ValueError('Missing complete native snapshot')
    return snapshots


def encode(snapshot):
    return {k: {'real': a.real.tolist(), 'imag': a.imag.tolist()}
            if np.iscomplexobj(a) else a.tolist() for k, a in snapshot.items()}


def analyze(s):
    v = float(np.hypot(s['vd'], s['vu']))
    cb, sb = s['vd']/v, s['vu']/v
    d = -sb*s['Ye'] + cb*s['Ye2']
    z, el, er = s['ZP'], s['ZEL'], s['ZER']
    sine = 2*z[1, 2]*(-sb*z[1, 0] + cb*z[1, 1])
    # e_mass = ZEL e_gauge, ec_mass = ZER ec_gauge. The external
    # neutrinos in UFv vertices are unrotated gauge-basis fields.
    f_mass = el.conj() @ s['fzee'] @ el.conj().T
    d_mass = er.conj() @ d @ el.conj().T
    analytic_mass = zee_matrix(f_mass, d_mass, s['MFe'], s['MHm'][1:], sine)
    native = s['mat1']*1e9
    native_mass = el.conj() @ native @ el.conj().T
    scale = max(float(np.max(abs(native))), 1e-30)
    masses = np.sort(np.linalg.svd(analytic_mass, compute_uv=False))
    u = s['UV_1L']
    diag = u.conj() @ native @ u.conj().T
    # Independent sum of both actual generated charged-scalar vertex arrays.
    # The arbitrary B0 scale cancels by orthogonality of charged mixing rows.
    ml, mh = s['MFe'], s['MHm']
    b0 = np.array([[1-np.log(h*h/160**2)
                    - m*m/(h*h-m*m)*np.log(h*h/(m*m)) for h in mh] for m in ml])
    sr = (np.einsum('aik,bik,i,ik->ab', s['cplUFvFecHmL'],
                    s['cplUFvFecHmR'].conj(), ml, b0)
          + np.einsum('iak,ibk,i,ik->ab', s['cplcFeUFvHmL'],
                      s['cplcFeUFvHmR'].conj(), ml, b0))
    vertex = -sr*1e9/(16*np.pi**2)
    m_yukawa = er.conj() @ (cb*s['Ye']+sb*s['Ye2']) @ el.conj().T
    return {
        'matched_analytic_masses_eV': masses.tolist(),
        'native_returned_masses_eV': (s['MFv_1L']*1e9).tolist(),
        'dmsq21_eV2': float(masses[1]**2-masses[0]**2),
        'dmsq31_eV2': float(masses[2]**2-masses[0]**2),
        'matrix_relative_residual': float(np.max(abs(native_mass-analytic_mass))/scale),
        'vertex_relative_residual': float(np.max(abs(vertex-native))/scale),
        'majorana_symmetry_relative_residual': float(np.max(abs(native-native.T))/scale),
        'sigSL_conjugate_sigSR_relative_residual':
            float(np.max(abs(s['sigSL'].conj()-s['sigSR']))*1e9/scale),
        'takagi_relative_residual': float(np.max(abs(diag-np.diag(s['MFv_1L']*1e9)))/scale),
        'charged_lepton_mass_matrix_relative_residual':
            float(np.max(abs(m_yukawa*v/np.sqrt(2)-np.diag(ml)))/max(ml)),
        'sin2gamma_ZP': float(sine),
        'sin2gamma_Mu4': float(-np.sqrt(2)*v*s['Mu4']/(mh[1]**2-mh[2]**2)),
        'D_diagonal': {'real': d.diagonal().real.tolist(), 'imag': d.diagonal().imag.tolist()},
        'native_Mnu_eV': encode({'matrix': native})['matrix'],
    }


def require_native_snapshots(directory):
    """SPheno sometimes exits 0 after an early spectrum/tadpole failure."""
    path=directory/'snapshots.dat'
    if path.is_file() and path.stat().st_size:
        return path
    log=directory/'run.log'
    tail="\n".join(log.read_text(errors='replace').splitlines()[-30:]) if log.exists() else "(run.log missing)"
    raise RuntimeError(
        "SPheno did not reach the instrumented OneLoopFv neutrino calculation; "
        "this is NOT a native neutrino-screen result. "
        f"Inspect {log} and {directory/'output.spc'}.\n"
        f"SPheno log tail:\n{tail}")


def run(spheno, card, directory, rotate_tree=False):
    spheno, card, directory = (p.resolve() for p in (spheno, card, directory))
    if sha(spheno/'bin/SPhenoZeeTypeII') != BINARY_SHA:
        raise ValueError('Unrecognized production executable')
    directory.mkdir(parents=True, exist_ok=False)
    model = spheno/'ZeeTypeII'
    sources = ('LoopMasses_ZeeTypeII.f90', 'SPhenoZeeTypeII.f90')
    for name in sources:
        content = (model/name).read_text()
        (directory/name).write_text(instrument(content, rotate_tree)
                                   if name.startswith('LoopMasses') else content)
    def command(args, log):
        result = subprocess.run(args, cwd=directory, capture_output=True, text=True, timeout=300)
        (directory/log).write_text(result.stdout+result.stderr)
        result.check_returncode()
    for name in sources:
        command(['gfortran', '-c', '-g', '-ffree-line-length-none', '-I'+str(directory),
                 '-I'+str(spheno/'include'), '-J'+str(directory), name,
                 '-o', name.replace('.f90', '.o')], name+'.build.log')
    command(['gfortran', '-o', 'SPhenoDiagnostic', 'SPhenoZeeTypeII.o',
             'LoopMasses_ZeeTypeII.o', str(spheno/'lib/libSPhenoZeeTypeII.a'),
             str(spheno/'lib/libSPheno.a')], 'link.log')
    (directory/'input.in').write_bytes(card.read_bytes())
    command([str(directory/'SPhenoDiagnostic'), str(directory/'input.in'),
             str(directory/'output.spc')], 'run.log')
    snapshots = read_snapshots(require_native_snapshots(directory))
    masses = dict(parse_entries((directory/'output.spc').read_text(), 'MASS'))
    pole = np.sort([masses[k]*1e9 for k in (12, 14, 16)])
    selected = [i for i, s in enumerate(snapshots)
                if np.allclose(s['MFv_1L']*1e9, pole, rtol=1e-7, atol=1e-15)]
    if len(selected) != 1:
        raise ValueError('Cannot uniquely associate native snapshot with SLHA masses')
    result = {'input_sha256': sha(card), 'rotated_tree_neutrinos': rotate_tree,
              'production_binary_sha256': BINARY_SHA,
              'diagnostic_binary_sha256': sha(directory/'SPhenoDiagnostic'),
              'sources_sha256': {n: sha(model/n) for n in sources},
              'libraries_sha256': {n: sha(spheno/'lib'/n)
                                   for n in ('libSPhenoZeeTypeII.a', 'libSPheno.a')},
              'spectrum_snapshot_index': selected[0], 'spheno_pole_masses_eV': pole.tolist(),
              'snapshots': [encode(s) for s in snapshots],
              'comparisons': [analyze(s) for s in snapshots]}
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spheno', required=True, type=Path)
    p.add_argument('--card', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--rotate-tree-neutrinos', action='store_true')
    a = p.parse_args()
    r = run(a.spheno, a.card, a.output_dir, a.rotate_tree_neutrinos)
    print(json.dumps(r['comparisons'][r['spectrum_snapshot_index']], indent=2))


if __name__ == '__main__':
    main()
