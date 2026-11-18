"""One SPheno spectrum: check lepton Higgs-basis matching, not the PMNS fit."""
import argparse
import json
import math
from pathlib import Path

from Model.check_native import blocks


def matrix(data, name):
    expected = {(i, j) for i in range(1, 4) for j in range(1, 4)}
    real = data.get(name, {})
    imag = data.get('IM' + name, {})
    if set(real) != expected or (imag and set(imag) != expected):
        raise ValueError('Missing or incomplete ' + name + '/IM' + name)
    return [[complex(real[i, j], imag.get((i, j), 0.))
             for j in range(1, 4)] for i in range(1, 4)]


def inspect(text):
    data = blocks(text)
    tb = data.get('MINPAR', {}).get((10,))
    if tb is None or not math.isfinite(tb) or tb <= 0:
        raise ValueError('Missing/invalid native tanbeta')
    beta = math.atan(tb)
    c, s = math.cos(beta), math.sin(beta)
    ye, ye2, f = (matrix(data, name) for name in ('YE', 'YEE2', 'FZEE'))
    mass_yukawa = [[c*ye[i][j] + s*ye2[i][j] for j in range(3)] for i in range(3)]
    d = [[-s*ye[i][j] + c*ye2[i][j] for j in range(3)] for i in range(3)]
    scale = max(abs(mass_yukawa[i][i]) for i in range(3))
    if scale <= 0:
        raise ValueError('Zero charged-lepton mass Yukawa scale')
    offdiag = max(abs(mass_yukawa[i][j]) for i in range(3) for j in range(3) if i != j)
    imag_diag = max(abs(mass_yukawa[i][i].imag) for i in range(3))
    antisym = max(abs(f[i][j] + f[j][i]) for i in range(3) for j in range(3))
    diagonal = [mass_yukawa[i][i].real for i in range(3)]
    passed = (offdiag/scale <= 5e-7 and imag_diag/scale <= 5e-7
              and 0 < diagonal[0] < diagonal[1] < diagonal[2]
              and antisym <= 5e-9)
    return {'status': 'TREE_LEPTON_BASIS_PASS' if passed else 'TREE_LEPTON_BASIS_FAIL',
            'beta_rad': beta, 'mass_yukawa_real_diagonal': diagonal,
            'mass_yukawa_offdiag_relative': offdiag/scale,
            'mass_yukawa_diagonal_imag_relative': imag_diag/scale,
            'f_antisymmetric_max_abs': antisym,
            'higgs_basis_D': {'real': [[z.real for z in row] for row in d],
                              'imag': [[z.imag for z in row] for row in d]},
            'pmns_from_tree_UV': 'NOT_ESTABLISHED',
            'radiative_neutrino_mass': 'NOT_EVALUATED'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectrum', required=True, type=Path)
    a = p.parse_args()
    result = inspect(a.spectrum.read_text())
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'TREE_LEPTON_BASIS_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
