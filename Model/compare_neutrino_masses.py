"""Compare one-loop Zee analytic singular values against SPheno pole masses.

This is a one-card diagnostic, not a PMNS or full Majorana self-energy proof.
NumPy is used for complex 3x3 singular values.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from Model.audit_neutrino_native import parse_entries, sha
from Model.check_basis import matrix
from Model.check_native import blocks


def finite_kernel(m1, m2, ml):
    """Exact finite integral int_0^1 log[(x ml^2+(1-x)m1^2)/(x ml^2+(1-x)m2^2)] dx."""
    if m1 == m2:
        return 0.
    if min(m1, m2) <= 0 or ml < 0:
        raise ValueError("Invalid loop masses")
    if ml == 0:
        return math.log((m1/m2)**2)
    z, a, b = ml**2, m1**2, m2**2
    if z == a or z == b:
        raise ValueError("Degenerate internal mass requires limiting kernel")
    return math.log(a/b) + z/(a-z)*math.log(a/z) - z/(b-z)*math.log(b/z)


def compare(input_card, baseline, loop_output):
    native = dict(parse_entries(loop_output.read_text(), "MASS"))
    tree = dict(parse_entries(baseline.read_text(), "MASS"))
    params = dict(parse_entries(input_card.read_text(), "MINPAR"))
    for pdg in (11, 13, 15, 37, 9000037):
        if pdg not in tree or tree[pdg] <= 0:
            raise ValueError(f"Missing tree-level mass PDG {pdg}")
    for pdg in (12, 14, 16):
        if pdg not in native or native[pdg] < 0:
            raise ValueError(f"Missing/negative native neutrino mass PDG {pdg}")
    if 10 not in params or 13 not in params:
        raise ValueError("Missing input tanbeta/Mu4 in MINPAR (10/13)")

    sb = math.sin(math.atan(params[10]))
    cb = math.cos(math.atan(params[10]))
    data = blocks(baseline.read_text())
    f = np.array(matrix(data, "FZEE"), dtype=complex)
    ye = np.array(matrix(data, "YE"), dtype=complex)
    ye2 = np.array(matrix(data, "YEE2"), dtype=complex)
    mass = np.array([tree[k] for k in (11, 13, 15)], dtype=float)
    sy = cb*ye + sb*ye2
    dy = -sb*ye + cb*ye2
    scale = np.max(abs(sy))
    if (np.max(abs(sy-np.diag(np.diag(sy)))) > 1e-7*scale
            or np.max(abs(sy.diagonal().imag)) > 1e-7*scale
            or np.max(abs(f+f.T)) > 5e-9):
        raise ValueError("Unvalidated charged-lepton or Fzee basis")
    v_each = math.sqrt(2)*mass/sy.diagonal().real
    if not np.all(np.isfinite(v_each)) or np.any(v_each <= 0):
        raise ValueError("Cannot reconstruct the common electroweak vev")
    v = float(np.median(v_each))
    if np.max(abs(v_each/v-1)) > 3e-4:
        raise ValueError("Charged-lepton masses and Yukawa vevs are inconsistent")
    charged = sorted([tree[37], tree[9000037]])
    sine2gamma = -math.sqrt(2)*v*params[13]/(charged[0]**2-charged[1]**2)
    if abs(sine2gamma) > 1.001:
        raise ValueError("Mu4/charged spectrum gives unphysical sin(2 gamma)")
    loop_w = (1e9 * sine2gamma/(32*math.pi**2)
              *mass*np.array([finite_kernel(*charged, ml) for ml in mass]))
    a = (f*loop_w[np.newaxis,:])@dy
    mnu = a + a.T
    analytic = np.sort(np.linalg.svd(mnu, compute_uv=False))
    pole = np.sort(np.array([native[k]*1e9 for k in (12, 14, 16)]))
    worst = float(np.max(abs(analytic-pole))/max(float(np.max(analytic)), 1e-12))
    # The SPheno loop-mass switch also alters other pole masses; a tree-vs-pole
    # comparison may not be at an identical mass/parameter scheme.
    scalar_drift = {str(pdg): (native[pdg]-tree[pdg])/tree[pdg]
                    for pdg in (37, 9000037) if pdg in native}
    lepton_drift = {str(pdg): (native[pdg]-tree[pdg])/tree[pdg]
                    for pdg in (11, 13, 15) if pdg in native}
    result = {
        "status": "SINGULAR_VALUES_AGREE_ONLY" if worst <= 0.01 else "NATIVE_ANALYTIC_MASS_MISMATCH",
        "one_loop_majorana_matrix_validated": False,
        "input_sha256": sha(input_card), "tree_spectrum_sha256": sha(baseline),
        "loop_spectrum_sha256": sha(loop_output),
        "charged_masses_tree_GeV": charged, "lepton_masses_tree_GeV": mass.tolist(),
        "charged_masses_loop_GeV": {str(pdg): native[pdg]
                                      for pdg in (37, 9000037) if pdg in native},
        "charged_mass_relative_shifts": scalar_drift,
        "lepton_mass_relative_shifts": lepton_drift,
        "matched_native_running_parameter_scheme": False,
        "vev_GeV": v, "sin_2gamma_from_Mu4": sine2gamma,
        "analytic_mass_singular_values_eV": analytic.tolist(),
        "spheno_pole_masses_eV": pole.tolist(),
        "analytic_dmsq21_eV2": float(analytic[1]**2-analytic[0]**2),
        "native_dmsq21_eV2": float(pole[1]**2-pole[0]**2),
        "analytic_dmsq31_eV2": float(analytic[2]**2-analytic[0]**2),
        "native_dmsq31_eV2": float(pole[2]**2-pole[0]**2),
        "max_abs_mass_difference_over_analytic_max": worst,
        "analytic_Mnu_eV": {
            "real": mnu.real.tolist(), "imag": mnu.imag.tolist()},
        "warning": "Compare only neutrino singular values. Native complex self-energy and UV basis remain unverified."
    }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--baseline", required=True, type=Path)
    p.add_argument("--loop", required=True, type=Path)
    a = p.parse_args()
    print(json.dumps(compare(a.input, a.baseline, a.loop), indent=2))


if __name__ == "__main__":
    main()
