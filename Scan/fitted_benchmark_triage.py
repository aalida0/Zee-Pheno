"""Source-bound, fixed-benchmark constraint triage; no physical acceptance.

Reads the two archived s00017 native controls, re-evaluates the validated
scalar one-loop LFV subset, and audits SPheno observable and decay blocks.
Partial calculations, model-dependent cosmology, and unvalidated FlavorKit
outputs are explicitly not promoted to exclusion/acceptance gates.
"""
import argparse
import json
import math
import re
from pathlib import Path

from LFV.pilot import evaluate as lfv_one_loop
from Model.audit_neutrino_native import parse_entries, sha

ROOT = Path(__file__).resolve().parent / "reference/native_benchmark"
MEG_II_2025_LIMIT_90CL = 1.5e-13
MEG_II_SOURCE = "MEG II 2021-2022 data, arXiv:2504.15711, 90% CL"


def mixing_entries(text, name):
    """Read one real 3x3 SARAH mixing block without assuming determinant +1."""
    entries = {}
    inside = False
    for line in text.splitlines():
        match = re.match(r"\s*Block\s+(\w+)", line, re.I)
        if match:
            inside = match[1].upper() == name.upper()
            continue
        if inside and re.match(r"\s*DECAY\b", line, re.I):
            inside = False
        if inside:
            cols = line.split("#", 1)[0].split()
            if len(cols) == 3:
                try:
                    i, j = int(cols[0]), int(cols[1])
                    entries[i, j] = float(cols[2].replace("D", "E"))
                except ValueError:
                    pass
    if set(entries) != {(i, j) for i in range(1, 4) for j in range(1, 4)}:
        raise ValueError(f"Incomplete {name}")
    return entries


def decay_block(text, pdg):
    """Parse one SLHA DECAY channel list; do not infer absent channels."""
    active, total_width, rows = False, None, []
    for line in text.splitlines():
        match = re.match(r"\s*DECAY\s+(\d+)\s+(\S+)", line, re.I)
        if match:
            active = int(match[1]) == pdg
            if active:
                total_width = float(match[2].replace("D", "E"))
            continue
        if not active or line.lstrip().startswith("#"):
            continue
        cols = line.split("#", 1)[0].split()
        if not cols:
            continue
        try:
            br, nda = float(cols[0]), int(cols[1])
            daughters = tuple(int(x) for x in cols[2:])
        except (IndexError, ValueError) as exc:
            raise ValueError("Invalid native decay table") from exc
        if len(daughters) != nda or not 0 <= br <= 1:
            raise ValueError("Invalid native decay width/BR")
        rows.append((br, daughters))
    if total_width is None or not math.isfinite(total_width) or total_width <= 0:
        raise ValueError(f"Missing or invalid DECAY {pdg}")
    return total_width, rows


def evaluate(order, root=ROOT):
    if order not in ("NO", "IO"):
        raise ValueError("NO/IO order required")
    p = root / f"s00017_{order}_benchmark"
    card, native, spectrum = p / "input.in", p / "native.json", p / "spectrum.slha"
    pilot = lfv_one_loop(native, card, spectrum)  # verifies card/native SHA
    text = spectrum.read_text()
    masses = dict(parse_entries(text, "MASS"))
    gauge = dict(parse_entries(text, "GAUGE"))
    flavorkit_lfv = dict(parse_entries(text, "FlavorKitLFV"))
    flavorkit_qfv = dict(parse_entries(text, "FlavorKitQFV"))
    ew = dict(parse_entries(text, "SPhenoLowEnergy"))
    zp = mixing_entries(text, "CHARGEMIX")
    heavy = masses[9000037]
    light = masses[37]
    mz = masses[23]
    width, channels = decay_block(text, 9000037)
    matches = [br for br, daughters in channels if daughters == (37, 23)]
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous charged cascade channel")
    cascade = matches[0]
    if heavy <= light + mz:
        raise ValueError("Expected open s00017 cascade")
    # The physical doublet/singlet mixing rows give 2 sin(gamma) cos(gamma).
    # g1 in GAUGE is SM hypercharge here, not GUT-normalized g1.
    sin2gamma = 2 * abs(zp[2, 3] * zp[3, 3])
    g1, g2 = gauge[1], gauge[2]
    cosw = g2 / math.hypot(g1, g2)
    gz12 = g2 * sin2gamma / (4 * cosw)
    lambda12 = (heavy**2 - (light + mz)**2) * (heavy**2 - (light - mz)**2)
    analytic_width = gz12**2 * lambda12**1.5 / (
        16 * math.pi * heavy**3 * mz**2)
    native_width = width * cascade
    residual = abs(analytic_width - native_width) / native_width
    if residual > 2e-5:
        raise ValueError("Independent scalar-Z width does not match SPheno")
    subset = {
        ch: record["BR_scalar_one_loop_subset"]
        for ch, record in pilot["processes"].items()
    }
    if any(x < 0 or not math.isfinite(x) for x in subset.values()):
        raise ValueError("Invalid partial LFV output")
    msum = sum(masses[i] for i in (12, 14, 16)) * 1e9
    if not 0 < msum < 2:
        raise ValueError("Unexpected native neutrino mass sum")
    return {
        "ordering": order,
        "status": "FIXED_BENCHMARK_PARTIAL_CONSTRAINT_EVIDENCE",
        "input_provenance": {
            "card_sha256": sha(card),
            "native_sha256": sha(native),
            "spectrum_sha256": sha(spectrum),
            "native_card_binding": "PASS",
        },
        "cascade": {
            "branching_fraction_spheno": cascade,
            "total_width_GeV_spheno": width,
            "partial_width_GeV_spheno": native_width,
            "partial_width_GeV_analytic_tree": analytic_width,
            "tree_width_relative_residual": residual,
            "listed_branching_fraction_sum": sum(br for br, _ in channels),
            "status": "TREE_TWO_BODY_WIDTH_CHECK_PASS_FULL_BR_UNCERTIFIED",
        },
        "neutrino": {
            "sum_masses_eV": msum,
            "cosmological_gate": "UNKNOWN_MODEL_DEPENDENT",
        },
        "lfv": {
            "mu_to_e_gamma": {
                "scalar_one_loop_subset_BR": subset["mu_to_e_gamma"],
                "spheno_flavorkit_BR": flavorkit_lfv[701],
                "flavorkit_to_scalar_subset_ratio":
                    flavorkit_lfv[701] / subset["mu_to_e_gamma"],
                "meg_ii_2025_upper_90CL": MEG_II_2025_LIMIT_90CL,
                "meg_ii_source": MEG_II_SOURCE,
                "scalar_subset_above_meg_ii_limit":
                    subset["mu_to_e_gamma"] > MEG_II_2025_LIMIT_90CL,
                "full_Zee_lfv_gate": "UNKNOWN",
                "reason": "Different amplitude subsets; missing terms unbounded",
            },
            "scalar_one_loop_subsets_BR": subset,
            "spheno_flavorkit": {
                "mu_to_3e": flavorkit_lfv[901],
                "mu_e_conversion_Al": flavorkit_lfv[800],
            },
            "flavorkit_model_and_loop_scope": "UNCERTIFIED",
        },
        "flavor": {
            "spheno_BR_B_to_Xs_gamma": flavorkit_qfv[200],
            "spheno_ratio_to_internal_SM": flavorkit_qfv[201],
            "physical_flavor_gate": "UNKNOWN",
        },
        "ewpo": {
            "S_spheno": ew[2], "T_spheno": ew[1], "U_spheno": ew[3],
            "correlated_fit_gate": "UNKNOWN",
        },
        "higgs": {
            "pole_mass_GeV": masses[25],
            "total_width_GeV_spheno": decay_block(text, 25)[0],
            "rates_and_searches_gate": "UNKNOWN",
        },
        "remaining_gates": {
            key: "UNKNOWN" for key in (
                "full_Zee_LFV", "flavor_validated", "Higgs_rates_and_searches",
                "unitarity", "global_vacuum", "EWPO_fit", "cosmology")
        },
        "globally_allowed": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "status": "FIXED_BENCHMARK_EVIDENCE_ONLY_NO_GLOBAL_ACCEPTANCE",
        "benchmarks": [evaluate(order) for order in ("NO", "IO")],
        "n_globally_allowed": None,
    }
    with args.output.open("x") as out:
        json.dump(report, out, indent=2, allow_nan=False)
        out.write("\n")
    for entry in report["benchmarks"]:
        print(entry["ordering"], "cascade BR:",
              f'{entry["cascade"]["branching_fraction_spheno"]:.7f}',
              "partial mu->e gamma BR:",
              f'{entry["lfv"]["mu_to_e_gamma"]["scalar_one_loop_subset_BR"]:.4e}')
    print("Physical acceptance: UNKNOWN")


if __name__ == "__main__":
    main()
