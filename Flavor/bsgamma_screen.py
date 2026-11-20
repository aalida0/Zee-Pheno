"""Conservative, explicitly NON-EXCLUSION B->Xs gamma screening.

Inputs: source-bound *executed* Zee HEPfit NLO result from Flavor/bsgamma.py.
Two diagnostic predictions:
  (A) published 2026 NNLO SM centre + same-scheme HEPfit NLO NP shift;
  (B) 2026 published leading-order NP formula (real LO Delta C7,8 at 160 GeV).

The LO-vs-NLO spread is NOT a statistical/theoretical uncertainty or a
certified remainder bound. One spread is used as a user-declared
*screening guard*, solely to prioritize scan candidates. All physical
95% CL confidence levels, and full-model acceptances, remain UNKNOWN.
No physical veto is authorized by any screening status.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
S00017 = HERE / "reference" / "hepfit_restored" / "s00017" / "bsgamma_nlo.json"
ARCHIVE = ROOT / "Scan" / "reference" / "native_benchmark"

REFERENCE = {
    "paper": "Misiak et al., arXiv:2608.29863v1 (2026), Eqs.(1),(3),(4)",
    "url": "https://arxiv.org/abs/2608.29863",
    "observable": "CP- and isospin-averaged BR(B->Xs gamma), Egamma>1.6 GeV",
    "units": "1e-4",
    "sm_central": 3.54,
    "sm_sigma": 0.14,
    "exp_central": 3.49,
    "exp_sigma": 0.19,
    "LO_delta_c7_coefficient": -8.51,
    "LO_delta_c8_coefficient": -2.16,
    "published_2HDM_II_bound_GeV": 670.,
    "published_2HDM_bound_reproduced_with_this_prescription": False,
}
# Only for a nominal Gaussian display, not the publication's 95% CL belt.
N_SIGMA = 1.96
# One non-probabilistic method spread as an engineering/triage guard.
# NOT an estimated 1-sigma NP theory error or confidence limit.
GUARD_SPREAD_MULTIPLIER = 1.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(name, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} not numeric")
    num = float(value)
    if not math.isfinite(num):
        raise ValueError(f"{name} not finite")
    return num


def screen_prediction(np_shift, c7, c8, *, label=None):
    """Evaluate conditional scan triage from already validated NP inputs.

    *np_shift*, c7 and c8 are dimensionless BR and ordinary LO BMU/NDR Wilson
    coefficients at mu0=160 GeV. Caller must verify that source contract.
    """
    shift4 = 1e4 * finite("np_shift", np_shift)
    dc7 = finite("Delta C7(0)", c7)
    dc8 = finite("Delta C8(0)", c8)
    if not label:
        raise ValueError("Missing point identifier")
    ref = REFERENCE
    hybrid = ref["sm_central"] + shift4
    lo_shift4 = ref["LO_delta_c7_coefficient"] * dc7 + ref["LO_delta_c8_coefficient"] * dc8
    lo = ref["sm_central"] + lo_shift4
    sigma = math.hypot(ref["exp_sigma"], ref["sm_sigma"])
    low, high = (ref["exp_central"] - N_SIGMA*sigma,
                 ref["exp_central"] + N_SIGMA*sigma)
    method_spread = abs(hybrid - lo)
    method_range = [min(hybrid, lo), max(hybrid, lo)]
    guard = [method_range[0] - GUARD_SPREAD_MULTIPLIER * method_spread,
             method_range[1] + GUARD_SPREAD_MULTIPLIER * method_spread]
    if guard[0] > high or guard[1] < low:
        decision = "B_SCREEN_FAIL"
    elif guard[0] >= low and guard[1] <= high:
        decision = "B_SCREEN_PASS"
    else:
        decision = "B_SCREEN_BORDERLINE"
    z_hybrid = (hybrid - ref["exp_central"])/sigma
    z_lo = (lo - ref["exp_central"])/sigma
    if hybrid > high:
        correction = high - hybrid
        relation = "decrease_to_upper_edge"
    elif hybrid < low:
        correction = low - hybrid
        relation = "increase_to_lower_edge"
    else:
        correction = None
        relation = "hybrid_inside_nominal_interval"
    return {
        "point_id": label,
        "screen_status": decision,
        "physical_constraint": "B_PHYSICAL_UNKNOWN",
        "physically_allowed": None,
        "experimental_exclusion_at_95CL": None,
        "NP_theory_sigma": None,
        "unit_for_rates": "1e-4",
        "HEPfit_NLO_delta_BR": np_shift,
        "modern_SM_plus_HEPfit_NLO_delta_BR": hybrid,
        "published_2026_linear_LO_delta_BR_in_1e-4": lo_shift4,
        "published_2026_linear_LO_total_in_1e-4": lo,
        "formal_statistical_interpretation": "NONE_FOR_SCREEN",
        "illustrative_gaussian_assumption": "independent exp and quoted SM sigmas; NP uncertainty omitted",
        "illustrative_gaussian_1sigma_in_1e-4": sigma,
        "illustrative_gaussian_95_interval_in_1e-4": [low, high],
        "illustrative_z_hybrid": z_hybrid,
        "illustrative_z_linear": z_lo,
        "method_spread_in_1e-4": method_spread,
        "method_spread_is_theory_uncertainty": False,
        "guard_definition": "method_range +/- 1*|hybrid - linear_LO|, triage-only",
        "method_range_in_1e-4": method_range,
        "guarded_range_in_1e-4": guard,
        "shift_needed_to_enter_nominal_95_band_in_1e-4": correction,
        "shift_to_band_direction": relation,
        "fraction_of_NP_shift_needed_to_enter_band":
            abs(correction / shift4) if (correction is not None and shift4 != 0.) else None,
        "reason": (
            "HEPfit NLO delta and 2026 linear LO delta are different "
            "perturbative approximations, not a statistically calibrated "
            "theory nuisance. The enlarged interval is only a priority "
            "ranking screen. No point may be removed from the physical "
            "scan solely because it fails this heuristic."
        ),
    }


def validate_archived_data(data, archive=ARCHIVE):
    if not isinstance(data, dict) or data.get("status") != "EXECUTED_CENTRAL_NLO_ONLY_NO_PHYSICAL_EXCLUSION":
        raise ValueError("Missing executed HEPfit NLO archive")
    for k, value in (("physical_bsgamma", "UNKNOWN"), ("theory_sigma", None),
                     ("experimental_exclusion", None), ("globally_allowed", None)):
        if k not in data or data[k] != value:
            raise ValueError("Archived constraints were promoted: " + k)
    if (data.get("matching_scale_GeV") != 160.0 or data.get("photon_cut_GeV") != 1.6):
        raise ValueError("Unsupported matching scale/cut")
    provenance = data.get("source", {})
    restored = provenance.get("restored_backend", {})
    if (restored.get("status") != "RESTORED_HEPFIT_NUMERICALLY_EQUIVALENT"
            or restored.get("validated_runs") != 37
            or restored.get("HEPfit_commit") != "322bd73a2dadd94ce94a6d71cb79d616e569120a"
            or provenance.get("driver_binary_sha256") != restored.get("binary_sha256")
            or provenance.get("matching_sha256") != restored.get("matching_sha256")):
        raise ValueError("Missing validated original HEPfit provenance")
    points = data.get("points")
    if not isinstance(points, list) or len(points) != 2 or {r.get("ordering") for r in points} != {"NO", "IO"}:
        raise ValueError("Expected two fixed s00017 native controls")
    for point in points:
        order = point["ordering"]
        scalar = point["scalar"]
        native_dir = archive / ("s00017_" + order + "_benchmark")
        for path, field in (("input.in", "card_sha256"), ("native.json", "native_sha256"),
                            ("spectrum.slha", "spectrum_sha256")):
            if sha(native_dir/path) != scalar[field]:
                raise ValueError("Mismatch with archived native benchmark: " + order + "/" + path)
        if (point.get("physical_bsgamma_gate") != "UNKNOWN" or
                point.get("uncertainty_sigma") is not None or
                point.get("full_model_allowed") is not None):
            raise ValueError("Point incorrectly promoted to a physical constraint")
        br = finite("HEPfit NLO", point["BR_NLO_Egamma_gt_1p6GeV"])
        sm = finite("HEPfit same scheme SM", point["BR_SM_same_scheme"])
        shift = finite("HEPfit delta", point["delta_BR_same_scheme"])
        if abs(br-sm-shift) > 1e-14:
            raise ValueError("HEPfit NP subtraction does not close")
        if point.get("value_checks", {}).get("high_scale_injected_coefficients") != "PASS":
            raise ValueError("High-scale NP matching not verified")
        # One ordering cannot covertly borrow the other's rate or scalars.
        m = point["matching"]
        for name in ("C7_LO","C8_LO","C4_NLO_raw","C7_NLO_raw","C8_NLO_raw"):
            finite("matching:"+name, m[name])


def evaluate_archive(path=S00017):
    archive = Path(path)
    data = json.loads(archive.read_text())
    validate_archived_data(data)
    results = [screen_prediction(
        p["delta_BR_same_scheme"], p["matching"]["C7_LO"],
        p["matching"]["C8_LO"], label="s00017_" + p["ordering"])
        for p in data["points"]]
    return {
        "status": "TWO_METHOD_CONSERVATIVE_EXPLORATORY_SCREEN_ONLY",
        "source_input": {"filename":str(archive), "sha256":sha(archive)},
        "reference": REFERENCE,
        "nominal_confidence_interval_is_validated_95CL": False,
        "reason_not_a_95CL": (
            "The Gaussian comparison uses only the experimental and quoted "
            "modern-SM uncertainties, with no validated model-specific NP "
            "error or correlated likelihood. No 670 GeV Type-II confidence "
            "belt has been reproduced for this prescription."
        ),
        "physical_bsgamma": "UNKNOWN",
        "theory_sigma": None,
        "experimental_exclusion": None,
        "globally_allowed": None,
        "points": results,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=S00017)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = evaluate_archive(args.input)
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    for entry in result["points"]:
        print(entry["point_id"], entry["screen_status"],
              f"hybrid={entry['modern_SM_plus_HEPfit_NLO_delta_BR']:.6f}e-4",
              f"linear={entry['published_2026_linear_LO_total_in_1e-4']:.6f}e-4",
              "physical=UNKNOWN")


if __name__ == "__main__":
    main()
