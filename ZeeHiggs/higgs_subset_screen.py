"""Published-method HiggsSignals reference *screen* for Zee model scans.

The HiggsTools 1.3 authors (arXiv:2608.05401v1, Sec 4.3) favor
Delta chi2 relative to a best-fit model, or (where justified) the SM
reference as an approximation. For a **properly profiled two-parameter**
inference, Delta chi2 > 5.991464547 is the asymptotic two-sided 95% CL
contour criterion (chi2_2), subject to Wilks regularity conditions.

This Zee pilot is NOT such a 2D profile fit; no Zee best-fit or
model-theory covariance was established. It is a robust exploratory
reference screen only and MUST NOT set a physical 95% CL exclusion.
DO NOT use chi2/159 or chi2-(SM chi2)>3.84 as a purported model 95% CL.

The four observed chi2 values are for the SAME official HiggsSignals
v1.1 dataset, native masses and model couplings, varying only the
prescribed h125 ggF normalization and BR method. Their spread is not
a probability distribution or an assigned systematic uncertainty.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/"Higgs/reference/s00017_HXSWG_ggF_width_validation_summary.json"
LOOP=ROOT/"Higgs/reference/s00017_h125_loop_preserved_summary.json"
HIGGS_BOUNDS=ROOT/"Higgs/reference/s00017_higgstools_v1_3_executed.json"

PRESCRIPTIONS=(
    "hxswg_xs_effective_BRs",
    "spheno_xs_effective_BRs",
    "hxswg_xs_explicit_BRs",
    "spheno_xs_explicit_BRs",
)
PRIMARY=PRESCRIPTIONS[0]
TWO_PARAMETER_95= -2.0*math.log(0.05)
# The authors note the SM proxy is useful if it reasonably describes
# the measurements. 159 is the official observableCount() of the
# exact frozen v1.1 dataset in the executed HS diagnostic. Not a fit ndf.
NOBS=159
EXPECTED_SOURCE={
    "artifact_id":14703179437,
    "original_reference_JSON_sha256":"38306e5b6ff40e21a32f72e9be349c705997cc6ea32ea7565a9659e24d5bb4ce",
    "original_width_diagnostic_JSON_sha256":"0917ab5b24feda3606e0c12c8ea38fa96c227327c25c8ed2e51b9eeebb30fa25",
    "higgsSignals_dataset":"v1.1",
    "higgstools_v":"1.3",
}
CITATION={
    "authors_method":"HiggsTools for LHC Run 3 and Beyond, arXiv:2608.05401, sec. 4.3",
    "url":"https://arxiv.org/abs/2608.05401",
    "criterion":"Delta chi2 = chi2(model) - chi2(best-fit); SM as proxy if adequate; 2D 95% 5.99",
    "assumption":"Two parameters of interest and all other relevant parameters fixed or properly profiled, plus Wilks regularity",
    "actual_Zee_pilot_satisfies_two_parameter_profile_assumption":False,
    "SM_proxy_is_verified_model_best_fit":False,
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(x, label):
    if isinstance(x,bool) or not isinstance(x,(float,int)):
        raise ValueError(f"{label} is not a numeric chi2")
    x=float(x)
    if not math.isfinite(x) or x < 0:
        raise ValueError(f"{label} is negative or nonfinite")
    return x


def screen_point(point_id, sm_chi2, corner_chi2, *, threshold=TWO_PARAMETER_95):
    """Pure numerical exploratory screen, usable by a future pilot scan.

    Never returns a physical pass/fail; only an SM-referenced diagnostic.
    For arbitrary scan points, the caller MUST check that the same data,
    particle mass treatment, loop conventions and provenance apply.
    """
    if not isinstance(point_id,str) or not point_id:
        raise ValueError("Empty point ID")
    if set(corner_chi2)!=set(PRESCRIPTIONS):
        raise ValueError("Exactly four predefined HiggsSignal prescriptions required")
    if not math.isclose(threshold,TWO_PARAMETER_95,rel_tol=0,abs_tol=1e-12):
        raise ValueError("No silent change to the declared two-parameter reference threshold")
    sm=finite(sm_chi2,"SM reference chi2")
    corners={k:finite(corner_chi2[k],k) for k in PRESCRIPTIONS}
    deltas={k:corners[k]-sm for k in PRESCRIPTIONS}
    dmin=min(deltas.values())
    dmax=max(deltas.values())
    # Guard on ALL four actual methods, rather than picking whichever
    # prescription makes a point pass or fail.
    if dmin>TWO_PARAMETER_95:
        flag="HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL"
    elif dmax<=TWO_PARAMETER_95:
        flag="HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_PASS"
    else:
        flag="HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_BORDERLINE"
    return {
        "point_id":point_id,
        "status":flag,
        "conditional_2d_SM_reference_threshold":TWO_PARAMETER_95,
        "chi2_SM_same_mass_and_dataset":sm,
        "chi2_by_prescription":corners,
        "delta_chi2_vs_SM_by_prescription":deltas,
        "primary_prescription":PRIMARY,
        "primary_delta_chi2":deltas[PRIMARY],
        "min_delta_chi2":dmin,
        "max_delta_chi2":dmax,
        "all_prescriptions_above_threshold":dmin>TWO_PARAMETER_95,
        "all_prescriptions_below_or_equal_threshold":dmax<=TWO_PARAMETER_95,
        "reference_chi2_gap_to_clear_screen":dmin-TWO_PARAMETER_95,
        "interpretation":(
            "CONDITIONAL two-dimensional SM-reference screening ONLY; not "
            "a profiled Zee best-fit test, not a p-value and not a "
            "validated 95 percent physical exclusion."
        ),
        "is_formal_95CL_model_exclusion":False,
        "higgsSignals_physical_95CL":None,
        "higgs_theory_covariance_validated":False,
        "physical_higgs":"UNKNOWN",
        "globally_allowed":None,
    }


def validate_receipts(reference, loop, hb):
    if (reference.get("status") !=
            "ACTUAL_OFFICIAL_HIGGSTOOLS_XSWG_SPHENO_FOUR_CORNER_VALIDATION_NO_PHYSICAL_EXCLUSION"):
        raise ValueError("Not the executed HXSWG four-prescription source")
    if not reference.get("all_chi2_same_HS_dataset_and_spectrum"):
        raise ValueError("HS measurements/particle inputs not common")
    if (reference.get("physical_higgs")!="UNKNOWN"
        or reference.get("theory_sigma") is not None
        or reference.get("higgsSignals_95CL_exclusion") is not None
        or reference.get("globally_allowed") is not None):
        raise ValueError("Previous source was promoted to physical certainty")
    for key,expected in EXPECTED_SOURCE.items():
        if reference.get("source",{}).get(key)!=expected:
            raise ValueError("Source provenance mismatch: "+key)
    if (loop.get("status") !=
            "ACTUAL_OFFICIAL_HIGGSTOOLS_V1_3_LOOP_PRESERVED_DIAGNOSTIC_NO_PHYSICAL_EXCLUSION"
        or loop.get("source",{}).get("flags")!={"calcggH":False,"calcHgamgam":False}
        or loop.get("phys_higgs")!="UNKNOWN"):
        raise ValueError("Missing validated full-Zee native hgg/hgammagamma loop policy")
    if (hb.get("status")!="HIGGSTOOLS_EXECUTED_BUT_PHYSICS_INCOMPLETE"
        or hb.get("physical_higgs")!="UNKNOWN"):
        raise ValueError("Invalid historic HiggsBounds diagnostic")
    if set(reference.get("results",{}))!={"NO","IO"}:
        raise ValueError("Expected both native orderings")
    loop_by_order={p["ordering"]:p for p in loop["points"]}
    hb_by_order={p["point_id"].split("_")[-1]:p for p in hb["points"]}
    if set(loop_by_order)!={"NO","IO"} or set(hb_by_order)!={"NO","IO"}:
        raise ValueError("Missing loop/HB control")
    for order,row in reference["results"].items():
        old=loop_by_order[order]
        if not math.isclose(row["HiggsSignals_chi2"]["spheno_xs_effective_BRs"],
                            old["corrected_effective_chi2"],rel_tol=0,abs_tol=1e-9):
            raise ValueError("Mismatch actual full-Zee-loop input "+order)
        if not math.isclose(row["HiggsSignals_chi2"]["spheno_xs_explicit_BRs"],
                            old["explicit_spheno_chi2"],rel_tol=0,abs_tol=1e-9):
            raise ValueError("Mismatch old SPheno explicit BR input "+order)
        if not math.isclose(row["higgsSignals_SM_chi2"],old["sm_chi2"],
                            rel_tol=0,abs_tol=1e-9):
            raise ValueError("SM reference changed "+order)
        if not math.isclose(row["HXSWG_normalized_Zee_ggF_pb"],
                            row["HXSWG_SM_ggF_pb"]*row["kappa_g2"],
                            rel_tol=5e-9):
            raise ValueError("Missing source-Zee loop closure in ggF "+order)
        if set(row["HiggsSignals_chi2"])!=set(PRESCRIPTIONS):
            raise ValueError("Incomplete four-corner HiggsSignals study")
        if not row["HiggsSignals_chi2_delta_vs_SM"] or set(
                row["HiggsSignals_chi2_delta_vs_SM"])!=set(PRESCRIPTIONS):
            raise ValueError("Missing archived relative chi2")
        for k in PRESCRIPTIONS:
            if not math.isclose(row["HiggsSignals_chi2"][k]-
                                row["higgsSignals_SM_chi2"],
                                row["HiggsSignals_chi2_delta_vs_SM"][k],
                                rel_tol=0,abs_tol=1e-9):
                raise ValueError("Saved chi2 difference does not close")
    return hb_by_order


def selected_higgsbounds_diagnostic(entry):
    """Historical raw HB neutral-only selected limits, OLD production scheme.

    Its input is NOT the preferred HXSWG-normalized h125 production rate,
    so this is explicitly NOT a synchronized HB subset PASS/FAIL.
    """
    if (entry.get("physical_higgs")!="UNKNOWN"
        or entry.get("higgsbounds_full_model_constraint")!="UNKNOWN"):
        raise ValueError("Historical HiggsBounds report promoted to full-model acceptance")
    branch=entry["branches"]["effective_coupling_h125_BR"]
    selected=branch["selectedLimits"]
    if not selected:
        raise ValueError("No HiggsBounds selected analyses recorded")
    ps={x["particle"] for x in selected}
    if not ps.issubset({"25","35","36"}):
        raise ValueError("Unexpected selected charged-Higgs limit; update validated subset")
    raw=bool(branch["HiggsBounds_allowed_with_incomplete_inputs"])
    if any(finite(x["obs_ratio"], "observed limit ratio")>1 for x in selected) and raw:
        raise ValueError("HB raw status inconsistent with a selected limit")
    return {
        "status":("HB_HISTORIC_NEUTRAL_SELECTED_NOT_EXCLUDED" if raw
                  else "HB_HISTORIC_SELECTED_EXCLUSION_DIAGNOSTIC"),
        "observed_ratio_max":max(x["obs_ratio"] for x in selected),
        "selected_particle_ids":sorted(ps),
        "charged_search_prediction_coverage":False,
        "selected_limit_input_matches_preferred_HXSWG_ggF":False,
        "is_full_Zee_HiggsBounds_pass":False,
        "physical_HiggsBounds_constraint":"UNKNOWN",
    }


def evaluate_archived(reference_path=REFERENCE, loop_path=LOOP, hb_path=HIGGS_BOUNDS):
    ref=json.loads(Path(reference_path).read_text())
    loop=json.loads(Path(loop_path).read_text())
    hb=json.loads(Path(hb_path).read_text())
    hb_by_order=validate_receipts(ref,loop,hb)
    results=[]
    for order in ("NO","IO"):
        data=ref["results"][order]
        row=screen_point("s00017_"+order,data["higgsSignals_SM_chi2"],
                         data["HiggsSignals_chi2"])
        row["HiggsBounds_historic_diagnostic"]=selected_higgsbounds_diagnostic(
            hb_by_order[order])
        row["measured_observable_count"]=NOBS
        row["chi2_SM_div_observable_count"]=data["higgsSignals_SM_chi2"]/NOBS
        row["SM_reference_is_approximated_best_fit"]=True
        row["SM_proxy_is_verified_Zee_best_fit"]=False
        row["profile_likelihood_for_full_Zee_performed"]=False
        row["selection_cannot_be_used_as_global_Higgs_cut"]=True
        results.append(row)
    return {
        "status":"HIGGS_SCOPED_REFERENCE_CHI2_SCREEN_VALIDATED_NO_FULL_95CL",
        "statistical_method":CITATION,
        "reference_data":{
            "input":str(reference_path),
            "input_sha256":sha256(reference_path),
            "loop_input":str(loop_path),
            "loop_input_sha256":sha256(loop_path),
            "higgsbounds_history":str(hb_path),
            "higgsbounds_history_sha256":sha256(hb_path),
            "HS_dataset_version":"v1.1",
            "HS_observable_count_from_original_executed_output":NOBS,
            "HS_observable_count_full_JSON_source_sha256":EXPECTED_SOURCE[
                "original_width_diagnostic_JSON_sha256"],
            "HiggsTools_version":"v1.3",
            "model_specific_theory_covariance":"NOT_VALIDATED",
            "Zee_best_fit":"NOT_YET_ESTABLISHED",
        },
        "two_dimensional_95pct_threshold":TWO_PARAMETER_95,
        "criteria_have_true_95pct_coverage_for_current_Zee_scan":False,
        "scanned_parameter_dimension_known":False,
        "model_nuisances_profiled":False,
        "points":results,
        "physical_higgs":"UNKNOWN",
        "higgsSignals_95CL_exclusion":None,
        "higgsBounds_full_model":"UNKNOWN",
        "globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,default=REFERENCE)
    p.add_argument("--loop",type=Path,default=LOOP)
    p.add_argument("--higgsbounds-history",type=Path,default=HIGGS_BOUNDS)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=evaluate_archived(a.reference,a.loop,a.higgsbounds_history)
    with a.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    for row in result["points"]:
        print(row["point_id"],row["status"],
              f"primary_delta_chi2={row['primary_delta_chi2']:.6f}",
              f"least_disfavored_delta_chi2={row['min_delta_chi2']:.6f}",
              f"2D_reference_criterion={TWO_PARAMETER_95:.6f}",
              "physical=UNKNOWN")


if __name__=="__main__":
    main()
