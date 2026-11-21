"""Per-native-SPheno-card scoped HiggsTools screens for future Zee scans.

The same schema_version=1 card/native/spectrum manifest as Scan.native_pipeline
is used. Source binding is delegated to the full independent Scan.theory_gate_native
checker; model-specific Higgs couplings, decays, rates and table range are always
taken from *this point's* spectrum, never from the s00017 reference.

Reports official HiggsSignals chi2 and the HiggsBounds selected subset using
HXSWG-normalized h125 ggF and official 13TeV Type-II charged production tables.
No invented model best-fit, theory covariance, 95% CL or global viability.
"""
import argparse
import json
import math
from pathlib import Path

from Scan.theory_gate_native import screen_one, ORIGINAL_BINARY_SHA256
from ZeeHiggs.spheno_audit import read_slha, sha256
from ZeeHiggs.higgstools_pilot import result_limits, dataset_receipt
from ZeeHiggs.diagnose_h125 import make_inputs,predictions_for
from ZeeHiggs.charged_effective_production import production_reference_official
from ZeeHiggs.higgs_subset_screen import screen_point


def audit_manifest_source(entry,root):
    """Authenticate raw card, spectrum and native mass table before Higgs input."""
    t=screen_one(entry,root,expected_binary_sha=ORIGINAL_BINARY_SHA256)
    p=(Path(root)/entry["spectrum"]).resolve()
    blocks,decays=read_slha(p.read_text())
    from ZeeHiggs.spheno_audit import simple_block
    masses=simple_block(blocks,"MASS")
    neutral={}
    for pdg in (25,35,36):
        key=(pdg,)
        if key not in masses or masses[key]<=0:
            raise ValueError("Missing neutral scalar mass")
        neutral[str(pdg)]=masses[key]
    charged={}
    for pdg in (37,9000037):
        key=(pdg,)
        if key not in masses or masses[key]<=0:
            raise ValueError("Missing charged scalar mass")
        charged[str(pdg)]=masses[key]
    for pdg in (6,25,35,36,37,9000037):
        record=decays.get(pdg)
        if record is None or record["total_width_GeV"]<=0:
            raise ValueError("Missing native top/scalar width "+str(pdg))
        br=sum(channel["br"] for channel in record["channels"])
        if abs(br-1)>1e-3:
            raise ValueError("Source native total BR closure failed "+str(pdg))
    xs={}
    for row,_ in blocks.get("HIGGSLHC13",[]):
        if len(row)<3:raise ValueError("Malformed HiggsLHC13 production row")
        if int(row[0])==1 and int(row[1]) in (25,35,36):
            pdg=int(row[1])
            if str(pdg) in xs:raise ValueError("Duplicate scalar gluon fusion rate")
            value=float(row[2].replace("D","E"))
            if not math.isfinite(value) or value<=0:
                raise ValueError("Invalid SPheno approximate neutral ggF rate")
            xs[str(pdg)]=value
    if set(xs)!={"25","35","36"}:
        raise ValueError("Incomplete native scalar ggF reference")
    return {
        "point_id":entry["point_id"],"ordering":entry["ordering"],
        "source":{
            "card_sha256":t["provenance"]["card_sha256"],
            "spectrum_sha256":t["provenance"]["spheno_spectrum_sha256"],
            "spheno_binary_sha256":ORIGINAL_BINARY_SHA256,
        },
        "spectrum_path":str(p),"masses_GeV":{**neutral,**charged},
        "neutral_ggF_LHC13_pb":xs,
        "theory_scope_pass":t["passed_in_declared_tree_and_one_loop_scalar_scope"],
        "theory_scope_status":t["status"],
    }


def evaluate_one(entry,root,bounds,signals,hinput):
    import Higgs.predictions as HP
    report=audit_manifest_source(entry,root)
    data,ids,removed,old=make_inputs(report,hinput)
    charged_status="TABULATED_13TEV_TYPE_II_INPUT"
    try:
        charged=production_reference_official(report)
        for pdg in (37,9000037):
            k=f"cxn_{ids[pdg]}_Hpmtb_LHC13"
            if k in data:
                raise ValueError("Refusing to overwrite actual charged production data")
            data[k]=charged["higgstools_cross_sections"][f"cxn_{pdg}_Hpmtb_LHC13"]
    except ValueError as exc:
        if "Outside justified HiggsTools Hpmtb tabulation/top decay domain" not in str(exc):
            raise
        charged=None
        charged_status="UNKNOWN_OUTSIDE_CHARGED_TABLE_DOMAIN"
    # Keep existing heavy H/A rates and Zee-induced hgg/hgamgam loop amplitudes.
    # Remove the historic h125 ggF so the HXSWG SM*Zee kappa_g^2 reference
    # is used by HiggsTools while all neutral/heavy states remain the same.
    hkey=f"cxn_{ids[25]}_ggH_LHC13"
    if hkey not in data:raise ValueError("Missing original 13TeV ggF h125 override")
    old_ggf=data[hkey]
    hxswg_data=dict(data)
    hxswg_data.pop(hkey)
    # Four predeclared cross-section/BR prescriptions: never select
    # whichever gives a convenient HiggsSignals decision.
    corners={
        "hxswg_xs_effective_BRs":(hxswg_data,False),
        "spheno_xs_effective_BRs":(data,False),
        "hxswg_xs_explicit_BRs":(hxswg_data,True),
        "spheno_xs_explicit_BRs":(data,True),
    }
    preds={name:predictions_for(vals,ids,hinput,explicit=explicit)
           for name,(vals,explicit) in corners.items()}
    hs_corners={name:float(signals(pred)) for name,pred in preds.items()}
    pred=preds["hxswg_xs_effective_BRs"]
    hs=hs_corners["hxswg_xs_effective_BRs"]
    sm=HP.Predictions()
    h=sm.addParticle(HP.BsmParticle("SM_h",charge="neutral",cp="even"))
    h.setMass(report["masses_GeV"]["25"])
    HP.effectiveCouplingInput(h,HP.smLikeEffCouplings)
    sm_chi=float(signals(sm))
    if not all(math.isfinite(x) for x in [sm_chi,*hs_corners.values()]):
        raise ValueError("Nonfinite official HiggsSignals chi2")
    reference_screen=screen_point(entry["point_id"],sm_chi,hs_corners)
    hb=bounds(pred)
    selected=result_limits(hb)
    observed=[float(row["obs_ratio"]) for row in selected]
    subset_status=("CENTRAL_RATE_EVALUATED_SUBSET_EXCLUSION_CANDIDATE"
                   if any(v>1 for v in observed)
                   else "NOT_EXCLUDED_WITHIN_EVALUATED_SUBSET" if observed
                   else "NO_SELECTED_LIMITS_FOR_INPUTS")
    return {
        "point_id":entry["point_id"],"ordering":entry["ordering"],
        "status":"SOURCE_BOUND_SCOPED_HIGGSTOOLS_SCREEN_NOT_FULL_ACCEPTANCE",
        "source":report["source"],
        "native_theory_gate_status":report["theory_scope_status"],
        "native_theory_gate_pass":report["theory_scope_pass"],
        "charged_reference_status":charged_status,
        "charged_production_reference":charged,
        "old_native_h125_ggF_pb":old_ggf,
        "HiggsSignals":{
            "model_chi2":hs,"SM_same_mass_chi2":sm_chi,
            "delta_chi2_vs_SM":hs-sm_chi,
            "four_predeclared_chi2_prescriptions":hs_corners,
            "tested_subset_reference_screen":reference_screen,
            "method":"Official HiggsSignals; HXSWG h125 ggF, Zee loop effective BR",
            "is_profiled_Zee_95CL_exclusion":False,
            "full_physical_HiggsSignals":"UNKNOWN",
        },
        "HiggsBounds":{
            "raw_allowed_with_incomplete_rates":bool(hb.allowed),
            "selectedLimits":selected,
            "tested_subset_central_input_status":subset_status,
            "tested_subset_is_full_Zee_exclusion":False,
            "all_analysis_topology_coverage_validated":False,
            "full_physical_HiggsBounds":"UNKNOWN",
        },
        "unchanged_native_HA_ggF_reference":{
            x:report["neutral_ggF_LHC13_pb"][x] for x in ("35","36")},
        "h125_loop_inputs_preserved":True,
        "uncalculated_Zgamma_keys_removed":removed,
        "full_higgs_status":"UNKNOWN",
        "globally_allowed":None,
    }


def evaluate_manifest(manifest_path,bounds,signals,hinput):
    manifest_path=Path(manifest_path).resolve()
    man=json.loads(manifest_path.read_text())
    if man.get("schema_version")!=1 or not isinstance(man.get("points"),list):
        raise ValueError("Expected schema_version=1 manifest")
    if not 1<=len(man["points"])<=1000:
        raise ValueError("Expected 1..1000 native manifest entries")
    ids=[r.get("point_id") for r in man["points"]]
    if len(ids)!=len(set(ids)):raise ValueError("Duplicate point_id")
    results=[]
    for entry in man["points"]:
        try:
            results.append(evaluate_one(entry,manifest_path.parent,bounds,signals,hinput))
        except (ValueError,KeyError,IndexError,TypeError,ArithmeticError,
                OSError,AssertionError) as exc:
            results.append({
                "point_id":entry.get("point_id"),"ordering":entry.get("ordering"),
                "status":"NOT_EVALUATED_SOURCE_OR_BACKEND_ERROR",
                "error":str(exc),"globally_allowed":None,
            })
    return {
        "status":"SCOPED_NATIVE_HIGGSTOOLS_HIGGS_SCAN_NOT_FULL_EXCLUSION",
        "manifest_sha256":sha256(manifest_path),
        "n_points":len(results),
        "n_scoped_evaluated":sum(r["status"]==
            "SOURCE_BOUND_SCOPED_HIGGSTOOLS_SCREEN_NOT_FULL_ACCEPTANCE"
            for r in results),
        "points":results,"globally_allowed":None,
    }


def main():
    import Higgs.bounds as HB
    import Higgs.signals as HS
    import Higgs.tools.Input as HI
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--hbdataset",type=Path,required=True)
    parser.add_argument("--hsdataset",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    bounds=HB.Bounds(str(args.hbdataset.resolve()))
    signals=HS.Signals(str(args.hsdataset.resolve()))
    out=evaluate_manifest(args.manifest,bounds,signals,HI)
    out["backend"]="Official HiggsTools v1.3"
    out["hbdataset"]=dataset_receipt(args.hbdataset)
    out["hsdataset"]=dataset_receipt(args.hsdataset)
    with args.output.open("x") as f:
        json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
    for p in out["points"]:
        print(p["point_id"],p["status"],
              p.get("charged_reference_status","UNKNOWN"),
              "physical=UNKNOWN")


if __name__=="__main__":
    main()
