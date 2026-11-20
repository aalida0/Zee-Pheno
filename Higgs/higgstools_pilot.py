"""Optional official HiggsTools 1.3 pilot; all results are DIAGNOSTIC.

Requires an independently installed HiggsTools >=v1.3 Python interface and
HB/HS JSON datasets. This module does not approximate the experimental data,
modify 2HDMC/ScannerS, or manufacture missing charged production. A raw
HiggsBounds verdict from incomplete model predictions is *not* a Zee exclusion.
"""
import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

from Higgs.spheno_audit import audit, sha256

NEUTRAL = (25, 35, 36)
CHARGED = (37, 9000037)


def infer_ids(data, report):
    """Use the actual official readHB5SLHA dictionary, not assumed h1/Hp1 IDs."""
    pairs = []
    for k, v in data.items():
        if not isinstance(k,str) or not k.startswith("m_") or not isinstance(v,(int,float)):
            continue
        if math.isfinite(v) and v > 0:
            pairs.append((k[2:],float(v)))
    ids = {}
    for pdg in NEUTRAL + CHARGED:
        mass = report["masses_GeV"][str(pdg)]
        matched = [key for key, v in pairs if abs(v - mass) <= max(1e-4, mass*1e-6)]
        if len(matched) != 1:
            raise ValueError(f"HiggsTools ambiguous/unmapped PDG {pdg}, mass {mass}: {matched}; m_ keys={pairs}")
        ids[pdg] = matched[0]
    if len(set(ids.values())) != len(ids):
        raise ValueError("HiggsTools collapsed distinct physical Zee scalars")
    return ids


def clean_missing_zero_couplings(data):
    """Do not insert NOT-CALCULATED SPheno hZgamma zeros into HiggsTools."""
    new = dict(data)
    removed = []
    for k in data:
        if isinstance(k,str) and re.search(r"zgamma|zgam|gamz",k,re.IGNORECASE):
            removed.append(k)
            new.pop(k)
    return new,removed


def dataset_receipt(path):
    path=Path(path).resolve()
    if not path.is_dir():
        raise FileNotFoundError("Missing official HiggsTools dataset: " + str(path))
    files = sorted(path.rglob("*.json"))
    if not files:
        raise ValueError("No HiggsTools JSON datasets at "+str(path))
    manifest = hashlib.sha256()
    for f in files:
        manifest.update(f.relative_to(path).as_posix().encode()+b"\0")
        manifest.update(bytes.fromhex(sha256(f)))
    result={"path":str(path),"n_json_files":len(files),"manifest_sha256":manifest.hexdigest()}
    if (path/".git").exists():
        try:
            p=subprocess.run(["git","-C",str(path),"rev-parse","HEAD"],
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                             text=True,check=True,timeout=12)
            result["git_head"]=p.stdout.strip()
        except (subprocess.CalledProcessError,subprocess.TimeoutExpired,FileNotFoundError):
            result["git_head"]=None
    return result


def result_limits(hb_result):
    out=[]
    for particle,limit in hb_result.selectedLimits.items():
        selected = limit.limit()
        record={
            "particle":str(particle),
            "obs_ratio":float(limit.obsRatio()),
            "exp_ratio":float(limit.expRatio()),
        }
        for attr in ("reference","citeKey","processDesc","experiment","collider"):
            try:
                fn=getattr(selected,attr)
                record[attr]=str(fn())
            except (AttributeError,TypeError,ValueError):
                record[attr]=None
        if not all(math.isfinite(record[k]) for k in ("obs_ratio","exp_ratio")):
            raise ValueError("Invalid HiggsBounds selected limit ratio")
        out.append(record)
    return out


def run_one(order, bounds, signals):
    report=audit(order)
    try:
        import Higgs.tools.Input as hinput
    except ImportError as e:
        raise RuntimeError("HiggsTools Python interface unavailable; install official HiggsTools v1.3") from e
    data=hinput.readHB5SLHA(
        report["spectrum_path"],
        neutralPDGs=list(NEUTRAL),
        chargedPDGs=list(CHARGED),
    )
    if not isinstance(data,dict):
        raise ValueError("HiggsTools SLHA importer did not return a dictionary")
    ids=infer_ids(data,report)
    data,missing=clean_missing_zero_couplings(data)
    # The historical SPheno native HiggsLHC13 contains ggF ONLY.
    # Use its actual 13 TeV numeric prediction as an explicit override:
    # this is approximate (SusHi1.5 rescaling) and NOT precision-accepted.
    overrides={}
    for pdg in NEUTRAL:
        k=f"cxn_{ids[pdg]}_ggH_LHC13"
        if k in data:
            observed=data[k]
            if not isinstance(observed,(int,float)) or not math.isfinite(observed):
                raise ValueError("Unexpected pre-existing ggH cross section")
        overrides[k]=report["neutral_ggF_LHC13_pb"][str(pdg)]
        data[k]=overrides[k]
    # Do not call it a full-production prediction: associated charged
    # tH+ and two-charged cascade process rates are NOT supplied.
    neutral_ids=[ids[p] for p in NEUTRAL]
    charged_ids=[ids[p] for p in CHARGED]
    # Compare two documented HiggsTools input prescriptions.
    # SM-like h effective-coupling BRs avoid SPheno width precision bias;
    # the explicit-all mode preserves native BRs as a diagnostic.
    all_explicit=hinput.predictionsFromDict(
        data, neutral_ids, charged_ids, [], useExplicitBr=True)
    effective_h125={x:True for x in neutral_ids+charged_ids}
    effective_h125[ids[25]]=False
    recommended=hinput.predictionsFromDict(
        data, neutral_ids, charged_ids, [], useExplicitBr=effective_h125)
    bounds_all=bounds(all_explicit)
    bounds_recommended=bounds(recommended)
    hs_explicit=float(signals(all_explicit))
    hs_recommended=float(signals(recommended))
    if not math.isfinite(hs_explicit) or not math.isfinite(hs_recommended):
        raise ValueError("Nonfinite HiggsSignals chi2")
    return {
        "point_id":report["point_id"],
        "status":"OFFICIAL_HIGGSTOOLS_DIAGNOSTIC_ONLY",
        "source":report["source"],
        "PDG_to_HiggsTools_ID":{str(k):v for k,v in ids.items()},
        "higgstools_imported_input_keys":sorted(str(k) for k in data),
        "uncalculated_Zgamma_keys_removed":missing,
        "ggH_SPheno_13TeV_pb_override":overrides,
        "branches":{
            "all_SPheno_explicit_BR":{
                "chi2_HiggsSignals":hs_explicit,
                "HiggsBounds_allowed_with_incomplete_inputs":bool(bounds_all.allowed),
                "selectedLimits":result_limits(bounds_all),
            },
            "effective_coupling_h125_BR":{
                "chi2_HiggsSignals":hs_recommended,
                "HiggsBounds_allowed_with_incomplete_inputs":bool(bounds_recommended.allowed),
                "selectedLimits":result_limits(bounds_recommended),
            }
        },
        "physical_higgs":"UNKNOWN",
        "higgssignals_confidence":"UNKNOWN",
        "higgsbounds_full_model_constraint":"UNKNOWN",
        "globally_allowed":None,
        "missing":[
            "Charged-Higgs production rates absent (must not default to physical zero)",
            "Source-validated h125 Zgamma and total-width precision missing",
            "Neutral ggF input rescaled from historical SusHi 1.5.0; QCD scheme not independently certified",
            "Not all nonstandard lepton Yukawa and exotic topologies are proven covered",
            "HB dataset's charged-cascade H2 -> H1 Z topology coverage not audited",
            "No validated HiggsSignals goodness-of-fit/p-value and correlated theory uncertainties",
        ],
        "note":"Raw HiggsTools verdict is not a full Zee 95%-CL exclusion when production and widths are incomplete."
    }


def calculate(orders,hbdataset,hsdataset):
    imports={
        "HiggsTools":"gitlab.com/higgsbounds/higgstools v1.3",
        "HBdataset":dataset_receipt(hbdataset),
        "HSdataset":dataset_receipt(hsdataset),
    }
    try:
        import Higgs.bounds as HB
        import Higgs.signals as HS
    except ImportError as exc:
        raise RuntimeError("HiggsTools Python interface unavailable; install v1.3") from exc
    hb=HB.Bounds(str(Path(hbdataset).resolve()))
    hs=HS.Signals(str(Path(hsdataset).resolve()))
    return {
        "status":"HIGGSTOOLS_EXECUTED_BUT_PHYSICS_INCOMPLETE",
        "backend":imports,
        "points":[run_one(o,hb,hs) for o in orders],
        "physical_higgs":"UNKNOWN",
        "globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hbdataset",type=Path,required=True)
    p.add_argument("--hsdataset",type=Path,required=True)
    p.add_argument("--order",choices=["NO","IO","both"],default="both")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    orders=("NO","IO") if args.order=="both" else (args.order,)
    result=calculate(orders,args.hbdataset,args.hsdataset)
    with args.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    for row in result["points"]:
        print(row["point_id"],"HS chi2 (effective-h125) =",
              row["branches"]["effective_coupling_h125_BR"]["chi2_HiggsSignals"],
              "HB diagnostic allowed =",
              row["branches"]["effective_coupling_h125_BR"]["HiggsBounds_allowed_with_incomplete_inputs"],
              "physical=UNKNOWN")


if __name__=="__main__":
    main()
