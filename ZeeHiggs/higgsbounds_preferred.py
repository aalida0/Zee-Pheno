"""Actual official HiggsBounds crosscheck at preferred HXSWG h125 ggF.

Run both preferred and historic 13TeV h125 ggF policies through the
SAME unmodified HiggsTools v1.3 HiggsBounds v6.1 and HBDataSet v1.7.

Native h125 effective couplings (including both Zee charged loop
contributions) are preserved. h/H/A and both H+ inputs are sourced
from audited native SPheno but heavy/charged production remains
incomplete; do not interpret allowed=True as full-model allowed.
"""
import argparse
import json
from pathlib import Path

from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.diagnose_h125 import make_inputs, predictions_for
from ZeeHiggs.higgstools_pilot import (
    dataset_receipt, result_limits
)
from ZeeHiggs.charged_production_bridge import (
    load_charged_rates, insert_charged_rates
)
from ZeeHiggs.charged_effective_production import production_reference_official


def run_one(order, bounds, charged_rates_dir=None, charged_reference=False):
    import Higgs.tools.Input as hinput
    source=audit(order)
    data,ids,removed,old_rates=make_inputs(source,hinput)
    if charged_rates_dir is not None and charged_reference:
        raise ValueError("Select exactly one charged production input policy")
    charged_receipt=None
    charged_reference_result=None
    if charged_rates_dir is not None:
        charged_receipt=load_charged_rates(
            Path(charged_rates_dir)/f"s00017_{order}.json",source)
        data=insert_charged_rates(data,ids,charged_receipt)
    if charged_reference:
        charged_reference_result=production_reference_official(source)
        for pdg in (37,9000037):
            key=f"cxn_{ids[pdg]}_Hpmtb_LHC13"
            if key in data:
                raise ValueError("Source charged Hpmtb already exists; no override permitted")
            data[key]=charged_reference_result["higgstools_cross_sections"][
                f"cxn_{pdg}_Hpmtb_LHC13"]
    hkey=f"cxn_{ids[25]}_ggH_LHC13"
    if hkey not in data:
        raise ValueError("Expected SPheno historic ggF override")
    old=predictions_for(data, ids, hinput, explicit=False)
    modern_data=dict(data)
    modern_data.pop(hkey)
    modern=predictions_for(modern_data, ids, hinput, explicit=False)

    def check(pred,label):
        output=bounds(pred)
        selected=result_limits(output)
        if not selected:
            raise ValueError("No official HiggsBounds selected analyses: "+label)
        return {
            "status":"OFFICIAL_HIGGSBOUNDS_SELECTED_LIMITS_DIAGNOSTIC_ONLY",
            "raw_HiggsBounds_allowed_for_incomplete_inputs":bool(output.allowed),
            "selectedLimits":selected,
            "largest_observed_ratio":max(v["obs_ratio"] for v in selected),
            "selected_neutral_particle_ids":sorted(set(
                v["particle"] for v in selected if v["particle"] in ("25","35","36"))),
            "selected_charged_particle_ids":sorted(set(
                v["particle"] for v in selected if v["particle"] in ("37","9000037"))),
            "is_full_Zee_model_Higgs_exclusion":False,
        }

    h=modern.particle(ids[25])
    modern_gg=float(h.cxn("LHC13","ggH"))
    old_gg=float(old.particle(ids[25]).cxn("LHC13","ggH"))
    if abs(old_gg-old_rates[hkey])>1e-7:
        raise ValueError("Old SPheno ggF input not reproduced")
    if abs(modern_gg-old_gg)<1:
        raise ValueError("Preferred HXSWG ggF change not reproduced")
    return {
        "point_id":source["point_id"],
        "status":"OFFICIAL_HIGGSBOUNDS_HXSWG_VS_SPHENO_RUN_COMPLETE",
        "source":source["source"],
        "old_spheno_ggF_13TeV_pb":old_gg,
        "preferred_HXSWG_ggF_13TeV_pb":modern_gg,
        "loop_coupling_policy":{"calcggH":False,"calcHgamgam":False},
        "old_spheno_ggF":check(old,"historic SPheno"),
        "preferred_HXSWG_ggF":check(modern,"preferred HXSWG"),
        "charged_production_complete":False,
        "charged_Hpmtb_13TeV_receipt":charged_receipt,
        "charged_Hpmtb_13TeV_reference":charged_reference_result,
        "charged_Hpmtb_13TeV_input_present":(
            charged_receipt is not None or charged_reference_result is not None),
        "charged_production_missing_or_unvalidated":(
            ["13 TeV charged Hpmtb missing"] if not (
                charged_receipt or charged_reference_result) else
            ["Generator raw rates/amplitudes not independently revalidated"]
            if charged_receipt is not None else
            ["Official HiggsTools tabulated reference, model matching requires validation"]
        )+["8 TeV and other process modes not source-matched"],
        "heavy_neutral_production_source_matched":False,
        "neutral_selected_limits_physical_acceptance":"UNKNOWN",
        "charged_selected_limits_physical_acceptance":"UNKNOWN",
        "physical_higgs":"UNKNOWN",
        "globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hbdataset",type=Path,required=True)
    p.add_argument("--charged-rates-dir",type=Path,
                   help="Optional exact-point s00017_NO/IO 13 TeV 5FS MG5 receipt JSON directory")
    p.add_argument("--charged-reference",action="store_true",
                   help="Use official HiggsTools 13 TeV tabulated tbH cross section for Zee Type-II doublet fractions")
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    import Higgs.bounds as HB
    hb=HB.Bounds(str(a.hbdataset.resolve()))
    out={
        "status":"OFFICIAL_HIGGSBOUNDS_PREFERRED_HXSWG_NATIVE_TWO_CHARGED_SCALARS",
        "backend":"HiggsTools v1.3 / HiggsBounds v6.1",
        "dataset":dataset_receipt(a.hbdataset),
        "points":[run_one(o,hb,a.charged_rates_dir,a.charged_reference) for o in ("NO","IO")],
        "physical_higgs":"UNKNOWN",
        "HiggsBounds_full_model_exclusion":None,
        "globally_allowed":None,
    }
    with a.output.open("x") as f:
        json.dump(out,f,indent=2,allow_nan=False)
        f.write("\n")
    for row in out["points"]:
        h=row["preferred_HXSWG_ggF"]
        print(row["point_id"],
              "HB_preferred_raw_allowed=",h["raw_HiggsBounds_allowed_for_incomplete_inputs"],
              "maxObsRatio=",h["largest_observed_ratio"],
              "chargedSelected=",h["selected_charged_particle_ids"],
              "physical=UNKNOWN")


if __name__=="__main__":
    main()
