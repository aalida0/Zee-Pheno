"""HiggsSignals/SPheno h125 BR and correlated-measurement audit.

HiggsTools is the ONLY Higgs experimental backend. We compare:
(A) all explicit SPheno branching ratios, (B) effective-coupling h125
branching ratios, against an official same-mass SM reference.

Individual Measurement.__call__ evaluations KEEP intra-measurement
correlations. Their sum is only a diagnostic: it may not reproduce the
global HiggsSignals chi2 if inter-measurement correlations are present.

A channel "splice" replaces ONE partial width in the effective-branch
particle with SPheno's Gamma_total * BR for that channel and then
recomputes the HiggsSignals chi2. Other BRs change via normalization:
these are finite-difference sensitivity probes, NOT additive chi2 terms.
"""
import argparse
import json
import math
from pathlib import Path

from ZeeHiggs.spheno_audit import audit, read_slha
from ZeeHiggs.higgstools_pilot import (
    NEUTRAL, CHARGED, infer_ids, clean_missing_zero_couplings,
    dataset_receipt, sm_higgs_reference,
)

CHANNELS = ("bb", "cc", "tautau", "mumu", "ss", "WW", "ZZ", "gg", "gamgam")
MAX_CHI2 = 1e12


def number(value, where):
    v = float(value)
    if not math.isfinite(v) or abs(v) > MAX_CHI2:
        raise ValueError(f"Nonfinite or excessive diagnostic {where}: {value}")
    return v


def make_inputs(report, hinput):
    data = hinput.readHB5SLHA(
        report["spectrum_path"], neutralPDGs=list(NEUTRAL),
        chargedPDGs=list(CHARGED))
    if not isinstance(data, dict):
        raise ValueError("Official readHB5SLHA did not return a dict")
    ids = infer_ids(data, report)
    data, removed = clean_missing_zero_couplings(data)
    production_overrides = {}
    for pdg in NEUTRAL:
        key = f"cxn_{ids[pdg]}_ggH_LHC13"
        production_overrides[key] = report["neutral_ggF_LHC13_pb"][str(pdg)]
        data[key] = production_overrides[key]
    return data, ids, removed, production_overrides


def predictions_for(data, ids, hinput, *, explicit):
    neutral = [ids[p] for p in NEUTRAL]
    charged = [ids[p] for p in CHARGED]
    if explicit:
        option = True
    else:
        option = {pid: (pid != ids[25]) for pid in neutral + charged}
    return hinput.predictionsFromDict(
        data, neutral, charged, [], useExplicitBr=option)


def particle_report(predictions, identifier):
    h = predictions.particle(identifier)
    br = {decay: number(h.br(decay), f"br_{identifier}_{decay}")
          for decay in CHANNELS}
    w = number(h.totalWidth(), f"width_{identifier}")
    if w <= 0 or any(b < 0 or b > 1.000001 for b in br.values()):
        raise ValueError("Invalid Higgs width or branching fraction")
    return {
        "width_GeV": w, "branching_ratios": br,
        "partial_widths_GeV": {c: b*w for c, b in br.items()},
        "listed_channel_sum": sum(br.values()),
        "note": ("Listed-channel sum excludes channels without an input, "
                 "including missing h->Zgamma. This is NOT a full-width "
                 "closure check.")
    }


def measure_all(signals, prediction, *, label):
    ref = signals.measurements()
    rows = []
    seen = set()
    for m in ref:
        m_id = str(m.id())
        if m_id in seen:
            raise ValueError("Duplicate measurement id " + m_id)
        seen.add(m_id)
        val = number(m(prediction), f"{label}:measurement:{m_id}")
        rows.append({
            "id": m_id,
            "chi2": val,
            "reference": str(m.reference()),
            "citeKey": str(m.citeKey()),
            "source": str(m.loadedFrom()),
            "n_subobservables": int(m.nSubMeasurements()),
            "included_in_measured_total": True,
        })
    return rows


def compare_measurements(a, b, c):
    indexed = [{str(r["id"]):r for r in x} for x in (a,b,c)]
    if len(indexed[0]) != len(a) or not set(indexed[0]) == set(indexed[1]) == set(indexed[2]):
        raise ValueError("Inconsistent per-measurement datasets")
    out = []
    for id_, row in indexed[0].items():
        row_b, row_c = indexed[1][id_], indexed[2][id_]
        out.append({
            "id": id_, "reference":row["reference"], "source":row["source"],
            "n_subobservables":row["n_subobservables"],
            "sm_chi2":row["chi2"],
            "effective_chi2":row_b["chi2"],
            "explicit_chi2":row_c["chi2"],
            "delta_effective_vs_SM":row_b["chi2"]-row["chi2"],
            "delta_explicit_vs_effective":row_c["chi2"]-row_b["chi2"],
            "delta_explicit_vs_SM":row_c["chi2"]-row["chi2"],
            "has_correlations_within_measurement": "included_by_official_measurement",
        })
    out.sort(key=lambda r:abs(r["delta_explicit_vs_effective"]), reverse=True)
    return out


def diagnose(order, signals):
    import Higgs.predictions as HP
    import Higgs.tools.Input as hinput

    source = audit(order)
    data, ids, missing, production = make_inputs(source, hinput)
    exp = predictions_for(data, ids, hinput, explicit=True)
    eff = predictions_for(data, ids, hinput, explicit=False)

    model_id = ids[25]
    mh = source["masses_GeV"]["25"]
    sm = HP.Predictions()
    h_sm = sm.addParticle(HP.BsmParticle("SM_h", charge="neutral", cp="even"))
    h_sm.setMass(mh)
    HP.effectiveCouplingInput(h_sm, HP.smLikeEffCouplings)
    chis = {"SM":number(signals(sm),"SM"),
            "effective":number(signals(eff),"effective"),
            "explicit":number(signals(exp),"explicit")}
    m_sm=measure_all(signals,sm,label="SM")
    m_eff=measure_all(signals,eff,label="effective")
    m_exp=measure_all(signals,exp,label="explicit")
    rows=compare_measurements(m_sm,m_eff,m_exp)

    exp_p=particle_report(exp,model_id)
    eff_p=particle_report(eff,model_id)
    sm_p=particle_report(sm,"SM_h")
    sp_width = source["total_widths_GeV"]["25"]
    slha_blocks, slha_decays = read_slha(Path(source["spectrum_path"]).read_text())
    daughters={
        (5,-5):"bb",(-5,5):"bb",
        (4,-4):"cc",(-4,4):"cc",
        (15,-15):"tautau",(-15,15):"tautau",
        (13,-13):"mumu",(-13,13):"mumu",
        (3,-3):"ss",(-3,3):"ss",
        (24,-24):"WW",(-24,24):"WW",
        (23,23):"ZZ", (21,21):"gg", (22,22):"gamgam",
    }
    native={k:0. for k in CHANNELS}
    for channel in slha_decays[25]["channels"]:
        if channel["pdgs"] in daughters:
            native[daughters[channel["pdgs"]]] += channel["br"]
    # The actual native BRs must match what HiggsTools imports to 1e-8;
    # check direct translation before interpreting any numerical difference.
    for mode,v in native.items():
        if v == 0:
            continue
        if abs(v-exp_p["branching_ratios"][mode])>2e-7:
            raise ValueError(f"Native SLHA to HiggsTools explicit BR mismatch {mode}")
    sp_native={
        "total_width_GeV": sp_width,
        "BR_by_channel": native,
        "partial_width_by_channel_GeV":{k:v*sp_width for k,v in native.items()}
    }

    finite_diff=[]
    for ch in CHANNELS:
        if native[ch] <= 0:
            continue
        variation=predictions_for(data, ids, hinput, explicit=False)
        p=variation.particle(model_id)
        old_g=p.totalWidth()*p.br(ch)
        actual=sp_width*native[ch]
        p.setDecayWidth(ch, actual)
        chi=number(signals(variation),"single_replaced_"+ch)
        finite_diff.append({
            "replaced_partial_width_channel":ch,
            "source_SPheno_Gamma_GeV":actual,
            "baseline_effective_Gamma_GeV":old_g,
            "ratio_SPheno_to_effective_Gamma":actual/old_g if old_g>0 else None,
            "modified_total_width_GeV":number(p.totalWidth(),"modified width"),
            "chi2_after_one_channel_replacement":chi,
            "delta_chi2_vs_effective":chi-chis["effective"],
            "note":"One-at-a-time partial-width splice changes all normalized BRs; not additive.",
        })
    finite_diff.sort(key=lambda r:abs(r["delta_chi2_vs_effective"]),reverse=True)

    # Isolate impact of explicit SPheno LHC13 ggH override, holding
    # all other effective-coupling and BR inputs fixed.
    without_overrides=dict(data)
    for key in production:
        without_overrides.pop(key)
    eff_no_gg_override=predictions_for(without_overrides,ids,hinput,explicit=False)
    gg_override_chi2=number(signals(eff_no_gg_override),"no_gg_override")

    # Compare *identical* process in same predictions object, against
    # the official SM-like reference, to detect normalization drift.
    native_ggf = float(data[f"cxn_{model_id}_ggH_LHC13"])
    eff_auto_ggf = number(eff_no_gg_override.particle(model_id).cxn("LHC13","ggH"),"auto ggH")
    sm_ggf = number(sm.particle("SM_h").cxn("LHC13","ggH"),"SM ggH")

    return {
        "point_id":source["point_id"],
        "status":"OFFICIAL_HIGGSTOOLS_WIDTH_AND_MEASUREMENT_DIAGNOSTIC",
        "physical_higgs":"UNKNOWN","higgssignals_exclusion":None,
        "source":source["source"], "h125_mass_GeV":mh,
        "higgsSignals_chi2":chis,
        "higgsSignals_delta": {
            "effective_vs_SM":chis["effective"]-chis["SM"],
            "explicit_vs_SM":chis["explicit"]-chis["SM"],
            "explicit_vs_effective":chis["explicit"]-chis["effective"],
        },
        "measurements_count":len(rows),
        "per_measurement_correlated_chi2":rows,
        "per_measurement_sum_diagnostic":{
            "SM":sum(r["chi2"] for r in m_sm),
            "effective":sum(r["chi2"] for r in m_eff),
            "explicit":sum(r["chi2"] for r in m_exp),
            "warning":("Per-measurement totals can differ from global Signals chi2 "
                       "because the global dataset may include cross-measurement correlations.")
        },
        "SPheno_raw":sp_native,
        "HiggsTools_explicit":exp_p,
        "HiggsTools_effective":eff_p,
        "HiggsTools_SM":sm_p,
        "partial_width_one_at_a_time_diagnostics":finite_diff,
        "production_ggH_13TeV_pb": {
            "SPheno_explicit_override":native_ggf,
            "HiggsTools_effective_prediction_without_override":eff_auto_ggf,
            "HiggsTools_SM_reference":sm_ggf,
            "chi2_without_ggH_override":gg_override_chi2,
            "delta_chi2_vs_baseline_effective":gg_override_chi2-chis["effective"],
        },
        "missing":{
            "SPheno_Zgamma":"NOT_CALCULATED, no physical zero inferred",
            "imported_uncalculated_key_removals":missing,
            "source_charged_production":"NOT_VALIDATED",
            "experimental_data_2026_coverage":"NOT_COMPLETE",
            "full_model_HiggsSignals_95CL":"UNKNOWN",
        },
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hsdataset",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    import Higgs.signals as HS
    data=dataset_receipt(args.hsdataset)
    signals=HS.Signals(str(args.hsdataset.resolve()))
    res={"status":"EXECUTED_OFFICIAL_HIGGSTOOLS_H125_DIAGNOSTIC",
         "backend":"HiggsTools v1.3, HiggsSignals v3.1",
         "hsdataset":data, "observable_count":int(signals.observableCount()),
         "points":[diagnose(order,signals) for order in ("NO","IO")],
         "physical_higgs":"UNKNOWN","globally_allowed":None}
    with args.output.open("x") as f:
        json.dump(res,f,indent=2,allow_nan=False)
        f.write("\n")
    for row in res["points"]:
        print(row["point_id"],"chi2(exp-eff)=",row["higgsSignals_delta"]["explicit_vs_effective"],
              "width(exp/eff)=",row["HiggsTools_explicit"]["width_GeV"]/row["HiggsTools_effective"]["width_GeV"],
              "highest measurement delta=",row["per_measurement_correlated_chi2"][0]["delta_explicit_vs_effective"],
              "full-model=UNKNOWN")


if __name__=="__main__":
    main()
