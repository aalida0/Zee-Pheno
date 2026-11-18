"""Reproducible *strict 2HDM* control sweep of omitted LFV amplitudes.

Compute the missing non-photonic / gauge terms by subtracting the same
author's 2017-like approximation from their full (erratum-aware) 2HDM sum.
This is NOT a bound on the general Zee F-dependent two-loop remainder.
Never change point parameters in response to any LFV limit.
"""
import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

from LFV.strict_2hdm_reference import (
    decoupled_input, evaluate as public_reference, load_public, SOURCE_BLOCKS,
)
from Model.model import scalar_parameters

# Only the subsets we implement independently for the Zee model.
# tHp/cHp/tauHp, Z, kinetic and W associated terms are *not* included.
RETAINED_SOURCE_BLOCKS = frozenset({
    "1loop", "th", "ch", "bh", "tauh", "HA"
})
OMITTED_2HDM_SOURCE_BLOCKS = frozenset(SOURCE_BLOCKS) - RETAINED_SOURCE_BLOCKS
AUTHOR_SHA256 = "9b4af7aa4bf7bf533a7b82c8c2603d118346424e3b7b5f9f06ab8cb41bc4e193"
DEFAULT_GRID = {
    "tan_beta": [1.5, 5., 12.],
    "beta_minus_alpha_rad": [-.15, 0., .15],
    "heavy_neutral_mass_GeV": [350., 700.],
}
CHANNELS = ("mu_to_e_gamma", "tau_to_e_gamma", "tau_to_mu_gamma")


def complex_of(value):
    return complex(float(value["real"]), float(value["imag"]))


def norm(coeff):
    return float(np.linalg.norm(np.asarray(coeff, complex)))


def assess_reference_channel(reference, channel, author_mw):
    dipoles = reference["reference_dipoles"][channel]
    ledger = dipoles.get("raw_public_source_components")
    if not isinstance(ledger, dict) or set(ledger) != set(SOURCE_BLOCKS):
        raise ValueError("Full source-chirality component coverage required")
    def vec(key):
        return np.array([
            complex_of(ledger[key]["C_L"]), complex_of(ledger[key]["C_R"])
        ], complex)
    approx = sum((vec(k) for k in RETAINED_SOURCE_BLOCKS), np.zeros(2,complex))
    omitted = sum((vec(k) for k in OMITTED_2HDM_SOURCE_BLOCKS),
                  np.zeros(2,complex))
    full = approx + omitted
    raw = np.array([complex_of(dipoles["C_L_GeV_minus2"]),
                    complex_of(dipoles["C_R_GeV_minus2"])])
    if not np.allclose(full,raw,rtol=2e-10,atol=1e-24):
        raise ArithmeticError("Retained + omitted source terms fail to reproduce full reference")
    # Paper-printed alternative for source HW extra M_W^-2 discrepancy:
    # keep both results, never silently "repair" the author source.
    group_norms={key:norm(vec(key)) for key in sorted(OMITTED_2HDM_SOURCE_BLOCKS)}
    largest=max(group_norms,key=group_norms.get)
    hw_adjustment=(author_mw**2-1.)*vec("HW")
    paper_full=full+hw_adjustment
    paper_omitted=omitted+hw_adjustment
    floor=1e-30
    def ratio(x,y):
        return norm(x)/max(norm(y),floor)
    return {
        "author_approx_norm_GeV_minus2":norm(approx),
        "author_omitted_norm_GeV_minus2":norm(omitted),
        "largest_omitted_source_group":largest,
        "largest_omitted_source_group_norm_GeV_minus2":group_norms[largest],
        "author_full_norm_GeV_minus2":norm(full),
        "author_omitted_over_approx":ratio(omitted,approx),
        "author_omitted_over_full":ratio(omitted,full),
        "author_full_BR_over_approx_BR":ratio(full,approx)**2,
        "paper_HW_variant_omitted_norm_GeV_minus2":norm(paper_omitted),
        "paper_HW_variant_full_norm_GeV_minus2":norm(paper_full),
        "paper_HW_variant_omitted_over_approx":ratio(paper_omitted,approx),
        "paper_HW_variant_full_BR_over_approx_BR":ratio(paper_full,approx)**2,
        "approx_zero":bool(norm(approx)<floor),
        "omitted_over_approx_above_10pct":
            bool(norm(omitted)>0.1*norm(approx)),
        "paper_variant_omitted_over_approx_above_10pct":
            bool(norm(paper_omitted)>0.1*norm(approx)),
        "genuine_Zee_F_gauge_remainder":"UNBOUNDED",
        "physical_Zee_LFV_gate":"UNKNOWN",
        "comparison_scope":"STRICT_DECOUPLED_2HDM_CONTROL_NOT_FULL_ZEE",
    }


def grid_points(template, grid=None):
    """Fixed-factor grid; no rejection/adaptation to favor LFV-compatible cards."""
    grid=DEFAULT_GRID if grid is None else grid
    if set(grid)!=set(DEFAULT_GRID):
        raise ValueError("Expected explicit three-factor 2HDM control grid")
    for key in DEFAULT_GRID:
        values=grid[key]
        if not isinstance(values,list) or not values or len(values)>15 or any(
            not isinstance(x,(int,float)) or not math.isfinite(x) for x in values
        ):
            raise ValueError(f"Invalid factor grid for {key}")
    for tanb,delta,mH in itertools.product(
        grid["tan_beta"],grid["beta_minus_alpha_rad"],
        grid["heavy_neutral_mass_GeV"]
    ):
        if not (0.4<tanb<25 and abs(delta)<.4 and mH>150):
            raise ValueError("Control factors outside documented grid domain")
        point=copy.deepcopy(template)
        beta=math.atan(tanb)
        scalar=point["scalar"]
        scalar.update(beta=beta,alpha=beta-delta,gamma=0.,
                      mh1=125.,mh2=float(mH),ma=float(mH+70),
                      mH1=float(mH+35),mH2=float(mH+535))
        source=scalar_parameters(**scalar,
            m12sq=(mH**2)*math.sin(beta)*math.cos(beta),
            lc=.2,k1=0.,k2=0.,k12=0.)
        point["potential"]={k:source[k] for k in (
            "Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
            "K1","K2","K12","Mu4"
        )}
        point["point_id"]=f"2hdm_t{tanb:g}_d{delta:+.2f}_M{int(mH)}"
        yield point, {
            "tan_beta":float(tanb),
            "beta_minus_alpha_rad":float(delta),
            "heavy_neutral_mass_GeV":float(mH),
        }


def sweep(template, sm, public_source, grid=None):
    # Fail early if stale or modified upstream source changes conventions.
    actual=hashlib.sha256(public_source.read_bytes()).hexdigest()
    if actual != AUTHOR_SHA256:
        raise ValueError("Public author thdm.py SHA-256 differs from independently reviewed source")
    rows=[]
    for point,coordinates in grid_points(template,grid):
        result=public_reference(point,sm,public_source)
        if result["Zee_full_dipoles"] is not None:
            raise ValueError("Reference wrongly promoted to a Zee prediction")
        for channel in CHANNELS:
            materiality=assess_reference_channel(result,channel,sm["MW"])
            rows.append(dict(point_id=point["point_id"],
                             channel=channel,**coordinates,**materiality))
    output={
        "status":"STRICT_2HDM_CONTROL_MATERIALITY_ONLY",
        "author_thdm_sha256":actual,
        "grid":DEFAULT_GRID if grid is None else grid,
        "n_points":len(rows)//len(CHANNELS),
        "n_channel_evaluations":len(rows),
        "retained_source_blocks":sorted(RETAINED_SOURCE_BLOCKS),
        "omitted_source_blocks":sorted(OMITTED_2HDM_SOURCE_BLOCKS),
        "true_Zee_full_LFV_gate":"UNKNOWN",
        "notes":[
            "No light neutrino fit, scalar viability, vacuum or collider gates applied.",
            "Only a strict F=gamma=Mu4=K1=K2=K12=0 2HDM control, not Zee LFV.",
            "The 2017-like approximation contains only author 1loop, neutral-fermion-photon and charged-scalar-photon blocks.",
            "Both raw public-HW and paper-printed-HW alternatives are reported; their difference is a source ambiguity, not a BSM error bar.",
            "An omitted/retained ratio is not a universal bound outside this exact decoupling limit.",
        ],
        "rows":rows,
    }
    output["summaries"]={}
    for channel in CHANNELS:
        z=[r for r in rows if r["channel"]==channel and not r["approx_zero"]]
        output["summaries"][channel]=dict(
            n_nonzero_control_points=len(z),
            fraction_raw_over_10pct=(
                sum(r["omitted_over_approx_above_10pct"] for r in z)/len(z)
                if z else None),
            fraction_paper_variant_over_10pct=(
                sum(r["paper_variant_omitted_over_approx_above_10pct"] for r in z)/len(z)
                if z else None),
            raw_ratio_min=min((r["author_omitted_over_approx"] for r in z),default=None),
            raw_ratio_median=float(np.median([
                r["author_omitted_over_approx"] for r in z])) if z else None,
            raw_ratio_max=max((r["author_omitted_over_approx"] for r in z),default=None),
            paper_variant_ratio_max=max((
                r["paper_HW_variant_omitted_over_approx"] for r in z),default=None),
            dominant_source_counts={
                group:sum(r["largest_omitted_source_group"]==group for r in z)
                for group in sorted(OMITTED_2HDM_SOURCE_BLOCKS)
            },
            physical_Zee_lfv_status="UNKNOWN",
        )
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("template","sm","public-source","output-dir"):
        ap.add_argument("--"+key,required=True,type=Path)
    ap.add_argument("--grid",type=Path,help="Optional explicit three-factor grid")
    a=ap.parse_args()
    data=sweep(json.loads(a.template.read_text()),json.loads(a.sm.read_text()),
               a.public_source,
               json.loads(a.grid.read_text()) if a.grid else None)
    a.output_dir.mkdir(parents=True,exist_ok=False)
    (a.output_dir/"materiality.json").write_text(
        json.dumps(data,indent=2,allow_nan=False)+"\n")
    with (a.output_dir/"materiality.csv").open("w",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    print(json.dumps({"status":data["status"],
                      "author_sha256":data["author_thdm_sha256"],
                      "n_points":data["n_points"],
                      "summary":data["summaries"]},indent=2))


if __name__=="__main__":
    main()
