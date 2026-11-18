"""A strictly labeled, *non-transferred* 2HDM electroweak control for Zee points.

The doublet diagonal charged mass before Zee singlet mixing is
m_D² = cos²(gamma) m(H1±)² + sin²(gamma) m(H2±)².

Construct an auxiliary pure 2HDM with exactly the *same doublet scalar
quartics and neutral scalar spectrum*, but with no charged singlet,
F=0, gamma=0. Evaluate the pinned complete public 2HDM source there.
Its omitted electroweak terms diagnose risk at related parameters.
They are NOT the actual Zee gauge correction or an upper/lower bound.
"""
import argparse
import copy
import csv
import json
import math
from pathlib import Path

import numpy as np

from LFV.approximation_study import (
    AUTHOR_SHA256, assess_reference_channel, CHANNELS,
)
from LFV.constraints import input_parameters
from LFV.nonzero_mixing_study import DEFAULT_GRID, grid_points
from LFV.strict_2hdm_reference import evaluate as public_2hdm
from Model.model import scalar_parameters


def doublet_proxy(actual):
    p,F,D,ml,pole,fm,qy,widths,alpha,pot=input_parameters(actual)
    if pot is None:
        raise ValueError("Need reconstructed physical restricted Zee potential")
    if not abs(p["gamma"])>1e-9 or not np.max(abs(F))>0:
        raise ValueError("Must start from genuinely mixed Zee point with F nonzero")
    c,s=math.cos(p["gamma"]),math.sin(p["gamma"])
    doublet_mass=math.sqrt(c*c*p["mH1"]**2+s*s*p["mH2"]**2)
    # This source input is an auxiliary pure 2HDM point. The noninteracting
    # singlet is an artificial spectator with no physical impact.
    reference=copy.deepcopy(actual)
    reference["scalar"]["gamma"]=0.
    reference["scalar"]["mH1"]=doublet_mass
    reference["scalar"]["mH2"]=max(p["mH2"],doublet_mass+1)
    reference["F"]={"real":np.zeros((3,3)).tolist(),
                    "imag":np.zeros((3,3)).tolist()}
    for key in ("K1","K2","K12","Mu4"):
        reference["potential"][key]=0.
    # The public 2HDM code uses one charged-lepton mass prescription;
    # this declaration changes the *auxiliary control*, not the Zee input.
    reference["internal_lepton_masses_GeV"]=actual["pole_lepton_masses_GeV"][:]
    reference["mass_yukawa_scheme"]="AUTHOR_2HDM_CONTROL_WITH_POLE_LEPTONS_NOT_ZEE"
    reference["point_id"]=actual["point_id"]+"__doublet_reference"
    # If actual Zee masses and quartics were not mutually consistent, the
    # doublet proxy is meaningless. Reconstruct 2HDM-only quartics from
    # its neutral/charged spectrum at identical soft m12sq.
    beta=p["beta"]
    soft=(p["ma"]**2+pot["Lambda5"]*p["v"]**2)
    m12sq=soft*math.sin(beta)*math.cos(beta)
    test=scalar_parameters(
        **reference["scalar"],m12sq=m12sq,lc=.2,k1=0.,k2=0.,k12=0.)
    residual=max(abs(pot[k]-test[k])/max(1.,abs(pot[k]))
                 for k in ("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"))
    if not residual<1e-9:
        raise ValueError(f"Physical Zee input does not reproduce a consistent 2HDM doublet sector (residual {residual:g})")
    return reference, dict(
        doublet_unmixed_mass_GeV=doublet_mass,
        doublet_mass_tree_consistency_residual=float(residual),
        original_gamma_rad=p["gamma"],
        original_physical_H1_GeV=p["mH1"],
        original_physical_H2_GeV=p["mH2"],
        note="AUXILIARY_2HDM_DOUBLET_REFERENCE_NOT_A_ZEE_AMPLITUDE",
    )


def run_proxy(template,inversion,grid,sm,public_source):
    nodes=list(grid_points(template,inversion,grid))
    rows=[]
    for actual,coords in nodes:
        ref,identity=doublet_proxy(actual)
        paper=public_2hdm(ref,sm,public_source)
        if paper["Zee_full_dipoles"] is not None:
            raise AssertionError("Wrongly promoted 2HDM control to Zee LFV")
        for name in CHANNELS:
            data=assess_reference_channel(paper,name,sm["MW"])
            rows.append({
                "Zee_point_id":actual["point_id"],"channel":name,**coords,
                **identity,
                "pure_2HDM_control_omitted_over_photon":
                    data["author_omitted_over_approx"],
                "pure_2HDM_control_paper_variant_omitted_over_photon":
                    data["paper_HW_variant_omitted_over_approx"],
                "pure_2HDM_control_omitted_above_10pct":
                    data["omitted_over_approx_above_10pct"],
                "true_Zee_missing_gauge_amplitude_upper_bound":None,
                "true_Zee_LFV_gate":"UNKNOWN",
            })
    if len(rows)!=3*len(nodes):
        raise AssertionError("Missing a channel in doublet reference")
    return {
        "status":"AUXILIARY_2HDM_DOUBLET_CONTROL_NOT_TRANSFERRED",
        "n_nonzero_Zee_input_points":len(nodes),
        "n_2HDM_reference_evaluations":len(rows),
        "public_author_sha256":AUTHOR_SHA256,
        "missing_Zee_F_and_mixed_gauge_corrections":"UNBOUNDED",
        "full_Zee_LFV_gate":"UNKNOWN",
        "notes":[
            "Full gauge-invariant 2HDM source evaluated at unmixed doublet mass, not at physical mixed Zee charged masses.",
            "The exact neutral doublet potential is held fixed; charged singlet/F removed only inside auxiliary reference.",
            "The public author code uses its own quark mass running and pole lepton masses. No source quantity is added to a physical Zee amplitude.",
            "Ratios in this file are source-internal pure-2HDM diagnostics, neither estimates nor bounds on Zee mixed charged gauge diagrams.",
        ],
        "rows":rows,
        "summaries":{
            name:{
                "n_control_points":len([r for r in rows if r["channel"]==name]),
                "fraction_2HDM_omitted_above_10pct":float(np.mean([
                    r["pure_2HDM_control_omitted_above_10pct"]
                    for r in rows if r["channel"]==name])),
                "min_ratio":min(r["pure_2HDM_control_omitted_over_photon"]
                                for r in rows if r["channel"]==name),
                "median_ratio":float(np.median([
                    r["pure_2HDM_control_omitted_over_photon"]
                    for r in rows if r["channel"]==name])),
                "max_ratio":max(r["pure_2HDM_control_omitted_over_photon"]
                                for r in rows if r["channel"]==name),
                "true_Zee_LFV_gate":"UNKNOWN",
            } for name in CHANNELS
        },
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("template","inversion","sm","public-source","output-dir"):
        parser.add_argument("--"+arg,type=Path,required=True)
    parser.add_argument("--grid",type=Path)
    a=parser.parse_args()
    template=json.loads(a.template.read_text())
    inversion=json.loads(a.inversion.read_text())
    sm=json.loads(a.sm.read_text())
    grid=json.loads(a.grid.read_text()) if a.grid else DEFAULT_GRID
    result=run_proxy(template,inversion,grid,sm,a.public_source)
    a.output_dir.mkdir(parents=True,exist_ok=False)
    (a.output_dir/"doublet_proxy.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")
    with (a.output_dir/"doublet_proxy.csv").open("w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(result["rows"][0]))
        writer.writeheader()
        writer.writerows(result["rows"])
    print(json.dumps({"status":result["status"],
                      "n_points":result["n_nonzero_Zee_input_points"],
                      "summaries":result["summaries"],
                      "Zee_physical_gate":result["full_Zee_LFV_gate"]},indent=2))


if __name__=="__main__":
    main()
