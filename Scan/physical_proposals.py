"""Reproducible scalar + full Zee neutrino-family proposal scan.

Generates *analytical proposals*, not accepted points. Each scalar is a
separate sampling unit; multiple complex F/D families share its scalar input.
Every saved SLHA card needs a fresh native SPheno run and source-bound gates.
"""
import argparse
import hashlib
import json
import math
import random
from pathlib import Path

import numpy as np

from Model.audit_neutrino_native import sha
from Model.check_native import blocks, fix_card
from Model.model import scalar_parameters
from Neutrino.fit import (REFERENCE_CARD_SHA, complex_matrix, encode,
                          inverse, target)

MZ=91.1876
SCALAR_KEYS=("M112","M222","M12Sq","MCS2","Mu4",
             "Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
             "LambdaC","K1","K2","K12")
RANGES=("tan_beta","gamma_rad","mH2_GeV","mH1_fraction",
        "mH_offset_GeV","mA_offset_GeV","alignment_delta_rad",
        "soft_mass_ratio","lambda_c","K1","K2","K12")
FAMILY_RANGES=("lightest_eV","Majorana_phase_rad","log10_F_scale",
               "log10_homogeneous_scale","null_t")
MAX_SCALARS=5000


def validate(config):
    if (not isinstance(config,dict) or
        set(config)!={"schema_version","seed","scalar_count",
                      "families_per_scalar","arms","neutrino"} or
        config["schema_version"]!=2 or
        type(config["seed"]) is not int or
        type(config["scalar_count"]) is not int or
        not 1<=config["scalar_count"]<=MAX_SCALARS or
        type(config["families_per_scalar"]) is not int or
        not 1<=config["families_per_scalar"]<=50):
        raise ValueError("Bad versioned scalar/neutrino scan configuration")
    arms=config["arms"]
    if not isinstance(arms,list) or not arms or len(arms)>10:
        raise ValueError("Require 1..10 named sampling arms")
    names=[]
    for arm in arms:
        if not isinstance(arm,dict) or set(arm)!={"name","weight","ranges"}:
            raise ValueError("Bad sampling arm")
        name=arm["name"]
        if not isinstance(name,str) or not name.isidentifier() or name in names:
            raise ValueError("Sampling-arm names must be unique identifiers")
        names.append(name)
        if type(arm["weight"]) is not int or arm["weight"]<=0:
            raise ValueError("Sampling arm weights must be positive integers")
        check_ranges(arm["ranges"],RANGES)
    nu=config["neutrino"]
    if not isinstance(nu,dict) or set(nu)!={"orderings","ranges"} or (
        not isinstance(nu["orderings"],list) or not nu["orderings"] or
        set(nu["orderings"])-{"NO","IO"} or len(set(nu["orderings"]))!=len(nu["orderings"])
    ):
        raise ValueError("Require explicitly named unique NO/IO targets")
    check_ranges(nu["ranges"],FAMILY_RANGES)
    return config


def check_ranges(ranges,required):
    if not isinstance(ranges,dict) or set(ranges)!=set(required):
        raise ValueError("Missing or unknown scan prior ranges")
    for key,pair in ranges.items():
        if (not isinstance(pair,list) or len(pair)!=2 or
            any(type(x) not in (float,int) or not math.isfinite(x) for x in pair) or
            pair[0]>pair[1]):
            raise ValueError(f"Invalid prior interval: {key}")


def cheap_scalar(p):
    """Necessary doublet BFB, sufficient full-Zee quartic BFB, |coupling| proxy."""
    l1,l2,l3,l4,l5=(p[k] for k in
                    ("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"))
    lc,k1,k2,k12=(p[k] for k in ("LambdaC","K1","K2","K12"))
    necessary=(l1>0 and l2>0 and
               l3+math.sqrt(l1*l2)>0 and
               l3+l4-abs(l5)+math.sqrt(l1*l2)>0)
    # Conservative sufficient condition via 2sqrt(xy)<=x+y and
    # |Phi1†Phi2|<=sqrt(xy); not a global-minimum certificate.
    sufficient=(necessary and lc>0 and
                k1>=abs(k12) and k2>=abs(k12))
    maxquartic=max(abs(p[k]) for k in SCALAR_KEYS if k not in
                   ("M112","M222","M12Sq","MCS2","Mu4"))
    return dict(
        doublet_BFB_necessary="PASS" if necessary else "FAIL",
        zee_quartic_BFB_sufficient="PASS" if sufficient else "UNKNOWN",
        quartic_lt_4pi_convention="PASS" if maxquartic<4*math.pi else "FAIL",
        max_abs_quartic=maxquartic,
        full_tree_unitarity="UNKNOWN",
        global_vacuum="UNKNOWN",
        EWPO="UNKNOWN",
        Higgs_rates="UNKNOWN",
        flavor="UNKNOWN",
    )


def scalar_point(arm,rng,v):
    r=arm["ranges"]
    sample=lambda key:rng.uniform(*r[key])
    tanb=sample("tan_beta")
    gamma=sample("gamma_rad")
    beta=math.atan(tanb)
    alpha=beta-sample("alignment_delta_rad")
    heavy=sample("mH2_GeV")
    light=heavy*sample("mH1_fraction")
    md=math.hypot(light*math.cos(gamma),heavy*math.sin(gamma))
    mH=md+sample("mH_offset_GeV")
    ma=md+sample("mA_offset_GeV")
    soft=md*sample("soft_mass_ratio")
    scalar=dict(v=v,beta=beta,alpha=alpha,gamma=gamma,mh1=125.25,
                mh2=mH,ma=ma,mH1=light,mH2=heavy)
    if not (heavy>light+MZ and mH>=125.25 and ma>0):
        return scalar,None,dict(status="FAIL_PREDECLARED_KINEMATIC_OR_MASS_PRIOR")
    try:
        lag=scalar_parameters(**scalar,
            m12sq=soft**2*math.sin(beta)*math.cos(beta),
            lc=sample("lambda_c"),k1=sample("K1"),
            k2=sample("K2"),k12=sample("K12"))
    except ValueError as exc:
        return scalar,None,dict(status="INVALID_SCALAR_INPUT",reason=str(exc))
    flags=cheap_scalar(lag)
    # Keep all points in the denominator; write expensive native cards
    # only if basic necessary and analysis-defined perturbative priors pass.
    flags["status"]=("ANALYTIC_SCALAR_PRESELECTION"
        if flags["doublet_BFB_necessary"]=="PASS" and
           flags["quartic_lt_4pi_convention"]=="PASS"
        else "SCREENED_AT_DECLARED_SCALAR_PRIORS")
    return scalar,lag,flags


def scalar_updates(p,tanbeta):
    from Model.check_native import SCALAR_BLOCKS
    assert SCALAR_BLOCKS==("ZEEPOT","ZEEQUARTIC")
    names=("M112","M222","M12Sq","MCS2","Mu4")
    quart=("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
           "LambdaC","K1","K2","K12")
    return {**{("ZEEPOT",(i,)):p[name] for i,name in enumerate(names,1)},
            **{("ZEEQUARTIC",(i,)):p[name] for i,name in enumerate(quart,1)},
            ("MINPAR",(3,)):tanbeta}


def prepare_card(template,scalar,lag,F,Ye2):
    updates=scalar_updates(lag,math.tan(scalar["beta"]))
    for name,mat in (("FZEEIN",F),("YEE2IN",Ye2)):
        for i in range(3):
            for j in range(3):
                updates[(name,(i+1,j+1))]=float(mat[i,j].real)
                updates[("IM"+name,(i+1,j+1))]=float(mat[i,j].imag)
    return fix_card(template,updates)


def generate(config_path,template_path,snapshot_path,output_dir):
    config=validate(json.loads(config_path.read_text()))
    if sha(template_path)!=REFERENCE_CARD_SHA:
        raise ValueError("Only the source-attested loop-enabled native card is accepted")
    original=template_path.read_text()
    data=blocks(original)
    required={
        "MINPAR":{(3,)},
        "ZEEPOT":{(i,) for i in range(1,6)},
        "ZEEQUARTIC":{(i,) for i in range(1,10)},
    }
    if any(not indices<=set(data.get(name,{}))
           for name,indices in required.items()):
        raise ValueError("Attested template lacks required scalar/tan-beta entries")
    snapjson=json.loads(snapshot_path.read_text())
    snap=(snapjson["snapshots"][snapjson["spectrum_snapshot_index"]]
          if "snapshots" in snapjson else snapjson["snapshot"])
    mass=np.asarray(snap["MFe"],float)
    if mass.shape!=(3,) or np.any(mass<=0):
        raise ValueError("Matched native lepton masses required")
    v=math.hypot(float(snap["vd"]),float(snap["vu"]))
    oldbeta=math.atan2(snap["vu"],snap["vd"])
    ye0=complex_matrix(snap["Ye"])
    ye20=complex_matrix(snap["Ye2"])
    sy=math.cos(oldbeta)*ye0+math.sin(oldbeta)*ye20
    el,er=(complex_matrix(snap[name]) for name in ("ZEL","ZER"))
    arms=[arm for arm in config["arms"] for _ in range(arm["weight"])]
    # All scalar attempts are independent; Yukawa families are conditional.
    rows=[]
    cards=[]
    for i in range(config["scalar_count"]):
        rng=random.Random(f"{config['seed']}:{i}")
        arm=arms[i%len(arms)]
        scalar,lag,gates=scalar_point(arm,rng,v)
        record=dict(scalar_id=f"s{i:05d}",arm=arm["name"],
                    scalar=scalar,lagrangian=lag,gates=gates,
                    families=[])
        rows.append(record)
        if gates["status"]!="ANALYTIC_SCALAR_PRESELECTION":
            continue
        charged=(scalar["mH1"],scalar["mH2"])
        sine=math.sin(2*scalar["gamma"])
        nu=config["neutrino"]
        r=nu["ranges"]
        def draw(key):return rng.uniform(*r[key])
        for k in range(config["families_per_scalar"]):
            ordering=nu["orderings"][k%len(nu["orderings"])]
            mlight=draw("lightest_eV")
            a21,a31=draw("Majorana_phase_rad"),draw("Majorana_phase_rad")
            root=(k//len(nu["orderings"]))%2
            chart=((0,1,2),(1,2,0),(2,0,1))[k%3]
            t=complex(draw("null_t"),draw("null_t"))
            fscale=10**draw("log10_F_scale")
            hscale=10**draw("log10_homogeneous_scale")
            z=np.zeros(4,complex) if k<len(nu["orderings"]) else hscale*(
                np.asarray([rng.uniform(-1,1)+1j*rng.uniform(-1,1)
                            for _ in range(4)]))
            try:
                M,_,_=target(ordering,mlight,a21,a31)
                F,D,res=inverse(M,mass,charged,sine,fscale,t,root=root,
                                homogeneous=z,f_phase=rng.uniform(-math.pi,math.pi),chart=chart)
                Fg=el.T@F@el
                Dg=er.T@D@el
                Ye2=math.sin(scalar["beta"])*sy+math.cos(scalar["beta"])*Dg
                family={"ordering":ordering,"mlight_eV":mlight,
                        "alpha21_rad":a21,"alpha31_rad":a31,
                        "root":root,"chart":chart,"t":(t.real,t.imag),
                        "homogeneous_coefficients":encode(z.reshape(1,4)),
                        "F_mass_basis":encode(F),"D_mass_basis":encode(D),
                        "analytical_closure":res,
                        "max_abs_yukawa":float(max(np.max(abs(Fg)),np.max(abs(Dg)),
                                                     np.max(abs(Ye2))))}
                card=prepare_card(original,scalar,lag,Fg,Ye2)
                family["status"]="NATIVE_REVALIDATION_REQUIRED"
                filename=f"{record['scalar_id']}_{ordering}_f{k:02d}.in"
                family["card"]=filename
                cards.append((filename,card,family))
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
                family={"ordering":ordering,"status":"ANALYTIC_INVERSE_FAILED",
                        "reason":str(exc)}
            record["families"].append(family)
    output_dir.mkdir(parents=True,exist_ok=False)
    (output_dir/"cards").mkdir()
    for filename,card,record in cards:
        path=output_dir/"cards"/filename
        with path.open("x") as out:out.write(card)
        record["card_sha256"]=sha(path)
    result=dict(schema_version=2,
                status="SCALAR_NEUTRINO_PROPOSALS_NATIVE_UNKNOWN",
                config_sha256=sha(config_path),
                template_sha256=sha(template_path),
                snapshot_sha256=sha(snapshot_path),
                n_scalar_attempted=len(rows),
                n_scalar_preselected=sum(r["gates"]["status"]==
                    "ANALYTIC_SCALAR_PRESELECTION" for r in rows),
                n_native_cards_written=len(cards),
                n_native_points_validated=0,
                n_globally_allowed=None,
                sampling_note="Arms are proposal priors; count scalar denominators separately from F/D families",
                native_FZEEIN_propagation="REQUIRES_NEW_NATIVE_SOURCE_PROOF",
                full_Zee_LFV="UNKNOWN",
                outputs=rows)
    with (output_dir/"proposals.json").open("x") as out:
        json.dump(result,out,indent=2,allow_nan=False)
        out.write("\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("config","template","snapshot","output-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    result=generate(a.config,a.template,a.snapshot,a.output_dir)
    print(json.dumps({k:v for k,v in result.items() if k!="outputs"},
                     indent=2,allow_nan=False))


if __name__=="__main__":
    main()
