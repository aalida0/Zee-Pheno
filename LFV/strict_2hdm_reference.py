"""Strict singlet-decoupling reference adapter for public general-2HDM code.

No replacement for a Zee calculation: F = gamma = Mu4 = 0 is mandatory.
Original public code: https://gitlab.com/jbrod/general-2hdm-pheno
Authors' cL/cR convention maps to ours as C=-c/(8*pi**2*v**2)
(see previous aalida0/zee LFV_COVERAGE.md). Public tool must be supplied
separately, unmodified; this adapter never vendors or patches it.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from LFV.barr_zee import physical_trilinears
from LFV.constraints import input_parameters

SM_KEYS={"alphaQED","alphasMZ","sw2","MZ","MW","MTpole","mbmb",
         "mcmc","ms_at_2GeV","md_at_2GeV","mu_at_2GeV"}
QUARK_INDEX={"u":"up","d":"down","c":"charm","s":"strange","t":"top","b":"bottom"}
LEPTON_MASS_NAMES=("melectron","mmuon","mtau")
CHANNELS={"mu_to_e_gamma":(2,1),"tau_to_e_gamma":(3,1),
          "tau_to_mu_gamma":(3,2)}
ZERO=1e-13


def cp_even_rotation(beta, alpha):
    """arXiv:2410.17313v4 Eq.(9) in Zee CP-conserving neutral convention."""
    theta=beta-alpha
    return dict(theta12=theta,theta13=0.,theta23=0.)


def decoupled_input(point, sm):
    """Build author's input dictionary; quark inputs use author's running scheme.

    Source mass/Yukawa matching must be checked numerically after the
    author's class has computed its own weak-scale running quark masses.
    """
    p,F,D,ml,pole,_,_,_,_,pot=input_parameters(point)
    if (np.max(abs(F))>ZERO or abs(p["gamma"])>ZERO or pot is None or
        abs(pot.get("Mu4",math.inf))>ZERO):
        raise ValueError("Strict F=0, gamma=0, Mu4=0 and complete potential required")
    if not np.allclose(ml,pole,rtol=1e-11,atol=1e-15):
        raise ValueError("Use one consistent charged-lepton mass scheme for 2HDM overlap")
    if not isinstance(sm,dict) or set(sm)!=SM_KEYS or any(
            not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0
            for v in sm.values()):
        raise ValueError("Explicit finite author-scheme SM reference inputs required")
    if not 0<sm["sw2"]<1:
        raise ValueError("Weak mixing angle must lie in (0,1)")
    delta=p["beta"]-p["alpha"]
    gh=physical_trilinears(p,pot)
    # gamma=0: physical charged eigenstate 1 is exactly the doublet H+.
    Z3=(math.cos(delta)*gh["h"][0]+math.sin(delta)*gh["H"][0])/p["v"]
    Z7=(-math.sin(delta)*gh["h"][0]+math.cos(delta)*gh["H"][0])/p["v"]
    rho_l=D.conj().T/math.sqrt(2)  # 2015/2026 2HDM rho_l^dagger = D/sqrt(2)
    data=dict(sm,**cp_even_rotation(p["beta"],p["alpha"]))
    data.update(dict(Mh1=p["mh1"],Mh2=p["mh2"],Mh3=p["ma"],
                     MHp=p["mH1"],vev=p["v"],Z3=float(Z3),
                     absZ7=float(abs(Z7)),argZ7=0. if Z7>=0 else math.pi))
    data.update(zip(LEPTON_MASS_NAMES,map(float,pole)))
    for family in ("u","d","l"):
        for i in range(1,4):
            for j in range(1,4):
                if family=="l":
                    val=rho_l[i-1,j-1]
                else:
                    val=0j  # filled from author-run quark masses below
                data[f"absrho{family}{i}{j}"]=float(abs(val))
                data[f"argrho{family}{i}{j}"]=float(np.angle(val)) if abs(val) else 0.
    return data


def fill_type_ii_quarks(author_input, reference):
    """Match rho_q to the SAME evolved running mass used inside public code."""
    point=dict(author_input)
    t=math.tan(reference["beta"])
    # Names from author public thdm.thdm initialization / class attributes.
    source={"u":reference["mumt"],"c":reference["mc"],"t":reference["mt"],
            "d":reference["mdmt"],"s":reference["msmt"],"b":reference["mb"]}
    for family,flavors,coef in (
            ("u",("u","c","t"),1/t),("d",("d","s","b"),-t)):
        for idx,label in enumerate(flavors,1):
            z=float(source[label]*coef/author_input["vev"])
            point[f"absrho{family}{idx}{idx}"]=abs(z)
            point[f"argrho{family}{idx}{idx}"]=0. if z>=0 else math.pi
    return point


def load_public(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    spec=importlib.util.spec_from_file_location("zee_public_2hdm_reference",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module,"thdm"):
        raise ValueError("Public file must expose thdm.thdm")
    return module


def evaluate(point,sm,public_source):
    from LFV.core import vertices
    from LFV.constraints import matrix
    kwargs=decoupled_input(point,sm)
    original=load_public(public_source)
    first=original.thdm(kwargs)
    masses={key:float(getattr(first,key))
            for key in ("mumt","mc","mt","mdmt","msmt","mb")}
    kwargs=fill_type_ii_quarks(kwargs,dict(masses,beta=point["scalar"]["beta"]))
    author=original.thdm(kwargs)
    # Exact neutral-vertex check; overall A rephase is a convention.
    D=matrix(point["D"],"D")
    ml=np.asarray(point["pole_lepton_masses_GeV"],float)
    neutral=vertices(D,point["scalar"],ml)
    q1=np.array([author.q1(i) for i in (1,2,3)])
    q2=np.array([author.q2(i) for i in (1,2,3)])
    vertex_residual={}
    for k,name in enumerate(("h","H","A")):
        model=q1[k]*np.diag(ml/point["scalar"]["v"])+q2[k]*D/math.sqrt(2)
        # A in original SARAH Zee convention is the negative of public a0.
        if name=="A":
            model=-model
        difference=float(np.max(abs(model-neutral[name][0])))
        vertex_residual[name]=difference
        if difference>5e-12:
            raise ValueError(f"2HDM neutral Yukawa convention mismatch: {name}")
    # The *entire* author result is a 2HDM reference, not a Zee prediction.
    common=-1/(8*math.pi**2*point["scalar"]["v"]**2)
    dipoles={}
    for label,(i,j) in CHANNELS.items():
        cl=complex(author.cL_llpg(i,j))*common
        cr=complex(author.cR_llpg(i,j))*common
        if not np.isfinite([cl,cr]).all():
            raise ArithmeticError("Nonfinite author reference dipole")
        dipoles[label]={
            "C_L_GeV_minus2":{"real":cl.real,"imag":cl.imag},
            "C_R_GeV_minus2":{"real":cr.real,"imag":cr.imag},
            "scope":"STRICT_2HDM_ALL_AUTHOR_TERMS_ONLY"}
    return {
        "status":"STRICT_2HDM_REFERENCE_EVALUATED_NOT_ZEE",
        "point_id":point["point_id"],
        "public_source_sha256":hashlib.sha256(public_source.read_bytes()).hexdigest(),
        "public_model":"jbrod/general-2hdm-pheno thdm.py",
        "dipole_conversion":"C_X=-c_X/(8*pi^2*v^2) in archived Zee convention",
        "neutral_vertex_max_residual":max(vertex_residual.values()),
        "neutral_vertex_residuals":vertex_residual,
        "rho_l":"rho_l^dagger=D/sqrt(2)",
        "theta12":kwargs["theta12"],
        "H_plus_mass_GeV":kwargs["MHp"],
        "Z3":kwargs["Z3"],"Z7_real":(kwargs["absZ7"]*math.cos(kwargs["argZ7"])),
        "reference_dipoles":dipoles,
        "Zee_full_dipoles":None,
        "Zee_gauge_extension":"UNVALIDATED",
        "note":"The source full 2HDM result is gauge invariant only in its original model; physical Zee gauge correction is missing. See LFV/REMAINDER_AUDIT.md.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("input","sm","public-source","output"):
        parser.add_argument("--"+arg,type=Path,required=True)
    a=parser.parse_args()
    result=evaluate(json.loads(a.input.read_text()),json.loads(a.sm.read_text()),a.public_source)
    with a.output.open("x") as fh:
        json.dump(result,fh,indent=2,allow_nan=False)
        fh.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k!="reference_dipoles"},indent=2))


if __name__=="__main__":
    main()
