"""General restricted-Zee LFV amplitude ledger; no neutrino fitting or scan logic.

Input is one explicit physical parameter point in the charged-lepton mass
basis, not a native SPheno result. Only independently audited photon subsets
are included. Physical LFV acceptance is UNKNOWN until missing contributions
are controlled. All masses GeV, angles radians, CL/CR GeV^-2.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from LFV.core import (BASIS, OPERATOR, vertices, dipole_ledger,
                      radiative_width, neutral_width, universality)
from LFV import tree
from LFV.barr_zee import (QUARKS, SPECIES, quark_coupling, fermion_photon,
                          charged_scalar_photon, physical_trilinears)

CHANNELS={"mu_to_e_gamma":(1,0),"tau_to_e_gamma":(2,0),"tau_to_mu_gamma":(2,1)}
SCALAR_KEYS={"v","beta","alpha","gamma","mh1","mh2","ma","mH1","mH2"}
MISSING=(
    "Unmatched full Zee W/Goldstone/ghost/charged-current LFV gauge aggregate (F-dependent and doublet sectors)",
    "Fermion Z-photon and physical mixed charged-scalar Z-photon Barr-Zee",
    "Other electroweak two-loop bosonic, WZ and gauge-kinetic pieces",
    "QED running, real/virtual radiative effects and scheme/scale uncertainty",
    "Exact external-invariant and finite-neutrino-mass remainders",
    "Coherent tree, dipole, nondipole photon, Z/scalar penguin and box contributions to 3-lepton decays",
    "Muon-electron conversion: semileptonic scalar/vector operators, photonic/Z penguins, boxes and nuclear matching",
    "Leptonic universality doublet-D and mixed F-D operators and higher powers",
    "Muonium-antimuonium exchange and separate lepton-flavor-preserving precision bounds",
)


def matrix(x, name):
    if not isinstance(x,dict) or set(x)!={"real","imag"}:
        raise ValueError(f"{name} requires complex real/imag arrays")
    a=np.asarray(x["real"],float)+1j*np.asarray(x["imag"],float)
    if a.shape!=(3,3) or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a finite complex 3x3 matrix")
    return a


def val(z):
    return {"real":float(z.real),"imag":float(z.imag)}


def positive_array(a, count, name):
    a=np.asarray(a,float)
    if a.shape!=(count,) or not np.isfinite(a).all() or np.min(a)<=0:
        raise ValueError(f"{name} requires {count} finite positive entries")
    return a


def input_parameters(point):
    if not isinstance(point,dict):
        raise ValueError("Expected point dictionary")
    required={"point_id","scalar","F","D","internal_lepton_masses_GeV",
              "pole_lepton_masses_GeV","fermion_loop_masses_GeV",
              "quark_yukawa_masses_GeV","alpha_em",
              "total_lepton_widths_GeV","mass_yukawa_scheme"}
    if not required<=point.keys() or point.keys()-required-{"potential","tree_quadrature_order","GF_GeV_minus2","radiative_limits"}:
        raise ValueError("Unexpected/missing point fields")
    if not isinstance(point["point_id"],str) or not point["point_id"].strip():
        raise ValueError("Nonempty point_id required")
    p=point["scalar"]
    if not isinstance(p,dict) or set(p)!=SCALAR_KEYS or not all(
            isinstance(v,(float,int)) and not isinstance(v,bool) and
            math.isfinite(v) for v in p.values()):
        raise ValueError("Scalar point must supply all finite physical inputs")
    if (min(p[k] for k in ("v","mh1","mh2","ma","mH1","mH2"))<=0 or
            p["mh1"]>p["mh2"] or p["mH1"]>p["mH2"] or
            not 1e-4<p["beta"]<math.pi/2-1e-4):
        raise ValueError("Invalid scalar spectrum/vev/beta")
    if not isinstance(point["mass_yukawa_scheme"],str) or not point["mass_yukawa_scheme"].strip():
        raise ValueError("Explicit mass/Yukawa scheme label is required")
    F,D=matrix(point["F"],"F"),matrix(point["D"],"D")
    if np.max(abs(F+F.T))>1e-10*max(1.,np.max(abs(F))):
        raise ValueError("F must be antisymmetric in the lepton mass basis")
    ml=positive_array(point["internal_lepton_masses_GeV"],3,"internal lepton masses")
    pole=positive_array(point["pole_lepton_masses_GeV"],3,"pole lepton masses")
    if not np.all(np.diff(ml)>0) or not np.all(np.diff(pole)>0):
        raise ValueError("Lepton flavor indices must be e,mu,tau")
    fm=point["fermion_loop_masses_GeV"]
    if not isinstance(fm,dict) or set(fm)!=set(SPECIES):
        raise ValueError("Provide all nine fermion loop masses, with explicit scheme")
    for k in fm:
        if not isinstance(fm[k],(int,float)) or not math.isfinite(fm[k]) or fm[k]<=0:
            raise ValueError("Invalid fermion loop mass")
    qy=point["quark_yukawa_masses_GeV"]
    if (not isinstance(qy,dict) or set(qy)!=set(QUARKS) or
            any(not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0
                for x in qy.values())):
        raise ValueError("Provide six positive quark Yukawa mass inputs in stated scheme")
    widths=point["total_lepton_widths_GeV"]
    if not isinstance(widths,dict) or set(widths)!={"mu","tau"} or any(
            not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0
            for x in widths.values()):
        raise ValueError("Muon and tau total widths must be explicit")
    alpha=point["alpha_em"]
    if not isinstance(alpha,(float,int)) or not math.isfinite(alpha) or not 0<alpha<1:
        raise ValueError("Positive finite alpha_em required")
    if min(p[k] for k in ("mh1","mh2","ma","mH1","mH2"))<20*ml[2]:
        raise ValueError("Outside validated external-invariant expansion")
    if "tree_quadrature_order" in point and (type(point["tree_quadrature_order"]) is not int
        or not 8<=point["tree_quadrature_order"]<=48):
        raise ValueError("tree_quadrature_order must be an integer from 8 to 48")
    if "GF_GeV_minus2" in point:
        gf=point["GF_GeV_minus2"]
        if not isinstance(gf,(int,float)) or not math.isfinite(gf) or gf<=0:
            raise ValueError("Fermi constant must be finite and positive")
    limits=point.get("radiative_limits",{})
    if not isinstance(limits,dict) or set(limits)-set(CHANNELS):
        raise ValueError("Unsupported radiative limits")
    for ch,item in limits.items():
        if not isinstance(item,dict) or set(item)-{
            "BR_90CL","source","missing_norm_upper_GeV_minus2","bound_source"}:
            raise ValueError(f"Invalid radiative limit record {ch}")
        if (not isinstance(item.get("BR_90CL"),(int,float)) or
            not math.isfinite(item["BR_90CL"]) or
            not 0<item["BR_90CL"]<1 or
            not isinstance(item.get("source"),str) or not item["source"].strip()):
            raise ValueError(f"Bad confidence level/source for {ch}")
        if "missing_norm_upper_GeV_minus2" in item:
            r=item["missing_norm_upper_GeV_minus2"]
            if (not isinstance(r,(int,float)) or not math.isfinite(r) or r<0
                or not isinstance(item.get("bound_source"),str)
                or not item["bound_source"].strip()):
                raise ValueError("Omitted-amplitude upper bound needs finite norm and source")
    potential=point.get("potential")
    if potential is not None and not isinstance(potential,dict):
        raise ValueError("Potential must be a dict when provided")
    return p,F,D,ml,pole,fm,qy,widths,alpha,potential



def amplitude_budget(norm, mi, mj, width, alpha, bound):
    """Rigorous two-chirality triangle inequality; conditional only."""
    factor=alpha*mi**5/4*(1-(mj/mi)**2)**3/width
    allowed=math.sqrt(bound["BR_90CL"]/factor)
    upper_missing=bound.get("missing_norm_upper_GeV_minus2")
    result={
        "experimental_limit_BR_90CL":bound["BR_90CL"],
        "experimental_source":bound["source"],
        "required_max_missing_norm_for_sufficient_pass_GeV_minus2":
            max(0.,allowed-norm),
        "retained_amplitude_norm_GeV_minus2":norm,
        "experimental_threshold_norm_GeV_minus2":allowed,
        "status":"UNKNOWN_UNBOUNDED_MISSING",
    }
    if upper_missing is not None:
        lowest=max(0.,norm-upper_missing)
        highest=norm+upper_missing
        decision=("PASS_IF_DECLARED_BOUND_VALID" if highest<allowed else
                  "FAIL_IF_DECLARED_BOUND_VALID" if lowest>allowed else
                  "UNKNOWN_EVEN_WITH_DECLARED_BOUND")
        result.update(status=decision,
                      omitted_norm_upper_GeV_minus2=upper_missing,
                      omitted_bound_source=bound["bound_source"],
                      allowed_BR_interval_if_bound_valid=[factor*lowest**2,
                                                          factor*highest**2],
                      caveat="The source's amplitude-norm bound has not been independently verified by this code.")
    return result


def evaluate(point):
    p,F,D,ml,pole,fm,qy,widths,alpha,potential=input_parameters(point)
    neutral=vertices(D,p,ml)
    scalar_G=physical_trilinears(p,potential) if potential is not None else None
    result={}
    for channel,(i,j) in CHANNELS.items():
        rows={}
        one=dipole_ledger(F,D,p,ml,i,j,external=pole)
        for label,row in one.items():
            rows["1loop_"+label]=dict(C_L=val(row["C_L"]),C_R=val(row["C_R"]),
                                      loop_order=1,source=row["provenance"])
        for s,(C,ms) in neutral.items():
            external=np.array([C[j,i],C[i,j].conjugate()])
            for name in SPECIES:
                inner=quark_coupling(name,qy[name],s,p) if name in QUARKS else C[{"electron":0,"muon":1,"tau":2}[name],{"electron":0,"muon":1,"tau":2}[name]]
                nc,q=QUARKS[name] if name in QUARKS else (1,-1)
                amp,err=fermion_photon(external,inner,fm[name],ms,pole[i],alpha,nc,q)
                rows[f"2loop_{s}_{name}_gamma"]=dict(
                    C_L=val(amp[0]),C_R=val(amp[1]),loop_order=2,
                    source="module4 MATCHING_AUDIT: fermion/photon; signed Lavoura+Ilisie/Chang",
                    numerical_error_GeV_minus2=float(err))
            if scalar_G is not None:
                for n,name in enumerate(("H1","H2")):
                    g=0. if s=="A" else float(scalar_G[s][n])
                    amp,err=charged_scalar_photon(external,g,p["mH1" if n==0 else "mH2"],ms,pole[i],alpha)
                    rows[f"2loop_{s}_{name}_charged_gamma"]=dict(
                        C_L=val(amp[0]),C_R=val(amp[1]),loop_order=2,
                        trilinear_GeV=g,numerical_error_GeV_minus2=float(err),
                        source="module4 CHARGED_SCALAR_AUDIT: two triangles + seagull; potential-derived G")
        total=np.array([sum(complex(v[key]["real"],v[key]["imag"]) for v in rows.values())
                        for key in ("C_L","C_R")])
        n=float(np.linalg.norm(total))
        sum_abs=float(sum(math.hypot(complex(v["C_L"]["real"],v["C_L"]["imag"]).__abs__(),
                                        complex(v["C_R"]["real"],v["C_R"]["imag"]).__abs__())
                          for v in rows.values()))
        br=radiative_width(total,pole[i],pole[j],alpha)/widths["mu" if i==1 else "tau"]
        bound=(amplitude_budget(n,pole[i],pole[j],widths["mu" if i==1 else "tau"],
                                alpha,point["radiative_limits"][channel])
               if channel in point.get("radiative_limits",{}) else None)
        result[channel]=dict(
            status="PARTIAL_COHERENT_PHOTON_SUBSETS",
            C_L_GeV_minus2=val(total[0]),C_R_GeV_minus2=val(total[1]),
            amplitude_norm_GeV_minus2=n,
            cancellation_ratio=n/sum_abs if sum_abs else None,
            branching_ratio_included_topologies_only=float(br),
            included_components=rows,
            conditional_missing_amplitude_budget=bound,
            physical_branching_ratio=None,experimental_gate="UNKNOWN")
    higgs_lfv={}
    for name,(C,mass) in neutral.items():
        higgs_lfv[name]={}
        for label,i,j in (("e_mu",0,1),("e_tau",0,2),("mu_tau",1,2)):
            higgs_lfv[name][label]={
                "partial_width_two_charges_GeV":2*neutral_width(C,mass,i,j,pole),
                "physical_branching_ratio":None,
                "status":"LO_NEUTRAL_SCALAR_PARTIAL_WIDTH_ONLY"}
    singlet_universality=(universality(F,p,point["GF_GeV_minus2"])
                          if "GF_GeV_minus2" in point else None)
    three_lepton={}
    if "tree_quadrature_order" in point:
        order=point["tree_quadrature_order"]
        for label,i,a,b in (("mu_to_3e",1,0,0),("tau_to_3e",2,0,0),
                             ("tau_to_3mu",2,1,1),("tau_to_mu_ee",2,1,0),
                             ("tau_to_e_mumu",2,0,1)):
            low=tree.width(neutral,i,a,b,pole,order)
            high=tree.width(neutral,i,a,b,pole,2*order)
            residual=abs(high-low)/max(abs(high),np.finfo(float).tiny)
            three_lepton[label]=dict(
                status="TREE_SUBSET_QUADRATURE_PASS" if residual<1e-4 else "UNKNOWN_NUMERICAL",
                BR_neutral_tree_only=high/widths["mu" if i==1 else "tau"],
                quadrature_relative_shift=residual,
                physical_branching_ratio=None,experimental_gate="UNKNOWN",
                note="No dipole, off-shell gamma/Z penguins, boxes or coherent interference included")
    return dict(point_id=point["point_id"],status="PARTIAL_PREDICTION",
                basis=BASIS,operator=OPERATOR,
                input_sha256=hashlib.sha256(json.dumps(point,sort_keys=True,separators=(",",":"),
                                                        allow_nan=False).encode()).hexdigest(),
                mass_yukawa_scheme=point["mass_yukawa_scheme"],
                included="seven scalar one-loop + 27 fermion-photon Barr-Zee per channel"
                         +(" + six physical charged-scalar-photon rows" if scalar_G is not None else ""),
                charged_scalar_photon="INCLUDED_FROM_RESTRICTED_POTENTIAL" if scalar_G is not None
                                       else "UNKNOWN_POTENTIAL_NOT_SUPPLIED",
                processes=result,three_lepton_tree=three_lepton,
                neutral_higgs_lfv_partial_widths=higgs_lfv,
                universality_F_subset=singlet_universality,
                missing_contributions=list(MISSING),
                full_lfv_gate="UNKNOWN",global_zee_acceptance=None,
                caveat="No experimental pass/exclusion is possible until omitted amplitudes and scheme uncertainties are controlled.")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True,help="Must not already exist")
    args=ap.parse_args()
    point=json.loads(args.input.read_text())
    report=evaluate(point)
    with args.output.open("x") as out:
        json.dump(report,out,indent=2,allow_nan=False)
        out.write("\n")
    print(json.dumps({"status":report["status"],"point_id":report["point_id"],
                      "charged_scalar_photon":report["charged_scalar_photon"],
                      "processes":{k:v["branching_ratio_included_topologies_only"]
                                   for k,v in report["processes"].items()},
                      "three_lepton_tree":{k:v["BR_neutral_tree_only"] for k,v in report["three_lepton_tree"].items()},
                      "full_lfv_gate":report["full_lfv_gate"]},indent=2))


if __name__=="__main__":
    main()
