"""Rigorous SUFFICIENT GLOBAL TREE vacuum criterion for restricted Zee.

Let M_C^2 > 0, lambda_C >= 0 and portal K = [[k1,-k12],[-k12,k2]]
be positive semidefinite.  The full Eq.(17) Zee scalar potential obeys

 V_Zee >= V_2HDM - (mu4^2/M_C^2) * |Phi1^T i sigma2 Phi2|^2,

because M_C^2 |chi|^2 - 2 |mu4 det| |chi| >= -mu4^2|det|^2/M_C^2.
K positive semidefinite and lambda_C>=0 add nonnegative terms.

Using |det(Phi1,Phi2)|^2 = |Phi1|^2|Phi2|^2-|Phi1^dag Phi2|^2,
the lower bound is EXACTLY a softly-broken Z2 CP-conserving 2HDM with

  lambda3_eff = lambda3 - d; lambda4_eff = lambda4 + d
  lambda1,lambda2,lambda5 unchanged; d=mu4^2/M_C^2.

At a neutral EW configuration det=0 and chi=0, the inequality
is saturated. Hence if this auxiliary *tree* 2HDM has the specified
normal EW vacuum as its global minimum, the full Zee tree potential
has that global minimum, including all charged and CP phases.

The 2HDM global-vacuum condition for a normal MINIMUM is D>0,
D = m12^2*(m11^2 - k^2*m22^2)*(tanbeta-k),
k=(lambda1/lambda2)^1/4, Barroso/Ferreira/Ivanov/Santos/Silva
2013 EPJC 73:2537. This is a SUFFICIENT Zee certificate because
some Zee potentials fail K PSD or the auxiliary condition but remain
globally stable. This is NOT a proof of a one-loop global vacuum.

The tree-level m11^2,m22^2 MUST satisfy normal tadpoles at the
EXACT stated vd,vu and lambda inputs. For SPheno native outputs,
the loop-tadpole derived ZEEPOT masses may FAIL this stationarity
condition. DO NOT certify the native running potential in that case:
report the certificate only for the explicitly reconstructed
tree-consistent reference potential, and keep the native gate UNKNOWN.
"""
import math

import numpy as np

from Scan.zee_scalar_theory import bfb,validate
from Scan.zee_chargebreaking_vacuum import normal_tree_mass_squared,VDIM


SOURCE={
  "Zee_potential":"Fontes/Romao/Silva, EPJC 81 (2021) 1087, Eq.(17)",
  "Zee_url":"https://link.springer.com/article/10.1140/epjc/s10052-021-09956-2",
  "2HDM_global_discriminant":"Barroso/Ferreira/Ivanov/Santos/Silva, EPJC 73 (2013) 2537",
  "2HDM_url":"https://orbi.uliege.be/bitstream/2268/136450/1/EPJC.pdf",
  "bound":"Completing complex charged singlet square + K positive-semidefinite; lower bound by auxiliary 2HDM; equality at neutral EW",
}


def auxiliary_tree_vacuum(q,vd,vu,m12sq,mc2,mu4,*,quad_masses=None,
                          tolerance=1e-7):
    """Return a mathematical one-sided global-minimum certificate.

    If quad_masses supplied, they must satisfy tree stationary tadpoles;
    a SPheno loop-corrected tadpole *will not* pass.
    """
    q=validate(q)
    vals=(vd,vu,m12sq,mc2,mu4)
    if any(type(z) not in (float,int) or not math.isfinite(z) for z in vals):
        raise ValueError("Nonfinite Zee vacuum input")
    vd,vu,m12sq,mc2,mu4=(float(z) for z in vals)
    if vd<=0 or vu<=0 or mc2<=0 or m12sq<=0:
        return {
            "status":"UNKNOWN_SUFFICIENT_PRECONDITIONS_NOT_MET",
            "reason":"vd,vu,m12^2 and M_C^2 must be positive for published-discriminant certificate",
            "certified_global_tree_EW":False,
            "native_potential_certified":False,
        }
    analytic=normal_tree_mass_squared(q,(vd,vu),m12sq)
    if quad_masses is None:
        used=analytic
        matched=True
        matching="ANALYTIC_TREE_TADPOLE_RECONSTRUCTION"
    else:
        used=np.asarray(quad_masses,dtype=float)
        if used.shape!=(2,) or not np.isfinite(used).all():
            raise ValueError("Nonfinite or wrongly shaped m11/m22 tree tadpole inputs")
        residual=np.asarray(used)-np.asarray(analytic)
        denom=max(1.,max(abs(z) for z in analytic))
        matched=bool(np.max(abs(residual))/denom<tolerance)
        matching=("MATCHED_SOURCE_TREE_TADPOLES" if matched else
                  "SOURCE_RUNNING_TADPOLES_NOT_TREE_STATIONARY")
        if not matched:
            return {
                "status":"UNKNOWN_SOURCE_TADPOLE_SCHEME_MISMATCH",
                "tree_tadpole_masses2_GeV2":list(analytic),
                "provided_mass_m11_m22_GeV2":list(map(float,used)),
                "stationarity_max_abs_residual_GeV2":float(np.max(abs(residual))),
                "relative_stationarity_max":float(np.max(abs(residual))/denom),
                "certified_global_tree_EW":False,
                "native_potential_certified":False,
                "reason":"No global theorem can certify a nonstationary EW point for those mass inputs",
            }
    cert_bfb=bfb(q)
    k1,k2,k12=(q[z] for z in ("K1","K2","K12"))
    portal_psd=k1>=0 and k2>=0 and k1*k2-k12*k12>=0
    c_positive=q["LambdaC"]>=0
    d=mu4*mu4/mc2
    effective=dict(q,Lambda3=q["Lambda3"]-d,
                   Lambda4=q["Lambda4"]+d)
    aux_bfb=bfb({**effective,"LambdaC":max(q["LambdaC"],1e-10),
                 "K1":0.,"K2":0.,"K12":0.})
    # Pure auxiliary 2HDM BFB necessary and sufficient at tree level.
    l1,l2,l3,l4,l5=(effective[z] for z in (
        "Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"))
    root=math.sqrt(max(0,l1*l2))
    aux_bfb_strict=(
        l1>0 and l2>0 and l3+root>0 and
        l3+l4-abs(l5)+root>0)
    v2=vd*vd+vu*vu
    tanbeta=vu/vd
    m2bar=m12sq*v2/(vd*vu)
    massmat=np.array([
        [m12sq*tanbeta+l1*vd*vd,
         -m12sq+(l3+l4+l5)*vd*vu],
        [-m12sq+(l3+l4+l5)*vd*vu,
         m12sq/tanbeta+l2*vu*vu],
    ])
    cp_even_eig=np.linalg.eigvalsh(massmat)
    mA2=m2bar-l5*v2
    mHp2=m2bar-.5*(l4+l5)*v2
    normal_local_minimum=bool(
        min(cp_even_eig)>0 and mA2>0 and mHp2>0)
    k=(l1/l2)**.25 if l1>0 and l2>0 else float("nan")
    discriminant=(m12sq*(used[0]-k*k*used[1])*(tanbeta-k)
                  if math.isfinite(k) else float("nan"))
    d_positive=math.isfinite(discriminant) and discriminant>0
    sufficient=bool(portal_psd and c_positive and aux_bfb_strict and
                    normal_local_minimum and d_positive)
    return {
       "status":("CERTIFIED_SUFFICIENT_GLOBAL_TREE_EW" if sufficient else
                 "UNKNOWN_SUFFICIENT_CONDITIONS_NOT_MET"),
       "certificate_source":SOURCE,
       "source_matching":matching,
       "mass_parameters_m11_m22_GeV2":list(map(float,used)),
       "tree_tadpole_masses2_GeV2":list(analytic),
       "charge_singlet_MC2_GeV2":mc2,
       "mu4_GeV":mu4,
       "complete_square_delta_lambda3_minus_delta_lambda4":d,
       "portal_K_PSD":bool(portal_psd),
       "singlet_lambdaC_nonnegative":bool(c_positive),
       "auxiliary_2HDM_lambda3":l3,
       "auxiliary_2HDM_lambda4":l4,
       "auxiliary_2HDM_BFB_strict":bool(aux_bfb_strict),
       "auxiliary_2HDM_necessary_BFB_status":aux_bfb["BFB"]["doublet_BFB_necessary"],
       "auxiliary_CP_even_masses_squared_GeV2":[float(z) for z in cp_even_eig],
       "auxiliary_CP_odd_mass_squared_GeV2":float(mA2),
       "auxiliary_charged_mass_squared_GeV2":float(mHp2),
       "auxiliary_normal_local_minimum":normal_local_minimum,
       "2HDM_k":float(k),
       "2HDM_global_discriminant_GeV4":float(discriminant),
       "2HDM_global_discriminant_positive":bool(d_positive),
       "certified_global_tree_EW":sufficient,
       "certificate_is_sufficient_not_necessary":True,
       "all_complex_charged_and_CP_violating_directions_covered_by_inequality":True,
       "one_loop_or_thermal_global_stability_certified":False,
       "native_potential_certified":bool(sufficient and quad_masses is not None),
    }


def from_native_point(p):
    """Compare two distinct physical questions, not interchangeable."""
    q=p["quartics"]
    args=(q,p["vevs"][0],p["vevs"][1],p["m12sq"],p["mCsq"],p["mu4"])
    native=auxiliary_tree_vacuum(*args,quad_masses=(p["m1sq"],p["m2sq"]))
    tree=auxiliary_tree_vacuum(*args)
    if native.get("native_potential_certified") and not tree.get("certified_global_tree_EW"):
        raise ValueError("Inconsistent source and tree certificates")
    return {
       "matched_native_potential":native,
       "defined_tree_tadpole_potential":tree,
       "physical_global_vacuum":"UNKNOWN" if not native.get("native_potential_certified")
                                    else "CERTIFIED_GLOBAL_TREE_ONLY",
       "RGE_and_loop_level_vacuum_global":"UNKNOWN",
    }
