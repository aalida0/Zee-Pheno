"""Conservative scalar-potential diagnostics for a restricted Zee point.

These are necessary doublet-only BFB inequalities and a conventional
absolute-quartic threshold, NOT full Zee BFB, global vacuum or unitarity.
An externally declared potential is not a native SPheno attestation.
"""
import math

from Model.model import scalar_parameters


SCALARS=("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
         "K1","K2","K12","Mu4")
INVERSION=("m12sq_GeV2","lambda_c","K1","K2","K12")


def check_potential(scalar,potential,inversion):
    if not isinstance(inversion,dict) or set(inversion)!=set(INVERSION):
        raise ValueError("Explicit scalar inversion controls are required")
    if (not isinstance(potential,dict) or set(potential)!=set(SCALARS)
        or any(type(v) not in (float,int) or not math.isfinite(v)
               for v in inversion.values()) or
        any(type(v) not in (float,int) or not math.isfinite(v)
            for v in potential.values())):
        raise ValueError("Incomplete/nonfinite restricted Zee scalar data")
    # Native SPheno may report gamma outside (-pi/2,pi/2).
    # gamma -> gamma +/- pi flips a physical charged eigenstate phase;
    # sin², cos², sin*cos, trilinears and reconstructed scalar potential
    # are invariant. Do NOT edit native evidence or LFV source parameters.
    canonical=dict(scalar)
    original_gamma=canonical["gamma"]
    canonical["gamma"]=(original_gamma+math.pi/2)%math.pi-math.pi/2
    expected=scalar_parameters(
        **canonical,m12sq=inversion["m12sq_GeV2"],
        lc=inversion["lambda_c"],
        k1=inversion["K1"],k2=inversion["K2"],k12=inversion["K12"])
    discrepancy=max(abs(expected[k]-potential[k])/max(1.,abs(expected[k]))
                    for k in SCALARS)
    # Require the physical spectrum/potential relation to be source-matched
    # internally before reporting even approximate scalar observables.
    if discrepancy>1e-8:
        raise ValueError("Declared Zee potential inconsistent with scalar masses and mixing")
    l1,l2,l3,l4,l5=(potential[k] for k in
                    ("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"))
    doublet_bfb=(l1>0 and l2>0 and
                 l3 > -math.sqrt(l1*l2) and
                 l3+l4-abs(l5) > -math.sqrt(l1*l2))
    quartics=[potential[k] for k in SCALARS if k!="Mu4"]
    quartics.append(inversion["lambda_c"])
    max_quartic=max(abs(v) for v in quartics)
    return {
        "internal_tree_scalar_inversion":"PASS",
        "original_native_gamma_rad":original_gamma,
        "principal_inversion_gamma_rad":canonical["gamma"],
        "gamma_mod_pi_equivalence":"CHARGED_FIELD_REPHASING_NOT_A_PHYSICS_RETUNE",
        "max_potential_relative_residual":discrepancy,
        "doublet_subpotential_BFB_necessary":
            "PASS" if doublet_bfb else "FAIL",
        "quartic_abs_lt_4pi_convention":
            "PASS" if max_quartic<4*math.pi else "FAIL",
        "max_abs_quartic":max_quartic,
        "full_Zee_BFB":"UNKNOWN",
        "full_Zee_unitarity":"UNKNOWN",
        "full_Zee_vacuum_globality":"UNKNOWN",
        "native_scalar_potential_provenance":"UNKNOWN",
    }
