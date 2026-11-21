"""HiggsTools tabulated 13 TeV associated charged-Higgs production in Zee.

For *Type-II quark Yukawas* the charged singlet has no tb Yukawa.
The physical H_i tb vertex is the Type-II charged-doublet coupling
times its coefficient in the native, full 3x3 CHARGEMIX (Goldstone,H1,H2).
This uses the official HiggsPredictions ppHpmtb 4FS/5FS matched reference
tables for a declared low-scale phenomenological prediction.

The reference is only valid for 13 TeV, approx 145--2000 GeV; here
restrict further to above on-shell top decay (m_H > m_t + m_b) so
BR(t->H+b)=0 can be source-justified. Out of domain -> UNKNOWN/error,
never a fabricated zero. No exhaustive Zee HiggsBounds coverage claimed.
"""
import math
from pathlib import Path

from ZeeHiggs.spheno_audit import CHARGED, read_slha
from Scan.zee_oblique_T_independent import read_matrix
from Scan.zee_ewpo import block_1d


def native_charged_couplings(source):
    """Derive Type-II kappa_t/b from measured SPheno vevs/mixing, NOT gamma guesses."""
    blocks,_=read_slha(Path(source["spectrum_path"]).read_text())
    hmix=block_1d(blocks["HMIX"],"HMIX")
    mass=block_1d(blocks["MASS"],"MASS")
    mixing=read_matrix(blocks,"CHARGEMIX",3)
    vd,vu=(hmix[k][0] for k in (102,103))
    if min(vd,vu)<=0 or not math.isfinite(vd+vu):
        raise ValueError("Non-positive or nonfinite source electroweak vevs")
    v=math.hypot(vd,vu)
    cb,sb=vd/v,vu/v
    # Rows = mass eigenstates; gauge columns = (Hd-,Hu-,singlet).
    # Row 0 is the Goldstone (cb,sb,0). H_D=(sb,-cb,0)
    # up to an unphysical common sign for each H_i.
    if (abs(mixing[0,0]-cb)>3e-3 or
        abs(mixing[0,1]-sb)>3e-3 or abs(mixing[0,2])>3e-3):
        raise ValueError("Native Goldstone/charged gauge basis mismatch")
    doublet=[float(mixing[row,0]*sb-mixing[row,1]*cb)
             for row in (1,2)]
    if not math.isclose(sum(d*d for d in doublet),1.,rel_tol=0,abs_tol=3e-3):
        raise ValueError("Native charged doublet strength does not sum to unity")
    out={}
    for i,pdg in enumerate(CHARGED):
        mh=source["masses_GeV"][str(pdg)]
        if not math.isclose(mass[pdg][0],mh,rel_tol=1e-8,abs_tol=1e-5):
            raise ValueError("Source charged mass/spectrum mismatch")
        mt=mass[6][0]
        mb=mass[5][0]
        if not (max(200.,mt+mb)<mh<=2000.):
            raise ValueError("Outside justified HiggsTools Hpmtb tabulation/top decay domain")
        d=doublet[i]
        out[pdg]={
            "mass_GeV":mh,
            "native_doublet_amplitude":d,
            "native_singlet_fraction_squared":float(mixing[i+1,2]**2),
            "kappa_t":d*vd/vu,
            "kappa_b":d*vu/vd,
            "br_top_to_Hpb":0.0,
        }
    return {
        "source":source["source"],"point_id":source["point_id"],
        "native_mix_basis":"CHARGEMIX rows Goldstone,H1,H2; columns Hd,Hu,charged singlet",
        "charged_doublet_projection":"H_D=(vu,-vd,0)/hypot(vd,vu)",
        "tan_beta_native":vu/vd,
        "mass_eigenstates":out,
        "status":"SOURCE_MATCHED_TYPE_II_CHARGED_COUPLING_PROJECTIONS",
    }


def production_reference(source, backend):
    """Numerical official HiggsTools reference cross sections (pb), not MG5."""
    couplings=native_charged_couplings(source)
    result={}
    for pdg, c in couplings["mass_eigenstates"].items():
        rate=float(backend(c["mass_GeV"],c["kappa_t"],c["kappa_b"],0.))
        tanb=couplings["tan_beta_native"]
        base=float(backend(c["mass_GeV"],1/tanb,tanb,0.))
        if (not math.isfinite(rate) or not math.isfinite(base) or
            rate<0 or base<=0):
            raise ValueError("HiggsTools tabulated Hpmtb rate nonfinite/invalid")
        expected=base*c["native_doublet_amplitude"]**2
        if not math.isclose(rate,expected,rel_tol=3e-5,abs_tol=1e-12):
            raise ValueError("HiggsTools reference violates shared-doublet quadratic scaling")
        result[pdg]={**c,
            "sigma_pp_Hpmtb_13TeV_pb":rate,
            "pure_doublet_reference_same_mass_pb":base,
            "source":"HiggsTools EffectiveCouplingCxns.ppHpmtb LHC13 tabulated 4FS/5FS",
            "theory_pdf_scale_uncertainty":"NOT_ASSIGNED",
            "production_model":"Type-II doublet projection; no charged-singlet Yukawas",
        }
    return {
        "point_id":source["point_id"],
        "status":"HIGGSTOOLS_TYPE_II_ZEE_CHARGED_13TEV_PRODUCTION_REFERENCE",
        "charged_mass_eigenstates":result,
        "doublet_projections":couplings,
        "higgstools_cross_sections":{
            f"cxn_{pdg}_Hpmtb_LHC13":row["sigma_pp_Hpmtb_13TeV_pb"]
            for pdg,row in result.items()},
        "full_Zee_charged_production_complete":False,
        "full_HiggsBounds_95CL":None,
        "globally_allowed":None,
    }


def production_reference_official(source):
    import Higgs.predictions as HP

    def table(mass,kt,kb,top_br):
        return HP.EffectiveCouplingCxns.ppHpmtb(
            HP.Collider.LHC13,mass,kt,kb,top_br)
    return production_reference(source,table)
