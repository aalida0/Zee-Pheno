"""Single reusable, source-bound Zee pre-scan *theory* gate for ANY native card.

Completes one workflow for the restricted CP-conserving Zee scalar
theory at the input's own running scale:
- complete all-field quartic BFB sufficient certificate/necessary fails
- all published high-energy tree-level scalar 2->2 unitarity eigenvalues
- all five complex running quark/lepton Yukawa spectral norms
- independent mixed-charged-scalar Grimus (2008) one-loop S,T,U +
  correlated PDG2025 3D ellipsoid
- exact all-complex-field global TREE EW vacuum sufficient lower-bound
  proof for a TREE-TADPOLE DEFINED potential, with native SPheno
  loop-tadpole matching reported SEPARATELY.

The defined tree-level consistency scheme is an explicit necessary
assumption, not an unverified statement that the native one-loop
effective potential has the same global minimum. No UV running
above native Q is required by the low-energy phenom scan's theory
scope. If users request a UV cutoff > native Q, this gate deliberately
marks the UV check UNKNOWN and blocks promotion until supplied.

The gate is re-evaluated *per new SPheno spectrum*. No s00017
masses, mixing matrices, couplings or neutrino-fit data are reused.
No data reported here constitutes full model viability: the B/LFV/
Higgs experimental constraints remain independent.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ZeeHiggs.spheno_audit import read_slha,sha256
from Scan.zee_scalar_theory import evaluate as scalar_evaluate,NAMES
from Scan.zee_ewpo import block_1d,correlated_pdg2025
from Scan.zee_oblique_T_independent import source_point_from_blocks,independent_t
from Scan.zee_oblique_STU_independent import independent_S_U,MAX_DELTA_ABSOLUTE
from Scan.zee_chargebreaking_vacuum import source_inputs_from_blocks
from Scan.zee_vacuum_global_sufficient import from_native_point
from Scan.zee_yukawa_perturbativity import analyze_matrices

REQUIRED=("point_id","ordering","card","native","spectrum")
ORIGINAL_BINARY_SHA256="2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52"


def _path(root,value):
    if not isinstance(value,str) or not value.strip():
        raise ValueError("Missing manifest file path")
    return (root/value).resolve()


def _finite(v):
    if not math.isfinite(v):
        raise ValueError("Nonfinite source scalar parameter")
    return float(v)


def _quartic_consistency(blocks,card):
    if not {"MINPAR","ZEEQUARTIC","ZEEPOT","SPHENO","HMIX"}<=set(blocks):
        raise ValueError("Missing original running scalar input/output blocks")
    a=block_1d(blocks["MINPAR"],"MINPAR")
    q=block_1d(blocks["ZEEQUARTIC"],"ZEEQUARTIC")
    z=block_1d(blocks["ZEEPOT"],"ZEEPOT")
    opts=block_1d(blocks["SPHENO"],"SPHENO")
    if opts[38][0]!=2:
        raise ValueError("Native source SPheno two-loop RGE flag not enabled")
    if opts[57][0]!=1:
        raise ValueError("Native SPheno EWPO computation not enabled")
    scale=opts[33][0]
    if scale<50 or scale>1e6:
        raise ValueError("Unexpected native running Q")
    quartics={}
    for idx,name in enumerate(NAMES,1):
        x,y=a[idx][0],q[idx][0]
        if not math.isclose(x,y,abs_tol=5e-8,rel_tol=3e-7):
            raise ValueError("Source MINPAR/ZEEQUARTIC mismatch "+name)
        quartics[name]=_finite(y)
    for minidx,pidx in ((11,3),(12,4),(13,5)):
        if not math.isclose(a[minidx][0],z[pidx][0],
                            rel_tol=3e-7,abs_tol=1e-4):
            raise ValueError("Source MINPAR/ZEEPOT mismatch "+str(minidx))
    return quartics,float(scale)


def screen_blocks(blocks,*,point_id="native",cutoff_Q_GeV=None):
    q,scale=_quartic_consistency(blocks,point_id)
    scalar=scalar_evaluate(q)
    p=source_point_from_blocks(blocks)
    t=independent_t(p,"running_gauge")
    su=independent_S_U(p,"running_gauge")
    actual={"S":p["native_S"],"T":p["native_T"],"U":p["native_U"]}
    ours={"S":su["S"],"T":t["independent_scalar_T"],"U":su["U"]}
    diffs={k:ours[k]-actual[k] for k in ours}
    theory_match=all(abs(diffs[k])<MAX_DELTA_ABSOLUTE[k] for k in diffs)
    experimental=correlated_pdg2025(actual)
    y=analyze_matrices(blocks,native_Q_GeV=scale,
                       cutoff_Q_GeV=cutoff_Q_GeV or scale)
    vacuum=from_native_point(source_inputs_from_blocks(blocks))
    qst=scalar["BFB"]["full_Zee_quartic_BFB"]
    unitarity=scalar["unitarity"]["tree_high_energy_full_Zee_unitarity"]
    native_vacuum=vacuum["matched_native_potential"]
    defined_tree=vacuum["defined_tree_tadpole_potential"]
    # The LOW-SCALE tree and one-loop-oblique theory scope is complete
    # in the declared hybrid perturbative prescription; full one-loop
    # effective potential, high-Q flow and nonoblique EW loops are
    # explicitly out of scope and NOT silently treated as PASSES.
    mandatory={
        "quartic_BFB_sufficient":qst=="CERTIFIED_SUFFICIENT_PASS",
        "scalar_high_energy_tree_unitarity":unitarity=="PASS",
        "local_scalar_4pi_prior":scalar["max_abs_quartic_lt_4pi_conventional_prior"]=="PASS",
        "local_Yukawa_singular_values":y["local_Yukawa_perturbativity"]=="PASS",
        "independent_mixed_charged_scalar_STU":theory_match,
        "correlated_PDG2025_STU":experimental["status"].endswith("_PASS"),
        "defined_tree_global_EW_sufficient":
             defined_tree.get("certified_global_tree_EW") is True,
    }
    required_uv=not y["RGE_to_cutoff_proven"]
    mandatory["optional_UV_RGE_cutoff_if_requested"]=not required_uv
    ready=all(mandatory.values())
    return {
        "point_id":point_id,
        "status":("DECLARED_LOW_SCALE_THEORY_SCOPE_PASS" if ready else
                  "DECLARED_LOW_SCALE_THEORY_SCOPE_FAIL_OR_UNPROVEN"),
        "passed_in_declared_tree_and_one_loop_scalar_scope":ready,
        "requirements":mandatory,
        "scale_GeV":scale,
        "requested_UV_cutoff_GeV":y["high_scale_cutoff_GeV"],
        "quartics":q,
        "BFB":scalar["BFB"],
        "scalar_unitarity":scalar["unitarity"],
        "Yukawa_perturbativity":y,
        "scalar_one_loop_STU":{
            "SPheno_native":actual,"independent_Grimus_2008":ours,
            "difference":diffs,"tolerance":MAX_DELTA_ABSOLUTE,
            "source_match_within_declared_absolute_tolerances":theory_match,
            "all_two_charged_states_included":True,
        },
        "PDG2025_correlated_STU":experimental,
        "global_tree_vacuum":{
            "defined_tree_potential_certificate":defined_tree,
            "original_native_running_potential_certificate":native_vacuum,
            "native_loop_tadpole_spheno_match":
               native_vacuum["status"]=="CERTIFIED_SUFFICIENT_GLOBAL_TREE_EW",
            "one_loop_effective_potential_global_minimum":"NOT_COMPUTED",
        },
        "physics_scope":{
            "tree_scalar_potential_defined_by_analytic_tadpoles":"YES",
            "native_loop_corrected_scalar_potential_globality":
                "CERTIFIED" if native_vacuum.get("certified_global_tree_EW") else "UNKNOWN",
            "finite_energy_unitarity":"NOT_COMPUTED",
            "Yukawa_RGE_above_source_Q":
                "NOT_REQUESTED" if not required_uv else
                "UNKNOWN_REQUIRES_MODEL_SOURCE_RGE_TRAJECTORY",
            "nonoblique_EW_observables":"NOT_COMPUTED",
            "Higgs_FLAVOR_LFV_experimental_gate":"NOT_PART_OF_THIS_THEORY_GATE",
        },
        "full_physical_model_viability":"UNKNOWN",
        "globally_allowed":None,
    }


def screen_one(entry,root,expected_binary_sha=ORIGINAL_BINARY_SHA256,
               cutoff_Q_GeV=None):
    if not isinstance(entry,dict) or not set(REQUIRED)<=set(entry):
        raise ValueError("Manifest entry missing card/native/spectrum/point/ordering")
    if entry["ordering"] not in ("NO","IO"):
        raise ValueError("Invalid neutrino ordering")
    paths={name:_path(root,entry[name]) for name in ("card","native","spectrum")}
    if not all(p.is_file() for p in paths.values()):
        raise ValueError("Missing source card, native or SPheno spectrum")
    native=json.loads(paths["native"].read_text())
    card_sha=sha256(paths["card"])
    if native.get("input_sha256")!=card_sha:
        raise ValueError("Source native/card SHA mismatch")
    if native.get("production_binary_sha256")!=expected_binary_sha:
        raise ValueError("Source original SPheno binary mismatch")
    if native.get("spectrum_status")!="PASS":
        raise ValueError("Source native spectrum is not passed")
    blocks,_=read_slha(paths["spectrum"].read_text())
    # A digest of the spectrum alone is not a proof it came from this
    # card. Require a SECOND, native-JSON mass table generated by the
    # authenticated SPheno execution to agree with every physical
    # scalar in the separate SLHA. New runs MUST preserve this table.
    native_mass=native.get("masses_GeV")
    if not isinstance(native_mass,dict):
        raise ValueError("Native source lacks independent physical mass table")
    parsed_mass=block_1d(blocks["MASS"],"MASS")
    for pdg in (25,35,36,37,9000037):
        sv=native_mass.get(str(pdg))
        if type(sv) not in (float,int) or not math.isfinite(sv):
            raise ValueError("Native separate scalar mass missing "+str(pdg))
        actual=parsed_mass[pdg][0]
        if not math.isclose(actual,sv,rel_tol=2e-7,abs_tol=3e-5):
            raise ValueError("Native source/spectrum scalar mass mismatch "+str(pdg))
    if native.get("spectrum_sha256") is not None and (
            native["spectrum_sha256"]!=sha256(paths["spectrum"])):
        raise ValueError("Native source/spectrum digest mismatch")
    r=screen_blocks(blocks,point_id=entry["point_id"],
                    cutoff_Q_GeV=cutoff_Q_GeV)
    r["point_id"]=entry["point_id"]
    r["ordering"]=entry["ordering"]
    r["provenance"]={
        "card_sha256":card_sha,
        "native_source_json_sha256":sha256(paths["native"]),
        "spheno_spectrum_sha256":sha256(paths["spectrum"]),
        "production_binary_sha256":expected_binary_sha,
        "source_card_matched_native":True,
        "spectrum_SHA256_computed":True,
        "source_native_independent_mass_concordance":True,
        "optional_native_JSON_spectrum_digest_matched":
            native.get("spectrum_sha256") is not None,
    }
    return r


def evaluate_manifest(manifest,cutoff_Q_GeV=None,
                      expected_binary_sha=ORIGINAL_BINARY_SHA256):
    path=Path(manifest).resolve()
    obj=json.loads(path.read_text())
    if obj.get("schema_version")!=1 or not isinstance(obj.get("points"),list):
        raise ValueError("Unsupported source-bound manifest")
    if not 1<=len(obj["points"])<=1000:
        raise ValueError("Need 1..1000 real native points")
    outputs=[]
    for row in obj["points"]:
        try:
            outputs.append(screen_one(row,path.parent,expected_binary_sha,
                                      cutoff_Q_GeV))
        except (ValueError,KeyError,IndexError,TypeError,
                ArithmeticError,AssertionError,OSError) as e:
            outputs.append({
                "point_id":row.get("point_id"),
                "status":"NOT_EVALUATED_SOURCE_OR_CONVENTION_ERROR",
                "error_type":type(e).__name__,"error":str(e),
                "passed_in_declared_tree_and_one_loop_scalar_scope":False,
                "globally_allowed":None,
            })
    return {
        "status":"SOURCE_BOUND_PRE_SCAN_THEORY_GATE_NO_GLOBAL_ACCEPTANCE",
        "input_manifest":str(path),
        "input_manifest_sha256":sha256(path),
        "native_binary_sha256":expected_binary_sha,
        "n_points":len(outputs),
        "n_theory_scope_pass":sum(x.get(
            "passed_in_declared_tree_and_one_loop_scalar_scope") is True
            for x in outputs),
        "n_source_or_convention_invalid":sum(x["status"]==
            "NOT_EVALUATED_SOURCE_OR_CONVENTION_ERROR" for x in outputs),
        "n_global_allowed":None,
        "points":outputs,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--cutoff-Gev",type=float)
    p.add_argument("--expected-binary-sha256",default=ORIGINAL_BINARY_SHA256)
    a=p.parse_args()
    result=evaluate_manifest(a.manifest,a.cutoff_Gev,
                             a.expected_binary_sha256)
    with a.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k!="points"},indent=2))


if __name__=="__main__":main()
