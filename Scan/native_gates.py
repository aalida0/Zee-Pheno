"""Native SPheno -> Neutrino and LFV evidence bridge (no retuning)."""
import json
from pathlib import Path

from Neutrino.fit import inspect_native
from LFV.native_point import from_native
from LFV.constraints import evaluate as lfv_evaluate
from Model.audit_neutrino_native import sha

MZ_GEV=91.1876
REQUIRED=("point_id","ordering","card","native","spectrum","quark_loop",
          "quark_yukawa","mass_yukawa_scheme")
OPTIONAL=("potential","potential_source","native_input_report")
UNKNOWN_SECTORS=(
    "full_Zee_LFV","full_three_lepton_LFV","mu_e_conversion",
    "B_flavor","Higgs_rates","EWPO","vacuum_global","scalar_BFB_unitarity",
    "direct_searches","collider_sensitivity","native_FZEEIN_propagation",
)

def gate(status,method,**evidence):
    assert status in ("PASS","FAIL","UNKNOWN")
    return dict(status=status,method=method,**evidence)

def check_entry(entry):
    if (not isinstance(entry,dict) or
        set(entry)-set(REQUIRED)-set(OPTIONAL) or
        not set(REQUIRED)<=set(entry)):
        raise ValueError("Unsupported/missing native point manifest fields")
    if entry["ordering"] not in ("NO","IO"):
        raise ValueError("Invalid ordering")
    for k in REQUIRED:
        if not isinstance(entry[k],str) or not entry[k].strip():
            raise ValueError(f"Empty or non-string {k}")
    if ("potential" in entry)!=("potential_source" in entry):
        raise ValueError("Potential requires a source declaration")

def evaluate_entry(entry,root):
    """A failed input binding aborts this point; no fabricated physics result."""
    check_entry(entry)
    paths={k:(root/entry[k]).resolve() for k in
           ("card","native","spectrum","quark_loop","quark_yukawa")}
    if any(not x.is_file() for x in paths.values()):
        raise FileNotFoundError("Source card/native/spectrum/quark file missing")
    potential=((root/entry["potential"]).resolve()
               if "potential" in entry else None)
    if potential and not potential.is_file():
        raise FileNotFoundError("Potential file absent")
    native=json.loads(paths["native"].read_text())
    binding=sha(paths["card"])
    if native.get("input_sha256")!=binding:
        raise ValueError("CARD_NATIVE_SHA_MISMATCH")
    propagation=gate("UNKNOWN","NO_SOURCE_MATCHED_FZEEIN_PROBE")
    if "native_input_report" in entry:
        path=(root/entry["native_input_report"]).resolve()
        if not path.is_file():
            raise FileNotFoundError("Declared native input-propagation report missing")
        proof=json.loads(path.read_text())
        if proof.get("input_sha256")!=binding:
            raise ValueError("INPUT_PROBE_REPORT_CARD_SHA_MISMATCH")
        if proof.get("status")!="NATIVE_INPUT_GATE_PASS":
            raise ValueError("INPUT_PROBE_NOT_NATIVE_PASS")
        checks=proof.get("checks",{})
        if (proof.get("native_executed") is not True or
            not isinstance(proof.get("binary_sha256"),str) or
            len(proof["binary_sha256"])!=64 or
            not isinstance(checks,dict) or
            not isinstance(proof.get("output_sha256"),dict) or
            set(proof["output_sha256"])!={"baseline","f_probe","ye2_probe"}):
            raise ValueError("Input report lacks attested binary and three native outputs")
        try:
            verified=(checks["baseline_F_relative"]<=5e-4 and
                      checks["baseline_Ye2_relative"]<=5e-4 and
                      checks["F_antisymmetric_max_abs"]<=5e-9 and
                      checks["F_probe_output_delta"]>=.5*proof["delta_f"] and
                      checks["Ye2_probe_output_delta"]>=.5*proof["delta_imye2"] and
                      proof["delta_f"]>0 and proof["delta_imye2"]>0)
        except (KeyError,TypeError) as exc:
            raise ValueError("Missing native input-propagation metrics") from exc
        if not verified:
            raise ValueError("Native F/Ye2 input propagation metrics fail published thresholds")
        propagation=gate("PASS","CARD_BOUND_FZEEIN_AND_YEE2_PROBE_NATIVE_GATE",
                         report_sha256=sha(path),
                         native_binary_sha256=proof.get("binary_sha256"),
                         probe_checks=checks)
    neutrino=inspect_native(paths["native"],entry["ordering"])
    if neutrino["status"] not in ("NATIVE_3SIGMA_SCREEN_PASS",
                                 "NATIVE_3SIGMA_SCREEN_FAIL"):
        raise ValueError("Native neutrino fit has no evaluated 3sigma status")
    point=from_native(
        paths["native"],paths["card"],paths["spectrum"],
        paths["quark_loop"],paths["quark_yukawa"],
        entry["mass_yukawa_scheme"],entry["point_id"],potential)
    scalar=point["scalar"]
    neutrino_ok=neutrino["status"]=="NATIVE_3SIGMA_SCREEN_PASS"
    # 125.25 +- 3 GeV is deliberately a loose necessary-only
    # restricted-h125 mass screen, NOT an actual Higgs likelihood.
    higgs_mass_ok=abs(scalar["mh1"]-125.25)<=3.
    cascade_ok=scalar["mH2"]>scalar["mH1"]+MZ_GEV
    gates={
        "card_native_binding":gate("PASS","SHA256_NATIVE_CARD",
                                 card_sha256=binding),
        "neutrino_1D_3sigma":gate(
            "PASS" if neutrino_ok else "FAIL",
            "NATIVE_TAKAGI_NUFIT_6_1_SIX_1D_INTERVALS",
            values=neutrino["parameters"],checks=neutrino["checks"]),
        "h125_mass_screen":gate(
            "PASS" if higgs_mass_ok else "FAIL",
            "RESTRICTED_H_IS_PDG25_WITHIN_3GEV_125_25",
            pole_mass_GeV=scalar["mh1"],full_Higgs_gate="UNKNOWN"),
        "cascade_on_shell":gate(
            "PASS" if cascade_ok else "FAIL",
            "H2_GT_H1_PLUS_MZ",mass_gap_GeV=scalar["mH2"]-scalar["mH1"]),
        "potential_provenance":gate(
            "UNKNOWN","NATIVE_POTENTIAL_MATCH_NOT_VERIFIED",
            user_declared_source=entry.get("potential_source"),
            potential_SHA256=sha(potential) if potential else None),
    }
    gates.update({name:gate("UNKNOWN","NOT_YET_VALIDATED")
                  for name in UNKNOWN_SECTORS})
    gates["native_FZEEIN_propagation"]=propagation
    lfv=lfv_evaluate(point)
    if lfv["full_lfv_gate"]!="UNKNOWN" or lfv["global_zee_acceptance"] is not None:
        raise AssertionError("Partial LFV promoted to physical pass")
    fail=[k for k in ("neutrino_1D_3sigma","h125_mass_screen",
                      "cascade_on_shell") if gates[k]["status"]=="FAIL"]
    return {
        "point_id":entry["point_id"],"ordering":entry["ordering"],
        "status":"NATIVE_SOURCE_BOUND_PARTIAL_EVALUATION",
        "verdict":"FAIL_NECESSARY_SCREEN" if fail else "UNKNOWN_PENDING_GATES",
        "globally_allowed":False if fail else None,
        "failed_necessary_screens":fail,"gates":gates,
        "native_neutrino_screen":neutrino,
        "native_scalar_parameters":scalar,
        "LFV_status":"APPROXIMATION_ONLY",
        "LFV_retained_BR":{
            name:item["branching_ratio_included_topologies_only"]
            for name,item in lfv["processes"].items()},
        "LFV_retained_amplitudes":{
            name:{"norm_GeV_minus2":item["amplitude_norm_GeV_minus2"],
                  "topology":item["retained_topology_materiality"]}
            for name,item in lfv["processes"].items()},
        "LFV_full_gate":"UNKNOWN",
        "charged_photon_topology":lfv["charged_scalar_photon"],
        "card_sha256":binding,
        "lfv_point_sha256":lfv["input_sha256"],
        "source_sha256":{name:sha(p) for name,p in paths.items()},
        "source_note":"No FZEEIN input propagation proof, no full Zee LFV or global-constraint certification.",
    }
