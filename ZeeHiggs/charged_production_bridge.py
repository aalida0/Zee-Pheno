"""Source-bound LO charged-Higgs rates for the official HiggsTools adapter.

No production physics is computed here. Rates must originate from a separate,
model-matched generator integration for the EXACT SPheno point. In particular,
13.6 TeV predictions cannot be silently entered as 13 TeV HiggsBounds inputs;
missing modes are not zero. The allowed process is 5FS pp -> tb H_i^+-,
summed over both charges, at LHC13, with nonoverlapping production modes.
"""
import json
import math
import re
from pathlib import Path

from ZeeHiggs.spheno_audit import CHARGED, sha256

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def validate_charged_rates(payload, source):
    """Return HiggsTools cxn keys and an explicitly *partial* coverage ledger.

    This validates data identity/format, not an independent MadGraph amplitude
    or QCD/theory-error calculation. The caller must source-match the actual
    model and generator run cards before treating any limits as physical.
    """
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Charged production receipt requires schema_version=1")
    if payload.get("point_id") != source["point_id"]:
        raise ValueError("Charged production receipt point ID mismatch")
    identity = payload.get("source")
    if not isinstance(identity, dict):
        raise ValueError("Missing SPheno source identity")
    for receipt_key, report_key in (
        ("card_sha256", "card_sha256"),
        ("spectrum_sha256", "spectrum_sha256"),
        ("spheno_binary_sha256", "spheno_binary_sha256"),
    ):
        if identity.get(receipt_key) != source["source"][report_key]:
            raise ValueError("Charged production receipt source mismatch: "+receipt_key)
    rows = payload.get("rates")
    if not isinstance(rows, list) or len(rows) != len(CHARGED):
        raise ValueError("Exactly two per-eigenstate charged rates required")
    data = {}
    evidence = []
    for row in rows:
        if not isinstance(row, dict) or type(row.get("pdg")) is not int:
            raise ValueError("Malformed charged rate row")
        pdg = row["pdg"]
        if pdg not in CHARGED or pdg in data:
            raise ValueError("Unexpected or duplicated charged particle")
        required = {
            "mode": "Hpmtb",
            "collider": "LHC13",
            "sqrt_s_TeV": 13,
            "flavor_scheme": "5FS",
            "perturbative_order": "LO",
            "charge_inclusive": True,
        }
        for k, expected in required.items():
            if type(expected) is bool:
                valid = row.get(k) is expected
            else:
                valid = row.get(k) == expected
            if not valid:
                raise ValueError("Incompatible charged rate: "+k)
        for k in ("sigma_pb", "integration_error_pb", "mass_GeV"):
            v = row.get(k)
            if type(v) not in (int,float) or not math.isfinite(v) or v < 0:
                raise ValueError("Invalid or nonfinite production "+k)
        mass = source["masses_GeV"][str(pdg)]
        if not math.isclose(row["mass_GeV"], mass, rel_tol=1e-7, abs_tol=1e-4):
            raise ValueError("Production and native charged scalar mass mismatch")
        if row["sigma_pb"] == 0 and not row.get("explicit_zero_reason"):
            raise ValueError("Zero rate requires evidence, never a missing-rate default")
        if row["sigma_pb"] > 0 and row["integration_error_pb"] >= row["sigma_pb"]:
            raise ValueError("Production numerical integration error too large")
        if not isinstance(row.get("generator"),str) or not row["generator"].strip():
            raise ValueError("Production generator not identified")
        for digest in ("param_card_sha256", "run_card_sha256",
                       "process_card_sha256", "generator_result_sha256"):
            if not isinstance(row.get(digest), str) or not SHA256_RE.fullmatch(row[digest]):
                raise ValueError("Missing or invalid generator receipt hash: "+digest)
        # These hashes identify external evidence but do not, by themselves,
        # prove that the rate was taken from that evidence.
        data[pdg] = row["sigma_pb"]
        evidence.append({
            "pdg": pdg, "mass_GeV":mass, "sigma_pb":row["sigma_pb"],
            "integration_error_pb":row["integration_error_pb"],
            "scheme":"13 TeV pp, LO 5FS, charge-inclusive Hpmtb",
            "generator": row["generator"],
            "generator_receipt_hashes_present":True,
            "generator_native_result_independently_checked":False,
            "pdf_scale_theory_uncertainty":"NOT_QUANTIFIED",
        })
    if set(data) != set(CHARGED):
        raise ValueError("Both charged scalar mass eigenstates required")
    return {
        "status":"CHARGED_LO_13TEV_INPUT_FORMAT_SOURCE_BOUND_NOT_PHYSICS_CERTIFIED",
        "higgstools_cross_sections":{f"cxn_{pdg}_Hpmtb_LHC13":data[pdg]
                                      for pdg in CHARGED},
        "evidence":evidence,
        "charged_13TeV_Hpmtb_input_present":True,
        "all_energy_and_topology_production_complete":False,
        "generator_amplitudes_and_raw_rates_independently_verified":False,
        "full_HiggsBounds_model_exclusion":None,
        "globally_allowed":None,
    }


def load_charged_rates(path, source):
    """Check file identity and retain the exact receipt file digest."""
    p = Path(path)
    result = validate_charged_rates(json.loads(p.read_text()),source)
    return {**result, "receipt_path":str(p.resolve()),
            "receipt_sha256":sha256(p)}


def insert_charged_rates(data, ids, validated):
    """Inject *only* Hpmtb; leave native top BRs and all other modes intact."""
    if validated.get("status") != "CHARGED_LO_13TEV_INPUT_FORMAT_SOURCE_BOUND_NOT_PHYSICS_CERTIFIED":
        raise ValueError("Unvalidated charged production receipt")
    merged = dict(data)
    if set(ids) < set(CHARGED):
        raise ValueError("Missing charged HiggsTools identifiers")
    for pdg in CHARGED:
        source_key=f"cxn_{pdg}_Hpmtb_LHC13"
        target_key=f"cxn_{ids[pdg]}_Hpmtb_LHC13"
        if source_key not in validated["higgstools_cross_sections"]:
            raise ValueError("Missing charged production rate")
        if target_key in merged:
            raise ValueError("Refusing to overwrite existing charged production")
        merged[target_key]=validated["higgstools_cross_sections"][source_key]
    return merged
