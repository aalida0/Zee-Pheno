"""Versioned manifest runner for native neutrino -> LFV -> tri-state gates.

Produces JSON and CSV, never globally accepts candidates. Relative paths
are resolved against the manifest's own directory. Source integrity
failures are recorded per candidate and never count as constraint passes.
"""
import argparse
import csv
import json
from pathlib import Path

from Model.audit_neutrino_native import sha
from Scan.native_gates import evaluate_entry,check_entry

def evaluate_manifest(manifest_path):
    manifest=json.loads(manifest_path.read_text())
    if not isinstance(manifest,dict) or set(manifest)!={"schema_version","points"}:
        raise ValueError("Manifest requires schema_version and points")
    if type(manifest["schema_version"]) is not int or manifest["schema_version"]!=1:
        raise ValueError("Unsupported schema_version")
    entries=manifest["points"]
    if not isinstance(entries,list) or not 1<=len(entries)<=1000:
        raise ValueError("Manifest needs 1..1000 points")
    for entry in entries:
        check_entry(entry)
    ids=[e["point_id"] for e in entries]
    if len(ids)!=len(set(ids)):
        raise ValueError("Duplicate point_id")
    reports=[]
    for entry in entries:
        try:
            reports.append(evaluate_entry(entry,manifest_path.parent))
        except (ValueError,TypeError,KeyError,IndexError,FileNotFoundError,
                ArithmeticError,AssertionError) as exc:
            # Failing source integrity is not a physics exclusion.
            reports.append({
                "point_id":entry["point_id"],"ordering":entry["ordering"],
                "status":"INVALID_SOURCE_OR_EVALUATION",
                "verdict":"NOT_EVALUATED_SOURCE_ERROR",
                "globally_allowed":None,
                "gates":{},"LFV_status":"NOT_EVALUATED",
                "source_error_type":type(exc).__name__,
                "source_error":str(exc),
            })
    return {
        "schema_version":1,
        "status":"PARTIAL_NATIVE_EVIDENCE_NO_GLOBAL_ACCEPTANCE",
        "manifest_sha256":sha(manifest_path),
        "n_entries":len(entries),
        "n_native_1D_pass":sum(
            r.get("gates",{}).get("neutrino_1D_3sigma",{}).get("status")=="PASS"
            for r in reports),
        "n_fail_necessary":sum(
            r["verdict"]=="FAIL_NECESSARY_SCREEN" for r in reports),
        "n_pending":sum(
            r["verdict"]=="UNKNOWN_PENDING_GATES" for r in reports),
        "n_invalid_source":sum(
            r["verdict"]=="NOT_EVALUATED_SOURCE_ERROR" for r in reports),
        "n_globally_allowed":None,
        "LFV_full_gate":"UNKNOWN",
        "reports":reports,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    result=evaluate_manifest(args.manifest)
    args.output_dir.mkdir(parents=True,exist_ok=False)
    with (args.output_dir/"points.json").open("x") as fh:
        json.dump(result,fh,indent=2,allow_nan=False)
        fh.write("\n")
    keys=("point_id","ordering","verdict","native_1D_nufit",
          "h125_mass","cascade_on_shell","LFV_full",
          "BR_mu_e_gamma_retained","BR_tau_e_gamma_retained",
          "BR_tau_mu_gamma_retained")
    with (args.output_dir/"points.csv").open("x",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=keys)
        writer.writeheader()
        for row in result["reports"]:
            gates=row.get("gates",{})
            part=row.get("LFV_retained_BR",{})
            writer.writerow({
                "point_id":row["point_id"],"ordering":row["ordering"],
                "verdict":row["verdict"],
                "native_1D_nufit":gates.get("neutrino_1D_3sigma",{}).get("status","UNKNOWN"),
                "h125_mass":gates.get("h125_mass_screen",{}).get("status","UNKNOWN"),
                "cascade_on_shell":gates.get("cascade_on_shell",{}).get("status","UNKNOWN"),
                "LFV_full":"UNKNOWN",
                "BR_mu_e_gamma_retained":part.get("mu_to_e_gamma"),
                "BR_tau_e_gamma_retained":part.get("tau_to_e_gamma"),
                "BR_tau_mu_gamma_retained":part.get("tau_to_mu_gamma"),
            })
    print(json.dumps({key:value for key,value in result.items()
                      if key!="reports"},indent=2))

if __name__=="__main__":
    main()
