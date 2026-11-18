"""Write a source-bound NO/IO manifest for the existing two native pilots.

These are NOT physics-accepted benchmarks: their PDG25 masses are
approximately 484 GeV, and the output is only a status/validation audit.
This command does not invent quark running/Yukawa masses.
"""
import argparse
import json
from pathlib import Path

from Model.audit_neutrino_native import sha


def prepare(work,output,quark_loop=None,quark_yukawa=None,scheme=None):
    quark=(quark_loop,quark_yukawa,scheme)
    if any(x is not None for x in quark) and not all(x is not None for x in quark):
        raise ValueError("Quark mass dictionaries and renormalization scheme must be supplied together")
    if all(x is not None for x in quark) and (
        not quark_loop.is_file() or not quark_yukawa.is_file()
        or not isinstance(scheme,str) or not scheme.strip()):
        raise ValueError("Missing quark dictionary files or scheme")
    rows=[]
    for ordering in ("NO","IO"):
        card=(work/f"step1_{ordering}_proposal.in").resolve()
        native=(work/f"step1_{ordering}_native/result.json").resolve()
        spectrum=(work/f"step1_{ordering}_native/output.spc").resolve()
        if not card.is_file() or not native.is_file() or not spectrum.is_file():
            raise FileNotFoundError(f"{ordering}: native card/diagnostic/spectrum not found")
        row={"point_id":f"{ordering}_reference_diagnostic",
             "ordering":ordering,
             "card":str(card),"native":str(native),"spectrum":str(spectrum)}
        if all(x is not None for x in quark):
            row.update(quark_loop=str(quark_loop.resolve()),
                       quark_yukawa=str(quark_yukawa.resolve()),
                       mass_yukawa_scheme=scheme)
        rows.append(row)
    manifest={"schema_version":1,"points":rows}
    with output.open("x") as f:
        json.dump(manifest,f,indent=2)
        f.write("\n")
    return {"status":"REFERENCE_MANIFEST_WRITTEN_NOT_PHYSICS_ACCEPTED",
            "n_points":len(rows),"manifest_sha256":sha(output),
            "LFV_scope":("SCALAR_AND_SELECTED_PHOTON_BZ"
                         if all(x is not None for x in quark)
                         else "SCALAR_ONE_LOOP_ONLY"),
            "native_F_input_propagation":"UNKNOWN_UNLESS_CARD_BOUND_PROBE_INCLUDED",
            "known_h125_mass_warning":"NATIVE_PDG25_AROUND_484_GEV"}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--quark-loop",type=Path)
    ap.add_argument("--quark-yukawa",type=Path)
    ap.add_argument("--mass-yukawa-scheme",type=str)
    args=ap.parse_args()
    print(json.dumps(prepare(args.work,args.output,args.quark_loop,
                             args.quark_yukawa,args.mass_yukawa_scheme),indent=2))

if __name__=="__main__":
    main()
