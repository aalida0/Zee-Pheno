"""Source-bound audit of the native SPheno Zee Higgs input.

This is NOT an implementation of HiggsBounds, HiggsSignals, or Higgs widths.
It authenticates and checks the native SPheno artifact before handing it to
the official HiggsTools interface. Never interpret a missing loop quantity as 0.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "Scan/reference/native_benchmark"
NEUTRAL = (25, 35, 36)
CHARGED = (37, 9000037)
PARTICLES = (6,) + NEUTRAL + CHARGED
SPHENO_ORIGINAL_SHA256 = "2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def numeric(s):
    f = float(s.replace("D", "E").replace("d", "e"))
    if not math.isfinite(f):
        raise ValueError("Nonfinite SLHA input")
    return f


def read_slha(text):
    """Conservative reader for SLHA BLOCK and DECAY/BR rows (not a SPheno replacement)."""
    blocks, decays = {}, {}
    active_block, active_decay = None, None
    for lineno, raw in enumerate(text.splitlines(), 1):
        body = raw.split("#", 1)[0].strip()
        if not body:
            continue
        fields = body.split()
        if fields[0].upper() == "BLOCK":
            if len(fields) < 2:
                raise ValueError(f"Malformed BLOCK at line {lineno}")
            active_block, active_decay = fields[1].upper(), None
            if active_block in blocks:
                raise ValueError("Duplicate SLHA BLOCK: " + active_block)
            blocks[active_block] = []
            continue
        if fields[0].upper() == "DECAY":
            if len(fields) != 3:
                raise ValueError("Malformed DECAY header")
            pid = int(fields[1])
            if pid in decays:
                raise ValueError("Duplicate DECAY header: " + str(pid))
            width = numeric(fields[2])
            if width < 0:
                raise ValueError("Negative DECAY width")
            decays[pid] = {"total_width_GeV": width, "channels": []}
            active_decay, active_block = pid, None
            continue
        if active_block is not None:
            blocks[active_block].append((fields, raw))
            continue
        if active_decay is not None:
            if len(fields) < 3:
                raise ValueError("Malformed decay channel")
            br, n = numeric(fields[0]), int(fields[1])
            ids = tuple(int(x) for x in fields[2:])
            if n != len(ids) or br < 0 or br > 1 + 1e-7:
                raise ValueError("Invalid branching fraction or daughter count")
            decays[active_decay]["channels"].append({"br": br, "pdgs": ids})
    return blocks, decays


def simple_block(blocks, name):
    if name not in blocks:
        raise ValueError("Missing SLHA block " + name)
    entries = {}
    for row, _ in blocks[name]:
        if len(row) < 2:
            raise ValueError("Invalid scalar block row: " + name)
        key = tuple(int(s) for s in row[:-1])
        if key in entries:
            raise ValueError("Duplicate scalar parameter: " + name)
        entries[key] = numeric(row[-1])
    return entries


def scalar_boson_coupling(blocks, pdg, v1, v2):
    if "HIGGSCOUPLINGSBOSONS" not in blocks:
        raise ValueError("Missing HiggsCouplingsBosons")
    matches = []
    for row, _ in blocks["HIGGSCOUPLINGSBOSONS"]:
        if len(row) < 5:
            raise ValueError("Malformed HiggsCouplingsBosons")
        count = int(row[1])
        if count < 2 or len(row) != 2 + count:
            raise ValueError("Bad boson-coupling leg count")
        legs = [int(s) for s in row[2:]]
        if pdg == legs[0] and len(legs) == 3 and set(legs[1:]) == {v1, v2}:
            matches.append(numeric(row[0]))
    if len(matches) != 1:
        raise ValueError("Missing or duplicated requested boson coupling")
    return matches[0]


def scalar_fermion_coupling(blocks, pdg, fermion):
    matches = []
    for row, _ in blocks.get("HIGGSCOUPLINGSFERMIONS", []):
        if len(row) < 6:
            raise ValueError("Malformed HiggsCouplingsFermions")
        r, im = numeric(row[0]), numeric(row[1])
        if int(row[2]) != 3:
            raise ValueError("Unexpected HiggsCouplingsFermions multiplicity")
        if (int(row[3]), int(row[4]), int(row[5])) == (pdg, fermion, fermion):
            matches.append((r, im))
    if len(matches) != 1:
        raise ValueError("Missing or duplicated neutral-fermion coupling")
    return matches[0]


def audit(order, native_root=NATIVE):
    if order not in ("NO", "IO"):
        raise ValueError("Only fixed native NO/IO controls are accepted")
    path = Path(native_root) / f"s00017_{order}_benchmark"
    input_path, native_path, spectrum_path = (path / f for f in
                                             ("input.in", "native.json", "spectrum.slha"))
    record = json.loads(native_path.read_text())
    if (sha256(input_path) != record["input_sha256"]
            or record.get("production_binary_sha256") != SPHENO_ORIGINAL_SHA256
            or record.get("spectrum_status") != "PASS"
            or record.get("neutrino_status") != "PASS"):
        raise ValueError("Unattested/invalid native SPheno spectrum")
    text = spectrum_path.read_text()
    blocks, decays = read_slha(text)
    masses = simple_block(blocks, "MASS")
    mass = {k: masses[(k,)] for k in PARTICLES if (k,) in masses}
    if len(mass) != len(PARTICLES):
        raise ValueError("Missing top or physical Zee scalar mass")
    if not (120 < mass[25] < 130 and
            0 < mass[25] < mass[35] < mass[36] and
            mass[37] < mass[9000037]):
        raise ValueError("Unexpected native scalar hierarchy")
    for pid in PARTICLES:
        if pid not in decays or decays[pid]["total_width_GeV"] <= 0:
            raise ValueError("Missing/zero total width: " + str(pid))
        channels = decays[pid]["channels"]
        if not channels:
            raise ValueError("Missing channels: " + str(pid))
        s = sum(x["br"] for x in channels)
        if abs(1 - s) > 1e-3:
            raise ValueError(f"DECAY BR sum fails for {pid}: {s}")
        decays[pid]["branching_sum"] = s
    if not {"HIGGSCOUPLINGSBOSONS", "HIGGSCOUPLINGSFERMIONS", "HIGGSLHC13"} <= blocks.keys():
        raise ValueError("Missing HiggsTools coupling or SPheno production blocks")
    scalar_couplings = {}
    for pid in NEUTRAL:
        fermions = {str(f): scalar_fermion_coupling(blocks, pid, f)
                    for f in (4, 5, 6, 13, 15)}
        gauge = {label:scalar_boson_coupling(blocks, pid, *legs)
                 for label, legs in (("WW", (24,24)),("ZZ",(23,23)),("gamgam",(22,22)),
                                     ("gg",(21,21)))}
        scalar_couplings[str(pid)] = {"fermion_re_im": fermions, "boson": gauge}
    production = {}
    for row, _ in blocks["HIGGSLHC13"]:
        if len(row) < 3:
            raise ValueError("Malformed HiggsLHC13 row")
        process, pid = int(row[0]), int(row[1])
        if process == 1 and pid in NEUTRAL:
            if pid in production:
                raise ValueError("Duplicate 13 TeV ggF production")
            production[pid] = numeric(row[2])
    if set(production) != set(NEUTRAL) or min(production.values()) <= 0:
        raise ValueError("Missing valid neutral ggF cross section")
    # Explicit SPheno comment: h->Zgamma is NOT calculated, not zero.
    for rows in (blocks.get("EFFHIGGSCOUPLINGS",[]),):
        for row, raw in rows:
            if len(row) >= 4 and row[0] in ("25","35","36") and row[1:3] == ["22","23"]:
                if not ("not yet calculated" in raw.lower() and numeric(row[3]) == 0):
                    raise ValueError("Unexpected provenance of Zgamma effective coupling")
                break
        else:
            raise ValueError("Missing explicit SPheno Zgamma not-calculated warning")
    h2 = decays[9000037]
    cascade = [ch["br"] for ch in h2["channels"] if ch["pdgs"] in ((37,23),(23,37))]
    if len(cascade) != 1 or not (.01 < cascade[0] < 1):
        raise ValueError("Missing/invalid charged-cascade channel")
    notes = [
        "SPheno explicitly says h/H/A -> Zgamma loop amplitudes NOT CALCULATED; zero is a placeholder",
        "Neutral Higgs 13 TeV ggF entries are SusHi-1.5.0-rescaled approximations, not certified production theory",
        "Charged H1/H2 t(b)H production cross sections not in the native SLHA; do not set to zero",
        "Check HiggsTools dataset topology coverage, narrow-width and production/decay normalization",
        "HiggsSignals effective-coupling normalization recommended for h125: unmodified SPheno h total-width precision not validated",
        "Full charged-scalar width, loop correction and HiggsTools prediction-contract validation pending",
    ]
    return {
        "status":"SOURCE_ATTESTED_NATIVE_SPHENO_HIGGS_INPUTS_AUDITED",
        "point_id":"s00017_"+order,
        "spheno_version":"4.0.5", "sarah_version":"4.15.4",
        "source":{"card_sha256":sha256(input_path),"native_sha256":sha256(native_path),
                  "spectrum_sha256":sha256(spectrum_path),
                  "spheno_binary_sha256":record["production_binary_sha256"]},
        "spectrum_path":str(spectrum_path),
        "neutral_pdgs":list(NEUTRAL),"charged_pdgs":list(CHARGED),
        "masses_GeV":{str(k):v for k,v in mass.items()},
        "total_widths_GeV":{str(k):v["total_width_GeV"] for k,v in decays.items()
                            if k in PARTICLES},
        "branching_sum":{str(k):decays[k]["branching_sum"] for k in PARTICLES},
        "charged_cascade_br":cascade[0],
        "neutral_ggF_LHC13_pb":{str(k):v for k,v in production.items()},
        "neutral_couplings":scalar_couplings,
        "input_limitations":notes,
        "run_higgstools_ready_for_diagnostic":True,
        "widths_precision_accepted":False,
        "charged_production_accepted":False,
        "higgsbounds_full_model_coverage_validated":False,
        "higgssignals_full_model_coverage_validated":False,
        "physical_higgs":"UNKNOWN","globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result={o:audit(o) for o in ("NO","IO")}
    with args.output.open("x") as dest:
        json.dump(result,dest,indent=2,allow_nan=False)
        dest.write("\n")
    for o,row in result.items():
        print(o,"audited",
              "mh=",row["masses_GeV"]["25"],
              "gamma_h=",row["total_widths_GeV"]["25"],
              "cascade=",row["charged_cascade_br"],
              "physical_higgs=UNKNOWN")


if __name__=="__main__":
    main()
