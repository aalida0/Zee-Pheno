"""Native Zee scalar spectrum -> source-locked BMU NLO -> HEPfit bridge.

Use the previous, validated HEPfit client in the separate ~/zee checkout.
No FlavorKit polynomial updates or physical exclusion from central rates.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import subprocess

from LFV.pilot import neutral_alpha
from Model.audit_neutrino_native import parse_entries

REPO_ROOT = Path(__file__).resolve().parents[1]
NATIVE_ROOT = REPO_ROOT / "Scan/reference/native_benchmark"
MATCHER_SHA = "0c2dfb32a7138b02fa4a52c163c63f942f313bb5d697545e138ef1915c200c4e"
HEPFIT_COMMIT = "322bd73a2dadd94ce94a6d71cb79d616e569120a"
SPHENO_SHA = "2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52"
MATCHING_KEYS = ("C7_LO", "C8_LO", "C4_NLO_raw", "C7_NLO_raw", "C8_NLO_raw")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def block_matrix(spectrum, name, size=3):
    values = {}
    active = False
    for line in spectrum.splitlines():
        pieces = line.split("#", 1)[0].split()
        if pieces and pieces[0].upper() == "BLOCK":
            active = len(pieces) > 1 and pieces[1].upper() == name.upper()
            continue
        if not active or not pieces:
            continue
        if len(pieces) == 3:
            i, j = int(pieces[0]), int(pieces[1])
            if (i, j) in values:
                raise ValueError(f"Duplicate {name} entry {i},{j}")
            values[i, j] = float(pieces[2].replace("D", "E"))
    expected = {(i, j) for i in range(1, size + 1) for j in range(1, size + 1)}
    if set(values) != expected or not all(math.isfinite(z) for z in values.values()):
        raise ValueError(f"Incomplete or nonfinite {name}")
    for i in range(1, size + 1):
        for j in range(1, size + 1):
            dot = sum(values[i, k] * values[j, k] for k in range(1, size + 1))
            if abs(dot - (1. if i == j else 0.)) > 1e-5:
                raise ValueError(f"Nonorthogonal {name}")
    return values


def native_scalar(order, root=NATIVE_ROOT):
    """Use pole charged masses and native charged rotations, never proposal masses."""
    if order not in ("NO", "IO"):
        raise ValueError("Ordering must be NO or IO")
    directory = root / f"s00017_{order}_benchmark"
    card, spectrum_path, native = (directory / x for x in
                                   ("input.in", "spectrum.slha", "native.json"))
    record = json.loads(native.read_text())
    if (record.get("production_binary_sha256") != SPHENO_SHA or
            record.get("spectrum_status") != "PASS" or
            record.get("neutrino_status") != "PASS"):
        raise ValueError("Unattested or invalid native production result")
    if sha(card) != record["input_sha256"]:
        raise ValueError("Card/native source mismatch")
    if "snapshots" not in record or "spectrum_snapshot_index" not in record:
        raise ValueError("Missing validated native snapshot")
    snapshot = record["snapshots"][record["spectrum_snapshot_index"]]
    vd, vu = float(snapshot["vd"]), float(snapshot["vu"])
    if not 0 < vd < vu or not all(map(math.isfinite, (vd, vu))):
        raise ValueError("Invalid s00017 EWSB")
    beta = math.atan2(vu, vd)
    tanb = vu / vd
    text = spectrum_path.read_text()
    masses = dict(parse_entries(text, "MASS"))
    minpar = dict(parse_entries(text, "MINPAR"))
    hm = dict(parse_entries(text, "HMIX"))
    # HMIX[10] is the scalar-sector betaH angle. It differs from
    # atan2(vu,vd) by ~4.8e-5 rad at this loop-matched benchmark.
    if abs(hm[10] - beta) > 2e-4 or abs(minpar[10] - tanb) > 2e-5:
        raise ValueError("Native beta mismatches SLHA")
    rot = block_matrix(text, "CHARGEMIX")
    # Directly measure the physical SU(2)-doublet norm of each charged row.
    # This avoids conflating betaH (the native scalar mixing angle) with
    # the Yukawa-running vev ratio after tadpole/loop iteration.
    weights = [sum(rot[row, j] ** 2 for j in (1, 2)) for row in (2, 3)]
    if any(w < 0 or w > 1 + 1e-6 for w in weights) or abs(sum(weights) - 1) > 1e-5:
        raise ValueError("Charged-state doublet projectors do not sum to 1")
    light, heavy = masses[37], masses[9000037]
    if not (light > 0 and heavy > light):
        raise ValueError("Bad pole charged masses")
    alpha = neutral_alpha(text)
    # HEPfit THDM bma is a host-only spectator for the *bridge* mode.
    # Map signed native hVV onto the conventional nonnegative sin(b-a).
    host_bma = math.asin(min(1., abs(math.cos(beta - alpha))))
    return {
        "ordering": order,
        "card_sha256": sha(card), "native_sha256": sha(native),
        "spectrum_sha256": sha(spectrum_path),
        "masses_GeV": {"H1": light, "H2": heavy, "h": masses[25],
                       "H": masses[35], "A": masses[36]},
        "tanbeta": tanb, "beta": beta,
        "native_hmix_beta_radians": hm[10],
        "native_hmix_vs_running_beta_radians": hm[10] - beta,
        "doublet_weights": weights,
        "soft_m12sq_GeV2": minpar[11],
        "hepfit_host_bma": host_bma,
        "legacy_flavorkit_BR": dict(parse_entries(text, "FlavorKitQFV"))[200],
    }


def backend(upstream, driver, hepfit_source, attestation=None):
    base = upstream / "research/full_constraints/bsgamma_hepfit"
    receipt_path = base / "execution/HEPFIT_EXECUTION.json"
    match_path = upstream / "scripts/zee_module5_bmu_nlo.py"
    source_lock_path = base / "SOURCE_LOCK.json"
    receipt, lock = json.loads(receipt_path.read_text()), json.loads(source_lock_path.read_text())
    if (receipt["HEPfit_commit"] != HEPFIT_COMMIT or lock["HEPfit_commit"] != HEPFIT_COMMIT
            or sha(match_path) != MATCHER_SHA
            or sha(base / "driver.cpp") != receipt["driver_sha256"]):
        raise ValueError("HEPfit or BMU runtime differs from source-attested validated execution")
    restored = None
    if attestation is not None:
        from Flavor.validate_hepfit import verify_attestation
        restored = verify_attestation(attestation,upstream,driver,hepfit_source)
    elif sha(driver) != receipt["binary_sha256"]:
        raise ValueError("Changed HEPfit binary requires a numerically validated build attestation")
    if not (hepfit_source / "THDM/tabs/bsgammatable.dat").is_file():
        raise FileNotFoundError("Pinned HEPfit source/tables not available")
    if (receipt.get("full_rate_status") != "B_FULL_RATE_BLOCKED"
            or receipt.get("physical_bsgamma") != "UNKNOWN"):
        raise ValueError("Unexpected historical rate-certification contract")
    match = runpy.run_path(str(match_path))["matching"]
    proof = {"upstream_receipt_sha256": sha(receipt_path),
             "matching_sha256": sha(match_path),
             "driver_binary_sha256": sha(driver),
             "driver_source_sha256": sha(base / "driver.cpp"),
             "HEPfit_commit": HEPFIT_COMMIT}
    if restored is not None:
        proof.update(restored_backend_attestation_sha256=sha(attestation),
                     restored_backend=restored)
    return receipt, match, proof


def ordinary(scalar, mt_msbar, matching, mu=160.):
    if not mt_msbar > 0 or not math.isfinite(mt_msbar):
        raise ValueError("Missing HEPfit running top mass")
    parts = [matching(scalar["masses_GeV"][label], scalar["tanbeta"], mt_msbar, mu=mu)
             for label in ("H1", "H2")]
    weights = scalar["doublet_weights"]
    combined = {k: sum(weights[i] * parts[i][k] for i in range(2))
                for k in MATCHING_KEYS}
    if not all(math.isfinite(x) for x in combined.values()):
        raise ValueError("Nonfinite NLO Wilson matching")
    np = [[0.] * 8 for _ in range(2)]
    np[0][6:8] = [combined["C7_LO"], combined["C8_LO"]]
    np[1][3] = combined["C4_NLO_raw"]
    np[1][6:8] = [combined["C7_NLO_raw"], combined["C8_NLO_raw"]]
    return np, combined


def parse_driver_output(text, matching_scale=160.):
    result = {}
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "UNKNOWN_INPUT":
            raise ValueError("HEPfit silently ignores an input: " + line)
        if parts[0] == "VALUE":
            if len(parts) != 3 or parts[1] in result:
                raise ValueError("Malformed or duplicate HEPfit VALUE")
            x = float(parts[2])
            if not math.isfinite(x):
                raise ValueError("Nonfinite HEPfit observable")
            result[parts[1]] = x
    needed = {"BR_NLO", "BR_stock_NNLO", "cut_GeV", "mt_MSbar_mut", "muw_GeV",
              "mu_b_GeV", "mu_c_GeV", "a_s_muw", "overall", "P_NLO", "N"}
    if not needed <= result.keys() or not result["BR_NLO"] > 0:
        raise ValueError("Incomplete executable HEPfit output")
    if abs(result["cut_GeV"] - 1.6) > 1e-12 or abs(result["muw_GeV"] - matching_scale) > 1e-10:
        raise ValueError("Wrong HEPfit photon cut or matching scale")
    if abs(result["BR_NLO"] - result["overall"] *
           (result["P_NLO"] + result["N"])) > 2e-15:
        raise ValueError("HEPfit NLO observable does not close")
    return result


def run_driver(label, mode, parameters, directory, driver, env):
    card = directory / (label + ".input")
    log = directory / (label + ".log")
    if card.exists() or log.exists():
        raise FileExistsError(label)
    if not all(math.isfinite(v) for v in parameters.values()):
        raise ValueError("Nonfinite HEPfit input")
    card.write_text("".join(f"{key} {value:.17g}\n" for key, value in sorted(parameters.items())))
    try:
        run = subprocess.run([str(driver), mode, str(card.resolve())],
                             env=env, text=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, timeout=120, check=False)
        output = run.stdout
        code = run.returncode
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        output += "\nTIMEOUT\n"
        code = None
    log.write_text(output)
    if code != 0:
        raise RuntimeError(f"HEPfit {label} failed; inspect {log}")
    return parse_driver_output(output), {"input_sha256": sha(card),
                                         "log_sha256": sha(log), "exit_code": code}


def injected_coeff_check(values, np):
    """Check actual *high-scale* order-tagged HEPfit bridge coefficients."""
    asmu = values["a_s_muw"]
    want = {
        "high_1_0_7": np[0][6], "high_1_0_8": np[0][7],
        "high_1_1_4": asmu * np[1][3],
        "high_1_1_7": asmu * (np[1][6] - 4 * np[1][3] / 9),
        "high_1_1_8": asmu * (np[1][7] - np[1][3] / 6),
    }
    for k, target in want.items():
        if k not in values or abs(values[k] - target) > 2e-12:
            raise ValueError("HEPfit BMU-to-effective-order bridge mismatch: " + k)
    return want


def calculate(upstream, driver, hepfit_source, compat_lib, root_lib, output_dir,
              attestation=None):
    """Runs native NO, IO, an SM reference and a zero-NP reference."""
    receipt, matcher, provenance = backend(upstream, driver, hepfit_source,attestation)
    scalars = [native_scalar(order) for order in ("NO", "IO")]
    sm_input = dict(receipt["runs"]["SM_reference"]["input"])
    # The archived baseline is a declared SM/hadronic scheme, not a
    # substitution of every running SPheno SMINPUTS value.
    if sm_input.get("muw") != 160. or sm_input.get("mut") != 160.:
        raise ValueError("Unexpected historical HEPfit scale")
    prototype = {k:v for k,v in receipt["runs"]["F13_heavy_0001330_central"]["input"].items()
                 if not k.startswith("NP")}
    output_dir.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, HEPFITTABS=str(hepfit_source.resolve()),
               LD_LIBRARY_PATH=f"{compat_lib.resolve()}:{root_lib.resolve()}")
    runs, sm_runs = [], {}
    for scalar in scalars:
        params = dict(prototype)
        masses = scalar["masses_GeV"]
        params.update(logtb=math.log10(scalar["tanbeta"]),
                      mHh2=masses["H"] ** 2, mA2=masses["A"] ** 2,
                      mHp2=masses["H1"] ** 2, m12_2=scalar["soft_m12sq_GeV2"],
                      mHl=masses["h"], bma=scalar["hepfit_host_bma"])
        zero = dict(params)
        for o in range(2):
            for i in range(8):
                zero[f"NP{o}_{i+1}"] = 0.
        name = scalar["ordering"]
        # Match each HEPfit SM control to the fitted point's Higgs mass.
        # Other SM/hadronic inputs deliberately retain the validated HEPfit prescription.
        current_sm_input = dict(sm_input, mHl=masses["h"])
        sm, sm_execution = run_driver(name + "_SM", "sm", current_sm_input,
                                      output_dir, driver, env)
        sm_runs[name] = {"BR_NLO": sm["BR_NLO"], "execution": sm_execution}
        no_np, zero_execution = run_driver(name + "_zero_NP", "bridge", zero,
                                           output_dir, driver, env)
        if abs(no_np["BR_NLO"] - sm["BR_NLO"]) > 2e-12:
            raise ValueError("Zero NP does not reproduce same-prescription SM")
        if abs(no_np["mt_MSbar_mut"] - sm["mt_MSbar_mut"]) > 1e-8:
            raise ValueError("Matching top mass is inconsistent")
        np, matching_coeffs = ordinary(scalar, no_np["mt_MSbar_mut"], matcher)
        p = dict(params)
        for o in range(2):
            for i in range(8):
                p[f"NP{o}_{i+1}"] = np[o][i]
        fitted, execution = run_driver(name + "_Zee_NLO", "bridge", p,
                                       output_dir, driver, env)
        injected_coeff_check(fitted, np)
        runs.append({
            "ordering": name, "scalar": scalar,
            "matching": matching_coeffs,
            "ordinary_NP_vector": np,
            "BR_NLO_Egamma_gt_1p6GeV": fitted["BR_NLO"],
            "BR_SM_same_scheme": sm["BR_NLO"],
            "delta_BR_same_scheme": fitted["BR_NLO"] - sm["BR_NLO"],
            "historical_FlavorKit_BR": scalar["legacy_flavorkit_BR"],
            "executions": {"zero_NP": zero_execution, "Zee_NLO": execution},
            "value_checks": {
                "zero_np_reproduces_sm": True,
                "high_scale_injected_coefficients": "PASS",
            },
            "note": "HEPfit BR_stock_NNLO is intentionally not used; NP matching NLO only.",
            "uncertainty_sigma": None,
            "physical_bsgamma_gate": "UNKNOWN",
            "full_model_allowed": None,
        })
    report = {
        "status": "EXECUTED_CENTRAL_NLO_ONLY_NO_PHYSICAL_EXCLUSION",
        "prediction_accuracy": "Real, top-only charged-Higgs NLO QCD NP matching; HEPfit NLO rate with native HQE/EM prescription.",
        "limits": "Missing model-specific primed, light quark, electroweak, resolved-photon and correlated uncertainty coverage.",
        "source": provenance, "SM_reference_runs": sm_runs,
        "matching_scale_GeV": 160., "photon_cut_GeV": 1.6,
        "SM_rates_same_HEPfit_scheme": {k: v["BR_NLO"] for k, v in sm_runs.items()},
        "points": runs, "theory_sigma": None,
        "physical_bsgamma": "UNKNOWN",
        "experimental_exclusion": None, "globally_allowed": None,
    }
    (output_dir / "bsgamma_nlo.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--upstream", type=Path, default=Path.home() / "zee",
                    help="Older source-locked ~/zee checkout, research/zee-full-constraints")
    ap.add_argument("--driver", type=Path, default=Path.home()/"zee-hepfit-build/driver")
    ap.add_argument("--hepfit-source", type=Path, default=Path.home()/
                    ("zee-hepfit-build/HEPfit-"+HEPFIT_COMMIT))
    ap.add_argument("--compat-lib-dir", type=Path, default=Path.home()/
                    "zee-hepfit-build/compat/usr/lib/x86_64-linux-gnu")
    ap.add_argument("--root-lib-dir", type=Path, default=Path.home() / "root/lib")
    ap.add_argument("--attestation",type=Path,
                    help="Fresh source/build/37-run equivalence receipt for a rebuilt binary")
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    result = calculate(args.upstream.resolve(), args.driver.resolve(),
                       args.hepfit_source.resolve(), args.compat_lib_dir.resolve(),
                       args.root_lib_dir.resolve(), args.output_dir.resolve(),
                       args.attestation.resolve() if args.attestation else None)
    print(json.dumps({
        "status": result["status"],
        "SM_rates": result["SM_rates_same_HEPfit_scheme"],
        "NO_IO_rates": {x["ordering"]: x["BR_NLO_Egamma_gt_1p6GeV"] for x in result["points"]},
        "physical_bsgamma_gate": "UNKNOWN",
        "result": str(args.output_dir / "bsgamma_nlo.json"),
    }, indent=2))


if __name__ == "__main__":
    main()
