"""Regression checks for the Zee SPheno input gate (synthetic, not physics)."""
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from Model.check_native import MATRICES, antisym, blocks, fix_card, mat


def sample():
    lines = ["Block MINPAR\n 1 0.2\n"]
    for name in MATRICES:
        lines.append(f"Block {name}\n")
        for i in range(1, 4):
            for j in range(1, 4):
                value = (1e-4 if (i, j) == (1, 2) else
                         -1e-4 if (i, j) == (2, 1) else 0.) if name == "FZEEIN" else (
                         3e-4 if i == j else 0.) if name == "YEE2IN" else (
                         1e-4 if i == j else 0.) if name == "IMYEE2IN" else 0.
                lines.append(f" {i} {j} {value:.17e}\n")
    return "".join(lines)


class NativeContracts(unittest.TestCase):
    def test_fourth_field_complex_matrix_and_no_input_mutation(self):
        source = sample()
        fixed = fix_card(source)
        self.assertEqual(source.count("#"), 0)
        self.assertEqual(fixed.count("# generated entry"), 37)
        for name in MATRICES:
            section = fixed.split(f"Block {name}\n", 1)[1].split("Block", 1)[0]
            for row in section.splitlines():
                self.assertEqual(row.split()[3], "#")
        self.assertEqual(antisym(mat(blocks(fixed), "FZEEIN")), 0.)

    def test_probe_and_invalid_antisimmetry(self):
        source = sample()
        new = fix_card(source, {("FZEEIN", (1, 2)): 2e-4,
                                ("FZEEIN", (2, 1)): -2e-4})
        self.assertAlmostEqual(mat(blocks(new), "FZEEIN")[1, 2].real, 2e-4)
        with self.assertRaises(ValueError):
            fix_card(source, {("FZEEIN", (1, 2)): 2e-4})
        with self.assertRaises(ValueError):
            fix_card(source.replace("Block IMFZEEIN", "Block LOST_IMFZEEIN"))

    def test_synthetic_native_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            card = p / "source.in"
            card.write_text(sample())
            fake = p / "mock_spheno.py"
            fake.write_text(
                '#!/usr/bin/env python3\n'
                'import pathlib,sys\n'
                'x=pathlib.Path(sys.argv[1]).read_text()\n'
                'for a,b in (("FZEEIN","FZEE"),("IMFZEEIN","IMFZEE"),'
                '("YEE2IN","YEE2"),("IMYEE2IN","IMYEE2")):\n'
                ' x=x.replace("Block "+a+"\\n","Block "+b+"\\n")\n'
                'pathlib.Path(sys.argv[2]).write_text(x)\n')
            fake.chmod(0o755)
            attestation = p / "attestation.json"
            attestation.write_text(json.dumps({
                "binary_sha256": hashlib.sha256(fake.read_bytes()).hexdigest(),
                "force_complete_rebuild": True, "build_exit_code": 0}))
            result = subprocess.run(
                [sys.executable, "-m", "Model.check_native",
                 "--input", str(card), "--output-dir", str(p / "run"),
                 "--exe", str(fake), "--attestation", str(attestation)],
                text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            saved = json.loads((p / "run" / "native-input-report.json").read_text())
            self.assertEqual(saved["status"], "NATIVE_INPUT_GATE_PASS")
            self.assertEqual(card.read_text(), sample())


if __name__ == "__main__":
    unittest.main()
