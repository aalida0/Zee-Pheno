"""Input mutation controls for the one-off attested SPheno loop probe."""
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
from Model.probe_neutrino_loop import enable_loop_masses, probe
from Model.audit_neutrino_native import parse_entries, sha


class LoopProbeTests(unittest.TestCase):
    def test_two_repeated_switches_only(self):
        original = ("Block SPHENOINPUT\n 55 0 # tree\n 11 0\n"
                    "Block FZEEIN\n 1 2 3e-4 # entry\n"
                    "Block SPHENOINPUT\n 55 0.0D+00 # tree\n 57 1 # LFV\n")
        updated = enable_loop_masses(original)
        self.assertEqual([v for k,v in parse_entries(updated,"SPHENOINPUT") if k==55], [1.,1.])
        self.assertIn(" 1 2 3e-4 # entry", updated)
        self.assertIn(" 57 1 # LFV", updated)
        changed = [(a, b) for a, b in zip(original.splitlines(), updated.splitlines()) if a != b]
        self.assertEqual(changed, [(" 55 0 # tree", " 55 1 # tree"),
                                   (" 55 0.0D+00 # tree", " 55 1 # tree")])

    def test_rejects_any_other_situation(self):
        for value in ("1", "2"):
            with self.assertRaises(ValueError):
                enable_loop_masses(f"Block SPHENOINPUT\n 55 {value}\nBlock SPHENOINPUT\n 55 0\n")
        with self.assertRaises(ValueError):
            enable_loop_masses("Block SPHENOINPUT\n 55 0\n")

    def test_one_card_native_mock_leaves_baseline_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            card, exe = root / "baseline.in", root / "mock.py"
            source = "Block SPHENOINPUT\n 55 0\nBlock SPHENOINPUT\n 55 0\n"
            card.write_text(source)
            exe.write_text("#!/usr/bin/env python3\nimport sys,pathlib\n"
                           "pathlib.Path(sys.argv[2]).write_text('Block MASS\\n 12 1.1e-11\\n 14 2.2e-11\\n 16 5.0e-11\\n')\n")
            exe.chmod(0o755)
            with patch("Model.probe_neutrino_loop.BASELINE_SHA", sha(card)), \
                 patch("Model.probe_neutrino_loop.BINARY_SHA", sha(exe)):
                result = probe(card, exe, root / "new")
            self.assertEqual(result["status"], "NEUTRINO_POLE_MASS_NONZERO_UNVALIDATED")
            self.assertFalse(result["native_one_loop_neutrino_mass_validated"])
            self.assertEqual(card.read_text(), source)
            self.assertTrue((root / "new" / "result.json").exists())


if __name__ == "__main__":
    unittest.main()
