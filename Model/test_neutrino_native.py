"""Tiny fail-closed tests; a source symbol occurrence never certifies a loop mass."""
import tempfile
import unittest
from pathlib import Path
from Model.audit_neutrino_native import parse_entries, audit


class AuditTests(unittest.TestCase):
    def test_flag_and_mass_sections(self):
        text = ("Block SPHENOINPUT\n 55 0 # tree\nBlock MASS\n"
                " 12 0.0\n 14 0.0\n 16 0.0\nDECAY 25 0.004\n")
        self.assertEqual(parse_entries(text, "SPHENOINPUT"), [(55, 0.)])
        self.assertEqual(dict(parse_entries(text, "MASS"))[12], 0.)

    def test_explicitly_disabled_loop_masses_cannot_pass(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "LoopMasses_ZeeTypeII.f90").write_text("subroutine SelfEnergyFv()\nend\n")
            card = root / "a.in"
            card.write_text("Block SPHENOINPUT\n 55 0\n")
            spec = root / "a.spc"
            spec.write_text("Block MASS\n 12 0.0\n 14 0.0\n 16 0.0\n")
            result = audit(root, card, spec)
            self.assertEqual(result["status"], "LOOP_MASSES_DISABLED")
            self.assertFalse(result["native_one_loop_neutrino_mass_validated"])
            self.assertTrue(result["source_lines_to_inspect"])


if __name__ == "__main__":
    unittest.main()
