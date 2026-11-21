"""Parser and physical-scope tests for pinned independent HDECAY SM control."""
import tempfile
from pathlib import Path
import unittest

from ZeeHiggs.hdecay_sm_reference import parse_table, parse_outputs


class HdecayParserTest(unittest.TestCase):
    def test_parse_ordered_sm_br_tables(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)
            (p/"br.sm1").write_text(
                "MHSM  BB TAU TAU MU MU SS CC TT\n"
                "125.250 0.58 0.062 0.0002 0.0002 0.029 0.0\n")
            (p/"br.sm2").write_text(
                "MHSM GG GAM GAM Z GAM WW ZZ WIDTH\n"
                "125.250 0.086 0.002 0.0016 0.213 0.026 0.0041\n")
            (p/"hdecay.in").write_text("HIGGS=5\n2HDM=0\nMHL = 125.250D0\n")
            (p/"run.stdout.log").write_text("OK\n")
            v=parse_outputs(p,125.25)
            self.assertAlmostEqual(v["width_GeV"],0.0041,delta=1e-12)
            self.assertAlmostEqual(v["BR"]["bb"],.58)
            self.assertAlmostEqual(v["partial_width_GeV"]["bb"],0.0041*.58)
            self.assertEqual(len(v["source_files_sha256"]),4)
            with self.assertRaisesRegex(ValueError,"Expected one"):
                parse_outputs(p,130.)

    def test_parse_fails_closed_for_nonfinite_negative_and_duplicate(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/"br.sm2"
            p.write_text("125.25 0.1 nan 0.0 0.5 0.3 0.004\n")
            with self.assertRaisesRegex(ValueError,"Expected one"):
                parse_table(p,("gg","gamgam","Zgam","WW","ZZ","total_width_GeV"),125.25)
            p.write_text("125.25 0.1 -0.01 0.0 0.5 0.3 0.004\n")
            with self.assertRaisesRegex(ValueError,"Negative"):
                parse_table(p,("gg","gamgam","Zgam","WW","ZZ","total_width_GeV"),125.25)
            p.write_text("125.25 0.1 0.01 0.0 0.5 0.3 0.004\n"*2)
            with self.assertRaisesRegex(ValueError,"Expected one"):
                parse_table(p,("gg","gamgam","Zgam","WW","ZZ","total_width_GeV"),125.25)


if __name__=="__main__":
    unittest.main()
