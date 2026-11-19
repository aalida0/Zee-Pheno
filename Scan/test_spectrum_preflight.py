"""Synthetic SPheno early-termination gate; not a physics validation."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Model.audit_neutrino_native import sha
from Scan.spectrum_preflight import preflight


class SpectrumPreflightTests(unittest.TestCase):
    def test_exit_zero_without_spectrum_is_failure_and_complete_spectrum_is_pending(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            cards=root/"cards"
            cards.mkdir()
            generated=[]
            for name in ("s00000_NO_f00.in","s00001_IO_f00.in"):
                card=cards/name
                card.write_text("dummy input\n")
                generated.append({"ordering":"NO" if "NO" in name else "IO",
                                  "status":"NATIVE_REVALIDATION_REQUIRED",
                                  "card":name,"card_sha256":sha(card)})
            catalog=root/"proposals.json"
            catalog.write_text(json.dumps({
                "schema_version":2,
                "status":"SCALAR_NEUTRINO_PROPOSALS_NATIVE_UNKNOWN",
                "outputs":[
                    {"scalar_id":"s00000","families":[generated[0]]},
                    {"scalar_id":"s00001","families":[generated[1]]}]}))
            executable=root/"fake-spheno"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "import pathlib,sys\n"
                "if 's00000' in sys.argv[1]:\n"
                " print('CalculateMhh: negative mass squared; TerminateProgram')\n"
                "else:\n"
                " pathlib.Path(sys.argv[2]).write_text('Block MASS\\n"
                " 25 125.0\\n 35 500.0\\n 36 550.0\\n"
                " 37 300.0\\n 9000037 650.0\\n')\n")
            executable.chmod(0o755)
            with patch("Scan.spectrum_preflight.BINARY_SHA",sha(executable)):
                result=preflight(catalog,executable,root/"out")
            self.assertEqual(result["n_cases"],2)
            self.assertEqual(result["n_failed"],1)
            self.assertEqual(result["n_spectrum_present"],1)
            self.assertEqual(result["cases"][0]["exit_code"],0)
            self.assertIn("negative mass squared",
                          (root/"out/s00000_NO_f00/run.log").read_text())
            self.assertEqual(result["cases"][1]["scalar_pole_masses_GeV"]["25"],125.0)
            self.assertEqual((cards/generated[0]["card"]).read_text(),"dummy input\n")
            with patch("Scan.spectrum_preflight.BINARY_SHA",sha(executable)):
                with self.assertRaises(FileExistsError):
                    preflight(catalog,executable,root/"out")


if __name__=="__main__":
    unittest.main()
