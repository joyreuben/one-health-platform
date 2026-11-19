import unittest
from unittest import mock

import engine
from engine import (AMRFinderClient, HL7Generator, predict_phenotype, phenotype_profile,
                    run_amr_genome, validate_genome)


def finding(symbol, klass, subclass, subtype="AMR", method="EXACTX", type_="AMR", seq_id="c1"):
    return {"seq_id": seq_id, "symbol": symbol, "name": symbol, "scope": "core", "type": type_,
            "subtype": subtype, "class": klass, "subclass": subclass, "method": method,
            "coverage": 100.0, "identity": 100.0, "accession": "WP_0", "ref_name": symbol,
            "start": "1", "stop": "900", "strand": "+"}


KPC = finding("blaKPC-2", "BETA-LACTAM", "CARBAPENEM", method="ALLELEX")
SHV12 = finding("blaSHV-12", "BETA-LACTAM", "CEFIDEROCOL/CEPHALOSPORIN", method="ALLELEX")
OMPK35 = finding("ompK35_E42RfsTer47", "BETA-LACTAM", "CARBAPENEM", subtype="POINT_DISRUPT", method="POINTX")
GYRA = finding("gyrA_S83I", "QUINOLONE", "QUINOLONE", subtype="POINT", method="POINTX")
OQXB = finding("oqxB5", "NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE",
               "NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE")
AAC = finding("aac(6')-Ib-cr5", "AMINOGLYCOSIDE/QUINOLONE", "AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN")
OXA_STOP = finding("blaOXA", "BETA-LACTAM", "CARBAPENEM", method="INTERNAL_STOP")


class PhenotypeRules(unittest.TestCase):
    def test_carbapenemase_predicts_carbapenem_and_cephalosporin_resistance(self):
        self.assertEqual(predict_phenotype("meropenem", [KPC]), (1, ["blaKPC-2"]))
        self.assertEqual(predict_phenotype("ceftriaxone", [KPC])[0], 1)

    def test_esbl_alone_is_not_carbapenem_resistance(self):
        self.assertEqual(predict_phenotype("meropenem", [SHV12])[0], 0)
        self.assertEqual(predict_phenotype("ceftazidime", [SHV12])[0], 1)

    def test_porin_disruption_needs_an_esbl(self):
        self.assertEqual(predict_phenotype("ertapenem", [OMPK35])[0], 0)
        p, why = predict_phenotype("ertapenem", [OMPK35, SHV12])
        self.assertEqual(p, 1)
        self.assertEqual(why, ["blaSHV-12", "ompK35_E42RfsTer47"])

    def test_porin_mutation_does_not_count_for_cephalosporins(self):
        self.assertEqual(predict_phenotype("cefepime", [OMPK35])[0], 0)

    def test_intrinsic_allele_variants_are_ignored(self):
        # Regression: v1 of the benchmark matched "oqxB" exactly and missed "oqxB5".
        self.assertEqual(predict_phenotype("ciprofloxacin", [OQXB])[0], 0)
        self.assertEqual(predict_phenotype("ciprofloxacin", [OQXB, GYRA]), (1, ["gyrA_S83I"]))

    def test_internal_stop_is_ignored(self):
        self.assertEqual(predict_phenotype("imipenem", [OXA_STOP])[0], 0)

    def test_subclass_tokens(self):
        self.assertEqual(predict_phenotype("amikacin", [AAC])[0], 1)
        self.assertEqual(predict_phenotype("tobramycin", [AAC])[0], 1)
        self.assertEqual(predict_phenotype("gentamicin", [AAC])[0], 0)

    def test_profile_covers_every_drug(self):
        prof = phenotype_profile([KPC])
        self.assertEqual({p["drug"] for p in prof}, set(engine.PHENOTYPE_DRUG_RULES))


class GenomeValidation(unittest.TestCase):
    def test_accepts_multi_contig_dna_above_single_sequence_limit(self):
        big = {"contig1": "ACGT" * 150_000, "contig2": "acgtn\nACGT " * 1000}
        ok, err, cleaned = validate_genome(big)
        self.assertTrue(ok, err)
        self.assertGreater(sum(map(len, cleaned.values())), engine.SequenceValidator.MAX_LENGTH)
        self.assertNotIn(" ", cleaned["contig2"])

    def test_rejects_protein(self):
        ok, err, _ = validate_genome({"p": "MKLEFPQ" * 300})
        self.assertFalse(ok)
        self.assertIn("non-DNA", err)

    def test_rejects_too_small_and_empty(self):
        self.assertFalse(validate_genome({"c": "ACGT"})[0])
        self.assertFalse(validate_genome({})[0])

    def test_rejects_too_large(self):
        with mock.patch.object(engine, "GENOME_MAX_BP", 1_000):
            self.assertFalse(validate_genome({"c": "A" * 2_000})[0])


class GenomePipeline(unittest.TestCase):
    INFO = {"ok": True, "software": "4.2.7", "database": "2026-10-29.1"}

    def _run(self, findings_by_contig, organism):
        with mock.patch.object(AMRFinderClient, "is_available", return_value=True), \
             mock.patch.object(AMRFinderClient, "info", return_value=self.INFO), \
             mock.patch.object(AMRFinderClient, "run", return_value=findings_by_contig) as run:
            res = run_amr_genome("ISO1", {"contigA": "ACGT" * 500, "contigB": "TTGA" * 500}, organism)
        return res, run

    def test_kpn_genome_reports_findings_phenotype_and_hl7(self):
        res, run = self._run({"c1": [dict(KPC)], "c2": [dict(GYRA)]}, "Klebsiella_pneumoniae")
        self.assertTrue(res["ok"])
        self.assertEqual(run.call_args.kwargs["organism"], "Klebsiella_pneumoniae")
        self.assertEqual(sorted(run.call_args.args[0]), ["c1", "c2"])   # contigs renamed, sent in one call
        self.assertEqual({f["seq_id"] for f in res["findings"]}, {"contigA", "contigB"})  # names restored
        self.assertEqual(res["contigs"], 2)
        self.assertEqual(res["genome_bp"], 4000)
        drugs = {p["drug"]: p["predicted_resistant"] for p in res["phenotype"]}
        self.assertTrue(drugs["meropenem"])
        self.assertTrue(drugs["ciprofloxacin"])
        self.assertFalse(drugs["colistin"])
        self.assertEqual(len(res["sequence_sha256"]), 64)
        self.assertIn("AMR_DETERMINANT", res["hl7"])
        self.assertIn("PREDICTED_SUSCEPTIBILITY^meropenem", res["hl7"])

    def test_no_phenotype_for_unvalidated_organism(self):
        res, _ = self._run({"c1": [dict(KPC)]}, "Escherichia")
        self.assertTrue(res["ok"])
        self.assertIsNone(res["phenotype"])

    def test_negative_genome(self):
        res, _ = self._run({"c1": [], "c2": []}, None)
        self.assertFalse(res["ok"])
        self.assertTrue(res["negative"])

    def test_unknown_organism_rejected_before_running(self):
        with mock.patch.object(AMRFinderClient, "is_available", return_value=True), \
             mock.patch.object(AMRFinderClient, "run") as run:
            res = run_amr_genome("ISO1", {"c": "ACGT" * 500}, "Not_a_species")
        self.assertFalse(res["ok"])
        run.assert_not_called()

    def test_needs_amrfinder(self):
        with mock.patch.object(AMRFinderClient, "is_available", return_value=False):
            res = run_amr_genome("ISO1", {"c": "ACGT" * 500}, None)
        self.assertFalse(res["ok"])
        self.assertIn("AMRFinderPlus", res["error"])


class OrganismFlag(unittest.TestCase):
    def test_organism_is_passed_and_validated(self):
        with mock.patch.object(AMRFinderClient, "_exec", return_value=(0, "", "")) as ex:
            AMRFinderClient.run({"s0": "ACGT" * 30}, protein=False, organism="Klebsiella_pneumoniae")
        self.assertIn("-O Klebsiella_pneumoniae", ex.call_args.args[0])
        with mock.patch.object(AMRFinderClient, "_exec", return_value=(0, "", "")) as ex:
            AMRFinderClient.run({"s0": "ACGT" * 30}, protein=False)
        self.assertNotIn("-O", ex.call_args.args[0])
        with self.assertRaises(ValueError):
            AMRFinderClient.run({"s0": "ACGT" * 30}, protein=False, organism="x; rm -rf /")


class GenomeHL7(unittest.TestCase):
    def test_escapes_and_numbers_segments(self):
        msg = HL7Generator.amr_genome_report("S|1", "Klebsiella pneumoniae", [KPC, GYRA],
                                             [{"drug": "meropenem", "predicted_resistant": True,
                                               "determinants": ["blaKPC-2"]}])
        lines = msg.split("\n")
        self.assertTrue(lines[0].startswith("MSH|"))
        self.assertIn("S\\F\\1", lines[1])
        obx = [l for l in lines if l.startswith("OBX|")]
        self.assertEqual([l.split("|")[1] for l in obx], ["1", "2", "3", "4"])
        self.assertIn("|blaKPC-2^CARBAPENEM^ALLELEX^100.00|", obx[1])
        odd = dict(KPC, symbol="gene^x|y")
        line = HL7Generator.amr_genome_report("S", None, [odd]).split("\n")[-1]
        self.assertIn("|gene\\S\\x\\F\\y^CARBAPENEM^", line)   # data escaped, separators kept


if __name__ == "__main__":
    unittest.main()
