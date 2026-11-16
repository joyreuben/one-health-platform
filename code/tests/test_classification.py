import unittest
from unittest import mock

import engine
from engine import AMRFinderClient, AMRValidationEngine

# AMRFinderPlus v4 column names (nucleotide run)
TSV_V4 = (
    "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\t"
    "Class\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t"
    "% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\t"
    "HMM accession\tHMM description\n"
    "NA\ts0\t1\t813\t+\tblaNDM-1\tsubclass B1 metallo-beta-lactamase NDM-1\tcore\tAMR\tAMR\t"
    "BETA-LACTAM\tCARBAPENEM\tALLELEX\t270\t270\t100.00\t100.00\t270\tWP_004201164.1\t"
    "subclass B1 metallo-beta-lactamase NDM-1\tNA\tNA\n"
    "NA\ts1\t10\t600\t-\ttet(A)\ttetracycline efflux MFS transporter Tet(A)\tcore\tAMR\tAMR\t"
    "TETRACYCLINE\tTETRACYCLINE\tPARTIALX\t197\t399\t49.37\t99.49\t197\tWP_000804064.1\t"
    "tetracycline efflux MFS transporter Tet(A)\tNA\tNA\n"
)

# AMRFinderPlus v3 column names (protein run)
TSV_V3 = (
    "Protein identifier\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\t"
    "Method\tTarget length\tReference sequence length\t% Coverage of reference sequence\t"
    "% Identity to reference sequence\tAlignment length\tAccession of closest sequence\t"
    "Name of closest sequence\tHMM id\tHMM description\n"
    "prot1\tmcr-1.1\tphosphoethanolamine--lipid A transferase MCR-1.1\tcore\tAMR\tAMR\tCOLISTIN\tCOLISTIN\t"
    "EXACTP\t541\t541\t100.00\t100.00\t541\tWP_049589868.1\tphosphoethanolamine--lipid A transferase MCR-1.1\tNA\tNA\n"
)


class AMRFinderParsing(unittest.TestCase):
    def test_v4_nucleotide_output(self):
        rows = AMRFinderClient.parse_tsv(TSV_V4)
        self.assertEqual(len(rows), 2)
        ndm = rows[0]
        self.assertEqual(ndm["seq_id"], "s0")           # "NA" protein id is skipped
        self.assertEqual(ndm["symbol"], "blaNDM-1")
        self.assertEqual(ndm["subclass"], "CARBAPENEM")
        self.assertEqual(ndm["method"], "ALLELEX")
        self.assertEqual(ndm["identity"], 100.0)
        self.assertEqual(ndm["accession"], "WP_004201164.1")
        self.assertAlmostEqual(rows[1]["coverage"], 49.37)

    def test_v3_protein_output(self):
        rows = AMRFinderClient.parse_tsv(TSV_V3)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seq_id"], "prot1")
        self.assertEqual(rows[0]["symbol"], "mcr-1.1")
        self.assertEqual(rows[0]["class"], "COLISTIN")

    def test_header_only_or_empty_means_no_findings(self):
        self.assertEqual(AMRFinderClient.parse_tsv(""), [])
        self.assertEqual(AMRFinderClient.parse_tsv(TSV_V4.splitlines()[0] + "\n"), [])


class AMRFinderInterpretation(unittest.TestCase):
    def test_last_resort_drug_classes_are_critical(self):
        self.assertEqual(AMRFinderClient.concern("BETA-LACTAM", "CARBAPENEM"), "CRITICAL")
        self.assertEqual(AMRFinderClient.concern("COLISTIN", "COLISTIN"), "CRITICAL")
        self.assertEqual(AMRFinderClient.concern("GLYCOPEPTIDE", "VANCOMYCIN"), "CRITICAL")
        self.assertEqual(AMRFinderClient.concern("BETA-LACTAM", "BETA-LACTAM"), "HIGH")
        self.assertEqual(AMRFinderClient.concern("MACROLIDE", "MACROLIDE"), "MODERATE")
        self.assertEqual(AMRFinderClient.concern("", ""), "MODERATE")

    def test_method_weight_ignores_search_type_suffix(self):
        self.assertEqual(AMRFinderClient.method_weight("ALLELEX"), 1.0)
        self.assertEqual(AMRFinderClient.method_weight("EXACTP"), 1.0)
        self.assertGreater(AMRFinderClient.method_weight("BLASTX"), AMRFinderClient.method_weight("PARTIALX"))
        self.assertGreater(AMRFinderClient.method_weight("PARTIALX"), AMRFinderClient.method_weight("INTERNAL_STOP"))
        self.assertEqual(AMRFinderClient.method_weight("HMM"), 0.5)

    def test_result_picks_best_gene_and_lists_all(self):
        findings = AMRFinderClient.parse_tsv(TSV_V4)
        with mock.patch.object(AMRFinderClient, "info",
                               return_value={"ok": True, "software": "4.0.0", "database": "2026-01-01.1"}):
            res = engine._amrfinder_result("ISO_1", "ATGC" * 30, findings)
        self.assertTrue(res["ok"])
        self.assertEqual(res["display_gene"], "blaNDM-1")
        self.assertEqual(res["gene_info"]["clinical_concern"], "CRITICAL")
        self.assertEqual(res["gene_info"]["source"], "AMRFINDERPLUS")
        self.assertEqual(len(res["findings"]), 2)
        self.assertGreater(res["confidence"], 0.95)
        self.assertIn("blaNDM-1", res["hl7"])
        self.assertIsNone(res["best"]["e_value"])


class TitleMatching(unittest.TestCase):
    """BLAST-title fallback classifier. Titles are normalised first, which
    turns '-' into a space, so patterns must tolerate that."""

    def match(self, title):
        hit = AMRValidationEngine._match_title(AMRValidationEngine._normalise(title))
        return hit["marker_found"] if hit else None

    def test_hyphenated_gene_names(self):
        self.assertEqual(self.match("Escherichia coli plasmid mcr-1 gene, complete cds"), "mcr")
        self.assertEqual(self.match("Klebsiella pneumoniae blaCTX-M-15 gene"), "ctx-m")
        self.assertEqual(self.match("Klebsiella pneumoniae KPC-2 carbapenemase"), "kpc")
        self.assertEqual(self.match("NDM-1 metallo-beta-lactamase"), "ndm")

    def test_descriptive_names(self):
        self.assertEqual(self.match("D-Ala-D-Ala ligase VanA [Enterococcus faecium]"), "vanA")
        self.assertEqual(self.match("penicillin-binding protein 2a"), "mecA")
        self.assertEqual(self.match("tetracycline resistance protein TetA"), "tetA")

    def test_unrelated_title_does_not_match(self):
        self.assertIsNone(self.match("hypothetical protein [Bacillus subtilis]"))
        self.assertIsNone(self.match("16S ribosomal RNA gene, partial sequence"))


if __name__ == "__main__":
    unittest.main()
