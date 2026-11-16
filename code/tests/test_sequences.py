import unittest

from engine import SequenceValidator as V

DNA = "ATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGC"


class ValidateAndClean(unittest.TestCase):
    def test_fasta_header_is_removed(self):
        ok, err, seq = V.validate_and_clean(f">sample1 some description\n{DNA}")
        self.assertTrue(ok, err)
        self.assertEqual(seq, DNA)

    def test_windows_line_endings_tabs_and_spaces_are_stripped(self):
        ok, err, seq = V.validate_and_clean(f">s\r\n{DNA[:30]}\r\n{DNA[30:40]} \t{DNA[40:]}\r\n")
        self.assertTrue(ok, err)
        self.assertEqual(seq, DNA)

    def test_lowercase_is_uppercased(self):
        ok, _, seq = V.validate_and_clean(DNA.lower())
        self.assertTrue(ok)
        self.assertEqual(seq, DNA)

    def test_trailing_protein_stop_symbol_is_dropped(self):
        ok, _, seq = V.validate_and_clean("MKKLLPLAVLAALLSGCQSTPEQ*")
        self.assertTrue(ok)
        self.assertFalse(seq.endswith("*"))

    def test_empty_and_header_only_are_rejected(self):
        self.assertFalse(V.validate_and_clean("")[0])
        self.assertFalse(V.validate_and_clean("   \n ")[0])
        self.assertFalse(V.validate_and_clean(">only a header")[0])

    def test_too_short_is_rejected(self):
        ok, err, _ = V.validate_and_clean("ATGC")
        self.assertFalse(ok)
        self.assertIn("short", err.lower())

    def test_invalid_characters_are_reported(self):
        ok, err, _ = V.validate_and_clean(DNA + "12#")
        self.assertFalse(ok)
        self.assertIn("Invalid characters", err)


class SequenceType(unittest.TestCase):
    def test_dna_is_not_protein(self):
        self.assertFalse(V.is_protein(DNA))
        self.assertFalse(V.is_protein("ATGCNNNNRYSWKM" * 3))   # IUPAC ambiguity codes are still DNA

    def test_protein_is_detected(self):
        self.assertTrue(V.is_protein("MELPNIMHPVAKLSTALAAALMLSGCMPGEIRPTIGQQ"))

    def test_likely_coding_needs_start_and_stop_codon(self):
        cds = "ATG" + "GCT" * 40 + "TAA"
        self.assertTrue(V.is_likely_coding(cds))
        self.assertFalse(V.is_likely_coding("GCT" * 40 + "TAA"))   # no start codon
        self.assertFalse(V.is_likely_coding("ATGGCTTAA"))          # too short


if __name__ == "__main__":
    unittest.main()
