import unittest

from engine import ConfidenceScorer as S, HL7Generator as H


class Confidence(unittest.TestCase):
    def test_blast_confidence_is_bounded(self):
        for args in [(100, 0.0, 100, "RESFINDER"), (0, 1.0, 0, "MANUAL_REVIEW"),
                     (250, 1e-400, 500, "LOCAL_DATABASE"), (-5, 5.0, -10, "nonsense")]:
            self.assertTrue(0.0 <= S.amr_confidence(*args) <= 1.0, args)

    def test_weaker_classification_source_scores_lower(self):
        strong = S.amr_confidence(99, 1e-150, 100, "LOCAL_DATABASE")
        motif = S.amr_confidence(99, 1e-150, 100, "SEQUENCE_MOTIF")
        manual = S.amr_confidence(99, 1e-150, 100, "MANUAL_REVIEW")
        self.assertGreater(strong, motif)
        self.assertGreater(motif, manual)

    def test_amrfinder_exact_full_length_match_is_certain(self):
        self.assertAlmostEqual(S.amrfinder_confidence(100, 100, "ALLELEX"), 1.0)

    def test_amrfinder_partial_or_divergent_hits_score_lower(self):
        exact = S.amrfinder_confidence(100, 100, "EXACTX")
        blast = S.amrfinder_confidence(96, 100, "BLASTX")
        partial = S.amrfinder_confidence(99, 50, "PARTIALX")
        self.assertGreater(exact, blast)
        self.assertGreater(blast, partial)
        self.assertTrue(0.0 <= partial <= 1.0)

    def test_risk_category_thresholds(self):
        self.assertEqual(S.risk_category(0.85)[0], "CRITICAL")
        self.assertEqual(S.risk_category(0.84)[0], "HIGH")
        self.assertEqual(S.risk_category(0.70)[0], "HIGH")
        self.assertEqual(S.risk_category(0.50)[0], "MODERATE")
        self.assertEqual(S.risk_category(0.30)[0], "LOW")
        self.assertEqual(S.risk_category(0.29)[0], "MINIMAL")


class HL7(unittest.TestCase):
    def test_amr_report_structure(self):
        msg = H.amr_report("ISO_1", "blaNDM-1", 0.97, 100.0, 1e-120, "CRITICAL")
        segs = msg.split("\n")
        self.assertEqual([s.split("|")[0] for s in segs], ["MSH", "PID", "OBR"] + ["OBX"] * 5)
        self.assertIn("ORU^R01", segs[0])
        self.assertTrue(segs[0].endswith("|2.5"))
        self.assertIn("blaNDM-1", segs[3])
        self.assertEqual(segs[3].split("|")[8], "H")       # abnormal flag set for a critical gene

    def test_delimiters_in_values_are_escaped(self):
        msg = H.amr_report("A|B^C&D\\E", "gene|x", 0.5, 90.0, 1e-5, "LOW")
        pid = msg.split("\n")[1]
        self.assertIn("A\\F\\B\\S\\C\\T\\D\\E\\E", pid)
        self.assertNotIn("A|B", msg)
        self.assertIn("gene\\F\\x", msg)

    def test_evalue_formatting(self):
        self.assertEqual(H._ef(None), "N/A")               # AMRFinderPlus calls have no E-value
        self.assertEqual(H._ef(0), "< 1e-300")
        self.assertEqual(H._ef(1e-320), "< 1e-300")
        self.assertEqual(H._ef(2.5e-30), "2.50e-30")
        self.assertIn("|N/A|", H.amr_report("S", "g", 0.9, 99.0, None, "HIGH"))

    def test_low_confidence_low_concern_is_not_flagged(self):
        msg = H.amr_report("S", "g", 0.4, 80.0, 1e-5, "MODERATE")
        self.assertEqual(msg.split("\n")[3].split("|")[8], "N")

    def test_epidemic_alert(self):
        hi = H.epidemic_alert("RPT_1", "Mpox", 0.82, 12, "rising", "Lagos, Nigeria").split("\n")
        lo = H.epidemic_alert("RPT_2", "Influenza", 0.15, 2, "stable", "Unknown").split("\n")
        self.assertEqual(len(hi), 8)
        self.assertEqual(hi[3].split("|")[8], "H")
        self.assertEqual(lo[3].split("|")[8], "N")
        self.assertIn("|12|count|", hi[5])
        self.assertIn("Lagos, Nigeria", hi[7])


if __name__ == "__main__":
    unittest.main()
