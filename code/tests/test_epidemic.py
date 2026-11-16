import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import engine
from engine import ArticleCredibilityChecker as C, EpidemicDetectionEngine as E
from storage import Store

DEMOS = E.DEMO_ARTICLES


class TempStore(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.store = Store(os.path.join(self._tmp.name, "test.db"))

    def tearDown(self):
        self._tmp.cleanup()


class CredibilityHeuristic(unittest.TestCase):
    def test_fake_demo_articles_are_blocked(self):
        for key in ("fake_1", "fake_2"):
            self.assertIn(C._heuristic_check(DEMOS[key])["verdict"], ("FAKE", "LIKELY_FAKE"), key)

    def test_real_demo_articles_are_not_blocked(self):
        for key in ("covid_2020", "mpox_2022", "normal"):
            self.assertNotIn(C._heuristic_check(DEMOS[key])["verdict"], ("FAKE", "LIKELY_FAKE"), key)

    def test_score_is_bounded_and_flags_are_listed(self):
        r = C._heuristic_check("SHOCKING!!! MIRACLE CURE!!! SHARE BEFORE DELETED")
        self.assertTrue(0.0 <= r["credibility_score"] <= 1.0)
        self.assertTrue(r["flags"])


class KeywordFallback(unittest.TestCase):
    def test_disease_keywords(self):
        self.assertEqual(E._demo_analysis("Suspected Ebola in the district")["disease"], "Ebola")
        self.assertEqual(E._demo_analysis("New SARS-CoV-2 variant")["disease"], "COVID-19")
        self.assertEqual(E._demo_analysis("monkeypox cluster")["disease"], "Mpox")
        routine = E._demo_analysis("two cases of seasonal flu")
        self.assertEqual(routine["disease"], "Unknown")
        self.assertLess(routine["risk_score"], 0.3)


class Coercion(unittest.TestCase):
    def test_to_int_handles_messy_model_output(self):
        self.assertEqual(engine._to_int("15+"), 15)
        self.assertEqual(engine._to_int("1,200"), 1200)
        self.assertEqual(engine._to_int(3.7), 3)
        self.assertEqual(engine._to_int(None), 0)
        self.assertEqual(engine._to_int("about 40 cases"), 40)

    def test_to_unit_float_clamps(self):
        self.assertEqual(engine._to_unit_float(1.7), 1.0)
        self.assertEqual(engine._to_unit_float(-2), 0.0)
        self.assertEqual(engine._to_unit_float("0.4"), 0.4)
        self.assertEqual(engine._to_unit_float("high", 0.5), 0.5)


class ClaudeJson(unittest.TestCase):
    """_claude_json against a fake client: parsing and retry behaviour."""

    def _client(self, *replies):
        msgs = [SimpleNamespace(content=[SimpleNamespace(type="text", text=t)]) for t in replies]
        return SimpleNamespace(messages=SimpleNamespace(create=mock.Mock(side_effect=msgs)))

    def _call(self, client):
        with mock.patch.object(engine, "anthropic_client", client), \
             mock.patch.object(engine, "_append_claude_log"), mock.patch.object(engine.time, "sleep"):
            return engine._claude_json("prompt", "TEST")

    def test_braces_inside_strings_do_not_break_parsing(self):
        out = self._call(self._client('Here you go: {"summary": "a } tricky { one", "cases": 3} thanks'))
        self.assertEqual(out, {"summary": "a } tricky { one", "cases": 3})

    def test_retries_after_empty_or_non_json_reply(self):
        client = self._client("", "no json here", '{"ok": true}')
        self.assertEqual(self._call(client), {"ok": True})
        self.assertEqual(client.messages.create.call_count, 3)

    def test_raises_when_every_attempt_fails(self):
        with self.assertRaises(ValueError):
            self._call(self._client("nope", "nope", "nope"))


class EpidemicPipeline(TempStore):
    def run_pipeline(self, article):
        """Run with Claude disconnected, so only the offline paths are used."""
        with mock.patch.object(engine, "anthropic_client", None):
            return engine.run_epidemic_pipeline("RPT_T", article, store=self.store)

    def test_real_report_is_scored_and_saved_as_keyword_fallback(self):
        out = self.run_pipeline(DEMOS["covid_2020"])
        self.assertTrue(out["ok"])
        self.assertFalse(out["blocked"])
        self.assertTrue(out["is_demo"])                      # no Claude -> never labelled "Real"
        self.assertIn("MSH|", out["hl7"])
        rec = self.store.get_analysis(out["analysis_id"])
        self.assertEqual(rec["type"], "EPIDEMIC")
        self.assertEqual(rec["is_demo"], 1)
        self.assertEqual(rec["location"], "")                # "Unknown" is never sent to the geocoder

    def test_fake_report_is_blocked_and_not_saved(self):
        out = self.run_pipeline(DEMOS["fake_1"])
        self.assertTrue(out["blocked"])
        self.assertNotIn("hl7", out)
        self.assertEqual(self.store.count_analyses(), 0)

    def test_unverifiable_dampening_never_raises_a_score(self):
        cred = {"verdict": "UNVERIFIABLE", "credibility_score": 0.6, "flags": [], "reasoning": "", "method": "TEST"}
        for original in (0.15, 0.40, 0.82):
            analysis = dict(E._demo_analysis("routine"), risk_score=original, source="CLAUDE_API")
            with mock.patch.object(C, "check", return_value=dict(cred)), \
                 mock.patch.object(E, "analyze_article", return_value=analysis):
                out = engine.run_epidemic_pipeline("RPT_D", "text", store=self.store)
            self.assertLessEqual(out["risk_score"], original + 1e-9, original)
            self.assertIsNotNone(out["dampened_note"])
        # a serious signal keeps a floor of 50%
        self.assertGreaterEqual(out["risk_score"], 0.5)

    def test_location_is_saved_for_the_map(self):
        analysis = dict(E._demo_analysis("mpox"), location="Lagos, Nigeria", source="CLAUDE_API")
        cred = {"verdict": "LIKELY_REAL", "credibility_score": 0.9, "flags": [], "reasoning": "", "method": "TEST"}
        with mock.patch.object(C, "check", return_value=cred), \
             mock.patch.object(E, "analyze_article", return_value=analysis):
            out = engine.run_epidemic_pipeline("RPT_L", "text", store=self.store)
        rec = self.store.get_analysis(out["analysis_id"])
        self.assertEqual(rec["location"], "Lagos, Nigeria")
        self.assertEqual(rec["is_demo"], 0)


if __name__ == "__main__":
    unittest.main()
