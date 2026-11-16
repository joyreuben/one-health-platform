import os
import tempfile
import unittest
from unittest import mock

import engine
import feeds
from storage import Store

WHO_PAYLOAD = {"value": [
    {"Id": "guid-1", "DonId": "2026-DON626", "UrlName": "2026-DON626",
     "Title": "Middle East respiratory syndrome (MERS) – Qatar", "UseOverrideTitle": False,
     "PublicationDate": "2026-11-13T13:48:47Z",
     "Summary": "On 29 October 2026, Qatar notified WHO of two cases.",
     "Overview": "<div><p paraid='1'>On 29 October 2026, Qatar notified WHO of two cases.</p>"
                 "<p>Both patients were hospitalised &amp; isolated.</p><script>alert(1)</script></div>",
     "Epidemiology": "<p>MERS is caused by a coronavirus.</p>"},
    {"Id": "guid-2", "UrlName": "2026-DON625", "Title": "Ignored title", "UseOverrideTitle": True,
     "OverrideTitle": "Cholera – Region X", "PublicationDate": "2026-11-01T09:00:00Z",
     "Summary": "", "Overview": "<p>" + "Long text. " * 2000 + "</p>", "Epidemiology": None},
    {"Title": "No id at all"},
]}


class Feeds(unittest.TestCase):
    def test_html_to_text(self):
        t = feeds.html_to_text("<p>One&nbsp;&amp; two</p><p>Three<br>Four</p><style>p{}</style>")
        self.assertEqual(t, "One & two\nThree\nFour")
        self.assertEqual(feeds.html_to_text(""), "")
        self.assertEqual(feeds.html_to_text(None), "")

    def test_parse_who_don(self):
        items = feeds.parse_who_don(WHO_PAYLOAD)
        self.assertEqual(len(items), 2)                       # the item without any id is skipped
        a, b = items
        self.assertEqual(a["id"], "2026-DON626")
        self.assertEqual(a["url"], "https://www.who.int/emergencies/disease-outbreak-news/item/2026-DON626")
        self.assertIn("Published 2026-11-13", a["text"])
        self.assertIn("hospitalised & isolated", a["text"])
        self.assertNotIn("<p", a["text"])
        self.assertNotIn("alert(1)", a["text"])
        self.assertFalse(a["truncated"])
        self.assertEqual(b["id"], "2026-DON625")
        self.assertEqual(b["title"], "Cholera – Region X")    # override title wins
        self.assertTrue(b["truncated"])
        self.assertLessEqual(len(b["text"]), feeds.MAX_TEXT_CHARS)

    def test_parse_tolerates_empty_payload(self):
        self.assertEqual(feeds.parse_who_don({}), [])
        self.assertEqual(feeds.parse_who_don(None), [])


class TempStore(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.path = os.path.join(self._tmp.name, "nested", "test.db")
        self.store = Store(self.path)

    def tearDown(self):
        self._tmp.cleanup()


class Storage(TempStore):
    def test_analysis_round_trip(self):
        self.assertTrue(self.store.get_logs().empty)
        i = self.store.log_analysis("ISO_1", "GENOMIC", "blaNDM-1", 0.973, "AMRFINDERPLUS",
                                    details={"gene_info": {"class": "Beta-Lactam"}}, hl7="MSH|x")
        j = self.store.log_analysis("RPT_1", "EPIDEMIC", "Mpox", 0.7, "KEYWORD_FALLBACK",
                                    is_demo=True, location="Lagos")
        df = self.store.get_logs()
        self.assertEqual(list(df["ID"]), [i, j])
        self.assertEqual(df.iloc[0]["Confidence"], "97.3%")
        self.assertEqual(df.iloc[0]["Status"], "✅ Real")
        self.assertEqual(df.iloc[1]["Status"], "⚠️ DEMO")
        self.assertEqual(df.iloc[1]["Location"], "Lagos")
        rec = self.store.get_analysis(i)
        self.assertEqual(rec["hl7"], "MSH|x")
        self.assertEqual(rec["details"]["gene_info"]["class"], "Beta-Lactam")
        self.assertIsNone(self.store.get_analysis(9999))

    def test_data_survives_reopening_the_database(self):
        self.store.log_analysis("ISO_1", "GENOMIC", "mcr-1.1", 1.0, "AMRFINDERPLUS")
        self.assertEqual(Store(self.path).count_analyses(), 1)

    def test_settings(self):
        self.assertEqual(self.store.get_setting("missing", "dflt"), "dflt")
        self.store.set_setting("k", "1")
        self.store.set_setting("k", "0")
        self.assertEqual(self.store.get_setting("k"), "0")

    def test_feed_items_are_deduplicated_and_tracked(self):
        items = feeds.parse_who_don(WHO_PAYLOAD)
        self.assertEqual(self.store.upsert_feed_items(items), 2)
        self.assertEqual(self.store.upsert_feed_items(items), 0)      # nothing new the second time
        listed = self.store.feed_items()
        self.assertEqual([i["id"] for i in listed], ["2026-DON626", "2026-DON625"])   # newest first
        self.store.set_feed_status("2026-DON626", "analyzed", analysis_id=7)
        self.assertEqual([i["id"] for i in self.store.feed_items(status="new")], ["2026-DON625"])

    def test_clearing_analyses_requeues_feed_items(self):
        self.store.upsert_feed_items(feeds.parse_who_don(WHO_PAYLOAD))
        aid = self.store.log_analysis("2026-DON626", "EPIDEMIC", "MERS", 0.6, "CLAUDE_API")
        self.store.set_feed_status("2026-DON626", "analyzed", analysis_id=aid)
        self.assertEqual(self.store.clear_analyses(), 1)
        self.assertEqual(len(self.store.feed_items(status="new")), 2)


class FeedSync(TempStore):
    def test_sync_stores_new_items_and_skips_analysis_without_claude(self):
        with mock.patch.object(feeds, "fetch_who_don", return_value=feeds.parse_who_don(WHO_PAYLOAD)), \
             mock.patch.object(engine, "anthropic_client", None):
            first = engine.sync_feed(self.store)
            second = engine.sync_feed(self.store)
        self.assertEqual((first["fetched"], first["new"], first["analyzed"], first["error"]), (2, 2, 0, ""))
        self.assertEqual(second["new"], 0)
        self.assertTrue(self.store.get_setting("feed_last_fetch"))
        self.assertEqual(self.store.count_analyses(), 0)

    def test_sync_reports_network_errors_instead_of_raising(self):
        with mock.patch.object(feeds, "fetch_who_don", side_effect=RuntimeError("network down")):
            out = engine.sync_feed(self.store)
        self.assertIn("network down", out["error"])
        self.assertIn("network down", self.store.get_setting("feed_last_error"))

    def test_auto_analysis_respects_the_setting_and_the_cap(self):
        many = [{"id": f"DON{i:02d}", "source": "WHO DON", "title": f"Report {i}", "url": "",
                 "published_at": f"2026-11-{i + 1:02d}T00:00:00Z", "text": engine.EpidemicDetectionEngine.DEMO_ARTICLES["mpox_2022"]}
                for i in range(8)]
        cred = {"verdict": "LIKELY_REAL", "credibility_score": 0.9, "flags": [], "reasoning": "", "method": "TEST"}
        with mock.patch.object(feeds, "fetch_who_don", return_value=many), \
             mock.patch.object(engine, "anthropic_client", object()), \
             mock.patch.object(engine.ArticleCredibilityChecker, "check", return_value=cred), \
             mock.patch.object(engine.EpidemicDetectionEngine, "analyze_article",
                               return_value=dict(engine.EpidemicDetectionEngine._demo_analysis("mpox"),
                                                 source="CLAUDE_API")):
            self.store.set_setting("feed_auto_analyze", "0")
            off = engine.sync_feed(self.store)
            self.store.set_setting("feed_auto_analyze", "1")
            on = engine.sync_feed(self.store)
        self.assertEqual(off["analyzed"], 0)
        self.assertEqual(on["analyzed"], engine.FEED_AUTO_ANALYZE_CAP)
        self.assertEqual(self.store.count_analyses(), engine.FEED_AUTO_ANALYZE_CAP)
        self.assertEqual(len(self.store.feed_items(status="new")), 8 - engine.FEED_AUTO_ANALYZE_CAP)


if __name__ == "__main__":
    unittest.main()
