"""One Health Intelligence Platform: persistence (SQLite).

Keeps analyses, the audit trail and outbreak-feed items across restarts.
Deliberately free of Streamlit so it can be unit-tested and used from
background threads. One short-lived connection per operation; WAL mode so
the feed worker and the web app can write concurrently.

What is stored: the finding, scores, the generated HL7 message and a JSON
"details" blob. Raw sequences and pasted article text are NOT stored.
"""
import json
import os
import sqlite3
from contextlib import contextmanager
from datetime import datetime
from typing import Any, Dict, List, Optional

import pandas as pd

DATA_DIR = os.environ.get("OHP_DATA_DIR") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
DB_PATH = os.path.join(DATA_DIR, "onehealth.db")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS analyses (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at  TEXT    NOT NULL,
    type        TEXT    NOT NULL,          -- GENOMIC | EPIDEMIC
    sample_id   TEXT    NOT NULL,
    finding     TEXT    NOT NULL,
    confidence  REAL    NOT NULL,          -- 0..1
    source      TEXT    NOT NULL,
    location    TEXT    NOT NULL DEFAULT '',
    is_demo     INTEGER NOT NULL DEFAULT 0,
    details     TEXT    NOT NULL DEFAULT '{}',
    hl7         TEXT    NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_analyses_created ON analyses(created_at);

CREATE TABLE IF NOT EXISTS feed_items (
    id           TEXT PRIMARY KEY,         -- e.g. WHO DonId
    source       TEXT NOT NULL,
    title        TEXT NOT NULL,
    url          TEXT NOT NULL DEFAULT '',
    published_at TEXT NOT NULL DEFAULT '',
    fetched_at   TEXT NOT NULL,
    text         TEXT NOT NULL DEFAULT '',
    status       TEXT NOT NULL DEFAULT 'new',   -- new | analyzed | blocked | error
    note         TEXT NOT NULL DEFAULT '',
    analysis_id  INTEGER
);

CREATE TABLE IF NOT EXISTS settings (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class Store:
    def __init__(self, path: str = DB_PATH):
        self.path = path
        folder = os.path.dirname(os.path.abspath(path))
        os.makedirs(folder, exist_ok=True)
        with self._conn() as con:
            con.execute("PRAGMA journal_mode=WAL")
            con.executescript(_SCHEMA)

    @contextmanager
    def _conn(self):
        con = sqlite3.connect(self.path, timeout=15)
        con.row_factory = sqlite3.Row
        try:
            yield con
            con.commit()
        finally:
            con.close()

    # ------------------------------------------------------------ analyses
    def log_analysis(self, sample_id: str, analysis_type: str, finding: str, confidence: float,
                     source: str, is_demo: bool = False, location: str = "",
                     details: Optional[Dict[str, Any]] = None, hl7: str = "") -> int:
        with self._conn() as con:
            cur = con.execute(
                "INSERT INTO analyses (created_at, type, sample_id, finding, confidence, source, "
                "location, is_demo, details, hl7) VALUES (?,?,?,?,?,?,?,?,?,?)",
                (_now(), analysis_type, str(sample_id), str(finding), float(confidence), str(source),
                 location or "", 1 if is_demo else 0,
                 json.dumps(details or {}, default=str, ensure_ascii=False), hl7 or ""))
            return int(cur.lastrowid)

    def get_logs(self) -> pd.DataFrame:
        """Audit table, oldest first, in the column layout the UI expects."""
        with self._conn() as con:
            rows = con.execute(
                "SELECT id, created_at, type, sample_id, finding, confidence, source, location, is_demo "
                "FROM analyses ORDER BY id").fetchall()
        if not rows:
            return pd.DataFrame()
        return pd.DataFrame([{
            "ID":            r["id"],
            "Timestamp":     r["created_at"],
            "Type":          r["type"],
            "Sample/Report": r["sample_id"],
            "Finding":       r["finding"],
            "Confidence":    f"{r['confidence'] * 100:.1f}%",
            "Source":        r["source"],
            "Location":      r["location"],
            "Status":        "⚠️ DEMO" if r["is_demo"] else "✅ Real",
        } for r in rows])

    def get_analysis(self, analysis_id: int) -> Optional[Dict[str, Any]]:
        with self._conn() as con:
            r = con.execute("SELECT * FROM analyses WHERE id = ?", (analysis_id,)).fetchone()
        if r is None:
            return None
        out = dict(r)
        try:
            out["details"] = json.loads(out.get("details") or "{}")
        except Exception:
            out["details"] = {}
        return out

    def count_analyses(self) -> int:
        with self._conn() as con:
            return int(con.execute("SELECT COUNT(*) FROM analyses").fetchone()[0])

    def clear_analyses(self) -> int:
        with self._conn() as con:
            n = con.execute("DELETE FROM analyses").rowcount
            con.execute("UPDATE feed_items SET status = 'new', analysis_id = NULL, note = '' "
                        "WHERE status IN ('analyzed', 'blocked', 'error')")
            return int(n)

    # ---------------------------------------------------------------- feed
    def upsert_feed_items(self, items: List[Dict[str, Any]]) -> int:
        """Insert items not seen before. Returns how many were new."""
        new = 0
        with self._conn() as con:
            for it in items:
                cur = con.execute(
                    "INSERT OR IGNORE INTO feed_items (id, source, title, url, published_at, fetched_at, text) "
                    "VALUES (?,?,?,?,?,?,?)",
                    (it["id"], it.get("source", ""), it.get("title", ""), it.get("url", ""),
                     it.get("published_at", ""), _now(), it.get("text", "")))
                new += cur.rowcount
        return new

    def feed_items(self, limit: int = 30, status: Optional[str] = None) -> List[Dict[str, Any]]:
        sql = "SELECT * FROM feed_items"
        args: list = []
        if status:
            sql += " WHERE status = ?"
            args.append(status)
        sql += " ORDER BY published_at DESC, id DESC LIMIT ?"
        args.append(limit)
        with self._conn() as con:
            return [dict(r) for r in con.execute(sql, args).fetchall()]

    def set_feed_status(self, item_id: str, status: str, analysis_id: Optional[int] = None, note: str = ""):
        with self._conn() as con:
            con.execute("UPDATE feed_items SET status = ?, analysis_id = ?, note = ? WHERE id = ?",
                        (status, analysis_id, note, item_id))

    # ------------------------------------------------------------ settings
    def get_setting(self, key: str, default: str = "") -> str:
        with self._conn() as con:
            r = con.execute("SELECT value FROM settings WHERE key = ?", (key,)).fetchone()
        return r["value"] if r else default

    def set_setting(self, key: str, value: str):
        with self._conn() as con:
            con.execute("INSERT INTO settings (key, value) VALUES (?, ?) "
                        "ON CONFLICT(key) DO UPDATE SET value = excluded.value", (key, str(value)))
