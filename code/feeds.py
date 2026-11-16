"""One Health Intelligence Platform: outbreak news feeds.

Currently one source: WHO Disease Outbreak News (DON), via WHO's public
JSON API. Parsing is separated from fetching so it can be unit-tested
offline. Streamlit-free.
"""
import html
import re
from typing import Any, Dict, List

import requests

WHO_DON_API = "https://www.who.int/api/news/diseaseoutbreaknews"
WHO_DON_ITEM_URL = "https://www.who.int/emergencies/disease-outbreak-news/item/"
MAX_TEXT_CHARS = 12_000   # per item sent for analysis; longer reports are cut and flagged
_UA = {"User-Agent": "OneHealthPlatform/1.0 (outbreak feed reader)"}


def html_to_text(fragment: str) -> str:
    """Strip tags from a WHO HTML fragment, keeping paragraph breaks."""
    if not fragment:
        return ""
    text = re.sub(r"(?is)<(script|style).*?>.*?</\1>", " ", fragment)
    text = re.sub(r"(?i)<\s*(br|/p|/div|/li|/h[1-6])\s*/?>", "\n", text)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text).replace("\xa0", " ")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r" *\n *", "\n", text)
    return re.sub(r"\n{3,}", "\n\n", text).strip()


def parse_who_don(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Turn the WHO DON API response into feed items."""
    items = []
    for it in (payload or {}).get("value", []) or []:
        item_id = it.get("DonId") or it.get("UrlName") or it.get("Id")
        if not item_id:
            continue
        title = (it.get("OverrideTitle") if it.get("UseOverrideTitle") and it.get("OverrideTitle")
                 else it.get("Title")) or "Untitled"
        published = (it.get("PublicationDate") or "")[:10]
        parts = [
            f"WHO Disease Outbreak News. Published {published}." if published else "WHO Disease Outbreak News.",
            title.strip(),
            html_to_text(it.get("Summary", "")),
            html_to_text(it.get("Overview", "")),
            html_to_text(it.get("Epidemiology", "")),
        ]
        # Summary and Overview often repeat each other; drop exact duplicates.
        seen, body = set(), []
        for p in parts:
            if p and p not in seen:
                seen.add(p)
                body.append(p)
        text = "\n\n".join(body)
        truncated = len(text) > MAX_TEXT_CHARS
        items.append({
            "id": str(item_id),
            "source": "WHO DON",
            "title": title.strip(),
            "url": WHO_DON_ITEM_URL + str(it.get("UrlName") or item_id),
            "published_at": it.get("PublicationDate") or "",
            "text": text[:MAX_TEXT_CHARS],
            "truncated": truncated,
        })
    return items


def fetch_who_don(limit: int = 10, timeout: int = 30) -> List[Dict[str, Any]]:
    """Newest WHO Disease Outbreak News items. Raises on network/HTTP errors."""
    r = requests.get(WHO_DON_API,
                     params={"$orderby": "PublicationDate desc", "$top": str(int(limit))},
                     headers=_UA, timeout=timeout)
    r.raise_for_status()
    return parse_who_don(r.json())
