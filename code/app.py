"""One Health Intelligence Platform: web app (Streamlit entry point).

    streamlit run code/app.py

Multi-page site with top navigation (no sidebar). Each page has its own
colour theme and animated header (see ui_theme.py); all analysis logic
lives in engine.py.
"""
import html
import re
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from engine import (
    __version__, anthropic_client, init_session_state, run_startup_tests, service_status,
    run_amr_pipeline, run_amr_batch, log_amr_result, run_epidemic_pipeline, parse_fasta_files,
    EpidemicDetectionEngine, ArticleCredibilityChecker, AMRFinderClient, ENGINE_MODES,
    LOCATION_CACHE, LOCATION_CACHE_LOCK, _save_location_cache, maybe_geocode_in_background,
    get_store, start_feed_worker, sync_feed, analyze_feed_item, FEED_INTERVAL_S, FEED_AUTO_ANALYZE_CAP,
)
import ui_theme as ui

PAGES = {}  # name -> st.Page, filled in main(); used by the home page buttons

_PLOT_LAYOUT = dict(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
                    plot_bgcolor="rgba(0,0,0,0)", font=dict(family="Inter, sans-serif"))


# =====================================================================
# HOME
# =====================================================================
def page_home():
    ui.home_hero()

    st.write("")
    c1, c2, c3, _ = st.columns([1.15, 1.15, 1, 1.2])
    if c1.button("🧬  Screen a sequence", type="primary", width="stretch"):
        st.switch_page(PAGES["genomic"])
    if c2.button("🦠  Analyze a health report", width="stretch"):
        st.switch_page(PAGES["epidemic"])
    if c3.button("🌍  Open dashboard", width="stretch"):
        st.switch_page(PAGES["dashboard"])

    ui.home_stats()

    ui.section_heading("home", "What it does", "Three tools, <span class='grad'>one health picture</span>",
                       "Human, animal and environmental signals in one place: from a bacterial gene to a city-wide outbreak.")
    ui.home_features()

    ui.section_heading("home", "How it works", "From raw input to an <span class='grad'>actionable report</span>",
                       "Every analysis follows the same five steps and leaves an auditable trail.")
    ui.home_pipeline()

    ui.home_disclaimer()


# =====================================================================
# GENOMIC SURVEILLANCE
# =====================================================================
def page_genomic():
    ui.page_hero("genomic", "Genomic AMR Surveillance",
                 "Antimicrobial <span class='grad'>resistance screening</span>",
                 "Paste a DNA or protein sequence (or upload FASTA files). NCBI AMRFinderPlus finds the "
                 "resistance genes it carries; we score the evidence and produce an HL7 lab report.")

    amr_info = AMRFinderClient.info()
    if amr_info.get("ok"):
        st.caption(f"🧬 Detection engine: **AMRFinderPlus {amr_info['software']}** · "
                   f"reference database **{amr_info['database']}**")
    else:
        st.warning(f"⚠️ AMRFinderPlus is unavailable ({amr_info.get('error', 'unknown')}). "
                   "Using web BLAST with keyword classification instead, which is slower and less reliable.")

    col1, col2 = st.columns([1, 1.6], gap="large")

    with col1:
        with st.container(border=True):
            st.markdown("#### Input")
            sample_id = st.text_input("Sample ID", value=st.session_state['default_sample_id'])
            sequence = st.text_area("DNA / Protein Sequence", height=200,
                                    placeholder=">header\nATGCGAC...\nor raw sequence")
            cr, cc = st.columns(2)
            trigger = cr.button("▶️ Analyze", type="primary", width="stretch")
            if cc.button("🔄 Reset", width="stretch"):
                st.session_state.pop('last_amr_result', None)
                st.session_state.pop('last_amr_batch', None)
                st.rerun()

        with st.expander("📁 Batch: upload FASTA file(s)"):
            uploaded_fastas = st.file_uploader(
                "FASTA files", type=["fa", "fasta", "fna", "txt"],
                accept_multiple_files=True, key="fasta_uploader")
            batch_trigger = st.button("▶️ Analyze Batch", width="stretch", disabled=not uploaded_fastas)

        with st.expander("⚙️ Search settings"):
            st.selectbox("Detection engine", list(ENGINE_MODES.keys()),
                         format_func=ENGINE_MODES.get, key="engine_mode")
            st.slider("BLASTN identity threshold (%)", 80, 100, 90, key="identity_threshold",
                      help="Only used by web BLAST. AMRFinderPlus applies its own curated cutoff per gene.")
            st.markdown("""
**Confidence, AMRFinderPlus results**

| Component | Weight |
|---|---|
| Identity to the reference gene | 45% |
| Coverage of the reference gene | 30% |
| Detection method | 25% |

Method weights: exact / allele match 1.00 · point mutation 0.90 · BLAST above
curated cutoff 0.85 · partial 0.60–0.65 · HMM only 0.50 · internal stop 0.35

**Confidence, web BLAST results**

| Component | Weight |
|---|---|
| BLAST identity | 40% |
| E-value (log) | 30% |
| Query coverage | 15% |
| Classification source reliability | 15% |
""")

    with col2:
        with st.container(border=True):
            st.markdown("#### Results")
            _genomic_run(trigger, batch_trigger, sample_id, sequence, uploaded_fastas,
                         st.session_state.get("identity_threshold", 90),
                         st.session_state.get("engine_mode", "amrfinder"))
            _genomic_render()


def _genomic_run(trigger, batch_trigger, sample_id, sequence, uploaded_fastas, identity_threshold, mode):
    # Results are stored in session_state so they survive unrelated reruns,
    # and successful ones are saved to the database.
    if trigger:
        if not sequence.strip():
            st.session_state['last_amr_result'] = {"ok": False, "error": "❌ Enter a sequence", "sample_id": sample_id}
        else:
            try:
                result = run_amr_pipeline(sample_id, sequence, identity_threshold, mode)
                log_amr_result(result)
                st.session_state['last_amr_result'] = result
            except Exception as exc:
                st.session_state['last_amr_result'] = {"ok": False, "error": f"❌ Error: {exc}", "sample_id": sample_id}

    if batch_trigger and uploaded_fastas:
        batch_sequences = parse_fasta_files(uploaded_fastas)
        if not batch_sequences:
            st.session_state['last_amr_batch'] = []
            st.warning("⚠️ No usable sequences found in the uploaded file(s)")
        else:
            progress = st.progress(0.0, text=f"Analyzing {len(batch_sequences)} sequence(s)…")
            batch_results = run_amr_batch(
                batch_sequences, identity_threshold, mode,
                progress=lambda i, total: progress.progress(i / total, text=f"Analyzed {i}/{total}"))
            for res in batch_results:
                log_amr_result(res)
            progress.empty()
            st.session_state['last_amr_batch'] = batch_results


def _genomic_render():
    result = st.session_state.get('last_amr_result')
    if result:
        if result.get("negative"):
            st.success(f"✅ {result['error']}")
            st.caption("A negative result means no match in the reference database. It does not prove "
                       "the organism is susceptible: point mutations and novel genes are not covered here.")
        elif not result.get("ok"):
            st.error(result.get("error", "Unknown error"))
        else:
            if result["is_demo"]:
                st.markdown('<div class="demo-alert">⚠️ DEMO MODE — Results randomly generated</div>',
                            unsafe_allow_html=True)
            best, gene_info = result["best"], result["gene_info"]
            via_amrfinder = result.get("engine") == "AMRFINDERPLUS"
            if via_amrfinder:
                st.caption(f"Called by AMRFinderPlus {result.get('tool_version', '?')} · database "
                           f"{result.get('db_version', '?')} · {gene_info.get('element_name', '')}")

            if not via_amrfinder and re.search(
                    r'chromosome\s+\d|complete\s+(genome|sequence)|whole.genome|plasmid\s+\w+,\s+complete',
                    best["title"], re.IGNORECASE):
                st.warning("⚠️ Top BLAST hit is a whole-genome assembly, not a named resistance gene entry. "
                           "The sequence likely integrated into a host chromosome. Classification will use "
                           "the sample ID and sequence context only.")

            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Gene", result["display_gene"])
            m2.metric("Confidence", f"{result['confidence']*100:.1f}%")
            m3.metric("Identity", f"{best['identity_percentage']:.1f}%")
            m4.metric("Risk", f"{result['emoji']} {result['risk']}")

            st.dataframe(pd.DataFrame([{
                "Sample":     result["sample_id"],
                "Gene Class": gene_info.get("class", "Unknown"),
                "Resistance": gene_info.get("resistance", "Unknown"),
                "Concern":    gene_info.get("clinical_concern", "Unknown"),
                "Identity %": f"{best['identity_percentage']:.2f}",
                ("Method" if via_amrfinder else "E-value"):
                    gene_info.get("method", "") if via_amrfinder else f"{best['e_value']:.2e}",
                "Coverage %": f"{best.get('query_coverage', 0):.1f}",
                "Source":     gene_info.get("source", "Unknown"),
            }]), width="stretch", hide_index=True)

            if via_amrfinder:
                findings = result.get("findings", [])
                if len(findings) > 1:
                    st.markdown(f"##### 🧬 All {len(findings)} resistance genes found in this sequence")
                with st.expander("🔍 AMRFinderPlus findings", expanded=len(findings) > 1):
                    st.dataframe(pd.DataFrame([{
                        "Gene": f["symbol"], "Name": f["name"], "Class": f["class"], "Subclass": f["subclass"],
                        "Identity %": f["identity"], "Coverage %": f["coverage"], "Method": f["method"],
                        "Reference": f["accession"],
                    } for f in findings]), width="stretch", hide_index=True)
            else:
                with st.expander("🔍 Raw BLAST hit details"):
                    st.write(f"**Raw title (hit #1):** `{best['title']}`")
                    st.write(f"**Accession:** `{best['accession']}`")
                    st.write(f"**Classification source:** `{gene_info.get('source', '?')}`")
                    st.write(f"**Total hits scanned:** `{len(result['hits'])}`")
                    if len(result['hits']) > 1:
                        st.write("**Top 5 hit titles:**")
                        for i, h in enumerate(result['hits'][:5], 1):
                            st.write(f"  {i}. `{h['title'][:100]}`")

            st.markdown("##### 📋 HL7 ORU^R01 (v2.5)")
            st.code(result["hl7"], language="text")
            st.download_button("📥 Download HL7", data=result["hl7"],
                               file_name=f"{result['sample_id']}_AMR.hl7", mime="text/plain")

    batch = st.session_state.get('last_amr_batch')
    if batch:
        st.markdown("##### 📦 Batch Results")
        rows = []
        for r in batch:
            if r.get("ok"):
                genes = ", ".join(f["symbol"] for f in r.get("findings", [])) or r["display_gene"]
                rows.append({"Sample": r["sample_id"], "Gene(s)": genes,
                             "Resistance": r["gene_info"].get("resistance", ""),
                             "Confidence": f"{r['confidence']*100:.1f}%",
                             "Identity %": f"{r['best']['identity_percentage']:.2f}",
                             "Risk": f"{r['emoji']} {r['risk']}",
                             "Source": r["gene_info"].get("source", "Unknown")})
            else:
                rows.append({"Sample": r.get("sample_id", "?"),
                             "Gene(s)": "none detected" if r.get("negative") else "—",
                             "Resistance": "—", "Confidence": "—", "Identity %": "—",
                             "Risk": "—" if r.get("negative") else r.get("error", "No match"),
                             "Source": "AMRFINDERPLUS" if r.get("negative") else "—"})
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        all_hl7 = "\n\n".join(r["hl7"] for r in batch if r.get("ok"))
        if all_hl7:
            st.download_button("📥 Download Combined HL7 (batch)", data=all_hl7,
                               file_name="batch_AMR.hl7", mime="text/plain")

    if not result and not batch:
        st.info("👈 Enter a sequence (or upload FASTA files) and click **Analyze**.")


# =====================================================================
# EPIDEMIC DETECTION
# =====================================================================
DEMO_ARTICLES = {
    "COVID-19 Outbreak (Real)": "covid_2020",
    "Mpox Cluster (Real)":      "mpox_2022",
    "Routine Report (Real)":    "normal",
    "5G Conspiracy (Fake)":     "fake_1",
    "Miracle Cure Hoax (Fake)": "fake_2",
}


def page_epidemic():
    ui.page_hero("epidemic", "Epidemic Early Warning",
                 "Outbreak <span class='grad'>early-warning radar</span>",
                 "Paste a news story, WHO alert or hospital report. We screen it for misinformation first, "
                 "then extract the disease, case count, trend and location, and rate the outbreak risk.")

    if not anthropic_client:
        st.info("ℹ️ Claude API is not connected, so analysis uses the keyword fallback. "
                "See the **Status** page for details.")

    col1, col2 = st.columns([1, 1.6], gap="large")

    with col1:
        with st.container(border=True):
            st.markdown("#### Input")
            report_id = st.text_input("Report ID", value=st.session_state['default_report_id'])
            input_mode = st.radio("Input type", ["Paste Article", "Use Demo"], horizontal=True)
            if input_mode == "Paste Article":
                article = st.text_area("Health article / alert text", height=220,
                                       placeholder="Paste a real health article, WHO alert, or hospital report...")
            else:
                demo_key = st.selectbox("Demo article", list(DEMO_ARTICLES.keys()))
                article = EpidemicDetectionEngine.DEMO_ARTICLES[DEMO_ARTICLES[demo_key]]
                st.text_area("Preview", value=article, height=140, disabled=True)
            cr, cc = st.columns(2)
            trigger = cr.button("▶️ Analyze", type="primary", width="stretch")
            if cc.button("🔄 Reset", width="stretch"):
                st.session_state.pop('last_epi_result', None)
                st.rerun()

    with col2:
        with st.container(border=True):
            st.markdown("#### Results")
            if trigger:
                _epidemic_run(report_id, article)
            _epidemic_render()


def _epidemic_run(report_id, article):
    if not article.strip():
        st.session_state['last_epi_result'] = {"ok": False, "error": "❌ Enter article text"}
        return
    try:
        with st.spinner("Screening and analyzing…"):
            st.session_state['last_epi_result'] = run_epidemic_pipeline(report_id, article)
    except Exception as exc:
        st.session_state['last_epi_result'] = {"ok": False, "error": f"❌ Error: {exc}"}


def _epidemic_render():
    epi = st.session_state.get('last_epi_result')
    if not epi:
        st.info("👈 Paste an article (or pick a demo) and click **Analyze**.")
        return
    if not epi.get("ok"):
        st.error(epi.get("error", "Unknown error"))
        return

    st.markdown("##### 🔍 Article Credibility Check")
    ArticleCredibilityChecker.render_verdict(epi["cred_result"])
    if epi["blocked"]:
        return

    st.markdown("##### 📊 Epidemic Risk Analysis")
    if epi["cred_result"]["verdict"] == "UNVERIFIABLE":
        st.warning("⚠️ This article could not be fully verified. Treat the risk score below with caution "
                   "and cross-reference with official sources before acting.")
    analysis = epi["analysis"]
    if epi.get("dampened_note"):
        st.info(epi["dampened_note"])
    if epi["is_demo"]:
        reason = ("Claude API not connected" if not anthropic_client
                  else "Claude call failed, keyword-based fallback "
                       f"({html.escape(str(analysis.get('error', 'unknown error')))})")
        st.markdown(f'<div class="demo-alert">⚠️ {analysis["source"]} — {reason}</div>',
                    unsafe_allow_html=True)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Disease", analysis["disease"])
    m2.metric("Risk Score", f"{epi['risk_score']*100:.1f}%")
    m3.metric("Trend", analysis["trend"].upper())
    m4.metric("Severity", f"{epi['emoji']} {epi['risk']}")

    summary = analysis["summary"]
    st.dataframe(pd.DataFrame([{
        "Report":      epi["report_id"],
        "Location":    analysis["location"],
        "Cases":       analysis["cases"],
        "Baseline":    analysis["baseline"],
        "Severity":    analysis["severity"],
        "Credibility": f"{epi['cred_result']['credibility_score']*100:.0f}%",
        "Summary":     summary if len(summary) <= 80 else summary[:80] + "…",
    }]), width="stretch", hide_index=True)

    st.markdown("##### 📋 Recommended Actions")
    rs = epi["risk_score"]
    if rs >= 0.85:
        st.error("🔴 CRITICAL → Activate emergency response, notify public health authorities immediately")
    elif rs >= 0.70:
        st.warning("🟠 HIGH → Increase surveillance, prepare response team, issue health alert")
    elif rs >= 0.50:
        st.info("🟡 MODERATE → Monitor closely, prepare communication plan")
    else:
        st.success("🟢 LOW → Routine monitoring, standard protocols")

    st.markdown("##### 📋 HL7 Epidemic Alert (v2.5)")
    st.code(epi["hl7"], language="text")
    st.download_button("📥 Download HL7 Alert", data=epi["hl7"],
                       file_name=f"{epi['report_id']}_epidemic_alert.hl7", mime="text/plain")


# =====================================================================
# OUTBREAK FEED
# =====================================================================
def page_feed():
    ui.page_hero("feed", "Outbreak Feed",
                 "WHO outbreak news, <span class='grad'>watched for you</span>",
                 "New WHO Disease Outbreak News reports are pulled in automatically, screened, and scored "
                 "with the same pipeline as the Epidemic page. Results land on the dashboard map and audit trail.")

    store = get_store()
    c1, c2 = st.columns([1, 2.6], vertical_alignment="center")
    if c1.button("↻ Fetch latest now", type="primary", width="stretch"):
        with st.spinner("Fetching WHO Disease Outbreak News…"):
            summary = sync_feed(store)
        if summary["error"]:
            st.error(f"Could not fetch the feed: {summary['error']}")
        else:
            st.toast(f"{summary['new']} new report(s), {summary['analyzed']} analyzed")

    saved_auto = store.get_setting("feed_auto_analyze", "1") == "1"
    auto = c2.toggle(f"Analyze new reports automatically with Claude (up to {FEED_AUTO_ANALYZE_CAP} per check)",
                     value=saved_auto, help="Each analysis makes two Claude API calls and is billed to your API key.")
    if auto != saved_auto:
        store.set_setting("feed_auto_analyze", "1" if auto else "0")

    last = store.get_setting("feed_last_fetch", "")
    note = (f"Last checked {last}" if last else "Not checked yet; the first check runs a few seconds after the app starts")
    st.caption(f"{note} · checks every {FEED_INTERVAL_S // 3600} hours while the app is running · "
               "source: who.int Disease Outbreak News")
    if store.get_setting("feed_last_error", ""):
        st.warning(f"The last automatic check failed: {store.get_setting('feed_last_error')}")
    if not anthropic_client:
        st.info("ℹ️ Claude API is not connected, so reports are collected but not analyzed. "
                "Set ANTHROPIC_API_KEY and restart to turn analysis on.")

    items = store.feed_items(limit=20)
    if not items:
        st.info("No reports collected yet. Click **Fetch latest now**.")
        return

    m1, m2, m3 = st.columns(3)
    m1.metric("Reports shown", len(items))
    m2.metric("Analyzed", sum(1 for i in items if i["status"] == "analyzed"))
    m3.metric("Waiting", sum(1 for i in items if i["status"] == "new"))

    for item in items:
        with st.container(border=True):
            rec = store.get_analysis(item["analysis_id"]) if item.get("analysis_id") else None
            status = item["status"] if (rec or item["status"] != "analyzed") else "new"
            ui.feed_item_header(item["title"], item["url"], item["published_at"][:10], item["source"], status)

            if rec:
                a = rec["details"].get("analysis", {})
                k1, k2, k3, k4 = st.columns(4)
                k1.metric("Disease", a.get("disease", "?"))
                k2.metric("Risk", f"{rec['confidence'] * 100:.0f}% · {rec['details'].get('risk', '')}")
                k3.metric("Trend", str(a.get("trend", "?")).upper())
                k4.metric("Location", a.get("location", "?"))
                with st.expander("Summary, credibility and HL7 alert"):
                    st.write(a.get("summary", ""))
                    cred = rec["details"].get("credibility", {})
                    st.caption(f"Cases reported: {a.get('cases', '?')} · Credibility: {cred.get('verdict', '?')} "
                               f"({float(cred.get('credibility_score') or 0) * 100:.0f}%) · Source: {rec['source']}")
                    st.code(rec["hl7"], language="text")
            else:
                if item.get("note"):
                    st.caption(item["note"])
                label = "Analyze" if status == "new" else "Retry analysis"
                if st.button(label, key=f"feed_analyze_{item['id']}", disabled=not anthropic_client,
                             help=None if anthropic_client else "Needs the Claude API"):
                    with st.spinner("Analyzing…"):
                        analyze_feed_item(item, store)
                    st.rerun()


# =====================================================================
# CORRELATION DASHBOARD
# =====================================================================
_KNOWN_LOCATIONS = {
    'new york': {'lon': -74.0060, 'lat': 40.7128}, 'ny': {'lon': -74.0060, 'lat': 40.7128},
    'lagos': {'lon': 3.3792, 'lat': 6.5244}, 'nairobi': {'lon': 36.8219, 'lat': -1.2921},
    'tokyo': {'lon': 139.6917, 'lat': 35.6895}, 'johannesburg': {'lon': 28.0473, 'lat': -26.2041},
    'london': {'lon': -0.1278, 'lat': 51.5074}, 'beijing': {'lon': 116.4074, 'lat': 39.9042},
    'nigeria': {'lon': 8.6753, 'lat': 9.0820},
}


def page_dashboard():
    ui.page_hero("dashboard", "One Health Correlation",
                 "Signals, <span class='grad'>side by side</span>",
                 "Every saved genomic finding and epidemic signal, plotted over time and on a world map.")

    audit_df = get_store().get_logs()
    if audit_df.empty:
        st.info("No analyses yet. Run something on the **Genomic** or **Epidemic** pages and it will appear here.")
    else:
        genomic_df = audit_df[audit_df["Type"] == "GENOMIC"]
        epidemic_df = audit_df[audit_df["Type"] == "EPIDEMIC"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Analyses", len(audit_df))
        c2.metric("Genomic Findings", len(genomic_df))
        c3.metric("Epidemic Signals", len(epidemic_df))
        c4.metric("Real Results", len(audit_df[audit_df["Status"] == "✅ Real"]))

        cols = ["Timestamp", "Sample/Report", "Finding", "Confidence", "Source", "Status"]
        t1, t2 = st.columns(2, gap="large")
        with t1:
            st.markdown("##### 🧬 Genomic Findings")
            if genomic_df.empty:
                st.caption("None yet.")
            else:
                st.dataframe(genomic_df[cols], width="stretch", hide_index=True)
        with t2:
            st.markdown("##### 🦠 Epidemic Signals")
            if epidemic_df.empty:
                st.caption("None yet.")
            else:
                st.dataframe(epidemic_df[cols], width="stretch", hide_index=True)

        if len(audit_df) >= 2:
            st.markdown("##### 📊 Confidence Timeline")
            conf = audit_df["Confidence"].str.replace("%", "").astype(float)
            colors = ["#F87171" if c >= 85 else "#FB923C" if c >= 70 else "#FBBF24" if c >= 50 else "#4ADE80"
                      for c in conf]
            fig = go.Figure(go.Bar(x=audit_df["Timestamp"], y=conf, marker_color=colors,
                                   text=audit_df["Finding"], textposition="outside"))
            fig.update_layout(yaxis_title="Confidence %", height=350, margin=dict(t=20, b=40), **_PLOT_LAYOUT)
            st.plotly_chart(fig, width="stretch")

    st.markdown("##### 🌍 Global Map")
    st.checkbox("Look up unknown locations with Nominatim (background, cached)", value=True,
                key="enable_geocoding")
    _render_map(audit_df)

    with st.expander("Location cache"):
        st.caption("Shared by every user of this app instance.")
        st.json(LOCATION_CACHE, expanded=False)
        cc1, cc2 = st.columns(2)
        if cc1.button("Clear entire location cache"):
            with LOCATION_CACHE_LOCK:
                LOCATION_CACHE.clear()
                _save_location_cache()
            st.success("Location cache cleared (disk cache updated).")
            st.rerun()
        if cc2.button("Retry failed lookups now"):
            with LOCATION_CACHE_LOCK:
                failed = [loc for loc, e in LOCATION_CACHE.items()
                          if isinstance(e, dict) and e.get('value') is None]
                for loc in failed:
                    del LOCATION_CACHE[loc]
                _save_location_cache()
            for loc in failed:
                maybe_geocode_in_background(loc)
            st.success(f"Re-queued {len(failed)} failed lookup(s) for retry.")
            st.rerun()


def _render_map(audit_df):
    points_by_type = {}
    for _, row in (audit_df.iterrows() if not audit_df.empty else []):
        loc = row.get('Location', '') if isinstance(row.get('Location', ''), str) else ''
        if not loc:
            continue
        with LOCATION_CACHE_LOCK:
            entry = LOCATION_CACHE.get(loc)
        latlon = entry.get('value') if entry else None
        if latlon is None:
            lk = loc.lower()
            for k, v in _KNOWN_LOCATIONS.items():
                if k in lk:
                    latlon = v
                    with LOCATION_CACHE_LOCK:
                        LOCATION_CACHE[loc] = {'value': v, 'attempts': 0,
                                               'last_try': datetime.now().timestamp(), 'status': 'ok'}
                        _save_location_cache()
                    break
            if latlon is None:
                try:
                    maybe_geocode_in_background(loc)
                except Exception:
                    pass
        if not latlon:
            continue
        try:
            confv = float(str(row.get('Confidence', '0%')).replace('%', ''))
        except Exception:
            confv = 50.0
        hover = (f"{row.get('Finding', '')} · {row.get('Sample/Report', '')}<br>"
                 f"Conf: {confv:.1f}% · Source: {row.get('Source', '')}")
        points_by_type.setdefault(row.get('Type', 'UNKNOWN'), []).append(
            (latlon['lon'], latlon['lat'], hover, 8 + (confv / 100.0) * 24))

    if not points_by_type:
        st.info("No mapped locations yet. Epidemic reports that name a place will appear here "
                "(new places are geocoded in the background; refresh in a few seconds).")
        return

    color_map = {'GENOMIC': '#2DD4BF', 'EPIDEMIC': '#F87171', 'UNKNOWN': '#94A3B8'}
    fig = go.Figure()
    for t, pts in points_by_type.items():
        fig.add_trace(go.Scattergeo(
            lon=[p[0] for p in pts], lat=[p[1] for p in pts], mode='markers', name=f"{t} ({len(pts)})",
            marker=dict(size=[p[3] for p in pts], color=color_map.get(t, '#94A3B8'),
                        line=dict(width=1, color='rgba(255,255,255,.6)'), opacity=.85),
            hoverinfo='text', hovertext=[p[2] for p in pts]))
    fig.update_layout(
        geo=dict(projection_type='natural earth', showland=True, landcolor='#1E293B', showocean=True,
                 oceancolor='#0B1222', showcountries=True, countrycolor='#334155',
                 coastlinecolor='#475569', bgcolor='rgba(0,0,0,0)'),
        height=460, margin=dict(t=10, b=10, l=0, r=0), legend=dict(title='Type (count)'), **_PLOT_LAYOUT)
    st.plotly_chart(fig, width="stretch")


# =====================================================================
# AUDIT TRAIL
# =====================================================================
def page_audit():
    ui.page_hero("audit", "Audit Trail",
                 "Every analysis, <span class='grad'>on the record</span>",
                 "A permanent, timestamped log of every finding, labelled real or demo. "
                 "Saved to the database, so it survives refreshes and restarts.")

    store = get_store()
    audit_df = store.get_logs()
    if audit_df.empty:
        st.info("No analyses yet.")
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Analyses", len(audit_df))
    c2.metric("Real Results", len(audit_df[audit_df["Status"] == "✅ Real"]))
    c3.metric("Demo Results", len(audit_df[audit_df["Status"] == "⚠️ DEMO"]))

    real_only = st.toggle("Show real results only", value=False)
    display_df = audit_df[audit_df["Status"] == "✅ Real"] if real_only else audit_df
    st.dataframe(display_df.iloc[::-1], width="stretch", hide_index=True)
    st.download_button("📥 Export Audit CSV", data=display_df.to_csv(index=False).encode(),
                       file_name=f"audit_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv")

    st.markdown("##### 📂 Reopen a saved report")
    id_by_label = {f"#{r['ID']} · {r['Timestamp']} · {r['Type']} · {r['Sample/Report']} · {r['Finding']}": int(r["ID"])
                   for _, r in display_df.iloc[::-1].iterrows()}
    chosen = st.selectbox("Saved analysis", list(id_by_label), index=None,
                          placeholder="Choose an analysis to view its details and HL7 message")
    if chosen is not None:
        rec = store.get_analysis(id_by_label[chosen])
        if rec:
            if rec.get("hl7"):
                st.code(rec["hl7"], language="text")
                st.download_button("📥 Download HL7", data=rec["hl7"], mime="text/plain",
                                   file_name=f"{rec['sample_id']}_{rec['type'].lower()}.hl7")
            with st.expander("Stored details"):
                st.json(rec["details"])

    with st.expander("🗑️ Delete all saved analyses"):
        st.caption("Removes every saved analysis for all users of this app. This cannot be undone. "
                   "Export the CSV first if you need a copy.")
        confirm = st.checkbox("I understand, delete everything")
        if st.button("Delete all saved analyses", disabled=not confirm):
            n = store.clear_analyses()
            st.success(f"Deleted {n} saved analyses.")
            st.rerun()


# =====================================================================
# SYSTEM STATUS
# =====================================================================
def page_status():
    ui.page_hero("status", "System Status",
                 "Service <span class='grad'>health</span>",
                 "Live checks of every external service the platform depends on. "
                 "Results are cached for 5 minutes.")

    if st.button("↻ Recheck services"):
        service_status.clear()
    with st.spinner("Checking services…"):
        status = service_status()
    ui.status_cards(status["results"])
    st.caption(f"Last checked {status['checked_at']}")

    st.markdown("##### 🧪 Unit tests")
    tests_ok = st.session_state.get('startup_tests_ok')
    if tests_ok is None:
        st.info("Not configured: no `tests/` folder next to app.py.")
    elif tests_ok:
        st.success("All startup unit tests passed.")
    else:
        st.error("Startup unit tests failed.")
        with st.expander("View test output"):
            st.code(st.session_state.get('startup_test_output', ''), language='text')


# =====================================================================
# APP SHELL
# =====================================================================
def setup():
    st.set_page_config(page_title="One Health Intelligence Platform", page_icon="🧬",
                       layout="wide", initial_sidebar_state="collapsed")
    ui.inject_global_css()
    init_session_state()
    st.session_state.setdefault('default_sample_id', f"ISO_{datetime.now().strftime('%H%M%S')}")
    st.session_state.setdefault('default_report_id', f"RPT_{datetime.now().strftime('%H%M%S')}")
    if 'startup_tests_ok' not in st.session_state:
        run_startup_tests()
    start_feed_worker()


def main():
    setup()
    PAGES.update({
        "home":      st.Page(page_home, title="Home", icon=":material/home:", default=True),
        "genomic":   st.Page(page_genomic, title="Genomic", icon=":material/genetics:", url_path="genomic"),
        "epidemic":  st.Page(page_epidemic, title="Epidemic", icon=":material/coronavirus:", url_path="epidemic"),
        "feed":      st.Page(page_feed, title="Outbreak Feed", icon=":material/rss_feed:", url_path="feed"),
        "dashboard": st.Page(page_dashboard, title="Dashboard", icon=":material/public:", url_path="dashboard"),
        "audit":     st.Page(page_audit, title="Audit Trail", icon=":material/receipt_long:", url_path="audit"),
        "status":    st.Page(page_status, title="Status", icon=":material/monitor_heart:", url_path="status"),
    })
    st.navigation(list(PAGES.values()), position="top").run()
    ui.footer(__version__)


if __name__ == "__main__":
    main()
