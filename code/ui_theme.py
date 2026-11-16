"""One Health Intelligence Platform: visual design.

Global stylesheet, per-page colour themes and pure-CSS/SVG animated motifs
(DNA helix, radar, orbit, scanner, heartbeat). No JavaScript: everything is
rendered with st.html, and all motion is disabled for users who prefer
reduced motion.
"""
import html as _html

import streamlit as st

# ---------------------------------------------------------------------
# Global stylesheet
# ---------------------------------------------------------------------
GLOBAL_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Space+Grotesk:wght@500;600;700&family=JetBrains+Mono:wght@400;600&display=swap');

:root{
  --bg:#060A13; --surface:rgba(15,23,42,.62); --surface-2:rgba(30,41,59,.55);
  --border:rgba(148,163,184,.16); --text:#E2E8F0; --muted:#94A3B8;
  --teal:#2DD4BF; --violet:#A78BFA; --amber:#FBBF24; --red:#F87171;
  --blue:#60A5FA; --green:#4ADE80; --slate:#94A3B8; --sky:#38BDF8;
  --display:'Space Grotesk',system-ui,sans-serif; --body:'Inter',system-ui,sans-serif;
  --mono:'JetBrains Mono',ui-monospace,monospace;
}
/* page themes: --a1 / --a2 are the two accent colours */
.t-home{--a1:var(--teal);--a2:var(--violet)}
.t-genomic{--a1:#2DD4BF;--a2:#34D399}
.t-epidemic{--a1:#FBBF24;--a2:#F87171}
.t-dashboard{--a1:#A78BFA;--a2:#60A5FA}
.t-audit{--a1:#94A3B8;--a2:#38BDF8}
.t-status{--a1:#4ADE80;--a2:#2DD4BF}
.t-feed{--a1:#FB7185;--a2:#FB923C}

.stApp{background:
   radial-gradient(1200px 600px at 85% -10%, rgba(167,139,250,.10), transparent 60%),
   radial-gradient(900px 500px at -10% 10%, rgba(45,212,191,.08), transparent 60%),
   var(--bg);}
[data-testid="stMainBlockContainer"], .block-container{max-width:1240px;padding-top:1.2rem;padding-bottom:3rem}
header[data-testid="stHeader"]{background:rgba(6,10,19,.72);backdrop-filter:blur(12px);
   border-bottom:1px solid var(--border)}
[data-testid="stSidebarCollapsedControl"],[data-testid="stSidebar"]{display:none}

/* Streamlit widgets restyled to match */
.stButton>button,.stDownloadButton>button{border-radius:999px;font-weight:600;
   border:1px solid var(--border);transition:transform .15s ease, box-shadow .15s ease}
.stButton>button:hover,.stDownloadButton>button:hover{transform:translateY(-1px);
   box-shadow:0 6px 22px rgba(45,212,191,.18)}
.stButton>button[kind="primary"],[data-testid="stBaseButton-primary"]{
   background:linear-gradient(135deg,#14B8A6,#8B5CF6);border:none;color:#fff}
[data-testid="stMetric"]{background:var(--surface);border:1px solid var(--border);
   border-radius:14px;padding:14px 16px}
[data-testid="stMetricLabel"]{color:var(--muted)}
[data-testid="stVerticalBlockBorderWrapper"]{border-radius:16px}
[data-testid="stExpander"] details{border-radius:12px;border-color:var(--border);background:var(--surface)}
textarea{font-family:var(--mono)!important;font-size:13px!important}
hr{border-color:var(--border)!important}

/* status / alert chips used by analysis results */
.demo-alert{background:#2a2a00;color:#FFD700;padding:10px 16px;border-radius:10px;
   border-left:4px solid #FFD700;margin-bottom:12px;font-weight:500}
.status-ok,.status-warn,.status-info{padding:8px 12px;border-radius:10px;font-size:13px;margin:4px 0}
.status-ok  {background:#0D2A1A;color:#4ADE80;border-left:3px solid #4ADE80}
.status-warn{background:#2A2500;color:#FBBF24;border-left:3px solid #FBBF24}
.status-info{background:#0D1A2A;color:#60A5FA;border-left:3px solid #60A5FA}

/* ---------- shared building blocks ---------- */
.ohp{font-family:var(--body);color:var(--text)}
.ohp h1,.ohp h2,.ohp h3{font-family:var(--display);letter-spacing:-.02em;margin:0}
.ohp .eyebrow{display:inline-flex;align-items:center;gap:8px;font:600 12px/1 var(--body);
   letter-spacing:.14em;text-transform:uppercase;color:var(--a1);
   background:color-mix(in srgb,var(--a1) 12%,transparent);
   border:1px solid color-mix(in srgb,var(--a1) 35%,transparent);padding:7px 12px;border-radius:999px}
.ohp .eyebrow i{width:7px;height:7px;border-radius:50%;background:var(--a1);
   box-shadow:0 0 10px var(--a1);animation:ohp-blink 1.6s ease-in-out infinite}
.ohp .grad{background:linear-gradient(100deg,var(--a1),var(--a2));-webkit-background-clip:text;
   background-clip:text;color:transparent}
.ohp .lead{color:var(--muted);font-size:17px;line-height:1.65;max-width:620px}
.ohp .section-title{font:700 30px/1.15 var(--display);margin:8px 0 6px}
.ohp .section-sub{color:var(--muted);font-size:15px;margin-bottom:22px}

/* ---------- home hero ---------- */
.ohp-hero{position:relative;overflow:hidden;border-radius:28px;border:1px solid var(--border);
   padding:56px 52px;display:grid;grid-template-columns:1.15fr .85fr;gap:24px;align-items:center;
   background:linear-gradient(160deg,rgba(15,23,42,.92),rgba(6,10,19,.92))}
.ohp-hero::before{content:"";position:absolute;inset:0;pointer-events:none;
   background-image:linear-gradient(rgba(148,163,184,.06) 1px,transparent 1px),
                    linear-gradient(90deg,rgba(148,163,184,.06) 1px,transparent 1px);
   background-size:42px 42px;mask-image:radial-gradient(ellipse at 70% 40%,#000 30%,transparent 75%)}
.ohp-blob{position:absolute;border-radius:50%;filter:blur(60px);opacity:.45;pointer-events:none;
   animation:ohp-drift 14s ease-in-out infinite alternate}
.ohp-blob.b1{width:340px;height:340px;background:var(--a1);top:-120px;right:8%}
.ohp-blob.b2{width:300px;height:300px;background:var(--a2);bottom:-140px;right:30%;animation-delay:-6s}
.ohp-hero h1{font-size:56px;line-height:1.04;margin:18px 0 18px}
.ohp-hero .chips{display:flex;flex-wrap:wrap;gap:8px;margin-top:22px}
.ohp-hero .chip{font:500 13px var(--body);color:var(--text);padding:7px 12px;border-radius:999px;
   background:var(--surface-2);border:1px solid var(--border)}
.ohp-hero .art{position:relative;display:flex;justify-content:center;min-height:360px}
.ohp-base{position:absolute;font:600 15px var(--mono);color:var(--a1);opacity:0;
   animation:ohp-rise 7s linear infinite}

/* ---------- DNA helix ---------- */
.ohp-helix{perspective:900px;display:flex;flex-direction:column;align-items:center;gap:7px;
   position:relative;z-index:1}
.ohp-rung{position:relative;width:var(--w,230px);height:10px;transform-style:preserve-3d;
   animation:ohp-spin-y 6s linear infinite}
.ohp-rung .bar{position:absolute;top:4px;left:12px;right:12px;height:2px;
   background:linear-gradient(90deg,var(--a1) 0 50%,var(--a2) 50% 100%);opacity:.5}
.ohp-rung .n{position:absolute;top:0;width:10px;height:10px;border-radius:50%}
.ohp-rung .n.l{left:0;background:var(--a1);box-shadow:0 0 14px var(--a1)}
.ohp-rung .n.r{right:0;background:var(--a2);box-shadow:0 0 14px var(--a2)}
.ohp-helix.sm{gap:4px}
.ohp-helix.sm .ohp-rung{height:8px}
.ohp-helix.sm .n{width:8px;height:8px}
.ohp-helix.sm .bar{top:3px}

/* ---------- feature cards ---------- */
.ohp-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}
.ohp-card{position:relative;background:var(--surface);border:1px solid var(--border);border-radius:20px;
   padding:26px;overflow:hidden;transition:transform .25s ease,border-color .25s ease}
.ohp-card:hover{transform:translateY(-4px);border-color:color-mix(in srgb,var(--a1) 55%,transparent)}
.ohp-card::after{content:"";position:absolute;left:0;right:0;top:0;height:3px;
   background:linear-gradient(90deg,var(--a1),var(--a2))}
.ohp-card .ic{width:48px;height:48px;border-radius:14px;display:grid;place-items:center;margin-bottom:16px;
   background:color-mix(in srgb,var(--a1) 14%,transparent);border:1px solid color-mix(in srgb,var(--a1) 35%,transparent)}
.ohp-card .ic svg{width:24px;height:24px;stroke:var(--a1);fill:none;stroke-width:1.8;stroke-linecap:round;stroke-linejoin:round}
.ohp-card h3{font-size:20px;margin-bottom:8px}
.ohp-card p{color:var(--muted);font-size:14.5px;line-height:1.6;margin:0 0 12px}
.ohp-card ul{margin:0;padding-left:18px;color:var(--text);font-size:13.5px;line-height:1.75}
.ohp-card li::marker{color:var(--a1)}

/* ---------- stats strip ---------- */
.ohp-stats{display:grid;grid-template-columns:repeat(4,1fr);gap:1px;border-radius:18px;overflow:hidden;
   background:var(--border);border:1px solid var(--border)}
.ohp-stat{background:rgba(10,15,28,.95);padding:22px 20px}
.ohp-stat b{display:block;font:700 30px var(--display);background:linear-gradient(100deg,var(--teal),var(--violet));
   -webkit-background-clip:text;background-clip:text;color:transparent}
.ohp-stat span{color:var(--muted);font-size:13.5px}

/* ---------- pipeline ---------- */
.ohp-flow{position:relative;display:grid;grid-template-columns:repeat(5,1fr);gap:14px;padding-top:8px}
.ohp-flow .track{position:absolute;left:10%;right:10%;top:36px;height:2px;
   background:linear-gradient(90deg,var(--a1),var(--a2));opacity:.35}
.ohp-flow .pulse{position:absolute;top:31px;left:10%;width:12px;height:12px;border-radius:50%;
   background:#fff;box-shadow:0 0 16px var(--a1),0 0 32px var(--a2);animation:ohp-travel 4.5s ease-in-out infinite}
.ohp-step{position:relative;text-align:center;padding:0 6px}
.ohp-step .dot{width:56px;height:56px;margin:0 auto 14px;border-radius:50%;display:grid;place-items:center;
   font:700 18px var(--display);background:#0B1222;border:2px solid var(--a1);color:var(--a1);
   box-shadow:0 0 0 6px rgba(45,212,191,.08)}
.ohp-step h4{font:600 15px var(--display);margin:0 0 6px;color:var(--text)}
.ohp-step p{font-size:13px;color:var(--muted);line-height:1.5;margin:0}

/* ---------- page hero (inner pages) ---------- */
.ohp-page-hero{position:relative;overflow:hidden;border-radius:24px;border:1px solid var(--border);
   padding:34px 40px;margin-bottom:22px;display:grid;grid-template-columns:1fr auto;gap:24px;align-items:center;
   background:linear-gradient(120deg,color-mix(in srgb,var(--a1) 16%,#0B1222) 0%,#0B1222 55%,
              color-mix(in srgb,var(--a2) 14%,#0B1222) 100%)}
.ohp-page-hero::before{content:"";position:absolute;inset:0;pointer-events:none;
   background-image:radial-gradient(color-mix(in srgb,var(--a1) 35%,transparent) 1px,transparent 1px);
   background-size:22px 22px;opacity:.25;mask-image:linear-gradient(90deg,transparent,#000 70%)}
.ohp-page-hero h1{font-size:38px;line-height:1.1;margin:14px 0 10px}
.ohp-page-hero .lead{font-size:15.5px}
.ohp-page-hero .motif{position:relative;width:170px;height:150px;display:grid;place-items:center}

/* radar (epidemic) */
.ohp-radar{position:relative;width:150px;height:150px;border-radius:50%;overflow:hidden;
   border:1px solid color-mix(in srgb,var(--a1) 45%,transparent);
   background:radial-gradient(circle,color-mix(in srgb,var(--a1) 12%,transparent),transparent 70%)}
.ohp-radar .grid{position:absolute;inset:25%;border-radius:50%;border:1px solid color-mix(in srgb,var(--a1) 25%,transparent)}
.ohp-radar .sweep{position:absolute;inset:0;border-radius:50%;
   background:conic-gradient(from 0deg,color-mix(in srgb,var(--a1) 55%,transparent),transparent 28%);
   animation:ohp-rot 3.2s linear infinite}
.ohp-radar .ring{position:absolute;inset:0;border-radius:50%;border:2px solid var(--a2);
   animation:ohp-ping 3.2s ease-out infinite}
.ohp-radar .blip{position:absolute;width:8px;height:8px;border-radius:50%;background:var(--a2);
   box-shadow:0 0 12px var(--a2);animation:ohp-blink 1.6s ease-in-out infinite}

/* orbit (dashboard) */
.ohp-orbit{position:relative;width:150px;height:150px}
.ohp-orbit .globe{position:absolute;inset:34px;border-radius:50%;
   background:radial-gradient(circle at 35% 30%,color-mix(in srgb,var(--a1) 60%,transparent),#1E1B4B 72%);
   border:1px solid color-mix(in srgb,var(--a1) 50%,transparent);box-shadow:0 0 40px color-mix(in srgb,var(--a1) 35%,transparent)}
.ohp-orbit .globe::before,.ohp-orbit .globe::after{content:"";position:absolute;inset:0;border-radius:50%;
   border:1px solid rgba(255,255,255,.18)}
.ohp-orbit .globe::before{transform:scaleX(.45)} .ohp-orbit .globe::after{transform:scaleY(.45)}
.ohp-orbit .path{position:absolute;inset:0;border-radius:50%;border:1px dashed color-mix(in srgb,var(--a1) 45%,transparent);
   animation:ohp-rot 9s linear infinite}
.ohp-orbit .path.p2{inset:16px;animation-duration:5.5s;animation-direction:reverse}
.ohp-orbit .path::after{content:"";position:absolute;top:-5px;left:calc(50% - 5px);width:10px;height:10px;border-radius:50%;
   background:var(--a2);box-shadow:0 0 14px var(--a2)}

/* scanner (audit) */
.ohp-scan{position:relative;width:150px;height:130px;border-radius:14px;overflow:hidden;padding:14px;
   background:rgba(2,6,23,.6);border:1px solid var(--border);display:flex;flex-direction:column;gap:9px}
.ohp-scan .ln{height:7px;border-radius:4px;background:color-mix(in srgb,var(--a1) 35%,transparent)}
.ohp-scan .ln.hl{background:color-mix(in srgb,var(--a2) 70%,transparent)}
.ohp-scan .beam{position:absolute;left:0;right:0;height:26px;top:-26px;
   background:linear-gradient(180deg,transparent,color-mix(in srgb,var(--a2) 45%,transparent),transparent);
   animation:ohp-scanline 2.8s ease-in-out infinite}

/* heartbeat (status) */
.ohp-ecg{width:170px;height:80px}
.ohp-ecg path{fill:none;stroke:var(--a1);stroke-width:2.5;stroke-linejoin:round;stroke-linecap:round;
   stroke-dasharray:320;stroke-dashoffset:320;animation:ohp-draw 2.4s linear infinite;
   filter:drop-shadow(0 0 6px var(--a1))}

/* broadcast (feed) */
.ohp-cast{position:relative;width:150px;height:150px;display:grid;place-items:center}
.ohp-cast .core{width:16px;height:16px;border-radius:50%;background:var(--a1);
   box-shadow:0 0 18px var(--a1);z-index:1}
.ohp-cast .wave{position:absolute;inset:0;border-radius:50%;border:2px solid var(--a1);opacity:0;
   animation:ohp-ping 3s ease-out infinite}
.ohp-cast .wave:nth-child(2){animation-delay:-1s;border-color:var(--a2)}
.ohp-cast .wave:nth-child(3){animation-delay:-2s}

/* feed items */
.ohp-chip{display:inline-block;font:600 11px/1 var(--body);letter-spacing:.08em;text-transform:uppercase;
   padding:5px 9px;border-radius:999px;border:1px solid var(--border);color:var(--muted)}
.ohp-chip.new{color:#60A5FA;border-color:rgba(96,165,250,.45);background:rgba(96,165,250,.10)}
.ohp-chip.analyzed{color:#4ADE80;border-color:rgba(74,222,128,.45);background:rgba(74,222,128,.10)}
.ohp-chip.blocked,.ohp-chip.error{color:#FBBF24;border-color:rgba(251,191,36,.45);background:rgba(251,191,36,.10)}
.ohp-item-title{font:600 17px/1.35 var(--display);color:var(--text);margin:8px 0 4px}
.ohp-item-title a{color:inherit;text-decoration:none}
.ohp-item-title a:hover{color:var(--a1)}
.ohp-item-meta{color:var(--muted);font-size:13px}

/* status cards */
.ohp-svc{display:grid;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));gap:14px}
.ohp-svc .c{background:var(--surface);border:1px solid var(--border);border-radius:16px;padding:18px}
.ohp-svc .top{display:flex;align-items:center;gap:10px;margin-bottom:6px}
.ohp-svc .led{width:10px;height:10px;border-radius:50%}
.ohp-svc .ok .led{background:#4ADE80;box-shadow:0 0 12px #4ADE80;animation:ohp-blink 2s ease-in-out infinite}
.ohp-svc .warn .led{background:#FBBF24;box-shadow:0 0 12px #FBBF24}
.ohp-svc .info .led{background:#60A5FA;box-shadow:0 0 12px #60A5FA}
.ohp-svc b{font:600 16px var(--display)}
.ohp-svc .purpose{color:var(--muted);font-size:13px;margin-bottom:10px}
.ohp-svc .detail{font:500 13px var(--mono)}
.ohp-svc .ok .detail{color:#4ADE80} .ohp-svc .warn .detail{color:#FBBF24} .ohp-svc .info .detail{color:#60A5FA}

/* disclaimer + footer */
.ohp-note{border-radius:18px;padding:22px 26px;border:1px solid rgba(251,191,36,.3);
   background:linear-gradient(90deg,rgba(251,191,36,.08),rgba(251,191,36,.02));color:var(--text);font-size:14px}
.ohp-footer{margin-top:48px;padding:26px 0 8px;border-top:1px solid var(--border);display:flex;
   flex-wrap:wrap;justify-content:space-between;gap:12px;color:#64748B;font:13px var(--body)}
.ohp-footer b{color:var(--text);font-family:var(--display)}

/* ---------- keyframes ---------- */
@keyframes ohp-spin-y{to{transform:rotateY(360deg)}}
@keyframes ohp-rot{to{transform:rotate(360deg)}}
@keyframes ohp-ping{0%{transform:scale(.15);opacity:1}100%{transform:scale(1.05);opacity:0}}
@keyframes ohp-blink{0%,100%{opacity:1}50%{opacity:.25}}
@keyframes ohp-drift{to{transform:translate(-40px,30px) scale(1.12)}}
@keyframes ohp-rise{0%{transform:translateY(40px);opacity:0}15%{opacity:.75}85%{opacity:.5}100%{transform:translateY(-320px);opacity:0}}
@keyframes ohp-travel{0%{left:10%}50%{left:calc(90% - 12px)}100%{left:10%}}
@keyframes ohp-scanline{0%{top:-26px}100%{top:130px}}
@keyframes ohp-draw{to{stroke-dashoffset:0}}

/* ---------- responsive ---------- */
@media (max-width:900px){
  .ohp-hero{grid-template-columns:1fr;padding:34px 22px}
  .ohp-hero h1{font-size:38px}
  .ohp-hero .art{min-height:280px}
  .ohp-grid{grid-template-columns:1fr}
  .ohp-stats{grid-template-columns:repeat(2,1fr)}
  .ohp-flow{grid-template-columns:1fr;gap:22px}
  .ohp-flow .track,.ohp-flow .pulse{display:none}
  .ohp-page-hero{grid-template-columns:1fr;padding:26px 22px}
  .ohp-page-hero .motif{display:none}
}
@media (prefers-reduced-motion:reduce){
  .ohp *,.ohp *::before,.ohp *::after{animation:none!important;transition:none!important}
  .ohp-base{opacity:.4}
}
</style>
"""


def inject_global_css():
    st.html(GLOBAL_CSS)


# ---------------------------------------------------------------------
# Motifs
# ---------------------------------------------------------------------
def helix(rungs: int = 18, width: int = 230, period: float = 6.0, small: bool = False) -> str:
    """Rotating 3D double helix. Each rung spins on the Y axis with a phase
    offset, which produces the twisted-ladder illusion."""
    step = period / rungs * 1.1
    out = []
    for i in range(rungs):
        out.append(
            f'<div class="ohp-rung" style="--w:{width}px;animation-delay:-{i * step:.2f}s;'
            f'animation-duration:{period}s"><span class="n l"></span><span class="bar"></span>'
            f'<span class="n r"></span></div>')
    return f'<div class="ohp-helix{" sm" if small else ""}">{"".join(out)}</div>'


def _radar() -> str:
    blips = "".join(
        f'<span class="blip" style="top:{t}%;left:{l}%;animation-delay:-{d}s"></span>'
        for t, l, d in [(28, 62, 0.2), (58, 30, 0.9), (66, 70, 1.3)])
    return (f'<div class="ohp-radar"><span class="grid"></span><span class="sweep"></span>'
            f'<span class="ring"></span><span class="ring" style="animation-delay:-1.6s"></span>{blips}</div>')


def _orbit() -> str:
    return '<div class="ohp-orbit"><span class="path"></span><span class="path p2"></span><span class="globe"></span></div>'


def _scan() -> str:
    widths = [92, 70, 84, 55, 78, 64, 88]
    lines = "".join(f'<span class="ln{" hl" if i == 3 else ""}" style="width:{w}%"></span>'
                    for i, w in enumerate(widths))
    return f'<div class="ohp-scan">{lines}<span class="beam"></span></div>'


def _ecg() -> str:
    return ('<svg class="ohp-ecg" viewBox="0 0 170 80" aria-hidden="true"><path d="M0 44 H42 L52 44 L60 14 '
            'L70 70 L78 44 H104 L112 32 L120 44 H170"/></svg>')


def _cast() -> str:
    return ('<div class="ohp-cast"><span class="wave"></span><span class="wave"></span>'
            '<span class="wave"></span><span class="core"></span></div>')


def feed_item_header(title: str, url: str, published: str, source: str, status: str):
    """Title, link, date and status chip for one outbreak-feed item."""
    safe_url = _html.escape(url, quote=True) if url.startswith("https://") else "#"
    st.html(f'<div class="ohp t-feed"><span class="ohp-chip {_html.escape(status)}">{_html.escape(status)}</span>'
            f'<div class="ohp-item-title"><a href="{safe_url}" target="_blank" rel="noopener noreferrer">'
            f'{_html.escape(title)}</a></div>'
            f'<div class="ohp-item-meta">{_html.escape(source)} · published {_html.escape(published)}</div></div>')


MOTIFS = {
    "feed": _cast,
    "genomic": lambda: helix(rungs=11, width=120, period=5.0, small=True),
    "epidemic": _radar,
    "dashboard": _orbit,
    "audit": _scan,
    "status": _ecg,
}


# ---------------------------------------------------------------------
# Page building blocks
# ---------------------------------------------------------------------
def page_hero(theme: str, eyebrow: str, title: str, subtitle: str):
    """Coloured header band with an animated motif; one theme per page."""
    motif = MOTIFS.get(theme, lambda: "")()
    st.html(f"""
<section class="ohp ohp-page-hero t-{theme}">
  <div>
    <span class="eyebrow"><i></i>{_html.escape(eyebrow)}</span>
    <h1>{title}</h1>
    <p class="lead">{_html.escape(subtitle)}</p>
  </div>
  <div class="motif" aria-hidden="true">{motif}</div>
</section>""")


def section_heading(theme: str, eyebrow: str, title: str, subtitle: str = ""):
    sub = f'<div class="section-sub">{_html.escape(subtitle)}</div>' if subtitle else ""
    st.html(f'<div class="ohp t-{theme}" style="margin:46px 0 4px"><span class="eyebrow"><i></i>'
            f'{_html.escape(eyebrow)}</span><div class="section-title">{title}</div>{sub}</div>')


def footer(version: str):
    st.html(f"""
<footer class="ohp ohp-footer">
  <div><b>🧬 One Health Intelligence Platform</b> · v{_html.escape(version)}</div>
  <div>NCBI BLAST · Curated AMR markers · Claude API · HL7 v2.5</div>
  <div>⚕️ Research &amp; surveillance use only — not for clinical diagnosis</div>
</footer>""")


# ---------------------------------------------------------------------
# Home page sections
# ---------------------------------------------------------------------
_ICONS = {
    "dna": '<svg viewBox="0 0 24 24"><path d="M7 3c0 6 10 6 10 12s-10 6-10 6M17 3c0 6-10 6-10 12s10 6 10 6M8.5 6h7M8.5 18h7M10 9.5h4M10 14.5h4"/></svg>',
    "pulse": '<svg viewBox="0 0 24 24"><path d="M3 12h4l2-6 4 12 2-6h6"/></svg>',
    "globe": '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="9"/><path d="M3 12h18M12 3c3 3.5 3 14.5 0 18M12 3c-3 3.5-3 14.5 0 18"/></svg>',
}


def home_hero():
    bases = "".join(
        f'<span class="ohp-base" style="left:{l}%;bottom:{b}%;animation-delay:-{d}s">{ch}</span>'
        for ch, l, b, d in [("A", 8, 0, 0), ("T", 86, 4, 1.4), ("G", 18, 10, 2.8), ("C", 76, 0, 4.1),
                            ("G", 92, 14, 5.3), ("A", 4, 22, 3.5), ("T", 60, 0, 6.2), ("C", 34, 6, 0.8)])
    st.html(f"""
<section class="ohp ohp-hero t-home">
  <span class="ohp-blob b1"></span><span class="ohp-blob b2"></span>
  <div style="position:relative;z-index:1">
    <span class="eyebrow"><i></i>One Health · Genomics × Epidemiology</span>
    <h1>Read the genome.<br><span class="grad">Catch the outbreak.</span></h1>
    <p class="lead">One workspace that screens bacterial DNA for antimicrobial-resistance genes,
    checks health reports for misinformation, scores outbreak risk, and turns every finding
    into a standards-ready HL7 report.</p>
    <div class="chips">
      <span class="chip">🧬 AMR gene detection</span><span class="chip">🦠 Outbreak early warning</span>
      <span class="chip">🛡️ Fake-news screening</span><span class="chip">📋 HL7 v2.5 output</span>
    </div>
  </div>
  <div class="art" aria-hidden="true">{bases}{helix(rungs=20, width=240)}</div>
</section>""")


def home_stats():
    st.html("""
<section class="ohp ohp-stats" style="margin-top:22px">
  <div class="ohp-stat"><b>NCBI-curated</b><span>resistance genes called with AMRFinderPlus, in seconds</span></div>
  <div class="ohp-stat"><b>WHO feed</b><span>Disease Outbreak News pulled in automatically</span></div>
  <div class="ohp-stat"><b>HL7 v2.5</b><span>ORU^R01 lab reports &amp; epidemic alerts</span></div>
  <div class="ohp-stat"><b>Saved</b><span>every analysis is stored in a permanent audit trail</span></div>
</section>""")


def home_features():
    cards = [
        ("genomic", "dna", "Genomic AMR Surveillance",
         "Paste a DNA or protein sequence, or upload FASTA files, and find out whether it carries a resistance gene.",
         ["NCBI AMRFinderPlus with its curated reference database",
          "Finds every resistance gene in the sequence, not just one",
          "Flags last-resort threats: carbapenem, colistin, vancomycin",
          "Downloadable HL7 lab report per sample"]),
        ("epidemic", "pulse", "Epidemic Early Warning",
         "Paste a news story, WHO alert or hospital report and get a structured outbreak risk assessment.",
         ["Watches WHO Disease Outbreak News automatically",
          "Blocks fake or sensational articles first",
          "Claude extracts disease, cases, trend & location",
          "Risk score, response level and HL7 epidemic alert"]),
        ("dashboard", "globe", "One Health Correlation",
         "See genomic and epidemic signals side by side, on a world map and over time.",
         ["Live map of geocoded outbreak reports",
          "Confidence timeline across all analyses",
          "Saved audit trail with real vs demo labels",
          "Reopen any past report; one-click CSV export"]),
    ]
    items = []
    for theme, icon, title, text, bullets in cards:
        lis = "".join(f"<li>{_html.escape(b)}</li>" for b in bullets)
        items.append(f'<article class="ohp-card t-{theme}"><div class="ic">{_ICONS[icon]}</div>'
                     f'<h3>{title}</h3><p>{_html.escape(text)}</p><ul>{lis}</ul></article>')
    st.html(f'<section class="ohp ohp-grid">{"".join(items)}</section>')


def home_pipeline():
    steps = [
        ("1", "Input", "A sequence, FASTA batch, or a health article"),
        ("2", "Validate", "Clean & check the sequence; screen the article for misinformation"),
        ("3", "Search", "AMRFinderPlus for genes · Claude for outbreak signals"),
        ("4", "Classify & score", "Match resistance markers, weigh evidence, rate risk"),
        ("5", "Report", "HL7 v2.5 message, audit log entry and map point"),
    ]
    html_steps = "".join(f'<div class="ohp-step"><div class="dot">{n}</div><h4>{t}</h4><p>{_html.escape(d)}</p></div>'
                         for n, t, d in steps)
    st.html(f'<section class="ohp ohp-flow t-home"><span class="track"></span><span class="pulse"></span>'
            f'{html_steps}</section>')


def home_disclaimer():
    st.html("""
<div class="ohp ohp-note" style="margin-top:40px">
  <b>⚕️ Research &amp; surveillance use only.</b> This platform is a decision-support prototype, not a validated
  diagnostic or official public-health reporting system. Confirm findings with accredited laboratory methods
  and official sources before acting.
</div>""")


def status_cards(results):
    cards = []
    for name, purpose, level, detail in results:
        cards.append(f'<div class="c {level}"><div class="top"><span class="led"></span><b>{_html.escape(name)}</b></div>'
                     f'<div class="purpose">{_html.escape(purpose)}</div>'
                     f'<div class="detail">{_html.escape(detail)}</div></div>')
    st.html(f'<section class="ohp ohp-svc">{"".join(cards)}</section>')
