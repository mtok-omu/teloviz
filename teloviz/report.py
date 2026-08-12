"""Standalone HTML telomere-cap report for teloviz (spec section 10).

Writes one self-contained HTML file: a header echoing the settings and the exact
command used, a summary count, and a per-chromosome table of end-region counts
and calls. Pure string building — no extra dependencies, no browser needed.
"""

from __future__ import annotations

import html
from pathlib import Path

from .calling import SHORT_FACTOR, Call
from .features import FeatureSet

_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
 margin:2rem;color:#1a1a1a;line-height:1.45}
h1{font-size:1.4rem;margin:0 0 .3rem} .sub{color:#666;margin:0 0 1.2rem}
table{border-collapse:collapse;font-size:.9rem} th,td{padding:.35rem .7rem;
 border-bottom:1px solid #e2e2e2;text-align:right} th{text-align:right;
 border-bottom:2px solid #bbb} td.l,th.l{text-align:left}
code{background:#f4f4f4;padding:.1rem .35rem;border-radius:3px}
.settings{background:#f7f7f7;border:1px solid #e2e2e2;border-radius:6px;
 padding:.8rem 1rem;margin:0 0 1.2rem;font-size:.9rem}
.settings dt{font-weight:600;color:#444} .settings dd{margin:0 0 .4rem}
.yes{color:#0a7a2f;font-weight:600} .no{color:#b00;font-weight:600}
.both{color:#0a7a2f;font-weight:600} .one{color:#b8860b;font-weight:600}
.none{color:#b00;font-weight:600}
td.note{color:#444;font-size:.82rem}
.warn{color:#b8860b;cursor:help} .foot{color:#666;font-size:.82rem;margin-top:.8rem}
"""


def _cell_call(ok: bool) -> str:
    return '<td class="yes">✓</td>' if ok else '<td class="no">·</td>'


def _end_note(features: FeatureSet, cid: str, length: int, end: str,
              proximity_bp: int) -> str:
    """Modifier text for one *uncapped* end (spec section 5.1).

    Names the nearest feature within ``proximity_bp`` of that end (an rDNA/NOR
    array is the reason a real telomere may be unreachable), or states plainly
    that no feature is nearby — the absence is itself evidence for a genuine gap.
    ``end`` is ``"5'"`` (distance measured from position 0) or ``"3'"`` (from the
    chromosome length).
    """
    best = None  # (distance_bp, feature)
    for f in features.for_chrom(cid):
        dist = max(0, f.start) if end == "5'" else max(0, length - f.end)
        if dist <= proximity_bp and (best is None or dist < best[0]):
            best = (dist, f)
    if best is None:
        return f"cap missing — no feature within {proximity_bp / 1000:g} kb"
    dist, f = best
    return (f"cap missing — {f.type} ({html.escape(f.name)}) "
            f"{dist / 1000:g} kb from the {end} end")


def _status_class(status: str) -> str:
    return {"both": "both", "none": "none"}.get(status, "one")


def write_report(path: str | Path, calls: list[Call], *, meta: dict,
                 features: FeatureSet | None = None,
                 proximity_bp: int = 200_000) -> Path:
    """Render the report to ``path``; return the Path written.

    ``meta`` echoes the run: keys ``command``, ``input``, ``mode``, ``dist_kb``,
    ``call_min``, ``min_count``, ``motif``, ``window_size`` (and, with features,
    ``features``). When ``features`` is given, uncapped ends gain a modifier note
    naming a nearby feature or stating that none is within ``proximity_bp``; the
    report is otherwise unchanged (subject stays telomere, no feature statistics).
    """
    e = html.escape
    n = len(calls)
    n_both = sum(c.both for c in calls)
    n_five = sum(c.five for c in calls)
    n_three = sum(c.three for c in calls)
    n_none = sum(c.status == "none" for c in calls)
    has_feat = features is not None

    settings = [
        ("Input", f"<code>{e(str(meta.get('input', '')))}</code>"),
        ("Command", f"<code>{e(str(meta.get('command', '')))}</code>"),
        ("Call distance from each end", f"{e(str(meta.get('dist_kb')))} kb"),
        ("Call threshold (forward+reverse in that region)", f"&ge; {e(str(meta.get('call_min')))}"),
        ("Noise floor (--min-count)", e(str(meta.get("min_count")))),
        ("Motif", e(str(meta.get("motif") or "all"))),
        ("tidk window size", f"{e(str(meta.get('window_size')))} bp"),
        ("Plot mode(s)", e(str(meta.get("mode")))),
    ]
    if has_feat:
        settings.append(("Feature BED", f"<code>{e(str(meta.get('features', '')))}</code>"))
        settings.append(("Feature proximity (annotation only)", f"{proximity_bp / 1000:g} kb"))
    settings_html = "\n".join(
        f"<dt>{k}</dt><dd>{v}</dd>" for k, v in settings
    )

    # A sequence shorter than 2x the call distance has its 5' and 3' end regions
    # overlap, so the same windows are counted for both ends and a "both" call can
    # be spurious. ``Call.short`` decides which rows are flagged (same flag the
    # plot ambers); meta only supplies the threshold in kb for the wording.
    try:
        dist_bp = round(float(meta.get("dist_kb", 0)) * 1000)
    except (TypeError, ValueError):
        dist_bp = 0
    short_thresh = SHORT_FACTOR * dist_bp
    thresh_txt = (f" ({short_thresh / 1000:g} kb)" if dist_bp else "")
    warn_title = (
        f"length &lt; 2&times; call distance{thresh_txt}: "
        "the 5' and 3' end regions overlap, so a &ldquo;both ends&rdquo; call "
        "may be unreliable"
    )

    rows = []
    any_short = False
    for c in calls:
        note_cell = ""
        if has_feat:
            notes = []
            if not c.five:
                notes.append("5' " + _end_note(features, c.id, c.length, "5'", proximity_bp))
            if not c.three:
                notes.append("3' " + _end_note(features, c.id, c.length, "3'", proximity_bp))
            note_cell = f'<td class="l note">{"<br>".join(notes)}</td>'
        any_short = any_short or c.short
        warn = f' <span class="warn" title="{warn_title}">&#9888;</span>' if c.short else ""
        rows.append(
            f'<tr><td class="l">{e(c.id)}{c.suffix}</td>'
            f"<td>{c.length / 1e6:.2f}</td>"
            f"<td>{c.five_count}</td>{_cell_call(c.five)}"
            f"<td>{c.three_count}</td>{_cell_call(c.three)}"
            f'<td class="{_status_class(c.status)}">{e(c.status)}{warn}</td>{note_cell}</tr>'
        )
    rows_html = "\n".join(rows)
    foot_html = (
        f'\n<p class="foot">&#9888; sequence shorter than 2&times; the call '
        f"distance ({short_thresh / 1000:g} kb): its 5' and 3' end regions "
        "overlap, so the same windows are counted for both ends and a "
        "&ldquo;both ends&rdquo; call can be spurious. teloviz assumes "
        "chromosome-level input.</p>"
        if any_short else ""
    )
    star_note = "; &quot;(*)&quot; = flagged below" if any_short else ""
    feat_th = '<th class="l">feature notes</th>' if has_feat else ""

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>teloviz telomere report</title><style>{_CSS}</style></head>
<body>
<h1>teloviz — telomere-cap report</h1>
<p class="sub">{n} sequences &middot; both ends capped: <b>{n_both}</b> &middot;
 5' capped: {n_five} &middot; 3' capped: {n_three} &middot; no telomere: {n_none}
 <span style="color:#999">(&quot;*&quot; = both ends{star_note})</span></p>
<dl class="settings">
{settings_html}
</dl>
<table>
<thead><tr>
 <th class="l">chromosome</th><th>length (Mb)</th>
 <th>5' count</th><th>5' call</th>
 <th>3' count</th><th>3' call</th><th>status</th>{feat_th}
</tr></thead>
<tbody>
{rows_html}
</tbody></table>{foot_html}
</body></html>
"""
    p = Path(path)
    if p.parent != Path(""):
        p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(doc, encoding="utf-8")
    return p
