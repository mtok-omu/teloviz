"""Standalone HTML telomere-cap report for teloviz (spec section 10).

Writes one self-contained HTML file: a header echoing the settings and the exact
command used, a summary count, and a per-chromosome table of end-region counts
and calls. Pure string building — no extra dependencies, no browser needed.
"""

from __future__ import annotations

import html
from pathlib import Path

from .calling import Call

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
"""


def _cell_call(ok: bool) -> str:
    return '<td class="yes">✓</td>' if ok else '<td class="no">·</td>'


def _status_class(status: str) -> str:
    return {"both": "both", "none": "none"}.get(status, "one")


def write_report(path: str | Path, calls: list[Call], *, meta: dict) -> Path:
    """Render the report to ``path``; return the Path written.

    ``meta`` echoes the run: keys ``command``, ``input``, ``mode``, ``dist_kb``,
    ``call_min``, ``min_count``, ``motif``, ``window_size``.
    """
    e = html.escape
    n = len(calls)
    n_both = sum(c.both for c in calls)
    n_five = sum(c.five for c in calls)
    n_three = sum(c.three for c in calls)
    n_none = sum(c.status == "none" for c in calls)

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
    settings_html = "\n".join(
        f"<dt>{k}</dt><dd>{v}</dd>" for k, v in settings
    )

    rows = []
    for c in calls:
        rows.append(
            f'<tr><td class="l">{e(c.id)}{" *" if c.both else ""}</td>'
            f"<td>{c.length / 1e6:.2f}</td>"
            f"<td>{c.five_count}</td>{_cell_call(c.five)}"
            f"<td>{c.three_count}</td>{_cell_call(c.three)}"
            f'<td class="{_status_class(c.status)}">{e(c.status)}</td></tr>'
        )
    rows_html = "\n".join(rows)

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>teloviz telomere report</title><style>{_CSS}</style></head>
<body>
<h1>teloviz — telomere-cap report</h1>
<p class="sub">{n} sequences &middot; both ends capped: <b>{n_both}</b> &middot;
 5' capped: {n_five} &middot; 3' capped: {n_three} &middot; no telomere: {n_none}
 <span style="color:#999">(&quot;*&quot; = both ends)</span></p>
<dl class="settings">
{settings_html}
</dl>
<table>
<thead><tr>
 <th class="l">chromosome</th><th>length (Mb)</th>
 <th>5' count</th><th>5' call</th>
 <th>3' count</th><th>3' call</th><th>status</th>
</tr></thead>
<tbody>
{rows_html}
</tbody></table>
</body></html>
"""
    p = Path(path)
    if p.parent != Path(""):
        p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(doc, encoding="utf-8")
    return p
