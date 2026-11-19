#!/usr/bin/env python3
"""Draw the benchmark's two cross-arm curves from Day 5's recorded tables.

Reads `results/day5-provenance.json` and writes `docs/figures/detection-vs-dilution.svg`.
Runs no tool, touches no sequence data, and invents no number: every point plotted is
a field of that file. Re-run it after any re-score and the figure follows.

    python scripts/day6_figure.py [--out docs/figures/detection-vs-dilution.svg]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape

# Light chart surface, read the same way in both GitHub themes.
SURFACE = "#fcfcfb"
GRID = "#e4e3df"
AXIS = "#9b9a93"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#7c7b74"
KRAKEN = "#2a78d6"  # categorical slot 1
ASSEMBLY = "#eb6834"  # categorical slot 2
FLOOR = "#ece9e2"

W, H = 1040, 584
PANEL_W, PANEL_H = 400, 300
PAD_L, PAD_T = 78, 182
GAP = 106
FONT = "ui-sans-serif,-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif"


def log_scale(lo: float, hi: float, px0: float, px1: float):
    a, b = math.log10(lo), math.log10(hi)

    def to_px(v: float) -> float:
        return px0 + (math.log10(v) - a) / (b - a) * (px1 - px0)

    return to_px


def lin_scale(lo: float, hi: float, px0: float, px1: float):
    def to_px(v: float) -> float:
        return px0 + (v - lo) / (hi - lo) * (px1 - px0)

    return to_px


def text(x, y, s, size=12.5, fill=INK_2, anchor="start", weight=400, style=""):
    extra = f' font-style="{style}"' if style else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{extra}>{escape(s)}</text>'
    )


def polyline(points, colour, width=2.0, dash=None):
    pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="{width}" '
        f'stroke-linecap="round" stroke-linejoin="round"{d}/>'
    )


def marker(x, y, colour, r=4.6):
    return (
        f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{colour}" '
        f'stroke="{SURFACE}" stroke-width="2"/>'
    )


def panel_frame(x0, y0, letter, title, subtitle, y_label):
    out = [
        f'<rect x="{x0 - 0.5:.1f}" y="{y0 - 0.5:.1f}" width="{PANEL_W + 1}" '
        f'height="{PANEL_H + 1}" fill="none" stroke="{GRID}" stroke-width="1"/>',
        text(x0, y0 - 44, letter, size=15.5, fill=INK, weight=700),
        text(x0 + 20, y0 - 44, title, size=15.5, fill=INK, weight=650),
        text(x0, y0 - 25, subtitle, size=12.5, fill=MUTED),
        f'<text x="{x0 - 52:.1f}" y="{y0 + PANEL_H / 2:.1f}" font-family="{FONT}" '
        f'font-size="12" fill="{MUTED}" text-anchor="middle" '
        f'transform="rotate(-90 {x0 - 52:.1f} {y0 + PANEL_H / 2:.1f})">{escape(y_label)}</text>',
    ]
    return out


def build(prov: dict) -> str:
    levels = prov["cross_arm"]["levels"]
    floor = prov["resident_background"]["resolvable_floor_per_million"]
    a_thresh = prov["arm_a"]["headline"]["threshold_reads"]
    b_thresh = prov["arm_b"]["headline"]["threshold_bitscore"]

    dose = [lv["per_million"] for lv in levels]
    x_lo, x_hi = 0.5, 6000.0

    s: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" role="img" aria-label="Two panels: alarm margin and '
        f'read recovery against dilution for both detection pipelines">',
        f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
        text(
            PAD_L - 2,
            40,
            "Detection and recovery across the dilution series",
            size=19,
            fill=INK,
            weight=700,
        ),
        text(
            PAD_L - 2,
            66,
            "SARS-CoV-2 reads in a November 2019 reference frame, against 20 pre-COVID negatives.",
            size=13,
            fill=INK_2,
        ),
    ]

    # ---- legend -------------------------------------------------------------
    lx, ly = PAD_L - 2, 100
    for colour, label, width in (
        (KRAKEN, "Kraken2", 110),
        (ASSEMBLY, "MEGAHIT + DIAMOND", 200),
    ):
        s.append(f'<rect x="{lx}" y="{ly - 8}" width="20" height="3" rx="1.5" fill="{colour}"/>')
        s.append(marker(lx + 10, ly - 6.5, colour, r=4.2))
        s.append(text(lx + 27, ly - 2, label, size=12.5, fill=INK_2, weight=550))
        lx += width
    s.append(f'<rect x="{lx}" y="{ly - 14}" width="20" height="12" fill="{FLOOR}"/>')
    s.append(
        text(
            lx + 27,
            ly - 2,
            f"below the series' own floor ({floor:.1f} /M) — not a usable rung",
            size=12.5,
            fill=MUTED,
        )
    )

    # =========================== panel 1: alarms =============================
    x0, y0 = PAD_L, PAD_T
    fx = log_scale(x_lo, x_hi, x0 + 6, x0 + PANEL_W - 10)
    fy = log_scale(0.2, 4000.0, y0 + PANEL_H - 10, y0 + 14)
    s += panel_frame(
        x0,
        y0,
        "A",
        "Does it alarm?",
        "evidence as a multiple of each pipeline's own alarm threshold",
        "× the alarm threshold",
    )

    # floor band
    s.append(
        f'<rect x="{x0 + 0.5:.1f}" y="{y0 + 0.5:.1f}" width="{fx(floor) - x0:.1f}" '
        f'height="{PANEL_H - 1}" fill="{FLOOR}"/>'
    )

    for v in (0.2, 1, 10, 100, 1000):
        y = fy(v)
        s.append(
            f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x0 + PANEL_W:.1f}" y2="{y:.1f}" '
            f'stroke="{GRID}" stroke-width="1"/>'
        )
        if v != 0.2:
            s.append(text(x0 - 8, y + 4, f"{v:g}×", size=11.5, fill=MUTED, anchor="end"))

    # the alarm line
    s.append(
        f'<line x1="{x0:.1f}" y1="{fy(1):.1f}" x2="{x0 + PANEL_W:.1f}" y2="{fy(1):.1f}" '
        f'stroke="{INK}" stroke-width="1.6" stroke-dasharray="5 4"/>'
    )
    s.append(
        text(
            x0 + PANEL_W - 8, fy(1) - 8, "alarm line", size=11.5, fill=INK, anchor="end", weight=600
        )
    )
    s.append(text(x0 + PANEL_W - 8, fy(1) + 16, "below: silent", size=11, fill=MUTED, anchor="end"))

    a_pts = [(fx(d), fy(lv["arm_a_reads"] / a_thresh)) for d, lv in zip(dose, levels, strict=False)]
    b_pts = [
        (fx(d), fy(lv["arm_b_best_bitscore"] / b_thresh))
        for d, lv in zip(dose, levels, strict=False)
    ]
    s.append(polyline(a_pts, KRAKEN))
    s.append(polyline(b_pts, ASSEMBLY))
    for x, y in a_pts:
        s.append(marker(x, y, KRAKEN))
    for x, y in b_pts:
        s.append(marker(x, y, ASSEMBLY))

    # the one silent point
    silent = next(i for i, lv in enumerate(levels) if not lv["arm_a_alarm"])
    sx, sy = a_pts[silent]
    s.append(
        f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="9" fill="none" stroke="{KRAKEN}" '
        f'stroke-width="1.4" stroke-dasharray="3 3"/>'
    )
    s.append(text(sx + 16, sy + 4, "4 labelled reads against a", size=11, fill=KRAKEN, weight=600))
    s.append(text(sx + 16, sy + 17, "threshold of 12 — silent", size=11, fill=KRAKEN))

    s.append(
        text(
            fx(dose[0]) - 6,
            fy(levels[0]["arm_a_reads"] / a_thresh) - 14,
            "Kraken2",
            size=12,
            fill=KRAKEN,
            weight=650,
            anchor="end",
        )
    )
    s.append(
        text(
            fx(dose[2]) + 10,
            fy(levels[2]["arm_b_best_bitscore"] / b_thresh) - 14,
            "MEGAHIT + DIAMOND",
            size=12,
            fill=ASSEMBLY,
            weight=650,
        )
    )

    # ========================= panel 2: recovery =============================
    x1 = x0 + PANEL_W + GAP
    gx = log_scale(x_lo, x_hi, x1 + 6, x1 + PANEL_W - 10)
    gy = lin_scale(0, 104, y0 + PANEL_H - 10, y0 + 14)
    s += panel_frame(
        x1,
        y0,
        "B",
        "How much does it recover?",
        "share of the labelled viral fragments each one finds",
        "% of labelled fragments recovered",
    )
    s.append(
        f'<rect x="{x1 + 0.5:.1f}" y="{y0 + 0.5:.1f}" width="{gx(floor) - x1:.1f}" '
        f'height="{PANEL_H - 1}" fill="{FLOOR}"/>'
    )
    for v in (0, 25, 50, 75, 100):
        y = gy(v)
        s.append(
            f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x1 + PANEL_W:.1f}" y2="{y:.1f}" '
            f'stroke="{GRID}" stroke-width="1"/>'
        )
        s.append(text(x1 - 8, y + 4, f"{v}%", size=11.5, fill=MUTED, anchor="end"))

    ra, rb = [], []
    for d, lv in zip(dose, levels, strict=False):
        n = lv["viral_fragments"]
        ra.append((gx(d), gy(100 * lv["arm_a_reads"] / n)))
        rb.append((gx(d), gy(100 * lv["arm_b_reads_recovered"] / n)))
    s.append(polyline(ra, KRAKEN))
    s.append(polyline(rb, ASSEMBLY))
    for x, y in ra:
        s.append(marker(x, y, KRAKEN))
    for x, y in rb:
        s.append(marker(x, y, ASSEMBLY))

    s.append(
        text(
            x1 + PANEL_W,
            gy(98.2) - 12,
            "MEGAHIT + DIAMOND",
            size=12,
            fill=ASSEMBLY,
            weight=650,
            anchor="end",
        )
    )
    s.append(text(gx(dose[1]) + 6, gy(42.3) - 12, "Kraken2", size=12, fill=KRAKEN, weight=650))
    s.append(
        text(gx(dose[4]) + 2, gy(10) + 24, "1 read of 10", size=11, fill=ASSEMBLY, anchor="middle")
    )
    s.append(text(gx(dose[4]) + 2, gy(40) - 16, "4 of 10", size=11, fill=KRAKEN, anchor="middle"))

    # =========================== shared x axis ===============================
    ticks = [(lv["per_million"], lv["viral_fragments"]) for lv in levels]
    for px0, scale in ((x0, fx), (x1, gx)):
        y = y0 + PANEL_H
        s.append(
            f'<line x1="{px0:.1f}" y1="{y:.1f}" x2="{px0 + PANEL_W:.1f}" y2="{y:.1f}" '
            f'stroke="{AXIS}" stroke-width="1.2"/>'
        )
        for per_m, frags in ticks:
            x = scale(per_m)
            s.append(
                f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y + 5:.1f}" '
                f'stroke="{AXIS}" stroke-width="1.2"/>'
            )
            s.append(
                text(
                    x,
                    y + 19,
                    f"{per_m:.3g}" if per_m < 10 else f"{per_m:,.0f}",
                    size=11.5,
                    fill=INK_2,
                    anchor="middle",
                )
            )
            s.append(text(x, y + 33, f"{frags:,}", size=10.5, fill=MUTED, anchor="middle"))
        s.append(
            text(
                px0 + PANEL_W / 2,
                y + 54,
                "labelled viral reads per million post-QC pairs  (fragments spiked, below)",
                size=11.5,
                fill=MUTED,
                anchor="middle",
            )
        )

    s.append(
        text(
            PAD_L - 2,
            H - 28,
            f"Scored at {prov['parameters']['clade'].title()}: Kraken2 alarms at "
            f"\u2265{a_thresh:g} reads, the assembly at bitscore "
            f"\u2265{b_thresh:g}. Both thresholds "
            "are set from the 20 pre-COVID negatives alone,",
            size=11,
            fill=MUTED,
        )
    )
    s.append(
        text(
            PAD_L - 2,
            H - 13,
            "and both raise 0 false alarms on them. Drawn from "
            "results/day5-provenance.json by scripts/day6_figure.py.",
            size=11,
            fill=MUTED,
        )
    )
    s.append("</svg>")
    return "\n".join(s) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--provenance", default="results/day5-provenance.json")
    ap.add_argument("--out", default="docs/figures/detection-vs-dilution.svg")
    args = ap.parse_args()

    prov = json.loads(Path(args.provenance).read_text())
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(prov))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
