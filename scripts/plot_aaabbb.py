#!/usr/bin/env python3
"""
Generate publication-grade vector SVG figures for the AAAbbb Shuttle-Exchange confirmation dossier:
  1. assets/aaabbb_stability_island.svg      — 2D mass-space phase diagram (m1=1, m2, m3)
  2. assets/aaabbb_corotating_vs_inertial.svg — Inertial vs. Co-Rotating (m2-m3 binary) 3+3 shuttle frame

Pure Python standard library (zero matplotlib/numpy dependency).
"""
import csv
import json
import math
import os


def viridis_like(t):
    """Perceptually clear teal -> emerald -> gold -> warm orange colormap for t in [0, 1]."""
    t = max(0.0, min(1.0, t))
    stops = [
        (0.00, (186, 230, 253)),  # light sky
        (0.25, (56, 189, 248)),   # bright cyan
        (0.50, (16, 185, 129)),   # emerald
        (0.75, (245, 158, 11)),   # amber gold
        (1.00, (220, 38, 38)),    # deep scarlet peak
    ]
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            u = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            r = int(round(c0[0] + u * (c1[0] - c0[0])))
            g = int(round(c0[1] + u * (c1[1] - c0[1])))
            b = int(round(c0[2] + u * (c1[2] - c0[2])))
            return f"#{r:02x}{g:02x}{b:02x}"
    return "#dc2626"


def build_stability_island_svg(csv_path, out_path):
    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            rows.append({
                "m2": float(r["m2"]),
                "m3": float(r["m3"]),
                "regime": r["regime"],
                "marg3": float(r["stab_margin_3d"]),
                "rmin": float(r["rmin"]),
            })
    # Determine grid spacing
    m2s = sorted({round(r["m2"], 4) for r in rows})
    step = round(min(m2s[i + 1] - m2s[i] for i in range(len(m2s) - 1)), 4)

    # Display window: [0.54, 1.06] x [0.54, 1.06]
    lo, hi = 0.54, 1.06
    W, H = 840, 680
    ml, mr, mt, mb = 88, 175, 72, 82
    pw, ph = W - ml - mr, H - mt - mb

    def sx(m2):
        return ml + (m2 - lo) / (hi - lo) * pw

    def sy(m3):
        return mt + ph - (m3 - lo) / (hi - lo) * ph

    out = []
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Inter, -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif">')
    out.append('<defs>')
    out.append('<pattern id="kreinHatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">')
    out.append('<rect width="6" height="6" fill="#fee2e2"/>')
    out.append('<line x1="0" y1="0" x2="0" y2="6" stroke="#fca5a5" stroke-width="1.6"/>')
    out.append('</pattern>')
    out.append('<linearGradient id="cbarGrad" x1="0" y1="1" x2="0" y2="0">')
    for k in range(11):
        t = k / 10.0
        out.append(f'<stop offset="{int(t*100)}%" stop-color="{viridis_like(t)}"/>')
    out.append('</linearGradient>')
    out.append('<clipPath id="plotClip"><rect x="{ml}" y="{mt}" width="{pw}" height="{ph}"/></clipPath>'.format(ml=ml, mt=mt, pw=pw, ph=ph))
    out.append('</defs>')

    # Background card
    out.append(f'<rect x="0" y="0" width="{W}" height="{H}" rx="12" fill="#ffffff" stroke="#e2e8f0" stroke-width="1.5"/>')

    # Title & subtitle
    out.append(f'<text x="{ml}" y="32" font-size="18" font-weight="700" fill="#0f172a">AAAbbb Shuttle-Exchange: 2D Mass-Space Full-3D Stability Island (m₁ = 1)</text>')
    out.append(f'<text x="{ml}" y="53" font-size="12.5" fill="#475569">Exact inertial orbits (α = 0, q = 1) · 80-bit &amp; 128-bit analytical variational monodromy · Symmetric under m₂ ↔ m₃</text>')

    # Plot background
    out.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.2"/>')

    # Grid cells (clipped)
    out.append('<g clip-path="url(#plotClip)">')
    cw = (step / (hi - lo)) * pw + 0.6
    ch = (step / (hi - lo)) * ph + 0.6
    max_marg = 1.60
    for r in rows:
        if r["m2"] < lo - step or r["m2"] > hi + step or r["m3"] < lo - step or r["m3"] > hi + step:
            continue
        cx = sx(r["m2"]) - cw / 2.0
        cy = sy(r["m3"]) - ch / 2.0
        if r["regime"] == "S3":
            col = viridis_like(r["marg3"] / max_marg)
            out.append(f'<rect x="{cx:.2f}" y="{cy:.2f}" width="{cw:.2f}" height="{ch:.2f}" fill="{col}" stroke="{col}" stroke-width="0.3"/>')
        elif r["regime"] == "K":
            out.append(f'<rect x="{cx:.2f}" y="{cy:.2f}" width="{cw:.2f}" height="{ch:.2f}" fill="url(#kreinHatch)"/>')
        else:
            out.append(f'<rect x="{cx:.2f}" y="{cy:.2f}" width="{cw:.2f}" height="{ch:.2f}" fill="#e2e8f0" fill-opacity="0.55"/>')

    # Symmetry diagonal m2 = m3
    out.append(f'<line x1="{sx(lo):.2f}" y1="{sy(lo):.2f}" x2="{sx(hi):.2f}" y2="{sy(hi):.2f}" stroke="#0f172a" stroke-width="1.2" stroke-dasharray="5,4" opacity="0.55"/>')
    # Equal-mass lines m2 = 1 and m3 = 1
    out.append(f'<line x1="{sx(1.0):.2f}" y1="{sy(lo):.2f}" x2="{sx(1.0):.2f}" y2="{sy(hi):.2f}" stroke="#64748b" stroke-width="1" stroke-dasharray="3,3" opacity="0.6"/>')
    out.append(f'<line x1="{sx(lo):.2f}" y1="{sy(1.0):.2f}" x2="{sx(hi):.2f}" y2="{sy(1.0):.2f}" stroke="#64748b" stroke-width="1" stroke-dasharray="3,3" opacity="0.6"/>')

    # Region text callouts
    out.append(f'<text x="{sx(0.62):.1f}" y="{sy(0.60):.1f}" font-size="11.5" font-weight="600" fill="#475569">Hyperbolic (H)</text>')
    out.append(f'<text x="{sx(0.60):.1f}" y="{sy(0.575):.1f}" font-size="10.5" fill="#64748b">k₁ &gt; +2 tangent exit</text>')
    out.append(f'<text x="{sx(0.975):.1f}" y="{sy(0.79):.1f}" font-size="11" font-weight="600" fill="#991b1b" text-anchor="middle">Krein (K)</text>')
    out.append(f'<text x="{sx(0.79):.1f}" y="{sy(0.985):.1f}" font-size="11" font-weight="600" fill="#991b1b" text-anchor="middle">Krein (K)</text>')

    # Key showcase markers
    markers = [
        (0.90, 0.80, "2D-DEM-002 (1.206)", "star", "#0f172a", "#fde047", -95, -24),
        (0.94, 0.71, "PEAK [1,0.94,0.71] (1.596)", "diamond", "#0f172a", "#ffffff", -148, 28),
        (0.89, 0.89, "ISO89 (1.091)", "circle", "#0f172a", "#ffffff", -102, -16),
        (0.70, 0.70, "WIDE (rmin=0.243)", "circle", "#0f172a", "#ffffff", -22, -22),
        (0.80, 0.60, "2D-DEM-019 (0.160)", "circle", "#0f172a", "#ffffff", 16, 16),
        (1.00, 1.00, "Equal-Mass [1,1,1] (UNSTABLE k₁=3.57)", "cross", "#991b1b", "#ef4444", -212, -18),
    ]
    for m2, m3, lbl, kind, stroke_c, fill_c, dx, dy in markers:
        x, y = sx(m2), sy(m3)
        tx, ty = x + dx, y + dy
        out.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{tx + (60 if dx < 0 else 0):.1f}" y2="{ty - 3:.1f}" stroke="#0f172a" stroke-width="0.9" opacity="0.7"/>')
        if kind == "star":
            pts = []
            for k in range(10):
                ang = -math.pi / 2 + k * math.pi / 5
                rad = 8.5 if k % 2 == 0 else 3.8
                pts.append(f"{x + rad*math.cos(ang):.1f},{y + rad*math.sin(ang):.1f}")
            out.append(f'<polygon points="{" ".join(pts)}" fill="{fill_c}" stroke="{stroke_c}" stroke-width="1.5"/>')
        elif kind == "diamond":
            out.append(f'<polygon points="{x:.1f},{y-7:.1f} {x+7:.1f},{y:.1f} {x:.1f},{y+7:.1f} {x-7:.1f},{y:.1f}" fill="{fill_c}" stroke="{stroke_c}" stroke-width="1.7"/>')
        elif kind == "cross":
            out.append(f'<line x1="{x-6:.1f}" y1="{y-6:.1f}" x2="{x+6:.1f}" y2="{y+6:.1f}" stroke="{stroke_c}" stroke-width="2.6"/>')
            out.append(f'<line x1="{x-6:.1f}" y1="{y+6:.1f}" x2="{x+6:.1f}" y2="{y-6:.1f}" stroke="{stroke_c}" stroke-width="2.6"/>')
        else:
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.8" fill="{fill_c}" stroke="{stroke_c}" stroke-width="1.6"/>')
        # Label pill
        bw = len(lbl) * 6.3 + 12
        out.append(f'<rect x="{tx - 4:.1f}" y="{ty - 13:.1f}" width="{bw:.1f}" height="17" rx="3.5" fill="#ffffff" fill-opacity="0.92" stroke="#cbd5e1" stroke-width="0.8"/>')
        out.append(f'<text x="{tx + 2:.1f}" y="{ty - 1:.1f}" font-size="10.8" font-weight="600" fill="{stroke_c}">{lbl}</text>')

    out.append('</g>')

    # Axes ticks & labels
    ticks = [0.6, 0.7, 0.8, 0.9, 1.0]
    for tv in ticks:
        tx, ty = sx(tv), sy(tv)
        out.append(f'<line x1="{tx:.1f}" y1="{mt+ph}" x2="{tx:.1f}" y2="{mt+ph+6}" stroke="#334155" stroke-width="1.2"/>')
        out.append(f'<text x="{tx:.1f}" y="{mt+ph+22}" font-size="12" fill="#334155" text-anchor="middle">{tv:.1f}</text>')
        out.append(f'<line x1="{ml-6}" y1="{ty:.1f}" x2="{ml}" y2="{ty:.1f}" stroke="#334155" stroke-width="1.2"/>')
        out.append(f'<text x="{ml-12}" y="{ty+4:.1f}" font-size="12" fill="#334155" text-anchor="end">{tv:.1f}</text>')

    out.append(f'<text x="{ml + pw/2:.1f}" y="{H - 22}" font-size="14" font-weight="600" fill="#0f172a" text-anchor="middle">Companion mass m₂  (with shuttle mass m₁ = 1)</text>')
    out.append(f'<text x="26" y="{mt + ph/2:.1f}" font-size="14" font-weight="600" fill="#0f172a" text-anchor="middle" transform="rotate(-90 26 {mt + ph/2:.1f})">Companion mass m₃  (with shuttle mass m₁ = 1)</text>')

    # Colorbar & legend on right
    cbx, cby, cbw, cbh = ml + pw + 28, mt + 20, 20, 290
    out.append(f'<text x="{cbx}" y="{mt + 8}" font-size="12" font-weight="700" fill="#0f172a">3D Stability Margin</text>')
    out.append(f'<rect x="{cbx}" y="{cby}" width="{cbw}" height="{cbh}" fill="url(#cbarGrad)" stroke="#94a3b8" stroke-width="1"/>')
    for val in [0.0, 0.4, 0.8, 1.2, 1.6]:
        yy = cby + cbh * (1.0 - val / max_marg)
        out.append(f'<line x1="{cbx+cbw}" y1="{yy:.1f}" x2="{cbx+cbw+5}" y2="{yy:.1f}" stroke="#334155" stroke-width="1"/>')
        out.append(f'<text x="{cbx+cbw+9}" y="{yy+4:.1f}" font-size="11.5" fill="#334155">{val:.1f}</text>')

    ly = cby + cbh + 42
    out.append(f'<text x="{cbx}" y="{ly}" font-size="12" font-weight="700" fill="#0f172a">Instability Regimes</text>')
    out.append(f'<rect x="{cbx}" y="{ly+12}" width="18" height="14" fill="url(#kreinHatch)" stroke="#fca5a5" stroke-width="1"/>')
    out.append(f'<text x="{cbx+25}" y="{ly+23}" font-size="11.5" fill="#334155">Krein Quartet (K)</text>')
    out.append(f'<rect x="{cbx}" y="{ly+34}" width="18" height="14" fill="#e2e8f0" stroke="#cbd5e1" stroke-width="1"/>')
    out.append(f'<text x="{cbx+25}" y="{ly+45}" font-size="11.5" fill="#334155">Hyperbolic (H)</text>')
    out.append(f'<text x="{cbx}" y="{ly+78}" font-size="10.8" fill="#475569">All 3D-stable points:</text>')
    out.append(f'<text x="{cbx}" y="{ly+94}" font-size="10.8" fill="#475569">|k₁|, |k₂|, |k_z| &lt; 2</text>')
    out.append(f'<text x="{cbx}" y="{ly+110}" font-size="10.8" fill="#475569">Mean (m₂+m₃)/2 &lt; m₁</text>')
    out.append(f'<text x="{cbx}" y="{ly+126}" font-size="10.3" fill="#64748b">(incl. [1, 1, 0.88])</text>')

    out.append('</svg>')
    with open(out_path, "w") as f:
        f.write("\n".join(out))
    print(f"Wrote {out_path}")


def rk4_integrate_2d(y0, m, T, steps=3200):
    def deriv(s):
        dx01, dy01 = s[2]-s[0], s[3]-s[1]
        dx02, dy02 = s[4]-s[0], s[5]-s[1]
        dx12, dy12 = s[4]-s[2], s[5]-s[3]
        ir01 = (dx01*dx01 + dy01*dy01) ** -1.5
        ir02 = (dx02*dx02 + dy02*dy02) ** -1.5
        ir12 = (dx12*dx12 + dy12*dy12) ** -1.5
        return [
            s[6], s[7], s[8], s[9], s[10], s[11],
            m[1]*dx01*ir01 + m[2]*dx02*ir02, m[1]*dy01*ir01 + m[2]*dy02*ir02,
            -m[0]*dx01*ir01 + m[2]*dx12*ir12, -m[0]*dy01*ir01 + m[2]*dy12*ir12,
            -m[0]*dx02*ir02 - m[1]*dx12*ir12, -m[0]*dy02*ir02 - m[1]*dy12*ir12,
        ]
    h = T / steps
    y = list(y0)
    traj = [list(y)]
    for _ in range(steps):
        k1 = deriv(y)
        k2 = deriv([y[i] + 0.5*h*k1[i] for i in range(12)])
        k3 = deriv([y[i] + 0.5*h*k2[i] for i in range(12)])
        k4 = deriv([y[i] + h*k3[i] for i in range(12)])
        y = [y[i] + (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(12)]
        traj.append(list(y))
    return traj


def build_corotating_svg(quad_json_path, out_path):
    data = json.load(open(quad_json_path))
    orb = next(o for o in data["orbits"] if o["id"] == "AAAbbb-DEM002")
    m = orb["masses"]
    T = float(orb["T"])
    pos = orb["pos_34d"]
    vel = orb["vel_34d"]
    y0 = [
        float(pos[0][0]), float(pos[0][1]),
        float(pos[1][0]), float(pos[1][1]),
        float(pos[2][0]), float(pos[2][1]),
        float(vel[0][0]), float(vel[0][1]),
        float(vel[1][0]), float(vel[1][1]),
        float(vel[2][0]), float(vel[2][1]),
    ]
    traj = rk4_integrate_2d(y0, m, T, steps=3600)

    # Compute co-rotating frame aligned with the m2-m3 binary axis (origin = midpoint of r2, r3; +x along r2 - r3)
    corot_b1 = []
    corot_b2 = []
    corot_b3 = []
    for s in traj:
        mx = 0.5 * (s[2] + s[4])
        my = 0.5 * (s[3] + s[5])
        ux = s[2] - s[4]
        uy = s[3] - s[5]
        r23 = math.hypot(ux, uy)
        ca, sa = ux / r23, uy / r23
        def rot_pt(px, py):
            dx, dy = px - mx, py - my
            return (ca*dx + sa*dy, -sa*dx + ca*dy)
        corot_b1.append(rot_pt(s[0], s[1]))
        corot_b2.append(rot_pt(s[2], s[3]))
        corot_b3.append(rot_pt(s[4], s[5]))

    W, H = 960, 490
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Inter, -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif">',
        f'<rect x="0" y="0" width="{W}" height="{H}" rx="12" fill="#ffffff" stroke="#e2e8f0" stroke-width="1.5"/>',
        f'<text x="36" y="34" font-size="18" font-weight="700" fill="#0f172a">AAAbbb Shuttle-Exchange (2D-DEM-002 · m = [1, 0.9, 0.8]): Inertial vs. Binary Co-Rotating Frame</text>',
        f'<text x="36" y="54" font-size="12.5" fill="#475569">T = 15.4951 · T* = 21.8013 · Λ = 1.8043 · Full-3D Margin = 1.20605 · Windings (w₁₂, w₁₃, w₂₃) = (−5, −5, −2)</text>',
    ]

    # Left panel: Inertial frame
    lx, ly, lw, lh = 36, 74, 428, 372
    out.append(f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.2"/>')
    out.append(f'<text x="{lx+14}" y="{ly+24}" font-size="13.5" font-weight="700" fill="#0f172a">A. Inertial Center-of-Mass Frame (1 Period)</text>')
    lim_l = 1.48
    def px_l(x, y):
        cx, cy = lx + lw / 2.0, ly + lh / 2.0 + 10
        sc = (min(lw, lh) - 52) / (2 * lim_l)
        return cx + x * sc, cy - y * sc

    # Subtle crosshair
    cx0, cy0 = px_l(0, 0)
    out.append(f'<line x1="{lx+18}" y1="{cy0:.1f}" x2="{lx+lw-18}" y2="{cy0:.1f}" stroke="#cbd5e1" stroke-width="0.9" stroke-dasharray="4,4"/>')
    out.append(f'<line x1="{cx0:.1f}" y1="{ly+34}" x2="{cx0:.1f}" y2="{ly+lh-16}" stroke="#cbd5e1" stroke-width="0.9" stroke-dasharray="4,4"/>')

    colors = ["#dc2626", "#2563eb", "#059669"]
    labels = ["m₁ = 1.0 (Shuttle)", "m₂ = 0.9", "m₃ = 0.8"]
    for b in [2, 1, 0]:
        pts = " ".join(f"{px_l(s[2*b], s[2*b+1])[0]:.1f},{px_l(s[2*b], s[2*b+1])[1]:.1f}" for s in traj[::4])
        sw = 2.0 if b == 0 else 1.5
        op = 0.92 if b == 0 else 0.78
        out.append(f'<polyline points="{pts}" fill="none" stroke="{colors[b]}" stroke-width="{sw}" opacity="{op}"/>')
    for b in range(3):
        x, y = px_l(y0[2*b], y0[2*b+1])
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{6.2 if b==0 else 5.2}" fill="{colors[b]}" stroke="#ffffff" stroke-width="1.6"/>')

    # Right panel: Co-rotating frame (m2-m3 binary axis)
    rx, ry, rw, rh = 496, 74, 428, 372
    out.append(f'<rect x="{rx}" y="{ry}" width="{rw}" height="{rh}" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.2"/>')
    out.append(f'<text x="{rx+14}" y="{ry+24}" font-size="13.5" font-weight="700" fill="#0f172a">B. Co-Rotating m₂–m₃ Binary Frame: 3 + 3 Exchange (AAAbbb)</text>')
    lim_r = 1.90
    def px_r(x, y):
        cx, cy = rx + rw / 2.0, ry + rh / 2.0 + 10
        sc = (min(rw, rh) - 44) / (2 * lim_r)
        return cx + x * sc, cy - y * sc

    cxr, cyr = px_r(0, 0)
    out.append(f'<line x1="{rx+18}" y1="{cyr:.1f}" x2="{rx+rw-18}" y2="{cyr:.1f}" stroke="#cbd5e1" stroke-width="0.9" stroke-dasharray="4,4"/>')
    out.append(f'<line x1="{cxr:.1f}" y1="{ry+34}" x2="{cxr:.1f}" y2="{ry+rh-16}" stroke="#cbd5e1" stroke-width="0.9" stroke-dasharray="4,4"/>')

    # Draw m2 and m3 radial oscillation bars along x-axis
    x2_lo = min(p[0] for p in corot_b2); x2_hi = max(p[0] for p in corot_b2)
    x3_lo = min(p[0] for p in corot_b3); x3_hi = max(p[0] for p in corot_b3)
    p2a, p2b = px_r(x2_lo, 0), px_r(x2_hi, 0)
    p3a, p3b = px_r(x3_lo, 0), px_r(x3_hi, 0)
    out.append(f'<line x1="{p2a[0]:.1f}" y1="{p2a[1]:.1f}" x2="{p2b[0]:.1f}" y2="{p2b[1]:.1f}" stroke="#2563eb" stroke-width="5" stroke-linecap="round" opacity="0.45"/>')
    out.append(f'<line x1="{p3a[0]:.1f}" y1="{p3a[1]:.1f}" x2="{p3b[0]:.1f}" y2="{p3b[1]:.1f}" stroke="#059669" stroke-width="5" stroke-linecap="round" opacity="0.45"/>')

    # Color first half [0, T/2] (3 loops around m2/m3) vs second half [T/2, T] so the exchange is crystal-clear!
    half = len(corot_b1) // 2
    pts1 = " ".join(f"{px_r(p[0], p[1])[0]:.1f},{px_r(p[0], p[1])[1]:.1f}" for p in corot_b1[:half+1:3])
    pts2 = " ".join(f"{px_r(p[0], p[1])[0]:.1f},{px_r(p[0], p[1])[1]:.1f}" for p in corot_b1[half::3])
    out.append(f'<polyline points="{pts1}" fill="none" stroke="#dc2626" stroke-width="2.1" opacity="0.92"/>')
    out.append(f'<polyline points="{pts2}" fill="none" stroke="#f59e0b" stroke-width="2.1" opacity="0.90" stroke-dasharray="4.5,2.5"/>')

    # Mean positions of m2 and m3 + starting marker of m1
    m2x, m2y = px_r(0.5*(x2_lo+x2_hi), 0)
    m3x, m3y = px_r(0.5*(x3_lo+x3_hi), 0)
    out.append(f'<circle cx="{m2x:.1f}" cy="{m2y:.1f}" r="6.5" fill="#2563eb" stroke="#ffffff" stroke-width="1.6"/>')
    out.append(f'<circle cx="{m3x:.1f}" cy="{m3y:.1f}" r="6.5" fill="#059669" stroke="#ffffff" stroke-width="1.6"/>')
    out.append(f'<text x="{m2x:.1f}" y="{m2y+104:.1f}" font-size="12" font-weight="700" fill="#1d4ed8" text-anchor="middle">m₂ (3 prograde loops)</text>')
    out.append(f'<text x="{m3x:.1f}" y="{m3y+104:.1f}" font-size="12" font-weight="700" fill="#047857" text-anchor="middle">m₃ (3 prograde loops)</text>')
    s1x, s1y = px_r(corot_b1[0][0], corot_b1[0][1])
    out.append(f'<circle cx="{s1x:.1f}" cy="{s1y:.1f}" r="5.8" fill="#dc2626" stroke="#ffffff" stroke-width="1.6"/>')
    out.append(f'<text x="{s1x+12:.1f}" y="{s1y-11:.1f}" font-size="11.5" font-weight="700" fill="#991b1b" text-anchor="middle">t = 0</text>')

    # Bottom legend
    out.append(f'<circle cx="52" cy="468" r="5" fill="#dc2626"/><text x="64" y="472" font-size="12" fill="#334155">m₁ = 1.0 (heavy shuttle · solid red = [0, T/2], dashed amber = [T/2, T] in Panel B)</text>')
    out.append(f'<circle cx="585" cy="468" r="5" fill="#2563eb"/><text x="597" y="472" font-size="12" fill="#334155">m₂ = 0.9</text>')
    out.append(f'<circle cx="685" cy="468" r="5" fill="#059669"/><text x="697" y="472" font-size="12" fill="#334155">m₃ = 0.8</text>')

    out.append('</svg>')
    with open(out_path, "w") as f:
        f.write("\n".join(out))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    grid_csv = "data/aaabbb/mass_grid_002.csv" if os.path.exists("data/aaabbb/mass_grid_002.csv") else "data/aaabbb/mass_grid_005.csv"
    build_stability_island_svg(grid_csv, "assets/aaabbb_stability_island.svg")
    if os.path.exists("data/aaabbb/aaabbb_quad32.json"):
        build_corotating_svg("data/aaabbb/aaabbb_quad32.json", "assets/aaabbb_corotating_vs_inertial.svg")
