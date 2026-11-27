#!/usr/bin/env python3
import glob
import json
import math
import os
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def load_jsonl(pattern):
    rows = []
    for path in sorted(glob.glob(os.path.join(DATA, pattern))):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except Exception:
                    pass
    return rows


def scale_invariants(o):
    E = o["E"]
    T = o["T"]
    L = abs(o["L"])
    inv1 = E * (T ** (2.0 / 3.0))
    inv2 = L * math.sqrt(max(0.0, -E))
    return inv1, inv2


def deduplicate(rows):
    # Sort by return error ascending so we keep the cleanest representative
    rows = sorted(rows, key=lambda o: o.get("ret_err", 1.0))
    kept = []
    for o in rows:
        # Strict quality filter
        if o.get("survived_100T", 0) < 100:
            continue
        if o.get("stab_margin", -1) <= 1e-4:
            continue
        if o.get("ret_err", 1.0) > 1e-7:
            continue
        if o["dim"] == 3 and o.get("thickness_3d", 0.0) < 0.05:
            continue

        inv1, inv2 = scale_invariants(o)
        m_key = tuple(round(x, 3) for x in o["masses"])
        is_dup = False
        for k in kept:
            if k["dim"] != o["dim"]:
                continue
            km = tuple(round(x, 3) for x in k["masses"])
            if km != m_key:
                continue
            kinv1, kinv2 = scale_invariants(k)
            if abs(inv1 - kinv1) < 0.03 * max(1.0, abs(inv1)) and abs(inv2 - kinv2) < 0.03 * max(0.2, abs(inv2)):
                is_dup = True
                break
        if not is_dup:
            kept.append(o)
    return kept


def load_continuation(prefix):
    pts = []
    for suffix in ("dn", "up"):
        p = os.path.join(DATA, "continuations", f"{prefix}_{suffix}.dat")
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 8:
                    continue
                pts.append({
                    "m3": float(parts[0]),
                    "vx": float(parts[1]),
                    "vy": float(parts[2]),
                    "T": float(parts[3]),
                    "k1": float(parts[4]),
                    "k2": float(parts[5]),
                    "stable": int(parts[6]) == 1,
                    "margin": float(parts[7]),
                })
    # Deduplicate by m3 and sort
    by_m = {}
    for pt in pts:
        by_m[round(pt["m3"], 6)] = pt
    return [by_m[k] for k in sorted(by_m.keys())]


def analyze_continuation(pts):
    if not pts:
        return None
    st = [p for p in pts if p["stable"]]
    if not st:
        return {"n_total": len(pts), "n_stable": 0}
    best = max(st, key=lambda p: p["margin"])
    return {
        "n_total": len(pts),
        "n_stable": len(st),
        "m3_min_traced": pts[0]["m3"],
        "m3_max_traced": pts[-1]["m3"],
        "m3_stable_min": st[0]["m3"],
        "m3_stable_max": st[-1]["m3"],
        "best": best,
    }


def main():
    raw_2d = load_jsonl("discoveries_2d*.jsonl")
    raw_3d = load_jsonl("discoveries_3d*.jsonl")
    all_raw = raw_2d + raw_3d
    unique = deduplicate(all_raw)

    # Sort: 3D first, then 2D democratic by stability margin, then 2D hierarchical by stability margin
    unique.sort(key=lambda o: (
        0 if o["dim"] == 3 else (1 if o["democratic"] else 2),
        -o["stab_margin"],
    ))

    # Assign clean IDs
    c3, c2d, c2h = 0, 0, 0
    for o in unique:
        if o["dim"] == 3:
            c3 += 1
            o["id"] = f"3D-{c3:03d}"
        elif o["democratic"]:
            c2d += 1
            o["id"] = f"2D-DEM-{c2d:03d}"
        else:
            c2h += 1
            o["id"] = f"2D-HIE-{c2h:03d}"

    # Load 3D vertical stability kz for 2D orbits if present
    kz_path = os.path.join(DATA, "kz_results.dat")
    kz_map = {}
    if os.path.exists(kz_path):
        with open(kz_path) as f:
            for line in f:
                p = line.split()
                if len(p) >= 3:
                    kz_map[p[0]] = (float(p[1]), int(p[2]) == 1)
    for o in unique:
        if o["dim"] == 2 and o["id"] in kz_map:
            kzv, s3 = kz_map[o["id"]]
            o["kz"] = kzv
            o["stable_3d"] = s3
            o["stab_margin_3d"] = (2.0 - max(abs(o["k1"]), abs(o["k2"]), abs(kzv))) if s3 else 0.0
        elif o["dim"] == 3:
            o["stable_3d"] = True
            o["stab_margin_3d"] = o["stab_margin"]

    # Load continuations
    conts = {
        "Braid-1.10 (m3 around 1.10)": analyze_continuation(load_continuation("braid110")),
        "Heavy-Center-3.00 (m3 around 3.00)": analyze_continuation(load_continuation("heavy300")),
        "Heavy-Braid-1.75 (m3 around 1.75)": analyze_continuation(load_continuation("dem175")),
        "Light-Braid-0.55 (m3 around 0.55)": analyze_continuation(load_continuation("dem055")),
        "Figure-8 (m3 around 1.00)": analyze_continuation(load_continuation("fig8")),
    }

    catalog = {
        "generated_at": datetime.now().isoformat(),
        "total_raw": len(all_raw),
        "total_unique": len(unique),
        "counts": {
            "3d": sum(1 for o in unique if o["dim"] == 3),
            "2d_democratic": sum(1 for o in unique if o["dim"] == 2 and o["democratic"]),
            "2d_hierarchical": sum(1 for o in unique if o["dim"] == 2 and not o["democratic"]),
            "2d_collinear_L0": sum(1 for o in unique if o["subspace"] == "collinear_L0"),
            "2d_perp_L_nonzero": sum(1 for o in unique if o["subspace"] == "perp_L_nonzero"),
        },
        "continuations": conts,
        "orbits": unique,
    }

    with open(os.path.join(ROOT, "catalog.json"), "w") as f:
        json.dump(catalog, f, indent=2)

    # Write CATALOG.md
    md = []
    md.append("# 3-Body Steady (Linearly Stable) Periodic Orbit Discovery Catalog\n")
    md.append(f"Generated: `{catalog['generated_at']}`  \n")
    md.append(f"**Total unique verified linearly stable orbits:** **{catalog['total_unique']}** (from {catalog['total_raw']} raw detections)\n")
    md.append(f"- **3D Non-Planar Linearly Stable Orbits:** {catalog['counts']['3d']}\n")
    md.append(f"- **2D Democratic (Partner-Swapping) Linearly Stable Orbits:** {catalog['counts']['2d_democratic']}\n")
    md.append(f"- **2D Hierarchical Linearly Stable Orbits:** {catalog['counts']['2d_hierarchical']}\n")
    md.append(f"- **By Subspace:** {catalog['counts']['2d_collinear_L0']} collinear ($L=0$), {catalog['counts']['2d_perp_L_nonzero']} perpendicular-crossing ($L\\neq 0$), {catalog['counts']['3d']} 3D non-planar\n")

    md.append("\n## Verification Criteria\n")
    md.append("Every orbit in this catalog passed **all four** rigorous criteria:\n")
    md.append("1. **Periodicity to high precision:** $|\\mathbf{z}(T) - \\mathbf{z}(0)| < 10^{-7}$ (most $< 10^{-10}$) using adaptive Dormand-Prince 5(4) at tolerance $10^{-13}$.\n")
    md.append("2. **Exact Floquet linear stability:** All non-trivial Broucke stability indices of the $12\\times 12$ (2D) or $18\\times 18$ (3D) monodromy matrix are real with $|k_j| < 2$ (stability margin $2 - \\max_j |k_j| > 0$, Lyapunov exponent $= 0$).\n")
    md.append("3. **100-period empirical perturbation survival:** Survives 100/100 full periods bounded after a momentum-preserving $10^{-5}$ velocity kick (in all 3 axes for 3D orbits).\n")
    md.append("4. **True 3D check (for 3D orbits):** Out-of-best-fit-plane RMS thickness $h_{\\text{3D}} > 0.05$ (strictly excluding tilted planar orbits).\n")

    # 3D section
    orbits_3d = [o for o in unique if o["dim"] == 3]
    if orbits_3d:
        md.append("\n## 3D Non-Planar Steady Orbits\n")
        md.append("| ID | Subspace | Masses | Period $T$ | $|\\mathbf{L}|$ | $r_{\\min}$ | 3D Thickness | Swaps | Stability $(k_1, k_2, k_3)$ | Margin | 100T Drift |\n")
        md.append("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for o in orbits_3d:
            md.append(
                f"| `{o['id']}` | `{o['subspace']}` | `{o['masses']}` | {o['T']:.4f} | {o['L']:.3f} | "
                f"{o['rmin']:.4f} | {o['thickness_3d']:.4f} | {o['partner_swaps']} | "
                f"({o['k1']:+.3f}, {o['k2']:+.3f}, {o['k3']:+.3f}) | **{o['stab_margin']:.4f}** | {o['drift_100T']:.5f} |\n"
            )

    # 2D Democratic section
    orbits_2d_dem = [o for o in unique if o["dim"] == 2 and o["democratic"]]
    md.append("\n## 2D Democratic (Partner-Swapping) Steady Orbits\n")
    md.append("Ranked by linear stability margin $2 - \\max(|k_1|, |k_2|)$ (higher = deeper inside the KAM stability island; compare Figure-8 margin $= 0.0028$):\n\n")
    md.append("| ID | Subspace | Masses | Period $T$ | $L$ | $r_{\\min}$ | Swaps | Stability $(k_1, k_2)$ | Margin | 3D Stable ($k_z$) | 100T Drift |\n")
    md.append("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
    for o in orbits_2d_dem:
        kz_str = f"YES ({o['kz']:+.2f})" if o.get("stable_3d") else ("L=0 (2.00)" if "kz" in o else "—")
        md.append(
            f"| `{o['id']}` | `{o['subspace']}` | `{o['masses']}` | {o['T']:.4f} | {o['L']:+.3f} | "
            f"{o['rmin']:.4f} | {o['partner_swaps']} | "
            f"({o['k1']:+.4f}, {o['k2']:+.4f}) | **{o['stab_margin']:.4f}** | {kz_str} | {o['drift_100T']:.5f} |\n"
        )

    # 2D Hierarchical section (top 20)
    orbits_2d_hie = [o for o in unique if o["dim"] == 2 and not o["democratic"]]
    if orbits_2d_hie:
        md.append(f"\n## 2D Hierarchical Steady Orbits (Top {min(20, len(orbits_2d_hie))} of {len(orbits_2d_hie)})\n")
        md.append("| ID | Masses | Period $T$ | $L$ | $r_{\\min}$ | $r_{\\max}$ | Stability $(k_1, k_2)$ | Margin | 100T Drift |\n")
        md.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for o in orbits_2d_hie[:20]:
            md.append(
                f"| `{o['id']}` | `{o['masses']}` | {o['T']:.4f} | {o['L']:+.3f} | "
                f"{o['rmin']:.4f} | {o['rmax']:.3f} | "
                f"({o['k1']:+.4f}, {o['k2']:+.4f}) | **{o['stab_margin']:.4f}** | {o['drift_100T']:.5f} |\n"
            )

    # Mass continuation section
    md.append("\n## Mass-Continuation Stability Windows ($m_1 = m_2 = 1$, varying $m_3$)\n")
    for name, c in conts.items():
        if not c or c.get("n_stable", 0) == 0:
            continue
        b = c["best"]
        md.append(f"### {name}\n")
        md.append(f"- **Traced range:** $m_3 \\in [{c['m3_min_traced']:.4f},\\; {c['m3_max_traced']:.4f}]$ ({c['n_total']} points)\n")
        md.append(f"- **Linearly stable window:** $m_3 \\in [{c['m3_stable_min']:.4f},\\; {c['m3_stable_max']:.4f}]$ ({c['n_stable']} stable points)\n")
        md.append(f"- **Peak stability:** at $m_3^* = {b['m3']:.4f}$ with margin **{b['margin']:.4f}** ($k_1 = {b['k1']:+.4f},\\; k_2 = {b['k2']:+.4f}$, $T = {b['T']:.6f}$)\n\n")

    # Exact initial conditions for top showcase orbits
    md.append("\n## Exact Initial Conditions (Top Showcase Orbits)\n")
    showcase = orbits_3d[:8] + orbits_2d_dem[:12] + orbits_2d_hie[:5]
    for o in showcase:
        md.append(f"### `{o['id']}` — `{o['subspace']}`, masses `{o['masses']}`, $T = {o['T']:.12f}$\n")
        md.append(f"- **Margin:** `{o['stab_margin']:.6f}`, **Return err:** `{o['ret_err']:.2e}`, **100T drift:** `{o['drift_100T']:.6f}`\n")
        md.append(f"- **`y0`:** `{json.dumps(o['y0'])}`\n\n")

    with open(os.path.join(ROOT, "CATALOG.md"), "w") as f:
        f.writelines(md)

    print(f"Wrote CATALOG.md and catalog.json: {catalog['total_unique']} unique steady orbits "
          f"({catalog['counts']['3d']} 3D, {catalog['counts']['2d_democratic']} 2D democratic, "
          f"{catalog['counts']['2d_hierarchical']} 2D hierarchical).")


if __name__ == "__main__":
    main()
