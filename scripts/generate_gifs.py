#!/usr/bin/env python3
"""Render smooth animated GIFs for the showcase 3D, 2D, 1PN relativistic, and 4-body steady orbit families."""
import json
import math
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(ROOT, "assets")
BIN = os.path.join(ROOT, "bin", "render_gif")

# Showcase list: (orbit_id, out_filename, frames, fps, T_ghost_override, T_anim_override, yaw0, yaw_span, pitch)
SHOWCASE = [
    # --- Track 1: Linearly Stable Free-Fall (Brake) Orbit (cleanest rmin=0.00445 representative) ---
    ("2D-DEM-015", "2d_brake_stable.gif", 100, 25, None, None, 0.0, 0.0, 0.0),
    # --- Track 2: Period-Doubled 3D BHH Vertical Bifurcation Family (3D-PD40/48) ---
    ("3D-012", "3d_vertical_pd_3fold.gif", 150, 25, None, None, 0.4, 2.0 * math.pi, 0.45),
    # --- Track 3: 2D & 3D Rotating-Frame Choreographies (CFM Gamma_z & Marchal P12 / Gamma_y) ---
    ("2D-DEM-001", "2d_choreo_5petal_5fold.gif", 130, 25, None, None, 0.0, 0.0, 0.0),
    ("2D-DEM-018", "2d_choreo_trefoil_3fold.gif", 120, 25, None, None, 0.0, 0.0, 0.0),
    ("3D-005", "3d_choreo_helical_11fold.gif", 160, 25, None, 62.1634, 0.3, 2.0 * math.pi, 0.42),
    # --- Track 4: Post-Newtonian (1PN ADM) Relativistic Closures (synchronized 80-bit T_inert IDs) ---
    ("1PN-021", "1pn_2d_dem001_10fold.gif", 150, 25, None, 128.3801, 0.0, 0.0, 0.0),
    ("1PN-001", "1pn_3d_delta14_12fold.gif", 150, 25, None, 111.1074, 0.4, 2.0 * math.pi, 0.44),
    # --- Track 5: 4-Body Full-3D Linearly Stable Planar 3+1 Hybrid (24x24 80-bit verified) ---
    ("4B-001", "4b_2d_hybrid_dem001.gif", 130, 25, None, None, 0.0, 0.0, 0.0),
    # --- Core 3D & 2D Showcase Families (synchronized 80-bit audited IDs) ---
    ("3D-002", "3d_delta14_5fold.gif", 150, 25, None, None, 0.4, 2.0 * math.pi, 0.42),
    ("3D-011", "3d_delta14_7fold.gif", 160, 25, None, None, 0.4, 2.0 * math.pi, 0.45),
    ("3D-016", "3d_rpo_distinct.gif", 120, 25, 33.0985, 33.0985, 0.5, 2.0 * math.pi, 0.48),
    ("2D-HIE-001", "2d_A7_deepest.gif", 110, 25, None, None, 0.0, 0.0, 0.0),
    ("2D-DEM-002", "2d_AAAbbb_shuttle.gif", 100, 25, None, None, 0.0, 0.0, 0.0),
    ("2D-DEM-016", "2d_heavy_center.gif", 100, 25, None, None, 0.0, 0.0, 0.0),
    ("2D-DEM-008", "2d_braid16_m055.gif", 120, 25, None, None, 0.0, 0.0, 0.0),
]


def main():
    os.makedirs(ASSETS, exist_ok=True)
    os.makedirs(os.path.dirname(BIN), exist_ok=True)

    subprocess.run(
        ["g++", "-O3", "-march=native", "-std=c++20",
         os.path.join(ROOT, "src", "render_gif.cpp"), "-o", BIN],
        check=True,
    )

    cat = json.load(open(os.path.join(ROOT, "catalog.json")))
    by_id = {}
    for bucket in ("orbits", "relativistic_1pn", "four_body"):
        for o in cat.get(bucket, []):
            by_id[o["id"]] = o

    W, H = 480, 360
    for oid, fname, frames, fps, tg_ov, ta_ov, yaw0, yaw_span, pitch in SHOWCASE:
        o = by_id[oid]
        nb = len(o["masses"])
        beta = float(o.get("beta_inv_c2", 0.0))
        state = o["z0_canonical"] if beta > 0.0 else o["y0"]
        T_ghost = tg_ov if tg_ov is not None else o["T"]
        T_anim = ta_ov if ta_ov is not None else o["T"]
        payload = " ".join(
            [str(nb), str(o["dim"]), f"{beta:.12f}", str(W), str(H), str(frames),
             f"{T_ghost:.12f}", f"{T_anim:.12f}",
             f"{yaw0:.6f}", f"{yaw_span:.6f}", f"{pitch:.6f}"]
            + [f"{x:.12f}" for x in o["masses"]]
            + [f"{x:.13f}" for x in state]
        ) + "\n"

        out_gif = os.path.join(ASSETS, fname)
        p_render = subprocess.Popen(
            [BIN],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_ff = subprocess.Popen(
            [
                "ffmpeg", "-y", "-loglevel", "error",
                "-f", "rawvideo", "-pix_fmt", "rgb24",
                "-s", f"{W}x{H}", "-r", str(fps), "-i", "pipe:0",
                "-vf", "split[s0][s1];[s0]palettegen=max_colors=128:stats_mode=diff[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3",
                out_gif,
            ],
            stdin=p_render.stdout,
        )
        p_render.stdout.close()
        p_render.stdin.write(payload.encode())
        p_render.stdin.close()
        p_ff.wait()
        rc = p_render.wait()
        if rc != 0:
            raise RuntimeError(f"render_gif failed for {oid} (rc={rc})")
        sz = os.path.getsize(out_gif) // 1024
        print(f"Rendered {fname:28s} ({oid:10s}, {nb}B {o['dim']}D beta={beta:.4f}, {frames}f) -> {sz:4d} KB")


if __name__ == "__main__":
    main()
