#!/usr/bin/env python3
"""Render smooth animated GIFs for the showcase 3D and 2D new orbit families."""
import json
import math
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(ROOT, "assets")
BIN = os.path.join(ROOT, "bin", "render_gif")

# Showcase list: (orbit_id, out_filename, frames, fps, T_ghost_override, T_anim_override, yaw0, yaw_span, pitch)
SHOWCASE = [
    # 3D-Δ14 Family (5-fold exact inertial closure: show full 5-fold rose ghost + 1 full 5-fold loop with 360° gentle 3D camera rotation!)
    ("3D-001", "3d_delta14_5fold.gif", 150, 25, None, None, 0.4, 2.0 * math.pi, 0.42),
    # 3D-Δ14 Family (7-fold exact inertial closure)
    ("3D-004", "3d_delta14_7fold.gif", 160, 25, None, None, 0.4, 2.0 * math.pi, 0.45),
    # 3D Distinct-Mass Trio (1, 0.7, 0.4) — RPO, show 3 periods with 360° 3D camera rotation
    ("3D-005", "3d_rpo_distinct.gif", 120, 25, 33.0985, 33.0985, 0.5, 2.0 * math.pi, 0.48),
    # 2D A^7 Deepest Stability Island (margin 1.8725, 668x deeper than Fig-8)
    ("2D-DEM-001", "2d_A7_deepest.gif", 110, 25, None, None, 0.0, 0.0, 0.0),
    # 2D AABaaB All-Distinct-Mass (1, 0.9, 0.8), full-3D stable
    ("2D-DEM-006", "2d_AABaaB_distinct.gif", 100, 25, None, None, 0.0, 0.0, 0.0),
    # 2D (Ab)^7 Wide-Separation (rmin=1.272), 14-letter alternating dance
    ("2D-DEM-013", "2d_Ab7_wide.gif", 120, 25, None, None, 0.0, 0.0, 0.0),
    # 2D AABABABabaB Heavy-Center L=0 (m=[1,1,3], stable m3 in [2.84, 3.42])
    ("2D-DEM-022", "2d_heavy_center.gif", 100, 25, None, None, 0.0, 0.0, 0.0),
    # 2D 20-letter Unequal-Mass L=0 Braid (m=[1,1,0.55], rmin=0.161)
    ("2D-DEM-015", "2d_braid20.gif", 120, 25, None, None, 0.0, 0.0, 0.0),
]


def main():
    os.makedirs(ASSETS, exist_ok=True)
    os.makedirs(os.path.dirname(BIN), exist_ok=True)

    # Compile renderer
    subprocess.run(
        ["g++", "-O3", "-march=native", "-std=c++20",
         os.path.join(ROOT, "src", "render_gif.cpp"), "-o", BIN],
        check=True,
    )

    cat = json.load(open(os.path.join(ROOT, "catalog.json")))
    by_id = {o["id"]: o for o in cat["orbits"]}

    W, H = 480, 360
    for oid, fname, frames, fps, tg_ov, ta_ov, yaw0, yaw_span, pitch in SHOWCASE:
        o = by_id[oid]
        T_ghost = tg_ov if tg_ov is not None else o["T"]
        T_anim = ta_ov if ta_ov is not None else o["T"]
        payload = " ".join(
            [str(o["dim"]), str(W), str(H), str(frames),
             f"{T_ghost:.12f}", f"{T_anim:.12f}",
             f"{yaw0:.6f}", f"{yaw_span:.6f}", f"{pitch:.6f}"]
            + [f"{x:.12f}" for x in o["masses"]]
            + [f"{x:.13f}" for x in o["y0"]]
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
        p_render.wait()
        sz = os.path.getsize(out_gif) // 1024
        print(f"Rendered {fname:26s} ({oid:10s}, dim={o['dim']}D, {frames}f @ {fps}fps) -> {sz:4d} KB")


if __name__ == "__main__":
    main()
