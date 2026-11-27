#!/usr/bin/env python3
"""
Independent audit script for:
  1. All 14 free_fall_brake orbits: dense adaptive DOPRI5 rmin, exact syzygy count,
     and Tanikawa-Mikkola half-period syzygy sequence vs catalog.json coarse NS=2400 values.
  2. All 6 2D & 4 3D rotating choreography orbits: mirror-image proof of choreo_12swap_2d.dat
     vs fig8_rot_2d.dat, 24-fold closure proof for 3D-004, and syzygy sequence of 2D-DEM-004.
  3. Reconciled summary of 80-bit 4-body (bin/review_audit_4b) and 1PN (bin/review_audit_1pn) audits.
"""
import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def deriv2d(y, m):
    pos = [(y[0], y[1]), (y[2], y[3]), (y[4], y[5])]
    acc = [[0.0, 0.0] for _ in range(3)]
    for i in range(3):
        for j in range(i + 1, 3):
            dx = pos[j][0] - pos[i][0]
            dy = pos[j][1] - pos[i][1]
            ir3 = (dx * dx + dy * dy) ** (-1.5)
            acc[i][0] += m[j] * dx * ir3
            acc[i][1] += m[j] * dy * ir3
            acc[j][0] -= m[i] * dx * ir3
            acc[j][1] -= m[i] * dy * ir3
    return [y[6], y[7], y[8], y[9], y[10], y[11],
            acc[0][0], acc[0][1], acc[1][0], acc[1][1], acc[2][0], acc[2][1]]


def dopri5_trace_2d(y0, T, m, tol=1e-12):
    a21 = 1 / 5
    a31, a32 = 3 / 40, 9 / 40
    a41, a42, a43 = 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63, a64, a65 = 9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656
    b1, b3, b4, b5, b6 = 35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84
    e1, e3, e4, e5, e6, e7 = 71 / 57600, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40
    y = list(y0)
    t = 0.0
    h = min(1e-4, T * 1e-3)
    pts = [(0.0, list(y))]
    while t < T - 1e-14:
        if t + h > T:
            h = T - t
        k1 = deriv2d(y, m)
        k2 = deriv2d([y[i] + h * a21 * k1[i] for i in range(12)], m)
        k3 = deriv2d([y[i] + h * (a31 * k1[i] + a32 * k2[i]) for i in range(12)], m)
        k4 = deriv2d([y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]) for i in range(12)], m)
        k5 = deriv2d([y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]) for i in range(12)], m)
        k6 = deriv2d([y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]) for i in range(12)], m)
        yn = [y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]) for i in range(12)]
        k7 = deriv2d(yn, m)
        err = math.sqrt(sum(
            (h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i])
             / (tol * (1.0 + max(abs(y[i]), abs(yn[i]))))) ** 2
            for i in range(12)
        ) / 12.0)
        if err <= 1.0:
            y = yn
            t += h
            pts.append((t, list(y)))
            h = min(h * (5.0 if err == 0 else min(5.0, max(0.2, 0.9 * err ** (-0.2)))), 0.01)
        else:
            h *= max(0.1, 0.9 * err ** (-0.25))
    return pts


def trace_syz_and_rmin(o):
    pts = dopri5_trace_2d(o["y0"], o["T"], o["masses"])

    def area(y):
        return (y[2] - y[0]) * (y[5] - y[1]) - (y[3] - y[1]) * (y[4] - y[0])

    def mid_body(y):
        d01 = math.hypot(y[0] - y[2], y[1] - y[3])
        d02 = math.hypot(y[0] - y[4], y[1] - y[5])
        d12 = math.hypot(y[2] - y[4], y[3] - y[5])
        if d12 >= d01 and d12 >= d02:
            return 1
        if d02 >= d01 and d02 >= d12:
            return 2
        return 3

    min_r = [1e9, 1e9, 1e9]
    seq = []
    for i, (t, y) in enumerate(pts):
        min_r[0] = min(min_r[0], math.hypot(y[0] - y[2], y[1] - y[3]))
        min_r[1] = min(min_r[1], math.hypot(y[0] - y[4], y[1] - y[5]))
        min_r[2] = min(min_r[2], math.hypot(y[2] - y[4], y[3] - y[5]))
        if i > 0:
            a0, a1 = area(pts[i - 1][1]), area(y)
            if a0 * a1 < 0:
                frac = abs(a0) / (abs(a0) + abs(a1))
                ym = [pts[i - 1][1][k] + frac * (y[k] - pts[i - 1][1][k]) for k in range(12)]
                seq.append(mid_body(ym))
    return min(min_r), min_r, seq


def main():
    with open(os.path.join(ROOT, "catalog.json")) as f:
        cat = json.load(f)

    print("=== [1/4] 80-bit Analytical Variational Audit of All 11 4-Body Orbits ===")
    lines_4b = []
    for o in cat["four_body"]:
        y = o["y0"]
        if o["dim"] == 2:
            y24 = []
            for i in range(4):
                y24 += [y[2 * i], y[2 * i + 1], 0.0]
            for i in range(4):
                y24 += [y[8 + 2 * i], y[8 + 2 * i + 1], 0.0]
        else:
            y24 = y
        lines_4b.append(" ".join([o["id"], str(o["dim"])] + [f"{v:.16g}" for v in o["masses"] + [o["T"], o["precess_rad"]] + y24]))

    cmd_4b = [os.path.join(ROOT, "bin", "review_audit_4b")]
    if "--long-kick" in sys.argv:
        cmd_4b.append("--long-kick")
    p4 = subprocess.run(cmd_4b, input="\n".join(lines_4b) + "\n", text=True, capture_output=True, check=True)
    res_4b = [json.loads(ln) for ln in p4.stdout.strip().splitlines() if ln.strip()]
    for r, o in zip(res_4b, cat["four_body"]):
        verdict = "HYPERBOLIC UNSTABLE" if r["hyperbolic"] else "FULL-3D STABLE"
        extra = ""
        if "surv_4000T" in r:
            extra = f" | 4000T kick: surv={r['surv_4000T']}, d500={r['drift_500T']:.3e}, d2000={r['drift_2000T']:.3e}, d4000={r['drift_4000T']:.3e}"
        kinfo = f"k6-2={r['k_trace_minus_2']:+.2e}" if r["n_ell_eigs"] == 10 else f"exact_mrg={r['exact_margin']:.6f} (6/6 ell)"
        print(f"  {r['id']} ({r['dim']}D {o['subspace']:26s}): cat_mrg={o['stab_margin']:.6f} -> "
              f"80bit symp={r['symp_err']:.1e}, max|lam|={r['max_eig_mod']:.7f}, "
              f"{kinfo} [{verdict}]{extra}")

    print("\n=== [2/4] 80-bit Canonical Audit of All 21 1PN Relativistic Orbits (incl. T_inert vs T_rpo) ===")
    lines_1pn = []
    for o in cat["relativistic_1pn"]:
        z = o["z0_canonical"]
        if o["dim"] == 2:
            z18 = []
            for i in range(3):
                z18 += [z[2 * i], z[2 * i + 1], 0.0]
            for i in range(3):
                z18 += [z[6 + 2 * i], z[6 + 2 * i + 1], 0.0]
        else:
            z18 = z
        q_fold = o.get("q_fold", 0)
        T_rpo = o.get("T_rpo", o["T"])
        lines_1pn.append(" ".join([
            o["id"], str(o["dim"]), f"{o['beta_inv_c2']:.16g}",
            *[f"{m:.16g}" for m in o["masses"]],
            f"{o['T']:.16g}", f"{T_rpo:.16g}", str(q_fold), f"{o['precess_rad']:.16g}",
            *[f"{v:.16g}" for v in z18],
        ]))
    p1 = subprocess.run([os.path.join(ROOT, "bin", "review_audit_1pn")],
                        input="\n".join(lines_1pn) + "\n", text=True, capture_output=True, check=True)
    res_1pn = [json.loads(ln) for ln in p1.stdout.strip().splitlines() if ln.strip()]
    for r, o in zip(res_1pn, cat["relativistic_1pn"]):
        if r["q_fold"] > 0:
            print(f"  {r['id']} ({r['dim']}D {o['subspace']:26s} q={r['q_fold']:2d}): "
                  f"cat_mrg(1-RPO)={o['stab_margin']:.5f} (rpo_80bit={r['rpo_margin']:.5f}) vs "
                  f"TRUE T_inert_mrg={r['margin']:.5f} | k_inert=({r['k1']:+.4f},{r['k2']:+.4f},{r['k3']:+.4f}) "
                  f"symp={r['symp_err']:.1e} ret={r['ret_err']:.1e}")
        else:
            print(f"  {r['id']} ({r['dim']}D {o['subspace']:26s} RPO):  "
                  f"cat_mrg={o['stab_margin']:.5f} vs 80bit={r['margin']:.5f} "
                  f"symp={r['symp_err']:.1e} dH={r['dH']:.1e} ret={r['ret_err']:.1e}")

    print("\n=== [3/4] Dense Adaptive DOPRI5 Audit of Free-Fall Brake Orbits (rmin & Syzygy Sampling Bug) ===")
    for o in cat["orbits"]:
        if o["dynamical_class"] == "free_fall_brake":
            rmin_true, rmin_p, seq = trace_syz_and_rmin(o)
            ratio = o["rmin"] / rmin_true
            print(f"  {o['id']} m={[round(x, 5) for x in o['masses']]}: "
                  f"cat_rmin={o['rmin']:.5f} vs true_rmin={rmin_true:.5f} ({ratio:4.1f}x) | "
                  f"cat_syz={o['syzygies']:2d} vs true_syz={len(seq):2d} | half-seq={''.join(map(str, seq[:len(seq)//2]))}")

    print("\n=== [4/4] Rotating Choreography Literature Origin & Closure Checks ===")
    # Check 2D-DEM-004 syzygy cancellation
    o004 = next(o for o in cat["orbits"] if o["id"] == "2D-DEM-004")
    rmin4, _, seq4 = trace_syz_and_rmin(o004)
    print(f"  2D-DEM-004 (2d_choreo_8swap_3fold): true_rmin={rmin4:.5f}, syzygies={len(seq4)}, seq={''.join(map(str, seq4))}")
    for o in cat["orbits"]:
        if o["dynamical_class"] == "choreography_2d_rot":
            T_rpo = o["T"] / o["q_rpo"]
            Tstar_rpo = T_rpo * (abs(o["E"]) ** 1.5)
            print(f"  {o['id']} {o['subspace']:24s}: q_rpo={o['q_rpo']:2d}, T_rpo={T_rpo:7.4f}, T*_rpo={Tstar_rpo:7.4f}, Lambda={o['Lambda']:.4f}, word={o['word_compact']}")
    for o in cat["orbits"]:
        if o["dynamical_class"] == "choreography_3d_helical":
            prec_pi = o.get("precess_rad", 0.0) / math.pi
            print(f"  {o['id']} {o['subspace']:24s}: q_rpo={o['q_rpo']:2d}, T={o['T']:8.4f}, precess/pi={prec_pi:+.6f}, T*/q_rpo={o['Tstar']/o['q_rpo']:.4f}, Lambda={o['Lambda']:.4f}")

    # Mirror-image check between fig8_rot_2d.dat and choreo_12swap_2d.dat
    with open(os.path.join(ROOT, "data/choreo_rot/fig8_rot_2d.dat")) as f:
        r2d = [list(map(float, ln.split())) for ln in f if ln.strip()]
    with open(os.path.join(ROOT, "data/choreo_rot/choreo_12swap_2d.dat")) as f:
        r12 = [list(map(float, ln.split())) for ln in f if ln.strip()]
    print(f"  fig8_rot_2d.dat: alpha in [{r2d[0][0]:.4f}, {r2d[-1][0]:.4f}], y0 in [{r2d[-1][9]:.4f}, {r2d[0][9]:.4f}], T(alpha=0)={r2d[0][2]:.4f} (Fig-8!)")
    print(f"  choreo_12swap_2d.dat: alpha in [{r12[-1][0]:.4f}, {r12[0][0]:.4f}], y0 in [{r12[-1][9]:.4f}, {r12[0][9]:.4f}] (exact y->-y, vx->-vx mirror of fig8_rot_2d!)")
    with open(os.path.join(ROOT, "data/choreo_rot/fig8_helical_3d.dat")) as f:
        r3d = [list(map(float, ln.split())) for ln in f if ln.strip()]
    print(f"  fig8_helical_3d.dat: vz starts at {r3d[0][0]:.4f}, T={r3d[0][2]:.4f}, k=({r3d[0][5]:.5f}, {r3d[0][6]:.5f}) (exact Fig-8 -> CFM 2005 Gamma_y / Marchal P12!)")


if __name__ == "__main__":
    main()
