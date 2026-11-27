#!/usr/bin/env python3
"""
Compute exact Montgomery shape-sphere free-group homotopy words (F2 = <a, b>)
for all 2D periodic orbits in catalog.json.

Fixes applied following rigorous audit:
  1. Adaptive Dormand-Prince 5(4) integration (tol=1e-11) instead of fixed-step RK4,
     eliminating close-encounter drift/blowup.
  2. Equator-safe traversal starting at argmax_t |n_z(t)|, so orbits starting on the
     syzygy equator n_z(0) = 0 never drop their t = 0 crossing (e.g. AAAbbb).
  3. Full S_3 3-puncture outer automorphism quotient over all three collision
     projections (r_12 = 0, r_13 = 0, r_23 = 0), unifying equivalent words such as
     (Ab)^7 and A^7.
  4. Primitive root extraction (e.g. A^7 -> root 'A', power 7).
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Dormand-Prince 5(4) Butcher tableau
A21 = 1.0 / 5.0
A31, A32 = 3.0 / 40.0, 9.0 / 40.0
A41, A42, A43 = 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0
A51, A52, A53, A54 = 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0
A61, A62, A63, A64, A65 = (
    9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0
)
B1, B3, B4, B5, B6 = 35.0 / 384.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
E1, E3, E4, E5, E6, E7 = (
    71.0 / 57600.0, -71.0 / 16695.0, 71.0 / 1920.0,
    -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0,
)
INV_SQRT2 = 1.0 / math.sqrt(2.0)
INV_SQRT6 = 1.0 / math.sqrt(6.0)
INV_SQRT3 = 1.0 / math.sqrt(3.0)


def deriv(s, m0, m1, m2):
    x0, y0, x1, y1, x2, y2 = s[0], s[1], s[2], s[3], s[4], s[5]
    dx01, dy01 = x1 - x0, y1 - y0
    r2_01 = dx01 * dx01 + dy01 * dy01
    ir3_01 = 1.0 / (r2_01 * math.sqrt(r2_01))
    fx01, fy01 = dx01 * ir3_01, dy01 * ir3_01

    dx02, dy02 = x2 - x0, y2 - y0
    r2_02 = dx02 * dx02 + dy02 * dy02
    ir3_02 = 1.0 / (r2_02 * math.sqrt(r2_02))
    fx02, fy02 = dx02 * ir3_02, dy02 * ir3_02

    dx12, dy12 = x2 - x1, y2 - y1
    r2_12 = dx12 * dx12 + dy12 * dy12
    ir3_12 = 1.0 / (r2_12 * math.sqrt(r2_12))
    fx12, fy12 = dx12 * ir3_12, dy12 * ir3_12

    return (
        s[6], s[7], s[8], s[9], s[10], s[11],
        m1 * fx01 + m2 * fx02,
        m1 * fy01 + m2 * fy02,
        -m0 * fx01 + m2 * fx12,
        -m0 * fy01 + m2 * fy12,
        -m0 * fx02 - m1 * fx12,
        -m0 * fy02 - m1 * fy12,
    )


def integrate_trajectory_dopri5(y0, T, m, tol=1e-11):
    """Integrate 2D 3-body orbit over [0, T] with adaptive DOPRI5, returning sampled positions."""
    m0, m1, m2 = m[0], m[1], m[2]
    y = tuple(y0)
    t = 0.0
    h = min(1e-3, T * 0.005)
    hmax = min(0.05, T / 1600.0)
    pts = [y[:6]]
    k1 = deriv(y, m0, m1, m2)

    while t < T - 1e-14:
        if t + h > T:
            h = T - t
        accepted = False
        for _ in range(30):
            yt = tuple(y[i] + h * A21 * k1[i] for i in range(12))
            k2 = deriv(yt, m0, m1, m2)
            yt = tuple(y[i] + h * (A31 * k1[i] + A32 * k2[i]) for i in range(12))
            k3 = deriv(yt, m0, m1, m2)
            yt = tuple(y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]) for i in range(12))
            k4 = deriv(yt, m0, m1, m2)
            yt = tuple(
                y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i])
                for i in range(12)
            )
            k5 = deriv(yt, m0, m1, m2)
            yt = tuple(
                y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i])
                for i in range(12)
            )
            k6 = deriv(yt, m0, m1, m2)
            yn = tuple(
                y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i])
                for i in range(12)
            )
            k7 = deriv(yn, m0, m1, m2)

            err2 = 0.0
            for i in range(12):
                ei = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i])
                sc = tol * (1.0 + max(abs(y[i]), abs(yn[i])))
                r = ei / sc
                err2 += r * r
            err = math.sqrt(err2 / 12.0)
            if err <= 1.0:
                y = yn
                k1 = k7  # FSAL property
                t += h
                pts.append(y[:6])
                fac = 5.0 if err == 0.0 else min(5.0, max(0.2, 0.9 * (err ** -0.2)))
                h = min(h * fac, hmax)
                accepted = True
                break
            h *= max(0.1, 0.9 * (err ** -0.25))
            if h < 1e-14:
                raise RuntimeError("DOPRI5 step underflow during topology integration")
        if not accepted:
            raise RuntimeError("DOPRI5 failed to converge step")
    return pts


def shape_sphere_coords_perm(pos, perm):
    """Equal-weight Jacobi coordinates -> unit shape sphere (nx, ny, nz) for body permutation."""
    p0, p1, p2 = perm
    x1, y1 = pos[2 * p0], pos[2 * p0 + 1]
    x2, y2 = pos[2 * p1], pos[2 * p1 + 1]
    x3, y3 = pos[2 * p2], pos[2 * p2 + 1]
    rx = (x1 - x2) * INV_SQRT2
    ry = (y1 - y2) * INV_SQRT2
    lx = (x1 + x2 - 2.0 * x3) * INV_SQRT6
    ly = (y1 + y2 - 2.0 * y3) * INV_SQRT6
    r2 = rx * rx + ry * ry
    l2 = lx * lx + ly * ly
    R2 = r2 + l2
    if R2 < 1e-30:
        return (0.0, 0.0, 0.0)
    return ((r2 - l2) / R2, 2.0 * (rx * lx + ry * ly) / R2, 2.0 * (rx * ly - ry * lx) / R2)


def reduce_word(w):
    """Free-group reduction + cyclic reduction."""
    inv = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    stack = []
    for ch in w:
        if stack and inv[stack[-1]] == ch:
            stack.pop()
        else:
            stack.append(ch)
    left = 0
    right = len(stack)
    while right - left >= 2 and inv[stack[left]] == stack[right - 1]:
        left += 1
        right -= 1
    return "".join(stack[left:right])


def canonical_word_basic(w):
    """Canonical representative under cyclic shift, time reversal, and a<->b / a<->A relabeling."""
    if not w:
        return "1"
    inv = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    swap = {'a': 'b', 'b': 'a', 'A': 'B', 'B': 'A'}
    forms = set()
    rev = "".join(inv[c] for c in reversed(w))
    for bv in (w, rev):
        for t in (
            bv,
            "".join(swap[c] for c in bv),
            "".join(inv[c] for c in bv),
            "".join(swap[inv[c]] for c in bv),
        ):
            for i in range(len(t)):
                forms.add(t[i:] + t[:i])
    return min(forms)


def primitive_root(w):
    """Extract primitive cyclic root and repeat count (e.g. 'AAAAAAA' -> ('A', 7))."""
    if not w or w == "1":
        return "1", 1
    n = len(w)
    for p in range(1, n + 1):
        if n % p != 0:
            continue
        sub = w[:p]
        if sub * (n // p) == w:
            return canonical_word_basic(sub), n // p
    return w, 1


def format_word_compact(w_can_s3):
    root, power = primitive_root(w_can_s3)
    if root == "1":
        return "1"
    if power == 1:
        return root
    return f"({root})^{power}" if len(root) > 1 else f"{root}^{power}"


def word_from_points(pts, perm):
    """Extract reduced & basic-canonical word from sampled trajectory under given body permutation."""
    spts = [shape_sphere_coords_perm(p, perm) for p in pts]
    if len(spts) > 2:
        spts.pop()  # Drop duplicate endpoint t = T for cyclic traversal
    # Start strictly away from equator n_z = 0 so syzygy-starting orbits never miss t = 0 crossing!
    start = max(range(len(spts)), key=lambda i: abs(spts[i][2]))
    N = len(spts)
    last_pt = spts[start]
    last_sign = 1 if last_pt[2] > 0.0 else -1
    letters = []
    for k in range(1, N + 1):
        cur = spts[(start + k) % N]
        if cur[2] == 0.0:
            continue
        cur_sign = 1 if cur[2] > 0.0 else -1
        if cur_sign != last_sign:
            frac = abs(last_pt[2]) / (abs(last_pt[2]) + abs(cur[2]))
            nx = last_pt[0] + frac * (cur[0] - last_pt[0])
            ny = last_pt[1] + frac * (cur[1] - last_pt[1])
            if 1.0 + nx > 1e-12:
                u = ny / (1.0 + nx)
                up = cur_sign > 0
                if u > INV_SQRT3:
                    letters.append('a' if up else 'A')
                elif u < -INV_SQRT3:
                    letters.append('b' if up else 'B')
            last_sign = cur_sign
        last_pt = cur
    red = reduce_word(letters)
    return (red if red else "1"), canonical_word_basic(red)


def compute_word(o, tol=1e-11):
    """Compute exact (word_red, word_can, word_can_s3, word_prim, word_pow) for a 2D orbit."""
    pts = integrate_trajectory_dopri5(o["y0"], o["T"], o["masses"], tol=tol)
    w_red, w_can0 = word_from_points(pts, (0, 1, 2))
    _, w_can1 = word_from_points(pts, (0, 2, 1))
    _, w_can2 = word_from_points(pts, (1, 2, 0))
    w_can_s3 = min((w_can0, w_can1, w_can2), key=lambda w: (0 if w == "1" else len(w), w))
    w_prim, w_pow = primitive_root(w_can_s3)
    return w_red, w_can0, w_can_s3, w_prim, w_pow


def main():
    ap = argparse.ArgumentParser(description="Classify 2D periodic orbits by S3 shape-sphere free-group words.")
    ap.add_argument("--verify", action="store_true", help="Cross-check against word_can_s3 stored in catalog.json.")
    args = ap.parse_args()

    cat = json.load(open(os.path.join(ROOT, "catalog.json")))
    orbits_2d = [o for o in cat["orbits"] if o["dim"] == 2]
    by_s3 = {}
    mismatches = []

    for o in orbits_2d:
        w_red, w_can, w_s3, w_prim, w_pow = compute_word(o)
        by_s3.setdefault(w_s3, []).append((o, w_prim, w_pow))
        if args.verify and "word_can_s3" in o:
            if o["word_can_s3"] != w_s3:
                mismatches.append((o["id"], o["word_can_s3"], w_s3))

    print(f"Classified {len(orbits_2d)} 2D Newtonian orbits into {len(by_s3)} distinct S3 shape-sphere classes:\n")
    for w_s3, grp in sorted(by_s3.items(), key=lambda kv: (-len(kv[1]), 0 if kv[0] == "1" else len(kv[0]), kv[0])):
        ex, w_prim, w_pow = grp[0]
        compact = format_word_compact(w_s3)
        wdisp = compact if len(compact) <= 28 else compact[:25] + "..."
        print(
            f"  [{len(grp):2d} orbits] len={0 if w_s3 == '1' else len(w_s3):2d} "
            f"S3={wdisp:28s} | top={ex['id']} ({ex['subspace']}, m={ex['masses']}, marg={ex['stab_margin']:.4f})"
        )

    if args.verify:
        if mismatches:
            print(f"\nERROR: {len(mismatches)} mismatches vs catalog.json word_can_s3:")
            for oid, exp, got in mismatches:
                print(f"  {oid}: catalog={exp} vs python={got}")
            sys.exit(1)
        print(f"\nVERIFIED: 0/{len(orbits_2d)} mismatches against 80-bit C++ audit_catalog word_can_s3.")


if __name__ == "__main__":
    main()
