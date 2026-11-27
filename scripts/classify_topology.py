#!/usr/bin/env python3
"""
Compute exact Montgomery shape-sphere free-group homotopy words (F2 = <a, b>)
for all 2D democratic periodic orbits in catalog.json.

On the shape sphere (Jacobi coordinates rho, lambda -> Hopf map to S^2),
the three binary collisions lie on the equator (n_z = 0) at 120-degree intervals.
Stereographic projection from one collision (say r_12 = 0) sends the other two
collisions to two punctures in the plane, and the equatorial crossings yield
the standard free-group word in {a, A, b, B} (Šuvakov & Dmitrašinović 2013).
"""
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def deriv(s, m):
    out = [s[6], s[7], s[8], s[9], s[10], s[11], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(3):
        for j in range(i + 1, 3):
            dx = s[2*j] - s[2*i]
            dy = s[2*j + 1] - s[2*i + 1]
            r2 = dx*dx + dy*dy
            ir3 = 1.0 / (r2 * math.sqrt(r2))
            out[6 + 2*i] += m[j] * dx * ir3
            out[7 + 2*i] += m[j] * dy * ir3
            out[6 + 2*j] -= m[i] * dx * ir3
            out[7 + 2*j] -= m[i] * dy * ir3
    return out


def rk4_step(y, h, m):
    k1 = deriv(y, m)
    y2 = [y[i] + 0.5*h*k1[i] for i in range(12)]
    k2 = deriv(y2, m)
    y3 = [y[i] + 0.5*h*k2[i] for i in range(12)]
    k3 = deriv(y3, m)
    y4 = [y[i] + h*k3[i] for i in range(12)]
    k4 = deriv(y4, m)
    return [y[i] + (h/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(12)]


def shape_sphere_coords(y):
    """Equal-weight Jacobi coordinates -> unit shape sphere (nx, ny, nz)."""
    x1, y1, x2, y2, x3, y3 = y[0], y[1], y[2], y[3], y[4], y[5]
    # rho = (r1 - r2)/sqrt(2), lam = (r1 + r2 - 2*r3)/sqrt(6)
    rx = (x1 - x2) / math.sqrt(2.0)
    ry = (y1 - y2) / math.sqrt(2.0)
    lx = (x1 + x2 - 2.0*x3) / math.sqrt(6.0)
    ly = (y1 + y2 - 2.0*y3) / math.sqrt(6.0)
    r2 = rx*rx + ry*ry
    l2 = lx*lx + ly*ly
    R2 = r2 + l2
    if R2 < 1e-15:
        return (0.0, 0.0, 0.0)
    nx = (r2 - l2) / R2
    ny = 2.0 * (rx*lx + ry*ly) / R2
    nz = 2.0 * (rx*ly - ry*lx) / R2
    return (nx, ny, nz)


def reduce_word(w):
    """Free-group reduction + cyclic reduction."""
    inv = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    stack = []
    for ch in w:
        if stack and inv[stack[-1]] == ch:
            stack.pop()
        else:
            stack.append(ch)
    # Cyclic reduction
    while len(stack) >= 2 and inv[stack[0]] == stack[-1]:
        stack.pop(0)
        stack.pop()
    return "".join(stack)


def canonical_word(w):
    """Canonical representative under cyclic shift, inversion, and a<->b relabeling."""
    if not w:
        return "1"
    inv = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    swap = {'a': 'b', 'b': 'a', 'A': 'B', 'B': 'A'}
    flip = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    forms = set()
    base_variants = [
        w,
        "".join(inv[c] for c in reversed(w)),
    ]
    for bv in base_variants:
        for t in (
            bv,
            "".join(swap[c] for c in bv),
            "".join(flip[c] for c in bv),
            "".join(swap[flip[c]] for c in bv),
        ):
            for i in range(len(t)):
                forms.add(t[i:] + t[:i])
    return min(forms)


def compute_word(o):
    m = o["masses"]
    y = list(o["y0"])
    T = o["T"]
    N = max(2000, int(T * 120))
    dt = T / N

    # Stereographic projection from collision r12=0 (which is nx=-1, ny=0, nz=0)
    # Wait: r12=0 => rho=0 => nx=-1, ny=0, nz=0!
    # The other two collisions on the equator nz=0 are at:
    # r13=0 => (nx=+1/2, ny=+sqrt(3)/2, nz=0)
    # r23=0 => (nx=+1/2, ny=-sqrt(3)/2, nz=0)
    # Stereographic projection from (-1, 0, 0) onto (u, v) = (ny/(1+nx), nz/(1+nx)):
    #   r13=0 -> u = +1/sqrt(3) = +0.57735, v = 0
    #   r23=0 -> u = -1/sqrt(3) = -0.57735, v = 0
    # Crossings of v = 0 (equator nz = 0):
    #   u > +1/sqrt(3): outside right -> no generator (or basepoint)
    #   -1/sqrt(3) < u < +1/sqrt(3): between punctures -> 'a' (if dv>0) / 'A' (if dv<0)
    #   u < -1/sqrt(3): left of both -> 'b' (if dv>0) / 'B' (if dv<0)
    # Wait: in F2 with two punctures at P1=(+c, 0) and P2=(-c, 0),
    # crossing the real axis v=0 with u in (-c, +c) separates P1 from P2,
    # and crossing u < -c goes around P2 (or u > +c goes around P1).
    # Standard branch cuts: (-inf, -c] for 'b'/'B', and [+c, +inf) for 'a'/'A'!
    c = 1.0 / math.sqrt(3.0)
    prev = shape_sphere_coords(y)
    letters = []

    for _ in range(N):
        # Substep if close encounter
        r12 = math.hypot(y[2]-y[0], y[3]-y[1])
        r13 = math.hypot(y[4]-y[0], y[5]-y[1])
        r23 = math.hypot(y[4]-y[2], y[5]-y[3])
        rmin = min(r12, r13, r23)
        sub = 1 if rmin > 0.1 else (8 if rmin > 0.02 else 32)
        h = dt / sub
        for _s in range(sub):
            y = rk4_step(y, h, m)
        cur = shape_sphere_coords(y)
        if prev[2] * cur[2] < 0.0:
            # Linear interpolation for equator crossing nz = 0
            frac = abs(prev[2]) / (abs(prev[2]) + abs(cur[2]))
            nx = prev[0] + frac * (cur[0] - prev[0])
            ny = prev[1] + frac * (cur[1] - prev[1])
            if 1.0 + nx > 1e-9:
                u = ny / (1.0 + nx)
                up = cur[2] > prev[2]
                if u > c:
                    letters.append('a' if up else 'A')
                elif u < -c:
                    letters.append('b' if up else 'B')
        prev = cur

    red = reduce_word(letters)
    canon = canonical_word(red)
    return red, canon


def main():
    cat = json.load(open(os.path.join(ROOT, "catalog.json")))
    dem = [o for o in cat["orbits"] if o["dim"] == 2 and o["democratic"]]
    by_canon = {}
    for o in dem:
        red, canon = compute_word(o)
        o["word_reduced"] = red
        o["word_canonical"] = canon
        by_canon.setdefault(canon, []).append(o)

    print(f"Classified {len(dem)} democratic 2D steady orbits into {len(by_canon)} distinct topological free-group classes:\n")
    for canon, grp in sorted(by_canon.items(), key=lambda kv: (-len(kv[1]), len(kv[0]))):
        ex = grp[0]
        wdisp = canon if len(canon) <= 28 else canon[:25] + "..."
        masses = sorted({tuple(x["masses"]) for x in grp})
        print(f"  [{len(grp):2d} orbits] len={len(canon):2d} word={wdisp:28s} | top={ex['id']} ({ex['subspace']}, m={ex['masses']}, marg={ex['stab_margin']:.3f})")


if __name__ == "__main__":
    main()
