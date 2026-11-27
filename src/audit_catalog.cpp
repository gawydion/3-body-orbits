#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Production 80-bit (long double) analytical variational monodromy, RPO/cover,
// and Montgomery S3 shape-sphere topology auditor for Newtonian 3-body orbits.
//
// Modes:
//   audit_catalog [--threads N]
//     Reads lines: id dim subspace m1 m2 m3 T precess_rad y0[0..17]
//     (Pass precess_rad = 0 for inertial closures; nonzero ONLY for rotating-frame RPOs.)
//   audit_catalog --collinear-ultra [--threads N]
//     Reads lines: id m3 vx vy T
//     Exploits half-period spatial involution P(theta) * y(T/2) = y(0) (P^2 = I)
//     at tol = 1e-17L so M(T) = (P * M(T/2))^2 and k_j = h_j^2 - 2.

using ld = long double;
using cplx = std::complex<ld>;
using Masses = std::array<ld, 3>;
using State3D = std::array<ld, 18>;
using State2D = std::array<ld, 12>;

static constexpr ld PI_L = 3.141592653589793238462643383279502884L;

// =========================================================================
// Hessenberg + Francis double-shift QR eigenvalues in 80-bit long double
// =========================================================================
static std::vector<cplx> eigenvalues_qr(std::vector<std::vector<ld>> A) {
    int n = (int)A.size();
    for (int k = 0; k < n - 2; ++k) {
        ld norm2 = 0;
        for (int i = k + 1; i < n; ++i) norm2 += A[i][k] * A[i][k];
        if (norm2 == 0) continue;
        ld alpha = -std::copysign(std::sqrt(norm2), A[k + 1][k]);
        std::vector<ld> v(n - k - 1, 0);
        v[0] = A[k + 1][k] - alpha;
        for (int i = 1; i < n - k - 1; ++i) v[i] = A[k + 1 + i][k];
        ld v2 = 0;
        for (ld x : v) v2 += x * x;
        if (v2 == 0) continue;
        ld beta = 2.0L / v2;
        for (int j = k; j < n; ++j) {
            ld dot = 0;
            for (int i = 0; i < n - k - 1; ++i) dot += v[i] * A[k + 1 + i][j];
            for (int i = 0; i < n - k - 1; ++i) A[k + 1 + i][j] -= beta * dot * v[i];
        }
        for (int i = 0; i < n; ++i) {
            ld dot = 0;
            for (int j = 0; j < n - k - 1; ++j) dot += A[i][k + 1 + j] * v[j];
            for (int j = 0; j < n - k - 1; ++j) A[i][k + 1 + j] -= beta * dot * v[j];
        }
        for (int i = k + 2; i < n; ++i) A[i][k] = 0;
    }

    std::vector<cplx> eigs;
    int m = n, iter = 0;
    const ld eps = 1e-17L;
    while (m > 0) {
        if (m == 1) { eigs.emplace_back(A[0][0], 0); break; }
        int l = m - 1;
        while (l > 0) {
            ld s = std::abs(A[l - 1][l - 1]) + std::abs(A[l][l]);
            if (s == 0) s = 1.0L;
            if (std::abs(A[l][l - 1]) <= eps * s) { A[l][l - 1] = 0; break; }
            --l;
        }
        if (l == m - 1) { eigs.emplace_back(A[m - 1][m - 1], 0); --m; iter = 0; continue; }
        if (l == m - 2) {
            ld a = A[m - 2][m - 2], b = A[m - 2][m - 1];
            ld c = A[m - 1][m - 2], d = A[m - 1][m - 1];
            ld tr = a + d, det = a * d - b * c;
            ld disc = tr * tr - 4.0L * det;
            if (disc >= 0) {
                ld sd = std::sqrt(disc);
                eigs.emplace_back(0.5L * (tr + sd), 0);
                eigs.emplace_back(0.5L * (tr - sd), 0);
            } else {
                ld sd = std::sqrt(-disc);
                eigs.emplace_back(0.5L * tr, 0.5L * sd);
                eigs.emplace_back(0.5L * tr, -0.5L * sd);
            }
            m -= 2; iter = 0; continue;
        }
        ++iter;
        if (iter > 2500) {
            for (int i = 0; i < m; ++i) eigs.emplace_back(A[i][i], 0);
            break;
        }
        ld s, t;
        if (iter % 15 == 0) {
            ld w = std::abs(A[m - 1][m - 2]) + std::abs(A[m - 2][m - 3]);
            s = 1.5L * w; t = w * w;
        } else {
            s = A[m - 2][m - 2] + A[m - 1][m - 1];
            t = A[m - 2][m - 2] * A[m - 1][m - 1] - A[m - 2][m - 1] * A[m - 1][m - 2];
        }
        ld x = A[l][l] * A[l][l] + A[l][l + 1] * A[l + 1][l] - s * A[l][l] + t;
        ld y = A[l + 1][l] * (A[l][l] + A[l + 1][l + 1] - s);
        ld z = (l + 2 < m) ? A[l + 1][l] * A[l + 2][l + 1] : 0;

        for (int k = l; k <= m - 2; ++k) {
            int sz = std::min(3, m - k);
            ld v[3] = {x, y, (sz == 3 ? z : 0.0L)};
            ld norm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (norm == 0) break;
            v[0] += std::copysign(norm, v[0]);
            ld v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if (v2 == 0) break;
            ld beta = 2.0L / v2;
            int col_start = std::max(l, k - 1);
            for (int j = col_start; j < n; ++j) {
                ld dot = 0;
                for (int i = 0; i < sz; ++i) dot += v[i] * A[k + i][j];
                for (int i = 0; i < sz; ++i) A[k + i][j] -= beta * dot * v[i];
            }
            int row_end = std::min(m, k + 4);
            for (int i = 0; i < row_end; ++i) {
                ld dot = 0;
                for (int j = 0; j < sz; ++j) dot += A[i][k + j] * v[j];
                for (int j = 0; j < sz; ++j) A[i][k + j] -= beta * dot * v[j];
            }
            if (k < m - 2) {
                x = A[k + 1][k]; y = A[k + 2][k];
                z = (k + 3 < m) ? A[k + 3][k] : 0;
            }
        }
    }
    return eigs;
}

static ld min_eig_sym_3x3(const ld C[3][3]) {
    ld p1 = C[0][1]*C[0][1] + C[0][2]*C[0][2] + C[1][2]*C[1][2];
    ld q = (C[0][0] + C[1][1] + C[2][2]) / 3.0L;
    if (p1 < 1e-32L) return std::min({C[0][0], C[1][1], C[2][2]});
    ld p2 = (C[0][0]-q)*(C[0][0]-q) + (C[1][1]-q)*(C[1][1]-q) + (C[2][2]-q)*(C[2][2]-q) + 2.0L*p1;
    ld p = std::sqrt(p2 / 6.0L);
    ld B[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        B[i][j] = (C[i][j] - (i == j ? q : 0.0L)) / p;
    ld detB = B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
            - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
            + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);
    ld r = std::clamp(0.5L * detB, -1.0L, 1.0L);
    ld phi = std::acos(r) / 3.0L;
    ld e0 = q + 2.0L * p * std::cos(phi);
    ld e1 = q + 2.0L * p * std::cos(phi + 2.0L * PI_L / 3.0L);
    ld e2 = 3.0L * q - e0 - e1;
    return std::max(0.0L, std::min({e0, e1, e2}));
}

// =========================================================================
// 3D Gravitational N=3 Dynamics + Analytical Variational Equations
// =========================================================================
static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, ld* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    ld a[9] = {0}, rmin2 = 1e300L;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            ld dx = s[3 * j] - s[3 * i];
            ld dy = s[3 * j + 1] - s[3 * i + 1];
            ld dz = s[3 * j + 2] - s[3 * i + 2];
            ld r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rmin2) rmin2 = r2;
            ld ir = 1.0L / std::sqrt(r2);
            ld ir3 = ir * ir * ir;
            a[3 * i] += m[j] * dx * ir3;
            a[3 * i + 1] += m[j] * dy * ir3;
            a[3 * i + 2] += m[j] * dz * ir3;
            a[3 * j] -= m[i] * dx * ir3;
            a[3 * j + 1] -= m[i] * dy * ir3;
            a[3 * j + 2] -= m[i] * dz * ir3;
        }
    }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline void compute_accel_jacobian(const State3D& s, const Masses& m, ld H[9][9]) {
    for (int i = 0; i < 9; ++i) for (int j = 0; j < 9; ++j) H[i][j] = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            ld d[3] = {s[3 * j] - s[3 * i], s[3 * j + 1] - s[3 * i + 1], s[3 * j + 2] - s[3 * i + 2]};
            ld r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            ld ir = 1.0L / std::sqrt(r2);
            ld ir3 = ir * ir * ir;
            ld ir5 = ir3 * ir * ir;
            for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                    ld Tpq = ((p == q) ? ir3 : 0.0L) - 3.0L * d[p] * d[q] * ir5;
                    H[3 * i + p][3 * j + q] += m[j] * Tpq;
                    H[3 * i + p][3 * i + q] -= m[j] * Tpq;
                    H[3 * j + p][3 * i + q] += m[i] * Tpq;
                    H[3 * j + p][3 * j + q] -= m[i] * Tpq;
                }
            }
        }
    }
}

struct VarState { State3D y; ld M[18][18]; };

static inline void deriv_var(const VarState& vs, const Masses& m, VarState& dvs, ld* rmin_out = nullptr) {
    deriv3d(vs.y, m, dvs.y, rmin_out);
    ld H[9][9];
    compute_accel_jacobian(vs.y, m, H);
    for (int col = 0; col < 18; ++col) {
        for (int r = 0; r < 9; ++r) dvs.M[r][col] = vs.M[9 + r][col];
        for (int r = 0; r < 9; ++r) {
            ld sum = 0;
            for (int k = 0; k < 9; ++k) sum += H[r][k] * vs.M[k][col];
            dvs.M[9 + r][col] = sum;
        }
    }
}

static bool dopri5_step_state(State3D& y, ld& t, ld& h, ld t_end, const Masses& m, ld tol, ld& rmin_track) {
    static constexpr ld a21 = 1.0L / 5.0L, a31 = 3.0L / 40.0L, a32 = 9.0L / 40.0L,
        a41 = 44.0L / 45.0L, a42 = -56.0L / 15.0L, a43 = 32.0L / 9.0L,
        a51 = 19372.0L / 6561.0L, a52 = -25360.0L / 2187.0L, a53 = 64448.0L / 6561.0L, a54 = -212.0L / 729.0L,
        a61 = 9017.0L / 3168.0L, a62 = -355.0L / 33.0L, a63 = 46732.0L / 5247.0L, a64 = 49.0L / 176.0L, a65 = -5103.0L / 18656.0L,
        b1 = 35.0L / 384.0L, b3 = 500.0L / 1113.0L, b4 = 125.0L / 192.0L, b5 = -2187.0L / 6784.0L, b6 = 11.0L / 84.0L,
        e1 = 71.0L / 57600.0L, e3 = -71.0L / 16695.0L, e4 = 71.0L / 1920.0L, e5 = -17253.0L / 339200.0L, e6 = 22.0L / 525.0L, e7 = -1.0L / 40.0L;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0) return true;
    State3D k1, k2, k3, k4, k5, k6, k7, yt, yn;
    ld rm = 1e300L;
    for (int att = 0; att < 30; ++att) {
        deriv3d(y, m, k1, &rm);
        if (rm < 1e-8L) return false;
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * a21 * k1[i];
        deriv3d(yt, m, k2);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        deriv3d(yt, m, k3);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        deriv3d(yt, m, k4);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        deriv3d(yt, m, k5);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        deriv3d(yt, m, k6);
        for (int i = 0; i < 18; ++i) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        deriv3d(yn, m, k7, &rm);
        if (rm < 1e-8L) return false;
        ld err = 0;
        for (int i = 0; i < 18; ++i) {
            ld ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            ld sc = tol * (1.0L + std::max(std::abs(y[i]), std::abs(yn[i])));
            ld r = ei / sc;
            err += r * r;
        }
        err = std::sqrt(err / 18.0L);
        if (err <= 1.0L) {
            y = yn; t += h;
            if (rm < rmin_track) rmin_track = rm;
            h = std::min(h * ((err == 0) ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.05L);
            return true;
        }
        h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
        if (h < 1e-16L) return false;
    }
    return false;
}

static bool integrate_state(State3D& y, ld T, const Masses& m, ld tol = 1e-15L, ld* rmin_out = nullptr) {
    ld t = 0, h = std::min(1e-3L, T * 0.01L), rmin_track = 1e300L;
    int steps = 0;
    while (t < T - 1e-17L && steps < 2000000) {
        if (!dopri5_step_state(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-14L;
}

static bool integrate_var(const State3D& y0, ld T, const Masses& m, State3D& yT, ld Mout[18][18], ld tol = 1e-14L) {
    static constexpr ld a21 = 1.0L / 5.0L, a31 = 3.0L / 40.0L, a32 = 9.0L / 40.0L,
        a41 = 44.0L / 45.0L, a42 = -56.0L / 15.0L, a43 = 32.0L / 9.0L,
        a51 = 19372.0L / 6561.0L, a52 = -25360.0L / 2187.0L, a53 = 64448.0L / 6561.0L, a54 = -212.0L / 729.0L,
        a61 = 9017.0L / 3168.0L, a62 = -355.0L / 33.0L, a63 = 46732.0L / 5247.0L, a64 = 49.0L / 176.0L, a65 = -5103.0L / 18656.0L,
        b1 = 35.0L / 384.0L, b3 = 500.0L / 1113.0L, b4 = 125.0L / 192.0L, b5 = -2187.0L / 6784.0L, b6 = 11.0L / 84.0L,
        e1 = 71.0L / 57600.0L, e3 = -71.0L / 16695.0L, e4 = 71.0L / 1920.0L, e5 = -17253.0L / 339200.0L, e6 = 22.0L / 525.0L, e7 = -1.0L / 40.0L;

    VarState vs{};
    vs.y = y0;
    for (int i = 0; i < 18; ++i) for (int j = 0; j < 18; ++j) vs.M[i][j] = (i == j ? 1.0L : 0.0L);

    auto axpy = [](const VarState& a, VarState& out, ld h, int nk, const VarState* ks[], const ld cs[]) {
        for (int i = 0; i < 18; ++i) {
            ld s = 0;
            for (int k = 0; k < nk; ++k) s += cs[k] * ks[k]->y[i];
            out.y[i] = a.y[i] + h * s;
        }
        for (int r = 0; r < 18; ++r) {
            for (int c = 0; c < 18; ++c) {
                ld s = 0;
                for (int k = 0; k < nk; ++k) s += cs[k] * ks[k]->M[r][c];
                out.M[r][c] = a.M[r][c] + h * s;
            }
        }
    };

    VarState k1, k2, k3, k4, k5, k6, k7, yt, yn;
    ld t = 0, h = std::min(1e-3L, T * 0.01L);
    int steps = 0;
    while (t < T - 1e-17L && steps < 2500000) {
        if (t + h > T) h = T - t;
        bool ok = false;
        for (int att = 0; att < 30; ++att) {
            ld rm = 1e300L;
            deriv_var(vs, m, k1, &rm);
            if (rm < 1e-8L) return false;
            { const VarState* ks[1] = {&k1}; ld cs[1] = {a21}; axpy(vs, yt, h, 1, ks, cs); }
            deriv_var(yt, m, k2);
            { const VarState* ks[2] = {&k1, &k2}; ld cs[2] = {a31, a32}; axpy(vs, yt, h, 2, ks, cs); }
            deriv_var(yt, m, k3);
            { const VarState* ks[3] = {&k1, &k2, &k3}; ld cs[3] = {a41, a42, a43}; axpy(vs, yt, h, 3, ks, cs); }
            deriv_var(yt, m, k4);
            { const VarState* ks[4] = {&k1, &k2, &k3, &k4}; ld cs[4] = {a51, a52, a53, a54}; axpy(vs, yt, h, 4, ks, cs); }
            deriv_var(yt, m, k5);
            { const VarState* ks[5] = {&k1, &k2, &k3, &k4, &k5}; ld cs[5] = {a61, a62, a63, a64, a65}; axpy(vs, yt, h, 5, ks, cs); }
            deriv_var(yt, m, k6);
            { const VarState* ks[5] = {&k1, &k3, &k4, &k5, &k6}; ld cs[5] = {b1, b3, b4, b5, b6}; axpy(vs, yn, h, 5, ks, cs); }
            deriv_var(yn, m, k7, &rm);
            if (rm < 1e-8L) return false;

            ld err = 0;
            for (int i = 0; i < 18; ++i) {
                ld ei = h * (e1 * k1.y[i] + e3 * k3.y[i] + e4 * k4.y[i] + e5 * k5.y[i] + e6 * k6.y[i] + e7 * k7.y[i]);
                ld sc = tol * (1.0L + std::max(std::abs(vs.y[i]), std::abs(yn.y[i])));
                ld r = ei / sc; err += r * r;
            }
            for (int r = 0; r < 18; ++r) {
                for (int c = 0; c < 18; ++c) {
                    ld ei = h * (e1 * k1.M[r][c] + e3 * k3.M[r][c] + e4 * k4.M[r][c] + e5 * k5.M[r][c] + e6 * k6.M[r][c] + e7 * k7.M[r][c]);
                    ld sc = tol * (1.0L + std::max(std::abs(vs.M[r][c]), std::abs(yn.M[r][c])));
                    ld rr = ei / sc; err += rr * rr;
                }
            }
            err = std::sqrt(err / 342.0L);
            if (err <= 1.0L) {
                vs = yn; t += h;
                h = std::min(h * ((err == 0) ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.04L);
                ok = true;
                break;
            }
            h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
            if (h < 1e-16L) return false;
        }
        if (!ok) return false;
        ++steps;
    }
    yT = vs.y;
    for (int i = 0; i < 18; ++i) for (int j = 0; j < 18; ++j) Mout[i][j] = vs.M[i][j];
    return t >= T - 1e-14L;
}

static inline void rot_vec(const ld u[3], ld ang, ld& x, ld& y, ld& z) {
    ld c = std::cos(ang), s = std::sin(ang), omc = 1.0L - c;
    ld dot = u[0] * x + u[1] * y + u[2] * z;
    ld cx = u[1] * z - u[2] * y, cy = u[2] * x - u[0] * z, cz = u[0] * y - u[1] * x;
    ld nx = x * c + cx * s + u[0] * dot * omc;
    ld ny = y * c + cy * s + u[1] * dot * omc;
    ld nz = z * c + cz * s + u[2] * dot * omc;
    x = nx; y = ny; z = nz;
}

static State3D rotate_state(const State3D& s, const ld u[3], ld ang) {
    State3D out = s;
    for (int b = 0; b < 6; ++b) rot_vec(u, ang, out[3 * b], out[3 * b + 1], out[3 * b + 2]);
    return out;
}

static State3D permute_state(const State3D& s, const int p[3]) {
    State3D out{};
    for (int b = 0; b < 3; ++b) {
        for (int d = 0; d < 3; ++d) {
            out[3 * b + d] = s[3 * p[b] + d];
            out[9 + 3 * b + d] = s[9 + 3 * p[b] + d];
        }
    }
    return out;
}

static State3D reflect_plane(const State3D& s, const ld u[3]) {
    State3D out = s;
    for (int b = 0; b < 6; ++b) {
        ld dot = u[0] * out[3 * b] + u[1] * out[3 * b + 1] + u[2] * out[3 * b + 2];
        out[3 * b] -= 2.0L * dot * u[0];
        out[3 * b + 1] -= 2.0L * dot * u[1];
        out[3 * b + 2] -= 2.0L * dot * u[2];
    }
    return out;
}

static ld optimal_rot_angle(const State3D& s, const State3D& target, const ld u[3], ld& min_dist) {
    ld A = 0, B = 0;
    for (int b = 0; b < 6; ++b) {
        ld sx = s[3 * b], sy = s[3 * b + 1], sz = s[3 * b + 2];
        ld tx = target[3 * b], ty = target[3 * b + 1], tz = target[3 * b + 2];
        ld us = u[0] * sx + u[1] * sy + u[2] * sz;
        ld ut = u[0] * tx + u[1] * ty + u[2] * tz;
        A += (tx * sx + ty * sy + tz * sz) - ut * us;
        ld cx = u[1] * sz - u[2] * sy;
        ld cy = u[2] * sx - u[0] * sz;
        ld cz = u[0] * sy - u[1] * sx;
        B -= (tx * cx + ty * cy + tz * cz);
    }
    ld theta = std::atan2(B, A);
    State3D rs = rotate_state(s, u, -theta);
    ld d2 = 0;
    for (int i = 0; i < 18; ++i) { ld d = rs[i] - target[i]; d2 += d * d; }
    min_dist = std::sqrt(d2);
    return theta;
}

static inline ld energy3d(const State3D& s, const Masses& m) {
    ld K = 0, U = 0;
    for (int i = 0; i < 3; ++i)
        K += 0.5L * m[i] * (s[9 + 3 * i] * s[9 + 3 * i] + s[10 + 3 * i] * s[10 + 3 * i] + s[11 + 3 * i] * s[11 + 3 * i]);
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            ld dx = s[3 * j] - s[3 * i], dy = s[3 * j + 1] - s[3 * i + 1], dz = s[3 * j + 2] - s[3 * i + 2];
            U -= m[i] * m[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    return K + U;
}

static inline std::array<ld, 3> angmom3d(const State3D& s, const Masses& m) {
    ld Lx = 0, Ly = 0, Lz = 0;
    for (int i = 0; i < 3; ++i) {
        Lx += m[i] * (s[3 * i + 1] * s[11 + 3 * i] - s[3 * i + 2] * s[10 + 3 * i]);
        Ly += m[i] * (s[3 * i + 2] * s[9 + 3 * i] - s[3 * i] * s[11 + 3 * i]);
        Lz += m[i] * (s[3 * i] * s[10 + 3 * i] - s[3 * i + 1] * s[9 + 3 * i]);
    }
    return {Lx, Ly, Lz};
}

static inline ld state_dist(const State3D& a, const State3D& b) {
    ld s = 0;
    for (int i = 0; i < 18; ++i) { ld d = a[i] - b[i]; s += d * d; }
    return std::sqrt(s);
}

// Extract exact Broucke stability indices via Newton trace invariants (immune to
// Jordan-block root-splitting at trivial eigenvalues!).
// For an RPO with rotation angle theta around u:
//   - 2D planar 12x12 has 4 trivial eigenvalues at +1 and 2 pairs at e^{+/-i*theta},
//     so tr_triv(Mp^m) = 4 + 4*cos(m*theta).
//   - 3D 18x18 has 6 trivial eigenvalues at +1 and 3 pairs at e^{+/-i*theta},
//     so tr_triv(M^m)  = 6 + 6*cos(m*theta).
static void broucke_2d_from_matrix(const std::vector<std::vector<ld>>& Mp, ld theta,
                                   bool& real_k, ld& k1, ld& k2, ld& margin) {
    ld tr1 = 0, tr2 = 0;
    for (int i = 0; i < 12; ++i) tr1 += Mp[i][i];
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j)
            tr2 += Mp[i][j] * Mp[j][i];
    ld s = tr1 - (4.0L + 4.0L * std::cos(theta));
    ld p2 = (tr2 - (4.0L + 4.0L * std::cos(2.0L * theta))) + 4.0L;
    ld p = 0.5L * (s * s - p2);
    ld disc = s * s - 4.0L * p;
    if (disc >= -1e-10L) {
        ld sd = std::sqrt(std::max(0.0L, disc));
        real_k = true;
        k1 = 0.5L * (s + sd);
        k2 = 0.5L * (s - sd);
        margin = 2.0L - std::max(std::abs(k1), std::abs(k2));
    } else {
        real_k = false;
        k1 = 0.5L * s;
        k2 = 0.5L * std::sqrt(-disc);
        margin = -10.0L;
    }
}

static void broucke_3d_from_matrix(const std::vector<std::vector<ld>>& M, ld theta,
                                   bool& real_k, ld& k1, ld& k2, ld& k3, ld& margin) {
    ld M2[18][18] = {{0}};
    ld tr1 = 0, tr2 = 0, tr3 = 0;
    for (int i = 0; i < 18; ++i) {
        tr1 += M[i][i];
        for (int j = 0; j < 18; ++j) {
            ld s = 0;
            for (int k = 0; k < 18; ++k) s += M[i][k] * M[k][j];
            M2[i][j] = s;
        }
        tr2 += M2[i][i];
    }
    for (int i = 0; i < 18; ++i)
        for (int j = 0; j < 18; ++j)
            tr3 += M2[i][j] * M[j][i];

    ld t1 = tr1 - (6.0L + 6.0L * std::cos(theta));
    ld t2 = tr2 - (6.0L + 6.0L * std::cos(2.0L * theta));
    ld t3 = tr3 - (6.0L + 6.0L * std::cos(3.0L * theta));
    ld p2 = t2 + 6.0L;
    ld p3 = t3 + 3.0L * t1;
    ld e1 = t1;
    ld e2 = 0.5L * (e1 * e1 - p2);
    ld e3 = (p3 - e1 * e1 * e1 + 3.0L * e1 * e2) / 3.0L;

    // Solve k^3 - e1*k^2 + e2*k - e3 = 0
    ld p = e2 - e1 * e1 / 3.0L;
    ld q = -2.0L * e1 * e1 * e1 / 27.0L + e1 * e2 / 3.0L - e3;
    ld D = 0.25L * q * q + (p / 3.0L) * (p / 3.0L) * (p / 3.0L);

    auto polish = [&](ld x) {
        for (int it = 0; it < 6; ++it) {
            ld f = ((x - e1) * x + e2) * x - e3;
            ld df = (3.0L * x - 2.0L * e1) * x + e2;
            if (std::abs(df) < 1e-22L) break;
            x -= f / df;
        }
        return x;
    };

    if (D <= 1e-14L) {
        real_k = true;
        ld r = std::sqrt(std::max(0.0L, -p / 3.0L));
        ld r3 = r * r * r;
        ld arg = (r3 > 1e-24L) ? std::clamp(-0.5L * q / r3, -1.0L, 1.0L) : 0.0L;
        ld phi = std::acos(arg) / 3.0L;
        std::array<ld, 3> ks = {
            polish(2.0L * r * std::cos(phi) + e1 / 3.0L),
            polish(2.0L * r * std::cos(phi + 2.0L * PI_L / 3.0L) + e1 / 3.0L),
            polish(2.0L * r * std::cos(phi + 4.0L * PI_L / 3.0L) + e1 / 3.0L)
        };
        std::sort(ks.begin(), ks.end(), [](ld a, ld b) { return a > b; });
        k1 = ks[0]; k2 = ks[1]; k3 = ks[2];
        margin = 2.0L - std::max({std::abs(k1), std::abs(k2), std::abs(k3)});
    } else {
        real_k = false;
        ld sd = std::sqrt(D);
        ld u = std::cbrt(-0.5L * q + sd);
        ld v = std::cbrt(-0.5L * q - sd);
        ld kr = polish(u + v + e1 / 3.0L);
        ld re = 0.5L * (e1 - kr);
        ld c0 = e2 + kr * (kr - e1);
        ld im = std::sqrt(std::max(0.0L, c0 - re * re));
        k1 = kr; k2 = re; k3 = im;
        margin = -10.0L;
    }
}

static ld max_nontriv_eig_mod(std::vector<cplx> eigs, int n_one, int n_rot_pairs, ld theta) {
    std::vector<cplx> targets;
    for (int i = 0; i < n_one; ++i) targets.emplace_back(1.0L, 0.0L);
    cplx ep(std::cos(theta), std::sin(theta));
    cplx em(std::cos(theta), -std::sin(theta));
    for (int i = 0; i < n_rot_pairs; ++i) {
        targets.push_back(ep);
        targets.push_back(em);
    }
    for (const auto& tgt : targets) {
        if (eigs.empty()) break;
        size_t best = 0;
        ld bd = std::abs(eigs[0] - tgt);
        for (size_t i = 1; i < eigs.size(); ++i) {
            ld d = std::abs(eigs[i] - tgt);
            if (d < bd) { bd = d; best = i; }
        }
        eigs.erase(eigs.begin() + best);
    }
    ld mx = 0.0L;
    for (const auto& ev : eigs) mx = std::max(mx, std::abs(ev));
    return mx;
}

// =========================================================================
// Exact Montgomery S3 shape-sphere free-group word (2D)
// =========================================================================
static std::array<ld, 3> shape_sphere_2d(const State3D& s) {
    ld x1 = s[0], y1 = s[1], x2 = s[3], y2 = s[4], x3 = s[6], y3 = s[7];
    ld rx = (x1 - x2) / std::sqrt(2.0L);
    ld ry = (y1 - y2) / std::sqrt(2.0L);
    ld lx = (x1 + x2 - 2.0L * x3) / std::sqrt(6.0L);
    ld ly = (y1 + y2 - 2.0L * y3) / std::sqrt(6.0L);
    ld r2 = rx * rx + ry * ry;
    ld l2 = lx * lx + ly * ly;
    ld R2 = r2 + l2;
    return {(r2 - l2) / R2, 2.0L * (rx * lx + ry * ly) / R2, 2.0L * (rx * ly - ry * lx) / R2};
}

static std::string reduce_word(const std::string& w) {
    auto inv = [](char c) -> char {
        if (c == 'a') return 'A';
        if (c == 'A') return 'a';
        if (c == 'b') return 'B';
        return 'b';
    };
    std::string st;
    for (char c : w) {
        if (!st.empty() && inv(st.back()) == c) st.pop_back();
        else st.push_back(c);
    }
    size_t l = 0, r = st.size();
    while (r - l >= 2 && inv(st[l]) == st[r - 1]) { ++l; --r; }
    return st.substr(l, r - l);
}

static std::string canonical_word_basic(const std::string& w) {
    if (w.empty()) return "1";
    auto inv = [](char c) { return c=='a'?'A':(c=='A'?'a':(c=='b'?'B':'b')); };
    auto swp = [](char c) { return c=='a'?'b':(c=='b'?'a':(c=='A'?'B':'A')); };
    std::string rev;
    for (auto it = w.rbegin(); it != w.rend(); ++it) rev.push_back(inv(*it));
    std::string best;
    for (const std::string& base : {w, rev}) {
        std::string t1 = base, t2, t3, t4;
        for (char c : base) { t2.push_back(swp(c)); t3.push_back(inv(c)); t4.push_back(swp(inv(c))); }
        for (const std::string& t : {t1, t2, t3, t4}) {
            for (size_t i = 0; i < t.size(); ++i) {
                std::string rot = t.substr(i) + t.substr(0, i);
                if (best.empty() || rot < best) best = rot;
            }
        }
    }
    return best;
}

static std::pair<std::string, int> primitive_root_word(const std::string& w) {
    if (w.empty() || w == "1") return {"1", 1};
    size_t n = w.size();
    for (size_t p = 1; p <= n; ++p) {
        if (n % p != 0) continue;
        bool ok = true;
        for (size_t i = p; i < n; ++i) {
            if (w[i] != w[i % p]) { ok = false; break; }
        }
        if (ok) return {canonical_word_basic(w.substr(0, p)), (int)(n / p)};
    }
    return {w, 1};
}

static std::pair<std::string, std::string> exact_word_2d(const State3D& y0, ld T, const Masses& m, int perm_idx = 0) {
    static constexpr int perms[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
    std::vector<std::array<ld, 3>> pts;
    State3D y = y0;
    ld t = 0, h = 1e-4L, rmin_track = 1e300L;
    pts.push_back(shape_sphere_2d(permute_state(y, perms[perm_idx])));
    while (t < T - 1e-16L) {
        ld hcap = std::min(h, T / 4000.0L);
        if (!dopri5_step_state(y, t, hcap, T, m, 1e-14L, rmin_track)) break;
        h = hcap;
        pts.push_back(shape_sphere_2d(permute_state(y, perms[perm_idx])));
    }
    if (pts.size() > 2) pts.pop_back();
    size_t start = 0;
    ld best_nz = -1;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::abs(pts[i][2]) > best_nz) { best_nz = std::abs(pts[i][2]); start = i; }
    }
    const ld c = 1.0L / std::sqrt(3.0L);
    std::string letters;
    size_t N = pts.size();
    int last_sign = (pts[start][2] > 0) ? 1 : -1;
    auto last_pt = pts[start];
    for (size_t k = 1; k <= N; ++k) {
        const auto& cur = pts[(start + k) % N];
        if (cur[2] == 0.0L) continue;
        int cur_sign = (cur[2] > 0) ? 1 : -1;
        if (cur_sign != last_sign) {
            ld frac = std::abs(last_pt[2]) / (std::abs(last_pt[2]) + std::abs(cur[2]));
            ld nx = last_pt[0] + frac * (cur[0] - last_pt[0]);
            ld ny = last_pt[1] + frac * (cur[1] - last_pt[1]);
            ld u = ny / (1.0L + nx);
            bool up = (cur_sign > 0);
            if (u > c) letters.push_back(up ? 'a' : 'A');
            else if (u < -c) letters.push_back(up ? 'b' : 'B');
            last_sign = cur_sign;
        }
        last_pt = cur;
    }
    std::string red = reduce_word(letters);
    return {red, canonical_word_basic(red)};
}

// =========================================================================
// Mode 2: Ultra-precision half-period symplectic involution verifier (2D L=0)
// =========================================================================
static inline void deriv2d(const State2D& s, const Masses& m, State2D& ds, ld* rm = nullptr) {
    for (int i = 0; i < 6; ++i) ds[i] = s[6 + i];
    ld a[6] = {0}, r2m = 1e300L;
    for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) {
        ld dx = s[2 * j] - s[2 * i], dy = s[2 * j + 1] - s[2 * i + 1];
        ld r2 = dx * dx + dy * dy; if (r2 < r2m) r2m = r2;
        ld ir = 1.0L / std::sqrt(r2), ir3 = ir * ir * ir;
        a[2 * i] += m[j] * dx * ir3; a[2 * i + 1] += m[j] * dy * ir3;
        a[2 * j] -= m[i] * dx * ir3; a[2 * j + 1] -= m[i] * dy * ir3;
    }
    for (int i = 0; i < 6; ++i) ds[6 + i] = a[i];
    if (rm) *rm = std::sqrt(r2m);
}

static inline void accel_jac2d(const State2D& s, const Masses& m, ld H[6][6]) {
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) H[i][j] = 0;
    for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) {
        ld d[2] = {s[2 * j] - s[2 * i], s[2 * j + 1] - s[2 * i + 1]};
        ld r2 = d[0] * d[0] + d[1] * d[1], ir = 1.0L / std::sqrt(r2), ir3 = ir * ir * ir, ir5 = ir3 * ir * ir;
        for (int p = 0; p < 2; ++p) for (int q = 0; q < 2; ++q) {
            ld Tpq = ((p == q) ? ir3 : 0.0L) - 3.0L * d[p] * d[q] * ir5;
            H[2 * i + p][2 * j + q] += m[j] * Tpq; H[2 * i + p][2 * i + q] -= m[j] * Tpq;
            H[2 * j + p][2 * i + q] += m[i] * Tpq; H[2 * j + p][2 * j + q] -= m[i] * Tpq;
        }
    }
}

struct Var2D { State2D y; ld M[12][12]; };

static inline void deriv_var2d(const Var2D& vs, const Masses& m, Var2D& dvs, ld* rm = nullptr) {
    deriv2d(vs.y, m, dvs.y, rm);
    ld H[6][6]; accel_jac2d(vs.y, m, H);
    for (int c = 0; c < 12; ++c) {
        for (int r = 0; r < 6; ++r) dvs.M[r][c] = vs.M[6 + r][c];
        for (int r = 0; r < 6; ++r) {
            ld s = 0; for (int k = 0; k < 6; ++k) s += H[r][k] * vs.M[k][c];
            dvs.M[6 + r][c] = s;
        }
    }
}

static bool integrate_var2d(const State2D& y0, ld T, const Masses& m, State2D& yT, ld Mout[12][12], ld tol = 1e-17L, ld mtol_scale = 1.0L) {
    static constexpr ld a21 = 1.0L / 5, a31 = 3.0L / 40, a32 = 9.0L / 40,
        a41 = 44.0L / 45, a42 = -56.0L / 15, a43 = 32.0L / 9,
        a51 = 19372.0L / 6561, a52 = -25360.0L / 2187, a53 = 64448.0L / 6561, a54 = -212.0L / 729,
        a61 = 9017.0L / 3168, a62 = -355.0L / 33, a63 = 46732.0L / 5247, a64 = 49.0L / 176, a65 = -5103.0L / 18656,
        b1 = 35.0L / 384, b3 = 500.0L / 1113, b4 = 125.0L / 192, b5 = -2187.0L / 6784, b6 = 11.0L / 84,
        e1 = 71.0L / 57600, e3 = -71.0L / 16695, e4 = 71.0L / 1920, e5 = -17253.0L / 339200, e6 = 22.0L / 525, e7 = -1.0L / 40;
    Var2D vs{}; vs.y = y0;
    for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) vs.M[i][j] = (i == j ? 1.0L : 0.0L);
    auto axpy = [](const Var2D& a, Var2D& out, ld h, int nk, const Var2D* ks[], const ld cs[]) {
        for (int i = 0; i < 12; ++i) { ld s = 0; for (int k = 0; k < nk; ++k) s += cs[k] * ks[k]->y[i]; out.y[i] = a.y[i] + h * s; }
        for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) { ld s = 0; for (int k = 0; k < nk; ++k) s += cs[k] * ks[k]->M[r][c]; out.M[r][c] = a.M[r][c] + h * s; }
    };
    Var2D k1, k2, k3, k4, k5, k6, k7, yt, yn;
    ld t = 0, h = std::min(1e-4L, T * 0.01L); int steps = 0;
    while (t < T - 1e-18L && steps < 4000000) {
        if (t + h > T) h = T - t;
        bool ok = false;
        for (int att = 0; att < 30; ++att) {
            ld rm = 1e300L; deriv_var2d(vs, m, k1, &rm); if (rm < 1e-6L) return false;
            {const Var2D* ks[1] = {&k1}; ld cs[1] = {a21}; axpy(vs, yt, h, 1, ks, cs);} deriv_var2d(yt, m, k2);
            {const Var2D* ks[2] = {&k1, &k2}; ld cs[2] = {a31, a32}; axpy(vs, yt, h, 2, ks, cs);} deriv_var2d(yt, m, k3);
            {const Var2D* ks[3] = {&k1, &k2, &k3}; ld cs[3] = {a41, a42, a43}; axpy(vs, yt, h, 3, ks, cs);} deriv_var2d(yt, m, k4);
            {const Var2D* ks[4] = {&k1, &k2, &k3, &k4}; ld cs[4] = {a51, a52, a53, a54}; axpy(vs, yt, h, 4, ks, cs);} deriv_var2d(yt, m, k5);
            {const Var2D* ks[5] = {&k1, &k2, &k3, &k4, &k5}; ld cs[5] = {a61, a62, a63, a64, a65}; axpy(vs, yt, h, 5, ks, cs);} deriv_var2d(yt, m, k6);
            {const Var2D* ks[5] = {&k1, &k3, &k4, &k5, &k6}; ld cs[5] = {b1, b3, b4, b5, b6}; axpy(vs, yn, h, 5, ks, cs);}
            deriv_var2d(yn, m, k7, &rm); if (rm < 1e-6L) return false;
            ld err = 0;
            for (int i = 0; i < 12; ++i) {
                ld ei = h * (e1 * k1.y[i] + e3 * k3.y[i] + e4 * k4.y[i] + e5 * k5.y[i] + e6 * k6.y[i] + e7 * k7.y[i]);
                ld sc = tol * (1.0L + std::max(std::abs(vs.y[i]), std::abs(yn.y[i])));
                ld r = ei / sc; err += r * r;
            }
            for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) {
                ld ei = h * (e1 * k1.M[r][c] + e3 * k3.M[r][c] + e4 * k4.M[r][c] + e5 * k5.M[r][c] + e6 * k6.M[r][c] + e7 * k7.M[r][c]);
                ld sc = (tol * mtol_scale) * (1.0L + std::max(std::abs(vs.M[r][c]), std::abs(yn.M[r][c])));
                ld rr = ei / sc; err += rr * rr;
            }
            err = std::sqrt(err / 156.0L);
            if (err <= 1.0L) {
                vs = yn; t += h;
                h = std::min(h * ((err == 0) ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.025L);
                ok = true; break;
            }
            h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
            if (h < 1e-17L) return false;
        }
        if (!ok) return false; ++steps;
    }
    yT = vs.y;
    for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) Mout[i][j] = vs.M[i][j];
    return t >= T - 1e-14L;
}

static void build_P(ld th, int mode, ld P[12][12]) {
    for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) P[i][j] = 0;
    ld c = std::cos(-th), s = std::sin(-th);
    if (mode == 0) {
        for (int b = 0; b < 6; ++b) {
            P[2 * b][2 * b] = c; P[2 * b][2 * b + 1] = -s;
            P[2 * b + 1][2 * b] = -s; P[2 * b + 1][2 * b + 1] = -c;
        }
    } else {
        int mapb[3] = {1, 0, 2};
        for (int off : {0, 6}) {
            for (int b = 0; b < 3; ++b) {
                int nb = mapb[b];
                P[off + 2 * nb][off + 2 * b] = -c; P[off + 2 * nb][off + 2 * b + 1] = s;
                P[off + 2 * nb + 1][off + 2 * b] = s; P[off + 2 * nb + 1][off + 2 * b + 1] = c;
            }
        }
    }
}

static State2D apply_P(const ld P[12][12], const State2D& y) {
    State2D r{}; for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) r[i] += P[i][j] * y[j]; return r;
}

static ld symp_err_12(const ld M[12][12], const Masses& m) {
    ld sc[12] = {1, 1, 1, 1, 1, 1, m[0], m[0], m[1], m[1], m[2], m[2]};
    ld symp = 0;
    for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) {
        ld s = 0;
        for (int k = 0; k < 6; ++k)
            s += (sc[k] * M[k][i] / sc[i]) * (sc[k + 6] * M[k + 6][j] / sc[j])
               - (sc[k + 6] * M[k + 6][i] / sc[i]) * (sc[k] * M[k][j] / sc[j]);
        symp = std::max(symp, std::abs(s - ((j == i + 6) ? 1.0L : ((i == j + 6) ? -1.0L : 0.0L))));
    }
    return symp;
}

static State2D make_col(ld vx, ld vy, ld m3) {
    return {-1, 0, 1, 0, 0, 0, vx, vy, vx, vy, -2 * vx / m3, -2 * vy / m3};
}

static int run_collinear_ultra(int nthreads) {
    struct ColJob { std::string id; ld m3, vx, vy, T; };
    std::vector<ColJob> jobs;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        ColJob jb;
        if (ss >> jb.id >> jb.m3 >> jb.vx >> jb.vy >> jb.T) jobs.push_back(jb);
    }
    std::vector<std::string> outs(jobs.size());
    std::atomic<size_t> next_idx{0};
    std::vector<std::thread> ths;
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= jobs.size()) break;
                auto jb = jobs[idx];
                ld m3 = jb.m3, vx = jb.vx, vy = jb.vy, T = jb.T;
                Masses m = {1.0L, 1.0L, m3};
                auto dist12 = [](const State2D& a, const State2D& b) {
                    ld s = 0; for (int i = 0; i < 12; ++i) s += (a[i] - b[i]) * (a[i] - b[i]); return std::sqrt(s);
                };

                State2D y0 = make_col(vx, vy, m3), yH;
                ld MH[12][12];
                bool okH = integrate_var2d(y0, T / 2.0L, m, yH, MH, 1e-15L, 10.0L);
                ld th = okH ? std::atan2(yH[3], yH[2]) : 0.0L;
                ld P0[12][12], P1[12][12];
                build_P(th, 0, P0); build_P(th, 1, P1);
                ld e0 = okH ? dist12(apply_P(P0, yH), y0) : 1e9L;
                ld e1 = okH ? dist12(apply_P(P1, yH), y0) : 1e9L;

                std::ostringstream os;
                os << std::setprecision(13);
                if (std::min(e0, e1) < 1e-3L) {
                    int mode = (e0 < e1) ? 0 : 1;
                    ld Th = T / 2.0L;
                    ld min_dk1 = 1e30L, max_dk1 = -1e30L, last_k1 = 0, last_k2 = 0, last_ret = 0, last_symp = 0;
                    for (int it = 0; it < 6; ++it) {
                        y0 = make_col(vx, vy, m3);
                        if (!integrate_var2d(y0, Th, m, yH, MH, 1e-17L, 1.0L)) break;
                        ld P[12][12], dP[12][12];
                        build_P(th, mode, P);
                        ld c = std::cos(-th), s = std::sin(-th);
                        for (int i = 0; i < 12; ++i) for (int k = 0; k < 12; ++k) dP[i][k] = 0;
                        if (mode == 0) {
                            for (int b = 0; b < 6; ++b) {
                                dP[2 * b][2 * b] = s; dP[2 * b][2 * b + 1] = c;
                                dP[2 * b + 1][2 * b] = c; dP[2 * b + 1][2 * b + 1] = -s;
                            }
                        } else {
                            int mapb[3] = {1, 0, 2};
                            for (int off : {0, 6}) for (int b = 0; b < 3; ++b) {
                                int nb = mapb[b];
                                dP[off + 2 * nb][off + 2 * b] = -s; dP[off + 2 * nb][off + 2 * b + 1] = -c;
                                dP[off + 2 * nb + 1][off + 2 * b] = -c; dP[off + 2 * nb + 1][off + 2 * b + 1] = s;
                            }
                        }
                        State2D PyH = apply_P(P, yH);
                        State2D F; ld fn = 0;
                        for (int i = 0; i < 12; ++i) { F[i] = PyH[i] - y0[i]; fn += F[i] * F[i]; }
                        ld retH = std::sqrt(fn);
                        ld A[12][12] = {{0}};
                        for (int i = 0; i < 12; ++i) for (int c_ = 0; c_ < 12; ++c_)
                            for (int k = 0; k < 12; ++k) A[i][c_] += P[i][k] * MH[k][c_];
                        ld tr1 = 0, tr2 = 0;
                        for (int i = 0; i < 12; ++i) {
                            tr1 += A[i][i];
                            for (int c_ = 0; c_ < 12; ++c_) tr2 += A[i][c_] * A[c_][i];
                        }
                        ld sA = tr1, sq = tr2 - 4.0L, prod = 0.5L * (sA * sA - sq), disc = sA * sA - 4.0L * prod;
                        ld h1 = 0.5L * (sA + std::sqrt(std::max(0.0L, disc)));
                        ld h2 = 0.5L * (sA - std::sqrt(std::max(0.0L, disc)));
                        ld k1 = h1 * h1 - 2.0L, k2 = h2 * h2 - 2.0L;
                        if (std::abs(k2 - 2.0L) < std::abs(k1 - 2.0L)) { std::swap(k1, k2); std::swap(h1, h2); }
                        if (it >= 1) {
                            min_dk1 = std::min(min_dk1, k1 - 2.0L);
                            max_dk1 = std::max(max_dk1, k1 - 2.0L);
                        }
                        last_k1 = k1; last_k2 = k2; last_ret = retH; last_symp = symp_err_12(A, m);
                        if (it == 5 || (it >= 2 && retH < 1e-16L)) break;
                        ld dvx[12] = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3, 0};
                        ld dvy[12] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3};
                        State2D fH; deriv2d(yH, m, fH);
                        State2D PfH = apply_P(P, fH), dPyH = apply_P(dP, yH);
                        ld J[12][4];
                        for (int r = 0; r < 12; ++r) {
                            ld s0 = -dvx[r], s1 = -dvy[r];
                            for (int k = 0; k < 12; ++k) { s0 += A[r][k] * dvx[k]; s1 += A[r][k] * dvy[k]; }
                            J[r][0] = s0; J[r][1] = s1; J[r][2] = PfH[r]; J[r][3] = dPyH[r];
                        }
                        ld G[4][5] = {{0}};
                        for (int r = 0; r < 4; ++r) {
                            for (int c_ = 0; c_ < 4; ++c_) for (int i = 0; i < 12; ++i) G[r][c_] += J[i][r] * J[i][c_];
                            for (int i = 0; i < 12; ++i) G[r][4] -= J[i][r] * F[i];
                        }
                        for (int c_ = 0; c_ < 4; ++c_) {
                            int piv = c_;
                            for (int r = c_ + 1; r < 4; ++r) if (std::abs(G[r][c_]) > std::abs(G[piv][c_])) piv = r;
                            for (int k = c_; k < 5; ++k) std::swap(G[c_][k], G[piv][k]);
                            for (int r = c_ + 1; r < 4; ++r) {
                                ld f = G[r][c_] / G[c_][c_];
                                for (int k = c_; k < 5; ++k) G[r][k] -= f * G[c_][k];
                            }
                        }
                        ld d[4] = {0};
                        for (int r = 3; r >= 0; --r) {
                            ld s_ = G[r][4];
                            for (int c_ = r + 1; c_ < 4; ++c_) s_ -= G[r][c_] * d[c_];
                            d[r] = s_ / G[r][r];
                        }
                        vx += d[0]; vy += d[1]; Th += d[2]; th += d[3];
                    }
                    ld mid_dk1 = 0.5L * (min_dk1 + max_dk1);
                    ld half_span = 0.5L * (max_dk1 - min_dk1);
                    ld k1_fin = 2.0L + mid_dk1;
                    ld margin = 2.0L - std::max(std::abs(k1_fin), std::abs(last_k2));
                    std::string verdict;
                    if (min_dk1 > 0.0L || std::abs(last_k2) >= 2.0L) verdict = "UNSTABLE";
                    else if (max_dk1 >= 0.0L) verdict = "NEAR-PARABOLIC";
                    else if (margin < 1e-4L) verdict = "SUB<1e-4";
                    else verdict = "YES";
                    os << "{\"id\":\"" << jb.id << "\",\"mode\":\"half_" << mode << "\""
                       << ",\"m3\":" << (double)m3 << ",\"vx\":" << (double)vx
                       << ",\"vy\":" << (double)vy << ",\"T\":" << (double)(2.0L * Th)
                       << ",\"ret_err\":" << (double)last_ret << ",\"symp_err\":" << (double)last_symp
                       << ",\"k1\":" << (double)k1_fin << ",\"k2\":" << (double)last_k2
                       << ",\"dk1_mid\":" << (double)mid_dk1 << ",\"dk1_half_span\":" << (double)half_span
                       << ",\"stab_margin\":" << (double)margin << ",\"verdict\":\"" << verdict << "\"}";
                } else {
                    State2D yT; ld M[12][12];
                    for (int it = 0; it < 6; ++it) {
                        y0 = make_col(vx, vy, m3);
                        if (!integrate_var2d(y0, T, m, yT, M, 1e-16L, 10.0L)) break;
                        State2D F; ld fn = 0;
                        for (int i = 0; i < 12; ++i) { F[i] = yT[i] - y0[i]; fn += F[i] * F[i]; }
                        if (std::sqrt(fn) < 1e-15L) break;
                        ld dvx[12] = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3, 0};
                        ld dvy[12] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3};
                        State2D fT; deriv2d(yT, m, fT);
                        ld J[12][3];
                        for (int r = 0; r < 12; ++r) {
                            ld s0 = -dvx[r], s1 = -dvy[r];
                            for (int k = 0; k < 12; ++k) { s0 += M[r][k] * dvx[k]; s1 += M[r][k] * dvy[k]; }
                            J[r][0] = s0; J[r][1] = s1; J[r][2] = fT[r];
                        }
                        ld A[3][4] = {{0}};
                        for (int r = 0; r < 3; ++r) {
                            for (int c_ = 0; c_ < 3; ++c_) for (int i = 0; i < 12; ++i) A[r][c_] += J[i][r] * J[i][c_];
                            for (int i = 0; i < 12; ++i) A[r][3] -= J[i][r] * F[i];
                            A[r][r] += 1e-22L;
                        }
                        for (int c_ = 0; c_ < 3; ++c_) {
                            int piv = c_;
                            for (int r = c_ + 1; r < 3; ++r) if (std::abs(A[r][c_]) > std::abs(A[piv][c_])) piv = r;
                            for (int k = c_; k < 4; ++k) std::swap(A[c_][k], A[piv][k]);
                            for (int r = c_ + 1; r < 3; ++r) {
                                ld f = A[r][c_] / A[c_][c_];
                                for (int k = c_; k < 4; ++k) A[r][k] -= f * A[c_][k];
                            }
                        }
                        ld d[3] = {0};
                        for (int r = 2; r >= 0; --r) {
                            ld s = A[r][3];
                            for (int c_ = r + 1; c_ < 3; ++c_) s -= A[r][c_] * d[c_];
                            d[r] = s / A[r][r];
                        }
                        ld sn = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                        if (sn > 1e-5L) for (int i = 0; i < 3; ++i) d[i] *= (1e-5L / sn);
                        vx += d[0]; vy += d[1]; T += d[2];
                    }
                    y0 = make_col(vx, vy, m3);
                    integrate_var2d(y0, T, m, yT, M, 1e-16L, 10.0L);
                    ld ret = dist12(yT, y0), symp = symp_err_12(M, m);
                    ld tr1 = 0, tr2 = 0;
                    for (int i = 0; i < 12; ++i) tr1 += M[i][i];
                    for (int i = 0; i < 12; ++i) for (int c_ = 0; c_ < 12; ++c_) tr2 += M[i][c_] * M[c_][i];
                    ld s = tr1 - 8.0L, p = 0.5L * (s * s - (tr2 - 4.0L)), disc = s * s - 4.0L * p;
                    bool real_k = (disc >= -1e-10L);
                    ld k1 = real_k ? 0.5L * (s + std::sqrt(std::max(0.0L, disc))) : 0.5L * s;
                    ld k2 = real_k ? 0.5L * (s - std::sqrt(std::max(0.0L, disc))) : 0.5L * std::sqrt(-disc);
                    bool st = real_k && (std::abs(k1) < 2.0L) && (std::abs(k2) < 2.0L);
                    ld marg = real_k ? (2.0L - std::max(std::abs(k1), std::abs(k2))) : -10.0L;
                    std::string verdict = st ? (marg > 1e-4L ? "YES" : "SUB<1e-4") : "UNSTABLE";
                    os << "{\"id\":\"" << jb.id << "\",\"mode\":\"full\""
                       << ",\"m3\":" << (double)m3 << ",\"vx\":" << (double)vx
                       << ",\"vy\":" << (double)vy << ",\"T\":" << (double)T
                       << ",\"ret_err\":" << (double)ret << ",\"symp_err\":" << (double)symp
                       << ",\"k1\":" << (double)k1 << ",\"k2\":" << (double)k2
                       << ",\"dk1_mid\":" << (double)(k1 - 2.0L) << ",\"dk1_half_span\":0"
                       << ",\"stab_margin\":" << (double)marg << ",\"verdict\":\"" << verdict << "\"}";
                }
                outs[idx] = os.str();
            }
        });
    }
    for (auto& th : ths) th.join();
    for (const auto& s : outs) std::cout << s << "\n";
    return 0;
}

// =========================================================================
// Mode 1: Full 3D/2D 80-bit Variational, RPO/Cover, Geometry & S3 Word Audit
// =========================================================================
struct AuditJob {
    std::string id, subspace;
    int dim;
    Masses m;
    ld T, precess_rad;
    State3D y0;
};

static std::string audit_one(const AuditJob& jb) {
    const Masses& m = jb.m;
    const State3D& y0 = jb.y0;
    ld T = jb.T, precess_rad = jb.precess_rad;
    int dim = jb.dim;

    ld Mtot = m[0] + m[1] + m[2];
    ld Rcm[3] = {0}, Vcm[3] = {0};
    for (int b = 0; b < 3; ++b)
        for (int d = 0; d < 3; ++d) {
            Rcm[d] += m[b] * y0[3 * b + d] / Mtot;
            Vcm[d] += m[b] * y0[9 + 3 * b + d] / Mtot;
        }
    ld rcm_err = std::sqrt(Rcm[0]*Rcm[0] + Rcm[1]*Rcm[1] + Rcm[2]*Rcm[2]);
    ld vcm_err = std::sqrt(Vcm[0]*Vcm[0] + Vcm[1]*Vcm[1] + Vcm[2]*Vcm[2]);
    ld E0 = energy3d(y0, m);
    auto L0 = angmom3d(y0, m);
    ld Lmag = std::sqrt(L0[0]*L0[0] + L0[1]*L0[1] + L0[2]*L0[2]);
    ld u[3] = {0, 0, 1};
    if (Lmag > 1e-10L) { u[0] = L0[0]/Lmag; u[1] = L0[1]/Lmag; u[2] = L0[2]/Lmag; }

    State3D yT{};
    ld Mraw[18][18];
    bool ok_var = integrate_var(y0, T, m, yT, Mraw, 1e-14L);
    ld ET = energy3d(yT, m);
    auto LT = angmom3d(yT, m);
    ld dE = std::abs(ET - E0);
    ld dL = std::sqrt((LT[0]-L0[0])*(LT[0]-L0[0]) + (LT[1]-L0[1])*(LT[1]-L0[1]) + (LT[2]-L0[2])*(LT[2]-L0[2]));

    bool is_rpo = (std::abs(precess_rad) > 1e-12L);
    State3D yT_chk = is_rpo ? rotate_state(yT, u, -precess_rad) : yT;
    ld ret_err = state_dist(yT_chk, y0);

    auto build_M_and_symp = [&](const ld raw[18][18], std::vector<std::vector<ld>>& M) -> ld {
        M.assign(18, std::vector<ld>(18, 0.0L));
        if (is_rpo) {
            for (int col = 0; col < 18; ++col) {
                State3D v{};
                for (int r = 0; r < 18; ++r) v[r] = raw[r][col];
                v = rotate_state(v, u, -precess_rad);
                for (int r = 0; r < 18; ++r) M[r][col] = v[r];
            }
        } else {
            for (int r = 0; r < 18; ++r) for (int c = 0; c < 18; ++c) M[r][c] = raw[r][c];
        }
        ld scale[18];
        for (int b = 0; b < 3; ++b) for (int d = 0; d < 3; ++d) { scale[3*b+d] = 1.0L; scale[9+3*b+d] = m[b]; }
        ld symp_err = 0;
        for (int i = 0; i < 18; ++i) {
            for (int j = 0; j < 18; ++j) {
                ld s = 0;
                for (int k = 0; k < 9; ++k) {
                    ld mki = scale[k] * M[k][i] / scale[i];
                    ld mkj = scale[k] * M[k][j] / scale[j];
                    ld m9i = scale[k + 9] * M[k + 9][i] / scale[i];
                    ld m9j = scale[k + 9] * M[k + 9][j] / scale[j];
                    s += mki * m9j - m9i * mkj;
                }
                ld expected = (j == i + 9) ? 1.0L : ((i == j + 9) ? -1.0L : 0.0L);
                symp_err = std::max(symp_err, std::abs(s - expected));
            }
        }
        return symp_err;
    };

    std::vector<std::vector<ld>> M;
    ld symp_err = build_M_and_symp(Mraw, M);
    // Automatic high-accuracy retry if symp_err > 1e-7
    if (ok_var && symp_err > 1e-7L) {
        State3D yT2{};
        ld Mraw2[18][18];
        if (integrate_var(y0, T, m, yT2, Mraw2, 1e-16L)) {
            std::vector<std::vector<ld>> M2;
            ld s2 = build_M_and_symp(Mraw2, M2);
            if (s2 < symp_err) {
                yT = yT2; M = M2; symp_err = s2;
                yT_chk = is_rpo ? rotate_state(yT, u, -precess_rad) : yT;
                ret_err = state_dist(yT_chk, y0);
                ET = energy3d(yT, m); LT = angmom3d(yT, m);
                dE = std::abs(ET - E0);
                dL = std::sqrt((LT[0]-L0[0])*(LT[0]-L0[0]) + (LT[1]-L0[1])*(LT[1]-L0[1]) + (LT[2]-L0[2])*(LT[2]-L0[2]));
            }
        }
    }

    // Eigenvalues + exact Broucke indices
    bool real_k = false;
    ld k1 = 0, k2 = 0, k3 = 0, kz = 2.0L, stab_margin = -10.0L;
    ld max_eig_planar = 0, max_eig_vert = 0, max_eig_3d = 0, max_eig_nontriv = 0;
    bool stable_3d = false;
    ld theta_rpo = is_rpo ? precess_rad : 0.0L;
    if (dim == 2) {
        const int pidx[12] = {0,1,3,4,6,7, 9,10,12,13,15,16};
        const int zidx[6] = {2,5,8, 11,14,17};
        std::vector<std::vector<ld>> Mp(12, std::vector<ld>(12)), Mz(6, std::vector<ld>(6));
        for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) Mp[r][c] = M[pidx[r]][pidx[c]];
        for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) Mz[r][c] = M[zidx[r]][zidx[c]];
        broucke_2d_from_matrix(Mp, theta_rpo, real_k, k1, k2, stab_margin);
        ld trz = 0;
        for (int i = 0; i < 6; ++i) trz += Mz[i][i];
        kz = trz - (2.0L + 2.0L * std::cos(theta_rpo));
        auto ep = eigenvalues_qr(Mp);
        auto ez = eigenvalues_qr(Mz);
        for (auto& ev : ep) max_eig_planar = std::max(max_eig_planar, std::abs(ev));
        for (auto& ev : ez) max_eig_vert = std::max(max_eig_vert, std::abs(ev));
        max_eig_3d = std::max(max_eig_planar, max_eig_vert);
        max_eig_nontriv = max_nontriv_eig_mod(ep, 4, 2, theta_rpo);
        stable_3d = real_k && (stab_margin > 1e-4L) && (Lmag > 1e-7L) && (std::abs(kz) < 2.0L - 1e-4L);
    } else {
        broucke_3d_from_matrix(M, theta_rpo, real_k, k1, k2, k3, stab_margin);
        auto e3 = eigenvalues_qr(M);
        for (auto& ev : e3) max_eig_3d = std::max(max_eig_3d, std::abs(ev));
        max_eig_nontriv = max_nontriv_eig_mod(e3, 6, 3, theta_rpo);
        stable_3d = real_k && (stab_margin > 1e-4L);
    }

    // Multiple-cover & RPO divisor scan q = 2..36
    std::vector<std::array<int, 3>> valid_perms;
    int all_p[6][3] = {{0,1,2},{1,0,2},{0,2,1},{2,1,0},{1,2,0},{2,0,1}};
    for (auto& p : all_p) {
        if (std::abs(m[p[0]] - m[0]) < 1e-9L && std::abs(m[p[1]] - m[1]) < 1e-9L && std::abs(m[p[2]] - m[2]) < 1e-9L) {
            valid_perms.push_back({p[0], p[1], p[2]});
        }
    }
    int best_q_inert = 1, best_q_perm = 1, best_q_rpo = 1, best_q_rpo_perm = 1, best_q_rpo_refl = 1;
    ld best_rpo_angle = 0;
    for (int q = 2; q <= 36; ++q) {
        State3D yq = y0;
        if (!integrate_state(yq, T / q, m, 1e-14L)) continue;
        if (state_dist(yq, y0) < 1e-4L) best_q_inert = std::max(best_q_inert, q);
        for (auto& pm : valid_perms) {
            State3D yp = permute_state(yq, pm.data());
            if (state_dist(yp, y0) < 1e-4L) best_q_perm = std::max(best_q_perm, q);
            ld md = 1e9L;
            ld th = optimal_rot_angle(yp, y0, u, md);
            if (pm[0]==0 && pm[1]==1 && pm[2]==2) {
                if (md < 1e-4L && q > best_q_rpo) { best_q_rpo = q; best_rpo_angle = th; }
            }
            if (md < 1e-4L && q > best_q_rpo_perm) { best_q_rpo_perm = q; }
            if (dim == 3) {
                State3D ypr = reflect_plane(yp, u);
                ld mdr = 1e9L;
                optimal_rot_angle(ypr, y0, u, mdr);
                if (mdr < 1e-4L && q > best_q_rpo_refl) best_q_rpo_refl = q;
            }
        }
    }

    // Dense 2400-step adaptive geometry, swaps, syzygies, and 3D thickness
    const int NS = 2400;
    ld dt = T / NS;
    State3D yc = y0;
    ld rmin_pair[3] = {1e300L, 1e300L, 1e300L};
    ld rmax_pair[3] = {0, 0, 0};
    ld Rmin_body[3] = {1e300L, 1e300L, 1e300L};
    ld Rmax_body[3] = {0, 0, 0};
    int closest_counts[3] = {0, 0, 0};
    int prev_cl = -1, exact_swaps = 0, syzygies = 0;
    ld prev_area = 0;
    ld C[3][3] = {{0}};
    int cov_cnt = 0;

    for (int k = 0; k <= NS; ++k) {
        auto dist3 = [&](int a, int b) {
            ld dx = yc[3*b]-yc[3*a], dy = yc[3*b+1]-yc[3*a+1], dz = yc[3*b+2]-yc[3*a+2];
            return std::sqrt(dx*dx + dy*dy + dz*dz);
        };
        ld r[3] = {dist3(0, 1), dist3(0, 2), dist3(1, 2)};
        for (int p = 0; p < 3; ++p) {
            rmin_pair[p] = std::min(rmin_pair[p], r[p]);
            rmax_pair[p] = std::max(rmax_pair[p], r[p]);
        }
        for (int b = 0; b < 3; ++b) {
            ld qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
            ld Rb = std::sqrt(qv[0]*qv[0] + qv[1]*qv[1] + qv[2]*qv[2]);
            Rmin_body[b] = std::min(Rmin_body[b], Rb);
            Rmax_body[b] = std::max(Rmax_body[b], Rb);
            if (dim == 3 && k < NS) {
                for (int r_ = 0; r_ < 3; ++r_) for (int c_ = 0; c_ < 3; ++c_) C[r_][c_] += qv[r_] * qv[c_];
                cov_cnt++;
            }
        }
        int cl = (r[0] <= r[1] && r[0] <= r[2]) ? 0 : ((r[1] <= r[2]) ? 1 : 2);
        if (k < NS) closest_counts[cl]++;
        if (prev_cl != -1 && cl != prev_cl) exact_swaps++;
        prev_cl = cl;

        if (dim == 2) {
            ld area = (yc[3] - yc[0]) * (yc[7] - yc[1]) - (yc[4] - yc[1]) * (yc[6] - yc[0]);
            if (k > 0 && ((prev_area > 0 && area < 0) || (prev_area < 0 && area > 0))) syzygies++;
            if (area != 0.0L) prev_area = area;
        }
        if (k < NS) integrate_state(yc, dt, m, 1e-13L);
    }

    ld thickness_3d = 0.0L;
    if (dim == 3 && cov_cnt > 0) {
        for (int r_ = 0; r_ < 3; ++r_) for (int c_ = 0; c_ < 3; ++c_) C[r_][c_] /= cov_cnt;
        thickness_3d = std::sqrt(min_eig_sym_3x3(C));
    }

    ld rmin_all = std::min({rmin_pair[0], rmin_pair[1], rmin_pair[2]});
    ld rmax_all = std::max({rmax_pair[0], rmax_pair[1], rmax_pair[2]});

    // Exact S3 shape-sphere word for 2D
    std::string w_red = "-", w_can = "-", w_can_s3 = "-", w_prim = "-";
    int w_pow = 1;
    if (dim == 2) {
        auto [wr0, wc0] = exact_word_2d(y0, T, m, 0);
        auto [wr1, wc1] = exact_word_2d(y0, T, m, 1);
        auto [wr2, wc2] = exact_word_2d(y0, T, m, 2);
        w_red = wr0.empty() ? "1" : wr0;
        w_can = wc0;
        w_can_s3 = wc0;
        auto better = [](const std::string& a, const std::string& b) {
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        };
        if (better(wc1, w_can_s3)) w_can_s3 = wc1;
        if (better(wc2, w_can_s3)) w_can_s3 = wc2;
        auto [wp, wpw] = primitive_root_word(w_can_s3);
        w_prim = wp;
        w_pow = wpw;
    }

    ld Tstar = T * std::pow(std::abs(E0), 1.5L);
    ld Lambda = Lmag * std::sqrt(std::max(0.0L, -E0));
    ld L_out = (dim == 2) ? L0[2] : Lmag;

    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"id\":\"" << jb.id << "\""
       << ",\"dim\":" << dim
       << ",\"subspace\":\"" << jb.subspace << "\""
       << ",\"ok_var\":" << (ok_var ? "true" : "false")
       << ",\"masses\":[" << (double)m[0] << "," << (double)m[1] << "," << (double)m[2] << "]"
       << ",\"T\":" << (double)T
       << ",\"E\":" << (double)E0
       << ",\"L\":" << (double)L_out
       << ",\"Lmag\":" << (double)Lmag
       << ",\"Lvec\":[" << (double)L0[0] << "," << (double)L0[1] << "," << (double)L0[2] << "]"
       << ",\"Tstar\":" << (double)Tstar
       << ",\"Lambda\":" << (double)Lambda
       << ",\"rcm_err\":" << (double)rcm_err
       << ",\"vcm_err\":" << (double)vcm_err
       << ",\"ret_err_ld\":" << (double)ret_err
       << ",\"dE\":" << (double)dE
       << ",\"dL\":" << (double)dL
       << ",\"symp_err\":" << (double)symp_err
       << ",\"real_k\":" << (real_k ? "true" : "false")
       << ",\"k1\":" << (double)k1
       << ",\"k2\":" << (double)k2;
    if (dim == 3) {
        os << ",\"k3\":" << (double)k3;
    } else {
        os << ",\"kz\":" << (double)kz;
    }
    os << ",\"stab_margin\":" << (double)stab_margin
       << ",\"stable_3d\":" << (stable_3d ? "true" : "false")
       << ",\"max_eig_planar\":" << (double)max_eig_planar
       << ",\"max_eig_vert\":" << (double)max_eig_vert
       << ",\"max_eig_3d\":" << (double)max_eig_3d
       << ",\"max_eig_nontriv\":" << (double)max_eig_nontriv
       << ",\"q_inert\":" << best_q_inert
       << ",\"q_perm\":" << best_q_perm
       << ",\"q_rpo\":" << best_q_rpo
       << ",\"q_rpo_perm\":" << best_q_rpo_perm
       << ",\"q_rpo_refl\":" << best_q_rpo_refl
       << ",\"rpo_angle\":" << (double)best_rpo_angle
       << ",\"rmin\":" << (double)rmin_all
       << ",\"rmax\":" << (double)rmax_all
       << ",\"rmin_pair\":[" << (double)rmin_pair[0] << "," << (double)rmin_pair[1] << "," << (double)rmin_pair[2] << "]"
       << ",\"rmax_pair\":[" << (double)rmax_pair[0] << "," << (double)rmax_pair[1] << "," << (double)rmax_pair[2] << "]"
       << ",\"Rmin_body\":[" << (double)Rmin_body[0] << "," << (double)Rmin_body[1] << "," << (double)Rmin_body[2] << "]"
       << ",\"Rmax_body\":[" << (double)Rmax_body[0] << "," << (double)Rmax_body[1] << "," << (double)Rmax_body[2] << "]"
       << ",\"closest_counts\":[" << closest_counts[0] << "," << closest_counts[1] << "," << closest_counts[2] << "]"
       << ",\"exact_swaps\":" << exact_swaps
       << ",\"syzygies\":" << syzygies
       << ",\"thickness_3d\":" << (double)thickness_3d
       << ",\"word_red\":\"" << w_red << "\""
       << ",\"word_can\":\"" << w_can << "\""
       << ",\"word_can_s3\":\"" << w_can_s3 << "\""
       << ",\"word_prim\":\"" << w_prim << "\""
       << ",\"word_pow\":" << w_pow
       << "}";
    return os.str();
}

int main(int argc, char** argv) {
    bool collinear_ultra = false;
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--collinear-ultra") collinear_ultra = true;
        else if (arg == "--threads" && i + 1 < argc) nthreads = std::max(1, std::atoi(argv[++i]));
    }
    if (collinear_ultra) return run_collinear_ultra(nthreads);

    std::vector<AuditJob> jobs;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        AuditJob jb{};
        if (!(ss >> jb.id >> jb.dim >> jb.subspace >> jb.m[0] >> jb.m[1] >> jb.m[2] >> jb.T >> jb.precess_rad)) continue;
        for (int i = 0; i < 18; ++i) ss >> jb.y0[i];
        jobs.push_back(jb);
    }

    std::vector<std::string> outs(jobs.size());
    std::atomic<size_t> next_idx{0};
    std::vector<std::thread> ths;
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= jobs.size()) break;
                outs[idx] = audit_one(jobs[idx]);
            }
        });
    }
    for (auto& th : ths) th.join();
    for (const auto& s : outs) std::cout << s << "\n";
    return 0;
}
