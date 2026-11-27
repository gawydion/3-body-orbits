#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using ld = long double;
using cplx = std::complex<ld>;

static constexpr ld PI_L = 3.141592653589793238462643383279502884L;

// =========================================================================
// General real matrix eigenvalues via Hessenberg + Francis double-shift QR
// =========================================================================
static std::vector<cplx> eigenvalues_qr(std::vector<std::vector<ld>> A) {
    int n = (int)A.size();
    // 1. Reduce to upper Hessenberg via Householder reflections
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
        // Left multiply: A[k+1:n, k:n] -= beta * v * (v^T * A[k+1:n, k:n])
        for (int j = k; j < n; ++j) {
            ld dot = 0;
            for (int i = 0; i < n - k - 1; ++i) dot += v[i] * A[k + 1 + i][j];
            for (int i = 0; i < n - k - 1; ++i) A[k + 1 + i][j] -= beta * dot * v[i];
        }
        // Right multiply: A[0:n, k+1:n] -= beta * (A[0:n, k+1:n] * v) * v^T
        for (int i = 0; i < n; ++i) {
            ld dot = 0;
            for (int j = 0; j < n - k - 1; ++j) dot += A[i][k + 1 + j] * v[j];
            for (int j = 0; j < n - k - 1; ++j) A[i][k + 1 + j] -= beta * dot * v[j];
        }
        for (int i = k + 2; i < n; ++i) A[i][k] = 0;
    }

    // 2. Francis double-shift QR iterations
    std::vector<cplx> eigs;
    int m = n;
    int iter = 0;
    const ld eps = 1e-17L;
    while (m > 0) {
        if (m == 1) {
            eigs.emplace_back(A[0][0], 0);
            break;
        }
        // Check deflation
        int l = m - 1;
        while (l > 0) {
            ld s = std::abs(A[l - 1][l - 1]) + std::abs(A[l][l]);
            if (s == 0) s = 1.0L;
            if (std::abs(A[l][l - 1]) <= eps * s) {
                A[l][l - 1] = 0;
                break;
            }
            --l;
        }
        if (l == m - 1) {
            eigs.emplace_back(A[m - 1][m - 1], 0);
            --m;
            iter = 0;
            continue;
        }
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
            m -= 2;
            iter = 0;
            continue;
        }
        ++iter;
        if (iter > 2000) {
            // Fallback: extract diagonal
            for (int i = 0; i < m; ++i) eigs.emplace_back(A[i][i], 0);
            break;
        }
        // Compute shifts from bottom 2x2 (or exceptional shift)
        ld s, t;
        if (iter % 15 == 0) {
            ld w = std::abs(A[m - 1][m - 2]) + std::abs(A[m - 2][m - 3]);
            s = 1.5L * w;
            t = w * w;
        } else {
            s = A[m - 2][m - 2] + A[m - 1][m - 1];
            t = A[m - 2][m - 2] * A[m - 1][m - 1] - A[m - 2][m - 1] * A[m - 1][m - 2];
        }
        // First column of (H - s1 I)(H - s2 I)
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
                x = A[k + 1][k];
                y = A[k + 2][k];
                z = (k + 3 < m) ? A[k + 3][k] : 0;
            }
        }
    }
    return eigs;
}

// =========================================================================
// 3D Gravitational N=3 Dynamics + Analytical Variational Equations (18 + 324)
// (2D orbits are embedded in 3D with z=0, vz=0, giving BOTH planar AND
//  vertical stability simultaneously!)
// =========================================================================
using Masses = std::array<ld, 3>;
using State3D = std::array<ld, 18>;

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

// Analytical Jacobian A(s) of f(s): ds_dot = A(s) * delta_s
// Upper right 9x9 = I, Lower left 9x9 = H (3x3 blocks of 3x3)
static inline void compute_accel_jacobian(const State3D& s, const Masses& m, ld H[9][9]) {
    for (int i = 0; i < 9; ++i) for (int j = 0; j < 9; ++j) H[i][j] = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            ld d[3] = {s[3 * j] - s[3 * i], s[3 * j + 1] - s[3 * i + 1], s[3 * j + 2] - s[3 * i + 2]};
            ld r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            ld ir = 1.0L / std::sqrt(r2);
            ld ir3 = ir * ir * ir;
            ld ir5 = ir3 * ir * ir;
            // d(a_i)/d(r_j) = m_j * (I / r^3 - 3 d d^T / r^5)
            // d(a_i)/d(r_i) = -d(a_i)/d(r_j)
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

// Combined state + variational matrix (18 + 18*18 = 342)
struct VarState {
    State3D y;
    ld M[18][18];
};

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

static bool integrate_state(State3D& y, ld T, const Masses& m, ld tol = 1e-16L, ld* rmin_out = nullptr) {
    ld t = 0, h = std::min(1e-3L, T * 0.01L), rmin_track = 1e300L;
    int steps = 0;
    while (t < T - 1e-17L && steps < 2000000) {
        if (!dopri5_step_state(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-14L;
}

static bool integrate_var(const State3D& y0, ld T, const Masses& m, State3D& yT, ld Mout[18][18], ld tol = 1e-15L) {
    static constexpr ld a21 = 1.0L / 5.0L, a31 = 3.0L / 40.0L, a32 = 9.0L / 40.0L,
        a41 = 44.0L / 45.0L, a42 = -56.0L / 15.0L, a43 = 32.0L / 9.0L,
        a51 = 19372.0L / 6561.0L, a52 = -25360.0L / 2187.0L, a53 = 64448.0L / 6561.0L, a54 = -212.0L / 729.0L,
        a61 = 9017.0L / 3168.0L, a62 = -355.0L / 33.0L, a63 = 46732.0L / 5247.0L, a64 = 49.0L / 176.0L, a65 = -5103.0L / 18656.0L,
        b1 = 35.0L / 384.0L, b3 = 500.0L / 1113.0L, b4 = 125.0L / 192.0L, b5 = -2187.0L / 6784.0L, b6 = 11.0L / 84.0L,
        e1 = 71.0L / 57600.0L, e3 = -71.0L / 16695.0L, e4 = 71.0L / 1920.0L, e5 = -17253.0L / 339200.0L, e6 = 22.0L / 525.0L, e7 = -1.0L / 40.0L;

    VarState vs{};
    vs.y = y0;
    for (int i = 0; i < 18; ++i) for (int j = 0; j < 18; ++j) vs.M[i][j] = (i == j ? 1.0L : 0.0L);

    auto axpy = [](const VarState& a, VarState& out, ld h,
                   int nk, const VarState* ks[], const ld cs[]) {
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
    while (t < T - 1e-17L && steps < 2000000) {
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

// Optimal rotation angle around u minimizing ||R_u(-theta) s - target||^2
static ld optimal_rot_angle(const State3D& s, const State3D& target, const ld u[3], ld& min_dist) {
    // R_u(-theta) v = v cos(theta) - (u x v) sin(theta) + u (u.v) (1 - cos(theta))
    // Maximize sum_k target_k . R_u(-theta) s_k
    // = const + A cos(theta) + B sin(theta)
    // where A = sum (t.s - (u.t)(u.s)), B = -sum t . (u x s)
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

// Exact Montgomery shape-sphere word for 2D orbit (adaptive DOPRI5, shifted start away from equator!)
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

// Compute exact word by integrating adaptive DOPRI5 and shifting start point to max |nz|
static std::pair<std::string, std::string> exact_word_2d(const State3D& y0, ld T, const Masses& m, int perm_idx = 0) {
    static constexpr int perms[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
    // Collect dense adaptive points over [0, T]
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
    // Note: pts.back() is at t=T, identical to pts[0]. Drop duplicate endpoint for cyclic traversal.
    if (pts.size() > 2) pts.pop_back();
    // Find point maximally away from equator nz = 0
    size_t start = 0;
    ld best_nz = -1;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::abs(pts[i][2]) > best_nz) { best_nz = std::abs(pts[i][2]); start = i; }
    }
    const ld c = 1.0L / std::sqrt(3.0L);
    std::string letters;
    size_t N = pts.size();
    // Filter out isolated exact zeros and compare last nonzero sign
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

int main(int argc, char** argv) {
    // Input lines:
    // id dim subspace m1 m2 m3 T precess_rad y0[0..17] (for 2D embedded as x,y,0)
    std::string line;
    std::cout << std::setprecision(12);
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string id, subspace;
        int dim;
        Masses m;
        ld T, precess_rad;
        State3D y0{};
        ss >> id >> dim >> subspace >> m[0] >> m[1] >> m[2] >> T >> precess_rad;
        for (int i = 0; i < 18; ++i) ss >> y0[i];

        // 1. Conservation laws & COM check
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

        // 2. Integrate with exact variational equations over T
        State3D yT{};
        ld Mraw[18][18];
        bool ok_var = integrate_var(y0, T, m, yT, Mraw, 1e-14L);
        ld ET = energy3d(yT, m);
        auto LT = angmom3d(yT, m);
        ld dE = std::abs(ET - E0);
        ld dL = std::sqrt((LT[0]-L0[0])*(LT[0]-L0[0]) + (LT[1]-L0[1])*(LT[1]-L0[1]) + (LT[2]-L0[2])*(LT[2]-L0[2]));

        bool is_rpo = (subspace == "3d_rpo_rotating");
        State3D yT_chk = is_rpo ? rotate_state(yT, u, -precess_rad) : yT;
        ld ret_err = state_dist(yT_chk, y0);

        // Build rotated monodromy matrix M = R_u(-precess_rad) * Mraw (or Mraw if inertial)
        std::vector<std::vector<ld>> M(18, std::vector<ld>(18, 0.0L));
        if (is_rpo) {
            for (int col = 0; col < 18; ++col) {
                State3D v{};
                for (int r = 0; r < 18; ++r) v[r] = Mraw[r][col];
                v = rotate_state(v, u, -precess_rad);
                for (int r = 0; r < 18; ++r) M[r][col] = v[r];
            }
        } else {
            for (int r = 0; r < 18; ++r) for (int c = 0; c < 18; ++c) M[r][c] = Mraw[r][c];
        }

        // Check symplecticity in canonical coords (q, p = m*v):
        // M_can = D M D^-1 where D = diag(1..1, m0,m0,m0,m1,m1,m1,m2,m2,m2)
        ld scale[18];
        for (int b = 0; b < 3; ++b) for (int d = 0; d < 3; ++d) { scale[3*b+d] = 1.0L; scale[9+3*b+d] = m[b]; }
        ld symp_err = 0;
        // (M_can^T J M_can)_{ij} = sum_{k=0..8} (M_can_{k,i} M_can_{k+9,j} - M_can_{k+9,i} M_can_{k,j})
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

        // Compute ALL eigenvalues!
        // For 2D: separate into 12x12 planar and 6x6 vertical blocks!
        ld max_eig_mod_planar = 0, max_eig_mod_vert = 0, max_eig_mod_3d = 0;
        std::vector<cplx> eigs_planar, eigs_vert, eigs_3d;
        if (dim == 2) {
            const int pidx[12] = {0,1,3,4,6,7, 9,10,12,13,15,16};
            const int zidx[6] = {2,5,8, 11,14,17};
            std::vector<std::vector<ld>> Mp(12, std::vector<ld>(12)), Mz(6, std::vector<ld>(6));
            for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) Mp[r][c] = M[pidx[r]][pidx[c]];
            for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) Mz[r][c] = M[zidx[r]][zidx[c]];
            eigs_planar = eigenvalues_qr(Mp);
            eigs_vert = eigenvalues_qr(Mz);
            for (auto& ev : eigs_planar) max_eig_mod_planar = std::max(max_eig_mod_planar, std::abs(ev));
            for (auto& ev : eigs_vert) max_eig_mod_vert = std::max(max_eig_mod_vert, std::abs(ev));
            max_eig_mod_3d = std::max(max_eig_mod_planar, max_eig_mod_vert);
        } else {
            eigs_3d = eigenvalues_qr(M);
            for (auto& ev : eigs_3d) max_eig_mod_3d = std::max(max_eig_mod_3d, std::abs(ev));
        }

        // 3. Check exact minimal period & RPO / permutation / reflection divisors & continuous scan
        // Valid permutations of equal masses:
        std::vector<std::array<int, 3>> valid_perms;
        int all_p[6][3] = {{0,1,2},{1,0,2},{0,2,1},{2,1,0},{1,2,0},{2,0,1}};
        for (auto& p : all_p) {
            if (std::abs(m[p[0]] - m[0]) < 1e-9L && std::abs(m[p[1]] - m[1]) < 1e-9L && std::abs(m[p[2]] - m[2]) < 1e-9L) {
                valid_perms.push_back({p[0], p[1], p[2]});
            }
        }

        // Dense trajectory sampling over [0, T] for geometry, swaps, closest pairs, and period detection!
        const int NS = 2400;
        ld dt = T / NS;
        State3D yc = y0;
        ld rmin_pair[3] = {1e300L, 1e300L, 1e300L};
        ld rmax_pair[3] = {0, 0, 0};
        ld Rmin_body[3] = {1e300L, 1e300L, 1e300L};
        ld Rmax_body[3] = {0, 0, 0};
        int closest_counts[3] = {0, 0, 0};
        int prev_cl = -1, exact_swaps = 0;

        // Track divisor matches at exact fractions T/q for q=2..30
        int best_q_inert = 1, best_q_perm = 1, best_q_rpo = 1, best_q_rpo_perm = 1, best_q_rpo_refl = 1;
        ld best_rpo_angle = 0, best_rpo_err = 1e9L;
        ld best_rpo_perm_angle = 0, best_rpo_perm_err = 1e9L;
        int best_rpo_perm_id[3] = {0,1,2};

        // Check exact divisors directly by integrating to T/q
        for (int q = 2; q <= 24; ++q) {
            State3D yq = y0;
            if (!integrate_state(yq, T / q, m, 1e-14L)) continue;
            // Inertial exact repeat
            if (state_dist(yq, y0) < 1e-4L) best_q_inert = std::max(best_q_inert, q);
            // Permuted inertial repeat
            for (auto& pm : valid_perms) {
                State3D yp = permute_state(yq, pm.data());
                if (state_dist(yp, y0) < 1e-4L) best_q_perm = std::max(best_q_perm, q);
                // Rotating-frame RPO repeat around u
                ld md = 1e9L;
                ld th = optimal_rot_angle(yp, y0, u, md);
                if (pm[0]==0 && pm[1]==1 && pm[2]==2) {
                    if (md < 1e-4L && q > best_q_rpo) {
                        best_q_rpo = q; best_rpo_angle = th; best_rpo_err = md;
                    }
                }
                if (md < 1e-4L && q > best_q_rpo_perm) {
                    best_q_rpo_perm = q; best_rpo_perm_angle = th; best_rpo_perm_err = md;
                    best_rpo_perm_id[0]=pm[0]; best_rpo_perm_id[1]=pm[1]; best_rpo_perm_id[2]=pm[2];
                }
                // With reflection across plane perp to u (for 3D)
                if (dim == 3) {
                    State3D ypr = reflect_plane(yp, u);
                    ld mdr = 1e9L;
                    optimal_rot_angle(ypr, y0, u, mdr);
                    if (mdr < 1e-4L && q > best_q_rpo_refl) best_q_rpo_refl = q;
                }
            }
        }

        // Also scan continuously along [0, T] for any non-divisor / shape-space repeats & geometry
        ld best_cont_rpo_t = T, best_cont_rpo_err = 1e9L;
        for (int k = 0; k <= NS; ++k) {
            ld r[3];
            auto dist3 = [&](int a, int b) {
                ld dx = yc[3*b]-yc[3*a], dy = yc[3*b+1]-yc[3*a+1], dz = yc[3*b+2]-yc[3*a+2];
                return std::sqrt(dx*dx + dy*dy + dz*dz);
            };
            r[0] = dist3(0, 1); r[1] = dist3(0, 2); r[2] = dist3(1, 2);
            for (int p = 0; p < 3; ++p) {
                rmin_pair[p] = std::min(rmin_pair[p], r[p]);
                rmax_pair[p] = std::max(rmax_pair[p], r[p]);
            }
            for (int b = 0; b < 3; ++b) {
                ld Rb = std::sqrt(yc[3*b]*yc[3*b] + yc[3*b+1]*yc[3*b+1] + yc[3*b+2]*yc[3*b+2]);
                Rmin_body[b] = std::min(Rmin_body[b], Rb);
                Rmax_body[b] = std::max(Rmax_body[b], Rb);
            }
            int cl = (r[0] <= r[1] && r[0] <= r[2]) ? 0 : ((r[1] <= r[2]) ? 1 : 2);
            if (k < NS) closest_counts[cl]++;
            if (prev_cl != -1 && cl != prev_cl) exact_swaps++;
            prev_cl = cl;

            if (k > 10 && k < NS - 10) {
                for (auto& pm : valid_perms) {
                    State3D yp = permute_state(yc, pm.data());
                    ld md = 1e9L;
                    optimal_rot_angle(yp, y0, u, md);
                    if (md < best_cont_rpo_err) {
                        best_cont_rpo_err = md;
                        best_cont_rpo_t = k * dt;
                    }
                }
            }
            if (k < NS) integrate_state(yc, dt, m, 1e-13L);
        }

        // Exact free-group word for 2D (all 3 puncture projections!)
        std::string w_red = "-", w_can = "-", w_can_s3 = "-";
        if (dim == 2) {
            auto [wr0, wc0] = exact_word_2d(y0, T, m, 0);
            w_red = wr0.empty() ? "1" : wr0;
            w_can = wc0;
            // Full S3 canonical word across all 3 projections (r12=0, r13=0, r23=0):
            auto [wr1, wc1] = exact_word_2d(y0, T, m, 1);
            auto [wr2, wc2] = exact_word_2d(y0, T, m, 2);
            w_can_s3 = wc0;
            auto better = [](const std::string& a, const std::string& b) {
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            };
            if (better(wc1, w_can_s3)) w_can_s3 = wc1;
            if (better(wc2, w_can_s3)) w_can_s3 = wc2;
        }

        // Emit JSON result
        ld Tstar = T * std::pow(std::abs(E0), 1.5L);
        ld Lambda = Lmag * std::sqrt(std::max(0.0L, -E0));
        std::cout << "{\"id\":\"" << id << "\""
                  << ",\"dim\":" << dim
                  << ",\"subspace\":\"" << subspace << "\""
                  << ",\"masses\":[" << (double)m[0] << "," << (double)m[1] << "," << (double)m[2] << "]"
                  << ",\"T\":" << (double)T
                  << ",\"E\":" << (double)E0
                  << ",\"L\":" << (double)Lmag
                  << ",\"Tstar\":" << (double)Tstar
                  << ",\"Lambda\":" << (double)Lambda
                  << ",\"rcm_err\":" << (double)rcm_err
                  << ",\"vcm_err\":" << (double)vcm_err
                  << ",\"ret_err_ld\":" << (double)ret_err
                  << ",\"dE\":" << (double)dE
                  << ",\"dL\":" << (double)dL
                  << ",\"symp_err\":" << (double)symp_err
                  << ",\"max_eig_planar\":" << (double)max_eig_mod_planar
                  << ",\"max_eig_vert\":" << (double)max_eig_mod_vert
                  << ",\"max_eig_3d\":" << (double)max_eig_mod_3d
                  << ",\"q_inert\":" << best_q_inert
                  << ",\"q_perm\":" << best_q_perm
                  << ",\"q_rpo\":" << best_q_rpo
                  << ",\"q_rpo_perm\":" << best_q_rpo_perm
                  << ",\"q_rpo_refl\":" << best_q_rpo_refl
                  << ",\"rpo_angle\":" << (double)best_rpo_angle
                  << ",\"rpo_perm_angle\":" << (double)best_rpo_perm_angle
                  << ",\"cont_rpo_t\":" << (double)best_cont_rpo_t
                  << ",\"cont_rpo_err\":" << (double)best_cont_rpo_err
                  << ",\"exact_swaps\":" << exact_swaps
                  << ",\"closest_counts\":[" << closest_counts[0] << "," << closest_counts[1] << "," << closest_counts[2] << "]"
                  << ",\"rmin_pair\":[" << (double)rmin_pair[0] << "," << (double)rmin_pair[1] << "," << (double)rmin_pair[2] << "]"
                  << ",\"rmax_pair\":[" << (double)rmax_pair[0] << "," << (double)rmax_pair[1] << "," << (double)rmax_pair[2] << "]"
                  << ",\"Rmin_body\":[" << (double)Rmin_body[0] << "," << (double)Rmin_body[1] << "," << (double)Rmin_body[2] << "]"
                  << ",\"Rmax_body\":[" << (double)Rmax_body[0] << "," << (double)Rmax_body[1] << "," << (double)Rmax_body[2] << "]"
                  << ",\"word_red\":\"" << w_red << "\""
                  << ",\"word_can\":\"" << w_can << "\""
                  << ",\"word_can_s3\":\"" << w_can_s3 << "\"";
        // Also dump nontrivial eigenvalues (sorted by modulus descending)
        auto dump_eigs = [](const char* name, const std::vector<cplx>& evs) {
            std::cout << ",\"" << name << "\":[";
            for (size_t i = 0; i < evs.size(); ++i) {
                std::cout << "[" << (double)evs[i].real() << "," << (double)evs[i].imag() << "]"
                          << (i + 1 < evs.size() ? "," : "");
            }
            std::cout << "]";
        };
        if (dim == 2) {
            dump_eigs("eigs_planar", eigs_planar);
            dump_eigs("eigs_vert", eigs_vert);
        } else {
            dump_eigs("eigs_3d", eigs_3d);
        }
        std::cout << "}\n" << std::flush;
    }
    return 0;
}
