// Independent 80-bit long double analytical variational 24x24 monodromy auditor
// for all 11 4-body orbits (4B-001..4B-011) in catalog.json, plus optional
// long-time empirical perturbation test for weak hyperbolic modes.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using ld = long double;
using cplx = std::complex<ld>;
using State = std::array<ld, 24>;
using Masses = std::array<ld, 4>;
using VarState = std::array<ld, 24 + 24 * 24>;

static void deriv_state(const State& s, const Masses& m, State& ds) {
    for (int i = 0; i < 12; ++i) ds[i] = s[12 + i];
    for (int i = 12; i < 24; ++i) ds[i] = 0.0L;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            ld dx = s[3 * j] - s[3 * i], dy = s[3 * j + 1] - s[3 * i + 1], dz = s[3 * j + 2] - s[3 * i + 2];
            ld r2 = dx * dx + dy * dy + dz * dz;
            ld ir = 1.0L / std::sqrt(r2), ir3 = ir * ir * ir;
            ds[12 + 3 * i]     += m[j] * dx * ir3;
            ds[12 + 3 * i + 1] += m[j] * dy * ir3;
            ds[12 + 3 * i + 2] += m[j] * dz * ir3;
            ds[12 + 3 * j]     -= m[i] * dx * ir3;
            ds[12 + 3 * j + 1] -= m[i] * dy * ir3;
            ds[12 + 3 * j + 2] -= m[i] * dz * ir3;
        }
    }
}

static void deriv_and_jac(const State& s, const Masses& m, State& ds, ld J[24][24]) {
    for (int i = 0; i < 12; ++i) ds[i] = s[12 + i];
    for (int i = 12; i < 24; ++i) ds[i] = 0.0L;
    for (int r = 0; r < 24; ++r)
        for (int c = 0; c < 24; ++c) J[r][c] = 0.0L;
    for (int i = 0; i < 12; ++i) J[i][12 + i] = 1.0L;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            ld dx = s[3 * j] - s[3 * i], dy = s[3 * j + 1] - s[3 * i + 1], dz = s[3 * j + 2] - s[3 * i + 2];
            ld r2 = dx * dx + dy * dy + dz * dz;
            ld ir = 1.0L / std::sqrt(r2), ir3 = ir * ir * ir, ir5 = ir3 * ir * ir;
            ld d[3] = {dx, dy, dz};
            for (int a = 0; a < 3; ++a) {
                ds[12 + 3 * i + a] += m[j] * d[a] * ir3;
                ds[12 + 3 * j + a] -= m[i] * d[a] * ir3;
            }
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    ld T = ((a == b) ? ir3 : 0.0L) - 3.0L * d[a] * d[b] * ir5;
                    J[12 + 3 * i + a][3 * j + b] += m[j] * T;
                    J[12 + 3 * i + a][3 * i + b] -= m[j] * T;
                    J[12 + 3 * j + a][3 * i + b] += m[i] * T;
                    J[12 + 3 * j + a][3 * j + b] -= m[i] * T;
                }
            }
        }
    }
}

static void deriv_var(const VarState& z, const Masses& m, VarState& dz) {
    State s, ds;
    for (int i = 0; i < 24; ++i) s[i] = z[i];
    ld J[24][24];
    deriv_and_jac(s, m, ds, J);
    for (int i = 0; i < 24; ++i) dz[i] = ds[i];
    for (int c = 0; c < 24; ++c) {
        for (int r = 0; r < 24; ++r) {
            ld sum = 0.0L;
            if (r < 12) sum = z[24 + (12 + r) * 24 + c];
            else for (int k = 0; k < 12; ++k) sum += J[r][k] * z[24 + k * 24 + c];
            dz[24 + r * 24 + c] = sum;
        }
    }
}

static bool integrate_state(State& y, ld T, const Masses& m, ld tol = 1e-13L, int max_steps = 800000) {
    constexpr ld a21 = 1.0L / 5, a31 = 3.0L / 40, a32 = 9.0L / 40,
        a41 = 44.0L / 45, a42 = -56.0L / 15, a43 = 32.0L / 9,
        a51 = 19372.0L / 6561, a52 = -25360.0L / 2187, a53 = 64448.0L / 6561, a54 = -212.0L / 729,
        a61 = 9017.0L / 3168, a62 = -355.0L / 33, a63 = 46732.0L / 5247, a64 = 49.0L / 176, a65 = -5103.0L / 18656,
        b1 = 35.0L / 384, b3 = 500.0L / 1113, b4 = 125.0L / 192, b5 = -2187.0L / 6784, b6 = 11.0L / 84,
        e1 = 71.0L / 57600, e3 = -71.0L / 16695, e4 = 71.0L / 1920, e5 = -17253.0L / 339200, e6 = 22.0L / 525, e7 = -1.0L / 40;
    ld t = 0.0L, h = std::min(1e-3L, T * 0.01L);
    State k1, k2, k3, k4, k5, k6, k7, yt, yn;
    for (int step = 0; t < T - 1e-15L && step < max_steps; ++step) {
        if (t + h > T) h = T - t;
        deriv_state(y, m, k1);
        for (int i = 0; i < 24; ++i) yt[i] = y[i] + h * a21 * k1[i]; deriv_state(yt, m, k2);
        for (int i = 0; i < 24; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]); deriv_state(yt, m, k3);
        for (int i = 0; i < 24; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); deriv_state(yt, m, k4);
        for (int i = 0; i < 24; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); deriv_state(yt, m, k5);
        for (int i = 0; i < 24; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); deriv_state(yt, m, k6);
        for (int i = 0; i < 24; ++i) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]); deriv_state(yn, m, k7);
        ld err = 0.0L;
        for (int i = 0; i < 24; ++i) {
            ld ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            ld sc = tol * (1.0L + std::max(std::abs(y[i]), std::abs(yn[i])));
            err += (ei / sc) * (ei / sc);
        }
        err = std::sqrt(err / 24.0L);
        if (err <= 1.0L) {
            y = yn; t += h;
            for (int i = 0; i < 12; ++i) if (std::abs(y[i]) > 100.0L) return false;
            h = std::min(h * (err == 0 ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.05L);
        } else {
            h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
            if (h < 1e-16L) return false;
        }
    }
    return t >= T - 1e-14L;
}

static bool integrate_var(VarState& y, ld T, const Masses& m, ld tol = 1e-16L) {
    constexpr ld a21 = 1.0L / 5, a31 = 3.0L / 40, a32 = 9.0L / 40,
        a41 = 44.0L / 45, a42 = -56.0L / 15, a43 = 32.0L / 9,
        a51 = 19372.0L / 6561, a52 = -25360.0L / 2187, a53 = 64448.0L / 6561, a54 = -212.0L / 729,
        a61 = 9017.0L / 3168, a62 = -355.0L / 33, a63 = 46732.0L / 5247, a64 = 49.0L / 176, a65 = -5103.0L / 18656,
        b1 = 35.0L / 384, b3 = 500.0L / 1113, b4 = 125.0L / 192, b5 = -2187.0L / 6784, b6 = 11.0L / 84,
        e1 = 71.0L / 57600, e3 = -71.0L / 16695, e4 = 71.0L / 1920, e5 = -17253.0L / 339200, e6 = 22.0L / 525, e7 = -1.0L / 40;
    constexpr int D = 600;
    ld t = 0.0L, h = 1e-3L;
    VarState k1, k2, k3, k4, k5, k6, k7, yt, yn;
    for (int steps = 0; t < T - 1e-17L && steps < 600000; ++steps) {
        if (t + h > T) h = T - t;
        deriv_var(y, m, k1);
        for (int i = 0; i < D; ++i) yt[i] = y[i] + h * a21 * k1[i]; deriv_var(yt, m, k2);
        for (int i = 0; i < D; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]); deriv_var(yt, m, k3);
        for (int i = 0; i < D; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); deriv_var(yt, m, k4);
        for (int i = 0; i < D; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); deriv_var(yt, m, k5);
        for (int i = 0; i < D; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); deriv_var(yt, m, k6);
        for (int i = 0; i < D; ++i) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]); deriv_var(yn, m, k7);
        ld err = 0.0L;
        for (int i = 0; i < 24; ++i) {
            ld ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            ld sc = tol * (1.0L + std::max(std::abs(y[i]), std::abs(yn[i])));
            err += (ei / sc) * (ei / sc);
        }
        err = std::sqrt(err / 24.0L);
        if (err <= 1.0L) {
            y = yn; t += h;
            h = std::min(h * (err == 0 ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.04L);
        } else {
            h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
            if (h < 1e-16L) return false;
        }
    }
    return t >= T - 1e-15L;
}

static std::vector<cplx> qr_eigs(std::vector<std::vector<ld>> H) {
    int n = (int)H.size();
    for (int k = 0; k < n - 2; ++k) {
        ld scale = 0;
        for (int i = k + 1; i < n; ++i) scale += std::abs(H[i][k]);
        if (scale < 1e-35L) continue;
        ld norm2 = 0;
        std::vector<ld> v(n - k - 1);
        for (int i = k + 1; i < n; ++i) { v[i - k - 1] = H[i][k] / scale; norm2 += v[i - k - 1] * v[i - k - 1]; }
        ld alpha = (v[0] >= 0 ? -std::sqrt(norm2) : std::sqrt(norm2));
        v[0] -= alpha;
        ld v2 = 0; for (ld x : v) v2 += x * x;
        if (v2 < 1e-35L) continue;
        ld beta = 2.0L / v2;
        for (int j = k; j < n; ++j) {
            ld dot = 0; for (int i = k + 1; i < n; ++i) dot += v[i - k - 1] * H[i][j];
            for (int i = k + 1; i < n; ++i) H[i][j] -= beta * dot * v[i - k - 1];
        }
        for (int i = 0; i < n; ++i) {
            ld dot = 0; for (int j = k + 1; j < n; ++j) dot += H[i][j] * v[j - k - 1];
            for (int j = k + 1; j < n; ++j) H[i][j] -= beta * dot * v[j - k - 1];
        }
        for (int i = k + 2; i < n; ++i) H[i][k] = 0;
    }
    std::vector<cplx> eigs;
    ld norm = 0;
    for (int i = 0; i < n; ++i)
        for (int j = std::max(0, i - 1); j < n; ++j) norm += std::abs(H[i][j]);
    const ld eps = 1.0842021724855044e-19L;
    int nn = n - 1, iter = 0;
    while (nn >= 0) {
        int l = nn;
        while (l > 0) {
            ld s = std::abs(H[l - 1][l - 1]) + std::abs(H[l][l]);
            if (s == 0) s = norm;
            if (std::abs(H[l][l - 1]) <= eps * s) { H[l][l - 1] = 0; break; }
            --l;
        }
        if (l == nn) { eigs.emplace_back(H[nn][nn], 0.0L); --nn; iter = 0; }
        else if (l == nn - 1) {
            ld p = 0.5L * (H[nn - 1][nn - 1] - H[nn][nn]);
            ld q = p * p + H[nn][nn - 1] * H[nn - 1][nn];
            ld mid = 0.5L * (H[nn - 1][nn - 1] + H[nn][nn]);
            if (q >= 0) { ld sq = std::sqrt(q); eigs.emplace_back(mid + sq, 0); eigs.emplace_back(mid - sq, 0); }
            else { ld sq = std::sqrt(-q); eigs.emplace_back(mid, sq); eigs.emplace_back(mid, -sq); }
            nn -= 2; iter = 0;
        } else {
            if (iter > 500) break;
            ld x = H[nn][nn], y = H[nn - 1][nn - 1], w = H[nn][nn - 1] * H[nn - 1][nn];
            if (iter == 10 || iter == 20 || iter == 40 || iter == 80 || iter == 160) {
                ld s = std::abs(H[nn][nn - 1]) + std::abs(H[nn - 1][nn - 2]);
                x = y = 0.75L * s; w = -0.4375L * s * s;
            }
            ++iter;
            int m = nn - 2; ld p = 0, q = 0, r = 0;
            for (; m >= l; --m) {
                ld z = H[m][m]; r = x - z; ld s = y - z;
                p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
                q = H[m + 1][m + 1] - z - r - s; r = H[m + 2][m + 1];
                ld sc = std::abs(p) + std::abs(q) + std::abs(r); p /= sc; q /= sc; r /= sc;
                if (m == l) break;
                if (std::abs(H[m][m - 1]) * (std::abs(q) + std::abs(r)) <=
                    eps * std::abs(p) * (std::abs(H[m - 1][m - 1]) + std::abs(z) + std::abs(H[m + 1][m + 1]))) break;
            }
            for (int i = m + 2; i <= nn; ++i) { H[i][i - 2] = 0; if (i != m + 2) H[i][i - 3] = 0; }
            for (int k = m; k <= nn - 1; ++k) {
                bool notlast = (k != nn - 1);
                if (k != m) {
                    p = H[k][k - 1]; q = H[k + 1][k - 1]; r = notlast ? H[k + 2][k - 1] : 0;
                    ld sc = std::abs(p) + std::abs(q) + std::abs(r); if (sc == 0) continue; p /= sc; q /= sc; r /= sc;
                }
                ld s = std::sqrt(p * p + q * q + r * r); if (p < 0) s = -s;
                if (s != 0) {
                    if (k != m) H[k][k - 1] = -s * (std::abs(H[k][k - 1]) + std::abs(H[k + 1][k - 1]) + (notlast ? std::abs(H[k + 2][k - 1]) : 0));
                    else if (l != m) H[k][k - 1] = -H[k][k - 1];
                    p += s; ld xp = p / s, yp = q / s, zp = r / s; q /= p; r /= p;
                    for (int j = k; j < n; ++j) {
                        ld pp = H[k][j] + q * H[k + 1][j];
                        if (notlast) { pp += r * H[k + 2][j]; H[k + 2][j] -= pp * zp; }
                        H[k][j] -= pp * xp; H[k + 1][j] -= pp * yp;
                    }
                    int rmax = std::min(nn, k + 3);
                    for (int i = 0; i <= rmax; ++i) {
                        ld pp = xp * H[i][k] + yp * H[i][k + 1];
                        if (notlast) { pp += zp * H[i][k + 2]; H[i][k + 2] -= pp * r; }
                        H[i][k] -= pp; H[i][k + 1] -= pp * q;
                    }
                }
            }
        }
    }
    return eigs;
}

static void rot_mat3(const ld u[3], ld ang, ld R[3][3]) {
    ld c = std::cos(ang), s = std::sin(ang), omc = 1.0L - c;
    R[0][0] = c + u[0] * u[0] * omc; R[0][1] = u[0] * u[1] * omc - u[2] * s; R[0][2] = u[0] * u[2] * omc + u[1] * s;
    R[1][0] = u[1] * u[0] * omc + u[2] * s; R[1][1] = c + u[1] * u[1] * omc; R[1][2] = u[1] * u[2] * omc - u[0] * s;
    R[2][0] = u[2] * u[0] * omc - u[1] * s; R[2][1] = u[2] * u[1] * omc + u[0] * s; R[2][2] = c + u[2] * u[2] * omc;
}

static State rotate_state(const State& s, const ld u[3], ld ang) {
    ld R[3][3]; rot_mat3(u, ang, R);
    State out{};
    for (int blk = 0; blk < 8; ++blk)
        for (int a = 0; a < 3; ++a)
            out[3 * blk + a] = R[a][0] * s[3 * blk] + R[a][1] * s[3 * blk + 1] + R[a][2] * s[3 * blk + 2];
    return out;
}

int main(int argc, char** argv) {
    bool run_long = (argc > 1 && std::string(argv[1]) == "--long-kick");
    std::string id; int dim; Masses m; ld T, prec; State y0;
    while (std::cin >> id >> dim >> m[0] >> m[1] >> m[2] >> m[3] >> T >> prec) {
        for (int i = 0; i < 24; ++i) std::cin >> y0[i];
        ld L[3] = {0, 0, 0};
        for (int b = 0; b < 4; ++b) {
            L[0] += m[b] * (y0[3 * b + 1] * y0[12 + 3 * b + 2] - y0[3 * b + 2] * y0[12 + 3 * b + 1]);
            L[1] += m[b] * (y0[3 * b + 2] * y0[12 + 3 * b + 0] - y0[3 * b + 0] * y0[12 + 3 * b + 2]);
            L[2] += m[b] * (y0[3 * b + 0] * y0[12 + 3 * b + 1] - y0[3 * b + 1] * y0[12 + 3 * b + 0]);
        }
        ld Lm = std::sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
        ld u[3] = {0, 0, 1};
        if (Lm > 1e-12L) { u[0] = L[0] / Lm; u[1] = L[1] / Lm; u[2] = L[2] / Lm; }

        VarState z{};
        for (int i = 0; i < 24; ++i) z[i] = y0[i];
        for (int i = 0; i < 24; ++i) z[24 + i * 24 + i] = 1.0L;
        if (!integrate_var(z, T, m, 1e-16L)) {
            std::cout << "{\"id\":\"" << id << "\",\"ok\":false}\n";
            continue;
        }
        ld R[3][3]; rot_mat3(u, -prec, R);
        ld ret = 0.0L;
        for (int blk = 0; blk < 8; ++blk) {
            ld v[3] = {z[3 * blk], z[3 * blk + 1], z[3 * blk + 2]};
            for (int a = 0; a < 3; ++a) {
                ld rv = R[a][0] * v[0] + R[a][1] * v[1] + R[a][2] * v[2];
                ret += (rv - y0[3 * blk + a]) * (rv - y0[3 * blk + a]);
            }
        }
        std::vector<std::vector<ld>> M(24, std::vector<ld>(24, 0.0L));
        for (int c = 0; c < 24; ++c)
            for (int blk = 0; blk < 8; ++blk) {
                ld v[3] = {z[24 + (3 * blk) * 24 + c], z[24 + (3 * blk + 1) * 24 + c], z[24 + (3 * blk + 2) * 24 + c]};
                for (int a = 0; a < 3; ++a) M[3 * blk + a][c] = R[a][0] * v[0] + R[a][1] * v[1] + R[a][2] * v[2];
            }
        ld symp = 0.0L;
        for (int i = 0; i < 24; ++i)
            for (int j = 0; j < 24; ++j) {
                ld s = 0.0L;
                ld sci = (i >= 12) ? (1.0L / m[(i - 12) / 3]) : 1.0L;
                ld scj = (j >= 12) ? (1.0L / m[(j - 12) / 3]) : 1.0L;
                for (int k = 0; k < 12; ++k) {
                    ld mk = m[k / 3];
                    s += (M[k][i] * sci) * (mk * M[k + 12][j] * scj) - (mk * M[k + 12][i] * sci) * (M[k][j] * scj);
                }
                ld tgt = (j == i + 12) ? 1.0L : ((i == j + 12) ? -1.0L : 0.0L);
                symp = std::max(symp, std::abs(s - tgt));
            }
        ld trM = 0.0L; for (int i = 0; i < 24; ++i) trM += M[i][i];
        ld tr_nontriv = trM - (6.0L + 6.0L * std::cos(prec));
        auto eigs = qr_eigs(M);

        std::vector<cplx> expected;
        for (int i = 0; i < 6; ++i) expected.emplace_back(1.0L, 0.0L);
        cplx rp(std::cos(prec), std::sin(prec)), rm(std::cos(prec), -std::sin(prec));
        for (int i = 0; i < 3; ++i) { expected.push_back(rp); expected.push_back(rm); }
        std::vector<bool> used(24, false);
        ld triv_err = 0.0L;
        for (auto& tgt : expected) {
            int bi = -1; ld bd = 1e30L;
            for (int i = 0; i < 24; ++i) if (!used[i]) { ld d = std::abs(eigs[i] - tgt); if (d < bd) { bd = d; bi = i; } }
            used[bi] = true; triv_err = std::max(triv_err, bd);
        }
        std::vector<cplx> rem;
        for (int i = 0; i < 24; ++i) if (!used[i]) rem.push_back(eigs[i]);
        ld max_mod = 0.0L;
        for (auto& zc : rem) max_mod = std::max(max_mod, std::abs(zc));

        // Pair remaining eigenvalues and compute trace-invariant 6th index when 5 clear elliptic pairs exist
        ld sum_ell = 0.0L; int n_ell = 0;
        for (auto& zc : rem) if (std::abs(zc.imag()) > 1e-2L) { sum_ell += zc.real(); n_ell++; }
        ld k_trace_last = (n_ell == 10) ? (tr_nontriv - sum_ell) : 0.0L;
        bool has_hyperbolic = (n_ell == 10 && k_trace_last > 2.0L) || (max_mod > 1.0L + 1e-5L);

        // Extract all 6 stability indices (2*Re for elliptic pairs, lambda+1/lambda for real pairs)
        std::vector<ld> ks;
        std::vector<bool> pused(rem.size(), false);
        for (size_t i = 0; i < rem.size(); ++i) {
            if (pused[i]) continue;
            pused[i] = true;
            int mate = -1; ld bd = 1e30L;
            for (size_t j = i + 1; j < rem.size(); ++j) {
                if (pused[j]) continue;
                ld d = std::min(std::abs(rem[j] - std::conj(rem[i])), std::abs(rem[j] * rem[i] - cplx(1.0L, 0.0L)));
                if (d < bd) { bd = d; mate = (int)j; }
            }
            if (mate >= 0) { pused[mate] = true; ks.push_back((rem[i] + rem[mate]).real()); }
        }
        if (n_ell == 10) {
            // Replace the near-2 pair with the exact trace-invariant k6 (immune to Jordan splitting)
            auto it = std::max_element(ks.begin(), ks.end());
            if (it != ks.end()) *it = k_trace_last;
        }
        std::sort(ks.begin(), ks.end(), [](ld a, ld b) { return std::abs(a) > std::abs(b); });
        ld exact_margin = 2.0L - std::abs(ks[0]);

        int surv_long = -1; ld drift_500 = -1, drift_2000 = -1, drift_4000 = -1;
        if (run_long && id == "4B-008") {
            State yp = y0;
            ld kick = 1e-6L;
            yp[12] += kick; yp[15] -= kick; yp[14] += kick; yp[17] -= kick;
            surv_long = 0;
            for (int k = 1; k <= 4000; ++k) {
                if (!integrate_state(yp, T, m, 1e-12L, 400000)) break;
                State yr = rotate_state(yp, u, -prec * k);
                ld d2 = 0; for (int i = 0; i < 24; ++i) d2 += (yr[i] - y0[i]) * (yr[i] - y0[i]);
                ld d = std::sqrt(d2);
                if (k == 500) drift_500 = d;
                if (k == 2000) drift_2000 = d;
                if (k == 4000) drift_4000 = d;
                if (d > 5.0L) break;
                surv_long = k;
            }
        }

        std::cout << std::setprecision(13)
                  << "{\"id\":\"" << id << "\",\"dim\":" << dim
                  << ",\"ret_err\":" << (double)std::sqrt(ret)
                  << ",\"symp_err\":" << (double)symp
                  << ",\"triv_err\":" << (double)triv_err
                  << ",\"max_eig_mod\":" << (double)max_mod
                  << ",\"n_ell_eigs\":" << n_ell
                  << ",\"k_trace_last\":" << (double)k_trace_last
                  << ",\"k_trace_minus_2\":" << (double)(k_trace_last - 2.0L)
                  << ",\"exact_margin\":" << (double)exact_margin
                  << ",\"hyperbolic\":" << (has_hyperbolic ? "true" : "false")
                  << ",\"k_indices\":[";
        for (size_t i = 0; i < ks.size(); ++i) std::cout << (double)ks[i] << (i + 1 < ks.size() ? "," : "");
        std::cout << "]";
        if (surv_long >= 0) {
            std::cout << ",\"surv_4000T\":" << surv_long
                      << ",\"drift_500T\":" << (double)drift_500
                      << ",\"drift_2000T\":" << (double)drift_2000
                      << ",\"drift_4000T\":" << (double)drift_4000;
        }
        std::cout << "}\n";
    }
    return 0;
}
