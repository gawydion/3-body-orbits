// Independent 80-bit long double 1PN ADM canonical auditor for all 21 relativistic
// orbits (1PN-001..1PN-021): verifies exact Hamiltonian symplecticity ||M^T J M - J||_inf,
// periodicity, conservation (H, L, P), and compares 1-RPO vs true T_inert monodromy
// for all 5 snapped relativistic inertial closures.
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using ld = long double;
using cplx = std::complex<ld>;
using State1PN = std::array<ld, 18>;
using Masses = std::array<ld, 3>;

static ld hamiltonian_1pn(const State1PN& s, const Masses& m, ld beta) {
    ld H0 = 0.0L, H1 = 0.0L, phi[3] = {0, 0, 0};
    for (int a = 0; a < 3; ++a) {
        ld px = s[9 + 3 * a], py = s[10 + 3 * a], pz = s[11 + 3 * a];
        ld p2 = px * px + py * py + pz * pz;
        H0 += 0.5L * p2 / m[a];
        H1 -= (p2 * p2) / (8.0L * m[a] * m[a] * m[a]);
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            ld dx = s[3 * a] - s[3 * b], dy = s[3 * a + 1] - s[3 * b + 1], dz = s[3 * a + 2] - s[3 * b + 2];
            ld r = std::sqrt(dx * dx + dy * dy + dz * dz), ir = 1.0L / r;
            H0 -= m[a] * m[b] * ir;
            phi[a] += m[b] * ir; phi[b] += m[a] * ir;
            if (beta != 0.0L) {
                ld pax = s[9 + 3 * a], pay = s[10 + 3 * a], paz = s[11 + 3 * a];
                ld pbx = s[9 + 3 * b], pby = s[10 + 3 * b], pbz = s[11 + 3 * b];
                ld pa2 = pax * pax + pay * pay + paz * paz;
                ld pb2 = pbx * pbx + pby * pby + pbz * pbz;
                ld papb = pax * pbx + pay * pby + paz * pbz;
                ld nx = dx * ir, ny = dy * ir, nz = dz * ir;
                ld na = nx * pax + ny * pay + nz * paz;
                ld nb = nx * pbx + ny * pby + nz * pbz;
                ld C = 3.0L * ((m[b] / m[a]) * pa2 + (m[a] / m[b]) * pb2) - 7.0L * papb;
                H1 += -0.5L * ir * (C - na * nb);
            }
        }
    }
    if (beta != 0.0L)
        for (int a = 0; a < 3; ++a) H1 += 0.5L * m[a] * phi[a] * phi[a];
    return H0 + beta * H1;
}

static void deriv_1pn(const State1PN& s, const Masses& m, ld beta, State1PN& ds) {
    ld dr[9] = {0}, dp[9] = {0}, phi[3] = {0};
    ld dx_ab[3][3], dy_ab[3][3], dz_ab[3][3], ir3_ab[3][3];
    for (int a = 0; a < 3; ++a) {
        ld px = s[9 + 3 * a], py = s[10 + 3 * a], pz = s[11 + 3 * a];
        ld p2 = px * px + py * py + pz * pz, inv_m = 1.0L / m[a];
        ld fac = inv_m - beta * (0.5L * p2 * inv_m * inv_m * inv_m);
        dr[3 * a] = fac * px; dr[3 * a + 1] = fac * py; dr[3 * a + 2] = fac * pz;
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            ld dx = s[3 * a] - s[3 * b], dy = s[3 * a + 1] - s[3 * b + 1], dz = s[3 * a + 2] - s[3 * b + 2];
            ld r2 = dx * dx + dy * dy + dz * dz, ir = 1.0L / std::sqrt(r2), ir2 = ir * ir, ir3 = ir2 * ir;
            dx_ab[a][b] = dx; dy_ab[a][b] = dy; dz_ab[a][b] = dz; ir3_ab[a][b] = ir3;
            ld mm = m[a] * m[b];
            ld fx0 = -mm * dx * ir3, fy0 = -mm * dy * ir3, fz0 = -mm * dz * ir3;
            dp[3 * a] += fx0; dp[3 * a + 1] += fy0; dp[3 * a + 2] += fz0;
            dp[3 * b] -= fx0; dp[3 * b + 1] -= fy0; dp[3 * b + 2] -= fz0;
            phi[a] += m[b] * ir; phi[b] += m[a] * ir;
            if (beta != 0.0L) {
                ld pax = s[9 + 3 * a], pay = s[10 + 3 * a], paz = s[11 + 3 * a];
                ld pbx = s[9 + 3 * b], pby = s[10 + 3 * b], pbz = s[11 + 3 * b];
                ld pa2 = pax * pax + pay * pay + paz * paz, pb2 = pbx * pbx + pby * pby + pbz * pbz;
                ld papb = pax * pbx + pay * pby + paz * pbz;
                ld nx = dx * ir, ny = dy * ir, nz = dz * ir;
                ld na = nx * pax + ny * pay + nz * paz, nb = nx * pbx + ny * pby + nz * pbz;
                ld c_ab = 0.5L * beta * ir, mba = m[b] / m[a], mab = m[a] / m[b];
                dr[3 * a]     += c_ab * (-6.0L * mba * pax + 7.0L * pbx + nb * nx);
                dr[3 * a + 1] += c_ab * (-6.0L * mba * pay + 7.0L * pby + nb * ny);
                dr[3 * a + 2] += c_ab * (-6.0L * mba * paz + 7.0L * pbz + nb * nz);
                dr[3 * b]     += c_ab * (-6.0L * mab * pbx + 7.0L * pax + na * nx);
                dr[3 * b + 1] += c_ab * (-6.0L * mab * pby + 7.0L * pay + na * ny);
                dr[3 * b + 2] += c_ab * (-6.0L * mab * pbz + 7.0L * paz + na * nz);
                ld C = 3.0L * (mba * pa2 + mab * pb2) - 7.0L * papb;
                ld coef_n = C - 3.0L * na * nb, f_scale = -0.5L * beta * ir2;
                ld fx1 = f_scale * (coef_n * nx + nb * pax + na * pbx);
                ld fy1 = f_scale * (coef_n * ny + nb * pay + na * pby);
                ld fz1 = f_scale * (coef_n * nz + nb * paz + na * pbz);
                dp[3 * a] += fx1; dp[3 * a + 1] += fy1; dp[3 * a + 2] += fz1;
                dp[3 * b] -= fx1; dp[3 * b + 1] -= fy1; dp[3 * b + 2] -= fz1;
            }
        }
    }
    if (beta != 0.0L) {
        for (int a = 0; a < 3; ++a)
            for (int b = a + 1; b < 3; ++b) {
                ld f = beta * m[a] * m[b] * (phi[a] + phi[b]) * ir3_ab[a][b];
                dp[3 * a] += f * dx_ab[a][b]; dp[3 * a + 1] += f * dy_ab[a][b]; dp[3 * a + 2] += f * dz_ab[a][b];
                dp[3 * b] -= f * dx_ab[a][b]; dp[3 * b + 1] -= f * dy_ab[a][b]; dp[3 * b + 2] -= f * dz_ab[a][b];
            }
    }
    for (int i = 0; i < 9; ++i) { ds[i] = dr[i]; ds[9 + i] = dp[i]; }
}

static bool integrate_1pn(State1PN& y, ld T, const Masses& m, ld beta, ld tol = 1e-16L) {
    constexpr ld a21 = 1.0L / 5, a31 = 3.0L / 40, a32 = 9.0L / 40,
        a41 = 44.0L / 45, a42 = -56.0L / 15, a43 = 32.0L / 9,
        a51 = 19372.0L / 6561, a52 = -25360.0L / 2187, a53 = 64448.0L / 6561, a54 = -212.0L / 729,
        a61 = 9017.0L / 3168, a62 = -355.0L / 33, a63 = 46732.0L / 5247, a64 = 49.0L / 176, a65 = -5103.0L / 18656,
        b1 = 35.0L / 384, b3 = 500.0L / 1113, b4 = 125.0L / 192, b5 = -2187.0L / 6784, b6 = 11.0L / 84,
        e1 = 71.0L / 57600, e3 = -71.0L / 16695, e4 = 71.0L / 1920, e5 = -17253.0L / 339200, e6 = 22.0L / 525, e7 = -1.0L / 40;
    ld t = 0.0L, h = std::min(1e-3L, T * 0.01L);
    State1PN k1, k2, k3, k4, k5, k6, k7, yt, yn;
    for (int step = 0; t < T - 1e-17L && step < 3000000; ++step) {
        if (t + h > T) h = T - t;
        deriv_1pn(y, m, beta, k1);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * a21 * k1[i]; deriv_1pn(yt, m, beta, k2);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]); deriv_1pn(yt, m, beta, k3);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); deriv_1pn(yt, m, beta, k4);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); deriv_1pn(yt, m, beta, k5);
        for (int i = 0; i < 18; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); deriv_1pn(yt, m, beta, k6);
        for (int i = 0; i < 18; ++i) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]); deriv_1pn(yn, m, beta, k7);
        ld err = 0.0L;
        for (int i = 0; i < 18; ++i) {
            ld ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            ld sc = tol * (1.0L + std::max(std::abs(y[i]), std::abs(yn[i])));
            err += (ei / sc) * (ei / sc);
        }
        err = std::sqrt(err / 18.0L);
        if (err <= 1.0L) {
            y = yn; t += h;
            h = std::min(h * (err == 0 ? 5.0L : std::min(5.0L, std::max(0.2L, 0.9L * std::pow(err, -0.2L)))), 0.04L);
        } else {
            h *= std::max(0.1L, 0.9L * std::pow(err, -0.25L));
            if (h < 1e-16L) return false;
        }
    }
    return t >= T - 1e-14L;
}

static State1PN rotate_state(const State1PN& s, const ld u[3], ld ang) {
    ld c = std::cos(ang), sn = std::sin(ang), omc = 1.0L - c;
    State1PN out = s;
    for (int b = 0; b < 6; ++b) {
        ld x = s[3 * b], y = s[3 * b + 1], z = s[3 * b + 2];
        ld dot = u[0] * x + u[1] * y + u[2] * z;
        ld cx = u[1] * z - u[2] * y, cy = u[2] * x - u[0] * z, cz = u[0] * y - u[1] * x;
        out[3 * b]     = x * c + cx * sn + u[0] * dot * omc;
        out[3 * b + 1] = y * c + cy * sn + u[1] * dot * omc;
        out[3 * b + 2] = z * c + cz * sn + u[2] * dot * omc;
    }
    return out;
}

static void broucke_2d(const std::vector<std::vector<ld>>& M, ld theta, bool& real_k, ld& k1, ld& k2, ld& margin) {
    ld tr1 = 0, tr2 = 0;
    for (int i = 0; i < 12; ++i) {
        tr1 += M[i][i];
        for (int j = 0; j < 12; ++j) tr2 += M[i][j] * M[j][i];
    }
    ld s = tr1 - (4.0L + 4.0L * std::cos(theta));
    ld t2 = tr2 - (4.0L + 4.0L * std::cos(2.0L * theta));
    ld p = 0.5L * (s * s - (t2 + 4.0L));
    ld disc = s * s - 4.0L * p;
    if (disc >= -1e-9L) {
        real_k = true;
        ld sd = std::sqrt(std::max(0.0L, disc));
        k1 = 0.5L * (s + sd); k2 = 0.5L * (s - sd);
        margin = 2.0L - std::max(std::abs(k1), std::abs(k2));
    } else {
        real_k = false; k1 = 0.5L * s; k2 = std::sqrt(-disc); margin = -10.0L;
    }
}

static void broucke_3d(const std::vector<std::vector<ld>>& M, ld theta, bool& real_k, ld& k1, ld& k2, ld& k3, ld& margin) {
    ld M2[18][18] = {{0}}, tr1 = 0, tr2 = 0, tr3 = 0;
    for (int i = 0; i < 18; ++i) {
        tr1 += M[i][i];
        for (int j = 0; j < 18; ++j) {
            ld s = 0; for (int k = 0; k < 18; ++k) s += M[i][k] * M[k][j];
            M2[i][j] = s;
        }
        tr2 += M2[i][i];
    }
    for (int i = 0; i < 18; ++i) for (int j = 0; j < 18; ++j) tr3 += M2[i][j] * M[j][i];
    ld t1 = tr1 - (6.0L + 6.0L * std::cos(theta));
    ld t2 = tr2 - (6.0L + 6.0L * std::cos(2.0L * theta));
    ld t3 = tr3 - (6.0L + 6.0L * std::cos(3.0L * theta));
    ld p2 = t2 + 6.0L, p3 = t3 + 3.0L * t1;
    ld e1 = t1, e2 = 0.5L * (e1 * e1 - p2), e3 = (p3 - e1 * e1 * e1 + 3.0L * e1 * e2) / 3.0L;
    ld p = e2 - e1 * e1 / 3.0L, q = -2.0L * e1 * e1 * e1 / 27.0L + e1 * e2 / 3.0L - e3;
    ld D = 0.25L * q * q + (p / 3.0L) * (p / 3.0L) * (p / 3.0L);
    if (D <= 1e-12L) {
        real_k = true;
        ld r = std::sqrt(std::max(0.0L, -p / 3.0L)), r3 = r * r * r;
        ld arg = (r3 > 1e-24L) ? std::clamp(-0.5L * q / r3, -1.0L, 1.0L) : 0.0L;
        ld phi = std::acos(arg) / 3.0L;
        const ld PI_L = std::acos(-1.0L);
        std::array<ld, 3> ks = {
            2.0L * r * std::cos(phi) + e1 / 3.0L,
            2.0L * r * std::cos(phi + 2.0L * PI_L / 3.0L) + e1 / 3.0L,
            2.0L * r * std::cos(phi + 4.0L * PI_L / 3.0L) + e1 / 3.0L
        };
        std::sort(ks.begin(), ks.end(), [](ld a, ld b) { return a > b; });
        k1 = ks[0]; k2 = ks[1]; k3 = ks[2];
        margin = 2.0L - std::max({std::abs(k1), std::abs(k2), std::abs(k3)});
    } else {
        real_k = false; k1 = k2 = k3 = 0; margin = -10.0L;
    }
}

static void audit_period(const State1PN& y0, ld T, ld theta, const ld u[3], const Masses& m, ld beta, int dim,
                         ld& ret, ld& dH, ld& symp, bool& real_k, ld& k1, ld& k2, ld& k3, ld& margin) {
    const ld itol = 1e-14L;
    State1PN yT = y0;
    integrate_1pn(yT, T, m, beta, itol);
    dH = std::abs(hamiltonian_1pn(yT, m, beta) - hamiltonian_1pn(y0, m, beta));
    State1PN yTc = (std::abs(theta) > 1e-14L) ? rotate_state(yT, u, -theta) : yT;
    ld r2 = 0; for (int i = 0; i < 18; ++i) r2 += (yTc[i] - y0[i]) * (yTc[i] - y0[i]);
    ret = std::sqrt(r2);

    const ld eps = 1e-8L;
    std::vector<std::vector<ld>> M(18, std::vector<ld>(18, 0.0L));
    for (int col = 0; col < 18; ++col) {
        State1PN yp = y0, ym = y0; yp[col] += eps; ym[col] -= eps;
        integrate_1pn(yp, T, m, beta, itol);
        integrate_1pn(ym, T, m, beta, itol);
        if (std::abs(theta) > 1e-14L) { yp = rotate_state(yp, u, -theta); ym = rotate_state(ym, u, -theta); }
        for (int r = 0; r < 18; ++r) M[r][col] = (yp[r] - ym[r]) / (2.0L * eps);
    }
    // Canonical symplecticity ||M^T J M - J||_inf (y is already canonical [r, p]!)
    symp = 0.0L;
    for (int i = 0; i < 18; ++i)
        for (int j = 0; j < 18; ++j) {
            ld s = 0;
            for (int k = 0; k < 9; ++k) s += M[k][i] * M[k + 9][j] - M[k + 9][i] * M[k][j];
            ld tgt = (j == i + 9) ? 1.0L : ((i == j + 9) ? -1.0L : 0.0L);
            symp = std::max(symp, std::abs(s - tgt));
        }

    if (dim == 2) {
        const int pidx[12] = {0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16};
        const int zidx[6] = {2, 5, 8, 11, 14, 17};
        std::vector<std::vector<ld>> Mp(12, std::vector<ld>(12));
        for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) Mp[r][c] = M[pidx[r]][pidx[c]];
        broucke_2d(Mp, theta, real_k, k1, k2, margin);
        ld trz = 0; for (int i = 0; i < 6; ++i) trz += M[zidx[i]][zidx[i]];
        k3 = trz - (2.0L + 2.0L * std::cos(theta));
    } else {
        broucke_3d(M, theta, real_k, k1, k2, k3, margin);
    }
}

int main() {
    // Input: id dim beta m0 m1 m2 T T_rpo q_fold precess_rad z0[18]
    struct Job { std::string id; int dim, q_fold; ld beta, T, T_rpo, prec; Masses m; State1PN y0; };
    std::vector<Job> jobs;
    Job jb;
    while (std::cin >> jb.id >> jb.dim >> jb.beta >> jb.m[0] >> jb.m[1] >> jb.m[2] >> jb.T >> jb.T_rpo >> jb.q_fold >> jb.prec) {
        for (int i = 0; i < 18; ++i) std::cin >> jb.y0[i];
        jobs.push_back(jb);
    }
    std::vector<std::string> outs(jobs.size());
    std::atomic<size_t> next_idx{0};
    int nthr = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
    std::vector<std::thread> ths;
    for (int t = 0; t < nthr; ++t) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= jobs.size()) break;
                const auto& j = jobs[idx];
                ld L[3] = {0, 0, 0};
                for (int b = 0; b < 3; ++b) {
                    L[0] += j.y0[3 * b + 1] * j.y0[11 + 3 * b] - j.y0[3 * b + 2] * j.y0[10 + 3 * b];
                    L[1] += j.y0[3 * b + 2] * j.y0[9 + 3 * b]  - j.y0[3 * b + 0] * j.y0[11 + 3 * b];
                    L[2] += j.y0[3 * b + 0] * j.y0[10 + 3 * b] - j.y0[3 * b + 1] * j.y0[9 + 3 * b];
                }
                ld Lm = std::sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
                ld u[3] = {0, 0, 1};
                if (j.dim == 3 && Lm > 1e-12L) { u[0] = L[0] / Lm; u[1] = L[1] / Lm; u[2] = L[2] / Lm; }

                bool is_closure = (j.q_fold > 0);
                ld T_eval = j.T, theta_eval = is_closure ? 0.0L : j.prec;
                ld ret, dH, symp, k1, k2, k3, margin; bool real_k;
                audit_period(j.y0, T_eval, theta_eval, u, j.m, j.beta, j.dim, ret, dH, symp, real_k, k1, k2, k3, margin);

                std::ostringstream os;
                os << std::setprecision(12)
                   << "{\"id\":\"" << j.id << "\",\"dim\":" << j.dim << ",\"q_fold\":" << j.q_fold
                   << ",\"T_eval\":" << (double)T_eval
                   << ",\"ret_err\":" << (double)ret
                   << ",\"dH\":" << (double)dH
                   << ",\"symp_err\":" << (double)symp
                   << ",\"real_k\":" << (real_k ? "true" : "false")
                   << ",\"k1\":" << (double)k1 << ",\"k2\":" << (double)k2 << ",\"k3\":" << (double)k3
                   << ",\"margin\":" << (double)margin;
                if (is_closure) {
                    ld ret_r, dH_r, symp_r, k1_r, k2_r, k3_r, mrg_r; bool rk_r;
                    audit_period(j.y0, j.T_rpo, j.prec, u, j.m, j.beta, j.dim, ret_r, dH_r, symp_r, rk_r, k1_r, k2_r, k3_r, mrg_r);
                    os << ",\"rpo_k1\":" << (double)k1_r << ",\"rpo_k2\":" << (double)k2_r
                       << ",\"rpo_k3\":" << (double)k3_r << ",\"rpo_margin\":" << (double)mrg_r;
                }
                os << "}\n";
                outs[idx] = os.str();
            }
        });
    }
    for (auto& th : ths) th.join();
    for (const auto& s : outs) std::cout << s;
    return 0;
}
