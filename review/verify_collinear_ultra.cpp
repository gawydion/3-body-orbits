#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Ultra-precision 80-bit long double analytical-variational monodromy verifier
// for collinear L=0 orbits. Exploits the exact half-period spatial involution
// symmetry P(theta) * y(T/2) = y(0) (P^2 = I, canonical symplectic) when present,
// so M(T) = (P * M(T/2))^2 and k_j = h_j^2 - 2, reducing the monodromy condition
// number by ~100x and achieving sympH ~ 3e-11 .. 1e-8. Falls back to full-period
// analytical variational refinement for orbits without half-period involution
// symmetry (such as 2D-DEM-037).

using ld = long double;
using State = std::array<ld, 12>;
using Masses = std::array<ld, 3>;

static inline void deriv2d(const State& s, const Masses& m, State& ds, ld* rm = nullptr) {
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

static inline void accel_jac2d(const State& s, const Masses& m, ld H[6][6]) {
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

struct Var2D { State y; ld M[12][12]; };

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

static bool integrate_var2d(const State& y0, ld T, const Masses& m, State& yT, ld Mout[12][12], ld tol = 1e-17L, ld mtol_scale = 1.0L) {
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
        // P = Refl_x * Rot(-th)
        for (int b = 0; b < 6; ++b) {
            P[2 * b][2 * b] = c; P[2 * b][2 * b + 1] = -s;
            P[2 * b + 1][2 * b] = -s; P[2 * b + 1][2 * b + 1] = -c;
        }
    } else {
        // P = Swap12 * Refl_y * Rot(-th)
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

static State apply_P(const ld P[12][12], const State& y) {
    State r{}; for (int i = 0; i < 12; ++i) for (int j = 0; j < 12; ++j) r[i] += P[i][j] * y[j]; return r;
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

static State make_col(ld vx, ld vy, ld m3) {
    return {-1, 0, 1, 0, 0, 0, vx, vy, vx, vy, -2 * vx / m3, -2 * vy / m3};
}

struct Job { std::string id; ld m3, vx, vy, T, ck1, ck2; };

int main() {
    std::vector<Job> jobs; Job j;
    while (std::cin >> j.id >> j.m3 >> j.vx >> j.vy >> j.T >> j.ck1 >> j.ck2) jobs.push_back(j);
    std::vector<std::string> outs(jobs.size());
    std::vector<std::thread> ths;
    for (size_t idx = 0; idx < jobs.size(); ++idx) ths.emplace_back([&, idx]() {
        auto jb = jobs[idx];
        ld m3 = jb.m3, vx = jb.vx, vy = jb.vy, T = jb.T;
        Masses m = {1, 1, m3};
        auto dist12 = [](const State& a, const State& b) {
            ld s = 0; for (int i = 0; i < 12; ++i) s += (a[i] - b[i]) * (a[i] - b[i]); return std::sqrt(s);
        };

        // Test for exact half-period spatial involution symmetry
        State y0 = make_col(vx, vy, m3), yH; ld MH[12][12];
        integrate_var2d(y0, T / 2.0L, m, yH, MH, 1e-15L, 10.0L);
        ld th = std::atan2(yH[3], yH[2]);
        ld P0[12][12], P1[12][12];
        build_P(th, 0, P0); build_P(th, 1, P1);
        ld e0 = dist12(apply_P(P0, yH), y0), e1 = dist12(apply_P(P1, yH), y0);

        std::ostringstream os;
        if (std::min(e0, e1) < 1e-3L) {
            int mode = (e0 < e1) ? 0 : 1;
            ld Th = T / 2.0L;
            ld min_dk1 = 1e30L, max_dk1 = -1e30L, last_k1 = 0, last_k2 = 0, last_ret = 0, last_symp = 0;
            int n_conv = 0;
            for (int it = 0; it < 6; ++it) {
                y0 = make_col(vx, vy, m3);
                if (!integrate_var2d(y0, Th, m, yH, MH, 1e-17L, 1.0L)) break;
                ld P[12][12], dP[12][12];
                build_P(th, mode, P);
                // Exact analytical dP/dth
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
                State PyH = apply_P(P, yH);
                State F; ld fn = 0; for (int i = 0; i < 12; ++i) { F[i] = PyH[i] - y0[i]; fn += F[i] * F[i]; }
                ld retH = std::sqrt(fn);
                ld A[12][12] = {{0}};
                for (int i = 0; i < 12; ++i) for (int c_ = 0; c_ < 12; ++c_) for (int k = 0; k < 12; ++k) A[i][c_] += P[i][k] * MH[k][c_];
                ld tr1 = 0, tr2 = 0;
                for (int i = 0; i < 12; ++i) { tr1 += A[i][i]; for (int c_ = 0; c_ < 12; ++c_) tr2 += A[i][c_] * A[c_][i]; }
                ld sA = tr1, sq = tr2 - 4.0L, prod = 0.5L * (sA * sA - sq), disc = sA * sA - 4.0L * prod;
                ld h1 = 0.5L * (sA + std::sqrt(std::max(0.0L, disc)));
                ld h2 = 0.5L * (sA - std::sqrt(std::max(0.0L, disc)));
                ld k1 = h1 * h1 - 2.0L, k2 = h2 * h2 - 2.0L;
                if (std::abs(k2 - 2.0L) < std::abs(k1 - 2.0L)) { std::swap(k1, k2); std::swap(h1, h2); }
                if (it >= 1) {
                    min_dk1 = std::min(min_dk1, k1 - 2.0L);
                    max_dk1 = std::max(max_dk1, k1 - 2.0L);
                    ++n_conv;
                }
                last_k1 = k1; last_k2 = k2; last_ret = retH; last_symp = symp_err_12(A, m);
                if (it == 5 || (it >= 2 && retH < 1e-16L)) break;
                ld dvx[12] = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3, 0};
                ld dvy[12] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3};
                State fH; deriv2d(yH, m, fH);
                State PfH = apply_P(P, fH), dPyH = apply_P(dP, yH);
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
                    int piv = c_; for (int r = c_ + 1; r < 4; ++r) if (std::abs(G[r][c_]) > std::abs(G[piv][c_])) piv = r;
                    for (int k = c_; k < 5; ++k) std::swap(G[c_][k], G[piv][k]);
                    for (int r = c_ + 1; r < 4; ++r) { ld f = G[r][c_] / G[c_][c_]; for (int k = c_; k < 5; ++k) G[r][k] -= f * G[c_][k]; }
                }
                ld d[4] = {0}; for (int r = 3; r >= 0; --r) { ld s_ = G[r][4]; for (int c_ = r + 1; c_ < 4; ++c_) s_ -= G[r][c_] * d[c_]; d[r] = s_ / G[r][r]; }
                vx += d[0]; vy += d[1]; Th += d[2]; th += d[3];
            }
            ld mid_dk1 = 0.5L * (min_dk1 + max_dk1);
            ld half_span = 0.5L * (max_dk1 - min_dk1);
            std::string verdict;
            if (min_dk1 > 0.0L) verdict = "UNSTABLE";
            else if (max_dk1 >= 0.0L) verdict = "NEAR-PARABOLIC";
            else if (-mid_dk1 < 1e-4L) verdict = "SUB<1e-4";
            else verdict = "YES";
            os << jb.id << " [half m=" << mode << "] m3=" << std::fixed << std::setprecision(2) << (double)m3
               << " retH=" << std::scientific << std::setprecision(2) << (double)last_ret
               << " sympH=" << (double)last_symp
               << " k1-2=" << std::showpos << std::setprecision(3) << (double)mid_dk1
               << std::noshowpos << "+/-" << std::setprecision(1) << (double)half_span
               << std::fixed << std::setprecision(7)
               << " K=(" << (double)(2.0L + mid_dk1) << ", " << (double)last_k2 << ")"
               << " VERDICT=" << verdict
               << " CAT=(" << std::setprecision(4) << (double)jb.ck1 << ", " << (double)jb.ck2 << ")";
        } else {
            // Full-period refinement fallback (e.g., 2D-DEM-037)
            State yT; ld M[12][12];
            for (int it = 0; it < 6; ++it) {
                y0 = make_col(vx, vy, m3);
                if (!integrate_var2d(y0, T, m, yT, M, 1e-16L, 10.0L)) break;
                State F; ld fn = 0; for (int i = 0; i < 12; ++i) { F[i] = yT[i] - y0[i]; fn += F[i] * F[i]; }
                if (std::sqrt(fn) < 1e-15L) break;
                ld dvx[12] = {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3, 0};
                ld dvy[12] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -2.0L / m3};
                State fT; deriv2d(yT, m, fT);
                ld J[12][3];
                for (int r = 0; r < 12; ++r) { ld s0 = -dvx[r], s1 = -dvy[r]; for (int k = 0; k < 12; ++k) { s0 += M[r][k] * dvx[k]; s1 += M[r][k] * dvy[k]; } J[r][0] = s0; J[r][1] = s1; J[r][2] = fT[r]; }
                ld A[3][4] = {{0}};
                for (int r = 0; r < 3; ++r) { for (int c = 0; c < 3; ++c) for (int i = 0; i < 12; ++i) A[r][c] += J[i][r] * J[i][c]; for (int i = 0; i < 12; ++i) A[r][3] -= J[i][r] * F[i]; A[r][r] += 1e-22L; }
                for (int c = 0; c < 3; ++c) { int piv = c; for (int r = c + 1; r < 3; ++r) if (std::abs(A[r][c]) > std::abs(A[piv][c])) piv = r; for (int k = c; k < 4; ++k) std::swap(A[c][k], A[piv][k]); for (int r = c + 1; r < 3; ++r) { ld f = A[r][c] / A[c][c]; for (int k = c; k < 4; ++k) A[r][k] -= f * A[c][k]; } }
                ld d[3] = {0}; for (int r = 2; r >= 0; --r) { ld s = A[r][3]; for (int c = r + 1; c < 3; ++c) s -= A[r][c] * d[c]; d[r] = s / A[r][r]; }
                ld sn = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]); if (sn > 1e-5L) for (int i = 0; i < 3; ++i) d[i] *= (1e-5L / sn);
                vx += d[0]; vy += d[1]; T += d[2];
            }
            y0 = make_col(vx, vy, m3);
            integrate_var2d(y0, T, m, yT, M, 1e-16L, 10.0L);
            ld ret = dist12(yT, y0), symp = symp_err_12(M, m);
            ld tr1 = 0, tr2 = 0; for (int i = 0; i < 12; ++i) tr1 += M[i][i]; for (int i = 0; i < 12; ++i) for (int c = 0; c < 12; ++c) tr2 += M[i][c] * M[c][i];
            ld s = tr1 - 8.0L, p = 0.5L * (s * s - (tr2 - 4.0L)), disc = s * s - 4.0L * p;
            ld k1 = 0.5L * (s + std::sqrt(std::max(0.0L, disc))), k2 = 0.5L * (s - std::sqrt(std::max(0.0L, disc)));
            bool st = (std::abs(k1) < 2.0L) && (std::abs(k2) < 2.0L);
            ld marg = 2.0L - std::max(std::abs(k1), std::abs(k2));
            os << jb.id << " [full    ] m3=" << std::fixed << std::setprecision(2) << (double)m3
               << " ret =" << std::scientific << std::setprecision(2) << (double)ret
               << " symp =" << (double)symp
               << " k1-2=" << std::showpos << std::setprecision(3) << (double)(k1 - 2.0L)
               << std::noshowpos << "        " << std::fixed << std::setprecision(7)
               << " K=(" << (double)k1 << ", " << (double)k2 << ")"
               << " VERDICT=" << (st ? (marg > 1e-4L ? "YES" : "SUB<1e-4") : "UNSTABLE")
               << " CAT=(" << std::setprecision(4) << (double)jb.ck1 << ", " << (double)jb.ck2 << ")";
        }
        outs[idx] = os.str();
    });
    for (auto& th : ths) th.join();
    for (auto& s : outs) std::cout << s << "\n";
    return 0;
}
