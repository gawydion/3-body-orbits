#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using State = std::array<double, 12>;
using Masses = std::array<double, 3>;

static inline void deriv(const State& s, const Masses& m, State& ds, double* rmin_out = nullptr) {
    ds[0] = s[6];  ds[1] = s[7];
    ds[2] = s[8];  ds[3] = s[9];
    ds[4] = s[10]; ds[5] = s[11];

    double ax[3] = {0.0, 0.0, 0.0};
    double ay[3] = {0.0, 0.0, 0.0};
    double rmin2 = 1e300;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[2 * j] - s[2 * i];
            double dy = s[2 * j + 1] - s[2 * i + 1];
            double r2 = dx * dx + dy * dy;
            if (r2 < rmin2) rmin2 = r2;
            double inv_r = 1.0 / std::sqrt(r2);
            double inv_r3 = inv_r * inv_r * inv_r;
            double fx = dx * inv_r3, fy = dy * inv_r3;
            ax[i] += m[j] * fx; ay[i] += m[j] * fy;
            ax[j] -= m[i] * fx; ay[j] -= m[i] * fy;
        }
    }
    ds[6] = ax[0]; ds[7] = ay[0];
    ds[8] = ax[1]; ds[9] = ay[1];
    ds[10] = ax[2]; ds[11] = ay[2];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double compute_energy(const State& s, const Masses& m) {
    double K = 0.0, U = 0.0;
    for (int i = 0; i < 3; ++i)
        K += 0.5 * m[i] * (s[6 + 2 * i] * s[6 + 2 * i] + s[7 + 2 * i] * s[7 + 2 * i]);
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            U -= m[i] * m[j] / std::hypot(s[2 * j] - s[2 * i], s[2 * j + 1] - s[2 * i + 1]);
    return K + U;
}

static inline double compute_L(const State& s, const Masses& m) {
    double L = 0.0;
    for (int i = 0; i < 3; ++i)
        L += m[i] * (s[2 * i] * s[7 + 2 * i] - s[2 * i + 1] * s[6 + 2 * i]);
    return L;
}

static inline double state_dist(const State& a, const State& b) {
    double sum = 0.0;
    for (int i = 0; i < 12; ++i) { double d = a[i] - b[i]; sum += d * d; }
    return std::sqrt(sum);
}

static bool dopri5_step(State& y, double& t, double& h, double t_end, const Masses& m,
                        double tol, double& rmin_track, double r_abort = 1e-6) {
    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0,
                            a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                            a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                            a65 = -5103.0 / 18656.0;
    static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0,
                            b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0,
                            b6 = 11.0 / 84.0;
    static constexpr double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0,
                            e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0,
                            e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;

    State k1, k2, k3, k4, k5, k6, k7, yt, ynew;
    double rmin_local = 1e300;
    for (int attempt = 0; attempt < 25; ++attempt) {
        deriv(y, m, k1, &rmin_local);
        if (rmin_local < r_abort) return false;
        for (int i = 0; i < 12; ++i) yt[i] = y[i] + h * a21 * k1[i];
        deriv(yt, m, k2);
        for (int i = 0; i < 12; ++i) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        deriv(yt, m, k3);
        for (int i = 0; i < 12; ++i) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        deriv(yt, m, k4);
        for (int i = 0; i < 12; ++i) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        deriv(yt, m, k5);
        for (int i = 0; i < 12; ++i) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        deriv(yt, m, k6);
        for (int i = 0; i < 12; ++i) ynew[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        deriv(ynew, m, k7, &rmin_local);
        if (rmin_local < r_abort) return false;

        double err = 0.0;
        for (int i = 0; i < 12; ++i) {
            double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            double sc = tol * (1.0 + std::max(std::abs(y[i]), std::abs(ynew[i])));
            double ratio = ei / sc;
            err += ratio * ratio;
        }
        err = std::sqrt(err / 12.0);
        if (err <= 1.0) {
            y = ynew;
            t += h;
            if (rmin_local < rmin_track) rmin_track = rmin_local;
            double fac = (err == 0.0) ? 5.0 : std::min(5.0, std::max(0.2, 0.9 * std::pow(err, -0.2)));
            h = std::min(h * fac, 0.05);
            return true;
        }
        h *= std::max(0.1, 0.9 * std::pow(err, -0.25));
        if (h < 1e-14) return false;
    }
    return false;
}

static bool integrate_to(State& y, double T, const Masses& m, double tol = 1e-12,
                         int max_steps = 250000, double* rmin_out = nullptr, double r_abort = 1e-6) {
    double t = 0.0, h = std::min(1e-3, T * 0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i = 0; i < 6; ++i) if (std::abs(y[i]) > 60.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

// Agekyan-Anosova free-fall initial state:
// r1 = (-0.5, 0), r2 = (0.5, 0), r3 = (x3, y3), shifted to CM = (0,0), v = 0.
static State make_brake(double x3, double y3, const Masses& m) {
    double M = m[0] + m[1] + m[2];
    double xcm = (-0.5 * m[0] + 0.5 * m[1] + m[2] * x3) / M;
    double ycm = (m[2] * y3) / M;
    return {
        -0.5 - xcm, -ycm,
         0.5 - xcm, -ycm,
         x3  - xcm,  y3 - ycm,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
}

static inline double vel_norm(const State& s) {
    double v2 = 0.0;
    for (int i = 6; i < 12; ++i) v2 += s[i] * s[i];
    return std::sqrt(v2);
}

// Refine (x3, y3, Th) so all velocities vanish at t = Th: F = (vx1, vy1, vx2, vy2, vx3, vy3) = 0
static bool refine_brake(double& x3, double& y3, double& Th, const Masses& m, double& final_verr,
                         double r_abort = 5e-4) {
    double lambda = 1e-4;
    auto eval_F = [&](double cx3, double cy3, double cTh, std::array<double, 6>& F, State& yend,
                      double tol, int max_st) -> bool {
        if (cTh < 0.2 || cTh > 45.0 || cy3 < 1e-4) return false;
        double d13 = std::hypot(cx3 + 0.5, cy3);
        double d23 = std::hypot(cx3 - 0.5, cy3);
        if (std::min(d13, d23) < 0.04) return false;
        State y = make_brake(cx3, cy3, m);
        double rm = 1e300;
        if (!integrate_to(y, cTh, m, tol, max_st, &rm, r_abort)) return false;
        yend = y;
        for (int i = 0; i < 6; ++i) F[i] = y[6 + i];
        return true;
    };

    std::array<double, 6> F;
    State yend;
    double init_tol = (r_abort < 1e-5) ? 1e-13 : 1e-11;
    int init_st = (r_abort < 1e-5) ? 200000 : 30000;
    if (!eval_F(x3, y3, Th, F, yend, init_tol, init_st)) return false;
    double cost = 0.0;
    for (double v : F) cost += v * v;

    for (int iter = 0; iter < 30; ++iter) {
        double sc = std::sqrt(cost);
        if (sc < 1e-11) break;
        if (iter == 7 && sc > 0.2) return false;
        double cur_tol = (r_abort < 1e-5 || sc < 1e-6) ? 1e-13 : 1e-11;
        int cur_st = (r_abort < 1e-5 || sc < 1e-6) ? 200000 : 30000;
        double eps = std::max(1e-8, std::min(1e-6, sc * 1e-2));
        double J[6][3];
        std::array<double, 6> Fp, Fm;
        State dummy;
        if (!eval_F(x3 + eps, y3, Th, Fp, dummy, cur_tol, cur_st) ||
            !eval_F(x3 - eps, y3, Th, Fm, dummy, cur_tol, cur_st)) return false;
        for (int i = 0; i < 6; ++i) J[i][0] = (Fp[i] - Fm[i]) / (2.0 * eps);
        if (!eval_F(x3, y3 + eps, Th, Fp, dummy, cur_tol, cur_st) ||
            !eval_F(x3, y3 - eps, Th, Fm, dummy, cur_tol, cur_st)) return false;
        for (int i = 0; i < 6; ++i) J[i][1] = (Fp[i] - Fm[i]) / (2.0 * eps);
        State dy;
        deriv(yend, m, dy);
        for (int i = 0; i < 6; ++i) J[i][2] = dy[6 + i];

        double A[3][4] = {{0}};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int i = 0; i < 6; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0.0;
            for (int i = 0; i < 6; ++i) s -= J[i][r] * F[i];
            A[r][3] = s;
        }

        bool sing = false;
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col + 1; r < 3; ++r)
                if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (std::abs(A[piv][col]) < 1e-18) { sing = true; break; }
            for (int c = col; c < 4; ++c) std::swap(A[col][c], A[piv][c]);
            for (int r = col + 1; r < 3; ++r) {
                double fac = A[r][col] / A[col][col];
                for (int c = col; c < 4; ++c) A[r][c] -= fac * A[col][c];
            }
        }
        if (sing) { lambda *= 10.0; continue; }
        double d[3] = {0};
        for (int r = 2; r >= 0; --r) {
            double s = A[r][3];
            for (int c = r + 1; c < 3; ++c) s -= A[r][c] * d[c];
            d[r] = s / A[r][r];
        }
        // Keep spatial step bounded inside the local AA basin so LM never leaps across collision curves!
        double sxy = std::hypot(d[0], d[1]);
        double scale = 1.0;
        if (sxy > 0.025) scale = std::min(scale, 0.025 / sxy);
        if (std::abs(d[2]) > 0.15) scale = std::min(scale, 0.15 / std::abs(d[2]));
        if (scale < 1.0) for (int i = 0; i < 3; ++i) d[i] *= scale;

        std::array<double, 6> nF;
        State nyend;
        if (!eval_F(x3 + d[0], y3 + d[1], Th + d[2], nF, nyend, cur_tol, cur_st)) { lambda *= 10.0; continue; }
        double ncost = 0.0;
        for (double v : nF) ncost += v * v;
        if (ncost < cost) {
            x3 += d[0]; y3 += d[1]; Th += d[2];
            F = nF; yend = nyend; cost = ncost;
            lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0;
            if (lambda > 1e6) break;
        }
    }
    final_verr = std::sqrt(cost);
    return (final_verr < 1e-9 && Th > 0.2);
}

struct TopoInfo {
    int partner_swaps;
    int syzygies;
    double rmin;
    double rmax;
    bool democratic;
};

static TopoInfo analyze_topology(const State& y0, double T, const Masses& m) {
    const int N = 400;
    double dt = T / N;
    State y = y0;
    int prev_closest = -1, swaps = 0, syz = 0;
    double prev_area = 0.0, rmin = 1e300, rmax = 0.0;
    for (int k = 0; k <= N; ++k) {
        double r12 = std::hypot(y[2] - y[0], y[3] - y[1]);
        double r13 = std::hypot(y[4] - y[0], y[5] - y[1]);
        double r23 = std::hypot(y[4] - y[2], y[5] - y[3]);
        rmin = std::min({rmin, r12, r13, r23});
        rmax = std::max({rmax, r12, r13, r23});
        int cl = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
        if (prev_closest != -1 && cl != prev_closest) swaps++;
        prev_closest = cl;
        double area = (y[2] - y[0]) * (y[5] - y[1]) - (y[3] - y[1]) * (y[4] - y[0]);
        if (k > 0 && ((prev_area > 0 && area < 0) || (prev_area < 0 && area > 0))) syz++;
        prev_area = area;
        if (k < N) integrate_to(y, dt, m, 1e-11, 80000, nullptr, 1e-6);
    }
    return {swaps, syz, rmin, rmax, swaps >= 2};
}

struct StabilityInfo {
    bool real_k;
    double k1, k2;
    bool linearly_stable;
    double stab_margin;
    int survived_100T;
    double drift_100T;
};

static StabilityInfo analyze_stability(const State& y0, double T, const Masses& m) {
    double M[12][12];
    const double eps = 1e-7;
    StabilityInfo info{};
    info.stab_margin = -999.0;
    for (int col = 0; col < 12; ++col) {
        State yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        if (!integrate_to(yp, T, m, 1e-13, 400000, nullptr, 1e-6) ||
            !integrate_to(ym, T, m, 1e-13, 400000, nullptr, 1e-6)) {
            return info;
        }
        for (int row = 0; row < 12; ++row)
            M[row][col] = (yp[row] - ym[row]) / (2.0 * eps);
    }
    double tr1 = 0.0, tr2 = 0.0;
    for (int i = 0; i < 12; ++i) tr1 += M[i][i];
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j)
            tr2 += M[i][j] * M[j][i];

    double s = tr1 - 8.0;
    double p = 0.5 * (s * s - (tr2 - 4.0));
    double disc = s * s - 4.0 * p;

    if (disc >= -1e-4) {
        double sd = std::sqrt(std::max(0.0, disc));
        info.real_k = true;
        info.k1 = 0.5 * (s + sd);
        info.k2 = 0.5 * (s - sd);
        double mk = std::max(std::abs(info.k1), std::abs(info.k2));
        info.stab_margin = 2.0 - mk;
        info.linearly_stable = (mk < 2.0 - 1e-4);
    } else {
        info.real_k = false;
        info.k1 = 0.5 * s;
        info.k2 = 0.5 * std::sqrt(-disc);
        info.linearly_stable = false;
        info.stab_margin = -10.0;
    }

    if (info.linearly_stable) {
        State yp = y0;
        yp[6] += 1e-5;
        yp[8] -= 1e-5 * m[0] / m[1];
        yp[7] += 1e-5;
        yp[9] -= 1e-5 * m[0] / m[1];
        double max_drift = 0.0;
        int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate_to(yp, T, m, 1e-11, 250000, nullptr, 1e-6)) break;
            double d = state_dist(yp, y0);
            if (d > max_drift) max_drift = d;
            if (d > 5.0) break;
            surv = k + 1;
        }
        info.survived_100T = surv;
        info.drift_100T = max_drift;
    }
    return info;
}

static void emit_json(const Masses& m, double x3, double y3, const State& y0, double T,
                      double ret_err, const StabilityInfo& st, const TopoInfo& tp) {
    double E = compute_energy(y0, m);
    double L = compute_L(y0, m);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"dim\":2,\"subspace\":\"free_fall_brake\""
       << ",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "]"
       << ",\"aa_x3\":" << x3 << ",\"aa_y3\":" << y3
       << ",\"y0\":[";
    for (int i = 0; i < 12; ++i) os << y0[i] << (i < 11 ? "," : "");
    os << "],\"T\":" << T
       << ",\"E\":" << E
       << ",\"L\":" << L
       << ",\"ret_err\":" << ret_err
       << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2
       << ",\"stab_margin\":" << st.stab_margin
       << ",\"survived_100T\":" << st.survived_100T
       << ",\"drift_100T\":" << st.drift_100T
       << ",\"rmin\":" << tp.rmin
       << ",\"rmax\":" << tp.rmax
       << ",\"partner_swaps\":" << tp.partner_swaps
       << ",\"syzygies\":" << tp.syzygies
       << ",\"democratic\":" << (tp.democratic ? "true" : "false")
       << "}\n";
    std::cout << os.str() << std::flush;
}

static void run_scan(const Masses& m, int nx, int ny, double Thmax,
                     double xmin, double xmax, double ymin, double ymax, bool strict_aa, bool log_near) {
    struct Cell { double x3, y3, Th, vn; };
    struct Seed { double x3, y3, Th, vn; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    // 3 time bands per grid point so short, medium, and longer brake orbits never mask each other!
    const double t_cuts[4] = {0.25, std::min(4.0, Thmax * 0.35), std::min(9.0, Thmax * 0.7), Thmax};
    std::vector<std::vector<std::array<Cell, 3>>> grid(nx, std::vector<std::array<Cell, 3>>(ny));
    for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy)
            for (int b = 0; b < 3; ++b)
                grid[ix][iy][b] = {0.0, 0.0, 0.0, 1e9};

    std::vector<std::thread> ths;
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            for (int ix = tid; ix < nx; ix += nthreads) {
                double x3 = (nx > 1) ? (xmin + (xmax - xmin) * ix / (nx - 1)) : xmin;
                for (int iy = 0; iy < ny; ++iy) {
                    double y3 = (ny > 1) ? (ymin + (ymax - ymin) * iy / (ny - 1)) : ymin;
                    if (strict_aa && (x3 + 0.5) * (x3 + 0.5) + y3 * y3 > 1.0001) continue;
                    double d13 = std::hypot(x3 + 0.5, y3);
                    double d23 = std::hypot(x3 - 0.5, y3);
                    if (std::min(d13, d23) < 0.06) continue;

                    State y = make_brake(x3, y3, m);
                    double Eabs = std::max(0.1, std::abs(compute_energy(y, m)));
                    double v_scale = std::sqrt(2.0 * Eabs / (m[0] + m[1] + m[2]));
                    double t = 0.0, h = 1e-3, rmin = 1e300;
                    double pv = 1e9, p2v = 1e9, pt = 0.0;
                    bool accel = false;
                    int steps = 0;
                    while (t < Thmax && steps < 25000) {
                        double hc = std::min(h, 0.02);
                        if (!dopri5_step(y, t, hc, Thmax, m, 1e-8, rmin, 1e-4)) break;
                        h = hc; ++steps;
                        double maxr = 0;
                        for (int i = 0; i < 6; ++i) maxr = std::max(maxr, std::abs(y[i]));
                        if (maxr > 22.0) break;
                        double K = 0.0;
                        for (int b = 0; b < 3; ++b)
                            K += 0.5 * m[b] * (y[6 + 2*b]*y[6 + 2*b] + y[7 + 2*b]*y[7 + 2*b]);
                        double v = std::sqrt(K / Eabs);
                        if (v > 0.25) accel = true;
                        if (accel && pt > t_cuts[0] && steps >= 3 && pv < p2v && pv < v) {
                            int band = (pt < t_cuts[1]) ? 0 : ((pt < t_cuts[2]) ? 1 : 2);
                            if (pv < grid[ix][iy][band].vn) {
                                grid[ix][iy][band] = {x3, y3, pt, pv};
                            }
                        }
                        p2v = pv; pv = v; pt = t;
                    }
                    (void)v_scale;
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    // Extract 2D local minima in each time band!
    std::vector<Seed> filt;
    for (int b = 0; b < 3; ++b) {
        for (int ix = 1; ix < nx - 1; ++ix) {
            for (int iy = 1; iy < ny - 1; ++iy) {
                double v = grid[ix][iy][b].vn;
                if (v > 0.55 || grid[ix][iy][b].Th < 0.25) continue;
                bool is_min = true;
                for (int dx = -1; dx <= 1 && is_min; ++dx)
                    for (int dy = -1; dy <= 1 && is_min; ++dy)
                        if ((dx || dy) && grid[ix + dx][iy + dy][b].vn < v) is_min = false;
                if (is_min) {
                    filt.push_back({grid[ix][iy][b].x3, grid[ix][iy][b].y3, grid[ix][iy][b].Th, v});
                }
            }
        }
    }
    std::sort(filt.begin(), filt.end(), [](const Seed& a, const Seed& b) { return a.vn < b.vn; });
    if (filt.size() > 1200) filt.resize(1200);
    std::cerr << "[brake m=(" << m[0] << "," << m[1] << "," << m[2] << ")] local_min_seeds=" << filt.size() << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 3>> dedup;
    std::atomic<size_t> next_idx{0};
    int n_closed = 0, n_stable = 0;
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= filt.size()) break;
                double x3 = filt[idx].x3, y3 = filt[idx].y3, Th = filt[idx].Th, verr = 0.0;
                if (!refine_brake(x3, y3, Th, m, verr)) continue;

                double T = 2.0 * Th;
                State y0 = make_brake(x3, y3, m);
                // Check if sub-multiple is already periodic
                for (int div : {2, 3}) {
                    State yt = y0;
                    if (integrate_to(yt, T / div, m, 1e-12, 180000, nullptr, 1e-5) &&
                        state_dist(yt, y0) < 1e-5) {
                        Th /= div;
                        T = 2.0 * Th;
                        refine_brake(x3, y3, Th, m, verr);
                        y0 = make_brake(x3, y3, m);
                        break;
                    }
                }

                State yT = y0;
                double rmin_full = 1e300;
                if (!integrate_to(yT, T, m, 1e-13, 350000, &rmin_full, 1e-5)) continue;
                double ret = state_dist(yT, y0);
                if (ret > 1e-7) continue;

                double E = compute_energy(y0, m);
                double inv1 = E * std::pow(T, 2.0 / 3.0);
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup) {
                        if (std::abs(f[0] - inv1) < 1e-3 &&
                            std::hypot(f[1] - x3, f[2] - y3) < 2e-2) { dup = true; break; }
                    }
                    if (dup) continue;
                    dedup.push_back({inv1, x3, y3});
                    n_closed++;
                }

                auto st = analyze_stability(y0, T, m);
                auto tp = analyze_topology(y0, T, m);
                if (st.linearly_stable && st.survived_100T >= 100) {
                    std::lock_guard<std::mutex> lk(mtx);
                    n_stable++;
                    emit_json(m, x3, y3, y0, T, ret, st, tp);
                } else if (log_near && st.real_k && std::max(std::abs(st.k1), std::abs(st.k2)) < 25.0) {
                    std::lock_guard<std::mutex> lk(mtx);
                    std::cerr << std::setprecision(10)
                              << "  [near-brake m=(" << m[0] << "," << m[1] << "," << m[2]
                              << ")] x3=" << x3 << " y3=" << y3 << " Th=" << Th
                              << " k=(" << st.k1 << "," << st.k2 << ") rmin=" << tp.rmin
                              << " swaps=" << tp.partner_swaps << "\n";
                }
            }
        });
    }
    for (auto& th : ths) th.join();
    std::cerr << "[brake m=(" << m[0] << "," << m[1] << "," << m[2] << ")] closed=" << n_closed
              << " linearly_stable=" << n_stable << "\n";
}

static void run_continue(Masses m, double x3, double y3, double Th, int mass_idx, double dm, int nsteps, bool emit_stable_json) {
    std::cout << std::setprecision(12);
    for (int step = 0; step <= nsteps; ++step) {
        Masses cm = m;
        cm[mass_idx] = m[mass_idx] + step * dm;
        if (cm[mass_idx] <= 0.01) break;
        double verr = 0.0;
        if (!refine_brake(x3, y3, Th, cm, verr, 1e-6)) break;
        double T = 2.0 * Th;
        State y0 = make_brake(x3, y3, cm);
        State yT = y0;
        if (!integrate_to(yT, T, cm, 1e-13, 400000, nullptr, 1e-6)) break;
        double ret = state_dist(yT, y0);
        if (ret > 1e-7) break;
        auto st = analyze_stability(y0, T, cm);
        auto tp = analyze_topology(y0, T, cm);
        if (emit_stable_json) {
            if (st.linearly_stable && st.survived_100T >= 100) {
                emit_json(cm, x3, y3, y0, T, ret, st, tp);
            }
            std::cerr << std::setprecision(10)
                      << "cont m=(" << cm[0] << "," << cm[1] << "," << cm[2] << ") x3=" << x3 << " y3=" << y3 << " T=" << T
                      << " k=(" << st.k1 << "," << st.k2 << ") stab=" << st.linearly_stable
                      << " margin=" << st.stab_margin << " rmin=" << tp.rmin
                      << " swaps=" << tp.partner_swaps << "\n";
        } else {
            std::cout << cm[mass_idx] << " " << x3 << " " << y3 << " " << T << " "
                      << st.k1 << " " << st.k2 << " " << (st.linearly_stable ? 1 : 0) << " "
                      << st.stab_margin << " " << tp.rmin << " " << tp.partner_swaps << "\n" << std::flush;
        }
    }
}

static void run_hunt(int nsamples, double Thmax, int seed_off) {
    struct Cand { Masses m; double x3, y3, Th, vn; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Cand>> tcands(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(441001 + seed_off * 10007 + tid * 9973);
            std::uniform_real_distribution<double> u01(0.0, 1.0);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                Masses m = {1.0, 1.0, 1.0};
                double mode_r = u01(rng);
                if (mode_r < 0.45) {
                    // Symmetric unequal m1=m2=1, log-uniform m3 in [0.05, 6.0]
                    m[2] = std::exp(std::log(0.05) + u01(rng) * std::log(6.0 / 0.05));
                } else if (mode_r < 0.85) {
                    // All-distinct masses
                    m[1] = std::exp(std::log(0.15) + u01(rng) * std::log(3.0 / 0.15));
                    m[2] = std::exp(std::log(0.05) + u01(rng) * std::log(4.0 / 0.05));
                } // else equal masses (1,1,1)

                bool sym12 = (std::abs(m[0] - m[1]) < 1e-12);
                double x3 = sym12 ? (u01(rng) * 1.1) : (-0.9 + u01(rng) * 1.8);
                double y3 = 0.04 + u01(rng) * 1.05;
                if (std::min(std::hypot(x3 + 0.5, y3), std::hypot(x3 - 0.5, y3)) < 0.08) continue;

                State y = make_brake(x3, y3, m);
                double Eabs = std::max(0.1, std::abs(compute_energy(y, m)));
                double t = 0.0, h = 1e-3, rm = 1e300;
                double pv = 1e9, p2v = 1e9, pt = 0.0, best_v = 1e9, best_t = 0.0;
                bool accel = false;
                int steps = 0;
                while (t < Thmax && steps < 12000) {
                    double hc = std::min(h, 0.02);
                    if (!dopri5_step(y, t, hc, Thmax, m, 1e-8, rm, 3e-4)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i=0;i<6;++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 18.0) break;
                    double K = 0.0;
                    for (int b=0;b<3;++b) K += 0.5*m[b]*(y[6+2*b]*y[6+2*b]+y[7+2*b]*y[7+2*b]);
                    double v = std::sqrt(K / Eabs);
                    if (v > 0.25) accel = true;
                    if (accel && pt > 0.25 && steps >= 3 && pv < p2v && pv < v && pv < best_v) {
                        best_v = pv; best_t = pt;
                    }
                    p2v = pv; pv = v; pt = t;
                }
                if (best_v < 0.18 && best_t > 0.25) {
                    tcands[tid].push_back({m, x3, y3, best_t, best_v});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Cand> all;
    for (auto& v : tcands) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), [](const Cand& a, const Cand& b){ return a.vn < b.vn; });
    size_t limit = std::min(all.size(), (size_t)1200);
    std::cerr << "[brake_hunt] samples=" << nsamples << " candidates=" << all.size()
              << " refining=" << limit << "\n";

    std::mutex mtx;
    std::atomic<size_t> next_idx{0};
    int n_closed = 0, n_stable = 0;
    double best_mk = 1e9;
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                Masses m = all[idx].m;
                double x3 = all[idx].x3, y3 = all[idx].y3, Th = all[idx].Th, verr = 0.0;
                if (!refine_brake(x3, y3, Th, m, verr)) continue;
                double T = 2.0 * Th;
                State y0 = make_brake(x3, y3, m), yT = y0;
                if (!integrate_to(yT, T, m, 1e-13, 350000, nullptr, 1e-5)) continue;
                double ret = state_dist(yT, y0);
                if (ret > 1e-7) continue;

                auto st = analyze_stability(y0, T, m);
                auto tp = analyze_topology(y0, T, m);
                double mk = st.real_k ? std::max(std::abs(st.k1), std::abs(st.k2)) : 1e9;
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    n_closed++;
                    if (mk < best_mk) best_mk = mk;
                }
                if (st.linearly_stable && st.survived_100T >= 100) {
                    std::lock_guard<std::mutex> lk(mtx);
                    n_stable++;
                    emit_json(m, x3, y3, y0, T, ret, st, tp);
                } else if (mk < 10.0) {
                    {
                        std::lock_guard<std::mutex> lk(mtx);
                        std::cerr << std::setprecision(9)
                                  << "  [hot-brake m=(" << m[0] << "," << m[1] << "," << m[2]
                                  << ")] x3=" << x3 << " y3=" << y3 << " Th=" << Th
                                  << " k=(" << st.k1 << "," << st.k2 << ") rmin=" << tp.rmin
                                  << " swaps=" << tp.partner_swaps << "\n";
                    }
                    // Automatic local mass continuation along m3 and m2 to dip into stability!
                    for (int axis : {2, 1}) {
                        for (double sgn : {-1.0, 1.0}) {
                            Masses cm = m;
                            double cx3 = x3, cy3 = y3, cTh = Th;
                            double step = 0.02 * sgn * std::max(0.2, m[axis]);
                            for (int stp = 1; stp <= 25; ++stp) {
                                cm[axis] = m[axis] + stp * step;
                                if (cm[axis] <= 0.01) break;
                                double ve = 0;
                                if (!refine_brake(cx3, cy3, cTh, cm, ve, 1e-6)) break;
                                double cT = 2.0 * cTh;
                                State cy0 = make_brake(cx3, cy3, cm), cyT = cy0;
                                if (!integrate_to(cyT, cT, cm, 1e-13, 350000, nullptr, 1e-6)) break;
                                double cret = state_dist(cyT, cy0);
                                if (cret > 1e-7) break;
                                auto cst = analyze_stability(cy0, cT, cm);
                                if (cst.linearly_stable && cst.survived_100T >= 100) {
                                    auto ctp = analyze_topology(cy0, cT, cm);
                                    std::lock_guard<std::mutex> lk(mtx);
                                    n_stable++;
                                    emit_json(cm, cx3, cy3, cy0, cT, cret, cst, ctp);
                                }
                                double cmk = cst.real_k ? std::max(std::abs(cst.k1), std::abs(cst.k2)) : 1e9;
                                if (cmk < 4.5) {
                                    std::lock_guard<std::mutex> lk(mtx);
                                    std::cerr << std::setprecision(9)
                                              << "    [VERY-HOT m=(" << cm[0] << "," << cm[1] << "," << cm[2]
                                              << ")] x3=" << cx3 << " y3=" << cy3 << " Th=" << cTh
                                              << " k=(" << cst.k1 << "," << cst.k2 << ") stab=" << cst.linearly_stable << "\n";
                                }
                                if (cmk > mk + 8.0) break;
                            }
                        }
                    }
                }
            }
        });
    }
    for (auto& th : ths) th.join();
    std::cerr << "[brake_hunt] closed=" << n_closed << " best_max_|k|=" << best_mk
              << " linearly_stable=" << n_stable << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  search_brake --scan <m1> <m2> <m3> <nx> <ny> <Thmax> [xmin xmax ymin ymax] [--near]\n"
                  << "  search_brake --continue <m1> <m2> <m3> <x3> <y3> <Th> <dm3> <nsteps> [--json]\n"
                  << "  search_brake --continue-m2 <m1> <m2> <m3> <x3> <y3> <Th> <dm2> <nsteps> [--json]\n"
                  << "  search_brake --hunt <nsamples> <Thmax> [seed]\n";
        return 1;
    }
    std::string mode = argv[1];
    if (mode == "--scan") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        int nx = (argc > 5) ? std::atoi(argv[5]) : 80;
        int ny = (argc > 6) ? std::atoi(argv[6]) : 80;
        double Thmax = (argc > 7) ? std::atof(argv[7]) : 10.0;
        bool equal_all = (std::abs(m[0] - m[1]) < 1e-12 && std::abs(m[1] - m[2]) < 1e-12);
        bool sym12 = (std::abs(m[0] - m[1]) < 1e-12);
        double xmin = (argc > 8) ? std::atof(argv[8]) : (sym12 ? 0.0 : -0.9);
        double xmax = (argc > 9) ? std::atof(argv[9]) : (equal_all ? 0.5 : 1.4);
        double ymin = (argc > 10) ? std::atof(argv[10]) : 0.03;
        double ymax = (argc > 11) ? std::atof(argv[11]) : (equal_all ? 0.87 : 1.4);
        bool log_near = false;
        for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--near") log_near = true;
        run_scan(m, nx, ny, Thmax, xmin, xmax, ymin, ymax, equal_all && argc <= 8, log_near);
    } else if (mode == "--continue" || mode == "--continue-m2") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double x3 = std::atof(argv[5]), y3 = std::atof(argv[6]), Th = std::atof(argv[7]);
        double dm = std::atof(argv[8]);
        int nsteps = std::atoi(argv[9]);
        bool emit_j = false;
        for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--json") emit_j = true;
        run_continue(m, x3, y3, Th, (mode == "--continue-m2") ? 1 : 2, dm, nsteps, emit_j);
    } else if (mode == "--hunt") {
        int n = (argc > 2) ? std::atoi(argv[2]) : 35000;
        double th = (argc > 3) ? std::atof(argv[3]) : 8.0;
        int sd = (argc > 4) ? std::atoi(argv[4]) : 0;
        run_hunt(n, th, sd);
    }
    return 0;
}
