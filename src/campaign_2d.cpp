#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
            h = std::min(h * fac, 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9 * std::pow(err, -0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate_to(State& y, double T, const Masses& m, double tol = 1e-12,
                         int max_steps = 120000, double* rmin_out = nullptr, double r_abort = 1e-6) {
    double t = 0.0, h = std::min(1e-3, T * 0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i = 0; i < 6; ++i) if (std::abs(y[i]) > 80.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

// Analyze topology & democracy over 1 period:
// - partner_swaps: number of times the identity of the closest pair changes
// - syzygies: number of collinear alignments (area of triangle changes sign)
// - rmin, rmax: min and max pairwise separation
// - is_choreography: whether bodies trace the same curve (approx check)
struct TopoInfo {
    int partner_swaps;
    int syzygies;
    double rmin;
    double rmax;
    double max_rad;
    bool democratic;
};

static TopoInfo analyze_topology(const State& y0, double T, const Masses& m) {
    const int N = 400;
    double dt = T / N;
    State y = y0;
    int prev_closest = -1;
    int swaps = 0;
    double prev_area = 0.0;
    int syz = 0;
    double rmin = 1e300, rmax = 0.0, max_rad = 0.0;

    for (int k = 0; k <= N; ++k) {
        double r12 = std::hypot(y[2] - y[0], y[3] - y[1]);
        double r13 = std::hypot(y[4] - y[0], y[5] - y[1]);
        double r23 = std::hypot(y[4] - y[2], y[5] - y[3]);
        rmin = std::min({rmin, r12, r13, r23});
        rmax = std::max({rmax, r12, r13, r23});
        for (int b = 0; b < 3; ++b) max_rad = std::max(max_rad, std::hypot(y[2*b], y[2*b+1]));

        int closest = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
        if (prev_closest != -1 && closest != prev_closest) swaps++;
        prev_closest = closest;

        double area = (y[2] - y[0]) * (y[5] - y[1]) - (y[3] - y[1]) * (y[4] - y[0]);
        if (k > 0 && ((prev_area > 0 && area < 0) || (prev_area < 0 && area > 0))) syz++;
        prev_area = area;

        if (k < N) integrate_to(y, dt, m, 1e-11, 50000);
    }
    return {swaps, syz, rmin, rmax, max_rad, swaps >= 2};
}

struct StabilityInfo {
    bool real_k;
    double k1, k2;
    bool linearly_stable;
    double stab_margin; // 2 - max(|k1|, |k2|)
    double max_lyap;
    int survived_100T;
    double drift_100T;
};

static StabilityInfo analyze_stability(const State& y0, double T, const Masses& m) {
    double M[12][12];
    const double eps = 1e-7;
    for (int col = 0; col < 12; ++col) {
        State yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        integrate_to(yp, T, m, 1e-13, 300000);
        integrate_to(ym, T, m, 1e-13, 300000);
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

    StabilityInfo info{};
    info.stab_margin = -1.0;
    if (disc >= -1e-5) {
        double sd = std::sqrt(std::max(0.0, disc));
        info.real_k = true;
        info.k1 = 0.5 * (s + sd);
        info.k2 = 0.5 * (s - sd);
        double mk = std::max(std::abs(info.k1), std::abs(info.k2));
        info.stab_margin = 2.0 - mk;
        info.linearly_stable = (mk < 2.0 - 1e-4);
        auto max_mag = [](double k) {
            if (std::abs(k) <= 2.0) return 1.0;
            double ak = std::abs(k);
            return 0.5 * (ak + std::sqrt(ak * ak - 4.0));
        };
        info.max_lyap = std::log(std::max(max_mag(info.k1), max_mag(info.k2))) / T;
    } else {
        info.real_k = false;
        info.k1 = 0.5 * s;
        info.k2 = 0.5 * std::sqrt(-disc);
        info.linearly_stable = false;
        info.max_lyap = 1.0;
    }

    if (info.linearly_stable) {
        State yp = y0;
        yp[6] += 1e-5;
        yp[8] -= 1e-5 * m[0] / m[1];
        double max_drift = 0.0;
        int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate_to(yp, T, m, 1e-11, 200000)) break;
            double d = state_dist(yp, y0);
            if (d > max_drift) max_drift = d;
            if (d > 5.0) break;
            surv = k + 1;
        }
        info.survived_100T = surv;
        info.drift_100T = max_drift;
    } else {
        info.survived_100T = 0;
        info.drift_100T = 999.0;
    }
    return info;
}

static void emit_json(const std::string& subspace, const Masses& m, const State& y0, double T,
                      double ret_err, const StabilityInfo& st, const TopoInfo& tp) {
    double E = compute_energy(y0, m);
    double L = compute_L(y0, m);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"dim\":2,\"subspace\":\"" << subspace << "\""
       << ",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "]"
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

// ============================================================
// MODE 1: Collinear L=0 (m1=m2=1, m3=m3)
// ============================================================
static State make_collinear(double vx, double vy, double m3) {
    return {-1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            vx, vy, vx, vy, -2.0 * vx / m3, -2.0 * vy / m3};
}

static bool refine_collinear(double& vx, double& vy, double& T, double m3, double& final_dist) {
    Masses m = {1.0, 1.0, m3};
    double lambda = 1e-4;
    auto eval_F = [&](double cvx, double cvy, double cT, State& diff) -> bool {
        if (cT < 1.0 || cT > 60.0) return false;
        State y0 = make_collinear(cvx, cvy, m3);
        State yT = y0;
        double rm = 1e300;
        if (!integrate_to(yT, cT, m, 1e-12, 120000, &rm, 5e-4)) return false;
        for (int i = 0; i < 12; ++i) diff[i] = yT[i] - y0[i];
        return true;
    };

    State F;
    if (!eval_F(vx, vy, T, F)) return false;
    double cost = 0.0;
    for (int i = 0; i < 12; ++i) cost += F[i] * F[i];

    for (int iter = 0; iter < 30; ++iter) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[12][3];
        State Fp, Fm;
        if (!eval_F(vx + eps, vy, T, Fp) || !eval_F(vx - eps, vy, T, Fm)) return false;
        for (int i = 0; i < 12; ++i) J[i][0] = (Fp[i] - Fm[i]) / (2.0 * eps);
        if (!eval_F(vx, vy + eps, T, Fp) || !eval_F(vx, vy - eps, T, Fm)) return false;
        for (int i = 0; i < 12; ++i) J[i][1] = (Fp[i] - Fm[i]) / (2.0 * eps);
        State y0 = make_collinear(vx, vy, m3), yT, fT;
        for (int i = 0; i < 12; ++i) yT[i] = y0[i] + F[i];
        deriv(yT, m, fT);
        for (int i = 0; i < 12; ++i) J[i][2] = fT[i];

        double M[3][4] = {{0}};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                double s = 0; for (int i = 0; i < 12; ++i) s += J[i][r] * J[i][c];
                M[r][c] = s;
            }
            M[r][r] += lambda * (1.0 + M[r][r]);
            double s = 0; for (int i = 0; i < 12; ++i) s -= J[i][r] * F[i];
            M[r][3] = s;
        }
        bool sing = false;
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col + 1; r < 3; ++r)
                if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
            if (std::abs(M[piv][col]) < 1e-18) { sing = true; break; }
            for (int c = col; c < 4; ++c) std::swap(M[col][c], M[piv][c]);
            for (int r = col + 1; r < 3; ++r) {
                double fac = M[r][col] / M[col][col];
                for (int c = col; c < 4; ++c) M[r][c] -= fac * M[col][c];
            }
        }
        if (sing) { lambda *= 10.0; continue; }
        double d[3] = {0};
        for (int r = 2; r >= 0; --r) {
            double s = M[r][3];
            for (int c = r + 1; c < 3; ++c) s -= M[r][c] * d[c];
            d[r] = s / M[r][r];
        }
        double nvx = vx + d[0], nvy = vy + d[1], nT = T + d[2];
        State nF;
        if (!eval_F(nvx, nvy, nT, nF)) { lambda *= 10.0; continue; }
        double ncost = 0; for (int i = 0; i < 12; ++i) ncost += nF[i] * nF[i];
        if (ncost < cost) {
            vx = nvx; vy = nvy; T = nT; F = nF; cost = ncost;
            lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0;
            if (lambda > 1e7) break;
        }
    }
    final_dist = std::sqrt(cost);
    return (final_dist < 1e-8 && T > 1.5 && vx > 1e-4 && vy > 1e-4);
}

static void run_collinear(double m3, int gridN, double Tmax) {
    Masses m = {1.0, 1.0, m3};
    double vmax_x = std::min(1.05, 0.85 * std::sqrt(std::min(m3, 1.5)));
    double vmax_y = std::min(1.15, 0.98 * std::sqrt(std::min(m3, 1.5)));
    struct Cell { double vx, vy, T, dist; };
    std::vector<std::vector<Cell>> grid(gridN, std::vector<Cell>(gridN));
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            for (int ix = tid; ix < gridN; ix += nthreads) {
                double vx = 0.03 + (vmax_x - 0.03) * ix / (gridN - 1);
                for (int iy = 0; iy < gridN; ++iy) {
                    double vy = 0.03 + (vmax_y - 0.03) * iy / (gridN - 1);
                    State y0 = make_collinear(vx, vy, m3);
                    if (compute_energy(y0, m) >= -0.05) { grid[ix][iy] = {vx, vy, 0, 1e9}; continue; }
                    State y = y0;
                    double t = 0, h = 1e-3, rmin = 1e300;
                    double best_d = 1e9, best_T = 0, pd = 0, p2d = 0, pt = 0;
                    bool left = false;
                    int steps = 0;
                    while (t < Tmax && steps < 25000) {
                        double hc = std::min(h, 0.025);
                        if (!dopri5_step(y, t, hc, Tmax, m, 1e-8, rmin, 1e-3)) break;
                        h = hc; ++steps;
                        double maxr = 0; for (int i = 0; i < 6; ++i) maxr = std::max(maxr, std::abs(y[i]));
                        if (maxr > 22.0) break;
                        double d = state_dist(y, y0);
                        if (t > 1.0 && d > 0.5) left = true;
                        if (left && steps >= 3 && pd < p2d && pd < d && pd < best_d) {
                            best_d = pd; best_T = pt;
                        }
                        p2d = pd; pd = d; pt = t;
                    }
                    grid[ix][iy] = {vx, vy, best_T, best_d};
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Cell> seeds;
    for (int ix = 1; ix < gridN - 1; ++ix) {
        for (int iy = 1; iy < gridN - 1; ++iy) {
            double d = grid[ix][iy].dist;
            if (d > 0.35 || grid[ix][iy].T < 1.5) continue;
            bool is_min = true;
            for (int dx = -1; dx <= 1 && is_min; ++dx)
                for (int dy = -1; dy <= 1 && is_min; ++dy)
                    if ((dx || dy) && grid[ix + dx][iy + dy].dist < d) is_min = false;
            if (is_min) seeds.push_back(grid[ix][iy]);
        }
    }
    std::sort(seeds.begin(), seeds.end(), [](const Cell& a, const Cell& b) { return a.dist < b.dist; });
    std::cerr << "[collinear m3=" << m3 << "] seeds=" << seeds.size() << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 3>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= seeds.size()) break;
                double vx = seeds[idx].vx, vy = seeds[idx].vy, T = seeds[idx].T, dist = 0;
                if (!refine_collinear(vx, vy, T, m3, dist)) continue;

                // Check if half-period or third-period is already periodic (primitive period!)
                State y0 = make_collinear(vx, vy, m3);
                for (int div : {2, 3, 4}) {
                    State yt = y0;
                    if (integrate_to(yt, T / div, m, 1e-12, 120000) && state_dist(yt, y0) < 1e-5) {
                        T /= div;
                        refine_collinear(vx, vy, T, m3, dist);
                        break;
                    }
                }
                y0 = make_collinear(vx, vy, m3);

                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup) {
                        if (std::hypot(f[0] - vx, f[1] - vy) < 1e-3 && std::abs(f[2] - T) < 0.05) {
                            dup = true; break;
                        }
                    }
                    if (dup) continue;
                    dedup.push_back({vx, vy, T});
                }

                auto st = analyze_stability(y0, T, m);
                if (!st.linearly_stable || st.survived_100T < 100) continue;
                auto tp = analyze_topology(y0, T, m);
                std::lock_guard<std::mutex> lk(mtx);
                emit_json("collinear_L0", m, y0, T, dist, st, tp);
            }
        });
    }
    for (auto& th : ths) th.join();
}

// ============================================================
// MODE 2: Perpendicular-crossing L != 0 (arbitrary masses)
// ============================================================
static State make_perp(double x2, double v1, double v2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0] * x1 + m[1] * x2) / m[2];
    double v3 = -(m[0] * v1 + m[1] * v2) / m[2];
    return {x1, 0.0, x2, 0.0, x3, 0.0, 0.0, v1, 0.0, v2, 0.0, v3};
}

static bool refine_perp(double& x2, double& v1, double& v2, double& Th, const Masses& m, double& ferr) {
    double lambda = 1e-3;
    auto eval_F = [&](double cx2, double cv1, double cv2, double cTh, std::array<double, 4>& F) -> bool {
        if (cTh < 0.4 || cTh > 25.0) return false;
        State y = make_perp(cx2, cv1, cv2, m);
        double d12 = std::abs(y[2] - y[0]), d13 = std::abs(y[4] - y[0]), d23 = std::abs(y[4] - y[2]);
        if (std::min({d12, d13, d23}) < 0.1) return false;
        double rm = 1e300;
        if (!integrate_to(y, cTh, m, 1e-11, 35000, &rm, 1e-2)) return false;
        F = {y[1], y[3], y[6], y[8]};
        return true;
    };
    std::array<double, 4> F;
    if (!eval_F(x2, v1, v2, Th, F)) return false;
    double cost = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];

    for (int iter = 0; iter < 28; ++iter) {
        if (std::sqrt(cost) < 1e-10) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[4][4];
        std::array<double, 4> Fp, Fm;
        double p[4] = {x2, v1, v2, Th};
        for (int col = 0; col < 4; ++col) {
            double pp[4] = {p[0], p[1], p[2], p[3]};
            double pm[4] = {p[0], p[1], p[2], p[3]};
            pp[col] += eps; pm[col] -= eps;
            if (!eval_F(pp[0], pp[1], pp[2], pp[3], Fp) || !eval_F(pm[0], pm[1], pm[2], pm[3], Fm)) return false;
            for (int r = 0; r < 4; ++r) J[r][col] = (Fp[r] - Fm[r]) / (2.0 * eps);
        }
        double A[4][5] = {{0}};
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                double s = 0; for (int i = 0; i < 4; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0; for (int i = 0; i < 4; ++i) s -= J[i][r] * F[i];
            A[r][4] = s;
        }
        bool sing = false;
        for (int col = 0; col < 4; ++col) {
            int piv = col;
            for (int r = col + 1; r < 4; ++r)
                if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (std::abs(A[piv][col]) < 1e-18) { sing = true; break; }
            for (int c = col; c < 5; ++c) std::swap(A[col][c], A[piv][c]);
            for (int r = col + 1; r < 4; ++r) {
                double fac = A[r][col] / A[col][col];
                for (int c = col; c < 5; ++c) A[r][c] -= fac * A[col][c];
            }
        }
        if (sing) { lambda *= 10.0; continue; }
        double d[4] = {0};
        for (int r = 3; r >= 0; --r) {
            double s = A[r][4];
            for (int c = r + 1; c < 4; ++c) s -= A[r][c] * d[c];
            d[r] = s / A[r][r];
        }
        double sn = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + d[3]*d[3]);
        if (sn > 0.4) for (int i = 0; i < 4; ++i) d[i] *= (0.4 / sn);
        std::array<double, 4> nF;
        if (!eval_F(x2 + d[0], v1 + d[1], v2 + d[2], Th + d[3], nF)) { lambda *= 10.0; continue; }
        double ncost = nF[0]*nF[0] + nF[1]*nF[1] + nF[2]*nF[2] + nF[3]*nF[3];
        if (ncost < cost) {
            x2 += d[0]; v1 += d[1]; v2 += d[2]; Th += d[3];
            F = nF; cost = ncost; lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0; if (lambda > 1e7) break;
        }
    }
    ferr = std::sqrt(cost);
    return (ferr < 1e-9 && Th > 0.5);
}

static void run_perp(const Masses& m, int nsamples, double Thmax, int seed_offset) {
    struct Seed { double x2, v1, v2, Th, err; int swaps; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(777001 + seed_offset * 10007 + tid * 9973);
            std::uniform_real_distribution<double> ux2(-0.7, 2.2);
            std::uniform_real_distribution<double> uv(-1.4, 1.4);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double x2 = ux2(rng), v1 = uv(rng), v2 = uv(rng);
                State y0 = make_perp(x2, v1, v2, m);
                double d12 = std::abs(y0[2] - y0[0]), d13 = std::abs(y0[4] - y0[0]), d23 = std::abs(y0[4] - y0[2]);
                if (std::min({d12, d13, d23}) < 0.2) continue;
                if (compute_energy(y0, m) >= -0.05) continue;

                State y = y0;
                double t = 0, h = 1e-3, rmin = 1e300;
                double best_e = 1e9, best_t = 0, pe = 1e9, p2e = 1e9, pt = 0;
                int prev_cl = -1, swaps = 0, best_swaps = 0;
                int steps = 0;
                while (t < Thmax && steps < 12000) {
                    double hc = std::min(h, 0.025);
                    if (!dopri5_step(y, t, hc, Thmax, m, 1e-8, rmin, 1e-2)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i = 0; i < 6; ++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 18.0) break;

                    double r12 = std::hypot(y[2]-y[0], y[3]-y[1]);
                    double r13 = std::hypot(y[4]-y[0], y[5]-y[1]);
                    double r23 = std::hypot(y[4]-y[2], y[5]-y[3]);
                    int cl = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
                    if (prev_cl != -1 && cl != prev_cl) swaps++;
                    prev_cl = cl;

                    double e = std::sqrt(y[1]*y[1] + y[3]*y[3] + y[6]*y[6] + y[8]*y[8]);
                    if (t > 0.5 && steps >= 3 && pe < p2e && pe < e && pe < best_e) {
                        best_e = pe; best_t = pt; best_swaps = swaps;
                    }
                    p2e = pe; pe = e; pt = t;
                }
                if (best_e < 0.25 && best_t > 0.5) {
                    tseeds[tid].push_back({x2, v1, v2, best_t, best_e, best_swaps});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all_seeds;
    for (auto& v : tseeds) all_seeds.insert(all_seeds.end(), v.begin(), v.end());
    // Prioritize democratic seeds (swaps > 0) first, then lowest error!
    std::sort(all_seeds.begin(), all_seeds.end(), [](const Seed& a, const Seed& b) {
        if ((a.swaps > 0) != (b.swaps > 0)) return (a.swaps > 0) > (b.swaps > 0);
        return a.err < b.err;
    });
    size_t limit = std::min(all_seeds.size(), (size_t)400);
    std::cerr << "[perp m=(" << m[0] << "," << m[1] << "," << m[2] << ")] seeds=" << all_seeds.size()
              << " refining=" << limit << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 2>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                double x2 = all_seeds[idx].x2, v1 = all_seeds[idx].v1, v2 = all_seeds[idx].v2, Th = all_seeds[idx].Th, ferr = 0;
                if (!refine_perp(x2, v1, v2, Th, m, ferr)) continue;
                double T = 2.0 * Th;
                State y0 = make_perp(x2, v1, v2, m);
                State yT = y0;
                double rmin = 1e300;
                if (!integrate_to(yT, T, m, 1e-13, 250000, &rmin, 5e-3)) continue;
                double ret = state_dist(yT, y0);
                if (ret > 1e-7) continue;

                double E = compute_energy(y0, m);
                double L = compute_L(y0, m);
                double inv1 = E * std::pow(T, 2.0 / 3.0);
                double inv2 = std::abs(L) * std::sqrt(std::max(0.0, -E));
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup) {
                        if (std::abs(f[0] - inv1) < 2e-2 && std::abs(f[1] - inv2) < 2e-2) { dup = true; break; }
                    }
                    if (dup) continue;
                    dedup.push_back({inv1, inv2});
                }

                auto st = analyze_stability(y0, T, m);
                if (!st.linearly_stable || st.survived_100T < 100) continue;
                auto tp = analyze_topology(y0, T, m);
                std::lock_guard<std::mutex> lk(mtx);
                emit_json("perp_L_nonzero", m, y0, T, ret, st, tp);
            }
        });
    }
    for (auto& th : ths) th.join();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n  campaign_2d --collinear <m3> <gridN> <Tmax>\n"
                  << "  campaign_2d --perp <m1> <m2> <m3> <nsamples> <Thmax> [seed]\n";
        return 1;
    }
    std::string mode = argv[1];
    if (mode == "--collinear") {
        double m3 = std::atof(argv[2]);
        int gridN = (argc > 3) ? std::atoi(argv[3]) : 180;
        double Tmax = (argc > 4) ? std::atof(argv[4]) : 32.0;
        run_collinear(m3, gridN, Tmax);
    } else if (mode == "--perp") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        int nsamples = (argc > 5) ? std::atoi(argv[5]) : 16000;
        double Thmax = (argc > 6) ? std::atof(argv[6]) : 12.0;
        int seed = (argc > 7) ? std::atoi(argv[7]) : 0;
        run_perp(m, nsamples, Thmax, seed);
    }
    return 0;
}
