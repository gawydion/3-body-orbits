#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
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
            double fx = dx * inv_r3;
            double fy = dy * inv_r3;
            ax[i] += m[j] * fx;
            ay[i] += m[j] * fy;
            ax[j] -= m[i] * fx;
            ay[j] -= m[i] * fy;
        }
    }
    ds[6] = ax[0]; ds[7] = ay[0];
    ds[8] = ax[1]; ds[9] = ay[1];
    ds[10] = ax[2]; ds[11] = ay[2];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double compute_energy(const State& s, const Masses& m) {
    double K = 0.0;
    for (int i = 0; i < 3; ++i)
        K += 0.5 * m[i] * (s[6 + 2 * i] * s[6 + 2 * i] + s[7 + 2 * i] * s[7 + 2 * i]);
    double U = 0.0;
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
                        double tol, double& rmin_track) {
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
        if (rmin_local < 1e-5) return false;

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
        if (rmin_local < 1e-5) return false;

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
            h = std::min(h * fac, 0.1);
            return true;
        }
        h *= std::max(0.1, 0.9 * std::pow(err, -0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate_to(State& y, double T, const Masses& m, double tol = 1e-12,
                         int max_steps = 200000, double* rmin_out = nullptr) {
    double t = 0.0;
    double h = std::min(1e-3, T * 0.01);
    double rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
        for (int i = 0; i < 6; ++i) if (std::abs(y[i]) > 100.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

// Perpendicular crossing state:
// x1 = -1, x2 = x2, x3 = -(m1*(-1) + m2*x2)/m3
// vy1 = v1, vy2 = v2, vy3 = -(m1*v1 + m2*v2)/m3
// All y_i = 0, vx_i = 0
static State make_perp_state(double x2, double v1, double v2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0] * x1 + m[1] * x2) / m[2];
    double v3 = -(m[0] * v1 + m[1] * v2) / m[2];
    return {x1, 0.0, x2, 0.0, x3, 0.0,
            0.0, v1, 0.0, v2, 0.0, v3};
}

// Perpendicular crossing residual at t = T/2:
// F = (y1, y2, vx1, vx2)
static inline double perp_err(const State& s) {
    return std::sqrt(s[1]*s[1] + s[3]*s[3] + s[6]*s[6] + s[8]*s[8]);
}

// Refine (x2, v1, v2, Thalf) so that F(Thalf) = (y1, y2, vx1, vx2) = 0
static bool refine_perp(double& x2, double& v1, double& v2, double& Thalf, const Masses& m,
                        double& final_err, double& rmin_out) {
    double lambda = 1e-3;
    auto eval_F = [&](double cx2, double cv1, double cv2, double cTh, std::array<double, 4>& F, double& rmin) -> bool {
        if (cTh < 0.3 || cTh > 30.0) return false;
        State y = make_perp_state(cx2, cv1, cv2, m);
        // Check initial separation
        double d12 = std::abs(y[2] - y[0]);
        double d13 = std::abs(y[4] - y[0]);
        double d23 = std::abs(y[4] - y[2]);
        if (std::min({d12, d13, d23}) < 0.12) return false;
        if (!integrate_to(y, cTh, m, 1e-11, 35000, &rmin)) return false;
        if (rmin < 0.03) return false;
        F = {y[1], y[3], y[6], y[8]};
        return true;
    };

    std::array<double, 4> F;
    double rmin = 1e300;
    if (!eval_F(x2, v1, v2, Thalf, F, rmin)) return false;
    double cost = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];

    for (int iter = 0; iter < 28; ++iter) {
        if (std::sqrt(cost) < 1e-10) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[4][4];
        std::array<double, 4> Fp, Fm;
        double rm_tmp;

        double p[4] = {x2, v1, v2, Thalf};
        for (int col = 0; col < 4; ++col) {
            double pp[4] = {p[0], p[1], p[2], p[3]};
            double pm[4] = {p[0], p[1], p[2], p[3]};
            pp[col] += eps;
            pm[col] -= eps;
            if (!eval_F(pp[0], pp[1], pp[2], pp[3], Fp, rm_tmp) ||
                !eval_F(pm[0], pm[1], pm[2], pm[3], Fm, rm_tmp)) return false;
            for (int row = 0; row < 4; ++row) J[row][col] = (Fp[row] - Fm[row]) / (2.0 * eps);
        }

        // (J^T J + lambda I) delta = -J^T F
        double A[4][5] = {{0}};
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                double s = 0.0;
                for (int i = 0; i < 4; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0.0;
            for (int i = 0; i < 4; ++i) s -= J[i][r] * F[i];
            A[r][4] = s;
        }

        bool singular = false;
        for (int col = 0; col < 4; ++col) {
            int piv = col;
            for (int r = col + 1; r < 4; ++r)
                if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (std::abs(A[piv][col]) < 1e-18) { singular = true; break; }
            for (int c = col; c < 5; ++c) std::swap(A[col][c], A[piv][c]);
            for (int r = col + 1; r < 4; ++r) {
                double fac = A[r][col] / A[col][col];
                for (int c = col; c < 5; ++c) A[r][c] -= fac * A[col][c];
            }
        }
        if (singular) { lambda *= 10.0; continue; }

        double delta[4] = {0};
        for (int r = 3; r >= 0; --r) {
            double s = A[r][4];
            for (int c = r + 1; c < 4; ++c) s -= A[r][c] * delta[c];
            delta[r] = s / A[r][r];
        }

        // Clamp step size
        double step_norm = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] + delta[3]*delta[3]);
        if (step_norm > 0.5) {
            for (int i = 0; i < 4; ++i) delta[i] *= (0.5 / step_norm);
        }

        double nx2 = x2 + delta[0], nv1 = v1 + delta[1], nv2 = v2 + delta[2], nTh = Thalf + delta[3];
        std::array<double, 4> nF;
        double nrmin = 1e300;
        if (!eval_F(nx2, nv1, nv2, nTh, nF, nrmin)) { lambda *= 10.0; continue; }
        double ncost = nF[0]*nF[0] + nF[1]*nF[1] + nF[2]*nF[2] + nF[3]*nF[3];

        if (ncost < cost) {
            x2 = nx2; v1 = nv1; v2 = nv2; Thalf = nTh;
            F = nF; cost = ncost; rmin = nrmin;
            lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0;
            if (lambda > 1e8) break;
        }
    }

    final_err = std::sqrt(cost);
    rmin_out = rmin;
    return (final_err < 1e-9 && Thalf > 0.5);
}

struct StabilityInfo {
    bool real_k;
    double k1, k2;
    bool linearly_stable;
    double max_lyap;
    double drift50T;
};

static StabilityInfo analyze_stability(const State& y0, double T, const Masses& m) {
    double M[12][12];
    const double eps = 1e-7;
    for (int col = 0; col < 12; ++col) {
        State yp = y0, ym = y0;
        yp[col] += eps;
        ym[col] -= eps;
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
    if (disc >= -1e-6) {
        double sd = std::sqrt(std::max(0.0, disc));
        info.real_k = true;
        info.k1 = 0.5 * (s + sd);
        info.k2 = 0.5 * (s - sd);
        info.linearly_stable = (std::abs(info.k1) < 2.0 - 1e-4) && (std::abs(info.k2) < 2.0 - 1e-4);
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
        double kr = info.k1, ki = info.k2;
        double dr = kr*kr - ki*ki - 4.0, di = 2.0*kr*ki;
        double mag = std::hypot(dr, di);
        double sr = std::sqrt(0.5*(mag + dr));
        double si = (di >= 0 ? 1.0 : -1.0) * std::sqrt(0.5*(mag - dr));
        double l1 = std::hypot(0.5*(kr + sr), 0.5*(ki + si));
        double l2 = std::hypot(0.5*(kr - sr), 0.5*(ki - si));
        info.max_lyap = std::log(std::max(l1, l2)) / T;
    }

    // Empirical 50-period stability test with 1e-5 perturbation preserving momentum
    State yp = y0;
    yp[6] += 1e-5;
    yp[8] -= 1e-5 * m[0] / m[1];
    double max_drift = 0.0;
    bool survived = true;
    for (int k = 0; k < 50; ++k) {
        if (!integrate_to(yp, T, m, 1e-11, 200000)) { survived = false; max_drift = 999.0; break; }
        double d = state_dist(yp, y0);
        if (d > max_drift) max_drift = d;
        if (d > 5.0) { survived = false; break; }
    }
    info.drift50T = survived ? max_drift : 999.0;
    return info;
}

int main(int argc, char** argv) {
    double m1 = (argc > 1) ? std::atof(argv[1]) : 1.0;
    double m2 = (argc > 2) ? std::atof(argv[2]) : 1.0;
    double m3 = (argc > 3) ? std::atof(argv[3]) : 1.0;
    int nsamples = (argc > 4) ? std::atoi(argv[4]) : 8000;
    double Thalf_max = (argc > 5) ? std::atof(argv[5]) : 12.0;
    Masses m = {m1, m2, m3};

    std::cerr << "Searching perpendicular-crossing (L != 0) orbits for masses=("
              << m1 << "," << m2 << "," << m3 << "), samples=" << nsamples << "...\n";

    struct Seed { double x2, v1, v2, Thalf, err; };
    std::vector<Seed> all_seeds;
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> threads;

    for (int tid = 0; tid < nthreads; ++tid) {
        threads.emplace_back([&, tid]() {
            std::mt19937_64 rng(421337 + tid * 9973 + (int)(m3 * 1000) + (int)(m2 * 100));
            std::uniform_real_distribution<double> ux2(-0.6, 1.8);
            std::uniform_real_distribution<double> uv(-1.3, 1.3);

            int count = nsamples / nthreads;
            for (int k = 0; k < count; ++k) {
                double x2 = ux2(rng);
                double v1 = uv(rng);
                double v2 = uv(rng);
                State y0 = make_perp_state(x2, v1, v2, m);
                double d12 = std::abs(y0[2] - y0[0]);
                double d13 = std::abs(y0[4] - y0[0]);
                double d23 = std::abs(y0[4] - y0[2]);
                if (std::min({d12, d13, d23}) < 0.25) continue;
                if (compute_energy(y0, m) >= -0.05) continue; // must be bound

                State y = y0;
                double t = 0.0, h = 1e-3, rmin = 1e300;
                double best_e = 1e300, best_t = 0.0;
                double prev_e = 1e300, prev2_e = 1e300, prev_t = 0.0;
                int steps = 0;
                while (t < Thalf_max && steps < 15000) {
                    double hc = std::min(h, 0.025);
                    if (!dopri5_step(y, t, hc, Thalf_max, m, 1e-8, rmin)) break;
                    h = hc;
                    ++steps;
                    double maxr = 0.0;
                    for (int i = 0; i < 6; ++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 20.0) break;
                    double e = perp_err(y);
                    if (t > 0.6 && steps >= 3) {
                        if (prev_e < prev2_e && prev_e < e && prev_e < best_e) {
                            best_e = prev_e;
                            best_t = prev_t;
                        }
                    }
                    prev2_e = prev_e; prev_e = e; prev_t = t;
                }
                if (best_e < 0.28 && best_t > 0.6) {
                    tseeds[tid].push_back({x2, v1, v2, best_t, best_e});
                }
            }
        });
    }
    for (auto& th : threads) th.join();
    for (auto& v : tseeds) all_seeds.insert(all_seeds.end(), v.begin(), v.end());
    std::sort(all_seeds.begin(), all_seeds.end(), [](const Seed& a, const Seed& b){ return a.err < b.err; });
    std::cerr << "Found " << all_seeds.size() << " perp seeds. Refining top 250...\n";

    struct Orbit {
        double x2, v1, v2, T, E, L, rmin, ret_dist;
        StabilityInfo st;
    };
    std::vector<Orbit> found;
    int limit = std::min((int)all_seeds.size(), 100);

    for (int idx = 0; idx < limit; ++idx) {
        double x2 = all_seeds[idx].x2, v1 = all_seeds[idx].v1, v2 = all_seeds[idx].v2, Th = all_seeds[idx].Thalf;
        double err = 0.0, rmin = 0.0;
        if (!refine_perp(x2, v1, v2, Th, m, err, rmin)) continue;
        if (rmin < 0.05) continue;

        double T = 2.0 * Th;
        State y0 = make_perp_state(x2, v1, v2, m);
        State yT = y0;
        if (!integrate_to(yT, T, m, 1e-13, 300000, &rmin)) continue;
        double ret_dist = state_dist(yT, y0);
        if (ret_dist > 1e-8) continue;

        double E = compute_energy(y0, m);
        double L = compute_L(y0, m);

        // Deduplicate by (E, |L|, T) scale-invariant combo: E*T^(2/3) and |L|*sqrt(-E)
        double inv1 = E * std::pow(T, 2.0 / 3.0);
        double inv2 = std::abs(L) * std::sqrt(std::max(0.0, -E));
        bool dup = false;
        for (const auto& f : found) {
            double finv1 = f.E * std::pow(f.T, 2.0 / 3.0);
            double finv2 = std::abs(f.L) * std::sqrt(std::max(0.0, -f.E));
            if (std::abs(inv1 - finv1) < 1e-2 && std::abs(inv2 - finv2) < 1e-2) {
                dup = true;
                break;
            }
        }
        if (dup) continue;

        auto st = analyze_stability(y0, T, m);
        found.push_back({x2, v1, v2, T, E, L, rmin, ret_dist, st});
    }

    std::sort(found.begin(), found.end(), [](const Orbit& a, const Orbit& b) {
        if (a.st.linearly_stable != b.st.linearly_stable) return a.st.linearly_stable > b.st.linearly_stable;
        if (a.st.drift50T != b.st.drift50T) return a.st.drift50T < b.st.drift50T;
        return a.T < b.T;
    });

    std::cout << std::setprecision(12);
    std::cout << "FOUND " << found.size() << " PERP ORBITS m=(" << m1 << "," << m2 << "," << m3 << "):\n";
    for (size_t i = 0; i < found.size(); ++i) {
        const auto& o = found[i];
        std::cout << "#" << (i + 1)
                  << " x2=" << o.x2 << " v1=" << o.v1 << " v2=" << o.v2
                  << " T=" << o.T << " E=" << o.E << " L=" << o.L
                  << " rmin=" << o.rmin << " ret=" << o.ret_dist
                  << " stable=" << (o.st.linearly_stable ? "YES" : "NO")
                  << " k1=" << o.st.k1 << " k2=" << o.st.k2
                  << (o.st.real_k ? " (real)" : " (cplx)")
                  << " lyap=" << o.st.max_lyap
                  << " drift50T=" << o.st.drift50T << "\n";
    }
    return 0;
}
