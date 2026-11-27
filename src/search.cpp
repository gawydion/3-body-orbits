#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Planar 3-body problem (G = 1)
// State vector y[12]:
//   y[0..5]  = x1, y1, x2, y2, x3, y3
//   y[6..11] = vx1, vy1, vx2, vy2, vx3, vy3

using State = std::array<double, 12>;
using Masses = std::array<double, 3>;

struct Diagnostics {
    double E;
    double L;
    double r_min;
};

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
    for (int i = 0; i < 3; ++i) {
        K += 0.5 * m[i] * (s[6 + 2 * i] * s[6 + 2 * i] + s[7 + 2 * i] * s[7 + 2 * i]);
    }
    double U = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[2 * j] - s[2 * i];
            double dy = s[2 * j + 1] - s[2 * i + 1];
            U -= m[i] * m[j] / std::sqrt(dx * dx + dy * dy);
        }
    }
    return K + U;
}

static inline double compute_L(const State& s, const Masses& m) {
    double L = 0.0;
    for (int i = 0; i < 3; ++i) {
        L += m[i] * (s[2 * i] * s[7 + 2 * i] - s[2 * i + 1] * s[6 + 2 * i]);
    }
    return L;
}

static inline double state_dist(const State& a, const State& b) {
    double sum = 0.0;
    for (int i = 0; i < 12; ++i) {
        double d = a[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

// Dormand-Prince 5(4) adaptive step
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
        if (rmin_local < 1e-6) return false; // collision

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
        if (rmin_local < 1e-6) return false;

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
        double fac = std::max(0.1, 0.9 * std::pow(err, -0.25));
        h *= fac;
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
        // Escape check
        for (int i = 0; i < 6; ++i) {
            if (std::abs(y[i]) > 100.0) return false;
        }
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

// Build symmetric collinear initial state:
// r1 = (-1, 0), r2 = (1, 0), r3 = (0, 0)
// v1 = v2 = (vx, vy), v3 = (-2*vx/m3, -2*vy/m3)
// m1 = m2 = 1, m3 = m3
static State make_collinear_state(double vx, double vy, double m3) {
    return {-1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            vx, vy, vx, vy, -2.0 * vx / m3, -2.0 * vy / m3};
}

// Return distance at period T for collinear state
struct ReturnCandidate {
    double vx, vy, T, dist, rmin, m3;
};

// Scan a trajectory for minimum return distance to y0
static ReturnCandidate scan_return_collinear(double vx, double vy, double m3, double Tmax,
                                             double tol = 1e-9, int max_steps = 40000) {
    Masses m = {1.0, 1.0, m3};
    State y0 = make_collinear_state(vx, vy, m3);
    State y = y0;
    double t = 0.0;
    double h = 1e-3;
    double rmin_track = 1e300;
    double best_d = 1e300;
    double best_T = 0.0;
    double prev_d = 0.0;
    double prev2_d = 0.0;
    double prev_t = 0.0;
    bool left_origin = false;

    int steps = 0;
    while (t < Tmax && steps < max_steps) {
        double h_cap = std::min(h, 0.02);
        if (!dopri5_step(y, t, h_cap, Tmax, m, tol, rmin_track)) break;
        h = h_cap;
        ++steps;

        double max_r = 0.0;
        for (int i = 0; i < 6; ++i) max_r = std::max(max_r, std::abs(y[i]));
        if (max_r > 25.0) break;

        double d = state_dist(y, y0);
        if (t > 1.0 && d > 0.5) left_origin = true;
        if (left_origin && steps >= 3) {
            // Check local minimum at prev_t
            if (prev_d < prev2_d && prev_d < d && prev_d < best_d) {
                best_d = prev_d;
                best_T = prev_t;
            }
        }
        prev2_d = prev_d;
        prev_d = d;
        prev_t = t;
    }
    return {vx, vy, best_T, best_d, rmin_track, m3};
}

// Gauss-Newton / Levenberg-Marquardt refinement of (vx, vy, T) for collinear subspace
// Minimizes || Phi_T(y0(vx, vy)) - y0(vx, vy) ||^2
static bool refine_collinear(double& vx, double& vy, double& T, double m3,
                             double& final_dist, double& rmin_out) {
    Masses m = {1.0, 1.0, m3};
    double lambda = 1e-4;

    auto eval_F = [&](double cvx, double cvy, double cT, State& diff, double& rmin) -> bool {
        State y0 = make_collinear_state(cvx, cvy, m3);
        State yT = y0;
        if (!integrate_to(yT, cT, m, 1e-13, 300000, &rmin)) return false;
        for (int i = 0; i < 12; ++i) diff[i] = yT[i] - y0[i];
        return true;
    };

    State F;
    double rmin = 1e300;
    if (!eval_F(vx, vy, T, F, rmin)) return false;
    double cost = 0.0;
    for (int i = 0; i < 12; ++i) cost += F[i] * F[i];

    for (int iter = 0; iter < 40; ++iter) {
        if (std::sqrt(cost) < 1e-12) break;

        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[12][3];
        State Fp, Fm;
        double rm_tmp;

        // Column 0: dF/dvx
        if (!eval_F(vx + eps, vy, T, Fp, rm_tmp) || !eval_F(vx - eps, vy, T, Fm, rm_tmp)) return false;
        for (int i = 0; i < 12; ++i) J[i][0] = (Fp[i] - Fm[i]) / (2.0 * eps);

        // Column 1: dF/dvy
        if (!eval_F(vx, vy + eps, T, Fp, rm_tmp) || !eval_F(vx, vy - eps, T, Fm, rm_tmp)) return false;
        for (int i = 0; i < 12; ++i) J[i][1] = (Fp[i] - Fm[i]) / (2.0 * eps);

        // Column 2: dF/dT = f(y(T))
        State y0 = make_collinear_state(vx, vy, m3);
        State yT = y0;
        for (int i = 0; i < 12; ++i) yT[i] = y0[i] + F[i];
        State fT;
        deriv(yT, m, fT);
        for (int i = 0; i < 12; ++i) J[i][2] = fT[i];

        // Normal equations (J^T J + lambda I) delta = -J^T F
        double A[3][3] = {{0}}, b[3] = {0};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int i = 0; i < 12; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0.0;
            for (int i = 0; i < 12; ++i) s -= J[i][r] * F[i];
            b[r] = s;
        }

        // Solve 3x3 linear system with partial pivoting
        double M[3][4];
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) M[r][c] = A[r][c];
            M[r][3] = b[r];
        }
        bool singular = false;
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col + 1; r < 3; ++r) {
                if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
            }
            if (std::abs(M[piv][col]) < 1e-18) { singular = true; break; }
            for (int c = col; c < 4; ++c) std::swap(M[col][c], M[piv][c]);
            for (int r = col + 1; r < 3; ++r) {
                double fac = M[r][col] / M[col][col];
                for (int c = col; c < 4; ++c) M[r][c] -= fac * M[col][c];
            }
        }
        if (singular) { lambda *= 10.0; continue; }
        double delta[3] = {0};
        for (int r = 2; r >= 0; --r) {
            double s = M[r][3];
            for (int c = r + 1; c < 3; ++c) s -= M[r][c] * delta[c];
            delta[r] = s / M[r][r];
        }

        double nvx = vx + delta[0];
        double nvy = vy + delta[1];
        double nT  = T  + delta[2];
        if (nT <= 0.5 || nT > 100.0) { lambda *= 10.0; continue; }

        State nF;
        double nrmin = 1e300;
        if (!eval_F(nvx, nvy, nT, nF, nrmin)) { lambda *= 10.0; continue; }
        double ncost = 0.0;
        for (int i = 0; i < 12; ++i) ncost += nF[i] * nF[i];

        if (ncost < cost) {
            vx = nvx; vy = nvy; T = nT;
            F = nF; cost = ncost; rmin = nrmin;
            lambda = std::max(1e-12, lambda * 0.2);
        } else {
            lambda *= 5.0;
            if (lambda > 1e8) break;
        }
    }

    final_dist = std::sqrt(cost);
    rmin_out = rmin;
    return (final_dist < 1e-9 && T > 1.0);
}

// Compute 12x12 monodromy matrix M = dPhi_T / dy0 via central differences,
// and extract Broucke stability indices k1, k2.
struct StabilityInfo {
    double tr1, tr2;
    bool real_k;
    double k1, k2;
    double k_real, k_imag; // if complex
    bool linearly_stable;
    double max_lyap;       // max Re(ln(lambda))/T
    double perturb_drift_50T; // empirical drift over 50 periods with 1e-5 perturbation
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
        for (int row = 0; row < 12; ++row) {
            M[row][col] = (yp[row] - ym[row]) / (2.0 * eps);
        }
    }

    double tr1 = 0.0;
    for (int i = 0; i < 12; ++i) tr1 += M[i][i];

    double tr2 = 0.0;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            tr2 += M[i][j] * M[j][i];
        }
    }

    double s = tr1 - 8.0;
    double p = 0.5 * (s * s - (tr2 - 4.0));
    double disc = s * s - 4.0 * p;

    StabilityInfo info{};
    info.tr1 = tr1;
    info.tr2 = tr2;
    info.max_lyap = 0.0;

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
        double lmax = std::max(max_mag(info.k1), max_mag(info.k2));
        info.max_lyap = std::log(lmax) / T;
    } else {
        info.real_k = false;
        info.k_real = 0.5 * s;
        info.k_imag = 0.5 * std::sqrt(-disc);
        info.linearly_stable = false;
        // Complex k => lambda + 1/lambda = k_r + i k_i => |lambda| > 1
        // Solve lambda^2 - k*lambda + 1 = 0
        double kr = info.k_real, ki = info.k_imag;
        // sqrt(k^2 - 4)
        double dr = kr * kr - ki * ki - 4.0;
        double di = 2.0 * kr * ki;
        double mag = std::hypot(dr, di);
        double sr = std::sqrt(0.5 * (mag + dr));
        double si = (di >= 0 ? 1.0 : -1.0) * std::sqrt(0.5 * (mag - dr));
        double l1 = std::hypot(0.5 * (kr + sr), 0.5 * (ki + si));
        double l2 = std::hypot(0.5 * (kr - sr), 0.5 * (ki - si));
        info.max_lyap = std::log(std::max(l1, l2)) / T;
    }

    // Empirical perturbation test: perturb vx1, vx2 by +1e-5, vx3 by -2e-5/m3, integrate 50 periods
    State yp = y0;
    yp[6] += 1e-5;
    yp[8] += 1e-5;
    yp[10] -= 2e-5 / m[2];
    double max_drift = 0.0;
    bool survived = true;
    for (int p_idx = 0; p_idx < 50; ++p_idx) {
        if (!integrate_to(yp, T, m, 1e-11, 200000)) {
            survived = false;
            max_drift = 999.0;
            break;
        }
        double d = state_dist(yp, y0);
        if (d > max_drift) max_drift = d;
        if (d > 5.0) { survived = false; break; }
    }
    info.perturb_drift_50T = survived ? max_drift : 999.0;
    return info;
}

int main(int argc, char** argv) {
    // Quick sanity test: Chenciner-Montgomery Figure-8 orbit
    if (argc > 1 && std::string(argv[1]) == "--test") {
        double vx = 0.347111, vy = 0.532728, T = 6.32591;
        double dist = 0.0, rmin = 0.0;
        bool ok = refine_collinear(vx, vy, T, 1.0, dist, rmin);
        State y0 = make_collinear_state(vx, vy, 1.0);
        Masses m = {1.0, 1.0, 1.0};
        State yT = y0;
        integrate_to(yT, T, m, 1e-13);
        double dE = std::abs(compute_energy(yT, m) - compute_energy(y0, m));
        double dL = std::abs(compute_L(yT, m) - compute_L(y0, m));
        auto st = analyze_stability(y0, T, m);

        std::cout << std::setprecision(14);
        std::cout << "Figure-8 test: ok=" << ok << "\n";
        std::cout << "  vx=" << vx << " vy=" << vy << " T=" << T << "\n";
        std::cout << "  return_dist=" << dist << " rmin=" << rmin << " dE=" << dE << " dL=" << dL << "\n";
        std::cout << "  tr1=" << st.tr1 << " tr2=" << st.tr2 << " k1=" << st.k1 << " k2=" << st.k2
                  << " stable=" << st.linearly_stable << " drift50T=" << st.perturb_drift_50T << "\n";
        return 0;
    }

    // Parse optional mass m3, grid size, Tmax
    double m3 = (argc > 1) ? std::atof(argv[1]) : 1.0;
    int gridN = (argc > 2) ? std::atoi(argv[2]) : 160;
    double Tmax = (argc > 3) ? std::atof(argv[3]) : 25.0;

    std::cerr << "Searching collinear subspace m3=" << m3 << " grid=" << gridN << "x" << gridN
              << " Tmax=" << Tmax << "...\n";

    // Velocity bounds scale roughly with sqrt(m3) for m3 < 1
    double vmax_x = (m3 < 1.0) ? 0.85 * std::sqrt(m3) : 0.95;
    double vmax_y = (m3 < 1.0) ? 0.95 * std::sqrt(m3) : 1.05;

    std::vector<std::vector<ReturnCandidate>> grid(gridN, std::vector<ReturnCandidate>(gridN));
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;

    for (int tid = 0; tid < nthreads; ++tid) {
        threads.emplace_back([&, tid]() {
            for (int ix = tid; ix < gridN; ix += nthreads) {
                double vx = 0.02 + (vmax_x - 0.02) * ix / (gridN - 1);
                for (int iy = 0; iy < gridN; ++iy) {
                    double vy = 0.02 + (vmax_y - 0.02) * iy / (gridN - 1);
                    grid[ix][iy] = scan_return_collinear(vx, vy, m3, Tmax, 1e-8, 30000);
                }
            }
        });
    }
    for (auto& th : threads) th.join();

    // Collect local minima on the 2D return-distance landscape
    std::vector<ReturnCandidate> seeds;
    for (int ix = 1; ix < gridN - 1; ++ix) {
        for (int iy = 1; iy < gridN - 1; ++iy) {
            double d = grid[ix][iy].dist;
            if (d > 0.35 || grid[ix][iy].T < 1.5) continue;
            bool is_min = true;
            for (int dx = -1; dx <= 1 && is_min; ++dx) {
                for (int dy = -1; dy <= 1 && is_min; ++dy) {
                    if (dx == 0 && dy == 0) continue;
                    if (grid[ix + dx][iy + dy].dist < d) is_min = false;
                }
            }
            if (is_min) seeds.push_back(grid[ix][iy]);
        }
    }
    std::sort(seeds.begin(), seeds.end(), [](const ReturnCandidate& a, const ReturnCandidate& b) {
        return a.dist < b.dist;
    });
    std::cerr << "Found " << seeds.size() << " local minima seeds. Refining...\n";

    struct DiscoveredOrbit {
        double m3, vx, vy, T, dist, rmin, E;
        StabilityInfo st;
    };
    std::vector<DiscoveredOrbit> found;

    for (const auto& s : seeds) {
        double vx = s.vx, vy = s.vy, T = s.T;
        double dist = 0.0, rmin = 0.0;
        if (!refine_collinear(vx, vy, T, m3, dist, rmin)) continue;
        if (vx <= 1e-4 || vy <= 1e-4 || rmin < 1e-3) continue;

        // Deduplicate (also check if T is a multiple of an already-found orbit at same vx, vy)
        bool dup = false;
        for (const auto& f : found) {
            if (std::hypot(f.vx - vx, f.vy - vy) < 1e-3) {
                dup = true;
                break;
            }
        }
        if (dup) continue;

        State y0 = make_collinear_state(vx, vy, m3);
        Masses m = {1.0, 1.0, m3};
        double E = compute_energy(y0, m);
        auto st = analyze_stability(y0, T, m);
        found.push_back({m3, vx, vy, T, dist, rmin, E, st});
    }

    // Sort: linearly stable first, then by shortest period T
    std::sort(found.begin(), found.end(), [](const DiscoveredOrbit& a, const DiscoveredOrbit& b) {
        if (a.st.linearly_stable != b.st.linearly_stable) return a.st.linearly_stable > b.st.linearly_stable;
        if (a.st.max_lyap != b.st.max_lyap) return a.st.max_lyap < b.st.max_lyap;
        return a.T < b.T;
    });

    std::cout << std::setprecision(12);
    std::cout << "FOUND " << found.size() << " PERIODIC ORBITS (m3=" << m3 << "):\n";
    for (size_t i = 0; i < found.size(); ++i) {
        const auto& o = found[i];
        std::cout << "#" << (i + 1)
                  << " m3=" << o.m3
                  << " vx=" << o.vx
                  << " vy=" << o.vy
                  << " T=" << o.T
                  << " E=" << o.E
                  << " dist=" << o.dist
                  << " rmin=" << o.rmin
                  << " stable=" << (o.st.linearly_stable ? "YES" : "NO")
                  << " k1=" << (o.st.real_k ? o.st.k1 : o.st.k_real)
                  << " k2=" << (o.st.real_k ? o.st.k2 : o.st.k_imag)
                  << (o.st.real_k ? " (real)" : " (cplx)")
                  << " lyap=" << o.st.max_lyap
                  << " drift50T=" << o.st.perturb_drift_50T
                  << "\n";
    }
    return 0;
}
