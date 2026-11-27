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

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0};
    double rmin2 = 1e300;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]     - s[3*i];
            double dy = s[3*j + 1] - s[3*i + 1];
            double dz = s[3*j + 2] - s[3*i + 2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double inv_r = 1.0 / std::sqrt(r2);
            double inv_r3 = inv_r * inv_r * inv_r;
            double fx = dx * inv_r3, fy = dy * inv_r3, fz = dz * inv_r3;
            a[3*i]     += m[j] * fx; a[3*i + 1] += m[j] * fy; a[3*i + 2] += m[j] * fz;
            a[3*j]     -= m[i] * fx; a[3*j + 1] -= m[i] * fy; a[3*j + 2] -= m[i] * fz;
        }
    }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double energy3d(const State3D& s, const Masses& m) {
    double K = 0.0, U = 0.0;
    for (int i = 0; i < 3; ++i) {
        double vx = s[9 + 3*i], vy = s[10 + 3*i], vz = s[11 + 3*i];
        K += 0.5 * m[i] * (vx*vx + vy*vy + vz*vz);
    }
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j] - s[3*i], dy = s[3*j+1] - s[3*i+1], dz = s[3*j+2] - s[3*i+2];
            U -= m[i] * m[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    return K + U;
}

static inline std::array<double, 3> angmom3d(const State3D& s, const Masses& m) {
    double Lx = 0, Ly = 0, Lz = 0;
    for (int i = 0; i < 3; ++i) {
        double x = s[3*i], y = s[3*i+1], z = s[3*i+2];
        double vx = s[9+3*i], vy = s[10+3*i], vz = s[11+3*i];
        Lx += m[i] * (y*vz - z*vy);
        Ly += m[i] * (z*vx - x*vz);
        Lz += m[i] * (x*vy - y*vx);
    }
    return {Lx, Ly, Lz};
}

static inline double state_dist3d(const State3D& a, const State3D& b) {
    double sum = 0.0;
    for (int i = 0; i < 18; ++i) { double d = a[i] - b[i]; sum += d * d; }
    return std::sqrt(sum);
}

static bool dopri5_step3d(State3D& y, double& t, double& h, double t_end, const Masses& m,
                          double tol, double& rmin_track, double r_abort = 1e-5) {
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

    State3D k1, k2, k3, k4, k5, k6, k7, yt, ynew;
    double rmin_local = 1e300;
    for (int attempt = 0; attempt < 25; ++attempt) {
        deriv3d(y, m, k1, &rmin_local);
        if (rmin_local < r_abort) return false;
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
        for (int i = 0; i < 18; ++i) ynew[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        deriv3d(ynew, m, k7, &rmin_local);
        if (rmin_local < r_abort) return false;

        double err = 0.0;
        for (int i = 0; i < 18; ++i) {
            double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            double sc = tol * (1.0 + std::max(std::abs(y[i]), std::abs(ynew[i])));
            double ratio = ei / sc;
            err += ratio * ratio;
        }
        err = std::sqrt(err / 18.0);
        if (err <= 1.0) {
            y = ynew; t += h;
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

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-12,
                        int max_steps = 120000, double* rmin_out = nullptr, double r_abort = 1e-5) {
    double t = 0.0, h = std::min(1e-3, T * 0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step3d(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i = 0; i < 9; ++i) if (std::abs(y[i]) > 80.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

// Compute smallest eigenvalue of 3x3 symmetric covariance matrix to measure true 3D thickness
static double min_eig_3x3(double C[3][3]) {
    // Characteristic polynomial det(C - lam I) = 0 -> lam^3 - I1 lam^2 + I2 lam - I3 = 0
    double I1 = C[0][0] + C[1][1] + C[2][2];
    double I2 = C[0][0]*C[1][1] + C[0][0]*C[2][2] + C[1][1]*C[2][2]
              - C[0][1]*C[0][1] - C[0][2]*C[0][2] - C[1][2]*C[1][2];
    double I3 = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[1][2])
              - C[0][1]*(C[0][1]*C[2][2] - C[1][2]*C[0][2])
              + C[0][2]*(C[0][1]*C[1][2] - C[1][1]*C[0][2]);
    double p = I2 - I1 * I1 / 3.0;
    double q = -2.0 * I1 * I1 * I1 / 27.0 + I1 * I2 / 3.0 - I3;
    double r = std::sqrt(std::max(0.0, -p / 3.0));
    if (r < 1e-15) return std::max(0.0, I1 / 3.0);
    double c = std::clamp(-q / (2.0 * r * r * r), -1.0, 1.0);
    double phi = std::acos(c) / 3.0;
    double l1 = I1 / 3.0 + 2.0 * r * std::cos(phi);
    double l2 = I1 / 3.0 + 2.0 * r * std::cos(phi + 2.0 * M_PI / 3.0);
    double l3 = I1 / 3.0 + 2.0 * r * std::cos(phi + 4.0 * M_PI / 3.0);
    return std::max(0.0, std::min({l1, l2, l3}));
}

struct Topo3D {
    int partner_swaps;
    double rmin, rmax;
    double thickness_3d; // RMS out-of-best-fit-plane distance
    bool democratic;
};

static Topo3D analyze_topo3d(const State3D& y0, double T, const Masses& m) {
    const int N = 320;
    double dt = T / N;
    State3D y = y0;
    int prev_cl = -1, swaps = 0;
    double rmin = 1e300, rmax = 0.0;
    double C[3][3] = {{0}};
    int count = 0;

    for (int k = 0; k <= N; ++k) {
        double r12 = std::sqrt((y[3]-y[0])*(y[3]-y[0]) + (y[4]-y[1])*(y[4]-y[1]) + (y[5]-y[2])*(y[5]-y[2]));
        double r13 = std::sqrt((y[6]-y[0])*(y[6]-y[0]) + (y[7]-y[1])*(y[7]-y[1]) + (y[8]-y[2])*(y[8]-y[2]));
        double r23 = std::sqrt((y[6]-y[3])*(y[6]-y[3]) + (y[7]-y[4])*(y[7]-y[4]) + (y[8]-y[5])*(y[8]-y[5]));
        rmin = std::min({rmin, r12, r13, r23});
        rmax = std::max({rmax, r12, r13, r23});
        int cl = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
        if (prev_cl != -1 && cl != prev_cl) swaps++;
        prev_cl = cl;

        for (int b = 0; b < 3; ++b) {
            double p[3] = {y[3*b], y[3*b+1], y[3*b+2]};
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    C[r][c] += p[r] * p[c];
            count++;
        }
        if (k < N) integrate3d(y, dt, m, 1e-11, 50000);
    }
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            C[r][c] /= count;
    double thick = std::sqrt(min_eig_3x3(C));
    return {swaps, rmin, rmax, thick, swaps >= 2};
}

struct Stability3D {
    bool linearly_stable;
    double k1, k2, k3;
    double stab_margin;
    int survived_100T;
    double drift_100T;
};

static Stability3D analyze_stability3d(const State3D& y0, double T, const Masses& m) {
    double M[18][18];
    const double eps = 1e-7;
    for (int col = 0; col < 18; ++col) {
        State3D yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        integrate3d(yp, T, m, 1e-13, 300000);
        integrate3d(ym, T, m, 1e-13, 300000);
        for (int row = 0; row < 18; ++row)
            M[row][col] = (yp[row] - ym[row]) / (2.0 * eps);
    }

    double M2[18][18] = {{0}};
    for (int i = 0; i < 18; ++i)
        for (int k = 0; k < 18; ++k)
            for (int j = 0; j < 18; ++j)
                M2[i][j] += M[i][k] * M[k][j];

    double tr1 = 0.0, tr2 = 0.0, tr3 = 0.0;
    for (int i = 0; i < 18; ++i) {
        tr1 += M[i][i];
        tr2 += M2[i][i];
        for (int j = 0; j < 18; ++j) tr3 += M2[i][j] * M[j][i];
    }

    // 12 trivial unit eigenvalues (for L != 0):
    // p1 = k1 + k2 + k3 = tr1 - 12
    // p2 = k1^2 + k2^2 + k3^2 = tr2 - 6
    // p3 = k1^3 + k2^3 + k3^3 = tr3 - 12 + 3*p1
    double p1 = tr1 - 12.0;
    double p2 = tr2 - 6.0;
    double p3 = tr3 - 12.0 + 3.0 * p1;

    double e1 = p1;
    double e2 = 0.5 * (e1 * p1 - p2);
    double e3 = (e2 * p1 - e1 * p2 + p3) / 3.0;

    // Solve k^3 - e1 k^2 + e2 k - e3 = 0
    double Q = (e1 * e1 - 3.0 * e2) / 9.0;
    double R = (2.0 * e1 * e1 * e1 - 9.0 * e1 * e2 + 27.0 * e3) / 54.0;
    double disc = Q * Q * Q - R * R;

    Stability3D st{};
    st.stab_margin = -1.0;
    if (disc >= -1e-4 && Q > 0.0) {
        double sqQ = std::sqrt(Q);
        double arg = std::clamp(R / (sqQ * sqQ * sqQ), -1.0, 1.0);
        double th = std::acos(arg);
        st.k1 = 2.0 * sqQ * std::cos(th / 3.0) + e1 / 3.0;
        st.k2 = 2.0 * sqQ * std::cos((th + 2.0 * M_PI) / 3.0) + e1 / 3.0;
        st.k3 = 2.0 * sqQ * std::cos((th + 4.0 * M_PI) / 3.0) + e1 / 3.0;
        double mk = std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
        st.stab_margin = 2.0 - mk;
        st.linearly_stable = (mk < 2.0 - 1e-4);
    } else {
        st.linearly_stable = false;
    }

    // Empirical 100-period full-3D perturbation test (kick in x, y, AND z!)
    if (st.linearly_stable) {
        State3D yp = y0;
        double kick = 1e-5;
        yp[9]  += kick; yp[12] -= kick * m[0] / m[1];
        yp[10] += kick; yp[13] -= kick * m[0] / m[1];
        yp[11] += kick; yp[14] -= kick * m[0] / m[1];
        double max_drift = 0.0;
        int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate3d(yp, T, m, 1e-11, 200000)) break;
            double d = state_dist3d(yp, y0);
            if (d > max_drift) max_drift = d;
            if (d > 5.0) break;
            surv = k + 1;
        }
        st.survived_100T = surv;
        st.drift_100T = max_drift;
    }
    return st;
}

static void emit_json3d(const std::string& subspace, const Masses& m, const State3D& y0, double T,
                        double ret_err, const Stability3D& st, const Topo3D& tp) {
    double E = energy3d(y0, m);
    auto L = angmom3d(y0, m);
    double Lmag = std::sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"dim\":3,\"subspace\":\"" << subspace << "\""
       << ",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "]"
       << ",\"y0\":[";
    for (int i = 0; i < 18; ++i) os << y0[i] << (i < 17 ? "," : "");
    os << "],\"T\":" << T
       << ",\"E\":" << E
       << ",\"L\":" << Lmag
       << ",\"Lvec\":[" << L[0] << "," << L[1] << "," << L[2] << "]"
       << ",\"ret_err\":" << ret_err
       << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2 << ",\"k3\":" << st.k3
       << ",\"stab_margin\":" << st.stab_margin
       << ",\"survived_100T\":" << st.survived_100T
       << ",\"drift_100T\":" << st.drift_100T
       << ",\"rmin\":" << tp.rmin
       << ",\"rmax\":" << tp.rmax
       << ",\"thickness_3d\":" << tp.thickness_3d
       << ",\"partner_swaps\":" << tp.partner_swaps
       << ",\"democratic\":" << (tp.democratic ? "true" : "false")
       << "}\n";
    std::cout << os.str() << std::flush;
}

// ============================================================
// MODE 3A: 3D Perpendicular-Line-Crossing (6 vars, 6 eqs)
// ============================================================
static State3D make_perp3d(double x2, double vy1, double vz1, double vy2, double vz2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2) / m[2];
    double vy3 = -(m[0]*vy1 + m[1]*vy2) / m[2];
    double vz3 = -(m[0]*vz1 + m[1]*vz2) / m[2];
    return {x1, 0, 0, x2, 0, 0, x3, 0, 0,
            0, vy1, vz1, 0, vy2, vz2, 0, vy3, vz3};
}

static bool refine_perp3d(double p[6], const Masses& m, double& ferr) {
    // p = {x2, vy1, vz1, vy2, vz2, Th}
    double lambda = 1e-3;
    auto eval_F = [&](const double q[6], std::array<double, 6>& F) -> bool {
        if (q[5] < 0.4 || q[5] > 22.0) return false;
        State3D y = make_perp3d(q[0], q[1], q[2], q[3], q[4], m);
        double d12 = std::abs(y[3]-y[0]), d13 = std::abs(y[6]-y[0]), d23 = std::abs(y[6]-y[3]);
        if (std::min({d12, d13, d23}) < 0.12) return false;
        double rm = 1e300;
        if (!integrate3d(y, q[5], m, 1e-11, 35000, &rm, 1e-2)) return false;
        // F = (y1, z1, y2, z2, vx1, vx2)
        F = {y[1], y[2], y[4], y[5], y[9], y[12]};
        return true;
    };

    std::array<double, 6> F;
    if (!eval_F(p, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;

    for (int iter = 0; iter < 30; ++iter) {
        if (std::sqrt(cost) < 1e-10) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[6][6];
        std::array<double, 6> Fp, Fm;
        for (int col = 0; col < 6; ++col) {
            double pp[6], pm[6];
            for (int i = 0; i < 6; ++i) { pp[i] = p[i]; pm[i] = p[i]; }
            pp[col] += eps; pm[col] -= eps;
            if (!eval_F(pp, Fp) || !eval_F(pm, Fm)) return false;
            for (int r = 0; r < 6; ++r) J[r][col] = (Fp[r] - Fm[r]) / (2.0 * eps);
        }
        double A[6][7] = {{0}};
        for (int r = 0; r < 6; ++r) {
            for (int c = 0; c < 6; ++c) {
                double s = 0; for (int i = 0; i < 6; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0; for (int i = 0; i < 6; ++i) s -= J[i][r] * F[i];
            A[r][6] = s;
        }
        bool sing = false;
        for (int col = 0; col < 6; ++col) {
            int piv = col;
            for (int r = col + 1; r < 6; ++r)
                if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (std::abs(A[piv][col]) < 1e-18) { sing = true; break; }
            for (int c = col; c < 7; ++c) std::swap(A[col][c], A[piv][c]);
            for (int r = col + 1; r < 6; ++r) {
                double fac = A[r][col] / A[col][col];
                for (int c = col; c < 7; ++c) A[r][c] -= fac * A[col][c];
            }
        }
        if (sing) { lambda *= 10.0; continue; }
        double d[6] = {0};
        for (int r = 5; r >= 0; --r) {
            double s = A[r][6];
            for (int c = r + 1; c < 6; ++c) s -= A[r][c] * d[c];
            d[r] = s / A[r][r];
        }
        double sn = 0; for (int i = 0; i < 6; ++i) sn += d[i]*d[i];
        sn = std::sqrt(sn);
        if (sn > 0.35) for (int i = 0; i < 6; ++i) d[i] *= (0.35 / sn);
        double np[6]; for (int i = 0; i < 6; ++i) np[i] = p[i] + d[i];
        std::array<double, 6> nF;
        if (!eval_F(np, nF)) { lambda *= 10.0; continue; }
        double ncost = 0; for (double v : nF) ncost += v*v;
        if (ncost < cost) {
            for (int i = 0; i < 6; ++i) p[i] = np[i];
            F = nF; cost = ncost; lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0; if (lambda > 1e7) break;
        }
    }
    ferr = std::sqrt(cost);
    return (ferr < 1e-9 && p[5] > 0.5);
}

static void run_perp3d(const Masses& m, int nsamples, double Thmax, int seed_offset) {
    struct Seed { double p[6]; double err; int swaps; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(991001 + seed_offset * 10007 + tid * 9973);
            std::uniform_real_distribution<double> ux2(-0.6, 2.0);
            std::uniform_real_distribution<double> uv(-1.2, 1.2);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double x2 = ux2(rng), vy1 = uv(rng), vz1 = uv(rng), vy2 = uv(rng), vz2 = uv(rng);
                // Require genuinely non-planar initial velocities!
                if (std::abs(vy1 * vz2 - vz1 * vy2) < 0.15) continue;
                State3D y0 = make_perp3d(x2, vy1, vz1, vy2, vz2, m);
                double d12 = std::abs(y0[3]-y0[0]), d13 = std::abs(y0[6]-y0[0]), d23 = std::abs(y0[6]-y0[3]);
                if (std::min({d12, d13, d23}) < 0.25) continue;
                if (energy3d(y0, m) >= -0.05) continue;

                State3D y = y0;
                double t = 0, h = 1e-3, rmin = 1e300;
                double best_e = 1e9, best_t = 0, pe = 1e9, p2e = 1e9, pt = 0;
                int prev_cl = -1, swaps = 0, best_swaps = 0, steps = 0;
                while (t < Thmax && steps < 12000) {
                    double hc = std::min(h, 0.03);
                    if (!dopri5_step3d(y, t, hc, Thmax, m, 1e-8, rmin, 1e-2)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i = 0; i < 9; ++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 18.0) break;

                    double r12 = std::sqrt((y[3]-y[0])*(y[3]-y[0])+(y[4]-y[1])*(y[4]-y[1])+(y[5]-y[2])*(y[5]-y[2]));
                    double r13 = std::sqrt((y[6]-y[0])*(y[6]-y[0])+(y[7]-y[1])*(y[7]-y[1])+(y[8]-y[2])*(y[8]-y[2]));
                    double r23 = std::sqrt((y[6]-y[3])*(y[6]-y[3])+(y[7]-y[4])*(y[7]-y[4])+(y[8]-y[5])*(y[8]-y[5]));
                    int cl = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
                    if (prev_cl != -1 && cl != prev_cl) swaps++;
                    prev_cl = cl;

                    double e = std::sqrt(y[1]*y[1] + y[2]*y[2] + y[4]*y[4] + y[5]*y[5] + y[9]*y[9] + y[12]*y[12]);
                    if (t > 0.5 && steps >= 3 && pe < p2e && pe < e && pe < best_e) {
                        best_e = pe; best_t = pt; best_swaps = swaps;
                    }
                    p2e = pe; pe = e; pt = t;
                }
                if (best_e < 0.35 && best_t > 0.5) {
                    tseeds[tid].push_back({{x2, vy1, vz1, vy2, vz2, best_t}, best_e, best_swaps});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all_seeds;
    for (auto& v : tseeds) all_seeds.insert(all_seeds.end(), v.begin(), v.end());
    std::sort(all_seeds.begin(), all_seeds.end(), [](const Seed& a, const Seed& b) {
        if ((a.swaps > 0) != (b.swaps > 0)) return (a.swaps > 0) > (b.swaps > 0);
        return a.err < b.err;
    });
    size_t limit = std::min(all_seeds.size(), (size_t)350);
    std::cerr << "[3d_perp m=(" << m[0] << "," << m[1] << "," << m[2] << ")] seeds=" << all_seeds.size()
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
                double p[6]; for (int i = 0; i < 6; ++i) p[i] = all_seeds[idx].p[i];
                double ferr = 0;
                if (!refine_perp3d(p, m, ferr)) continue;
                double T = 2.0 * p[5];
                State3D y0 = make_perp3d(p[0], p[1], p[2], p[3], p[4], m);
                State3D yT = y0;
                double rmin = 1e300;
                if (!integrate3d(yT, T, m, 1e-13, 250000, &rmin, 5e-3)) continue;
                double ret = state_dist3d(yT, y0);
                if (ret > 1e-7) continue;

                // Reject tilted planar orbits: velocities in yz-plane must be non-collinear
                if (std::abs(p[1] * p[4] - p[2] * p[3]) < 1e-2) continue;
                auto tp = analyze_topo3d(y0, T, m);
                if (tp.thickness_3d < 0.05) continue; // Reject flat 2D orbits!

                double E = energy3d(y0, m);
                auto Lv = angmom3d(y0, m);
                double Lmag = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
                double inv1 = E * std::pow(T, 2.0 / 3.0);
                double inv2 = Lmag * std::sqrt(std::max(0.0, -E));
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup)
                        if (std::abs(f[0] - inv1) < 2e-2 && std::abs(f[1] - inv2) < 2e-2) { dup = true; break; }
                    if (dup) continue;
                    dedup.push_back({inv1, inv2});
                }

                auto st = analyze_stability3d(y0, T, m);
                if (!st.linearly_stable || st.survived_100T < 100) continue;
                std::lock_guard<std::mutex> lk(mtx);
                emit_json3d("3d_perp", m, y0, T, ret, st, tp);
            }
        });
    }
    for (auto& th : ths) th.join();
}

// ============================================================
// MODE 3B: 3D Symmetric Collinear (m1=m2=1, m3=m)
// r1=(-1,0,0), r2=(1,0,0), r3=(0,0,0)
// v1=(vx,vy,vz), v2=(vx,vy,-vz), v3=(-2vx/m, -2vy/m, 0)
// ============================================================
static State3D make_col3d(double vx, double vy, double vz, double m3) {
    return {-1, 0, 0, 1, 0, 0, 0, 0, 0,
            vx, vy, vz, vx, vy, -vz, -2*vx/m3, -2*vy/m3, 0};
}

static bool refine_col3d(double p[4], double m3, double& ferr) {
    // p = {vx, vy, vz, T}
    Masses m = {1.0, 1.0, m3};
    double lambda = 1e-4;
    auto eval_F = [&](const double q[4], State3D& diff) -> bool {
        if (q[3] < 1.0 || q[3] > 40.0) return false;
        State3D y0 = make_col3d(q[0], q[1], q[2], m3);
        State3D yT = y0;
        double rm = 1e300;
        if (!integrate3d(yT, q[3], m, 1e-12, 100000, &rm, 1e-3)) return false;
        for (int i = 0; i < 18; ++i) diff[i] = yT[i] - y0[i];
        return true;
    };
    State3D F;
    if (!eval_F(p, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;

    for (int iter = 0; iter < 30; ++iter) {
        if (std::sqrt(cost) < 1e-10) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[18][4];
        State3D Fp, Fm;
        for (int col = 0; col < 3; ++col) {
            double pp[4] = {p[0], p[1], p[2], p[3]}, pm[4] = {p[0], p[1], p[2], p[3]};
            pp[col] += eps; pm[col] -= eps;
            if (!eval_F(pp, Fp) || !eval_F(pm, Fm)) return false;
            for (int r = 0; r < 18; ++r) J[r][col] = (Fp[r] - Fm[r]) / (2.0 * eps);
        }
        State3D y0 = make_col3d(p[0], p[1], p[2], m3), yT, fT;
        for (int i = 0; i < 18; ++i) yT[i] = y0[i] + F[i];
        deriv3d(yT, m, fT);
        for (int r = 0; r < 18; ++r) J[r][3] = fT[r];

        double A[4][5] = {{0}};
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                double s = 0; for (int i = 0; i < 18; ++i) s += J[i][r] * J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda * (1.0 + A[r][r]);
            double s = 0; for (int i = 0; i < 18; ++i) s -= J[i][r] * F[i];
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
        double np[4] = {p[0]+d[0], p[1]+d[1], p[2]+d[2], p[3]+d[3]};
        State3D nF;
        if (!eval_F(np, nF)) { lambda *= 10.0; continue; }
        double ncost = 0; for (double v : nF) ncost += v*v;
        if (ncost < cost) {
            for (int i = 0; i < 4; ++i) p[i] = np[i];
            F = nF; cost = ncost; lambda = std::max(1e-12, lambda * 0.25);
        } else {
            lambda *= 5.0; if (lambda > 1e7) break;
        }
    }
    ferr = std::sqrt(cost);
    return (ferr < 1e-8 && p[3] > 1.5);
}

static void run_col3d(double m3, int nsamples, double Tmax, int seed_offset) {
    Masses m = {1.0, 1.0, m3};
    struct Seed { double p[4]; double err; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(551001 + seed_offset * 10007 + tid * 9973);
            std::uniform_real_distribution<double> uvx(0.0, 0.85);
            std::uniform_real_distribution<double> uvy(0.1, 0.85);
            std::uniform_real_distribution<double> uvz(0.1, 0.85);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double vx = uvx(rng), vy = uvy(rng), vz = uvz(rng);
                State3D y0 = make_col3d(vx, vy, vz, m3);
                if (energy3d(y0, m) >= -0.05) continue;
                State3D y = y0;
                double t = 0, h = 1e-3, rmin = 1e300;
                double best_d = 1e9, best_t = 0, pd = 0, p2d = 0, pt = 0;
                bool left = false;
                int steps = 0;
                while (t < Tmax && steps < 18000) {
                    double hc = std::min(h, 0.03);
                    if (!dopri5_step3d(y, t, hc, Tmax, m, 1e-8, rmin, 5e-3)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i = 0; i < 9; ++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 20.0) break;
                    double d = state_dist3d(y, y0);
                    if (t > 1.0 && d > 0.5) left = true;
                    if (left && steps >= 3 && pd < p2d && pd < d && pd < best_d) {
                        best_d = pd; best_t = pt;
                    }
                    p2d = pd; pd = d; pt = t;
                }
                if (best_d < 0.4 && best_t > 1.5) {
                    tseeds[tid].push_back({{vx, vy, vz, best_t}, best_d});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all_seeds;
    for (auto& v : tseeds) all_seeds.insert(all_seeds.end(), v.begin(), v.end());
    std::sort(all_seeds.begin(), all_seeds.end(), [](const Seed& a, const Seed& b){ return a.err < b.err; });
    size_t limit = std::min(all_seeds.size(), (size_t)300);
    std::cerr << "[3d_collinear m3=" << m3 << "] seeds=" << all_seeds.size() << " refining=" << limit << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 4>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                double p[4] = {all_seeds[idx].p[0], all_seeds[idx].p[1], all_seeds[idx].p[2], all_seeds[idx].p[3]};
                double ferr = 0;
                if (!refine_col3d(p, m3, ferr)) continue;
                State3D y0 = make_col3d(p[0], p[1], p[2], m3);
                auto tp = analyze_topo3d(y0, p[3], m);
                if (tp.thickness_3d < 0.05) continue;

                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup) {
                        double dv = std::sqrt((f[0]-p[0])*(f[0]-p[0]) + (f[1]-p[1])*(f[1]-p[1]) + (f[2]-p[2])*(f[2]-p[2]));
                        if (dv < 5e-3 && std::abs(f[3]-p[3]) < 0.1) { dup = true; break; }
                    }
                    if (dup) continue;
                    dedup.push_back({p[0], p[1], p[2], p[3]});
                }

                auto st = analyze_stability3d(y0, p[3], m);
                if (!st.linearly_stable || st.survived_100T < 100) continue;
                std::lock_guard<std::mutex> lk(mtx);
                emit_json3d("3d_collinear", m, y0, p[3], ferr, st, tp);
            }
        });
    }
    for (auto& th : ths) th.join();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n  campaign_3d --perp <m1> <m2> <m3> <nsamples> <Thmax> [seed]\n"
                  << "  campaign_3d --collinear <m3> <nsamples> <Tmax> [seed]\n";
        return 1;
    }
    std::string mode = argv[1];
    if (mode == "--perp") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        int nsamples = (argc > 5) ? std::atoi(argv[5]) : 25000;
        double Thmax = (argc > 6) ? std::atof(argv[6]) : 11.0;
        int seed = (argc > 7) ? std::atoi(argv[7]) : 0;
        run_perp3d(m, nsamples, Thmax, seed);
    } else if (mode == "--collinear") {
        double m3 = std::atof(argv[2]);
        int nsamples = (argc > 3) ? std::atoi(argv[3]) : 30000;
        double Tmax = (argc > 4) ? std::atof(argv[4]) : 22.0;
        int seed = (argc > 5) ? std::atoi(argv[5]) : 0;
        run_col3d(m3, nsamples, Tmax, seed);
    }
    return 0;
}
