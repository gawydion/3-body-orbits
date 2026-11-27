#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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
        if (rmin_local < 1e-6) return false;
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
            h = std::min(h * fac, 0.05);
            return true;
        }
        h *= std::max(0.1, 0.9 * std::pow(err, -0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate_to(State& y, double T, const Masses& m, double tol = 1e-13,
                         int max_steps = 400000, double* rmin_out = nullptr) {
    double t = 0.0, h = std::min(1e-3, T * 0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
        for (int i = 0; i < 6; ++i) if (std::abs(y[i]) > 100.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return (t >= T - 1e-13);
}

struct OrbitSpec {
    std::string id;
    std::string name;
    std::string family;
    Masses m;
    State y0;
    double T;
};

int main(int argc, char** argv) {
    // Reads orbit specs from stdin, one per line:
    // id|name|family|m1 m2 m3 x1 y1 x2 y2 x3 y3 vx1 vy1 vx2 vy2 vx3 vy3 T
    std::vector<OrbitSpec> specs;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        OrbitSpec sp;
        std::string rest;
        std::getline(ss, sp.id, '|');
        std::getline(ss, sp.name, '|');
        std::getline(ss, sp.family, '|');
        std::getline(ss, rest);
        std::stringstream rs(rest);
        rs >> sp.m[0] >> sp.m[1] >> sp.m[2];
        for (int i = 0; i < 12; ++i) rs >> sp.y0[i];
        rs >> sp.T;
        specs.push_back(sp);
    }

    int N_samples = 480;
    std::cout << std::setprecision(10) << "[\n";
    for (size_t idx = 0; idx < specs.size(); ++idx) {
        const auto& sp = specs[idx];
        State yT = sp.y0;
        double rmin = 1e300;
        integrate_to(yT, sp.T, sp.m, 1e-13, 500000, &rmin);
        double ret_dist = state_dist(yT, sp.y0);
        double E = compute_energy(sp.y0, sp.m);
        double L = compute_L(sp.y0, sp.m);
        double dE = std::abs(compute_energy(yT, sp.m) - E);
        double dL = std::abs(compute_L(yT, sp.m) - L);

        // Monodromy matrix & Floquet stability
        double M[12][12];
        const double eps = 1e-7;
        for (int col = 0; col < 12; ++col) {
            State yp = sp.y0, ym = sp.y0;
            yp[col] += eps; ym[col] -= eps;
            integrate_to(yp, sp.T, sp.m, 1e-13, 500000);
            integrate_to(ym, sp.T, sp.m, 1e-13, 500000);
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
        double k1 = 0.0, k2 = 0.0;
        bool stable = false;
        if (disc >= -1e-5) {
            double sd = std::sqrt(std::max(0.0, disc));
            k1 = 0.5 * (s + sd);
            k2 = 0.5 * (s - sd);
            stable = (std::abs(k1) < 2.0 - 1e-4) && (std::abs(k2) < 2.0 - 1e-4);
        } else {
            k1 = 0.5 * s;
            k2 = 0.5 * std::sqrt(-disc);
        }

        // 100-period empirical perturbation test (1e-5 velocity kick preserving momentum)
        State yp = sp.y0;
        yp[6] += 1e-5;
        yp[8] -= 1e-5 * sp.m[0] / sp.m[1];
        double max_drift = 0.0;
        int survived_periods = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate_to(yp, sp.T, sp.m, 1e-11, 250000)) break;
            double d = state_dist(yp, sp.y0);
            if (d > max_drift) max_drift = d;
            if (d > 5.0) break;
            survived_periods = k + 1;
        }

        // Sample full period trajectory
        std::vector<std::array<double, 6>> pts;
        pts.reserve(N_samples + 1);
        State yc = sp.y0;
        pts.push_back({yc[0], yc[1], yc[2], yc[3], yc[4], yc[5]});
        double dt = sp.T / N_samples;
        for (int k = 1; k <= N_samples; ++k) {
            integrate_to(yc, dt, sp.m, 1e-12, 100000);
            pts.push_back({yc[0], yc[1], yc[2], yc[3], yc[4], yc[5]});
        }

        std::cout << "  {\n";
        std::cout << "    \"id\": \"" << sp.id << "\",\n";
        std::cout << "    \"name\": \"" << sp.name << "\",\n";
        std::cout << "    \"family\": \"" << sp.family << "\",\n";
        std::cout << "    \"masses\": [" << sp.m[0] << ", " << sp.m[1] << ", " << sp.m[2] << "],\n";
        std::cout << "    \"y0\": [";
        for (int i = 0; i < 12; ++i) std::cout << sp.y0[i] << (i < 11 ? ", " : "");
        std::cout << "],\n";
        std::cout << "    \"T\": " << sp.T << ",\n";
        std::cout << "    \"E\": " << E << ",\n";
        std::cout << "    \"L\": " << L << ",\n";
        std::cout << "    \"rmin\": " << rmin << ",\n";
        std::cout << "    \"return_err\": " << ret_dist << ",\n";
        std::cout << "    \"dE\": " << dE << ", \"dL\": " << dL << ",\n";
        std::cout << "    \"k1\": " << k1 << ", \"k2\": " << k2 << ",\n";
        std::cout << "    \"stable\": " << (stable ? "true" : "false") << ",\n";
        std::cout << "    \"survived_100T\": " << survived_periods << ",\n";
        std::cout << "    \"max_drift_100T\": " << max_drift << ",\n";
        std::cout << "    \"pts\": [";
        std::cout << std::setprecision(5);
        for (size_t k = 0; k < pts.size(); ++k) {
            std::cout << "[" << pts[k][0] << "," << pts[k][1] << ","
                      << pts[k][2] << "," << pts[k][3] << ","
                      << pts[k][4] << "," << pts[k][5] << "]"
                      << (k + 1 < pts.size() ? "," : "");
        }
        std::cout << std::setprecision(10) << "]\n";
        std::cout << "  }" << (idx + 1 < specs.size() ? "," : "") << "\n";
    }
    std::cout << "]\n";
    return 0;
}
