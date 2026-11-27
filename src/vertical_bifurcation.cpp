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

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0}, rmin2 = 1e300;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]   += m[j]*dx*ir3; a[3*i+1] += m[j]*dy*ir3; a[3*i+2] += m[j]*dz*ir3;
            a[3*j]   -= m[i]*dx*ir3; a[3*j+1] -= m[i]*dy*ir3; a[3*j+2] -= m[i]*dz*ir3;
        }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static bool dopri5_step3d(State3D& y, double& t, double& h, double t_end, const Masses& m, double tol) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;
    State3D k1,k2,k3,k4,k5,k6,k7,yt,yn;
    double rm = 1e300;
    for (int att = 0; att < 25; ++att) {
        deriv3d(y, m, k1, &rm); if (rm < 1e-6) return false;
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt, m, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt, m, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt, m, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt, m, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt, m, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv3d(yn, m, k7, &rm); if (rm < 1e-6) return false;
        double err = 0;
        for (int i=0;i<18;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err/18.0);
        if (err <= 1.0) {
            y = yn; t += h;
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13) {
    double t = 0, h = std::min(1e-3, T*0.01);
    int steps = 0;
    while (t < T - 1e-15 && steps < 400000) {
        if (!dopri5_step3d(y, t, h, T, m, tol)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 80.0) return false;
    }
    return t >= T - 1e-13;
}

// Compute vertical (out-of-plane z) stability index kz of a planar orbit!
// The 6 z-variables are indices {2, 5, 8, 11, 14, 17}.
// 4 of the 6 eigenvalues are trivial +1 (z_cm, P_z, R_x, R_y), so tr(M_z) = 4 + kz!
static double compute_kz(const State3D& y0, double T, const Masses& m) {
    const int zidx[6] = {2, 5, 8, 11, 14, 17};
    double tr = 0.0;
    const double eps = 1e-7;
    for (int k = 0; k < 6; ++k) {
        int col = zidx[k];
        State3D yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        integrate3d(yp, T, m, 1e-13);
        integrate3d(ym, T, m, 1e-13);
        tr += (yp[col] - ym[col]) / (2.0 * eps);
    }
    return tr - 4.0;
}

int main() {
    // Reads lines: id m1 m2 m3 x1 y1 x2 y2 x3 y3 vx1 vy1 vx2 vy2 vx3 vy3 T
    std::string line;
    std::cout << std::setprecision(10);
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string id;
        Masses m;
        double y2d[12], T;
        ss >> id >> m[0] >> m[1] >> m[2];
        for (int i = 0; i < 12; ++i) ss >> y2d[i];
        ss >> T;
        State3D y0 = {
            y2d[0], y2d[1], 0,  y2d[2], y2d[3], 0,  y2d[4], y2d[5], 0,
            y2d[6], y2d[7], 0,  y2d[8], y2d[9], 0,  y2d[10], y2d[11], 0
        };
        double kz = compute_kz(y0, T, m);
        bool stable_3d = (std::abs(kz) < 2.0 - 1e-4);
        std::cout << id << " " << kz << " " << (stable_3d ? 1 : 0) << "\n" << std::flush;
    }
    return 0;
}
