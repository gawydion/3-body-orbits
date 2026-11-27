#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using State2D = std::array<double, 12>;
using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0}, rmin2 = 1e300;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j] - s[3*i], dy = s[3*j+1] - s[3*i+1], dz = s[3*j+2] - s[3*i+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]   += m[j]*dx*ir3; a[3*i+1] += m[j]*dy*ir3; a[3*i+2] += m[j]*dz*ir3;
            a[3*j]   -= m[i]*dx*ir3; a[3*j+1] -= m[i]*dy*ir3; a[3*j+2] -= m[i]*dz*ir3;
        }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double energy3d(const State3D& s, const Masses& m) {
    double K = 0, U = 0;
    for (int i = 0; i < 3; ++i)
        K += 0.5*m[i]*(s[9+3*i]*s[9+3*i] + s[10+3*i]*s[10+3*i] + s[11+3*i]*s[11+3*i]);
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            U -= m[i]*m[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    return K + U;
}

static inline std::array<double, 3> angmom3d(const State3D& s, const Masses& m) {
    double Lx = 0, Ly = 0, Lz = 0;
    for (int i = 0; i < 3; ++i) {
        Lx += m[i]*(s[3*i+1]*s[11+3*i] - s[3*i+2]*s[10+3*i]);
        Ly += m[i]*(s[3*i+2]*s[9+3*i]  - s[3*i]*s[11+3*i]);
        Lz += m[i]*(s[3*i]*s[10+3*i]   - s[3*i+1]*s[9+3*i]);
    }
    return {Lx, Ly, Lz};
}

static inline double state_dist3d(const State3D& a, const State3D& b) {
    double s = 0;
    for (int i = 0; i < 18; ++i) { double d = a[i] - b[i]; s += d*d; }
    return std::sqrt(s);
}

static inline void rot_vec(const double u[3], double ang, double& x, double& y, double& z) {
    double c = std::cos(ang), s = std::sin(ang), omc = 1.0 - c;
    double dot = u[0]*x + u[1]*y + u[2]*z;
    double cx = u[1]*z - u[2]*y, cy = u[2]*x - u[0]*z, cz = u[0]*y - u[1]*x;
    double nx = x*c + cx*s + u[0]*dot*omc;
    double ny = y*c + cy*s + u[1]*dot*omc;
    double nz = z*c + cz*s + u[2]*dot*omc;
    x = nx; y = ny; z = nz;
}

static State3D rotate_state(const State3D& s, const double u[3], double ang) {
    State3D out = s;
    for (int b = 0; b < 6; ++b) rot_vec(u, ang, out[3*b], out[3*b+1], out[3*b+2]);
    return out;
}

static bool dopri5_step3d(State3D& y, double& t, double& h, double t_end, const Masses& m,
                          double tol, double& rmin_track, double r_abort = 1e-4) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;
    State3D k1, k2, k3, k4, k5, k6, k7, yt, yn;
    double rm = 1e300;
    for (int att = 0; att < 25; ++att) {
        deriv3d(y, m, k1, &rm); if (rm < r_abort) return false;
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt, m, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt, m, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt, m, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt, m, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt, m, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv3d(yn, m, k7, &rm); if (rm < r_abort) return false;
        double err = 0;
        for (int i=0;i<18;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err/18.0);
        if (err <= 1.0) {
            y = yn; t += h;
            if (rm < rmin_track) rmin_track = rm;
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13,
                        int max_steps = 800000, double* rmin_out = nullptr, double r_abort = 1e-4) {
    double t = 0, h = std::min(1e-3, T*0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step3d(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 100.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

static State3D make_rpo_state(double x2, double vy1, double vy2, double vz2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2) / m[2];
    double vy3 = -(m[0]*vy1 + m[1]*vy2) / m[2];
    double vz3 = -(m[1]*vz2) / m[2];
    return {x1, 0, 0, x2, 0, 0, x3, 0, 0,
            0, vy1, 0, 0, vy2, vz2, 0, vy3, vz3};
}

// Generic LM solver for m_eq equations in n_var unknowns
template <int NEQ, int NVAR>
static bool solve_lm_step(const double J[NEQ][NVAR], const double F[NEQ], double lambda, double d[NVAR]) {
    double A[NVAR][NVAR + 1] = {{0}};
    for (int r = 0; r < NVAR; ++r) {
        for (int c = 0; c < NVAR; ++c) {
            double s = 0; for (int i = 0; i < NEQ; ++i) s += J[i][r] * J[i][c];
            A[r][c] = s;
        }
        A[r][r] += lambda * (1.0 + A[r][r]);
        double s = 0; for (int i = 0; i < NEQ; ++i) s -= J[i][r] * F[i];
        A[r][NVAR] = s;
    }
    for (int col = 0; col < NVAR; ++col) {
        int piv = col;
        for (int r = col + 1; r < NVAR; ++r)
            if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
        if (std::abs(A[piv][col]) < 1e-20) return false;
        for (int c = col; c <= NVAR; ++c) std::swap(A[col][c], A[piv][c]);
        for (int r = col + 1; r < NVAR; ++r) {
            double f = A[r][col] / A[col][col];
            for (int c = col; c <= NVAR; ++c) A[r][c] -= f * A[col][c];
        }
    }
    for (int r = NVAR - 1; r >= 0; --r) {
        double s = A[r][NVAR];
        for (int c = r + 1; c < NVAR; ++c) s -= A[r][c] * d[c];
        d[r] = s / A[r][r];
    }
    return true;
}

// Refine 2D perp_L_nonzero (4 vars: x2, vy1, vy2, Th) with alpha = 0
static bool refine_2d_perp(double p[4], const Masses& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[4], double F[4]) -> bool {
        if (q[3] < 0.3 || q[3] > 60.0) return false;
        State3D y = make_rpo_state(q[0], q[1], q[2], 0.0, m);
        if (!integrate3d(y, q[3], m, 1e-13, 300000)) return false;
        F[0] = y[1]; F[1] = y[4]; F[2] = y[9]; F[3] = y[12];
        return true;
    };
    double F[4];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<4;++i) cost += F[i]*F[i];
    for (int it = 0; it < 28; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[4][4], Fp[4], Fm[4];
        for (int col = 0; col < 4; ++col) {
            double pp[4], pm[4]; for (int i=0;i<4;++i) { pp[i]=p[i]; pm[i]=p[i]; }
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<4;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[4] = {0};
        if (!solve_lm_step<4, 4>(J, F, lambda, d)) { lambda *= 10; continue; }
        double sn = 0; for (int i=0;i<4;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.25) for (int i=0;i<4;++i) d[i] *= (0.25/sn);
        double np[4], nF[4]; for (int i=0;i<4;++i) np[i] = p[i]+d[i];
        if (!eval(np, nF)) { lambda *= 10; continue; }
        double nc = 0; for (int i=0;i<4;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for (int i=0;i<4;++i) p[i]=np[i]; for (int i=0;i<4;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13, lambda*0.2); }
        else { lambda *= 5; if (lambda > 1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Refine 2D perp RPO at fixed x2 (4 vars: vy1, vy2, Th, alpha)
static bool refine_2d_rpo_fixed_x2(double p[4], double x2, const Masses& m, double& ferr) {
    double lambda = 1e-5;
    const double uz[3] = {0, 0, 1};
    auto eval = [&](const double q[4], double F[4]) -> bool {
        if (q[2] < 0.3 || q[2] > 60.0) return false;
        State3D y = make_rpo_state(x2, q[0], q[1], 0.0, m);
        if (!integrate3d(y, q[2], m, 1e-13, 300000)) return false;
        State3D yt = rotate_state(y, uz, -q[3]);
        F[0] = yt[1]; F[1] = yt[4]; F[2] = yt[9]; F[3] = yt[12];
        return true;
    };
    double F[4];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<4;++i) cost += F[i]*F[i];
    for (int it = 0; it < 28; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[4][4], Fp[4], Fm[4];
        for (int col = 0; col < 4; ++col) {
            double pp[4], pm[4]; for (int i=0;i<4;++i) { pp[i]=p[i]; pm[i]=p[i]; }
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<4;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[4] = {0};
        if (!solve_lm_step<4, 4>(J, F, lambda, d)) { lambda *= 10; continue; }
        double sn = 0; for (int i=0;i<4;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.25) for (int i=0;i<4;++i) d[i] *= (0.25/sn);
        double np[4], nF[4]; for (int i=0;i<4;++i) np[i] = p[i]+d[i];
        if (!eval(np, nF)) { lambda *= 10; continue; }
        double nc = 0; for (int i=0;i<4;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for (int i=0;i<4;++i) p[i]=np[i]; for (int i=0;i<4;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13, lambda*0.2); }
        else { lambda *= 5; if (lambda > 1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Refine 3D RPO at fixed (vz2, m): 5 vars p = {x2, vy1, vy2, Th, alpha}
static bool refine_3d_rpo_fixed_vz(double p[5], double vz2, const Masses& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[5], double F[6]) -> bool {
        if (q[3] < 0.5 || q[3] > 80.0) return false;
        State3D y0 = make_rpo_state(q[0], q[1], q[2], vz2, m);
        double d12 = std::abs(y0[3]-y0[0]), d13 = std::abs(y0[6]-y0[0]), d23 = std::abs(y0[6]-y0[3]);
        if (std::min({d12, d13, d23}) < 0.08) return false;
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[3], m, 1e-13, 600000)) return false;
        State3D yt = rotate_state(y, u, -q[4]);
        F[0] = yt[1]; F[1] = yt[2]; F[2] = yt[4]; F[3] = yt[5]; F[4] = yt[9]; F[5] = yt[12];
        return true;
    };
    double F[6];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<6;++i) cost += F[i]*F[i];
    for (int it = 0; it < 35; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = 1e-7, J[6][5], Fp[6], Fm[6];
        for (int col = 0; col < 5; ++col) {
            double pp[5], pm[5]; for (int i=0;i<5;++i) { pp[i]=p[i]; pm[i]=p[i]; }
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5] = {0};
        if (!solve_lm_step<6, 5>(J, F, lambda, d)) { lambda *= 10; continue; }
        double sn = 0; for (int i=0;i<5;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.25) for (int i=0;i<5;++i) d[i] *= (0.25/sn);
        double np[5], nF[6]; for (int i=0;i<5;++i) np[i] = p[i]+d[i];
        if (!eval(np, nF)) { lambda *= 10; continue; }
        double nc = 0; for (int i=0;i<6;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for (int i=0;i<5;++i) p[i]=np[i]; for (int i=0;i<6;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13, lambda*0.2); }
        else { lambda *= 5; if (lambda > 1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Refine 3D RPO at fixed alpha (snapping exact rational inertial closure!): 5 vars p = {x2, vy1, vy2, vz2, Th}
static bool refine_3d_rpo_fixed_alpha(double p[5], double target_a, const Masses& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[5], double F[6]) -> bool {
        if (q[4] < 0.5 || q[4] > 80.0) return false;
        State3D y0 = make_rpo_state(q[0], q[1], q[2], q[3], m);
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[4], m, 1e-13, 600000)) return false;
        State3D yt = rotate_state(y, u, -target_a);
        F[0] = yt[1]; F[1] = yt[2]; F[2] = yt[4]; F[3] = yt[5]; F[4] = yt[9]; F[5] = yt[12];
        return true;
    };
    double F[6];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<6;++i) cost += F[i]*F[i];
    for (int it = 0; it < 35; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[6][5], Fp[6], Fm[6];
        for (int col = 0; col < 5; ++col) {
            double pp[5], pm[5]; for (int i=0;i<5;++i) { pp[i]=p[i]; pm[i]=p[i]; }
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5] = {0};
        if (!solve_lm_step<6, 5>(J, F, lambda, d)) { lambda *= 10; continue; }
        double sn = 0; for (int i=0;i<5;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.2) for (int i=0;i<5;++i) d[i] *= (0.2/sn);
        double np[5], nF[6]; for (int i=0;i<5;++i) np[i] = p[i]+d[i];
        if (!eval(np, nF)) { lambda *= 10; continue; }
        double nc = 0; for (int i=0;i<6;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for (int i=0;i<5;++i) p[i]=np[i]; for (int i=0;i<6;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13, lambda*0.2); }
        else { lambda *= 5; if (lambda > 1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

static double compute_kz_2d(const State3D& y0, double T, double two_alpha, const Masses& m) {
    const int zidx[6] = {2, 5, 8, 11, 14, 17};
    double tr = 0.0;
    const double eps = 1e-7;
    for (int k = 0; k < 6; ++k) {
        int col = zidx[k];
        State3D yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        integrate3d(yp, T, m, 1e-13, 400000);
        integrate3d(ym, T, m, 1e-13, 400000);
        tr += (yp[col] - ym[col]) / (2.0 * eps);
    }
    // 2 trivial (+1) from CM_z/P_z and 2 trivial (e^{+/-i 2*alpha}) from SO(3) tilt
    return tr - (2.0 + 2.0 * std::cos(two_alpha));
}

static void floquet_2d_rpo(const State3D& y0, double T, double two_alpha, const Masses& m,
                           double& k1, double& k2, bool& stable, double& margin) {
    const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
    const double uz[3] = {0, 0, 1};
    double M[12][12], eps = 1e-7;
    for (int c = 0; c < 12; ++c) {
        int col = xyidx[c];
        State3D yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        integrate3d(yp, T, m, 1e-13, 400000);
        integrate3d(ym, T, m, 1e-13, 400000);
        if (std::abs(two_alpha) > 1e-14) {
            yp = rotate_state(yp, uz, -two_alpha);
            ym = rotate_state(ym, uz, -two_alpha);
        }
        for (int r = 0; r < 12; ++r) M[r][c] = (yp[xyidx[r]] - ym[xyidx[r]]) / (2.0 * eps);
    }
    double tr1 = 0, tr2 = 0;
    for (int i = 0; i < 12; ++i) {
        tr1 += M[i][i];
        for (int j = 0; j < 12; ++j) tr2 += M[i][j] * M[j][i];
    }
    // In 2D rotating by 2a: 4 CM eigenvalues = 2 pairs of e^{+/- i 2a}, plus 4 (+1) from E/t and L/rot
    double c1 = std::cos(two_alpha), c2 = std::cos(2.0 * two_alpha);
    double s = tr1 - (4.0 + 4.0 * c1);
    double sum_k2 = tr2 - (4.0 + 4.0 * c2) + 4.0;
    double p = 0.5 * (s * s - sum_k2);
    double disc = s * s - 4.0 * p;
    stable = false; margin = -1.0; k1 = k2 = 0;
    if (disc >= -1e-4) {
        double sd = std::sqrt(std::max(0.0, disc));
        k1 = 0.5 * (s + sd); k2 = 0.5 * (s - sd);
        double mk = std::max(std::abs(k1), std::abs(k2));
        margin = 2.0 - mk;
        stable = (mk < 2.0 - 1e-4);
    }
}

static double min_eig_3x3(double C[3][3]) {
    double I1 = C[0][0]+C[1][1]+C[2][2];
    double I2 = C[0][0]*C[1][1]+C[0][0]*C[2][2]+C[1][1]*C[2][2]
              - C[0][1]*C[0][1]-C[0][2]*C[0][2]-C[1][2]*C[1][2];
    double I3 = C[0][0]*(C[1][1]*C[2][2]-C[1][2]*C[1][2])
              - C[0][1]*(C[0][1]*C[2][2]-C[1][2]*C[0][2])
              + C[0][2]*(C[0][1]*C[1][2]-C[1][1]*C[0][2]);
    double p = I2 - I1*I1/3.0, q = -2.0*I1*I1*I1/27.0 + I1*I2/3.0 - I3;
    double r = std::sqrt(std::max(0.0, -p/3.0));
    if (r < 1e-15) return std::max(0.0, I1/3.0);
    double phi = std::acos(std::clamp(-q/(2.0*r*r*r), -1.0, 1.0)) / 3.0;
    return std::max(0.0, std::min({I1/3.0+2*r*std::cos(phi),
                                   I1/3.0+2*r*std::cos(phi+2*M_PI/3),
                                   I1/3.0+2*r*std::cos(phi+4*M_PI/3)}));
}

struct Stab3D { bool stable; double k1, k2, k3, margin; int surv100; double drift100; };

static Stab3D analyze_stab3d_rpo(const State3D& y0, double T, const double u[3], double two_alpha,
                                 const Masses& m, bool run_100T = true) {
    double M[18][18], eps = 1e-7;
    Stab3D st{}; st.margin = -1.0;
    for (int col = 0; col < 18; ++col) {
        State3D yp = y0, ym = y0;
        yp[col] += eps; ym[col] -= eps;
        if (!integrate3d(yp, T, m, 1e-13, 800000) || !integrate3d(ym, T, m, 1e-13, 800000)) return st;
        yp = rotate_state(yp, u, -two_alpha);
        ym = rotate_state(ym, u, -two_alpha);
        for (int r = 0; r < 18; ++r) M[r][col] = (yp[r] - ym[r]) / (2.0 * eps);
    }
    double M2[18][18] = {{0}};
    for (int i=0;i<18;++i) for (int k=0;k<18;++k) for (int j=0;j<18;++j) M2[i][j] += M[i][k]*M[k][j];
    double tr1 = 0, tr2 = 0, tr3 = 0;
    for (int i=0;i<18;++i) {
        tr1 += M[i][i]; tr2 += M2[i][i];
        for (int j=0;j<18;++j) tr3 += M2[i][j]*M[j][i];
    }
    double c1 = std::cos(two_alpha), c2 = std::cos(2.0 * two_alpha), c3 = std::cos(3.0 * two_alpha);
    double p1 = tr1 - (6.0 + 6.0 * c1);
    double p2 = tr2 - (6.0 + 6.0 * c2) + 6.0;
    double p3 = tr3 - (6.0 + 6.0 * c3) + 3.0 * p1;
    double e1 = p1, e2 = 0.5*(e1*p1 - p2), e3 = (e2*p1 - e1*p2 + p3)/3.0;
    double Q = (e1*e1 - 3.0*e2)/9.0, R = (2.0*e1*e1*e1 - 9.0*e1*e2 + 27.0*e3)/54.0, disc = Q*Q*Q - R*R;
    if (disc >= -1e-4 && Q > 0) {
        double sqQ = std::sqrt(Q), th = std::acos(std::clamp(R/(sqQ*sqQ*sqQ), -1.0, 1.0));
        st.k1 = 2*sqQ*std::cos(th/3) + e1/3;
        st.k2 = 2*sqQ*std::cos((th + 2*M_PI)/3) + e1/3;
        st.k3 = 2*sqQ*std::cos((th + 4*M_PI)/3) + e1/3;
        double mk = std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
        st.margin = 2.0 - mk;
        st.stable = (mk < 2.0 - 1e-4);
    }
    if (st.stable && run_100T) {
        State3D yp = y0;
        double kick = 1e-5;
        yp[9]  += kick; yp[12] -= kick*m[0]/m[1];
        yp[10] += kick; yp[13] -= kick*m[0]/m[1];
        yp[11] += kick; yp[14] -= kick*m[0]/m[1];
        double md = 0; int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate3d(yp, T, m, 1e-11, 500000)) break;
            State3D yrot = rotate_state(yp, u, -two_alpha * (k + 1));
            double d = state_dist3d(yrot, y0);
            if (d > md) md = d;
            if (d > 5.0) break;
            surv = k + 1;
        }
        st.surv100 = surv; st.drift100 = md;
    }
    return st;
}

struct Geom3D { double rmin, rmax, thick; int swaps; };
static Geom3D measure_geom3d(const State3D& y0, double T, const Masses& m) {
    const int N = 300;
    double dt = T / N, C[3][3] = {{0}}, rmin = 1e300, rmax = 0;
    int prev_cl = -1, swaps = 0, cnt = 0;
    State3D yc = y0;
    for (int k = 0; k <= N; ++k) {
        double r12 = std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
        double r13 = std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
        double r23 = std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
        rmin = std::min({rmin, r12, r13, r23});
        rmax = std::max({rmax, r12, r13, r23});
        int cl = (r12<=r13 && r12<=r23) ? 0 : ((r13<=r23) ? 1 : 2);
        if (prev_cl != -1 && cl != prev_cl) swaps++;
        prev_cl = cl;
        for (int b = 0; b < 3; ++b) {
            double qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
            for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += qv[r]*qv[c];
            cnt++;
        }
        if (k < N) integrate3d(yc, dt, m, 1e-11, 80000);
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
    return {rmin, rmax, std::sqrt(min_eig_3x3(C)), swaps};
}

static void emit_3d_json(const std::string& subspace, const Masses& m, const State3D& y0,
                         double T, double two_alpha, double ret_err, const Stab3D& st,
                         const Geom3D& gm, int q_fold = 0, double T_rpo = 0.0) {
    double E = energy3d(y0, m);
    auto Lv = angmom3d(y0, m);
    double Lmag = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"dim\":3,\"subspace\":\"" << subspace << "\",\"masses\":["
       << m[0] << "," << m[1] << "," << m[2] << "],\"y0\":[";
    for (int i=0;i<18;++i) os << y0[i] << (i<17?",":"");
    os << "],\"T\":" << T;
    if (q_fold > 0) os << ",\"T_rpo\":" << T_rpo << ",\"q_fold\":" << q_fold;
    os << ",\"precess_rad\":" << two_alpha << ",\"E\":" << E << ",\"L\":" << Lmag
       << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]"
       << ",\"ret_err\":" << ret_err << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2 << ",\"k3\":" << st.k3
       << ",\"stab_margin\":" << st.margin << ",\"survived_100T\":" << st.surv100
       << ",\"drift_100T\":" << st.drift100 << ",\"rmin\":" << gm.rmin << ",\"rmax\":" << gm.rmax
       << ",\"thickness_3d\":" << gm.thick << ",\"partner_swaps\":" << gm.swaps
       << ",\"democratic\":" << (gm.swaps >= 2 ? "true" : "false") << "}\n";
    std::cout << os.str() << std::flush;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  branch_3d_vertical --trace-2d <m1> <m2> <m3> <x2> <vy1> <vy2> <Th> <dm3> <nsteps>\n"
                  << "  branch_3d_vertical --trace-2d-x2 <m1> <m2> <m3> <x2> <vy1> <vy2> <Th> <dx2> <nsteps>\n"
                  << "  branch_3d_vertical --branch-pd <m1> <m2> <m3> <x2> <vy1> <vy2> <Th_2d> <vz2_0> <dvz2> <nsteps> [--shift-half] [--json]\n"
                  << "  branch_3d_vertical --trace-3d <m1> <m2> <m3> <x2> <vy1> <vy2> <vz2> <Th> <alpha> <vz|m3> <step> <nsteps> [--json]\n"
                  << "  branch_3d_vertical --snap <m1> <m2> <m3> <x2> <vy1> <vy2> <vz2> <Th> <num_2a> <den_2a> <q_fold> <label>\n";
        return 1;
    }
    std::string mode = argv[1];
    std::cout << std::setprecision(12);

    if (mode == "--trace-2d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double p[4] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8])};
        double dm3 = std::atof(argv[9]);
        int nsteps = std::atoi(argv[10]);
        double prev_kz = 0, prev_m3 = m[2];
        double prev_p[4] = {p[0], p[1], p[2], p[3]};
        bool have_prev = false;

        for (int s = 0; s <= nsteps; ++s) {
            Masses cm = {m[0], m[1], m[2] + s * dm3};
            if (cm[2] <= 0.05) break;
            double ferr = 0;
            if (!refine_2d_perp(p, cm, ferr)) break;
            State3D y0 = make_rpo_state(p[0], p[1], p[2], 0.0, cm);
            double T = 2.0 * p[3];
            double k1, k2, marg; bool stab;
            floquet_2d_rpo(y0, T, 0.0, cm, k1, k2, stab, marg);
            double kz = compute_kz_2d(y0, T, 0.0, cm);
            std::cout << cm[2] << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
                      << T << " " << k1 << " " << k2 << " " << kz << " " << (stab ? 1 : 0) << " " << marg << "\n" << std::flush;

            // Detect kz crossing -2 or +2 and bisect!
            if (have_prev && ((prev_kz > -2.0 && kz <= -2.0) || (prev_kz < -2.0 && kz >= -2.0))) {
                double lo = prev_m3, hi = cm[2];
                double plo[4] = {prev_p[0], prev_p[1], prev_p[2], prev_p[3]};
                double pmid[4] = {prev_p[0], prev_p[1], prev_p[2], prev_p[3]};
                double kz_lo = prev_kz;
                for (int bi = 0; bi < 28; ++bi) {
                    double mid = 0.5 * (lo + hi);
                    Masses mm = {m[0], m[1], mid};
                    for (int i=0;i<4;++i) pmid[i] = plo[i];
                    double fe = 0;
                    if (!refine_2d_perp(pmid, mm, fe)) break;
                    double kzm = compute_kz_2d(make_rpo_state(pmid[0], pmid[1], pmid[2], 0.0, mm), 2.0*pmid[3], 0.0, mm);
                    if ((kz_lo + 2.0) * (kzm + 2.0) > 0) { lo = mid; kz_lo = kzm; for(int i=0;i<4;++i) plo[i]=pmid[i]; }
                    else hi = mid;
                }
                double mstar = 0.5 * (lo + hi);
                Masses ms = {m[0], m[1], mstar};
                double fe = 0; refine_2d_perp(pmid, ms, fe);
                double kzs = compute_kz_2d(make_rpo_state(pmid[0], pmid[1], pmid[2], 0.0, ms), 2.0*pmid[3], 0.0, ms);
                std::cerr << std::setprecision(13)
                          << "# BIFURCATION kz=-2 at m3*=" << mstar << " kz=" << kzs
                          << " x2=" << pmid[0] << " vy1=" << pmid[1] << " vy2=" << pmid[2]
                          << " Th=" << pmid[3] << " T=" << 2.0*pmid[3] << "\n";
            }
            prev_kz = kz; prev_m3 = cm[2];
            for (int i=0;i<4;++i) prev_p[i] = p[i];
            have_prev = true;
        }
    } else if (mode == "--trace-2d-x2") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double x2_0 = std::atof(argv[5]);
        double p[4] = {std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), 0.0};
        double dx2 = std::atof(argv[9]);
        int nsteps = std::atoi(argv[10]);
        for (int s = 0; s <= nsteps; ++s) {
            double cx2 = x2_0 + s * dx2;
            double ferr = 0;
            if (!refine_2d_rpo_fixed_x2(p, cx2, m, ferr)) break;
            State3D y0 = make_rpo_state(cx2, p[0], p[1], 0.0, m);
            double T = 2.0 * p[2], two_a = 2.0 * p[3];
            double k1, k2, marg; bool stab;
            floquet_2d_rpo(y0, T, two_a, m, k1, k2, stab, marg);
            double kz = compute_kz_2d(y0, T, two_a, m);
            std::cout << cx2 << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " "
                      << T << " " << two_a << " " << k1 << " " << k2 << " " << kz << " "
                      << (stab ? 1 : 0) << " " << marg << "\n" << std::flush;
        }
    } else if (mode == "--branch-pd") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double x2 = std::atof(argv[5]), vy1 = std::atof(argv[6]), vy2 = std::atof(argv[7]), Th2d = std::atof(argv[8]);
        double vz0 = std::atof(argv[9]), dvz = std::atof(argv[10]);
        int nsteps = std::atoi(argv[11]);
        bool shift_half = false, emit_j = false;
        for (int i = 12; i < argc; ++i) {
            if (std::string(argv[i]) == "--shift-half") shift_half = true;
            if (std::string(argv[i]) == "--json") emit_j = true;
        }
        if (shift_half) {
            State3D ys = make_rpo_state(x2, vy1, vy2, 0.0, m);
            integrate3d(ys, Th2d, m, 1e-13, 300000);
            // Rescale so x1 = -1 (flip sign if x1 > 0)
            double sgn = (ys[0] < 0) ? 1.0 : -1.0;
            double sc = std::abs(ys[0]);
            x2 = sgn * ys[3] / sc;
            vy1 = sgn * ys[10] * std::sqrt(sc);
            vy2 = sgn * ys[13] * std::sqrt(sc);
            Th2d = Th2d / std::pow(sc, 1.5);
            std::cerr << "# shifted half-period seed: x2=" << x2 << " vy1=" << vy1 << " vy2=" << vy2 << " Th2d=" << Th2d << "\n";
        }
        // Period-doubled 3D RPO has half-period Th_3d = 2 * Th2d (= T_2d)!!
        double p[5] = {x2, vy1, vy2, 2.0 * Th2d, 0.0};
        for (int s = 0; s <= nsteps; ++s) {
            double cvz = vz0 + s * dvz;
            if (cvz < 1e-4) break;
            double ferr = 0;
            if (!refine_3d_rpo_fixed_vz(p, cvz, m, ferr)) {
                std::cerr << "# branch failed at vz2=" << cvz << " ferr=" << ferr << "\n";
                break;
            }
            State3D y0 = make_rpo_state(p[0], p[1], p[2], cvz, m);
            double T = 2.0 * p[3], two_a = 2.0 * p[4];
            auto Lv = angmom3d(y0, m);
            double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
            double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
            State3D yT = y0;
            integrate3d(yT, T, m, 1e-13, 800000);
            double ret = state_dist3d(rotate_state(yT, u, -two_a), y0);
            auto gm = measure_geom3d(y0, T, m);
            auto st = analyze_stab3d_rpo(y0, T, u, two_a, m, emit_j);
            if (emit_j) {
                if (st.stable && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7) {
                    emit_3d_json("3d_rpo_vertical_pd", m, y0, T, two_a, ret, st, gm);
                }
                std::cerr << std::setprecision(10)
                          << "pd vz2=" << cvz << " T=" << T << " 2a=" << two_a
                          << " thick=" << gm.thick << " rmin=" << gm.rmin << " swaps=" << gm.swaps
                          << " k=(" << st.k1 << "," << st.k2 << "," << st.k3 << ") stab=" << st.stable
                          << " marg=" << st.margin << " surv=" << st.surv100 << " ret=" << ret << "\n";
            } else {
                std::cout << m[2] << " " << cvz << " " << T << " " << two_a << " " << gm.thick << " "
                          << gm.rmin << " " << gm.swaps << " " << st.k1 << " " << st.k2 << " " << st.k3 << " "
                          << (st.stable ? 1 : 0) << " " << st.margin << " "
                          << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << "\n" << std::flush;
            }
        }
    } else if (mode == "--trace-3d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double p[5] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[9]), std::atof(argv[10])};
        double vz0 = std::atof(argv[8]);
        std::string axis = argv[11];
        double step = std::atof(argv[12]);
        int nsteps = std::atoi(argv[13]);
        bool emit_j = false;
        for (int i = 14; i < argc; ++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int s = 0; s <= nsteps; ++s) {
            Masses cm = m;
            double cvz = vz0;
            if (axis == "m3") cm[2] = m[2] + s * step;
            else if (axis == "vz") cvz = vz0 + s * step;
            if (cm[2] < 0.05 || cvz < 1e-4) break;
            double ferr = 0;
            if (!refine_3d_rpo_fixed_vz(p, cvz, cm, ferr)) break;
            State3D y0 = make_rpo_state(p[0], p[1], p[2], cvz, cm);
            double T = 2.0 * p[3], two_a = 2.0 * p[4];
            auto Lv = angmom3d(y0, cm);
            double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
            double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
            State3D yT = y0;
            integrate3d(yT, T, cm, 1e-13, 800000);
            double ret = state_dist3d(rotate_state(yT, u, -two_a), y0);
            auto gm = measure_geom3d(y0, T, cm);
            auto st = analyze_stab3d_rpo(y0, T, u, two_a, cm, emit_j);
            if (emit_j) {
                if (st.stable && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7) {
                    emit_3d_json("3d_rpo_vertical_pd", cm, y0, T, two_a, ret, st, gm);
                }
                std::cerr << std::setprecision(10)
                          << "tr m3=" << cm[2] << " vz2=" << cvz << " T=" << T << " 2a=" << two_a
                          << " thick=" << gm.thick << " rmin=" << gm.rmin << " swaps=" << gm.swaps
                          << " k=(" << st.k1 << "," << st.k2 << "," << st.k3 << ") stab=" << st.stable
                          << " marg=" << st.margin << " surv=" << st.surv100 << " ret=" << ret << "\n";
            } else {
                std::cout << cm[2] << " " << cvz << " " << T << " " << two_a << " " << gm.thick << " "
                          << gm.rmin << " " << gm.swaps << " " << st.k1 << " " << st.k2 << " " << st.k3 << " "
                          << (st.stable ? 1 : 0) << " " << st.margin << " "
                          << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << "\n" << std::flush;
            }
        }
    } else if (mode == "--snap") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double p[5] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9])};
        double num = std::atof(argv[10]), den = std::atof(argv[11]);
        int q_fold = std::atoi(argv[12]);
        std::string label = (argc > 13) ? argv[13] : ("3d_inertial_" + std::to_string(q_fold) + "fold");
        // target 2*alpha = 2*pi * num / den  =>  alpha = pi * num / den
        double target_a = M_PI * num / den;
        double ferr = 0;
        if (!refine_3d_rpo_fixed_alpha(p, target_a, m, ferr)) {
            std::cerr << "Snap failed, ferr=" << ferr << "\n";
            return 2;
        }
        State3D y0 = make_rpo_state(p[0], p[1], p[2], p[3], m);
        double T_rpo = 2.0 * p[4];
        double T_inert = q_fold * T_rpo;
        double two_a = 2.0 * target_a;
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};

        State3D yq = y0;
        integrate3d(yq, T_inert, m, 1e-13, 2000000);
        double ret_inert = state_dist3d(yq, y0);
        auto gm = measure_geom3d(y0, T_rpo, m);
        gm.swaps *= q_fold;
        auto st = analyze_stab3d_rpo(y0, T_rpo, u, two_a, m, true);
        std::cerr << std::setprecision(13)
                  << "# SNAPPED " << label << " q=" << q_fold << " 2a=" << num << "*2pi/" << den
                  << " half_err=" << ferr << " T_rpo=" << T_rpo << " T_inert=" << T_inert
                  << " ret_inert=" << ret_inert << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                  << ") marg=" << st.margin << " surv=" << st.surv100 << " thick=" << gm.thick << "\n";
        if (st.stable && st.surv100 >= 100 && ret_inert < 1e-7 && gm.thick >= 0.05) {
            emit_3d_json(label, m, y0, T_inert, two_a, ret_inert, st, gm, q_fold, T_rpo);
        }
    }
    return 0;
}
