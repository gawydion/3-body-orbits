#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static constexpr Masses EQ_MASSES = {1.0, 1.0, 1.0};

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

static inline double energy3d(const State3D& s, const Masses& m) {
    double K = 0, U = 0;
    for (int i = 0; i < 3; ++i)
        K += 0.5*m[i]*(s[9+3*i]*s[9+3*i]+s[10+3*i]*s[10+3*i]+s[11+3*i]*s[11+3*i]);
    for (int i = 0; i < 3; ++i)
        for (int j = i+1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            U -= m[i]*m[j] / std::sqrt(dx*dx+dy*dy+dz*dz);
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
    double s = 0; for (int i = 0; i < 18; ++i) { double d = a[i]-b[i]; s += d*d; }
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
    State3D k1,k2,k3,k4,k5,k6,k7,yt,yn;
    double rm = 1e300;
    for (int att=0;att<25;++att) {
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
            y = yn; t += h; if (rm < rmin_track) rmin_track = rm;
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13,
                        int max_steps = 600000, double* rmin_out = nullptr, double r_abort = 1e-4) {
    double t = 0, h = std::min(1e-3, T*0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step3d(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 80.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

// Gauge-fixed symmetric choreography start (invariant under tau_1 = (t->-t) o (1 2) o R_y(pi)):
// r1=(-1, y0, 0), r2=(1, y0, 0), r3=(0, -2y0, 0)
// v1=(vx, vy, vz), v2=(vx, -vy, vz), v3=(-2vx, 0, -2vz)
// Exact L = (6*y0*vz, 0, -2*vy - 6*y0*vx) perpendicular to y-axis!
static State3D make_choreo_state(double y0, double vx, double vy, double vz) {
    return {-1.0, y0, 0.0,   1.0, y0, 0.0,   0.0, -2.0*y0, 0.0,
            vx, vy, vz,      vx, -vy, vz,    -2.0*vx, 0.0, -2.0*vz};
}

// Check choreography condition at t = T/3: r1(T/3)=R(2a)r2(0), r2(T/3)=R(2a)r3(0), r3(T/3)=R(2a)r1(0)
static double choreo_third_err(const State3D& y0, double T, const double u[3], double phi_third) {
    State3D y = y0;
    if (!integrate3d(y, T / 3.0, EQ_MASSES, 1e-13, 400000)) return 1e9;
    State3D yr = rotate_state(y, u, -phi_third);
    // Expected cyclic shift: body 0->1, 1->2, 2->0
    double e1 = 0, e2 = 0;
    for (int b = 0; b < 3; ++b) {
        int nb1 = (b + 1) % 3;
        int nb2 = (b + 2) % 3;
        for (int c = 0; c < 3; ++c) {
            double dr1 = yr[3*b + c] - y0[3*nb1 + c];
            double dv1 = yr[9 + 3*b + c] - y0[9 + 3*nb1 + c];
            e1 += dr1*dr1 + dv1*dv1;
            double dr2 = yr[3*b + c] - y0[3*nb2 + c];
            double dv2 = yr[9 + 3*b + c] - y0[9 + 3*nb2 + c];
            e2 += dr2*dr2 + dv2*dv2;
        }
    }
    return std::sqrt(std::min(e1, e2));
}

template <int NEQ, int NVAR>
static bool solve_lm(const double J[NEQ][NVAR], const double F[NEQ], double lambda, double d[NVAR]) {
    double A[NVAR][NVAR + 1] = {{0}};
    for (int r = 0; r < NVAR; ++r) {
        for (int c = 0; c < NVAR; ++c) {
            double s = 0; for (int i=0;i<NEQ;++i) s += J[i][r]*J[i][c];
            A[r][c] = s;
        }
        A[r][r] += lambda*(1.0 + A[r][r]);
        double s = 0; for (int i=0;i<NEQ;++i) s -= J[i][r]*F[i];
        A[r][NVAR] = s;
    }
    for (int col = 0; col < NVAR; ++col) {
        int piv = col;
        for (int r = col+1; r < NVAR; ++r) if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
        if (std::abs(A[piv][col]) < 1e-20) return false;
        for (int c = col; c <= NVAR; ++c) std::swap(A[col][c], A[piv][c]);
        for (int r = col+1; r < NVAR; ++r) {
            double f = A[r][col] / A[col][col];
            for (int c = col; c <= NVAR; ++c) A[r][c] -= f*A[col][c];
        }
    }
    for (int r = NVAR-1; r >= 0; --r) {
        double s = A[r][NVAR];
        for (int c = r+1; c < NVAR; ++c) s -= A[r][c]*d[c];
        d[r] = s / A[r][r];
    }
    return true;
}

// 2D rotating choreography refine at fixed alpha: p = {y0, vx, vy, T6}
static bool refine_2d_fixed_alpha(double p[4], double alpha, double& ferr) {
    double lambda = 1e-4;
    const double uz[3] = {0, 0, 1};
    auto eval = [&](const double q[4], double F[4]) -> bool {
        if (q[3] < 0.25 || q[3] > 20.0) return false;
        State3D y = make_choreo_state(q[0], q[1], q[2], 0.0);
        double d13 = std::hypot(-1.0, 3.0*q[0]);
        if (d13 < 0.15) return false;
        if (!integrate3d(y, q[3], EQ_MASSES, 1e-12, 120000, nullptr, 5e-3)) return false;
        State3D yt = rotate_state(y, uz, -alpha);
        F[0] = yt[0]; F[1] = yt[4] - yt[7]; F[2] = yt[10]; F[3] = yt[12] - yt[15];
        return true;
    };
    double F[4];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<4;++i) cost += F[i]*F[i];
    for (int it = 0; it < 32; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[4][4], Fp[4], Fm[4];
        for (int col=0;col<4;++col) {
            double pp[4], pm[4]; for(int i=0;i<4;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<4;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[4]={0};
        if (!solve_lm<4,4>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn = 0; for(int i=0;i<4;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<4;++i) d[i]*=(0.25/sn);
        double np[4], nF[4]; for(int i=0;i<4;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc = 0; for(int i=0;i<4;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<4;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// 2D rotating choreography refine at fixed y0: p = {vx, vy, T6, alpha}
static bool refine_2d_fixed_y0(double p[4], double y0, double& ferr) {
    double lambda = 1e-4;
    const double uz[3] = {0, 0, 1};
    auto eval = [&](const double q[4], double F[4]) -> bool {
        if (q[2] < 0.25 || q[2] > 20.0) return false;
        State3D y = make_choreo_state(y0, q[0], q[1], 0.0);
        if (!integrate3d(y, q[2], EQ_MASSES, 1e-12, 120000, nullptr, 5e-3)) return false;
        State3D yt = rotate_state(y, uz, -q[3]);
        F[0] = yt[0]; F[1] = yt[4] - yt[7]; F[2] = yt[10]; F[3] = yt[12] - yt[15];
        return true;
    };
    double F[4];
    if (!eval(p, F)) return false;
    double cost = 0; for (int i=0;i<4;++i) cost += F[i]*F[i];
    for (int it = 0; it < 30; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        if (it == 8 && std::sqrt(cost) > 0.25) return false;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[4][4], Fp[4], Fm[4];
        for (int col=0;col<4;++col) {
            double pp[4], pm[4]; for(int i=0;i<4;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<4;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[4]={0};
        if (!solve_lm<4,4>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn = 0; for(int i=0;i<4;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<4;++i) d[i]*=(0.25/sn);
        double np[4], nF[4]; for(int i=0;i<4;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc = 0; for(int i=0;i<4;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<4;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// 3D helical rotating choreography refine at fixed vz: p = {y0, vx, vy, T6, alpha}
static bool refine_3d_fixed_vz(double p[5], double vz, double& ferr) {
    double lambda = 1e-4;
    auto eval = [&](const double q[5], double F[6]) -> bool {
        if (q[3] < 0.25 || q[3] > 20.0) return false;
        State3D y0 = make_choreo_state(q[0], q[1], q[2], vz);
        if (std::hypot(-1.0, 3.0*q[0]) < 0.15) return false;
        auto Lv = angmom3d(y0, EQ_MASSES);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-3) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[3], EQ_MASSES, 1e-12, 150000, nullptr, 5e-3)) return false;
        State3D yt = rotate_state(y, u, -q[4]);
        F[0] = yt[0]; F[1] = yt[2]; F[2] = yt[4] - yt[7];
        F[3] = yt[10]; F[4] = yt[12] - yt[15]; F[5] = yt[14] - yt[17];
        return true;
    };
    double F[6];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<6;++i) cost += F[i]*F[i];
    for (int it = 0; it < 34; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        if (it == 9 && std::sqrt(cost) > 0.25) return false;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[6][5], Fp[6], Fm[6];
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5]={0};
        if (!solve_lm<6,5>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<5;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<5;++i) d[i]*=(0.25/sn);
        double np[5], nF[6]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<6;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<5;++i) p[i]=np[i]; for(int i=0;i<6;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// 3D helical rotating choreography snap at fixed alpha: p = {y0, vx, vy, vz, T6}
static bool refine_3d_fixed_alpha(double p[5], double target_a, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[5], double F[6]) -> bool {
        if (q[4] < 0.25 || q[4] > 20.0) return false;
        State3D y0 = make_choreo_state(q[0], q[1], q[2], q[3]);
        auto Lv = angmom3d(y0, EQ_MASSES);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-3) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[4], EQ_MASSES, 1e-13, 300000, nullptr, 5e-3)) return false;
        State3D yt = rotate_state(y, u, -target_a);
        F[0] = yt[0]; F[1] = yt[2]; F[2] = yt[4] - yt[7];
        F[3] = yt[10]; F[4] = yt[12] - yt[15]; F[5] = yt[14] - yt[17];
        return true;
    };
    double F[6];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<6;++i) cost += F[i]*F[i];
    for (int it = 0; it < 35; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[6][5], Fp[6], Fm[6];
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5]={0};
        if (!solve_lm<6,5>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<5;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.2) for(int i=0;i<5;++i) d[i]*=(0.2/sn);
        double np[5], nF[6]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<6;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<5;++i) p[i]=np[i]; for(int i=0;i<6;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
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

struct StabFull {
    bool stable_2d, stable_3d;
    double k1, k2, k3; // for 2D, k3 is vertical kz!
    double margin_2d, margin_3d;
    int surv100;
    double drift100;
};

static StabFull analyze_stab_choreo(const State3D& y0, double T, const double u[3], double phi_full,
                                    int dim, bool run_100T = true) {
    StabFull st{}; st.margin_2d = -1; st.margin_3d = -1;
    const double eps = 1e-7;
    if (dim == 2) {
        const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        double M[12][12];
        for (int c = 0; c < 12; ++c) {
            int col = xyidx[c];
            State3D yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            if (!integrate3d(yp, T, EQ_MASSES, 1e-13, 500000) ||
                !integrate3d(ym, T, EQ_MASSES, 1e-13, 500000)) return st;
            yp = rotate_state(yp, u, -phi_full);
            ym = rotate_state(ym, u, -phi_full);
            for (int r = 0; r < 12; ++r) M[r][c] = (yp[xyidx[r]] - ym[xyidx[r]]) / (2.0 * eps);
        }
        double tr1 = 0, tr2 = 0;
        for (int i=0;i<12;++i) { tr1 += M[i][i]; for (int j=0;j<12;++j) tr2 += M[i][j]*M[j][i]; }
        double c1 = std::cos(phi_full), c2 = std::cos(2.0 * phi_full);
        double s = tr1 - (4.0 + 4.0 * c1);
        double sum_k2 = tr2 - (4.0 + 4.0 * c2) + 4.0;
        double p = 0.5 * (s * s - sum_k2);
        double disc = s * s - 4.0 * p;
        if (disc >= -1e-4) {
            double sd = std::sqrt(std::max(0.0, disc));
            st.k1 = 0.5*(s + sd); st.k2 = 0.5*(s - sd);
            double mk = std::max(std::abs(st.k1), std::abs(st.k2));
            st.margin_2d = 2.0 - mk;
            st.stable_2d = (mk < 2.0 - 1e-4);
        }
        // Vertical kz
        const int zidx[6] = {2, 5, 8, 11, 14, 17};
        double trz = 0;
        for (int k=0;k<6;++k) {
            int col = zidx[k];
            State3D yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            integrate3d(yp, T, EQ_MASSES, 1e-13, 500000);
            integrate3d(ym, T, EQ_MASSES, 1e-13, 500000);
            trz += (yp[col] - ym[col]) / (2.0 * eps);
        }
        st.k3 = trz - (2.0 + 2.0 * std::cos(phi_full));
        st.stable_3d = st.stable_2d && (std::abs(st.k3) < 2.0 - 1e-4);
        st.margin_3d = st.stable_3d ? (2.0 - std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)})) : -1.0;

        if (st.stable_2d && run_100T) {
            State3D yp = y0;
            double kick = 1e-5;
            yp[9] += kick; yp[12] -= kick;
            yp[10] += kick; yp[13] -= kick;
            if (st.stable_3d) { yp[11] += kick; yp[14] -= kick; }
            double md = 0; int surv = 0;
            for (int k = 0; k < 100; ++k) {
                if (!integrate3d(yp, T, EQ_MASSES, 1e-11, 400000)) break;
                State3D yr = rotate_state(yp, u, -phi_full * (k + 1));
                double d = state_dist3d(yr, y0);
                if (d > md) md = d; if (d > 5.0) break;
                surv = k + 1;
            }
            st.surv100 = surv; st.drift100 = md;
        }
    } else {
        double M[18][18];
        for (int col=0;col<18;++col) {
            State3D yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            if (!integrate3d(yp, T, EQ_MASSES, 1e-13, 600000) ||
                !integrate3d(ym, T, EQ_MASSES, 1e-13, 600000)) return st;
            yp = rotate_state(yp, u, -phi_full);
            ym = rotate_state(ym, u, -phi_full);
            for (int r=0;r<18;++r) M[r][col] = (yp[r]-ym[r])/(2.0*eps);
        }
        double M2[18][18] = {{0}};
        for (int i=0;i<18;++i) for (int k=0;k<18;++k) for (int j=0;j<18;++j) M2[i][j] += M[i][k]*M[k][j];
        double tr1=0, tr2=0, tr3=0;
        for (int i=0;i<18;++i) { tr1+=M[i][i]; tr2+=M2[i][i]; for(int j=0;j<18;++j) tr3+=M2[i][j]*M[j][i]; }
        double c1=std::cos(phi_full), c2=std::cos(2*phi_full), c3=std::cos(3*phi_full);
        double p1 = tr1 - (6.0 + 6.0*c1);
        double p2 = tr2 - (6.0 + 6.0*c2) + 6.0;
        double p3 = tr3 - (6.0 + 6.0*c3) + 3.0*p1;
        double e1=p1, e2=0.5*(e1*p1-p2), e3=(e2*p1-e1*p2+p3)/3.0;
        double Q=(e1*e1-3*e2)/9.0, R=(2*e1*e1*e1-9*e1*e2+27*e3)/54.0, disc=Q*Q*Q-R*R;
        if (disc >= -1e-4 && Q > 0) {
            double sqQ=std::sqrt(Q), th=std::acos(std::clamp(R/(sqQ*sqQ*sqQ), -1.0, 1.0));
            st.k1 = 2*sqQ*std::cos(th/3) + e1/3;
            st.k2 = 2*sqQ*std::cos((th+2*M_PI)/3) + e1/3;
            st.k3 = 2*sqQ*std::cos((th+4*M_PI)/3) + e1/3;
            double mk = std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
            st.margin_3d = 2.0 - mk; st.margin_2d = st.margin_3d;
            st.stable_3d = (mk < 2.0 - 1e-4); st.stable_2d = st.stable_3d;
        }
        if (st.stable_3d && run_100T) {
            State3D yp = y0;
            double kick = 1e-5;
            yp[9]+=kick; yp[12]-=kick; yp[10]+=kick; yp[13]-=kick; yp[11]+=kick; yp[14]-=kick;
            double md = 0; int surv = 0;
            for (int k=0;k<100;++k) {
                if (!integrate3d(yp, T, EQ_MASSES, 1e-11, 400000)) break;
                State3D yr = rotate_state(yp, u, -phi_full * (k + 1));
                double d = state_dist3d(yr, y0);
                if (d > md) md = d; if (d > 5.0) break;
                surv = k + 1;
            }
            st.surv100 = surv; st.drift100 = md;
        }
    }
    return st;
}

struct Geom { double rmin, rmax, thick; int swaps; };
static Geom measure_geom(const State3D& y0, double T) {
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
        int cl = (r12<=r13 && r12<=r23)?0:((r13<=r23)?1:2);
        if (prev_cl != -1 && cl != prev_cl) swaps++;
        prev_cl = cl;
        for (int b=0;b<3;++b) {
            double qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
            for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += qv[r]*qv[c];
            cnt++;
        }
        if (k < N) integrate3d(yc, dt, EQ_MASSES, 1e-11, 80000);
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
    return {rmin, rmax, std::sqrt(min_eig_3x3(C)), swaps};
}

static void emit_choreo_json(int dim, const std::string& subspace, const State3D& y0, double T,
                             double phi_full, double choreo_err, double ret_err,
                             const StabFull& st, const Geom& gm, int q_fold = 0, double T_rpo = 0.0) {
    double E = energy3d(y0, EQ_MASSES);
    auto Lv = angmom3d(y0, EQ_MASSES);
    double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"dim\":" << dim << ",\"subspace\":\"" << subspace << "\",\"masses\":[1,1,1],\"y0\":[";
    if (dim == 2) {
        int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        for (int i=0;i<12;++i) os << y0[xyidx[i]] << (i<11?",":"");
    } else {
        for (int i=0;i<18;++i) os << y0[i] << (i<17?",":"");
    }
    os << "],\"T\":" << T;
    if (q_fold > 0) os << ",\"T_rpo\":" << T_rpo << ",\"q_fold\":" << q_fold;
    os << ",\"precess_rad\":" << phi_full
       << ",\"omega\":" << (phi_full / (q_fold > 0 ? T_rpo : T))
       << ",\"choreo_err\":" << choreo_err
       << ",\"E\":" << E << ",\"L\":" << (dim == 2 ? Lv[2] : Lmag);
    if (dim == 3) os << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]";
    os << ",\"ret_err\":" << ret_err
       << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2;
    if (dim == 3) {
        os << ",\"k3\":" << st.k3 << ",\"stab_margin\":" << st.margin_3d
           << ",\"thickness_3d\":" << gm.thick;
    } else {
        os << ",\"kz\":" << st.k3 << ",\"stable_3d\":" << (st.stable_3d ? "true" : "false")
           << ",\"stab_margin\":" << st.margin_2d;
    }
    os << ",\"survived_100T\":" << st.surv100 << ",\"drift_100T\":" << st.drift100
       << ",\"rmin\":" << gm.rmin << ",\"rmax\":" << gm.rmax
       << ",\"partner_swaps\":" << gm.swaps << ",\"democratic\":true}\n";
    std::cout << os.str() << std::flush;
}

static void run_scan_2d(int nsamples, double T6max, int seed_off) {
    struct Seed { double y0, vx, vy, T6, alpha, err; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;
    const double uz[3] = {0, 0, 1};

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(551001 + seed_off*10007 + tid*9973);
            std::uniform_real_distribution<double> u_y0(-1.2, 1.2);
            std::uniform_real_distribution<double> u_v(-1.3, 1.3);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double y0 = u_y0(rng), vx = u_v(rng), vy = u_v(rng);
                if (std::hypot(-1.0, 3.0*y0) < 0.3) continue;
                State3D s0 = make_choreo_state(y0, vx, vy, 0.0);
                if (energy3d(s0, EQ_MASSES) >= -0.08) continue;
                double Lz = angmom3d(s0, EQ_MASSES)[2];
                if (std::abs(Lz) < 0.05) continue;

                State3D s = s0;
                double t = 0, h = 1e-3, rm = 1e300;
                double best_e = 1e9, best_t = 0, best_a = 0, pe = 1e9, p2e = 1e9, pt = 0, pa = 0;
                int steps = 0;
                while (t < T6max && steps < 10000) {
                    double hc = std::min(h, 0.025);
                    if (!dopri5_step3d(s, t, hc, T6max, EQ_MASSES, 1e-8, rm, 1e-2)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i=0;i<9;++i) maxr = std::max(maxr, std::abs(s[i]));
                    if (maxr > 15.0) break;
                    double dx23 = s[3] - s[6], dy23 = s[4] - s[7];
                    double a0 = std::atan2(dy23, dx23);
                    double e = 1e9, ba = a0;
                    for (int sg = 0; sg < 2; ++sg) {
                        double a = a0 + sg * M_PI;
                        State3D st = rotate_state(s, uz, -a);
                        double ce = std::sqrt(st[0]*st[0] + st[10]*st[10] + (st[12]-st[15])*(st[12]-st[15]));
                        if (ce < e) { e = ce; ba = a; }
                    }
                    if (t > 0.35 && steps >= 3 && pe < p2e && pe < e && pe < best_e) {
                        best_e = pe; best_t = pt; best_a = pa;
                    }
                    p2e = pe; pe = e; pt = t; pa = ba;
                }
                if (best_e < 0.25 && best_t > 0.35) {
                    tseeds[tid].push_back({y0, vx, vy, best_t, best_a, best_e});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all;
    for (auto& v : tseeds) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), [](const Seed& a, const Seed& b){ return a.err < b.err; });
    size_t limit = std::min(all.size(), (size_t)600);
    std::cerr << "[choreo_2d] seeds=" << all.size() << " refining=" << limit << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 2>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                double y0 = all[idx].y0;
                double p[4] = {all[idx].vx, all[idx].vy, all[idx].T6, all[idx].alpha};
                double ferr = 0;
                if (!refine_2d_fixed_y0(p, y0, ferr)) continue;
                State3D s0 = make_choreo_state(y0, p[0], p[1], 0.0);
                double T = 6.0 * p[2], phi_full = 6.0 * p[3];
                double Lz = angmom3d(s0, EQ_MASSES)[2];
                if (std::abs(Lz) < 1e-3) continue;
                // Exclude trivial Lagrange equilateral triangle (constant r12=r13=r23)
                double cerr = choreo_third_err(s0, T, uz, 2.0 * p[3]);
                if (cerr > 1e-7) continue;
                State3D sT = s0;
                if (!integrate3d(sT, T, EQ_MASSES, 1e-13, 400000)) continue;
                double ret = state_dist3d(rotate_state(sT, uz, -phi_full), s0);
                if (ret > 1e-7) continue;
                auto gm = measure_geom(s0, T);
                if (gm.rmax - gm.rmin < 0.05) continue; // skip circular Lagrange equilateral!

                double E = energy3d(s0, EQ_MASSES);
                double inv1 = E * std::pow(T, 2.0/3.0);
                double inv2 = std::abs(Lz) * std::sqrt(std::max(0.0, -E));
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup)
                        if (std::abs(f[0]-inv1) < 2e-2 && std::abs(f[1]-inv2) < 2e-2) { dup = true; break; }
                    if (dup) continue;
                    dedup.push_back({inv1, inv2});
                }
                auto st = analyze_stab_choreo(s0, T, uz, phi_full, 2, true);
                if (!st.stable_2d || st.surv100 < 100) continue;
                std::lock_guard<std::mutex> lk(mtx);
                std::cerr << std::setprecision(10)
                          << "  [STABLE 2D CHOREO] y0=" << y0 << " vx=" << p[0] << " vy=" << p[1]
                          << " T6=" << p[2] << " a=" << p[3] << " T=" << T << " wT=" << phi_full
                          << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                          << " marg=" << st.margin_2d << " swaps=" << gm.swaps << "\n";
                emit_choreo_json(2, "2d_choreo_rotating", s0, T, phi_full, cerr, ret, st, gm);
            }
        });
    }
    for (auto& th : ths) th.join();
}

static void run_scan_3d(int nsamples, double T6max, int seed_off) {
    struct Seed { double y0, vx, vy, vz, T6, alpha, err; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(661001 + seed_off*10007 + tid*9973);
            std::uniform_real_distribution<double> u_y0(-1.1, 1.1);
            std::uniform_real_distribution<double> u_v(-1.1, 1.1);
            std::uniform_real_distribution<double> u_vz(0.15, 0.95);
            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double y0 = u_y0(rng), vx = u_v(rng), vy = u_v(rng), vz = u_vz(rng);
                if (std::abs(y0) < 0.08 || std::hypot(-1.0, 3.0*y0) < 0.3) continue;
                State3D s0 = make_choreo_state(y0, vx, vy, vz);
                if (energy3d(s0, EQ_MASSES) >= -0.08) continue;
                auto Lv = angmom3d(s0, EQ_MASSES);
                double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
                if (Lm < 0.1) continue;
                double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};

                State3D s = s0;
                double t = 0, h = 1e-3, rm = 1e300;
                double best_e = 1e9, best_t = 0, best_a = 0, pe = 1e9, p2e = 1e9, pt = 0, pa = 0;
                int steps = 0;
                while (t < T6max && steps < 10000) {
                    double hc = std::min(h, 0.025);
                    if (!dopri5_step3d(s, t, hc, T6max, EQ_MASSES, 1e-8, rm, 1e-2)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i=0;i<9;++i) maxr = std::max(maxr, std::abs(s[i]));
                    if (maxr > 15.0) break;
                    double dx = s[3]-s[6], dy = s[4]-s[7], dz = s[5]-s[8];
                    double a0 = std::atan2(dy, u[2]*dx - u[0]*dz);
                    double e = 1e9, ba = a0;
                    for (int sg = 0; sg < 2; ++sg) {
                        double a = a0 + sg * M_PI;
                        State3D st = rotate_state(s, u, -a);
                        double ce = std::sqrt(st[0]*st[0] + st[2]*st[2] + st[10]*st[10] +
                                              (st[12]-st[15])*(st[12]-st[15]) + (st[14]-st[17])*(st[14]-st[17]));
                        if (ce < e) { e = ce; ba = a; }
                    }
                    if (t > 0.35 && steps >= 3 && pe < p2e && pe < e && pe < best_e) {
                        best_e = pe; best_t = pt; best_a = pa;
                    }
                    p2e = pe; pe = e; pt = t; pa = ba;
                }
                if (best_e < 0.32 && best_t > 0.35) {
                    tseeds[tid].push_back({y0, vx, vy, vz, best_t, best_a, best_e});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all;
    for (auto& v : tseeds) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), [](const Seed& a, const Seed& b){ return a.err < b.err; });
    size_t limit = std::min(all.size(), (size_t)500);
    std::cerr << "[choreo_3d] seeds=" << all.size() << " refining=" << limit << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 2>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                double vz = all[idx].vz;
                double p[5] = {all[idx].y0, all[idx].vx, all[idx].vy, all[idx].T6, all[idx].alpha};
                double ferr = 0;
                if (!refine_3d_fixed_vz(p, vz, ferr)) continue;
                State3D s0 = make_choreo_state(p[0], p[1], p[2], vz);
                double T = 6.0 * p[3], phi_full = 6.0 * p[4];
                auto Lv = angmom3d(s0, EQ_MASSES);
                double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
                if (Lm < 1e-3) continue;
                double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};

                double cerr = choreo_third_err(s0, T, u, 2.0 * p[4]);
                if (cerr > 1e-7) continue;
                State3D sT = s0;
                if (!integrate3d(sT, T, EQ_MASSES, 1e-13, 500000)) continue;
                double ret = state_dist3d(rotate_state(sT, u, -phi_full), s0);
                if (ret > 1e-7) continue;
                auto gm = measure_geom(s0, T);
                if (gm.thick < 0.05) continue;

                double E = energy3d(s0, EQ_MASSES);
                double inv1 = E * std::pow(T, 2.0/3.0);
                double inv2 = Lm * std::sqrt(std::max(0.0, -E));
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup)
                        if (std::abs(f[0]-inv1) < 2e-2 && std::abs(f[1]-inv2) < 2e-2) { dup = true; break; }
                    if (dup) continue;
                    dedup.push_back({inv1, inv2});
                }
                auto st = analyze_stab_choreo(s0, T, u, phi_full, 3, true);
                if (!st.stable_3d || st.surv100 < 100) continue;
                std::lock_guard<std::mutex> lk(mtx);
                std::cerr << std::setprecision(10)
                          << "  [STABLE 3D HELICAL CHOREO] y0=" << p[0] << " vx=" << p[1] << " vy=" << p[2]
                          << " vz=" << vz << " T6=" << p[3] << " a=" << p[4] << " T=" << T << " wT=" << phi_full
                          << " thick=" << gm.thick << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                          << ") marg=" << st.margin_3d << "\n";
                emit_choreo_json(3, "3d_choreo_helical", s0, T, phi_full, cerr, ret, st, gm);
            }
        });
    }
    for (auto& th : ths) th.join();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  search_choreo_rot --scan-2d <nsamples> <T6max> [seed]\n"
                  << "  search_choreo_rot --scan-3d <nsamples> <T6max> [seed]\n"
                  << "  search_choreo_rot --continue-2d <y0> <vx> <vy> <T6> <a0> <da> <nsteps> [--json]\n"
                  << "  search_choreo_rot --branch-3d <y0> <vx> <vy> <T6> <a0> <vz0> <dvz> <nsteps> [--json]\n"
                  << "  search_choreo_rot --snap-2d <y0> <vx> <vy> <T6> <p_num> <q_den> <label>\n"
                  << "  search_choreo_rot --snap-3d <y0> <vx> <vy> <vz> <T6> <p_num> <q_den> <label>\n"
                  << "  search_choreo_rot --fig8-rot <da> <nsteps> [--json]\n";
        return 1;
    }
    std::string mode = argv[1];
    std::cout << std::setprecision(12);
    const double uz[3] = {0, 0, 1};

    if (mode == "--scan-2d") {
        int n = (argc > 2) ? std::atoi(argv[2]) : 30000;
        double t6 = (argc > 3) ? std::atof(argv[3]) : 5.0;
        int sd = (argc > 4) ? std::atoi(argv[4]) : 0;
        run_scan_2d(n, t6, sd);
    } else if (mode == "--scan-3d") {
        int n = (argc > 2) ? std::atoi(argv[2]) : 35000;
        double t6 = (argc > 3) ? std::atof(argv[3]) : 5.0;
        int sd = (argc > 4) ? std::atoi(argv[4]) : 0;
        run_scan_3d(n, t6, sd);
    } else if (mode == "--fig8-rot") {
        // Compute exact isosceles configuration of Figure-8 at t = T/12!
        // At t = 0: r1=(-1,0), r2=(1,0), r3=(0,0), v1=v2=(vx8,vy8), v3=(-2vx8,-2vy8).
        // At t = T/12: find the apex body (equidistant from the other two), permute it to body 3!
        double vx8 = 0.34711688812, vy8 = 0.532724945387, T8 = 6.32591398294;
        State3D f8 = {-1,0,0, 1,0,0, 0,0,0, vx8,vy8,0, vx8,vy8,0, -2*vx8,-2*vy8,0};
        integrate3d(f8, T8 / 12.0, EQ_MASSES, 1e-13, 200000);
        // Distances between pairs: find which body c has d(a,c) == d(b,c)
        double r01 = std::hypot(f8[3]-f8[0], f8[4]-f8[1]);
        double r02 = std::hypot(f8[6]-f8[0], f8[7]-f8[1]);
        double r12 = std::hypot(f8[6]-f8[3], f8[7]-f8[4]);
        int apex = 2, b0 = 0, b1 = 1;
        double e0 = std::abs(r01 - r02), e1 = std::abs(r01 - r12), e2 = std::abs(r02 - r12);
        if (e0 <= e1 && e0 <= e2) { apex = 0; b0 = 1; b1 = 2; }
        else if (e1 <= e0 && e1 <= e2) { apex = 1; b0 = 2; b1 = 0; }
        State3D fp{};
        int ord[3] = {b0, b1, apex};
        for (int k=0;k<3;++k) for (int c=0;c<3;++c) { fp[3*k+c]=f8[3*ord[k]+c]; fp[9+3*k+c]=f8[9+3*ord[k]+c]; }
        // Rotate so base (b1 - b0) lies along +x-axis
        double ang = std::atan2(fp[4] - fp[1], fp[3] - fp[0]);
        State3D s = rotate_state(fp, uz, -ang);
        double sc = std::abs(s[0]);
        double y0 = s[1] / sc;
        double vx = s[9] * std::sqrt(sc);
        double vy = s[10] * std::sqrt(sc);
        double T6 = (T8 / 6.0) / std::pow(sc, 1.5);
        // Check which sign of (vx, vy) puts body 1 (index 0) on the y-axis at t = T6
        State3D chk = make_choreo_state(y0, vx, vy, 0.0);
        integrate3d(chk, T6, EQ_MASSES, 1e-13, 200000);
        double ang_chk = std::atan2(chk[4]-chk[7], chk[3]-chk[6]);
        State3D cr = rotate_state(chk, uz, -ang_chk);
        if (std::abs(cr[0]) > 1e-3) { vx = -vx; vy = -vy; }
        std::cerr << std::setprecision(13) << "# Fig-8 isosceles gauge (apex=" << apex << "): y0=" << y0 << " vx=" << vx
                  << " vy=" << vy << " T6=" << T6 << "\n";
        double da = (argc > 2) ? std::atof(argv[2]) : 0.02;
        int nsteps = (argc > 3) ? std::atoi(argv[3]) : 40;
        bool emit_j = false;
        for (int i = 4; i < argc; ++i) if (std::string(argv[i]) == "--json") emit_j = true;
        double p[4] = {y0, vx, vy, T6};
        for (int k = 0; k <= nsteps; ++k) {
            double a = k * da;
            double ferr = 0;
            if (!refine_2d_fixed_alpha(p, a, ferr)) break;
            State3D s0 = make_choreo_state(p[0], p[1], p[2], 0.0);
            double T = 6.0 * p[3], phi = 6.0 * a;
            double cerr = choreo_third_err(s0, T, uz, 2.0 * a);
            State3D sT = s0; integrate3d(sT, T, EQ_MASSES, 1e-13, 400000);
            double ret = state_dist3d(rotate_state(sT, uz, -phi), s0);
            auto gm = measure_geom(s0, T);
            auto st = analyze_stab_choreo(s0, T, uz, phi, 2, emit_j);
            if (emit_j) {
                if (st.stable_2d && st.surv100 >= 100 && std::abs(a) > 1e-4 && ret < 1e-7) {
                    emit_choreo_json(2, "2d_choreo_rotating_fig8", s0, T, phi, cerr, ret, st, gm);
                }
                std::cerr << std::setprecision(10) << "fig8 a=" << a << " wT=" << phi << " T=" << T
                          << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                          << " marg2=" << st.margin_2d << " stab3=" << st.stable_3d
                          << " y0=" << p[0] << " vx=" << p[1] << " vy=" << p[2] << " T6=" << p[3] << "\n";
            } else {
                std::cout << a << " " << phi << " " << T << " " << st.k1 << " " << st.k2 << " "
                          << st.k3 << " " << (st.stable_2d?1:0) << " " << (st.stable_3d?1:0) << " "
                          << st.margin_2d << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << "\n" << std::flush;
            }
        }
    } else if (mode == "--continue-2d") {
        double p[4] = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5])};
        double a0 = std::atof(argv[6]), da = std::atof(argv[7]);
        int nsteps = std::atoi(argv[8]);
        bool emit_j = false;
        for (int i = 9; i < argc; ++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int k = 0; k <= nsteps; ++k) {
            double a = a0 + k * da;
            double ferr = 0;
            if (!refine_2d_fixed_alpha(p, a, ferr)) break;
            State3D s0 = make_choreo_state(p[0], p[1], p[2], 0.0);
            double T = 6.0 * p[3], phi = 6.0 * a;
            double cerr = choreo_third_err(s0, T, uz, 2.0 * a);
            State3D sT = s0; integrate3d(sT, T, EQ_MASSES, 1e-13, 400000);
            double ret = state_dist3d(rotate_state(sT, uz, -phi), s0);
            auto gm = measure_geom(s0, T);
            auto st = analyze_stab_choreo(s0, T, uz, phi, 2, emit_j);
            if (emit_j) {
                if (st.stable_2d && st.surv100 >= 100 && ret < 1e-7)
                    emit_choreo_json(2, "2d_choreo_rotating", s0, T, phi, cerr, ret, st, gm);
                std::cerr << std::setprecision(10) << "c2d a=" << a << " wT=" << phi << " T=" << T
                          << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                          << " marg2=" << st.margin_2d << " y0=" << p[0] << " vx=" << p[1]
                          << " vy=" << p[2] << " T6=" << p[3] << "\n";
            } else {
                std::cout << a << " " << phi << " " << T << " " << st.k1 << " " << st.k2 << " "
                          << st.k3 << " " << (st.stable_2d?1:0) << " " << (st.stable_3d?1:0) << " "
                          << st.margin_2d << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << "\n" << std::flush;
            }
        }
    } else if (mode == "--branch-3d") {
        double p[5] = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
        double vz0 = std::atof(argv[7]), dvz = std::atof(argv[8]);
        int nsteps = std::atoi(argv[9]);
        bool emit_j = false;
        for (int i = 10; i < argc; ++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int k = 0; k <= nsteps; ++k) {
            double vz = vz0 + k * dvz;
            if (vz < 1e-4) break;
            double ferr = 0;
            if (!refine_3d_fixed_vz(p, vz, ferr)) break;
            State3D s0 = make_choreo_state(p[0], p[1], p[2], vz);
            double T = 6.0 * p[3], phi = 6.0 * p[4];
            auto Lv = angmom3d(s0, EQ_MASSES);
            double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
            double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
            double cerr = choreo_third_err(s0, T, u, 2.0 * p[4]);
            State3D sT = s0; integrate3d(sT, T, EQ_MASSES, 1e-13, 500000);
            double ret = state_dist3d(rotate_state(sT, u, -phi), s0);
            auto gm = measure_geom(s0, T);
            auto st = analyze_stab_choreo(s0, T, u, phi, 3, emit_j);
            if (emit_j) {
                if (st.stable_3d && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7)
                    emit_choreo_json(3, "3d_choreo_helical", s0, T, phi, cerr, ret, st, gm);
                std::cerr << std::setprecision(10) << "b3d vz=" << vz << " wT=" << phi << " T=" << T
                          << " thick=" << gm.thick << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                          << ") stab=" << st.stable_3d << " marg=" << st.margin_3d
                          << " y0=" << p[0] << " vx=" << p[1] << " vy=" << p[2] << " T6=" << p[3] << " a=" << p[4] << "\n";
            } else {
                std::cout << vz << " " << phi << " " << T << " " << gm.thick << " " << gm.rmin << " "
                          << st.k1 << " " << st.k2 << " " << st.k3 << " " << (st.stable_3d?1:0) << " "
                          << st.margin_3d << " " << p[0] << " " << p[1] << " " << p[2] << " "
                          << p[3] << " " << p[4] << "\n" << std::flush;
            }
        }
    } else if (mode == "--snap-2d") {
        double p[4] = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5])};
        int pnum = std::atoi(argv[6]), qden = std::atoi(argv[7]);
        std::string label = (argc > 8) ? argv[8] : ("2d_choreo_inert_" + std::to_string(qden) + "fold");
        // wT = 6*a = 2*pi * pnum / qden  =>  a = pi * pnum / (3 * qden)
        double target_a = M_PI * pnum / (3.0 * qden);
        double ferr = 0;
        if (!refine_2d_fixed_alpha(p, target_a, ferr)) {
            std::cerr << "Snap 2D failed ferr=" << ferr << "\n"; return 2;
        }
        State3D s0 = make_choreo_state(p[0], p[1], p[2], 0.0);
        double T_rpo = 6.0 * p[3], T_inert = qden * T_rpo, phi = 6.0 * target_a;
        double cerr = choreo_third_err(s0, T_rpo, uz, 2.0 * target_a);
        State3D yq = s0; integrate3d(yq, T_inert, EQ_MASSES, 1e-13, 2000000);
        double ret_inert = state_dist3d(yq, s0);
        auto gm = measure_geom(s0, T_rpo); gm.swaps *= qden;
        auto st = analyze_stab_choreo(s0, T_rpo, uz, phi, 2, true);
        std::cerr << std::setprecision(13) << "# SNAPPED 2D " << label << " wT=" << pnum << "*2pi/" << qden
                  << " T_rpo=" << T_rpo << " T_inert=" << T_inert << " ret_inert=" << ret_inert
                  << " cerr=" << cerr << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                  << " marg=" << st.margin_2d << " surv=" << st.surv100 << "\n";
        if (st.stable_2d && st.surv100 >= 100 && ret_inert < 1e-7) {
            emit_choreo_json(2, label, s0, T_inert, phi, cerr, ret_inert, st, gm, qden, T_rpo);
        }
    } else if (mode == "--snap-3d") {
        double p[5] = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
        int pnum = std::atoi(argv[7]), qden = std::atoi(argv[8]);
        std::string label = (argc > 9) ? argv[9] : ("3d_choreo_inert_" + std::to_string(qden) + "fold");
        double target_a = M_PI * pnum / (3.0 * qden);
        double ferr = 0;
        if (!refine_3d_fixed_alpha(p, target_a, ferr)) {
            std::cerr << "Snap 3D failed ferr=" << ferr << "\n"; return 2;
        }
        State3D s0 = make_choreo_state(p[0], p[1], p[2], p[3]);
        double T_rpo = 6.0 * p[4], T_inert = qden * T_rpo, phi = 6.0 * target_a;
        auto Lv = angmom3d(s0, EQ_MASSES);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        double cerr = choreo_third_err(s0, T_rpo, u, 2.0 * target_a);
        State3D yq = s0; integrate3d(yq, T_inert, EQ_MASSES, 1e-13, 2000000);
        double ret_inert = state_dist3d(yq, s0);
        auto gm = measure_geom(s0, T_rpo); gm.swaps *= qden;
        auto st = analyze_stab_choreo(s0, T_rpo, u, phi, 3, true);
        std::cerr << std::setprecision(13) << "# SNAPPED 3D " << label << " wT=" << pnum << "*2pi/" << qden
                  << " T_rpo=" << T_rpo << " T_inert=" << T_inert << " ret_inert=" << ret_inert
                  << " cerr=" << cerr << " thick=" << gm.thick << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                  << ") marg=" << st.margin_3d << " surv=" << st.surv100 << "\n";
        if (st.stable_3d && st.surv100 >= 100 && gm.thick >= 0.05 && ret_inert < 1e-7) {
            emit_choreo_json(3, label, s0, T_inert, phi, cerr, ret_inert, st, gm, qden, T_rpo);
        }
    }
    return 0;
}
