#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using State = std::array<double, 12>;
using Masses = std::array<double, 3>;

static inline void deriv(const State& s, const Masses& m, State& ds, double* rmin_out = nullptr) {
    ds[0] = s[6];  ds[1] = s[7];
    ds[2] = s[8];  ds[3] = s[9];
    ds[4] = s[10]; ds[5] = s[11];
    double ax[3] = {0}, ay[3] = {0}, rmin2 = 1e300;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[2*j] - s[2*i], dy = s[2*j+1] - s[2*i+1];
            double r2 = dx*dx + dy*dy;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            ax[i] += m[j]*dx*ir3; ay[i] += m[j]*dy*ir3;
            ax[j] -= m[i]*dx*ir3; ay[j] -= m[i]*dy*ir3;
        }
    ds[6]=ax[0]; ds[7]=ay[0]; ds[8]=ax[1]; ds[9]=ay[1]; ds[10]=ax[2]; ds[11]=ay[2];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static bool dopri5_step(State& y, double& t, double& h, double t_end, const Masses& m, double tol) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;
    State k1,k2,k3,k4,k5,k6,k7,yt,yn;
    double rm = 1e300;
    for (int att = 0; att < 25; ++att) {
        deriv(y, m, k1, &rm); if (rm < 1e-6) return false;
        for (int i=0;i<12;++i) yt[i]=y[i]+h*a21*k1[i]; deriv(yt, m, k2);
        for (int i=0;i<12;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv(yt, m, k3);
        for (int i=0;i<12;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(yt, m, k4);
        for (int i=0;i<12;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(yt, m, k5);
        for (int i=0;i<12;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(yt, m, k6);
        for (int i=0;i<12;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv(yn, m, k7, &rm); if (rm < 1e-6) return false;
        double err = 0;
        for (int i=0;i<12;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err/12.0);
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

static bool integrate_to(State& y, double T, const Masses& m, double tol = 1e-13) {
    double t = 0, h = std::min(1e-3, T*0.01);
    int steps = 0;
    while (t < T - 1e-15 && steps < 300000) {
        if (!dopri5_step(y, t, h, T, m, tol)) return false;
        ++steps;
        for (int i=0;i<6;++i) if (std::abs(y[i]) > 80.0) return false;
    }
    return t >= T - 1e-13;
}

static State make_col(double vx, double vy, double m3) {
    return {-1, 0, 1, 0, 0, 0, vx, vy, vx, vy, -2*vx/m3, -2*vy/m3};
}

static bool refine_col(double& vx, double& vy, double& T, double m3) {
    Masses m = {1.0, 1.0, m3};
    double lambda = 1e-5;
    auto eval = [&](double cvx, double cvy, double cT, State& F) -> bool {
        State y0 = make_col(cvx, cvy, m3), yT = y0;
        if (!integrate_to(yT, cT, m, 1e-13)) return false;
        for (int i=0;i<12;++i) F[i] = yT[i] - y0[i];
        return true;
    };
    State F; if (!eval(vx, vy, T, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;
    for (int it = 0; it < 25; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7;
        double J[12][3];
        State Fp, Fm;
        if (!eval(vx+eps, vy, T, Fp) || !eval(vx-eps, vy, T, Fm)) return false;
        for (int i=0;i<12;++i) J[i][0] = (Fp[i]-Fm[i])/(2*eps);
        if (!eval(vx, vy+eps, T, Fp) || !eval(vx, vy-eps, T, Fm)) return false;
        for (int i=0;i<12;++i) J[i][1] = (Fp[i]-Fm[i])/(2*eps);
        State y0 = make_col(vx, vy, m3), yT, fT;
        for (int i=0;i<12;++i) yT[i] = y0[i] + F[i];
        deriv(yT, m, fT);
        for (int i=0;i<12;++i) J[i][2] = fT[i];
        double M[3][4]={{0}};
        for (int r=0;r<3;++r) {
            for (int c=0;c<3;++c) { double s=0; for(int i=0;i<12;++i) s+=J[i][r]*J[i][c]; M[r][c]=s; }
            M[r][r] += lambda*(1.0+M[r][r]);
            double s=0; for(int i=0;i<12;++i) s-=J[i][r]*F[i]; M[r][3]=s;
        }
        for (int col=0;col<3;++col) {
            int piv=col;
            for (int r=col+1;r<3;++r) if (std::abs(M[r][col])>std::abs(M[piv][col])) piv=r;
            if (std::abs(M[piv][col])<1e-18) return false;
            for (int c=col;c<4;++c) std::swap(M[col][c], M[piv][c]);
            for (int r=col+1;r<3;++r) {
                double f=M[r][col]/M[col][col];
                for (int c=col;c<4;++c) M[r][c]-=f*M[col][c];
            }
        }
        double d[3]={0};
        for (int r=2;r>=0;--r) { double s=M[r][3]; for(int c=r+1;c<3;++c) s-=M[r][c]*d[c]; d[r]=s/M[r][r]; }
        State nF;
        if (!eval(vx+d[0], vy+d[1], T+d[2], nF)) { lambda*=10; continue; }
        double nc=0; for(double v:nF) nc+=v*v;
        if (nc < cost) { vx+=d[0]; vy+=d[1]; T+=d[2]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.2); }
        else { lambda*=5; if (lambda>1e6) break; }
    }
    return std::sqrt(cost) < 1e-10;
}

static void floquet(const State& y0, double T, const Masses& m, double& k1, double& k2, bool& stable) {
    double M[12][12], eps=1e-7;
    for (int col=0;col<12;++col) {
        State yp=y0, ym=y0; yp[col]+=eps; ym[col]-=eps;
        integrate_to(yp, T, m); integrate_to(ym, T, m);
        for (int r=0;r<12;++r) M[r][col] = (yp[r]-ym[r])/(2*eps);
    }
    double tr1=0, tr2=0;
    for (int i=0;i<12;++i) tr1 += M[i][i];
    for (int i=0;i<12;++i) for (int j=0;j<12;++j) tr2 += M[i][j]*M[j][i];
    double s = tr1 - 8.0, p = 0.5*(s*s - (tr2 - 4.0)), disc = s*s - 4*p;
    if (disc >= -1e-5) {
        double sd = std::sqrt(std::max(0.0, disc));
        k1 = 0.5*(s+sd); k2 = 0.5*(s-sd);
        stable = (std::abs(k1) < 2.0-1e-4) && (std::abs(k2) < 2.0-1e-4);
    } else {
        k1 = 0.5*s; k2 = 0.5*std::sqrt(-disc); stable = false;
    }
}

int main(int argc, char** argv) {
    if (argc < 7) {
        std::cerr << "Usage: continuation <vx> <vy> <T> <m3_start> <dm3> <nsteps>\n";
        return 1;
    }
    double vx = std::atof(argv[1]), vy = std::atof(argv[2]), T = std::atof(argv[3]);
    double m3 = std::atof(argv[4]), dm3 = std::atof(argv[5]);
    int nsteps = std::atoi(argv[6]);
    std::cout << std::setprecision(12);
    for (int step = 0; step <= nsteps; ++step) {
        double cm3 = m3 + step * dm3;
        if (cm3 <= 0.05) break;
        if (!refine_col(vx, vy, T, cm3)) break;
        double k1=0, k2=0; bool stab=false;
        floquet(make_col(vx, vy, cm3), T, {1.0, 1.0, cm3}, k1, k2, stab);
        double margin = stab ? (2.0 - std::max(std::abs(k1), std::abs(k2))) : -1.0;
        std::cout << cm3 << " " << vx << " " << vy << " " << T << " " << k1 << " " << k2 << " "
                  << (stab ? 1 : 0) << " " << margin << "\n" << std::flush;
    }
    return 0;
}
