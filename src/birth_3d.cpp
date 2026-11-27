#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

static bool dopri5_step3d(State3D& y, double& t, double& h, double t_end, const Masses& m,
                          double tol, double& rmin_track) {
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
        deriv3d(y, m, k1, &rm); if (rm < 1e-5) return false;
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt, m, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt, m, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt, m, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt, m, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt, m, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv3d(yn, m, k7, &rm); if (rm < 1e-5) return false;
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

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13, double* rmin_out = nullptr) {
    double t = 0, h = std::min(1e-3, T*0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < 400000) {
        if (!dopri5_step3d(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 80.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

// Gauge-fixed 3D perp state: x1 = -1, vz1 = 0, vz2 = A (fixed non-zero!), free = {x2, vy1, vy2, Th}
static State3D make_gauge_perp3d(double x2, double vy1, double vy2, double vz2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2) / m[2];
    double vy3 = -(m[0]*vy1 + m[1]*vy2) / m[2];
    double vz3 = -(m[1]*vz2) / m[2];
    return {x1, 0, 0, x2, 0, 0, x3, 0, 0,
            0, vy1, 0, 0, vy2, vz2, 0, vy3, vz3};
}

// Refine 4 variables {x2, vy1, vy2, Th} with fixed vz2 to minimize 6 residuals {y1, z1, y2, z2, vx1, vx2}
static bool refine_branch(double p[4], double vz2, const Masses& m, double& ferr) {
    double lambda = 1e-3;
    auto eval = [&](const double q[4], std::array<double, 6>& F) -> bool {
        if (q[3] < 0.5 || q[3] > 60.0) return false;
        State3D y = make_gauge_perp3d(q[0], q[1], q[2], vz2, m);
        double rm = 1e300;
        if (!integrate3d(y, q[3], m, 1e-12, &rm)) return false;
        F = {y[1], y[2], y[4], y[5], y[9], y[12]};
        return true;
    };
    std::array<double, 6> F;
    if (!eval(p, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;

    for (int it = 0; it < 40; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[6][4];
        std::array<double, 6> Fp, Fm;
        for (int col = 0; col < 4; ++col) {
            double pp[4]={p[0],p[1],p[2],p[3]}, pm[4]={p[0],p[1],p[2],p[3]};
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double A[4][5] = {{0}};
        for (int r=0;r<4;++r) {
            for (int c=0;c<4;++c) { double s=0; for(int i=0;i<6;++i) s+=J[i][r]*J[i][c]; A[r][c]=s; }
            A[r][r] += lambda*(1.0 + A[r][r]);
            double s=0; for(int i=0;i<6;++i) s-=J[i][r]*F[i]; A[r][4]=s;
        }
        bool sing = false;
        for (int col=0;col<4;++col) {
            int piv=col;
            for (int r=col+1;r<4;++r) if (std::abs(A[r][col])>std::abs(A[piv][col])) piv=r;
            if (std::abs(A[piv][col])<1e-18) { sing=true; break; }
            for (int c=col;c<5;++c) std::swap(A[col][c], A[piv][c]);
            for (int r=col+1;r<4;++r) {
                double f=A[r][col]/A[col][col];
                for (int c=col;c<5;++c) A[r][c]-=f*A[col][c];
            }
        }
        if (sing) { lambda*=10; continue; }
        double d[4]={0};
        for (int r=3;r>=0;--r) { double s=A[r][4]; for(int c=r+1;c<4;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r]; }
        double sn=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
        if (sn > 0.3) for(int i=0;i<4;++i) d[i]*=(0.3/sn);
        double np[4]={p[0]+d[0], p[1]+d[1], p[2]+d[2], p[3]+d[3]};
        std::array<double, 6> nF;
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(double v:nF) nc+=v*v;
        if (nc < cost) { for(int i=0;i<4;++i) p[i]=np[i]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

static double min_eig_3x3(double C[3][3]) {
    double I1 = C[0][0]+C[1][1]+C[2][2];
    double I2 = C[0][0]*C[1][1]+C[0][0]*C[2][2]+C[1][1]*C[2][2] - C[0][1]*C[0][1]-C[0][2]*C[0][2]-C[1][2]*C[1][2];
    double I3 = C[0][0]*(C[1][1]*C[2][2]-C[1][2]*C[1][2]) - C[0][1]*(C[0][1]*C[2][2]-C[1][2]*C[0][2]) + C[0][2]*(C[0][1]*C[1][2]-C[1][1]*C[0][2]);
    double p = I2 - I1*I1/3.0, q = -2.0*I1*I1*I1/27.0 + I1*I2/3.0 - I3;
    double r = std::sqrt(std::max(0.0, -p/3.0));
    if (r < 1e-15) return std::max(0.0, I1/3.0);
    double phi = std::acos(std::clamp(-q/(2.0*r*r*r), -1.0, 1.0)) / 3.0;
    return std::max(0.0, std::min({I1/3.0+2*r*std::cos(phi), I1/3.0+2*r*std::cos(phi+2*M_PI/3), I1/3.0+2*r*std::cos(phi+4*M_PI/3)}));
}

int main(int argc, char** argv) {
    // Args: m1 m2 m3 x2 vy1 vy2 Thalf mult(1 or 2) vz2_amp
    Masses m = {std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
    double p[4] = {std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]) * std::atof(argv[8])};
    double vz2 = std::atof(argv[9]);
    double ferr = 0;
    bool ok = refine_branch(p, vz2, m, ferr);
    std::cerr << "vz2=" << vz2 << " ok=" << ok << " ferr=" << ferr << " Th=" << p[3] << "\n";
    if (!ok) return 0;

    double T = 2.0 * p[3];
    State3D y0 = make_gauge_perp3d(p[0], p[1], p[2], vz2, m);
    State3D yT = y0;
    double rmin = 1e300;
    integrate3d(yT, T, m, 1e-13, &rmin);
    double ret = state_dist3d(yT, y0);

    // 3D thickness & swaps
    const int N = 320;
    double dt = T/N, C[3][3]={{0}}, rmax=0;
    int prev_cl=-1, swaps=0, cnt=0;
    State3D yc = y0;
    for (int k=0;k<=N;++k) {
        double r12=std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
        double r13=std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
        double r23=std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
        rmax = std::max({rmax, r12, r13, r23});
        int cl = (r12<=r13 && r12<=r23)?0:((r13<=r23)?1:2);
        if (prev_cl!=-1 && cl!=prev_cl) swaps++;
        prev_cl = cl;
        for (int b=0;b<3;++b) {
            double q[3]={yc[3*b], yc[3*b+1], yc[3*b+2]};
            for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c]+=q[r]*q[c];
            cnt++;
        }
        if (k<N) integrate3d(yc, dt, m, 1e-11);
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c]/=cnt;
    double thick = std::sqrt(min_eig_3x3(C));

    // Full 18x18 Floquet
    double M[18][18], eps=1e-7;
    for (int col=0;col<18;++col) {
        State3D yp=y0, ym=y0; yp[col]+=eps; ym[col]-=eps;
        integrate3d(yp, T, m, 1e-13); integrate3d(ym, T, m, 1e-13);
        for (int r=0;r<18;++r) M[r][col]=(yp[r]-ym[r])/(2*eps);
    }
    double M2[18][18]={{0}};
    for (int i=0;i<18;++i) for (int k=0;k<18;++k) for (int j=0;j<18;++j) M2[i][j]+=M[i][k]*M[k][j];
    double tr1=0, tr2=0, tr3=0;
    for (int i=0;i<18;++i) { tr1+=M[i][i]; tr2+=M2[i][i]; for(int j=0;j<18;++j) tr3+=M2[i][j]*M[j][i]; }
    double p1=tr1-12.0, p2=tr2-6.0, p3=tr3-12.0+3.0*p1;
    double e1=p1, e2=0.5*(e1*p1-p2), e3=(e2*p1-e1*p2+p3)/3.0;
    double Q=(e1*e1-3.0*e2)/9.0, R=(2.0*e1*e1*e1-9.0*e1*e2+27.0*e3)/54.0, disc=Q*Q*Q-R*R;
    double k1=0, k2=0, k3=0, margin=-1;
    bool stable=false;
    if (disc >= -1e-4 && Q > 0) {
        double sqQ=std::sqrt(Q), th=std::acos(std::clamp(R/(sqQ*sqQ*sqQ),-1.0,1.0));
        k1=2*sqQ*std::cos(th/3)+e1/3; k2=2*sqQ*std::cos((th+2*M_PI)/3)+e1/3; k3=2*sqQ*std::cos((th+4*M_PI)/3)+e1/3;
        double mk=std::max({std::abs(k1), std::abs(k2), std::abs(k3)});
        margin=2.0-mk; stable=(mk < 2.0-1e-4);
    }

    // 100T empirical 3D kick test
    State3D yp = y0;
    double kick = 1e-5;
    yp[9]+=kick; yp[12]-=kick*m[0]/m[1];
    yp[10]+=kick; yp[13]-=kick*m[0]/m[1];
    yp[11]+=kick; yp[14]-=kick*m[0]/m[1];
    double md=0; int surv=0;
    for (int k=0;k<100;++k) {
        if (!integrate3d(yp, T, m, 1e-11)) break;
        double d=state_dist3d(yp, y0); if (d>md) md=d; if (d>5.0) break;
        surv=k+1;
    }

    std::cerr << "  3D BIRTH! thick=" << thick << " swaps=" << swaps << " ret=" << ret
              << " k=(" << k1 << "," << k2 << "," << k3 << ") stable=" << stable
              << " margin=" << margin << " surv100=" << surv << " drift=" << md << "\n";

    if (stable && surv >= 100 && thick > 0.05) {
        double E = energy3d(y0, m);
        auto Lv = angmom3d(y0, m);
        double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        std::cout << std::setprecision(13)
                  << "{\"dim\":3,\"subspace\":\"3d_vertical_bifurcation\",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "],\"y0\":[";
        for (int i=0;i<18;++i) std::cout << y0[i] << (i<17?",":"");
        std::cout << "],\"T\":" << T << ",\"E\":" << E << ",\"L\":" << Lmag
                  << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]"
                  << ",\"ret_err\":" << ret << ",\"k1\":" << k1 << ",\"k2\":" << k2 << ",\"k3\":" << k3
                  << ",\"stab_margin\":" << margin << ",\"survived_100T\":" << surv
                  << ",\"drift_100T\":" << md << ",\"rmin\":" << rmin << ",\"rmax\":" << rmax
                  << ",\"thickness_3d\":" << thick << ",\"partner_swaps\":" << swaps
                  << ",\"democratic\":" << (swaps >= 2 ? "true" : "false") << "}\n";
    }
    return 0;
}
