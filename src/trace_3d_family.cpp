#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0}, rmin2 = 1e300;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            double r2 = dx*dx+dy*dy+dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0/std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]+=m[j]*dx*ir3; a[3*i+1]+=m[j]*dy*ir3; a[3*i+2]+=m[j]*dz*ir3;
            a[3*j]-=m[i]*dx*ir3; a[3*j+1]-=m[i]*dy*ir3; a[3*j+2]-=m[i]*dz*ir3;
        }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline std::array<double, 3> angmom3d(const State3D& s, const Masses& m) {
    double Lx=0, Ly=0, Lz=0;
    for (int i=0;i<3;++i) {
        Lx += m[i]*(s[3*i+1]*s[11+3*i] - s[3*i+2]*s[10+3*i]);
        Ly += m[i]*(s[3*i+2]*s[9+3*i]  - s[3*i]*s[11+3*i]);
        Lz += m[i]*(s[3*i]*s[10+3*i]   - s[3*i+1]*s[9+3*i]);
    }
    return {Lx, Ly, Lz};
}

static inline double energy3d(const State3D& s, const Masses& m) {
    double K=0, U=0;
    for (int i=0;i<3;++i) K += 0.5*m[i]*(s[9+3*i]*s[9+3*i]+s[10+3*i]*s[10+3*i]+s[11+3*i]*s[11+3*i]);
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
        double dx=s[3*j]-s[3*i], dy=s[3*j+1]-s[3*i+1], dz=s[3*j+2]-s[3*i+2];
        U -= m[i]*m[j]/std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    return K + U;
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
    for (int att=0;att<25;++att) {
        deriv3d(y, m, k1, &rm); if (rm < 1e-4) return false;
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt, m, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt, m, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt, m, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt, m, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt, m, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv3d(yn, m, k7, &rm); if (rm < 1e-4) return false;
        double err=0;
        for (int i=0;i<18;++i) {
            double ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc=tol*(1.0+std::max(std::abs(y[i]), std::abs(yn[i])));
            double r=ei/sc; err+=r*r;
        }
        err=std::sqrt(err/18.0);
        if (err<=1.0) {
            y=yn; t+=h; if (rm < rmin_track) rmin_track = rm;
            h=std::min(h*((err==0)?5.0:std::min(5.0,std::max(0.2,0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h*=std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h<1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13, double* rmin_out = nullptr) {
    double t=0, h=std::min(1e-3, T*0.01), rmin_track=1e300;
    int steps=0;
    while (t < T - 1e-15 && steps < 400000) {
        if (!dopri5_step3d(y, t, h, T, m, tol, rmin_track)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 80.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

static State3D make_state(double x2, double vy1, double vy2, double vz2, const Masses& m) {
    double x1=-1.0, x3=-(m[0]*x1+m[1]*x2)/m[2];
    double vy3=-(m[0]*vy1+m[1]*vy2)/m[2], vz3=-(m[1]*vz2)/m[2];
    return {x1,0,0, x2,0,0, x3,0,0, 0,vy1,0, 0,vy2,vz2, 0,vy3,vz3};
}

// Refine 5 vars p = {x2, vy1, vy2, Th, alpha} at fixed vz2 and m
static bool refine_fixed_vz(double p[5], double vz2, const Masses& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[5], std::array<double, 6>& F) -> bool {
        if (q[3] < 0.5 || q[3] > 30.0) return false;
        State3D y0 = make_state(q[0], q[1], q[2], vz2, m);
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[3], m, 1e-12)) return false;
        State3D yt = rotate_state(y, u, -q[4]);
        F = {yt[1], yt[2], yt[4], yt[5], yt[9], yt[12]};
        return true;
    };
    std::array<double, 6> F;
    if (!eval(p, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;
    for (int it = 0; it < 25; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = 1e-7, J[6][5];
        std::array<double, 6> Fp, Fm;
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double A[5][6]={{0}};
        for (int r=0;r<5;++r) {
            for (int c=0;c<5;++c) { double s=0; for(int i=0;i<6;++i) s+=J[i][r]*J[i][c]; A[r][c]=s; }
            A[r][r] += lambda*(1.0+A[r][r]);
            double s=0; for(int i=0;i<6;++i) s-=J[i][r]*F[i]; A[r][5]=s;
        }
        for (int col=0;col<5;++col) {
            int piv=col;
            for (int r=col+1;r<5;++r) if (std::abs(A[r][col])>std::abs(A[piv][col])) piv=r;
            if (std::abs(A[piv][col]) < 1e-18) return false;
            for (int c=col;c<6;++c) std::swap(A[col][c], A[piv][c]);
            for (int r=col+1;r<5;++r) { double f=A[r][col]/A[col][col]; for(int c=col;c<6;++c) A[r][c]-=f*A[col][c]; }
        }
        double d[5]={0};
        for (int r=4;r>=0;--r) { double s=A[r][5]; for(int c=r+1;c<5;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r]; }
        double sn=0; for(int i=0;i<5;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<5;++i) d[i]*=(0.25/sn);
        double np[5]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        std::array<double, 6> nF;
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(double v:nF) nc+=v*v;
        if (nc < cost) { for(int i=0;i<5;++i) p[i]=np[i]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.25); }
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

static void floquet_rpo(const State3D& y0, double T, const double u[3], double two_alpha, const Masses& m,
                        double& k1, double& k2, double& k3, bool& stable, double& margin) {
    double M[18][18], eps = 1e-7;
    for (int col=0;col<18;++col) {
        State3D yp=y0, ym=y0; yp[col]+=eps; ym[col]-=eps;
        integrate3d(yp, T, m, 1e-13); integrate3d(ym, T, m, 1e-13);
        yp = rotate_state(yp, u, -two_alpha);
        ym = rotate_state(ym, u, -two_alpha);
        for (int r=0;r<18;++r) M[r][col]=(yp[r]-ym[r])/(2*eps);
    }
    double M2[18][18]={{0}};
    for (int i=0;i<18;++i) for (int k=0;k<18;++k) for (int j=0;j<18;++j) M2[i][j]+=M[i][k]*M[k][j];
    double tr1=0, tr2=0, tr3=0;
    for (int i=0;i<18;++i) { tr1+=M[i][i]; tr2+=M2[i][i]; for(int j=0;j<18;++j) tr3+=M2[i][j]*M[j][i]; }
    double c1=std::cos(two_alpha), c2=std::cos(2*two_alpha), c3=std::cos(3*two_alpha);
    double p1=tr1-(6+6*c1), p2=tr2-(6+6*c2)+6, p3=tr3-(6+6*c3)+3*p1;
    double e1=p1, e2=0.5*(e1*p1-p2), e3=(e2*p1-e1*p2+p3)/3.0;
    double Q=(e1*e1-3*e2)/9.0, R=(2*e1*e1*e1-9*e1*e2+27*e3)/54.0, disc=Q*Q*Q-R*R;
    stable = false; margin = -1.0; k1=k2=k3=0;
    if (disc >= -1e-4 && Q > 0) {
        double sqQ=std::sqrt(Q), th=std::acos(std::clamp(R/(sqQ*sqQ*sqQ),-1.0,1.0));
        k1=2*sqQ*std::cos(th/3)+e1/3; k2=2*sqQ*std::cos((th+2*M_PI)/3)+e1/3; k3=2*sqQ*std::cos((th+4*M_PI)/3)+e1/3;
        double mk=std::max({std::abs(k1), std::abs(k2), std::abs(k3)});
        margin = 2.0 - mk; stable = (mk < 2.0 - 1e-4);
    }
}

int main(int argc, char** argv) {
    // Args: mode ("vz" or "m3") step_val nsteps
    std::string mode = (argc > 1) ? argv[1] : "vz";
    double step_val = (argc > 2) ? std::atof(argv[2]) : 0.02;
    int nsteps = (argc > 3) ? std::atoi(argv[3]) : 30;

    // Seed from 3D-001:
    double m3 = 0.5;
    double vz2 = 0.6732282680926;
    double p[5] = {0.9479320264499, 0.2294279990126, -0.4557531444849, 9.146435729715, -3.769216675908 / 2.0};

    std::cout << std::setprecision(10);
    for (int s = 0; s <= nsteps; ++s) {
        double cm3 = (mode == "m3") ? (m3 + s * step_val) : m3;
        double cvz = (mode == "vz") ? (vz2 + s * step_val) : vz2;
        if (cm3 < 0.08 || cvz < 0.02) break;
        Masses m = {1.0, 1.0, cm3};
        double ferr = 0;
        if (!refine_fixed_vz(p, cvz, m, ferr)) break;

        State3D y0 = make_state(p[0], p[1], p[2], cvz, m);
        double T = 2.0 * p[3], two_a = 2.0 * p[4];
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};

        // 3D thickness & rmin
        const int N = 200;
        double dt = T / N, C[3][3]={{0}}, rmin = 1e300;
        int prev_cl=-1, swaps=0, cnt=0;
        State3D yc = y0;
        for (int k=0;k<=N;++k) {
            double r12=std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
            double r13=std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
            double r23=std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
            rmin = std::min({rmin, r12, r13, r23});
            int cl = (r12<=r13 && r12<=r23)?0:((r13<=r23)?1:2);
            if (prev_cl!=-1 && cl!=prev_cl) swaps++;
            prev_cl = cl;
            for (int b=0;b<3;++b) {
                double qv[3]={yc[3*b], yc[3*b+1], yc[3*b+2]};
                for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c]+=qv[r]*qv[c];
                cnt++;
            }
            if (k < N) integrate3d(yc, dt, m, 1e-11);
        }
        for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c]/=cnt;
        double thick = std::sqrt(min_eig_3x3(C));

        double k1, k2, k3, margin; bool stab;
        floquet_rpo(y0, T, u, two_a, m, k1, k2, k3, stab, margin);

        std::cout << cm3 << " " << cvz << " " << T << " " << two_a << " " << thick << " "
                  << rmin << " " << swaps << " " << k1 << " " << k2 << " " << k3 << " "
                  << (stab ? 1 : 0) << " " << margin << " "
                  << p[0] << " " << p[1] << " " << p[2] << "\n" << std::flush;
    }
    return 0;
}
