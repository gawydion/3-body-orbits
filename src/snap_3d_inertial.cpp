#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0};
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            double r2 = dx*dx+dy*dy+dz*dz, ir = 1.0/std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]+=m[j]*dx*ir3; a[3*i+1]+=m[j]*dy*ir3; a[3*i+2]+=m[j]*dz*ir3;
            a[3*j]-=m[i]*dx*ir3; a[3*j+1]-=m[i]*dy*ir3; a[3*j+2]-=m[i]*dz*ir3;
        }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
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
    for (int att=0;att<25;++att) {
        deriv3d(y, m, k1);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt, m, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt, m, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt, m, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt, m, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt, m, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); deriv3d(yn, m, k7);
        double err=0;
        for (int i=0;i<18;++i) {
            double ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc=tol*(1.0+std::max(std::abs(y[i]), std::abs(yn[i])));
            double r=ei/sc; err+=r*r;
        }
        err=std::sqrt(err/18.0);
        if (err<=1.0) {
            y=yn; t+=h;
            h=std::min(h*((err==0)?5.0:std::min(5.0,std::max(0.2,0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h*=std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h<1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-13) {
    double t=0, h=std::min(1e-3, T*0.01);
    int steps=0;
    while (t < T - 1e-15 && steps < 800000) {
        if (!dopri5_step3d(y, t, h, T, m, tol)) return false;
        ++steps;
    }
    return t >= T - 1e-13;
}

static State3D make_state(double x2, double vy1, double vy2, double vz2, const Masses& m) {
    double x1=-1.0, x3=-(m[0]*x1+m[1]*x2)/m[2];
    double vy3=-(m[0]*vy1+m[1]*vy2)/m[2], vz3=-(m[1]*vz2)/m[2];
    return {x1,0,0, x2,0,0, x3,0,0, 0,vy1,0, 0,vy2,vz2, 0,vy3,vz3};
}

int main(int argc, char** argv) {
    // Args: target_2a_num target_2a_den q_fold m3 x2 vy1 vy2 vz2 Th
    // Example 7-fold: -8 7 7 0.5 ...  (so 2*alpha = -8*pi/7 => alpha = -4*pi/7, 7*2*alpha = -8*pi)
    double num = std::atof(argv[1]), den = std::atof(argv[2]);
    int q_fold = std::atoi(argv[3]);
    Masses m = {1.0, 1.0, std::atof(argv[4])};
    double p[5] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9])};
    double target_a = 0.5 * M_PI * num / den;
    double lambda = 1e-5;

    auto eval = [&](const double q[5], std::array<double, 6>& F) -> bool {
        State3D y0 = make_state(q[0], q[1], q[2], q[3], m);
        auto Lv = angmom3d(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State3D y = y0;
        if (!integrate3d(y, q[4], m, 1e-13)) return false;
        State3D yt = rotate_state(y, u, -target_a);
        F = {yt[1], yt[2], yt[4], yt[5], yt[9], yt[12]};
        return true;
    };

    std::array<double, 6> F;
    eval(p, F);
    double cost = 0; for (double v : F) cost += v*v;
    for (int it = 0; it < 25; ++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[6][5];
        std::array<double, 6> Fp, Fm;
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            eval(pp, Fp); eval(pm, Fm);
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
            for (int c=col;c<6;++c) std::swap(A[col][c], A[piv][c]);
            for (int r=col+1;r<5;++r) { double f=A[r][col]/A[col][col]; for(int c=col;c<6;++c) A[r][c]-=f*A[col][c]; }
        }
        double d[5]={0};
        for (int r=4;r>=0;--r) { double s=A[r][5]; for(int c=r+1;c<5;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r]; }
        double np[5]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        std::array<double, 6> nF; eval(np, nF);
        double nc=0; for(double v:nF) nc+=v*v;
        if (nc < cost) { for(int i=0;i<5;++i) p[i]=np[i]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.2); }
        else { lambda*=5; }
    }

    State3D y0 = make_state(p[0], p[1], p[2], p[3], m);
    double T_rpo = 2.0 * p[4];
    double T_inert = q_fold * T_rpo;
    State3D yq = y0;
    integrate3d(yq, T_inert, m, 1e-13);
    double errq = 0; for(int i=0;i<18;++i) errq += (yq[i]-y0[i])*(yq[i]-y0[i]);
    errq = std::sqrt(errq);

    std::cout << std::setprecision(13);
    std::cout << q_fold << "-FOLD (2a=" << num << "*pi/" << den << ", m3=" << m[2] << "):"
              << " half_err=" << std::sqrt(cost)
              << " T_rpo=" << T_rpo << " T_inert=" << T_inert
              << " |y(T_inert)-y(0)|=" << errq
              << " vz2=" << p[3] << "\n";
    std::cout << "  y0=[";
    for (int i=0;i<18;++i) std::cout << y0[i] << (i<17?",":"");
    std::cout << "]\n";
    return 0;
}
