#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using ld = long double;
using State3D = std::array<ld, 18>;
using Masses = std::array<ld, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds) {
    for (int i=0;i<9;++i) ds[i] = s[9+i];
    ld a[9]={0};
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
        ld dx=s[3*j]-s[3*i], dy=s[3*j+1]-s[3*i+1], dz=s[3*j+2]-s[3*i+2];
        ld r2=dx*dx+dy*dy+dz*dz, ir=1.0L/std::sqrt(r2), ir3=ir*ir*ir;
        a[3*i]+=m[j]*dx*ir3; a[3*i+1]+=m[j]*dy*ir3; a[3*i+2]+=m[j]*dz*ir3;
        a[3*j]-=m[i]*dx*ir3; a[3*j+1]-=m[i]*dy*ir3; a[3*j+2]-=m[i]*dz*ir3;
    }
    for (int i=0;i<9;++i) ds[9+i]=a[i];
}

static inline void rot_z(State3D& s, ld ang) {
    ld c=std::cos(ang), sn=std::sin(ang);
    for (int b=0;b<6;++b) {
        ld x=s[3*b], y=s[3*b+1];
        s[3*b] = x*c - y*sn;
        s[3*b+1] = x*sn + y*c;
    }
}

static bool integrate3d(State3D& y, ld T, const Masses& m, ld tol=1e-13L) {
    static constexpr ld a21=1.0L/5, a31=3.0L/40, a32=9.0L/40,
        a41=44.0L/45, a42=-56.0L/15, a43=32.0L/9,
        a51=19372.0L/6561, a52=-25360.0L/2187, a53=64448.0L/6561, a54=-212.0L/729,
        a61=9017.0L/3168, a62=-355.0L/33, a63=46732.0L/5247, a64=49.0L/176, a65=-5103.0L/18656,
        b1=35.0L/384, b3=500.0L/1113, b4=125.0L/192, b5=-2187.0L/6784, b6=11.0L/84,
        e1=71.0L/57600, e3=-71.0L/16695, e4=71.0L/1920, e5=-17253.0L/339200, e6=22.0L/525, e7=-1.0L/40;
    ld t=0, h=std::min(1e-3L, T*0.01L);
    int steps=0;
    State3D k1,k2,k3,k4,k5,k6,k7,yt,yn;
    while (t < T - 1e-16L && steps < 400000) {
        if (t+h>T) h=T-t;
        bool ok=false;
        for (int att=0;att<25;++att) {
            deriv3d(y,m,k1);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt,m,k2);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt,m,k3);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt,m,k4);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt,m,k5);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt,m,k6);
            for(int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); deriv3d(yn,m,k7);
            ld err=0;
            for(int i=0;i<18;++i) {
                ld ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
                ld sc=tol*(1.0L+std::max(std::abs(y[i]),std::abs(yn[i])));
                ld r=ei/sc; err+=r*r;
            }
            err=std::sqrt(err/18.0L);
            if (err<=1.0L) { y=yn; t+=h; h=std::min(h*((err==0)?5.0L:std::min(5.0L,std::max(0.2L,0.9L*std::pow(err,-0.2L)))), 0.05L); ok=true; break; }
            h*=std::max(0.1L, 0.9L*std::pow(err,-0.25L));
            if (h<1e-14L) return false;
        }
        if (!ok) return false;
        ++steps;
    }
    return t >= T - 1e-13L;
}

static State3D make_planar(ld x2, ld vy1, ld vy2, const Masses& m) {
    ld x1 = -1.0L, x3 = -(m[0]*x1 + m[1]*x2)/m[2];
    ld vy3 = -(m[0]*vy1 + m[1]*vy2)/m[2];
    return {x1,0,0, x2,0,0, x3,0,0, 0,vy1,0, 0,vy2,0, 0,vy3,0};
}

static bool refine_planar_rpo(ld& x2, ld& vy1, ld vy2_fixed, ld& Th, ld& alpha, const Masses& m) {
    // 4 vars: x2, vy1, Th, alpha; 4 residuals at t=Th after R_z(-alpha): y1, y2, vx1, vx2
    ld p[4] = {x2, vy1, Th, alpha};
    auto eval = [&](const ld q[4], ld F[4]) {
        State3D y = make_planar(q[0], q[1], vy2_fixed, m);
        integrate3d(y, q[2], m, 1e-13L);
        rot_z(y, -q[3]);
        F[0]=y[1]; F[1]=y[4]; F[2]=y[9]; F[3]=y[12];
    };
    ld F[4]; eval(p, F);
    for (int it=0; it<20; ++it) {
        ld fn = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]+F[3]*F[3]);
        if (fn < 1e-13L) break;
        ld eps = 1e-7L, J[4][4], Fp[4], Fm[4];
        for (int c=0;c<4;++c) {
            ld pp[4]={p[0],p[1],p[2],p[3]}, pm[4]={p[0],p[1],p[2],p[3]};
            pp[c]+=eps; pm[c]-=eps;
            eval(pp,Fp); eval(pm,Fm);
            for (int r=0;r<4;++r) J[r][c]=(Fp[r]-Fm[r])/(2*eps);
        }
        ld A[4][5];
        for (int r=0;r<4;++r) { for(int c=0;c<4;++c) A[r][c]=J[r][c]; A[r][4]=-F[r]; }
        for (int c=0;c<4;++c) {
            int piv=c; for(int r=c+1;r<4;++r) if(std::abs(A[r][c])>std::abs(A[piv][c])) piv=r;
            if (std::abs(A[piv][c])<1e-18L) return false;
            for (int k=c;k<5;++k) std::swap(A[c][k], A[piv][k]);
            for (int r=c+1;r<4;++r) { ld f=A[r][c]/A[c][c]; for(int k=c;k<5;++k) A[r][k]-=f*A[c][k]; }
        }
        ld d[4]={0};
        for (int r=3;r>=0;--r) { ld s=A[r][4]; for(int c=r+1;c<4;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r]; }
        for (int i=0;i<4;++i) p[i]+=d[i];
        eval(p, F);
    }
    x2=p[0]; vy1=p[1]; Th=p[2]; alpha=p[3];
    return std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]+F[3]*F[3]) < 1e-10L;
}

int main() {
    // Seed from 2D-DEM-001 (m=[1,1,0.5], sign-flipped so L > 0!)
    // Original 2D-DEM-001: x2=0.9941458062522, vy1=0.07446872307688, vy2=-0.7818621047497, T=24.3591072214 (7-fold!)
    // Flip sign of vy so L > 0: vy1 = -0.07446872307688, vy2 = +0.7818621047497
    // T_1 = 24.3591072214 / 7 = 3.4798724602, Th = T_1 / 2 = 1.7399362301, 2*alpha = +3*2*pi/7 => alpha = 3*pi/7 = 1.34639685
    Masses m = {1.0L, 1.0L, 0.5L};
    ld x2 = 0.9941458062522L, vy1 = -0.07446872307688L, vy2 = 0.7818621047497L;
    ld Th = 24.3591072214L / 14.0L;
    ld alpha = 3.0L * 3.141592653589793238L / 7.0L;
    bool ok = refine_planar_rpo(x2, vy1, vy2, Th, alpha, m);
    std::cout << std::setprecision(8) << "Planar BHH seed ok=" << ok << " x2=" << (double)x2 << " vy1=" << (double)vy1
              << " vy2=" << (double)vy2 << " T1=" << (double)(2*Th) << " 2a=" << (double)(2*alpha) << "\n\n";
    std::cout << "vy2       x2        vy1       T1        5*T1      2a(rad)   5*2a-4pi  Lam       T1*       kz1        nu_z      5*nu_z\n";
    // Sweep vy2 from 0.78 up to 0.86 (and down to 0.72) to find exact vertical resonance!
    ld bx2=x2, bvy1=vy1, bTh=Th, balpha=alpha;
    for (int step = 0; step <= 40; ++step) {
        ld cvy2 = 0.7818621047497L + step * 0.002L;
        if (!refine_planar_rpo(bx2, bvy1, cvy2, bTh, balpha, m)) break;
        State3D y0 = make_planar(bx2, bvy1, cvy2, m);
        ld T1 = 2.0L * bTh, two_a = 2.0L * balpha;
        // Energy & L
        ld K=0, U=0, Lz=0;
        for(int i=0;i<3;++i) { K+=0.5L*m[i]*(y0[10+3*i]*y0[10+3*i]); Lz+=m[i]*y0[3*i]*y0[10+3*i]; }
        for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) U -= m[i]*m[j]/std::abs(y0[3*j]-y0[3*i]);
        ld E = K + U, Lam = std::abs(Lz)*std::sqrt(-E), T1star = T1*std::pow(-E, 1.5L);
        // Vertical trace at T1:
        const int zidx[6] = {2,5,8, 11,14,17};
        ld tr_z = 0, eps = 1e-7L;
        for (int k=0;k<6;++k) {
            int col = zidx[k];
            State3D yp=y0, ym=y0; yp[col]+=eps; ym[col]-=eps;
            integrate3d(yp, T1, m, 1e-13L); integrate3d(ym, T1, m, 1e-13L);
            tr_z += (yp[col]-ym[col])/(2*eps);
        }
        ld kz1 = tr_z - (2.0L + 2.0L*std::cos(two_a));
        ld nu_z = std::acos(std::clamp(kz1/2.0L, -1.0L, 1.0L)) / (2.0L * 3.141592653589793238L);
        ld five_2a_wrap = 5.0L * two_a - 4.0L * 3.141592653589793238L; // compare with 2a of 3D-Delta14 ~ -3.55..-3.77 (or +2.51..+2.73)!
        // Note: 5 * 2a = 5 * 2.69 ~ 13.46, mod 2pi => +0.90 or -5.38; wait 5*(2a - 2pi) = 5*(-3.59) => mod 2pi
        if (step % 2 == 0) {
            std::cout << std::fixed << std::setprecision(5)
                      << (double)cvy2 << "  " << (double)bx2 << "  " << (double)bvy1 << "  "
                      << (double)T1 << "  " << (double)(5*T1) << "  " << (double)two_a << "  "
                      << (double)(5*two_a - 6*3.141592653589793238L) << "  "
                      << (double)Lam << "  " << (double)T1star << "  " << (double)kz1 << "  "
                      << (double)nu_z << "  " << (double)(5*nu_z) << "\n";
        }
    }
    return 0;
}
