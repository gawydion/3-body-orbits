#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
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
    for (int b=0;b<6;++b) { ld x=s[3*b], y=s[3*b+1]; s[3*b]=x*c-y*sn; s[3*b+1]=x*sn+y*c; }
}
static bool integrate3d(State3D& y, ld T, const Masses& m, ld tol=1e-13L) {
    static constexpr ld a21=1.0L/5, a31=3.0L/40, a32=9.0L/40,
        a41=44.0L/45, a42=-56.0L/15, a43=32.0L/9,
        a51=19372.0L/6561, a52=-25360.0L/2187, a53=64448.0L/6561, a54=-212.0L/729,
        a61=9017.0L/3168, a62=-355.0L/33, a63=46732.0L/5247, a64=49.0L/176, a65=-5103.0L/18656,
        b1=35.0L/384, b3=500.0L/1113, b4=125.0L/192, b5=-2187.0L/6784, b6=11.0L/84,
        e1=71.0L/57600, e3=-71.0L/16695, e4=71.0L/1920, e5=-17253.0L/339200, e6=22.0L/525, e7=-1.0L/40;
    ld t=0, h=std::min(1e-3L, T*0.01L); int steps=0; State3D k1,k2,k3,k4,k5,k6,k7,yt,yn;
    while (t < T - 1e-16L && steps < 400000) {
        if (t+h>T) h=T-t; bool ok=false;
        for (int att=0;att<25;++att) {
            deriv3d(y,m,k1);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv3d(yt,m,k2);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv3d(yt,m,k3);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv3d(yt,m,k4);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv3d(yt,m,k5);
            for(int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv3d(yt,m,k6);
            for(int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); deriv3d(yn,m,k7);
            ld err=0; for(int i=0;i<18;++i) { ld ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]); ld sc=tol*(1+std::max(std::abs(y[i]),std::abs(yn[i]))); ld r=ei/sc; err+=r*r; }
            err=std::sqrt(err/18.0L);
            if(err<=1){y=yn;t+=h;h=std::min(h*((err==0)?5.0L:std::min(5.0L,std::max(0.2L,0.9L*std::pow(err,-0.2L)))),0.05L);ok=true;break;}
            h*=std::max(0.1L,0.9L*std::pow(err,-0.25L));
        }
        if(!ok) return false; ++steps;
    }
    return true;
}
int main() {
    // At planar BHH 5:2 vertical resonance (vy2=0.8156): refine planar RPO, compute ALL 3 Broucke indices (k1, k2, kz) for T1 AND 5*T1!
    Masses m = {1.0L, 1.0L, 0.5L};
    // From bhh_vertical_tune output: vy2=0.8156, x2=0.9540, vy1=-0.1002, T1=3.669, 2a=3.0715
    ld x2=0.95395L, vy1=-0.10020L, vy2=0.81560L, Th=3.6692L/2.0L, alpha=3.0715L/2.0L;
    auto make_pl = [&](ld q0, ld q1) -> State3D {
        ld x3 = -(m[0]*(-1.0L)+m[1]*q0)/m[2], vy3 = -(m[0]*q1+m[1]*vy2)/m[2];
        return {-1,0,0, q0,0,0, x3,0,0, 0,q1,0, 0,vy2,0, 0,vy3,0};
    };
    ld p[4]={x2,vy1,Th,alpha};
    auto eval = [&](const ld q[4], ld F[4]) {
        State3D y = make_pl(q[0], q[1]); integrate3d(y, q[2], m, 1e-14L); rot_z(y, -q[3]);
        F[0]=y[1]; F[1]=y[4]; F[2]=y[9]; F[3]=y[12];
    };
    ld F[4]; eval(p, F);
    for (int it=0;it<15;++it) {
        ld eps=1e-7L, J[4][4], Fp[4], Fm[4];
        for (int c=0;c<4;++c) { ld pp[4]={p[0],p[1],p[2],p[3]}, pm[4]={p[0],p[1],p[2],p[3]}; pp[c]+=eps; pm[c]-=eps; eval(pp,Fp); eval(pm,Fm); for(int r=0;r<4;++r) J[r][c]=(Fp[r]-Fm[r])/(2*eps); }
        ld A[4][5]; for(int r=0;r<4;++r){for(int c=0;c<4;++c) A[r][c]=J[r][c]; A[r][4]=-F[r];}
        for(int c=0;c<4;++c){int piv=c;for(int r=c+1;r<4;++r)if(std::abs(A[r][c])>std::abs(A[piv][c]))piv=r; for(int k=c;k<5;++k)std::swap(A[c][k],A[piv][k]); for(int r=c+1;r<4;++r){ld f=A[r][c]/A[c][c]; for(int k=c;k<5;++k)A[r][k]-=f*A[c][k];}}
        ld d[4]={0}; for(int r=3;r>=0;--r){ld s=A[r][4]; for(int c=r+1;c<4;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r];}
        for(int i=0;i<4;++i) p[i]+=d[i]; eval(p,F);
    }
    State3D y0 = make_pl(p[0], p[1]); ld T1 = 2*p[2], two_a = 2*p[3];
    // Compute 12x12 planar RPO monodromy + 6x6 vertical monodromy at T1:
    const int pidx[12]={0,1,3,4,6,7, 9,10,12,13,15,16};
    const int zidx[6]={2,5,8, 11,14,17};
    ld Mp[12][12], Mz[6][6], eps=1e-7L;
    for(int c=0;c<12;++c) {
        State3D yp=y0, ym=y0; yp[pidx[c]]+=eps; ym[pidx[c]]-=eps;
        integrate3d(yp, T1, m, 1e-14L); integrate3d(ym, T1, m, 1e-14L);
        rot_z(yp, -two_a); rot_z(ym, -two_a);
        for(int r=0;r<12;++r) Mp[r][c]=(yp[pidx[r]]-ym[pidx[r]])/(2*eps);
    }
    ld tr1=0, tr2=0;
    for(int i=0;i<12;++i) { tr1+=Mp[i][i]; for(int j=0;j<12;++j) tr2+=Mp[i][j]*Mp[j][i]; }
    // Trivial eigenvalues of 12x12 planar RPO (rotation two_a):
    // COM: 4 eigenvalues e^{+-i 2a} (2 pairs) -> sum = 4 cos(2a), sum sq = 4 cos(4a)
    // Flow/E + Rot/L: 4 eigenvalues +1 -> sum = 4, sum sq = 4
    // 2 nontrivial pairs: k1, k2 -> sum = k1+k2, sum sq = k1^2 + k2^2 - 4
    ld s = tr1 - (4.0L + 4.0L*std::cos(two_a));
    ld p2 = tr2 - (4.0L + 4.0L*std::cos(2*two_a)) + 4.0L;
    ld prod = 0.5L*(s*s - p2);
    ld disc = s*s - 4*prod;
    ld k1 = 0.5L*(s + std::sqrt(std::max(0.0L, disc))), k2 = 0.5L*(s - std::sqrt(std::max(0.0L, disc)));
    ld tr_z=0;
    for(int c=0;c<6;++c) { State3D yp=y0, ym=y0; yp[zidx[c]]+=eps; ym[zidx[c]]-=eps; integrate3d(yp,T1,m,1e-14L); integrate3d(ym,T1,m,1e-14L); tr_z+=(yp[zidx[c]]-ym[zidx[c]])/(2*eps); }
    ld kz = tr_z - (2.0L + 2.0L*std::cos(two_a));
    auto tune = [](ld k) { return std::acos(std::clamp(k/2.0L,-1.0L,1.0L))/(2*3.141592653589793238L); };
    auto k5 = [&](ld k) { return 2.0L*std::cos(5.0L*std::acos(std::clamp(k/2.0L,-1.0L,1.0L))); };
    std::cout << std::setprecision(8);
    std::cout << "Planar BHH at vy2=0.8156: T1=" << (double)T1 << " 5*T1=" << (double)(5*T1) << " 5*2a-6pi=" << (double)(5*two_a-6*3.141592653589793238L) << "\n";
    std::cout << "  1-period: k1=" << (double)k1 << " (nu=" << (double)tune(k1) << "), k2=" << (double)k2 << " (nu=" << (double)tune(k2) << "), kz=" << (double)kz << " (nu=" << (double)tune(kz) << ")\n";
    std::cout << "  5-period: k1^(5)=" << (double)k5(k1) << ", k2^(5)=" << (double)k5(k2) << ", kz^(5)=" << (double)k5(kz) << "\n";
    return 0;
}
