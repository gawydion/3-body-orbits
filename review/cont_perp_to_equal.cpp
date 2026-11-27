#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using ld = long double;
using State = std::array<ld, 12>;
using Masses = std::array<ld, 3>;

static inline void deriv(const State& s, const Masses& m, State& ds) {
    for (int i=0;i<6;++i) ds[i] = s[6+i];
    ld a[6]={0};
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
        ld dx=s[2*j]-s[2*i], dy=s[2*j+1]-s[2*i+1];
        ld r2=dx*dx+dy*dy, ir=1.0L/std::sqrt(r2), ir3=ir*ir*ir;
        a[2*i]+=m[j]*dx*ir3; a[2*i+1]+=m[j]*dy*ir3;
        a[2*j]-=m[i]*dx*ir3; a[2*j+1]-=m[i]*dy*ir3;
    }
    for (int i=0;i<6;++i) ds[6+i]=a[i];
}

static bool integrate(State& y, ld T, const Masses& m, ld tol=1e-13L) {
    static constexpr ld a21=1.0L/5, a31=3.0L/40, a32=9.0L/40,
        a41=44.0L/45, a42=-56.0L/15, a43=32.0L/9,
        a51=19372.0L/6561, a52=-25360.0L/2187, a53=64448.0L/6561, a54=-212.0L/729,
        a61=9017.0L/3168, a62=-355.0L/33, a63=46732.0L/5247, a64=49.0L/176, a65=-5103.0L/18656,
        b1=35.0L/384, b3=500.0L/1113, b4=125.0L/192, b5=-2187.0L/6784, b6=11.0L/84,
        e1=71.0L/57600, e3=-71.0L/16695, e4=71.0L/1920, e5=-17253.0L/339200, e6=22.0L/525, e7=-1.0L/40;
    ld t=0, h=std::min(1e-3L, T*0.01L); int steps=0; State k1,k2,k3,k4,k5,k6,k7,yt,yn;
    while (t < T - 1e-16L && steps < 300000) {
        if (t+h>T) h=T-t; bool ok=false;
        for (int att=0;att<25;++att) {
            deriv(y,m,k1);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*a21*k1[i]; deriv(yt,m,k2);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv(yt,m,k3);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(yt,m,k4);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(yt,m,k5);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(yt,m,k6);
            for(int i=0;i<12;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); deriv(yn,m,k7);
            ld err=0; for(int i=0;i<12;++i) { ld ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]); ld sc=tol*(1+std::max(std::abs(y[i]),std::abs(yn[i]))); ld r=ei/sc; err+=r*r; }
            err=std::sqrt(err/12.0L);
            if(err<=1){y=yn;t+=h;h=std::min(h*((err==0)?5.0L:std::min(5.0L,std::max(0.2L,0.9L*std::pow(err,-0.2L)))),0.05L);ok=true;break;}
            h*=std::max(0.1L,0.9L*std::pow(err,-0.25L));
            if(h<1e-14L) return false;
        }
        if(!ok) return false; ++steps;
    }
    return t >= T - 1e-13L;
}

static State make_perp(ld x2, ld v1, ld v2, const Masses& m) {
    ld x1=-1.0L, x3=-(m[0]*x1+m[1]*x2)/m[2], v3=-(m[0]*v1+m[1]*v2)/m[2];
    return {x1,0, x2,0, x3,0, 0,v1, 0,v2, 0,v3};
}

static bool refine(ld p[4], const Masses& m) {
    ld lam = 1e-5L;
    auto eval = [&](const ld q[4], ld F[4]) -> bool {
        State y = make_perp(q[0], q[1], q[2], m);
        if (!integrate(y, q[3], m, 1e-12L)) return false;
        F[0]=y[1]; F[1]=y[3]; F[2]=y[6]; F[3]=y[8];
        return true;
    };
    ld F[4]; if (!eval(p, F)) return false;
    ld cost = F[0]*F[0]+F[1]*F[1]+F[2]*F[2]+F[3]*F[3];
    for (int it=0;it<25;++it) {
        if (std::sqrt(cost) < 1e-11L) break;
        ld eps=1e-7L, J[4][4], Fp[4], Fm[4];
        for (int c=0;c<4;++c) {
            ld pp[4]={p[0],p[1],p[2],p[3]}, pm[4]={p[0],p[1],p[2],p[3]};
            pp[c]+=eps; pm[c]-=eps;
            if (!eval(pp,Fp) || !eval(pm,Fm)) return false;
            for (int r=0;r<4;++r) J[r][c]=(Fp[r]-Fm[r])/(2*eps);
        }
        ld A[4][5]={{0}};
        for (int r=0;r<4;++r) {
            for (int c=0;c<4;++c) for(int i=0;i<4;++i) A[r][c]+=J[i][r]*J[i][c];
            A[r][r]+=lam*(1+A[r][r]);
            for (int i=0;i<4;++i) A[r][4]-=J[i][r]*F[i];
        }
        for (int c=0;c<4;++c) {
            int piv=c; for(int r=c+1;r<4;++r) if(std::abs(A[r][c])>std::abs(A[piv][c])) piv=r;
            if (std::abs(A[piv][c])<1e-18L) return false;
            for(int k=c;k<5;++k) std::swap(A[c][k], A[piv][k]);
            for(int r=c+1;r<4;++r){ld f=A[r][c]/A[c][c]; for(int k=c;k<5;++k) A[r][k]-=f*A[c][k];}
        }
        ld d[4]={0}; for(int r=3;r>=0;--r){ld s=A[r][4]; for(int c=r+1;c<4;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r];}
        ld sn=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
        if(sn>0.2L) for(int i=0;i<4;++i) d[i]*=(0.2L/sn);
        ld np[4]={p[0]+d[0],p[1]+d[1],p[2]+d[2],p[3]+d[3]}, nF[4];
        if (!eval(np, nF)) { lam*=10; continue; }
        ld nc=nF[0]*nF[0]+nF[1]*nF[1]+nF[2]*nF[2]+nF[3]*nF[3];
        if (nc < cost) { for(int i=0;i<4;++i) p[i]=np[i]; for(int i=0;i<4;++i) F[i]=nF[i]; cost=nc; lam=std::max(1e-13L, lam*0.2L); }
        else { lam*=5; if(lam>1e7L) break; }
    }
    return std::sqrt(cost) < 1e-9L;
}

static void trace_to_equal(const char* name, Masses m0, ld p[4], int steps) {
    std::cout << "=== Continuing " << name << " from m=(" << (double)m0[0] << "," << (double)m0[1] << "," << (double)m0[2] << ") to (1,1,1) ===\n";
    for (int s=0; s<=steps; ++s) {
        ld frac = (ld)s / steps;
        Masses m = {m0[0] + frac*(1.0L-m0[0]), m0[1] + frac*(1.0L-m0[1]), m0[2] + frac*(1.0L-m0[2])};
        if (!refine(p, m)) { std::cout << "  Stopped at s=" << s << "\n"; break; }
        if (s == 0 || s == steps || s % (steps/4) == 0) {
            State y0 = make_perp(p[0], p[1], p[2], m);
            ld K=0, U=0, L=0;
            for(int i=0;i<3;++i){K+=0.5L*m[i]*y0[7+2*i]*y0[7+2*i]; L+=m[i]*y0[2*i]*y0[7+2*i];}
            for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) U-=m[i]*m[j]/std::abs(y0[2*j]-y0[2*i]);
            ld E=K+U, T=2*p[3], Ts=T*std::pow(-E,1.5L), Lm=std::abs(L)*std::sqrt(-E);
            std::cout << std::setprecision(8) << "  s=" << s << " m=(" << (double)m[0] << "," << (double)m[1] << "," << (double)m[2]
                      << ") x2=" << (double)p[0] << " v1=" << (double)p[1] << " v2=" << (double)p[2]
                      << " T=" << (double)T << " E=" << (double)E << " L=" << (double)L
                      << " T*=" << (double)Ts << " Lam=" << (double)Lm << "\n";
        }
    }
}

int main() {
    // 2D-DEM-006: m=[1, 0.9, 0.8], x2=-0.05951905193104, v1=0.7812388652309, v2=-0.1613813787143, T=15.49514153178
    ld p6[4] = {-0.05951905193104L, 0.7812388652309L, -0.1613813787143L, 15.49514153178L/2.0L};
    trace_to_equal("2D-DEM-006", {1.0L, 0.9L, 0.8L}, p6, 20);
    // 2D-DEM-007: m=[1, 0.6, 0.3], x2=-0.2924901591259, v1=-0.7862827491696, v2=1.190621532647, T=16.88230020881
    ld p7[4] = {-0.2924901591259L, -0.7862827491696L, 1.190621532647L, 16.88230020881L/2.0L};
    trace_to_equal("2D-DEM-007", {1.0L, 0.6L, 0.3L}, p7, 20);
    // 2D-DEM-029: m=[1, 1, 0.5]
    return 0;
}
