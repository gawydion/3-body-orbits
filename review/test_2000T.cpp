#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using State = std::array<double, 12>;
using Masses = std::array<double, 3>;

static inline void deriv(const State& s, const Masses& m, State& ds, double* rm = nullptr) {
    for (int i=0;i<6;++i) ds[i] = s[6+i];
    double a[6]={0}, r2m=1e300;
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
        double dx=s[2*j]-s[2*i], dy=s[2*j+1]-s[2*i+1];
        double r2=dx*dx+dy*dy; if(r2<r2m) r2m=r2;
        double ir=1.0/std::sqrt(r2), ir3=ir*ir*ir;
        a[2*i]+=m[j]*dx*ir3; a[2*i+1]+=m[j]*dy*ir3;
        a[2*j]-=m[i]*dx*ir3; a[2*j+1]-=m[i]*dy*ir3;
    }
    for (int i=0;i<6;++i) ds[6+i]=a[i];
    if(rm) *rm = std::sqrt(r2m);
}
static bool integrate(State& y, double T, const Masses& m, double tol=1e-12) {
    static constexpr double a21=1.0/5, a31=3.0/40, a32=9.0/40,
        a41=44.0/45, a42=-56.0/15, a43=32.0/9,
        a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729,
        a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656,
        b1=35.0/384, b3=500.0/1113, b4=125.0/192, b5=-2187.0/6784, b6=11.0/84,
        e1=71.0/57600, e3=-71.0/16695, e4=71.0/1920, e5=-17253.0/339200, e6=22.0/525, e7=-1.0/40;
    double t=0, h=std::min(1e-3, T*0.01); int steps=0; State k1,k2,k3,k4,k5,k6,k7,yt,yn;
    while (t < T - 1e-15 && steps < 400000) {
        if(t+h>T) h=T-t; bool ok=false;
        for(int att=0;att<25;++att) {
            double rm=1e300; deriv(y,m,k1,&rm); if(rm<1e-6) return false;
            for(int i=0;i<12;++i) yt[i]=y[i]+h*a21*k1[i]; deriv(yt,m,k2);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv(yt,m,k3);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(yt,m,k4);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(yt,m,k5);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(yt,m,k6);
            for(int i=0;i<12;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
            deriv(yn,m,k7,&rm); if(rm<1e-6) return false;
            double err=0; for(int i=0;i<12;++i){double ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]); double sc=tol*(1+std::max(std::abs(y[i]),std::abs(yn[i]))); double r=ei/sc; err+=r*r;}
            err=std::sqrt(err/12.0);
            if(err<=1){y=yn;t+=h;h=std::min(h*((err==0)?5.0:std::min(5.0,std::max(0.2,0.9*std::pow(err,-0.2)))),0.05);ok=true;break;}
            h*=std::max(0.1,0.9*std::pow(err,-0.25)); if(h<1e-13) return false;
        }
        if(!ok) return false; ++steps;
        for(int i=0;i<6;++i) if(std::abs(y[i])>50.0) return false;
    }
    return t >= T - 1e-13;
}
struct Job { std::string id; double m3, vx, vy, T; };
int main() {
    std::vector<Job> jobs; Job j;
    while (std::cin >> j.id >> j.m3 >> j.vx >> j.vy >> j.T) jobs.push_back(j);
    std::vector<std::string> outs(jobs.size());
    std::vector<std::thread> ths;
    for(size_t idx=0; idx<jobs.size(); ++idx) ths.emplace_back([&,idx]() {
        auto jb = jobs[idx];
        Masses m = {1.0, 1.0, jb.m3};
        State y0 = {-1,0, 1,0, 0,0, jb.vx,jb.vy, jb.vx,jb.vy, -2*jb.vx/jb.m3, -2*jb.vy/jb.m3};
        State yp = y0;
        // Same 1e-5 kick as campaign_2d.cpp (plus test an orthogonal vy kick!)
        yp[6] += 1e-5; yp[8] -= 1e-5;
        int surv = 0; double md100=0, md500=0, md2000=0;
        for (int k=1; k<=2000; ++k) {
            if (!integrate(yp, jb.T, m, 1e-12)) break;
            double d=0; for(int i=0;i<12;++i) d+=(yp[i]-y0[i])*(yp[i]-y0[i]); d=std::sqrt(d);
            if (k<=100) md100=std::max(md100, d);
            if (k<=500) md500=std::max(md500, d);
            md2000=std::max(md2000, d);
            if (d > 5.0) break;
            surv = k;
        }
        // Also test small kick 1e-4 in vy1/vy2 (symmetric break)
        State yp2 = y0; yp2[7] += 1e-5; yp2[9] -= 1e-5;
        int surv2 = 0; double md2 = 0;
        for (int k=1; k<=2000; ++k) {
            if (!integrate(yp2, jb.T, m, 1e-12)) break;
            double d=0; for(int i=0;i<12;++i) d+=(yp2[i]-y0[i])*(yp2[i]-y0[i]); d=std::sqrt(d);
            md2=std::max(md2, d); if (d > 5.0) break; surv2 = k;
        }
        std::ostringstream os;
        os << jb.id << " kick_vx: surv=" << surv << "/2000 (d100=" << std::setprecision(4) << md100
           << " d500=" << md500 << " d2000=" << md2000 << ") | kick_vy: surv=" << surv2 << "/2000 (d=" << md2 << ")";
        outs[idx] = os.str();
    });
    for(auto& th : ths) th.join();
    for(auto& s : outs) std::cout << s << "\n";
    return 0;
}
