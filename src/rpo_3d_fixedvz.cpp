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
    double cx = u[1]*z - u[2]*y;
    double cy = u[2]*x - u[0]*z;
    double cz = u[0]*y - u[1]*x;
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
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.1);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-12,
                        int max_steps = 200000, double* rmin_out = nullptr, double r_abort = 1e-4) {
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

// Gauge-fixed 3D perp state: x1=-1, vz1=0, free={x2, vy1, vy2, vz2, Th, alpha}
static State3D make_rpo_state(double x2, double vy1, double vy2, double vz2, const Masses& m) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2) / m[2];
    double vy3 = -(m[0]*vy1 + m[1]*vy2) / m[2];
    double vz3 = -(m[1]*vz2) / m[2];
    return {x1, 0, 0, x2, 0, 0, x3, 0, 0,
            0, vy1, 0, 0, vy2, vz2, 0, vy3, vz3};
}

// Compute optimal rotation angle alpha around u=L/|L| that aligns r1(t) onto the x-axis?
// Better: at t=Th, if y(Th) = R_u(alpha) y_perp, then r1(Th) = R_u(alpha) (x1', 0, 0).
// Since u = (0, uy, uz) is perpendicular to x-axis (because Lx = 0 at t=0!),
// R_u(alpha) e_x = (cos a, uz sin a, -uy sin a)!
// So alpha = atan2(r1_y * uz - r1_z * uy, r1_x)!!
static inline double infer_alpha(const State3D& y, const double u[3]) {
    double rx = y[0], ry = y[1], rz = y[2];
    double sin_part = ry * u[2] - rz * u[1];
    double cos_part = rx;
    // Flip sign if rx > 0 so r1 aligns with -x (or +x)
    return std::atan2(sin_part, -cos_part);
}

// Refine 6 vars p = {x2, vy1, vy2, vz2, Th, alpha} (or fix alpha if fix_alpha=true!)
// Residuals F = (yt_1, zt_1, yt_2, zt_2, vxt_1, vxt_2) where yt = R_u(-alpha) y(Th)
static bool refine_rpo(double p[6], const Masses& m, bool fix_alpha, double& ferr) {
    int nvars = fix_alpha ? 5 : 6;
    double lambda = 1e-3;
    auto eval_F = [&](const double q[6], std::array<double, 6>& F) -> bool {
        if (q[4] < 0.5 || q[4] > 30.0) return false;
        State3D y0 = make_rpo_state(q[0], q[1], q[2], q[3], m);
        double d12=std::abs(y0[3]-y0[0]), d13=std::abs(y0[6]-y0[0]), d23=std::abs(y0[6]-y0[3]);
        if (std::min({d12, d13, d23}) < 0.1) return false;
        auto Lv = angmom3d(y0, m);
        double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lmag < 1e-3) return false;
        double u[3] = {Lv[0]/Lmag, Lv[1]/Lmag, Lv[2]/Lmag};

        State3D y = y0;
        double rm = 1e300;
        if (!integrate3d(y, q[4], m, 1e-12, 80000, &rm, 5e-3)) return false;
        State3D yt = rotate_state(y, u, -q[5]);
        F = {yt[1], yt[2], yt[4], yt[5], yt[9], yt[12]};
        return true;
    };

    std::array<double, 6> F;
    if (!eval_F(p, F)) return false;
    double cost = 0; for (double v : F) cost += v*v;

    for (int iter = 0; iter < 35; ++iter) {
        if (std::sqrt(cost) < 1e-10) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[6][6] = {{0}};
        std::array<double, 6> Fp, Fm;
        for (int col = 0; col < nvars; ++col) { if (col == 3) { for (int r=0;r<6;++r) J[r][3]=0; continue; }
            double pp[6], pm[6];
            for (int i=0;i<6;++i) { pp[i]=p[i]; pm[i]=p[i]; }
            pp[col]+=eps; pm[col]-=eps;
            if (!eval_F(pp, Fp) || !eval_F(pm, Fm)) return false;
            for (int r=0;r<6;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double A[6][7] = {{0}};
        for (int r=0;r<nvars;++r) {
            for (int c=0;c<nvars;++c) { double s=0; for(int i=0;i<6;++i) s+=J[i][r]*J[i][c]; A[r][c]=s; }
            A[r][r] += lambda*(1.0 + A[r][r]) + (r==3 ? 1e12 : 0.0);
            double s=0; for(int i=0;i<6;++i) s-=J[i][r]*F[i]; A[r][nvars]=s;
        }
        bool sing = false;
        for (int col=0;col<nvars;++col) {
            int piv=col;
            for (int r=col+1;r<nvars;++r) if (std::abs(A[r][col])>std::abs(A[piv][col])) piv=r;
            if (std::abs(A[piv][col])<1e-18) { sing=true; break; }
            for (int c=col;c<=nvars;++c) std::swap(A[col][c], A[piv][c]);
            for (int r=col+1;r<nvars;++r) {
                double f=A[r][col]/A[col][col];
                for (int c=col;c<=nvars;++c) A[r][c]-=f*A[col][c];
            }
        }
        if (sing) { lambda*=10; continue; }
        double d[6]={0};
        for (int r=nvars-1;r>=0;--r) { double s=A[r][nvars]; for(int c=r+1;c<nvars;++c) s-=A[r][c]*d[c]; d[r]=s/A[r][r]; }
        double sn=0; for(int i=0;i<nvars;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.3) for(int i=0;i<nvars;++i) d[i]*=(0.3/sn);
        double np[6]; for(int i=0;i<6;++i) np[i]=p[i]+d[i];
        std::array<double, 6> nF;
        if (!eval_F(np, nF)) { lambda*=10; continue; }
        double nc=0; for(double v:nF) nc+=v*v;
        if (nc < cost) { for(int i=0;i<6;++i) p[i]=np[i]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return (ferr < 1e-9 && p[4] > 0.5);
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

struct Stab3D { bool stable; double k1, k2, k3, margin; int surv100; double drift100; };

static Stab3D analyze_stab3d_rpo(const State3D& y0, double T, const double u[3], double two_alpha, const Masses& m) {
    double M[18][18], eps = 1e-7;
    for (int col=0;col<18;++col) {
        State3D yp=y0, ym=y0; yp[col]+=eps; ym[col]-=eps;
        integrate3d(yp, T, m, 1e-13, 350000);
        integrate3d(ym, T, m, 1e-13, 350000);
        yp = rotate_state(yp, u, -two_alpha);
        ym = rotate_state(ym, u, -two_alpha);
        for (int r=0;r<18;++r) M[r][col]=(yp[r]-ym[r])/(2*eps);
    }
    double M2[18][18]={{0}};
    for (int i=0;i<18;++i) for (int k=0;k<18;++k) for (int j=0;j<18;++j) M2[i][j]+=M[i][k]*M[k][j];
    double tr1=0, tr2=0, tr3=0;
    for (int i=0;i<18;++i) { tr1+=M[i][i]; tr2+=M2[i][i]; for(int j=0;j<18;++j) tr3+=M2[i][j]*M[j][i]; }
    double c1 = std::cos(two_alpha), c2 = std::cos(2.0 * two_alpha), c3 = std::cos(3.0 * two_alpha);
    double p1 = tr1 - (6.0 + 6.0 * c1);
    double p2 = tr2 - (6.0 + 6.0 * c2) + 6.0;
    double p3 = tr3 - (6.0 + 6.0 * c3) + 3.0 * p1;
    double e1=p1, e2=0.5*(e1*p1-p2), e3=(e2*p1-e1*p2+p3)/3.0;
    double Q=(e1*e1-3.0*e2)/9.0, R=(2.0*e1*e1*e1-9.0*e1*e2+27.0*e3)/54.0, disc=Q*Q*Q-R*R;
    Stab3D st{}; st.margin=-1;
    if (disc >= -1e-4 && Q > 0) {
        double sqQ=std::sqrt(Q), th=std::acos(std::clamp(R/(sqQ*sqQ*sqQ),-1.0,1.0));
        st.k1=2*sqQ*std::cos(th/3)+e1/3; st.k2=2*sqQ*std::cos((th+2*M_PI)/3)+e1/3; st.k3=2*sqQ*std::cos((th+4*M_PI)/3)+e1/3;
        double mk=std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
        st.margin=2.0-mk; st.stable=(mk < 2.0-1e-4);
    }
    if (st.stable) {
        State3D yp = y0;
        double kick = 1e-5;
        yp[9]+=kick; yp[12]-=kick*m[0]/m[1];
        yp[10]+=kick; yp[13]-=kick*m[0]/m[1];
        yp[11]+=kick; yp[14]-=kick*m[0]/m[1];
        double md=0; int surv=0;
        for (int k=0;k<100;++k) {
            if (!integrate3d(yp, T, m, 1e-11, 250000)) break;
            State3D yrot = rotate_state(yp, u, -two_alpha * (k + 1));
            double d=state_dist3d(yrot, y0); if (d>md) md=d; if (d>5.0) break;
            surv=k+1;
        }
        st.surv100=surv; st.drift100=md;
    }
    return st;
}

int main(int argc, char** argv) {
    Masses m = {
        (argc > 1) ? std::atof(argv[1]) : 1.0,
        (argc > 2) ? std::atof(argv[2]) : 1.0,
        (argc > 3) ? std::atof(argv[3]) : 0.5
    };
    int nsamples = (argc > 4) ? std::atoi(argv[4]) : 12000;
    int seed_off = (argc > 5) ? std::atoi(argv[5]) : 0;

    struct Seed { double p[6]; double err; };
    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::vector<Seed>> tseeds(nthreads);
    std::vector<std::thread> ths;

    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            std::mt19937_64 rng(881001 + seed_off*10007 + tid*9973);
            std::uniform_real_distribution<double> u_x2(-0.6, 1.8);
            std::uniform_real_distribution<double> u_vy(-1.2, 1.2);
            std::uniform_real_distribution<double> u_vz(0.2, 1.2);

            int cnt = nsamples / nthreads;
            for (int k = 0; k < cnt; ++k) {
                double x2 = u_x2(rng), vy1 = u_vy(rng), vy2 = u_vy(rng), vz2 = u_vz(rng);
                if (std::abs(vy1) < 0.15) continue;
                State3D y0 = make_rpo_state(x2, vy1, vy2, vz2, m);
                double d12=std::abs(y0[3]-y0[0]), d13=std::abs(y0[6]-y0[0]), d23=std::abs(y0[6]-y0[3]);
                if (std::min({d12, d13, d23}) < 0.25) continue;
                if (energy3d(y0, m) >= -0.05) continue;
                auto Lv = angmom3d(y0, m);
                double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
                if (Lmag < 0.1) continue;
                double u[3] = {Lv[0]/Lmag, Lv[1]/Lmag, Lv[2]/Lmag};

                State3D y = y0;
                double t = 0, h = 1e-3, rmin = 1e300;
                double best_e = 1e9, best_t = 0, best_a = 0, pe = 1e9, p2e = 1e9, pt = 0, pa = 0;
                int steps = 0;
                while (t < 11.0 && steps < 10000) {
                    double hc = std::min(h, 0.03);
                    if (!dopri5_step3d(y, t, hc, 11.0, m, 1e-8, rmin, 1e-2)) break;
                    h = hc; ++steps;
                    double maxr = 0; for (int i=0;i<9;++i) maxr = std::max(maxr, std::abs(y[i]));
                    if (maxr > 18.0) break;

                    double a = infer_alpha(y, u);
                    State3D yt = rotate_state(y, u, -a);
                    double e = std::sqrt(yt[1]*yt[1]+yt[2]*yt[2]+yt[4]*yt[4]+yt[5]*yt[5]+yt[9]*yt[9]+yt[12]*yt[12]);
                    if (t > 0.6 && steps >= 3 && pe < p2e && pe < e && pe < best_e) {
                        best_e = pe; best_t = pt; best_a = pa;
                    }
                    p2e = pe; pe = e; pt = t; pa = a;
                }
                if (best_e < 0.3 && best_t > 0.6) {
                    tseeds[tid].push_back({{x2, vy1, vy2, vz2, best_t, best_a}, best_e});
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Seed> all;
    for (auto& v : tseeds) all.insert(all.end(), v.begin(), v.end());
    std::sort(all.begin(), all.end(), [](const Seed& a, const Seed& b){ return a.err < b.err; });
    size_t limit = std::min(all.size(), (size_t)500);
    std::cerr << "[rpo_3d m=(" << m[0] << "," << m[1] << "," << m[2] << ")] seeds=" << all.size()
              << " refining=" << limit << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 2>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= limit) break;
                double p[6]; for (int i=0;i<6;++i) p[i] = all[idx].p[i];
                double ferr = 0;
                if (!refine_rpo(p, m, false, ferr)) continue;
                if (std::abs(p[3]) < 0.1 || std::abs(p[1]) < 0.1) continue;

                double T = 2.0 * p[4];
                double two_alpha = 2.0 * p[5];
                State3D y0 = make_rpo_state(p[0], p[1], p[2], p[3], m);
                auto Lv = angmom3d(y0, m);
                double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
                if (Lmag < 1e-3) continue;
                double u[3] = {Lv[0]/Lmag, Lv[1]/Lmag, Lv[2]/Lmag};

                State3D yT = y0;
                double rmin = 1e300;
                if (!integrate3d(yT, T, m, 1e-13, 350000, &rmin, 5e-3)) continue;
                State3D yT_rot = rotate_state(yT, u, -two_alpha);
                double ret = state_dist3d(yT_rot, y0);
                if (ret > 1e-7) continue;

                // Check true 3D thickness & partner swaps
                const int N = 320;
                double dt = T / N, C[3][3] = {{0}}, rmax = 0;
                int prev_cl = -1, swaps = 0, cnt = 0;
                State3D yc = y0;
                for (int k = 0; k <= N; ++k) {
                    double r12=std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
                    double r13=std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
                    double r23=std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
                    rmax = std::max({rmax, r12, r13, r23});
                    int cl = (r12<=r13 && r12<=r23)?0:((r13<=r23)?1:2);
                    if (prev_cl != -1 && cl != prev_cl) swaps++;
                    prev_cl = cl;
                    for (int b=0;b<3;++b) {
                        double qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
                        for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += qv[r]*qv[c];
                        cnt++;
                    }
                    if (k < N) integrate3d(yc, dt, m, 1e-11, 50000);
                }
                for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
                double thick = std::sqrt(min_eig_3x3(C));
                if (thick < 0.08) continue;

                double E = energy3d(y0, m);
                double inv1 = E * std::pow(T, 2.0/3.0);
                double inv2 = Lmag * std::sqrt(std::max(0.0, -E));
                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup)
                        if (std::abs(f[0]-inv1) < 2e-2 && std::abs(f[1]-inv2) < 2e-2) { dup = true; break; }
                    if (dup) continue;
                    dedup.push_back({inv1, inv2});
                }

                auto st = analyze_stab3d_rpo(y0, T, u, two_alpha, m);
                if (!st.stable || st.surv100 < 100) continue;

                std::ostringstream os;
                os << std::setprecision(13);
                os << "{\"dim\":3,\"subspace\":\"3d_rpo_rotating\",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "],\"y0\":[";
                for (int i=0;i<18;++i) os << y0[i] << (i<17?",":"");
                os << "],\"T\":" << T << ",\"precess_rad\":" << two_alpha << ",\"E\":" << E << ",\"L\":" << Lmag
                   << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]"
                   << ",\"ret_err\":" << ret << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2 << ",\"k3\":" << st.k3
                   << ",\"stab_margin\":" << st.margin << ",\"survived_100T\":" << st.surv100
                   << ",\"drift_100T\":" << st.drift100 << ",\"rmin\":" << rmin << ",\"rmax\":" << rmax
                   << ",\"thickness_3d\":" << thick << ",\"partner_swaps\":" << swaps
                   << ",\"democratic\":" << (swaps >= 2 ? "true" : "false") << "}\n";
                std::lock_guard<std::mutex> lk(mtx);
                std::cout << os.str() << std::flush;
            }
        });
    }
    for (auto& th : ths) th.join();
    return 0;
}
