#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using State3D = std::array<double, 18>;
using Masses = std::array<double, 3>;

static inline void deriv3d(const State3D& s, const Masses& m, State3D& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 9; ++i) ds[i] = s[9 + i];
    double a[9] = {0};
    double rmin2 = 1e300;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j] - s[3*i], dy = s[3*j+1] - s[3*i+1], dz = s[3*j+2] - s[3*i+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]   += m[j]*dx*ir3; a[3*i+1] += m[j]*dy*ir3; a[3*i+2] += m[j]*dz*ir3;
            a[3*j]   -= m[i]*dx*ir3; a[3*j+1] -= m[i]*dy*ir3; a[3*j+2] -= m[i]*dz*ir3;
        }
    }
    for (int i = 0; i < 9; ++i) ds[9 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double energy3d(const State3D& s, const Masses& m) {
    double K = 0.0, U = 0.0;
    for (int i = 0; i < 3; ++i)
        K += 0.5 * m[i] * (s[9+3*i]*s[9+3*i] + s[10+3*i]*s[10+3*i] + s[11+3*i]*s[11+3*i]);
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            U -= m[i]*m[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
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
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate3d(State3D& y, double T, const Masses& m, double tol = 1e-12,
                        int max_steps = 120000, double* rmin_out = nullptr, double r_abort = 1e-4) {
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

// 3D C2-symmetric Sitnikov/Helical initial state:
// r1 = (1, 0, 0), r2 = (-1, 0, 0), r3 = (0, 0, 0)
// v1 = (0, vy, vz), v2 = (0, -vy, vz), v3 = (0, 0, -2*vz/m3)
static State3D make_sitnikov3d(double vy, double vz, double m3) {
    return {1, 0, 0,  -1, 0, 0,  0, 0, 0,
            0, vy, vz,  0, -vy, vz,  0, 0, -2.0 * vz / m3};
}

// Half-period condition type 0 (cross x-axis perpendicularly): y1=0, z1=0, vx1=0
// Half-period condition type 1 (turning point in z at x-axis): y1=0, vz1=0, vx1=0
static bool refine_sitnikov(double& vy, double& vz, double& Th, double m3, int cond_type, double& ferr) {
    Masses m = {1.0, 1.0, m3};
    double lambda = 1e-4;
    auto eval_F = [&](double cvy, double cvz, double cTh, std::array<double, 3>& F) -> bool {
        if (cTh < 0.5 || cTh > 25.0 || std::abs(cvy) < 0.05 || std::abs(cvz) < 0.05) return false;
        State3D y = make_sitnikov3d(cvy, cvz, m3);
        double rm = 1e300;
        if (!integrate3d(y, cTh, m, 1e-12, 80000, &rm, 5e-3)) return false;
        if (cond_type == 0) F = {y[1], y[2], y[9]};   // y1, z1, vx1
        else                F = {y[1], y[11], y[9]};  // y1, vz1, vx1
        return true;
    };
    std::array<double, 3> F;
    if (!eval_F(vy, vz, Th, F)) return false;
    double cost = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];

    for (int it = 0; it < 30; ++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost) * 1e-2));
        double J[3][3];
        std::array<double, 3> Fp, Fm;
        double p[3] = {vy, vz, Th};
        for (int col = 0; col < 3; ++col) {
            double pp[3] = {p[0], p[1], p[2]}, pm[3] = {p[0], p[1], p[2]};
            pp[col] += eps; pm[col] -= eps;
            if (!eval_F(pp[0], pp[1], pp[2], Fp) || !eval_F(pm[0], pm[1], pm[2], Fm)) return false;
            for (int r = 0; r < 3; ++r) J[r][col] = (Fp[r] - Fm[r]) / (2.0 * eps);
        }
        double A[3][4] = {{0}};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                double s = 0; for (int i = 0; i < 3; ++i) s += J[i][r]*J[i][c];
                A[r][c] = s;
            }
            A[r][r] += lambda*(1.0 + A[r][r]);
            double s = 0; for (int i = 0; i < 3; ++i) s -= J[i][r]*F[i];
            A[r][3] = s;
        }
        bool sing = false;
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col+1; r < 3; ++r) if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
            if (std::abs(A[piv][col]) < 1e-18) { sing = true; break; }
            for (int c = col; c < 4; ++c) std::swap(A[col][c], A[piv][c]);
            for (int r = col+1; r < 3; ++r) {
                double f = A[r][col]/A[col][col];
                for (int c = col; c < 4; ++c) A[r][c] -= f*A[col][c];
            }
        }
        if (sing) { lambda *= 10.0; continue; }
        double d[3] = {0};
        for (int r = 2; r >= 0; --r) {
            double s = A[r][3]; for (int c = r+1; c < 3; ++c) s -= A[r][c]*d[c];
            d[r] = s / A[r][r];
        }
        double sn = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        if (sn > 0.3) for (int i = 0; i < 3; ++i) d[i] *= (0.3 / sn);
        std::array<double, 3> nF;
        if (!eval_F(vy+d[0], vz+d[1], Th+d[2], nF)) { lambda *= 10.0; continue; }
        double nc = nF[0]*nF[0] + nF[1]*nF[1] + nF[2]*nF[2];
        if (nc < cost) { vy+=d[0]; vz+=d[1]; Th+=d[2]; F=nF; cost=nc; lambda=std::max(1e-12, lambda*0.25); }
        else { lambda *= 5.0; if (lambda > 1e7) break; }
    }
    ferr = std::sqrt(cost);
    return (ferr < 1e-10 && Th > 0.5);
}

static double min_eig_3x3(double C[3][3]) {
    double I1 = C[0][0] + C[1][1] + C[2][2];
    double I2 = C[0][0]*C[1][1] + C[0][0]*C[2][2] + C[1][1]*C[2][2]
              - C[0][1]*C[0][1] - C[0][2]*C[0][2] - C[1][2]*C[1][2];
    double I3 = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[1][2])
              - C[0][1]*(C[0][1]*C[2][2] - C[1][2]*C[0][2])
              + C[0][2]*(C[0][1]*C[1][2] - C[1][1]*C[0][2]);
    double p = I2 - I1*I1/3.0;
    double q = -2.0*I1*I1*I1/27.0 + I1*I2/3.0 - I3;
    double r = std::sqrt(std::max(0.0, -p/3.0));
    if (r < 1e-15) return std::max(0.0, I1/3.0);
    double c = std::clamp(-q / (2.0*r*r*r), -1.0, 1.0);
    double phi = std::acos(c) / 3.0;
    double l1 = I1/3.0 + 2.0*r*std::cos(phi);
    double l2 = I1/3.0 + 2.0*r*std::cos(phi + 2.0*M_PI/3.0);
    double l3 = I1/3.0 + 2.0*r*std::cos(phi + 4.0*M_PI/3.0);
    return std::max(0.0, std::min({l1, l2, l3}));
}

struct Stab3D { bool stable; double k1, k2, k3, margin; int surv100; double drift100; };

static Stab3D analyze_stab3d(const State3D& y0, double T, const Masses& m) {
    double M[18][18], eps = 1e-7;
    for (int col = 0; col < 18; ++col) {
        State3D yp = y0, ym = y0; yp[col] += eps; ym[col] -= eps;
        integrate3d(yp, T, m, 1e-13, 300000);
        integrate3d(ym, T, m, 1e-13, 300000);
        for (int r = 0; r < 18; ++r) M[r][col] = (yp[r] - ym[r]) / (2.0 * eps);
    }
    double M2[18][18] = {{0}};
    for (int i = 0; i < 18; ++i)
        for (int k = 0; k < 18; ++k)
            for (int j = 0; j < 18; ++j)
                M2[i][j] += M[i][k] * M[k][j];
    double tr1 = 0, tr2 = 0, tr3 = 0;
    for (int i = 0; i < 18; ++i) {
        tr1 += M[i][i]; tr2 += M2[i][i];
        for (int j = 0; j < 18; ++j) tr3 += M2[i][j] * M[j][i];
    }
    double p1 = tr1 - 12.0;
    double p2 = tr2 - 6.0;
    double p3 = tr3 - 12.0 + 3.0 * p1;
    double e1 = p1, e2 = 0.5*(e1*p1 - p2), e3 = (e2*p1 - e1*p2 + p3)/3.0;
    double Q = (e1*e1 - 3.0*e2)/9.0;
    double R = (2.0*e1*e1*e1 - 9.0*e1*e2 + 27.0*e3)/54.0;
    double disc = Q*Q*Q - R*R;

    Stab3D st{}; st.margin = -1.0;
    if (disc >= -1e-4 && Q > 0.0) {
        double sqQ = std::sqrt(Q);
        double th = std::acos(std::clamp(R / (sqQ*sqQ*sqQ), -1.0, 1.0));
        st.k1 = 2.0*sqQ*std::cos(th/3.0) + e1/3.0;
        st.k2 = 2.0*sqQ*std::cos((th + 2.0*M_PI)/3.0) + e1/3.0;
        st.k3 = 2.0*sqQ*std::cos((th + 4.0*M_PI)/3.0) + e1/3.0;
        double mk = std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
        st.margin = 2.0 - mk;
        st.stable = (mk < 2.0 - 1e-4);
    }

    if (st.stable) {
        // Symmetry-breaking 3D kick in ALL 3 axes!
        State3D yp = y0;
        double kick = 1e-5;
        yp[9]  += kick; yp[12] -= kick * m[0] / m[1];
        yp[10] += kick; yp[13] -= kick * m[0] / m[1];
        yp[11] += kick; yp[17] -= kick * m[0] / m[2];
        double md = 0; int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate3d(yp, T, m, 1e-11, 200000)) break;
            double d = state_dist3d(yp, y0);
            if (d > md) md = d;
            if (d > 5.0) break;
            surv = k + 1;
        }
        st.surv100 = surv; st.drift100 = md;
    }
    return st;
}

int main(int argc, char** argv) {
    double m3 = (argc > 1) ? std::atof(argv[1]) : 0.5;
    int gridN = (argc > 2) ? std::atoi(argv[2]) : 130;
    double Thmax = (argc > 3) ? std::atof(argv[3]) : 14.0;
    Masses m = {1.0, 1.0, m3};

    double vcirc = std::sqrt(m3 + 0.25);
    double vymin = 0.45 * vcirc, vymax = 1.35 * vcirc;
    double vzmin = 0.08 * vcirc, vzmax = 0.85 * vcirc;
    struct Cell { double vy, vz, Th, err; int type; };
    std::vector<std::vector<Cell>> g0(gridN, std::vector<Cell>(gridN));
    std::vector<std::vector<Cell>> g1(gridN, std::vector<Cell>(gridN));

    int nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> ths;
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&, tid]() {
            for (int ix = tid; ix < gridN; ix += nthreads) {
                double vy = vymin + (vymax - vymin) * ix / (gridN - 1);
                for (int iy = 0; iy < gridN; ++iy) {
                    double vz = vzmin + (vzmax - vzmin) * iy / (gridN - 1);
                    State3D y0 = make_sitnikov3d(vy, vz, m3);
                    if (energy3d(y0, m) >= -0.05) { g0[ix][iy]={vy,vz,0,1e9,0}; g1[ix][iy]={vy,vz,0,1e9,1}; continue; }
                    State3D y = y0;
                    double t = 0, h = 1e-3, rmin = 1e300;
                    double b0=1e9, t0=0, p0=1e9, pp0=1e9;
                    double b1=1e9, t1=0, p1=1e9, pp1=1e9;
                    double pt = 0;
                    int steps = 0;
                    while (t < Thmax && steps < 12000) {
                        double hc = std::min(h, 0.03);
                        if (!dopri5_step3d(y, t, hc, Thmax, m, 1e-8, rmin, 1e-2)) break;
                        h = hc; ++steps;
                        double maxr = 0; for (int i=0;i<9;++i) maxr = std::max(maxr, std::abs(y[i]));
                        if (maxr > 20.0) break;
                        double e0 = std::sqrt(y[1]*y[1] + y[2]*y[2] + y[9]*y[9]);
                        double e1 = std::sqrt(y[1]*y[1] + y[11]*y[11] + y[9]*y[9]);
                        if (t > 0.6 && steps >= 3) {
                            if (p0 < pp0 && p0 < e0 && p0 < b0) { b0 = p0; t0 = pt; }
                            if (p1 < pp1 && p1 < e1 && p1 < b1) { b1 = p1; t1 = pt; }
                        }
                        pp0 = p0; p0 = e0; pp1 = p1; p1 = e1; pt = t;
                    }
                    g0[ix][iy] = {vy, vz, t0, b0, 0};
                    g1[ix][iy] = {vy, vz, t1, b1, 1};
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::vector<Cell> seeds;
    auto harvest = [&](const std::vector<std::vector<Cell>>& g) {
        for (int ix = 1; ix < gridN - 1; ++ix)
            for (int iy = 1; iy < gridN - 1; ++iy) {
                double e = g[ix][iy].err;
                if (e > 0.25 || g[ix][iy].Th < 0.6) continue;
                bool ok = true;
                for (int dx=-1;dx<=1&&ok;++dx)
                    for (int dy=-1;dy<=1&&ok;++dy)
                        if ((dx||dy) && g[ix+dx][iy+dy].err < e) ok = false;
                if (ok) seeds.push_back(g[ix][iy]);
            }
    };
    harvest(g0); harvest(g1);
    std::sort(seeds.begin(), seeds.end(), [](const Cell& a, const Cell& b){ return a.err < b.err; });
    std::cerr << "[3d_sitnikov m3=" << m3 << "] seeds=" << seeds.size() << "\n";

    std::mutex mtx;
    std::vector<std::array<double, 3>> dedup;
    std::atomic<size_t> next_idx{0};
    ths.clear();
    for (int tid = 0; tid < nthreads; ++tid) {
        ths.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1);
                if (idx >= seeds.size()) break;
                double vy = seeds[idx].vy, vz = seeds[idx].vz, Th = seeds[idx].Th, ferr = 0;
                int type = seeds[idx].type;
                if (!refine_sitnikov(vy, vz, Th, m3, type, ferr)) continue;

                // Full period is 2*Th (for type 0) or 4*Th (check both 2*Th and 4*Th!)
                State3D y0 = make_sitnikov3d(vy, vz, m3);
                double T = 2.0 * Th;
                State3D yT = y0;
                double rmin = 1e300;
                if (!integrate3d(yT, T, m, 1e-13, 250000, &rmin, 5e-3)) continue;
                double ret = state_dist3d(yT, y0);
                if (ret > 1e-6) {
                    T = 4.0 * Th;
                    yT = y0;
                    if (!integrate3d(yT, T, m, 1e-13, 250000, &rmin, 5e-3)) continue;
                    ret = state_dist3d(yT, y0);
                    if (ret > 1e-6) continue;
                }

                // Topology & true 3D thickness
                const int N = 320;
                double dt = T / N, C[3][3] = {{0}}, rmax = 0;
                int prev_cl = -1, swaps = 0, cnt = 0;
                State3D yc = y0;
                for (int k = 0; k <= N; ++k) {
                    double r12 = std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
                    double r13 = std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
                    double r23 = std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
                    rmax = std::max({rmax, r12, r13, r23});
                    int cl = (r12 <= r13 && r12 <= r23) ? 0 : ((r13 <= r23) ? 1 : 2);
                    if (prev_cl != -1 && cl != prev_cl) swaps++;
                    prev_cl = cl;
                    for (int b = 0; b < 3; ++b) {
                        double p[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
                        for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += p[r]*p[c];
                        cnt++;
                    }
                    if (k < N) integrate3d(yc, dt, m, 1e-11, 50000);
                }
                for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
                double thick = std::sqrt(min_eig_3x3(C));
                if (thick < 0.08) continue;

                {
                    std::lock_guard<std::mutex> lk(mtx);
                    bool dup = false;
                    for (const auto& f : dedup)
                        if (std::hypot(f[0]-vy, f[1]-vz) < 2e-3 && std::abs(f[2]-T) < 0.1) { dup = true; break; }
                    if (dup) continue;
                    dedup.push_back({vy, vz, T});
                }

                auto st = analyze_stab3d(y0, T, m);
                if (!st.stable || st.surv100 < 100) continue;

                double E = energy3d(y0, m);
                auto Lv = angmom3d(y0, m);
                double Lmag = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
                std::ostringstream os;
                os << std::setprecision(13);
                os << "{\"dim\":3,\"subspace\":\"3d_sitnikov_C2\",\"masses\":[1,1," << m3 << "],\"y0\":[";
                for (int i=0;i<18;++i) os << y0[i] << (i<17?",":"");
                os << "],\"T\":" << T << ",\"E\":" << E << ",\"L\":" << Lmag
                   << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]"
                   << ",\"ret_err\":" << ret
                   << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2 << ",\"k3\":" << st.k3
                   << ",\"stab_margin\":" << st.margin
                   << ",\"survived_100T\":" << st.surv100
                   << ",\"drift_100T\":" << st.drift100
                   << ",\"rmin\":" << rmin << ",\"rmax\":" << rmax
                   << ",\"thickness_3d\":" << thick
                   << ",\"partner_swaps\":" << swaps
                   << ",\"democratic\":" << (swaps >= 2 ? "true" : "false")
                   << "}\n";
                std::lock_guard<std::mutex> lk(mtx);
                std::cout << os.str() << std::flush;
            }
        });
    }
    for (auto& th : ths) th.join();
    return 0;
}
