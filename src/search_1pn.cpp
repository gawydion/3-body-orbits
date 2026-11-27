#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// Canonical phase-space state for 3 bodies in 3D:
// s[0..8]  = r1(x,y,z), r2(x,y,z), r3(x,y,z)
// s[9..17] = p1(x,y,z), p2(x,y,z), p3(x,y,z)  (canonical momenta!)
using State1PN = std::array<double, 18>;
using Masses = std::array<double, 3>;

// Exact 1PN ADM / Einstein-Infeld-Hoffmann Hamiltonian (G = 1, beta = 1/c^2):
// H = sum p_A^2/(2 m_A) - sum_{A<B} m_A m_B / r_AB
//   + beta * [ -sum p_A^4/(8 m_A^3)
//              - 1/2 sum_{A<B} 1/r_AB * (3*(m_B/m_A p_A^2 + m_A/m_B p_B^2) - 7 p_A.p_B - (n_AB.p_A)(n_AB.p_B))
//              + 1/2 sum_A m_A * (sum_{B!=A} m_B / r_AB)^2 ]
static double hamiltonian_1pn(const State1PN& s, const Masses& m, double beta) {
    double H0 = 0.0, H1 = 0.0;
    double phi[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a) {
        double px = s[9+3*a], py = s[10+3*a], pz = s[11+3*a];
        double p2 = px*px + py*py + pz*pz;
        H0 += 0.5 * p2 / m[a];
        H1 -= (p2 * p2) / (8.0 * m[a] * m[a] * m[a]);
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            double dx = s[3*a] - s[3*b], dy = s[3*a+1] - s[3*b+1], dz = s[3*a+2] - s[3*b+2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            double ir = 1.0 / r;
            H0 -= m[a] * m[b] * ir;
            phi[a] += m[b] * ir;
            phi[b] += m[a] * ir;
            if (beta != 0.0) {
                double pax = s[9+3*a], pay = s[10+3*a], paz = s[11+3*a];
                double pbx = s[9+3*b], pby = s[10+3*b], pbz = s[11+3*b];
                double pa2 = pax*pax + pay*pay + paz*paz;
                double pb2 = pbx*pbx + pby*pby + pbz*pbz;
                double papb = pax*pbx + pay*pby + paz*pbz;
                double nx = dx * ir, ny = dy * ir, nz = dz * ir;
                double na = nx*pax + ny*pay + nz*paz;
                double nb = nx*pbx + ny*pby + nz*pbz;
                double C = 3.0 * ((m[b]/m[a])*pa2 + (m[a]/m[b])*pb2) - 7.0 * papb;
                H1 += -0.5 * ir * (C - na * nb);
            }
        }
    }
    if (beta != 0.0) {
        for (int a = 0; a < 3; ++a) H1 += 0.5 * m[a] * phi[a] * phi[a];
    }
    return H0 + beta * H1;
}

static inline void deriv_1pn(const State1PN& s, const Masses& m, double beta,
                             State1PN& ds, double* rmin_out = nullptr) {
    double dr[9] = {0}, dp[9] = {0};
    double phi[3] = {0, 0, 0};
    double dx_ab[3][3], dy_ab[3][3], dz_ab[3][3], ir_ab[3][3], ir3_ab[3][3];
    double rmin2 = 1e300;

    // 1. Kinetic terms (Newtonian + 1PN T1)
    for (int a = 0; a < 3; ++a) {
        double px = s[9+3*a], py = s[10+3*a], pz = s[11+3*a];
        double p2 = px*px + py*py + pz*pz;
        double inv_m = 1.0 / m[a];
        double fac = inv_m - beta * (0.5 * p2 * inv_m * inv_m * inv_m);
        dr[3*a]   = fac * px;
        dr[3*a+1] = fac * py;
        dr[3*a+2] = fac * pz;
    }

    // 2. Pairwise terms: Newtonian potential + 1PN V1 + accumulate phi_A
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            double dx = s[3*a] - s[3*b], dy = s[3*a+1] - s[3*b+1], dz = s[3*a+2] - s[3*b+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2);
            double ir2 = ir * ir, ir3 = ir2 * ir;
            dx_ab[a][b] = dx; dy_ab[a][b] = dy; dz_ab[a][b] = dz;
            ir_ab[a][b] = ir; ir3_ab[a][b] = ir3;

            double mm = m[a] * m[b];
            // Newtonian dp/dt = -dH0/dr
            double fx0 = -mm * dx * ir3, fy0 = -mm * dy * ir3, fz0 = -mm * dz * ir3;
            dp[3*a] += fx0; dp[3*a+1] += fy0; dp[3*a+2] += fz0;
            dp[3*b] -= fx0; dp[3*b+1] -= fy0; dp[3*b+2] -= fz0;

            phi[a] += m[b] * ir;
            phi[b] += m[a] * ir;

            if (beta != 0.0) {
                double pax = s[9+3*a], pay = s[10+3*a], paz = s[11+3*a];
                double pbx = s[9+3*b], pby = s[10+3*b], pbz = s[11+3*b];
                double pa2 = pax*pax + pay*pay + paz*paz;
                double pb2 = pbx*pbx + pby*pby + pbz*pbz;
                double papb = pax*pbx + pay*pby + paz*pbz;
                double nx = dx * ir, ny = dy * ir, nz = dz * ir;
                double na = nx*pax + ny*pay + nz*paz;
                double nb = nx*pbx + ny*pby + nz*pbz;

                // dr_A/dt += beta * dV1/dp_A
                double c_ab = 0.5 * beta * ir;
                double mba = m[b] / m[a], mab = m[a] / m[b];
                dr[3*a]   += c_ab * (-6.0*mba*pax + 7.0*pbx + nb*nx);
                dr[3*a+1] += c_ab * (-6.0*mba*pay + 7.0*pby + nb*ny);
                dr[3*a+2] += c_ab * (-6.0*mba*paz + 7.0*pbz + nb*nz);

                dr[3*b]   += c_ab * (-6.0*mab*pbx + 7.0*pax + na*nx);
                dr[3*b+1] += c_ab * (-6.0*mab*pby + 7.0*pay + na*ny);
                dr[3*b+2] += c_ab * (-6.0*mab*pbz + 7.0*paz + na*nz);

                // dp_A/dt += -beta * dV1/dr_A
                double C = 3.0 * (mba*pa2 + mab*pb2) - 7.0 * papb;
                double coef_n = C - 3.0 * na * nb;
                double f_scale = -0.5 * beta * ir2;
                double fx1 = f_scale * (coef_n * nx + nb * pax + na * pbx);
                double fy1 = f_scale * (coef_n * ny + nb * pay + na * pby);
                double fz1 = f_scale * (coef_n * nz + nb * paz + na * pbz);
                dp[3*a] += fx1; dp[3*a+1] += fy1; dp[3*a+2] += fz1;
                dp[3*b] -= fx1; dp[3*b+1] -= fy1; dp[3*b+2] -= fz1;
            }
        }
    }

    // 3. Static 3-body 1PN potential U1 = 1/2 sum_A m_A phi_A^2
    if (beta != 0.0) {
        for (int a = 0; a < 3; ++a) {
            for (int b = a + 1; b < 3; ++b) {
                double f = beta * m[a] * m[b] * (phi[a] + phi[b]) * ir3_ab[a][b];
                double fx = f * dx_ab[a][b], fy = f * dy_ab[a][b], fz = f * dz_ab[a][b];
                dp[3*a] += fx; dp[3*a+1] += fy; dp[3*a+2] += fz;
                dp[3*b] -= fx; dp[3*b+1] -= fy; dp[3*b+2] -= fz;
            }
        }
    }

    for (int i = 0; i < 9; ++i) { ds[i] = dr[i]; ds[9+i] = dp[i]; }
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline std::array<double, 3> angmom_1pn(const State1PN& s) {
    double Lx = 0, Ly = 0, Lz = 0;
    for (int i = 0; i < 3; ++i) {
        Lx += s[3*i+1]*s[11+3*i] - s[3*i+2]*s[10+3*i];
        Ly += s[3*i+2]*s[9+3*i]  - s[3*i]*s[11+3*i];
        Lz += s[3*i]*s[10+3*i]   - s[3*i+1]*s[9+3*i];
    }
    return {Lx, Ly, Lz};
}

static inline double state_dist(const State1PN& a, const State1PN& b) {
    double s = 0; for (int i=0;i<18;++i) { double d = a[i]-b[i]; s += d*d; }
    return std::sqrt(s);
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

static State1PN rotate_state(const State1PN& s, const double u[3], double ang) {
    State1PN out = s;
    for (int b = 0; b < 6; ++b) rot_vec(u, ang, out[3*b], out[3*b+1], out[3*b+2]);
    return out;
}

static bool dopri5_step(State1PN& y, double& t, double& h, double t_end, const Masses& m,
                        double beta, double tol, double& rmin_track, double r_abort = 1e-4) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;
    State1PN k1,k2,k3,k4,k5,k6,k7,yt,yn;
    double rm = 1e300;
    for (int att=0;att<25;++att) {
        deriv_1pn(y, m, beta, k1, &rm); if (rm < r_abort) return false;
        for (int i=0;i<18;++i) yt[i]=y[i]+h*a21*k1[i]; deriv_1pn(yt, m, beta, k2);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv_1pn(yt, m, beta, k3);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv_1pn(yt, m, beta, k4);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv_1pn(yt, m, beta, k5);
        for (int i=0;i<18;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv_1pn(yt, m, beta, k6);
        for (int i=0;i<18;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv_1pn(yn, m, beta, k7, &rm); if (rm < r_abort) return false;
        double err = 0;
        for (int i=0;i<18;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err/18.0);
        if (err <= 1.0) {
            y = yn; t += h; if (rm < rmin_track) rmin_track = rm;
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.08);
            return true;
        }
        h *= std::max(0.1, 0.9*std::pow(err,-0.25));
        if (h < 1e-13) return false;
    }
    return false;
}

static bool integrate_1pn(State1PN& y, double T, const Masses& m, double beta,
                          double tol = 1e-13, int max_steps = 800000, double* rmin_out = nullptr) {
    double t = 0, h = std::min(1e-3, T*0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step(y, t, h, T, m, beta, tol, rmin_track)) return false;
        ++steps;
        for (int i=0;i<9;++i) if (std::abs(y[i]) > 100.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

template <int NEQ, int NVAR>
static bool solve_lm(const double J[NEQ][NVAR], const double F[NEQ], double lambda, double d[NVAR]) {
    double A[NVAR][NVAR + 1] = {{0}};
    for (int r = 0; r < NVAR; ++r) {
        for (int c = 0; c < NVAR; ++c) {
            double s = 0; for (int i=0;i<NEQ;++i) s += J[i][r]*J[i][c];
            A[r][c] = s;
        }
        A[r][r] += lambda*(1.0 + A[r][r]);
        double s = 0; for (int i=0;i<NEQ;++i) s -= J[i][r]*F[i];
        A[r][NVAR] = s;
    }
    for (int col = 0; col < NVAR; ++col) {
        int piv = col;
        for (int r = col+1; r < NVAR; ++r) if (std::abs(A[r][col]) > std::abs(A[piv][col])) piv = r;
        if (std::abs(A[piv][col]) < 1e-20) return false;
        for (int c = col; c <= NVAR; ++c) std::swap(A[col][c], A[piv][c]);
        for (int r = col+1; r < NVAR; ++r) {
            double f = A[r][col] / A[col][col];
            for (int c = col; c <= NVAR; ++c) A[r][c] -= f*A[col][c];
        }
    }
    for (int r = NVAR-1; r >= 0; --r) {
        double s = A[r][NVAR];
        for (int c = r+1; c < NVAR; ++c) s -= A[r][c]*d[c];
        d[r] = s / A[r][r];
    }
    return true;
}

// Convert canonical state to physical velocities v = dH/dp (for exporting/visualization & v^2/c^2 measurement)
static std::array<double, 9> canonical_velocities(const State1PN& s, const Masses& m, double beta) {
    State1PN ds;
    deriv_1pn(s, m, beta, ds);
    std::array<double, 9> v;
    for (int i = 0; i < 9; ++i) v[i] = ds[i];
    return v;
}

// Canonical collinear L=0 state (m1=m2=1, m3=m3): r1=(-1,0,0), r2=(1,0,0), r3=(0,0,0)
// Optional (Px, Py) adds zero-Delta-L total canonical momentum P_A += (m_A/M)*(Px,Py) to cancel 1PN CM drift!
static State1PN make_col_1pn(double px, double py, double Px = 0.0, double Py = 0.0, double m3 = 1.0) {
    double M = 2.0 + m3;
    double w1 = 1.0 / M, w3 = m3 / M;
    return {-1,0,0, 1,0,0, 0,0,0,
            px + w1*Px, py + w1*Py, 0,
            px + w1*Px, py + w1*Py, 0,
            -2*px + w3*Px, -2*py + w3*Py, 0};
}

// Canonical perp RPO state (x1=-1+xc, x2=x2+xc, x3=x3_newt+xc, pz1=0)
// Optional P_par adds canonical momentum strictly along u = L/|L| with ZERO Delta L!
static State1PN make_rpo_1pn(double x2, double py1, double py2, double pz2,
                             const Masses& m, double xc = 0.0, double P_par = 0.0) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2) / m[2];
    double py3 = -(py1 + py2);
    double pz3 = -pz2;
    State1PN s = {x1+xc,0,0, x2+xc,0,0, x3+xc,0,0, 0,py1,0, 0,py2,pz2, 0,py3,pz3};
    if (std::abs(P_par) > 0.0) {
        auto Lv = angmom_1pn(s);
        double Lm = std::sqrt(Lv[0]*Lv[0] + Lv[1]*Lv[1] + Lv[2]*Lv[2]);
        if (Lm > 1e-12) {
            double uy = Lv[1] / Lm, uz = Lv[2] / Lm;
            double xn[3] = {x1, x2, x3};
            double M = m[0] + m[1] + m[2];
            double I0 = m[0]*x1*x1 + m[1]*x2*x2 + m[2]*x3*x3;
            for (int a = 0; a < 3; ++a) {
                double w = m[a] / M - (xc / I0) * m[a] * xn[a];
                s[10 + 3*a] += w * P_par * uy;
                s[11 + 3*a] += w * P_par * uz;
            }
        }
    }
    return s;
}

// Compute exact 1PN center-of-mass shift xc so sum_A E_A * x_A = 0
static double compute_1pn_xc(double x2, double py1, double py2, double pz2, const Masses& m, double beta) {
    if (beta == 0.0) return 0.0;
    State1PN s = make_rpo_1pn(x2, py1, py2, pz2, m, 0.0);
    double phi[3] = {0, 0, 0};
    for (int a=0;a<3;++a) for (int b=a+1;b<3;++b) {
        double r = std::abs(s[3*a] - s[3*b]);
        phi[a] += m[b]/r; phi[b] += m[a]/r;
    }
    double num = 0.0, den = 0.0;
    for (int a=0;a<3;++a) {
        double p2 = s[9+3*a]*s[9+3*a] + s[10+3*a]*s[10+3*a] + s[11+3*a]*s[11+3*a];
        double ma_eff = m[a] + beta * (0.5 * p2 / m[a] - 0.5 * m[a] * phi[a]);
        num += ma_eff * s[3*a];
        den += ma_eff;
    }
    return -num / den;
}

static bool verify_physics() {
    std::cout << std::setprecision(12);
    // Test 1: Analytical dH/dr, dH/dp vs central finite difference on asymmetric 3D state
    Masses m = {1.0, 0.7, 0.4};
    double beta = 0.02;
    State1PN s = {-1.0, 0.15, -0.22,  0.85, -0.31, 0.19,  0.21, 0.44, -0.07,
                   0.31, 0.42, -0.18, -0.52, 0.11, 0.37,  0.21, -0.53, -0.19};
    State1PN ds;
    deriv_1pn(s, m, beta, ds);
    double eps = 1e-7, max_grad_err = 0.0;
    for (int i = 0; i < 9; ++i) {
        State1PN sp = s, sm = s; sp[9+i] += eps; sm[9+i] -= eps;
        double num_dr = (hamiltonian_1pn(sp, m, beta) - hamiltonian_1pn(sm, m, beta)) / (2*eps);
        max_grad_err = std::max(max_grad_err, std::abs(ds[i] - num_dr));

        sp = s; sm = s; sp[i] += eps; sm[i] -= eps;
        double num_dp = -(hamiltonian_1pn(sp, m, beta) - hamiltonian_1pn(sm, m, beta)) / (2*eps);
        max_grad_err = std::max(max_grad_err, std::abs(ds[9+i] - num_dp));
    }

    // Test 2: Exact conservation of H_1PN, P, L over T = 20.0
    State1PN st = s;
    double H0 = hamiltonian_1pn(st, m, beta);
    auto L0 = angmom_1pn(st);
    integrate_1pn(st, 20.0, m, beta, 1e-13, 400000);
    double dH = std::abs(hamiltonian_1pn(st, m, beta) - H0);
    auto Lt = angmom_1pn(st);
    double dL = std::sqrt((Lt[0]-L0[0])*(Lt[0]-L0[0]) + (Lt[1]-L0[1])*(Lt[1]-L0[1]) + (Lt[2]-L0[2])*(Lt[2]-L0[2]));
    double dP = std::sqrt(std::pow(st[9]+st[12]+st[15], 2) + std::pow(st[10]+st[13]+st[16], 2) + std::pow(st[11]+st[14]+st[17], 2));

    // Test 3: 2-body GR perihelion precession Delta_phi = 6*pi*M*beta / (a*(1-e^2))
    auto test_precess = [&](double m1, double m2, double btest) -> double {
        Masses m2b = {m1, m2, 1e-18};
        double M = m1 + m2, mu = m1*m2/M;
        double a = 1.0, e = 0.4;
        double r_peri = a * (1.0 - e);
        double v_peri = std::sqrt(M * (1.0 + e) / r_peri);
        double p_peri = mu * v_peri;
        State1PN y = {-m2/M * r_peri, 0, 0,  m1/M * r_peri, 0, 0,  20.0, 0, 0,
                      0, -p_peri, 0,         0, p_peri, 0,         0, 0, 0};
        // Find first pericenter return (dr/dt crossing from negative to positive)
        double t = 0, h = 1e-3, rm = 1e300;
        double prev_rdot = 0;
        bool left = false;
        double phi_peri = 0;
        for (int step = 0; step < 200000; ++step) {
            State1PN yprev = y;
            double tprev = t;
            dopri5_step(y, t, h, 20.0, m2b, btest, 1e-14, rm);
            double rx = y[3]-y[0], ry = y[4]-y[1];
            State1PN dy; deriv_1pn(y, m2b, btest, dy);
            double rdot = rx*(dy[3]-dy[0]) + ry*(dy[4]-dy[1]);
            if (t > 1.0 && rdot < -0.05) left = true;
            if (left && prev_rdot < 0 && rdot >= 0) {
                // Linear/Newton root on rdot = 0
                double lo = 0, hi = t - tprev;
                for (int bi = 0; bi < 40; ++bi) {
                    double mid = 0.5*(lo + hi);
                    State1PN ym = yprev;
                    integrate_1pn(ym, mid, m2b, btest, 1e-14, 10000);
                    State1PN dym; deriv_1pn(ym, m2b, btest, dym);
                    double rm_dot = (ym[3]-ym[0])*(dym[3]-dym[0]) + (ym[4]-ym[1])*(dym[4]-dym[1]);
                    if (rm_dot < 0) lo = mid; else hi = mid;
                }
                State1PN ym = yprev;
                integrate_1pn(ym, 0.5*(lo+hi), m2b, btest, 1e-14, 10000);
                phi_peri = std::atan2(ym[4]-ym[1], ym[3]-ym[0]);
                break;
            }
            prev_rdot = rdot;
        }
        double expected = 6.0 * M_PI * M * btest / (a * (1.0 - e*e));
        return std::abs(phi_peri - expected) / expected;
    };
    double rel_tp = test_precess(1.0, 1e-6, 1e-4);
    double rel_eq = test_precess(0.5, 0.5, 1e-4);

    std::cout << "1PN Physics Self-Test:\n"
              << "  Hamiltonian gradient max err = " << max_grad_err << "\n"
              << "  Energy conservation |dH|     = " << dH << "\n"
              << "  AngMom conservation |dL|     = " << dL << "\n"
              << "  LinMom conservation |dP|     = " << dP << "\n"
              << "  2-body GR precession rel err (test particle) = " << rel_tp << "\n"
              << "  2-body GR precession rel err (equal mass)    = " << rel_eq << "\n";
    return (max_grad_err < 1e-8 && dH < 1e-11 && dL < 1e-11 && dP < 1e-12 && rel_tp < 5e-3 && rel_eq < 5e-3);
}

// Compute effective orbital beta_eff = <v^2> * beta and geometry
static double min_eig_3x3(double C[3][3]) {
    double I1 = C[0][0]+C[1][1]+C[2][2];
    double I2 = C[0][0]*C[1][1]+C[0][0]*C[2][2]+C[1][1]*C[2][2]
              - C[0][1]*C[0][1]-C[0][2]*C[0][2]-C[1][2]*C[1][2];
    double I3 = C[0][0]*(C[1][1]*C[2][2]-C[1][2]*C[1][2])
              - C[0][1]*(C[0][1]*C[2][2]-C[1][2]*C[0][2])
              + C[0][2]*(C[0][1]*C[1][2]-C[1][1]*C[0][2]);
    double p = I2 - I1*I1/3.0, q = -2.0*I1*I1*I1/27.0 + I1*I2/3.0 - I3;
    double r = std::sqrt(std::max(0.0, -p/3.0));
    if (r < 1e-15) return std::max(0.0, I1/3.0);
    double phi = std::acos(std::clamp(-q/(2.0*r*r*r), -1.0, 1.0)) / 3.0;
    return std::max(0.0, std::min({I1/3.0+2*r*std::cos(phi),
                                   I1/3.0+2*r*std::cos(phi+2*M_PI/3),
                                   I1/3.0+2*r*std::cos(phi+4*M_PI/3)}));
}

struct Geom1PN { double rmin, rmax, thick, v2_rms, v2_max; int swaps; };
static Geom1PN measure_geom(const State1PN& y0, double T, const Masses& m, double beta) {
    const int N = 260;
    double dt = T / N, C[3][3] = {{0}}, rmin = 1e300, rmax = 0, v2_sum = 0, v2_max = 0;
    int prev_cl = -1, swaps = 0, cnt = 0;
    State1PN yc = y0;
    for (int k = 0; k <= N; ++k) {
        double r12 = std::sqrt((yc[3]-yc[0])*(yc[3]-yc[0])+(yc[4]-yc[1])*(yc[4]-yc[1])+(yc[5]-yc[2])*(yc[5]-yc[2]));
        double r13 = std::sqrt((yc[6]-yc[0])*(yc[6]-yc[0])+(yc[7]-yc[1])*(yc[7]-yc[1])+(yc[8]-yc[2])*(yc[8]-yc[2]));
        double r23 = std::sqrt((yc[6]-yc[3])*(yc[6]-yc[3])+(yc[7]-yc[4])*(yc[7]-yc[4])+(yc[8]-yc[5])*(yc[8]-yc[5]));
        rmin = std::min({rmin, r12, r13, r23});
        rmax = std::max({rmax, r12, r13, r23});
        int cl = (r12<=r13 && r12<=r23)?0:((r13<=r23)?1:2);
        if (prev_cl != -1 && cl != prev_cl) swaps++;
        prev_cl = cl;
        auto v = canonical_velocities(yc, m, beta);
        for (int b=0;b<3;++b) {
            double v2 = v[3*b]*v[3*b] + v[3*b+1]*v[3*b+1] + v[3*b+2]*v[3*b+2];
            v2_sum += v2; if (v2 > v2_max) v2_max = v2;
            double qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
            for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += qv[r]*qv[c];
            cnt++;
        }
        if (k < N) integrate_1pn(yc, dt, m, beta, 1e-11, 80000);
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
    return {rmin, rmax, std::sqrt(min_eig_3x3(C)), v2_sum / cnt, v2_max, swaps};
}

// Exact Real Schur / Francis Double-Shift QR Eigenvalue Solver (General NxN)
static std::vector<std::complex<double>> real_eigenvalues(std::vector<std::vector<double>> H) {
    int n = (int)H.size();
    for (int k = 0; k < n - 2; ++k) {
        double scale = 0.0;
        for (int i = k + 1; i < n; ++i) scale += std::abs(H[i][k]);
        if (scale < 1e-30) continue;
        double norm2 = 0.0;
        std::vector<double> v(n - k - 1);
        for (int i = k + 1; i < n; ++i) { v[i - k - 1] = H[i][k] / scale; norm2 += v[i - k - 1]*v[i - k - 1]; }
        double alpha = (v[0] >= 0 ? -std::sqrt(norm2) : std::sqrt(norm2));
        v[0] -= alpha;
        double v2 = 0.0; for (double x : v) v2 += x*x;
        if (v2 < 1e-30) continue;
        double beta = 2.0 / v2;
        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k + 1; i < n; ++i) dot += v[i - k - 1] * H[i][j];
            for (int i = k + 1; i < n; ++i) H[i][j] -= beta * dot * v[i - k - 1];
        }
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;
            for (int j = k + 1; j < n; ++j) dot += H[i][j] * v[j - k - 1];
            for (int j = k + 1; j < n; ++j) H[i][j] -= beta * dot * v[j - k - 1];
        }
        for (int i = k + 2; i < n; ++i) H[i][k] = 0.0;
    }
    std::vector<std::complex<double>> eigs;
    double norm = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = std::max(0, i - 1); j < n; ++j) norm += std::abs(H[i][j]);
    const double eps = 2.220446049250313e-16;
    int nn = n - 1, iter = 0;
    double exshift = 0.0;
    while (nn >= 0) {
        int l = nn;
        while (l > 0) {
            double s = std::abs(H[l - 1][l - 1]) + std::abs(H[l][l]);
            if (s == 0.0) s = norm;
            if (std::abs(H[l][l - 1]) <= eps * s) { H[l][l - 1] = 0.0; break; }
            --l;
        }
        if (l == nn) {
            eigs.emplace_back(H[nn][nn], 0.0);
            --nn; iter = 0;
        } else if (l == nn - 1) {
            double p = 0.5 * (H[nn - 1][nn - 1] - H[nn][nn]);
            double q = p * p + H[nn][nn - 1] * H[nn - 1][nn];
            double mid = 0.5 * (H[nn - 1][nn - 1] + H[nn][nn]);
            if (q >= 0.0) {
                double sq = std::sqrt(q);
                eigs.emplace_back(mid + sq, 0.0);
                eigs.emplace_back(mid - sq, 0.0);
            } else {
                double sq = std::sqrt(-q);
                eigs.emplace_back(mid, sq);
                eigs.emplace_back(mid, -sq);
            }
            nn -= 2; iter = 0;
        } else {
            if (iter > 400) break;
            double x = H[nn][nn], y = H[nn - 1][nn - 1];
            double w = H[nn][nn - 1] * H[nn - 1][nn];
            if (iter == 10 || iter == 20 || iter == 40 || iter == 80 || iter == 160) {
                exshift += x;
                double s = std::abs(H[nn][nn - 1]) + std::abs(H[nn - 1][nn - 2]);
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }
            ++iter;
            int m = nn - 2;
            double p = 0, q = 0, r = 0;
            for (; m >= l; --m) {
                double z = H[m][m];
                r = x - z; double s = y - z;
                p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
                q = H[m + 1][m + 1] - z - r - s;
                r = H[m + 2][m + 1];
                double sc = std::abs(p) + std::abs(q) + std::abs(r);
                p /= sc; q /= sc; r /= sc;
                if (m == l) break;
                double lhs = std::abs(H[m][m - 1]) * (std::abs(q) + std::abs(r));
                double rhs = std::abs(p) * (std::abs(H[m - 1][m - 1]) + std::abs(z) + std::abs(H[m + 1][m + 1]));
                if (lhs <= eps * rhs) break;
            }
            for (int i = m + 2; i <= nn; ++i) {
                H[i][i - 2] = 0.0;
                if (i != m + 2) H[i][i - 3] = 0.0;
            }
            for (int k = m; k <= nn - 1; ++k) {
                bool notlast = (k != nn - 1);
                if (k != m) {
                    p = H[k][k - 1]; q = H[k + 1][k - 1];
                    r = notlast ? H[k + 2][k - 1] : 0.0;
                    double sc = std::abs(p) + std::abs(q) + std::abs(r);
                    if (sc == 0.0) continue;
                    p /= sc; q /= sc; r /= sc;
                }
                double s = std::sqrt(p*p + q*q + r*r);
                if (p < 0) s = -s;
                if (s != 0.0) {
                    if (k != m) H[k][k - 1] = -s * (std::abs(H[k][k - 1]) + std::abs(H[k + 1][k - 1]) + (notlast ? std::abs(H[k + 2][k - 1]) : 0.0));
                    else if (l != m) H[k][k - 1] = -H[k][k - 1];
                    p += s;
                    double xp = p / s, yp = q / s, zp = r / s;
                    q /= p; r /= p;
                    for (int j = k; j < n; ++j) {
                        double pp = H[k][j] + q * H[k + 1][j];
                        if (notlast) { pp += r * H[k + 2][j]; H[k + 2][j] -= pp * zp; }
                        H[k][j] -= pp * xp;
                        H[k + 1][j] -= pp * yp;
                    }
                    int rmax = std::min(nn, k + 3);
                    for (int i = 0; i <= rmax; ++i) {
                        double pp = xp * H[i][k] + yp * H[i][k + 1];
                        if (notlast) { pp += zp * H[i][k + 2]; H[i][k + 2] -= pp * r; }
                        H[i][k] -= pp;
                        H[i][k + 1] -= pp * q;
                    }
                }
            }
        }
    }
    return eigs;
}

struct Stab1PN {
    bool stable_2d, stable_3d;
    double k1, k2, k3; // for 2D, k3 = kz
    double margin_2d, margin_3d;
    double max_rad_err;
    double symp_err;
    int surv100;
    double drift100;
};

static bool extract_spectrum_1pn(const std::vector<std::vector<double>>& M, double two_alpha, int dim,
                                 std::vector<double>& ks, double& max_rad_err) {
    auto eigs = real_eigenvalues(M);
    int N = (int)M.size();
    if ((int)eigs.size() != N) return false;
    std::vector<std::complex<double>> expected;
    std::complex<double> rp(std::cos(two_alpha), std::sin(two_alpha));
    std::complex<double> rm(std::cos(two_alpha), -std::sin(two_alpha));
    if (dim == 2) {
        for (int i = 0; i < 4; ++i) expected.emplace_back(1.0, 0.0);
        for (int i = 0; i < 2; ++i) { expected.push_back(rp); expected.push_back(rm); }
    } else {
        for (int i = 0; i < 6; ++i) expected.emplace_back(1.0, 0.0);
        for (int i = 0; i < 3; ++i) { expected.push_back(rp); expected.push_back(rm); }
    }
    std::vector<bool> used(N, false);
    for (const auto& t : expected) {
        int best = -1; double bd = 1e300;
        for (int i = 0; i < N; ++i) {
            if (used[i]) continue;
            double d = std::abs(eigs[i] - t);
            if (d < bd) { bd = d; best = i; }
        }
        if (best >= 0) used[best] = true;
    }
    std::vector<std::complex<double>> phys;
    for (int i = 0; i < N; ++i) if (!used[i]) phys.push_back(eigs[i]);
    max_rad_err = 0.0;
    for (const auto& e : phys) max_rad_err = std::max(max_rad_err, std::abs(std::abs(e) - 1.0));
    ks.clear();
    std::vector<bool> pused(phys.size(), false);
    for (size_t i = 0; i < phys.size(); ++i) {
        if (pused[i]) continue;
        pused[i] = true;
        int mate = -1; double bd = 1e300;
        for (size_t j = i + 1; j < phys.size(); ++j) {
            if (pused[j]) continue;
            double d = std::min(std::abs(phys[j] - std::conj(phys[i])),
                                std::abs(phys[j] * phys[i] - std::complex<double>(1.0, 0.0)));
            if (d < bd) { bd = d; mate = (int)j; }
        }
        if (mate >= 0) {
            pused[mate] = true;
            ks.push_back((phys[i] + phys[mate]).real());
        }
    }
    std::sort(ks.begin(), ks.end(), std::greater<double>());
    return true;
}

static Stab1PN analyze_stab_1pn(const State1PN& y0, double T, const double u[3], double two_alpha,
                                const Masses& m, double beta, int dim, bool run_100T = true) {
    Stab1PN st{}; st.margin_2d = -1; st.margin_3d = -1; st.max_rad_err = 999.0; st.symp_err = 999.0;
    const double eps = 1e-7;
    const int max_steps = 3000000;
    auto check_symp = [](const std::vector<std::vector<double>>& M) -> double {
        int N = (int)M.size(), half = N / 2;
        double err = 0.0;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                double s = 0.0;
                for (int k = 0; k < half; ++k) s += M[k][i] * M[k + half][j] - M[k + half][i] * M[k][j];
                double tgt = (j == i + half) ? 1.0 : ((i == j + half) ? -1.0 : 0.0);
                err = std::max(err, std::abs(s - tgt));
            }
        return err;
    };
    if (dim == 2) {
        const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        std::vector<std::vector<double>> M(12, std::vector<double>(12, 0.0));
        for (int c = 0; c < 12; ++c) {
            int col = xyidx[c];
            State1PN yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            if (!integrate_1pn(yp, T, m, beta, 1e-13, max_steps) ||
                !integrate_1pn(ym, T, m, beta, 1e-13, max_steps)) return st;
            if (std::abs(two_alpha) > 1e-14) {
                yp = rotate_state(yp, u, -two_alpha);
                ym = rotate_state(ym, u, -two_alpha);
            }
            for (int r = 0; r < 12; ++r) M[r][c] = (yp[xyidx[r]] - ym[xyidx[r]]) / (2.0 * eps);
        }
        st.symp_err = check_symp(M);
        std::vector<double> ks;
        if (extract_spectrum_1pn(M, two_alpha, 2, ks, st.max_rad_err) && ks.size() >= 2) {
            st.k1 = ks[0]; st.k2 = ks[1];
            double mk = std::max(std::abs(st.k1), std::abs(st.k2));
            if (st.symp_err < 1e-4 && st.max_rad_err < 2e-3 && mk < 2.0 - 1e-4) {
                st.stable_2d = true;
                st.margin_2d = 2.0 - mk;
            } else {
                st.stable_2d = false;
                st.margin_2d = (st.max_rad_err < 2e-3) ? (2.0 - mk) : -st.max_rad_err;
            }
        }
        const int zidx[6] = {2, 5, 8, 11, 14, 17};
        double trz = 0;
        for (int k=0;k<6;++k) {
            int col = zidx[k];
            State1PN yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            integrate_1pn(yp, T, m, beta, 1e-13, max_steps);
            integrate_1pn(ym, T, m, beta, 1e-13, max_steps);
            trz += (yp[col] - ym[col]) / (2.0 * eps);
        }
        st.k3 = trz - (2.0 + 2.0 * std::cos(two_alpha));
        st.stable_3d = st.stable_2d && (std::abs(st.k3) < 2.0 - 1e-4);
        st.margin_3d = st.stable_3d ? (2.0 - std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)})) : -1.0;

        if (st.stable_2d && run_100T) {
            State1PN yp = y0;
            double kick = 1e-5;
            yp[9] += kick; yp[12] -= kick; // canonical P conserved!
            yp[10] += kick; yp[13] -= kick;
            double md = 0; int surv = 0;
            for (int k=0;k<100;++k) {
                if (!integrate_1pn(yp, T, m, beta, 1e-11, max_steps)) break;
                State1PN yr = (std::abs(two_alpha) > 1e-14) ? rotate_state(yp, u, -two_alpha*(k+1)) : yp;
                double d = state_dist(yr, y0);
                if (d > md) md = d; if (d > 5.0) break;
                surv = k + 1;
            }
            st.surv100 = surv; st.drift100 = md;
        }
    } else {
        std::vector<std::vector<double>> M(18, std::vector<double>(18, 0.0));
        for (int col=0;col<18;++col) {
            State1PN yp = y0, ym = y0; yp[col]+=eps; ym[col]-=eps;
            if (!integrate_1pn(yp, T, m, beta, 1e-13, max_steps) ||
                !integrate_1pn(ym, T, m, beta, 1e-13, max_steps)) return st;
            if (std::abs(two_alpha) > 1e-14) {
                yp = rotate_state(yp, u, -two_alpha);
                ym = rotate_state(ym, u, -two_alpha);
            }
            for (int r=0;r<18;++r) M[r][col] = (yp[r]-ym[r])/(2.0*eps);
        }
        st.symp_err = check_symp(M);
        std::vector<double> ks;
        if (extract_spectrum_1pn(M, two_alpha, 3, ks, st.max_rad_err) && ks.size() >= 3) {
            st.k1 = ks[0]; st.k2 = ks[1]; st.k3 = ks[2];
            double mk = std::max({std::abs(st.k1), std::abs(st.k2), std::abs(st.k3)});
            if (st.symp_err < 1e-4 && st.max_rad_err < 2e-3 && mk < 2.0 - 1e-4) {
                st.stable_3d = true; st.stable_2d = true;
                st.margin_3d = 2.0 - mk; st.margin_2d = st.margin_3d;
            } else {
                st.stable_3d = false; st.stable_2d = false;
                st.margin_3d = (st.max_rad_err < 2e-3) ? (2.0 - mk) : -st.max_rad_err;
                st.margin_2d = st.margin_3d;
            }
        }
        if (st.stable_3d && run_100T) {
            State1PN yp = y0;
            double kick = 1e-5;
            yp[9]+=kick; yp[12]-=kick; yp[10]+=kick; yp[13]-=kick; yp[11]+=kick; yp[14]-=kick;
            double md = 0; int surv = 0;
            for (int k=0;k<100;++k) {
                if (!integrate_1pn(yp, T, m, beta, 1e-11, max_steps)) break;
                State1PN yr = (std::abs(two_alpha) > 1e-14) ? rotate_state(yp, u, -two_alpha*(k+1)) : yp;
                double d = state_dist(yr, y0);
                if (d > md) md = d; if (d > 5.0) break;
                surv = k + 1;
            }
            st.surv100 = surv; st.drift100 = md;
        }
    }
    return st;
}

// Refine collinear L=0 1PN periodic/precessing orbit: p = {px, py, T, two_alpha, Px, Py}
static bool refine_col_1pn(double p[6], double m3, double beta, double& ferr) {
    Masses m = {1.0, 1.0, m3};
    const double uz[3] = {0, 0, 1};
    double lambda = 1e-5;
    auto eval = [&](const double q[6], double F[12]) -> bool {
        if (q[2] < 1.0 || q[2] > 50.0) return false;
        State1PN y0 = make_col_1pn(q[0], q[1], q[4], q[5], m3), yT = y0;
        if (!integrate_1pn(yT, q[2], m, beta, 1e-13, 350000)) return false;
        State1PN yr = (std::abs(q[3]) > 1e-15) ? rotate_state(yT, uz, -q[3]) : yT;
        const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        for (int i=0;i<12;++i) F[i] = yr[xyidx[i]] - y0[xyidx[i]];
        return true;
    };
    double F[12];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<12;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[12][6], Fp[12], Fm[12];
        for (int col=0;col<6;++col) {
            double pp[6], pm[6]; for(int i=0;i<6;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int i=0;i<12;++i) J[i][col] = (Fp[i]-Fm[i])/(2*eps);
        }
        double d[6]={0};
        if (!solve_lm<12,6>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn = 0; for(int i=0;i<6;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.15) for(int i=0;i<6;++i) d[i]*=(0.15/sn);
        double np[6], nF[12]; for(int i=0;i<6;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc = 0; for(int i=0;i<12;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<6;++i) p[i]=np[i]; for(int i=0;i<12;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    if (ferr >= 1e-9) std::cerr << "# refine_col_1pn failed beta=" << beta << " ferr=" << ferr << "\n";
    return ferr < 1e-9;
}

// Snap collinear L=0 1PN orbit to exact rational precession target_2a: p = {px, py, T, beta, Px, Py}
static bool snap_col_1pn(double p[6], double target_2a, double m3, double& ferr) {
    Masses m = {1.0, 1.0, m3};
    const double uz[3] = {0, 0, 1};
    double lambda = 1e-5;
    auto eval = [&](const double q[6], double F[12]) -> bool {
        if (q[2] < 1.0 || q[2] > 50.0 || q[3] < 0.0 || q[3] > 0.1) return false;
        State1PN y0 = make_col_1pn(q[0], q[1], q[4], q[5], m3), yT = y0;
        if (!integrate_1pn(yT, q[2], m, q[3], 1e-13, 400000)) return false;
        State1PN yr = rotate_state(yT, uz, -target_2a);
        const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        for (int i=0;i<12;++i) F[i] = yr[xyidx[i]] - y0[xyidx[i]];
        return true;
    };
    double F[12];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<12;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[12][6], Fp[12], Fm[12];
        for (int col=0;col<6;++col) {
            double hcol = (col == 3) ? 1e-8 : eps;
            double pp[6], pm[6]; for(int i=0;i<6;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=hcol; pm[col]-=hcol;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int i=0;i<12;++i) J[i][col] = (Fp[i]-Fm[i])/(2*hcol);
        }
        double d[6]={0};
        if (!solve_lm<12,6>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn = 0; for(int i=0;i<6;++i) sn += d[i]*d[i]; sn = std::sqrt(sn);
        if (sn > 0.1) for(int i=0;i<6;++i) d[i]*=(0.1/sn);
        double np[6], nF[12]; for(int i=0;i<6;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc = 0; for(int i=0;i<12;++i) nc += nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<6;++i) p[i]=np[i]; for(int i=0;i<12;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Refine 2D 1PN perp RPO at fixed (x2, beta): p = {py1, py2, Th, alpha, xc}
static bool refine_2d_rpo_1pn(double p[5], double x2, const Masses& m, double beta, double& ferr) {
    double lambda = 1e-5;
    const double uz[3] = {0, 0, 1};
    auto eval = [&](const double q[5], double F[5]) -> bool {
        if (q[2] < 0.5 || q[2] > 60.0) return false;
        State1PN y = make_rpo_1pn(x2, q[0], q[1], 0.0, m, q[4]);
        if (!integrate_1pn(y, q[2], m, beta, 1e-13, 400000)) return false;
        State1PN yt = rotate_state(y, uz, -q[3]);
        F[0] = yt[1]; F[1] = yt[4]; F[2] = yt[7]; F[3] = yt[9]; F[4] = yt[12];
        return true;
    };
    double F[5];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<5;++i) cost += F[i]*F[i];
    for (int it=0;it<32;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[5][5], Fp[5], Fm[5];
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<5;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5]={0};
        if (!solve_lm<5,5>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<5;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.2) for(int i=0;i<5;++i) d[i]*=(0.2/sn);
        double np[5], nF[5]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<5;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<5;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    if (ferr >= 1e-9) std::cerr << "# refine_2d_rpo_1pn failed beta=" << beta << " ferr=" << ferr << "\n";
    return ferr < 1e-9;
}

// Snap 2D 1PN rational closure at fixed (x2, target_alpha): p = {py1, py2, Th, beta, xc}
static bool snap_2d_1pn(double p[5], double x2, double target_a, const Masses& m, double& ferr) {
    double lambda = 1e-5;
    const double uz[3] = {0, 0, 1};
    auto eval = [&](const double q[5], double F[5]) -> bool {
        if (q[2] < 0.5 || q[2] > 60.0 || q[3] < 0.0 || q[3] > 0.2) return false;
        State1PN y = make_rpo_1pn(x2, q[0], q[1], 0.0, m, q[4]);
        if (!integrate_1pn(y, q[2], m, q[3], 1e-13, 400000)) return false;
        State1PN yt = rotate_state(y, uz, -target_a);
        F[0] = yt[1]; F[1] = yt[4]; F[2] = yt[7]; F[3] = yt[9]; F[4] = yt[12];
        return true;
    };
    double F[5];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<5;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[5][5], Fp[5], Fm[5];
        for (int col=0;col<5;++col) {
            double hcol = (col == 3) ? 1e-8 : eps;
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=hcol; pm[col]-=hcol;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<5;++r) J[r][col] = (Fp[r]-Fm[r])/(2*hcol);
        }
        double d[5]={0};
        if (!solve_lm<5,5>(J, F, lambda, d)) { lambda*=10; continue; }
        double np[5], nF[5]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<5;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<5;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Refine 3D 1PN perp RPO at fixed (pz2, m, beta): p = {x2, py1, py2, Th, alpha, xc, P_par}
static bool refine_3d_rpo_1pn(double p[7], double pz2, const Masses& m, double beta, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[7], double F[8]) -> bool {
        if (q[3] < 0.5 || q[3] > 60.0) return false;
        State1PN y0 = make_rpo_1pn(q[0], q[1], q[2], pz2, m, q[5], q[6]);
        auto Lv = angmom_1pn(y0);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State1PN y = y0;
        if (!integrate_1pn(y, q[3], m, beta, 1e-13, 500000)) return false;
        State1PN yt = rotate_state(y, u, -q[4]);
        F[0] = yt[1]; F[1] = yt[2]; F[2] = yt[4]; F[3] = yt[5];
        F[4] = yt[7]; F[5] = yt[8]; F[6] = yt[9]; F[7] = yt[12];
        return true;
    };
    double F[8];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<8;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[8][7], Fp[8], Fm[8];
        for (int col=0;col<7;++col) {
            double pp[7], pm[7]; for(int i=0;i<7;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<8;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[7]={0};
        if (!solve_lm<8,7>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<7;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.2) for(int i=0;i<7;++i) d[i]*=(0.2/sn);
        double np[7], nF[8]; for(int i=0;i<7;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<8;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<7;++i) p[i]=np[i]; for(int i=0;i<8;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    if (ferr >= 1e-9) std::cerr << "# refine_3d_rpo_1pn failed beta=" << beta << " ferr=" << ferr << "\n";
    return ferr < 1e-9;
}

// Snap 3D 1PN rational closure at fixed (m, beta, target_a): p = {x2, py1, py2, pz2, Th, xc, P_par}
static bool snap_3d_1pn(double p[7], double target_a, const Masses& m, double beta, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[7], double F[8]) -> bool {
        if (q[4] < 0.5 || q[4] > 60.0) return false;
        State1PN y0 = make_rpo_1pn(q[0], q[1], q[2], q[3], m, q[5], q[6]);
        auto Lv = angmom_1pn(y0);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State1PN y = y0;
        if (!integrate_1pn(y, q[4], m, beta, 1e-13, 500000)) return false;
        State1PN yt = rotate_state(y, u, -target_a);
        F[0] = yt[1]; F[1] = yt[2]; F[2] = yt[4]; F[3] = yt[5];
        F[4] = yt[7]; F[5] = yt[8]; F[6] = yt[9]; F[7] = yt[12];
        return true;
    };
    double F[8];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<8;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[8][7], Fp[8], Fm[8];
        for (int col=0;col<7;++col) {
            double pp[7], pm[7]; for(int i=0;i<7;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<8;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[7]={0};
        if (!solve_lm<8,7>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<7;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.2) for(int i=0;i<7;++i) d[i]*=(0.2/sn);
        double np[7], nF[8]; for(int i=0;i<7;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<8;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<7;++i) p[i]=np[i]; for(int i=0;i<8;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

static void emit_1pn_json(int dim, const std::string& label, const Masses& m, double beta,
                          const State1PN& s0, double T, double two_alpha, double ret_err,
                          const Stab1PN& st, const Geom1PN& gm, int q_fold = 0, double T_rpo = 0.0) {
    double E = hamiltonian_1pn(s0, m, beta);
    auto Lv = angmom_1pn(s0);
    double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"physics\":\"1pn_adm\",\"beta_inv_c2\":" << beta
       << ",\"beta_eff_v2c2\":" << (beta * gm.v2_rms)
       << ",\"beta_max_v2c2\":" << (beta * gm.v2_max)
       << ",\"dim\":" << dim << ",\"subspace\":\"" << label << "\""
       << ",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "]"
       << ",\"z0_canonical\":[";
    if (dim == 2) {
        const int xyidx[12] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16};
        for (int i=0;i<12;++i) os << s0[xyidx[i]] << (i<11?",":"");
    } else {
        for (int i=0;i<18;++i) os << s0[i] << (i<17?",":"");
    }
    os << "],\"T\":" << T;
    if (q_fold > 0) os << ",\"T_rpo\":" << T_rpo << ",\"q_fold\":" << q_fold;
    os << ",\"precess_rad\":" << two_alpha << ",\"E\":" << E
       << ",\"L\":" << (dim == 2 ? Lv[2] : Lmag)
       << ",\"ret_err\":" << ret_err
       << ",\"symp_err\":" << st.symp_err
       << ",\"k1\":" << st.k1 << ",\"k2\":" << st.k2;
    if (dim == 3) {
        os << ",\"k3\":" << st.k3 << ",\"stab_margin\":" << st.margin_3d
           << ",\"thickness_3d\":" << gm.thick;
    } else {
        os << ",\"kz\":" << st.k3 << ",\"stable_3d\":" << (st.stable_3d ? "true" : "false")
           << ",\"stab_margin\":" << st.margin_2d;
    }
    os << ",\"survived_100T\":" << st.surv100 << ",\"drift_100T\":" << st.drift100
       << ",\"rmin\":" << gm.rmin << ",\"rmax\":" << gm.rmax
       << ",\"partner_swaps\":" << gm.swaps
       << ",\"democratic\":" << (gm.swaps >= 2 ? "true" : "false") << "}\n";
    std::cout << os.str() << std::flush;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  search_1pn --verify-physics\n"
                  << "  search_1pn --trace-col <m3> <px> <py> <T> <dbeta> <nsteps> <label> [--json]\n"
                  << "  search_1pn --snap-col <m3> <px> <py> <T> <beta0> <num_2a> <den_2a> <q_fold> <label>\n"
                  << "  search_1pn --trace-2d <m1> <m2> <m3> <x2> <py1> <py2> <Th> <dbeta> <nsteps> <label> [--json]\n"
                  << "  search_1pn --snap-2d <m1> <m2> <m3> <x2> <py1> <py2> <Th> <beta0> <num_2a> <den_2a> <q_fold> <label>\n"
                  << "  search_1pn --trace-3d <m1> <m2> <m3> <x2> <py1> <py2> <pz2> <Th> <alpha> <dbeta> <nsteps> <label> [--json]\n"
                  << "  search_1pn --snap-3d <m1> <m2> <m3> <x2> <py1> <py2> <pz2> <Th> <beta> <num_2a> <den_2a> <q_fold> <label>\n";
        return 1;
    }
    std::string mode = argv[1];
    std::cout << std::setprecision(12);
    const double uz[3] = {0, 0, 1};

    if (mode == "--verify-physics") {
        return verify_physics() ? 0 : 2;
    } else if (mode == "--trace-col") {
        double m3 = std::atof(argv[2]);
        double p[6] = {std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]), 0.0, 0.0, 0.0};
        double dbeta = std::atof(argv[6]); int nsteps = std::atoi(argv[7]);
        std::string label = (argc > 8) ? argv[8] : "1pn_collinear";
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        Masses m = {1.0, 1.0, m3};
        for (int s = 0; s <= nsteps; ++s) {
            double beta = s * dbeta;
            double ferr = 0;
            if (!refine_col_1pn(p, m3, beta, ferr)) break;
            double T = p[2], two_a = p[3];
            State1PN s0 = make_col_1pn(p[0], p[1], p[4], p[5], m3), sT = s0;
            integrate_1pn(sT, T, m, beta, 1e-13, 400000);
            double ret = state_dist(rotate_state(sT, uz, -two_a), s0);
            auto gm = measure_geom(s0, T, m, beta);
            auto st = analyze_stab_1pn(s0, T, uz, two_a, m, beta, 2, emit_j);
            if (emit_j) {
                if (st.stable_2d && st.surv100 >= 100 && beta > 0.0 && ret < 1e-7)
                    emit_1pn_json(2, label, m, beta, s0, T, two_a, ret, st, gm);
                std::cerr << std::setprecision(10) << label << " beta=" << beta
                          << " beff=" << beta*gm.v2_rms << " T=" << T << " 2a=" << two_a
                          << " ret=" << ret << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                          << " stab=" << st.stable_2d << " marg=" << st.margin_2d << " surv=" << st.surv100 << "\n";
            } else {
                std::cout << beta << " " << (beta*gm.v2_rms) << " " << T << " " << two_a << " "
                          << st.k1 << " " << st.k2 << " " << st.k3 << " "
                          << (st.stable_2d?1:0) << " " << (st.stable_3d?1:0) << " " << st.margin_2d << " "
                          << p[0] << " " << p[1] << " " << p[4] << " " << p[5] << " " << gm.rmin << " " << ret << "\n" << std::flush;
            }
        }
    } else if (mode == "--snap-col") {
        double m3 = std::atof(argv[2]);
        double p[6] = {std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6]), 0.0, 0.0};
        double num = std::atof(argv[7]), den = std::atof(argv[8]);
        int q_fold = std::atoi(argv[9]);
        std::string label = (argc > 10) ? argv[10] : ("1pn_col_inert_" + std::to_string(q_fold) + "fold");
        double target_2a = 2.0 * M_PI * num / den;
        double ferr = 0;
        if (!snap_col_1pn(p, target_2a, m3, ferr)) {
            std::cerr << "Snap 1PN col failed ferr=" << ferr << "\n"; return 2;
        }
        Masses m = {1.0, 1.0, m3};
        double beta = p[3], T_rpo = p[2], T_inert = q_fold * T_rpo;
        State1PN s0 = make_col_1pn(p[0], p[1], p[4], p[5], m3), yq = s0;
        integrate_1pn(yq, T_inert, m, beta, 1e-13, 2500000);
        double ret_inert = state_dist(yq, s0);
        auto gm = measure_geom(s0, T_inert, m, beta);
        auto st = analyze_stab_1pn(s0, T_inert, uz, 0.0, m, beta, 2, true);
        std::cerr << std::setprecision(13) << "# SNAPPED 1PN COL " << label
                  << " beta=" << beta << " beff=" << beta*gm.v2_rms
                  << " 2a=" << num << "*2pi/" << den << " T_rpo=" << T_rpo << " T_inert=" << T_inert
                  << " ret_inert=" << ret_inert << " symp=" << st.symp_err
                  << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                  << " marg=" << st.margin_2d << " surv=" << st.surv100 << "\n";
        if (st.stable_2d && st.surv100 >= 100 && ret_inert < 1e-7) {
            emit_1pn_json(2, label, m, beta, s0, T_inert, target_2a, ret_inert, st, gm, q_fold, T_rpo);
        }
    } else if (mode == "--trace-2d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double x2 = std::atof(argv[5]);
        double p[5] = {std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), 0.0, 0.0};
        double dbeta = std::atof(argv[9]); int nsteps = std::atoi(argv[10]);
        std::string label = (argc > 11) ? argv[11] : "1pn_2d_rpo";
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int s = 0; s <= nsteps; ++s) {
            double beta = s * dbeta;
            if (s == 1) p[4] = compute_1pn_xc(x2, p[0], p[1], 0.0, m, beta);
            double ferr = 0;
            if (!refine_2d_rpo_1pn(p, x2, m, beta, ferr)) break;
            State1PN s0 = make_rpo_1pn(x2, p[0], p[1], 0.0, m, p[4]);
            double T = 2.0 * p[2], two_a = 2.0 * p[3];
            State1PN sT = s0; integrate_1pn(sT, T, m, beta, 1e-13, 500000);
            double ret = state_dist(rotate_state(sT, uz, -two_a), s0);
            auto gm = measure_geom(s0, T, m, beta);
            auto st = analyze_stab_1pn(s0, T, uz, two_a, m, beta, 2, emit_j);
            if (emit_j) {
                if (st.stable_2d && st.surv100 >= 100 && beta > 0.0 && ret < 1e-7)
                    emit_1pn_json(2, label, m, beta, s0, T, two_a, ret, st, gm);
                std::cerr << std::setprecision(10) << label << " beta=" << beta
                          << " beff=" << beta*gm.v2_rms << " T=" << T << " 2a=" << two_a
                          << " ret=" << ret << " xc=" << p[4]
                          << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                          << " marg2=" << st.margin_2d << " stab3=" << st.stable_3d
                          << " py1=" << p[0] << " py2=" << p[1] << " Th=" << p[2] << "\n";
            } else {
                std::cout << beta << " " << (beta*gm.v2_rms) << " " << T << " " << two_a << " "
                          << st.k1 << " " << st.k2 << " " << st.k3 << " "
                          << (st.stable_2d?1:0) << " " << (st.stable_3d?1:0) << " " << st.margin_2d
                          << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << ret << "\n" << std::flush;
            }
        }
    } else if (mode == "--snap-2d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double x2 = std::atof(argv[5]);
        double p[5] = {std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9]), 0.0};
        p[4] = compute_1pn_xc(x2, p[0], p[1], 0.0, m, p[3]);
        double num = std::atof(argv[10]), den = std::atof(argv[11]);
        int q_fold = std::atoi(argv[12]);
        std::string label = (argc > 13) ? argv[13] : ("1pn_2d_inert_" + std::to_string(q_fold) + "fold");
        double target_a = M_PI * num / den;
        double ferr = 0;
        if (!snap_2d_1pn(p, x2, target_a, m, ferr)) {
            std::cerr << "Snap 1PN 2D failed ferr=" << ferr << "\n"; return 2;
        }
        double beta = p[3];
        State1PN s0 = make_rpo_1pn(x2, p[0], p[1], 0.0, m, p[4]);
        double T_rpo = 2.0 * p[2], T_inert = q_fold * T_rpo, two_a = 2.0 * target_a;
        State1PN yq = s0; integrate_1pn(yq, T_inert, m, beta, 1e-13, 2500000);
        double ret_inert = state_dist(yq, s0);
        auto gm = measure_geom(s0, T_inert, m, beta);
        auto st = analyze_stab_1pn(s0, T_inert, uz, 0.0, m, beta, 2, true);
        std::cerr << std::setprecision(13) << "# SNAPPED 1PN 2D " << label
                  << " beta=" << beta << " beff=" << beta*gm.v2_rms
                  << " 2a=" << num << "*2pi/" << den << " T_rpo=" << T_rpo << " T_inert=" << T_inert
                  << " ret_inert=" << ret_inert << " symp=" << st.symp_err
                  << " k=(" << st.k1 << "," << st.k2 << ") kz=" << st.k3
                  << " marg=" << st.margin_2d << " surv=" << st.surv100 << "\n";
        if (st.stable_2d && st.surv100 >= 100 && ret_inert < 1e-7) {
            emit_1pn_json(2, label, m, beta, s0, T_inert, two_a, ret_inert, st, gm, q_fold, T_rpo);
        }
    } else if (mode == "--trace-3d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double p[7] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]),
                       std::atof(argv[9]), std::atof(argv[10]), 0.0, 0.0};
        double pz2 = std::atof(argv[8]);
        double dbeta = std::atof(argv[11]); int nsteps = std::atoi(argv[12]);
        std::string label = (argc > 13) ? argv[13] : "1pn_3d_rpo";
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int s = 0; s <= nsteps; ++s) {
            double beta = s * dbeta;
            if (s == 1) p[5] = compute_1pn_xc(p[0], p[1], p[2], pz2, m, beta);
            double ferr = 0;
            if (!refine_3d_rpo_1pn(p, pz2, m, beta, ferr)) break;
            State1PN s0 = make_rpo_1pn(p[0], p[1], p[2], pz2, m, p[5], p[6]);
            double T = 2.0 * p[3], two_a = 2.0 * p[4];
            auto Lv = angmom_1pn(s0);
            double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
            double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
            State1PN sT = s0; integrate_1pn(sT, T, m, beta, 1e-13, 600000);
            double ret = state_dist(rotate_state(sT, u, -two_a), s0);
            auto gm = measure_geom(s0, T, m, beta);
            auto st = analyze_stab_1pn(s0, T, u, two_a, m, beta, 3, emit_j);
            if (emit_j) {
                if (st.stable_3d && st.surv100 >= 100 && beta > 0.0 && gm.thick >= 0.05 && ret < 1e-7)
                    emit_1pn_json(3, label, m, beta, s0, T, two_a, ret, st, gm);
                std::cerr << std::setprecision(10) << label << " beta=" << beta
                          << " beff=" << beta*gm.v2_rms << " T=" << T << " 2a=" << two_a
                          << " ret=" << ret << " thick=" << gm.thick << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                          << ") stab=" << st.stable_3d << " marg=" << st.margin_3d
                          << " x2=" << p[0] << " py1=" << p[1] << " py2=" << p[2]
                          << " Th=" << p[3] << " xc=" << p[5] << " Ppar=" << p[6] << "\n";
            } else {
                std::cout << beta << " " << (beta*gm.v2_rms) << " " << T << " " << two_a << " "
                          << gm.thick << " " << st.k1 << " " << st.k2 << " " << st.k3 << " "
                          << (st.stable_3d?1:0) << " " << st.margin_3d << " "
                          << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4]
                          << " " << p[5] << " " << p[6] << " " << ret << "\n" << std::flush;
            }
        }
    } else if (mode == "--snap-3d") {
        Masses m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4])};
        double p[7] = {std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]),
                       std::atof(argv[8]), std::atof(argv[9]), 0.0, 0.0};
        double beta = std::atof(argv[10]);
        p[5] = compute_1pn_xc(p[0], p[1], p[2], p[3], m, beta);
        double num = std::atof(argv[11]), den = std::atof(argv[12]);
        int q_fold = std::atoi(argv[13]);
        std::string label = (argc > 14) ? argv[14] : ("1pn_3d_inert_" + std::to_string(q_fold) + "fold");
        double target_a = M_PI * num / den;
        double ferr = 0;
        if (!snap_3d_1pn(p, target_a, m, beta, ferr)) {
            std::cerr << "Snap 1PN 3D failed ferr=" << ferr << "\n"; return 2;
        }
        State1PN s0 = make_rpo_1pn(p[0], p[1], p[2], p[3], m, p[5], p[6]);
        double T_rpo = 2.0 * p[4], T_inert = q_fold * T_rpo, two_a = 2.0 * target_a;
        auto Lv = angmom_1pn(s0);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State1PN yq = s0; integrate_1pn(yq, T_inert, m, beta, 1e-13, 2500000);
        double ret_inert = state_dist(yq, s0);
        auto gm = measure_geom(s0, T_inert, m, beta);
        auto st = analyze_stab_1pn(s0, T_inert, u, 0.0, m, beta, 3, true);
        std::cerr << std::setprecision(13) << "# SNAPPED 1PN 3D " << label
                  << " beta=" << beta << " beff=" << beta*gm.v2_rms
                  << " 2a=" << num << "*2pi/" << den << " T_rpo=" << T_rpo << " T_inert=" << T_inert
                  << " ret_inert=" << ret_inert << " symp=" << st.symp_err << " thick=" << gm.thick
                  << " k=(" << st.k1 << "," << st.k2 << "," << st.k3
                  << ") marg=" << st.margin_3d << " surv=" << st.surv100 << "\n";
        if (st.stable_3d && st.surv100 >= 100 && gm.thick >= 0.05 && ret_inert < 1e-7) {
            emit_1pn_json(3, label, m, beta, s0, T_inert, two_a, ret_inert, st, gm, q_fold, T_rpo);
        }
    }
    return 0;
}
