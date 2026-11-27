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

using State4B = std::array<double, 24>;
using Masses4B = std::array<double, 4>;

static inline void deriv4b(const State4B& s, const Masses4B& m, State4B& ds, double* rmin_out = nullptr) {
    for (int i = 0; i < 12; ++i) ds[i] = s[12 + i];
    double a[12] = {0}, rmin2 = 1e300;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) rmin2 = r2;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            a[3*i]   += m[j]*dx*ir3; a[3*i+1] += m[j]*dy*ir3; a[3*i+2] += m[j]*dz*ir3;
            a[3*j]   -= m[i]*dx*ir3; a[3*j+1] -= m[i]*dy*ir3; a[3*j+2] -= m[i]*dz*ir3;
        }
    for (int i = 0; i < 12; ++i) ds[12 + i] = a[i];
    if (rmin_out) *rmin_out = std::sqrt(rmin2);
}

static inline double energy4b(const State4B& s, const Masses4B& m) {
    double K = 0, U = 0;
    for (int i = 0; i < 4; ++i)
        K += 0.5*m[i]*(s[12+3*i]*s[12+3*i] + s[13+3*i]*s[13+3*i] + s[14+3*i]*s[14+3*i]);
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            double dx = s[3*j]-s[3*i], dy = s[3*j+1]-s[3*i+1], dz = s[3*j+2]-s[3*i+2];
            U -= m[i]*m[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    return K + U;
}

static inline std::array<double, 3> angmom4b(const State4B& s, const Masses4B& m) {
    double Lx = 0, Ly = 0, Lz = 0;
    for (int i = 0; i < 4; ++i) {
        Lx += m[i]*(s[3*i+1]*s[14+3*i] - s[3*i+2]*s[13+3*i]);
        Ly += m[i]*(s[3*i+2]*s[12+3*i] - s[3*i]*s[14+3*i]);
        Lz += m[i]*(s[3*i]*s[13+3*i]   - s[3*i+1]*s[12+3*i]);
    }
    return {Lx, Ly, Lz};
}

static inline double state_dist4b(const State4B& a, const State4B& b) {
    double s = 0; for (int i = 0; i < 24; ++i) { double d = a[i]-b[i]; s += d*d; }
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

static State4B rotate_state4b(const State4B& s, const double u[3], double ang) {
    State4B out = s;
    for (int b = 0; b < 8; ++b) rot_vec(u, ang, out[3*b], out[3*b+1], out[3*b+2]);
    return out;
}

static bool dopri5_step4b(State4B& y, double& t, double& h, double t_end, const Masses4B& m,
                          double tol, double& rmin_track, double r_abort = 1e-4) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    if (t + h > t_end) h = t_end - t;
    if (h <= 0.0) return true;
    State4B k1,k2,k3,k4,k5,k6,k7,yt,yn;
    double rm = 1e300;
    for (int att=0;att<25;++att) {
        deriv4b(y, m, k1, &rm); if (rm < r_abort) return false;
        for (int i=0;i<24;++i) yt[i]=y[i]+h*a21*k1[i]; deriv4b(yt, m, k2);
        for (int i=0;i<24;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv4b(yt, m, k3);
        for (int i=0;i<24;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv4b(yt, m, k4);
        for (int i=0;i<24;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv4b(yt, m, k5);
        for (int i=0;i<24;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv4b(yt, m, k6);
        for (int i=0;i<24;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        deriv4b(yn, m, k7, &rm); if (rm < r_abort) return false;
        double err = 0;
        for (int i=0;i<24;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err/24.0);
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

static bool integrate4b(State4B& y, double T, const Masses4B& m, double tol = 1e-13,
                        int max_steps = 800000, double* rmin_out = nullptr, double r_abort = 1e-4) {
    double t = 0, h = std::min(1e-3, T*0.01), rmin_track = 1e300;
    int steps = 0;
    while (t < T - 1e-15 && steps < max_steps) {
        if (!dopri5_step4b(y, t, h, T, m, tol, rmin_track, r_abort)) return false;
        ++steps;
        for (int i=0;i<12;++i) if (std::abs(y[i]) > 150.0) return false;
    }
    if (rmin_out) *rmin_out = rmin_track;
    return t >= T - 1e-13;
}

// ============================================================================
// Exact Real Schur / Francis Double-Shift QR Eigenvalue Solver (General NxN)
// ============================================================================
static std::vector<std::complex<double>> real_eigenvalues(std::vector<std::vector<double>> H) {
    int n = (int)H.size();
    // 1. Reduce to upper Hessenberg via Householder reflections
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
        // Left multiply H[k+1..n-1][k..n-1]
        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k + 1; i < n; ++i) dot += v[i - k - 1] * H[i][j];
            for (int i = k + 1; i < n; ++i) H[i][j] -= beta * dot * v[i - k - 1];
        }
        // Right multiply H[0..n-1][k+1..n-1]
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;
            for (int j = k + 1; j < n; ++j) dot += H[i][j] * v[j - k - 1];
            for (int j = k + 1; j < n; ++j) H[i][j] -= beta * dot * v[j - k - 1];
        }
        for (int i = k + 2; i < n; ++i) H[i][k] = 0.0;
    }

    // 2. Francis double-shift QR iterations (Wilkinson/EISPACK HQR)
    std::vector<std::complex<double>> eigs;
    double norm = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = std::max(0, i - 1); j < n; ++j) norm += std::abs(H[i][j]);
    const double eps = 2.220446049250313e-16;
    int nn = n - 1;
    int iter = 0;
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
            // Look for two consecutive small subdiagonal elements
            int m = nn - 2;
            double p = 0, q = 0, r = 0;
            for (; m >= l; --m) {
                double z = H[m][m];
                r = x - z;
                double s = y - z;
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
            // Double QR step involving rows l..nn and columns m..nn
            for (int k = m; k <= nn - 1; ++k) {
                bool notlast = (k != nn - 1);
                if (k != m) {
                    p = H[k][k - 1];
                    q = H[k + 1][k - 1];
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
                    // Row modification
                    for (int j = k; j < n; ++j) {
                        double pp = H[k][j] + q * H[k + 1][j];
                        if (notlast) { pp += r * H[k + 2][j]; H[k + 2][j] -= pp * zp; }
                        H[k][j] -= pp * xp;
                        H[k + 1][j] -= pp * yp;
                    }
                    // Column modification
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

struct Stab4B {
    bool stable;
    double max_rad_err;
    double trivial_err;
    double margin;
    std::vector<double> k_indices;
    int surv100;
    double drift100;
};

static Stab4B analyze_stab_4b(const State4B& y0, double T, const double u[3], double two_alpha,
                              const Masses4B& m, int dim, bool run_100T = true) {
    Stab4B st{}; st.margin = -1.0; st.max_rad_err = 999.0; st.trivial_err = 999.0;
    const double eps = 1e-7;
    int N = (dim == 3) ? 24 : 16;
    const int xyidx[16] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16, 18,19, 21,22};
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));

    for (int c = 0; c < N; ++c) {
        int col = (dim == 3) ? c : xyidx[c];
        State4B yp = y0, ym = y0; yp[col] += eps; ym[col] -= eps;
        if (!integrate4b(yp, T, m, 1e-13, 800000) || !integrate4b(ym, T, m, 1e-13, 800000)) return st;
        if (std::abs(two_alpha) > 1e-14) {
            yp = rotate_state4b(yp, u, -two_alpha);
            ym = rotate_state4b(ym, u, -two_alpha);
        }
        for (int r = 0; r < N; ++r) {
            int row = (dim == 3) ? r : xyidx[r];
            M[r][c] = (yp[row] - ym[row]) / (2.0 * eps);
        }
    }

    auto eigs = real_eigenvalues(M);
    if ((int)eigs.size() != N) return st;

    // Build expected trivial eigenvalues and match greedily
    std::vector<std::complex<double>> expected;
    std::complex<double> rot_p(std::cos(two_alpha), std::sin(two_alpha));
    std::complex<double> rot_m(std::cos(two_alpha), -std::sin(two_alpha));
    if (dim == 2) {
        for (int i=0;i<4;++i) expected.emplace_back(1.0, 0.0);
        for (int i=0;i<2;++i) { expected.push_back(rot_p); expected.push_back(rot_m); }
    } else {
        for (int i=0;i<6;++i) expected.emplace_back(1.0, 0.0);
        for (int i=0;i<3;++i) { expected.push_back(rot_p); expected.push_back(rot_m); }
    }

    std::vector<bool> used(N, false);
    double max_triv = 0.0;
    for (const auto& targ : expected) {
        int best = -1; double bd = 1e9;
        for (int i = 0; i < N; ++i) {
            if (used[i]) continue;
            double d = std::abs(eigs[i] - targ);
            if (d < bd) { bd = d; best = i; }
        }
        if (best >= 0) { used[best] = true; max_triv = std::max(max_triv, bd); }
    }
    st.trivial_err = max_triv;

    // Collect remaining non-trivial eigenvalues
    std::vector<std::complex<double>> rem;
    for (int i = 0; i < N; ++i) if (!used[i]) rem.push_back(eigs[i]);

    double max_re = 0.0;
    for (const auto& z : rem) max_re = std::max(max_re, std::abs(std::abs(z) - 1.0));
    st.max_rad_err = max_re;

    // Extract stability indices k_j = 2 * Re(lambda) for conjugate pairs on the unit circle
    std::vector<double> ks;
    std::vector<bool> paired(rem.size(), false);
    for (size_t i = 0; i < rem.size(); ++i) {
        if (paired[i]) continue;
        paired[i] = true;
        // Find its conjugate or reciprocal partner
        int best = -1; double bd = 1e9;
        for (size_t j = i + 1; j < rem.size(); ++j) {
            if (paired[j]) continue;
            double d = std::min(std::abs(rem[j] - std::conj(rem[i])),
                                std::abs(rem[j] - 1.0 / rem[i]));
            if (d < bd) { bd = d; best = (int)j; }
        }
        if (best >= 0) {
            paired[best] = true;
            ks.push_back((rem[i] + rem[best]).real());
        } else {
            ks.push_back(2.0 * rem[i].real());
        }
    }
    std::sort(ks.begin(), ks.end(), [](double a, double b){ return std::abs(a) > std::abs(b); });
    st.k_indices = ks;

    double max_k = 0.0;
    for (double k : ks) max_k = std::max(max_k, std::abs(k));
    st.margin = (max_re < 2e-3) ? (2.0 - max_k) : -max_re;
    st.stable = (max_re < 2e-3 && max_triv < 5e-3 && max_k < 2.0 - 1e-6);

    if (st.stable && run_100T) {
        State4B yp = y0;
        double kick = 1e-5;
        // Momentum-preserving kick across bodies 0,1 and 2,3
        yp[12] += kick; yp[15] -= kick * m[0] / m[1];
        yp[13] += kick; yp[16] -= kick * m[0] / m[1];
        if (dim == 3) {
            yp[14] += kick; yp[17] -= kick * m[0] / m[1];
            yp[20] += 0.5*kick; yp[23] -= 0.5*kick * m[2] / m[3];
        }
        double md = 0.0; int surv = 0;
        for (int k = 0; k < 100; ++k) {
            if (!integrate4b(yp, T, m, 1e-11, 500000)) break;
            State4B yr = (std::abs(two_alpha) > 1e-14) ? rotate_state4b(yp, u, -two_alpha*(k+1)) : yp;
            double d = state_dist4b(yr, y0);
            if (d > md) md = d; if (d > 5.0) break;
            surv = k + 1;
        }
        st.surv100 = surv; st.drift100 = md;
    }
    return st;
}

struct Geom4B { double rmin, rmax, thick; int swaps; };
static Geom4B measure_geom4b(const State4B& y0, double T, const Masses4B& m) {
    const int N = 260;
    double dt = T / N, C[3][3] = {{0}}, rmin = 1e300, rmax = 0;
    int prev_cl = -1, swaps = 0, cnt = 0;
    State4B yc = y0;
    for (int k = 0; k <= N; ++k) {
        double best_r = 1e300; int cl = -1, pid = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j, ++pid) {
                double dx = yc[3*j]-yc[3*i], dy = yc[3*j+1]-yc[3*i+1], dz = yc[3*j+2]-yc[3*i+2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (r < rmin) rmin = r;
                if (r > rmax) rmax = r;
                if (r < best_r) { best_r = r; cl = pid; }
            }
        if (prev_cl != -1 && cl != prev_cl) swaps++;
        prev_cl = cl;
        for (int b = 0; b < 4; ++b) {
            double qv[3] = {yc[3*b], yc[3*b+1], yc[3*b+2]};
            for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] += qv[r]*qv[c];
            cnt++;
        }
        if (k < N) integrate4b(yc, dt, m, 1e-11, 80000);
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) C[r][c] /= cnt;
    return {rmin, rmax, std::sqrt(min_eig_3x3(C)), swaps};
}

static void emit_4b_json(int dim, const std::string& subspace, const Masses4B& m,
                         const State4B& y0, double T, double two_alpha, double ret_err,
                         const Stab4B& st, const Geom4B& gm, int q_fold = 0, double T_rpo = 0.0) {
    double E = energy4b(y0, m);
    auto Lv = angmom4b(y0, m);
    double Lmag = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
    std::ostringstream os;
    os << std::setprecision(13);
    os << "{\"n_bodies\":4,\"dim\":" << dim << ",\"subspace\":\"" << subspace << "\""
       << ",\"masses\":[" << m[0] << "," << m[1] << "," << m[2] << "," << m[3] << "],\"y0\":[";
    if (dim == 2) {
        const int xyidx[16] = {0,1, 3,4, 6,7, 9,10, 12,13, 15,16, 18,19, 21,22};
        for (int i=0;i<16;++i) os << y0[xyidx[i]] << (i<15?",":"");
    } else {
        for (int i=0;i<24;++i) os << y0[i] << (i<23?",":"");
    }
    os << "],\"T\":" << T;
    if (q_fold > 0) os << ",\"T_rpo\":" << T_rpo << ",\"q_fold\":" << q_fold;
    os << ",\"precess_rad\":" << two_alpha << ",\"E\":" << E
       << ",\"L\":" << (dim == 2 ? Lv[2] : Lmag);
    if (dim == 3) os << ",\"Lvec\":[" << Lv[0] << "," << Lv[1] << "," << Lv[2] << "]";
    os << ",\"ret_err\":" << ret_err << ",\"k_indices\":[";
    for (size_t i=0;i<st.k_indices.size();++i) os << st.k_indices[i] << (i+1<st.k_indices.size()?",":"");
    os << "],\"max_eig_rad_err\":" << st.max_rad_err
       << ",\"trivial_eig_err\":" << st.trivial_err
       << ",\"stab_margin\":" << st.margin
       << ",\"survived_100T\":" << st.surv100 << ",\"drift_100T\":" << st.drift100
       << ",\"rmin\":" << gm.rmin << ",\"rmax\":" << gm.rmax
       << ",\"thickness_3d\":" << gm.thick << ",\"partner_swaps\":" << gm.swaps << "}\n";
    std::cout << os.str() << std::flush;
}

// ============================================================================
// Subspace 1: Generalized D2d / S4 3D Hip-Hop
// Pair A (0,2) mass 1.0 at (+/-1, 0, 0), v=(0, +/-va, +vza)
// Pair B (1,3) mass mb  at (0, +/-rb, 0), v=(-/+vb, 0, -vza/mb)
// ============================================================================
static State4B make_hiphop(double mb, double rb, double va, double vb, double vza) {
    double vzb = vza / mb;
    return {
         1.0, 0.0, 0.0,   0.0,  rb, 0.0,  -1.0, 0.0, 0.0,   0.0, -rb, 0.0,
         0.0,  va, vza,   -vb, 0.0,-vzb,   0.0, -va, vza,    vb, 0.0,-vzb
    };
}

// Refine Hip-Hop at fixed (mb, vza, va_or_rb):
// Note: in D2d Hip-Hop, there is a 2-parameter family at fixed mb: (va, vza)!
// Wait: fix (mb, va, vza) -> 4 unknowns p = {rb, vb, Th, alpha} for 5 eqs (y0', z0', vx0', x1', vy1') = 0!
// Let's allow 5 unknowns p = {rb, va, vb, Th, alpha} at fixed (mb, vza) with 5 equations!
static bool refine_hiphop(double p[5], double mb, double vza, double& ferr) {
    Masses4B m = {1.0, mb, 1.0, mb};
    const double uz[3] = {0, 0, 1};
    double lambda = 1e-4;
    auto eval = [&](const double q[5], double F[5]) -> bool {
        if (q[0] < 0.1 || q[0] > 8.0 || q[3] < 0.2 || q[3] > 30.0) return false;
        State4B y = make_hiphop(mb, q[0], q[1], q[2], vza);
        if (!integrate4b(y, q[3], m, 1e-12, 200000, nullptr, 1e-3)) return false;
        State4B yt = rotate_state4b(y, uz, -q[4]);
        F[0] = yt[1];  // y0' = 0
        F[1] = yt[2];  // z0' = 0
        F[2] = yt[12]; // vx0' = 0
        F[3] = yt[3];  // x1' = 0
        F[4] = yt[16]; // vy1' = 0
        return true;
    };
    double F[5];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<5;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = std::max(1e-8, std::min(1e-6, std::sqrt(cost)*1e-2));
        double J[5][5], Fp[5], Fm[5];
        for (int col=0;col<5;++col) {
            double pp[5], pm[5]; for(int i=0;i<5;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<5;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[5]={0};
        if (!solve_lm<5,5>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<5;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<5;++i) d[i]*=(0.25/sn);
        double np[5], nF[5]; for(int i=0;i<5;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<5;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<5;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.25); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// ============================================================================
// Subspace 2: 3+1 Hybrid 4-Body 3D Perp-RPO (x1=-1, vz1=0)
// All 4 bodies start on x-axis with vx=0, CM=0, P=0:
// Unknowns (8): p = {x2, x4, vy1, vy2, vy4, vz2, Th, alpha} at fixed (m, vz4)
// Residuals (9, 8 independent via Lx=0): y0',z0', y1',z1', y2',z2', vx0',vx1',vx2'
// ============================================================================
static State4B make_perp_4b(double x2, double x4, double vy1, double vy2, double vy4,
                            double vz2, double vz4, const Masses4B& m) {
    double x1 = -1.0;
    double x3 = -(m[0]*x1 + m[1]*x2 + m[3]*x4) / m[2];
    double vy3 = -(m[0]*vy1 + m[1]*vy2 + m[3]*vy4) / m[2];
    double vz3 = -(m[1]*vz2 + m[3]*vz4) / m[2];
    return {
        x1, 0, 0,   x2, 0, 0,   x3, 0, 0,   x4, 0, 0,
        0, vy1, 0,  0, vy2, vz2, 0, vy3, vz3, 0, vy4, vz4
    };
}

static bool refine_perp4b_3d(double p[8], double vz4, const Masses4B& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[8], double F[11]) -> bool {
        if (q[6] < 0.5 || q[6] > 80.0) return false;
        State4B y0 = make_perp_4b(q[0], q[1], q[2], q[3], q[4], q[5], vz4, m);
        auto Lv = angmom4b(y0, m);
        double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
        if (Lm < 1e-4) return false;
        double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
        State4B y = y0;
        if (!integrate4b(y, q[6], m, 1e-13, 600000)) return false;
        State4B yt = rotate_state4b(y, u, -q[7]);
        // Include body 4 (index 3) explicitly so Jacobian w.r.t. (x4, vy4) is O(1) even when m4 << 1!
        F[0]=yt[1];  F[1]=yt[2];   // y0', z0'
        F[2]=yt[4];  F[3]=yt[5];   // y1', z1'
        F[4]=yt[7];  F[5]=yt[8];   // y2', z2'
        F[6]=yt[10]; F[7]=yt[11];  // y3', z3'
        F[8]=yt[12]; F[9]=yt[15]; F[10]=yt[21]; // vx0', vx1', vx3'
        return true;
    };
    double F[11];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<11;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = 1e-7, J[11][8], Fp[11], Fm[11];
        for (int col=0;col<8;++col) {
            double pp[8], pm[8]; for(int i=0;i<8;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<11;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[8]={0};
        if (!solve_lm<11,8>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<8;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<8;++i) d[i]*=(0.25/sn);
        double np[8], nF[11]; for(int i=0;i<8;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<11;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<8;++i) p[i]=np[i]; for(int i=0;i<11;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// 2D 3+1 Perp (vz2=vz4=0, alpha=0 fixed inertial closure):
// 6 vars p = {x2, x4, vy1, vy2, vy4, Th} for 8 eqs {y0, y1, y2, y3, vx0, vx1, vx2, vx3} (6 independent!)
static bool refine_perp4b_2d(double p[6], const Masses4B& m, double& ferr) {
    double lambda = 1e-5;
    auto eval = [&](const double q[6], double F[8]) -> bool {
        if (q[5] < 0.5 || q[5] > 80.0) return false;
        State4B y = make_perp_4b(q[0], q[1], q[2], q[3], q[4], 0.0, 0.0, m);
        if (!integrate4b(y, q[5], m, 1e-13, 600000)) return false;
        F[0]=y[1]; F[1]=y[4]; F[2]=y[7]; F[3]=y[10];
        F[4]=y[12]; F[5]=y[15]; F[6]=y[18]; F[7]=y[21];
        return true;
    };
    double F[8];
    if (!eval(p, F)) return false;
    double cost = 0; for(int i=0;i<8;++i) cost += F[i]*F[i];
    for (int it=0;it<35;++it) {
        if (std::sqrt(cost) < 1e-11) break;
        double eps = 1e-7, J[8][6], Fp[8], Fm[8];
        for (int col=0;col<6;++col) {
            double pp[6], pm[6]; for(int i=0;i<6;++i){pp[i]=p[i];pm[i]=p[i];}
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<8;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[6]={0};
        if (!solve_lm<8,6>(J, F, lambda, d)) { lambda*=10; continue; }
        double sn=0; for(int i=0;i<6;++i) sn+=d[i]*d[i]; sn=std::sqrt(sn);
        if (sn > 0.25) for(int i=0;i<6;++i) d[i]*=(0.25/sn);
        double np[6], nF[8]; for(int i=0;i<6;++i) np[i]=p[i]+d[i];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc=0; for(int i=0;i<8;++i) nc+=nF[i]*nF[i];
        if (nc < cost) { for(int i=0;i<6;++i) p[i]=np[i]; for(int i=0;i<8;++i) F[i]=nF[i]; cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

// Central-star / axial-body 1+3 (body 3 = heavy central mass or light axial Sitnikov body)
// 3 equal outer bodies (0,1,2) in equilateral / symmetric 3D ring + body 3 on z-axis:
// r_k = (R cos(2pi k/3), R sin(2pi k/3), z_ring), v_k = (-vt sin(2pi k/3), vt cos(2pi k/3), vz_ring)
// body 3 at (0, 0, -3*m_ring/m4 * z_ring), v3 = (0, 0, -3*m_ring/m4 * vz_ring)
// At t=0: z_ring = 0, vz_ring = vz0 (so body 3 starts at origin with vz3 = -3/m4 * vz0!)
// Periodicity condition at t=Th: z_ring(Th) = 0 and dr_ring/dt(Th) = 0! (Exact 2D reduced system!)
static State4B make_axial_3plus1(double m4, double vt, double vz_ring) {
    State4B s{};
    for (int k = 0; k < 3; ++k) {
        double th = 2.0 * M_PI * k / 3.0;
        double c = std::cos(th), sn = std::sin(th);
        s[3*k]   = c;  s[3*k+1] = sn; s[3*k+2] = 0.0;
        s[12+3*k] = -vt*sn; s[13+3*k] = vt*c; s[14+3*k] = vz_ring;
    }
    s[9] = 0; s[10] = 0; s[11] = 0;
    s[21] = 0; s[22] = 0; s[23] = -3.0 * vz_ring / m4;
    return s;
}

static bool refine_axial_3plus1(double p[3], double m4, double vz_ring, double& ferr) {
    // p = {vt, Th, alpha}: at t=Th, z0=0, radial velocity vr0=0, and ring rotated by alpha!
    Masses4B m = {1.0, 1.0, 1.0, m4};
    const double uz[3] = {0, 0, 1};
    double lambda = 1e-4;
    auto eval = [&](const double q[3], double F[3]) -> bool {
        if (q[1] < 0.2 || q[1] > 40.0) return false;
        State4B y = make_axial_3plus1(m4, q[0], vz_ring);
        if (!integrate4b(y, q[1], m, 1e-13, 300000)) return false;
        State4B yt = rotate_state4b(y, uz, -q[2]);
        F[0] = yt[1];  // y0' = 0
        F[1] = yt[2];  // z0' = 0
        F[2] = yt[12]; // vx0' = 0 (vr=0!)
        return true;
    };
    double F[3];
    if (!eval(p, F)) return false;
    double cost = F[0]*F[0]+F[1]*F[1]+F[2]*F[2];
    for (int it=0;it<30;++it) {
        if (std::sqrt(cost) < 1e-12) break;
        double eps = 1e-7, J[3][3], Fp[3], Fm[3];
        for (int col=0;col<3;++col) {
            double pp[3]={p[0],p[1],p[2]}, pm[3]={p[0],p[1],p[2]};
            pp[col]+=eps; pm[col]-=eps;
            if (!eval(pp, Fp) || !eval(pm, Fm)) return false;
            for (int r=0;r<3;++r) J[r][col] = (Fp[r]-Fm[r])/(2*eps);
        }
        double d[3]={0};
        if (!solve_lm<3,3>(J, F, lambda, d)) { lambda*=10; continue; }
        double np[3]={p[0]+d[0], p[1]+d[1], p[2]+d[2]}, nF[3];
        if (!eval(np, nF)) { lambda*=10; continue; }
        double nc = nF[0]*nF[0]+nF[1]*nF[1]+nF[2]*nF[2];
        if (nc < cost) { for(int i=0;i<3;++i){p[i]=np[i];F[i]=nF[i];} cost=nc; lambda=std::max(1e-13,lambda*0.2); }
        else { lambda*=5; if (lambda>1e7) break; }
    }
    ferr = std::sqrt(cost);
    return ferr < 1e-9;
}

static bool verify_eig() {
    // Test Francis QR eigensolver on a known 6x6 block-diagonal matrix with real & complex eigenvalues
    std::vector<std::vector<double>> A(6, std::vector<double>(6, 0.0));
    A[0][0] =  0.6; A[0][1] = -0.8; A[1][0] =  0.8; A[1][1] =  0.6; // e^{+-i atan2(4,3)}
    A[2][2] = -0.5; A[2][3] = -0.8660254037844386; A[3][2] = 0.8660254037844386; A[3][3] = -0.5;
    A[4][4] =  2.0; A[5][5] =  0.5;
    // Similarity transform by a non-trivial invertible matrix S
    double S[6][6] = {
        {1,0.3,0.1,0,0.2,0}, {0,1,0.4,0.1,0,0.2}, {0.2,0,1,0.3,0.1,0},
        {0,0.2,0,1,0.3,0.1}, {0.1,0,0.2,0,1,0.4}, {0,0.1,0,0.2,0,1}
    };
    // Build B = S * A * S^{-1} using Gauss-Jordan for S^{-1}
    double aug[6][12] = {{0}};
    for (int i=0;i<6;++i) { for(int j=0;j<6;++j) aug[i][j]=S[i][j]; aug[i][6+i]=1.0; }
    for (int col=0;col<6;++col) {
        double pv = aug[col][col];
        for (int c=0;c<12;++c) aug[col][c] /= pv;
        for (int r=0;r<6;++r) if (r!=col) { double f=aug[r][col]; for(int c=0;c<12;++c) aug[r][c]-=f*aug[col][c]; }
    }
    std::vector<std::vector<double>> B(6, std::vector<double>(6, 0.0)), tmp(6, std::vector<double>(6, 0.0));
    for (int i=0;i<6;++i) for(int k=0;k<6;++k) for(int j=0;j<6;++j) tmp[i][j] += S[i][k]*A[k][j];
    for (int i=0;i<6;++i) for(int k=0;k<6;++k) for(int j=0;j<6;++j) B[i][j] += tmp[i][k]*aug[k][6+j];
    auto eigs = real_eigenvalues(B);
    std::vector<std::complex<double>> exact = {
        {0.6, 0.8}, {0.6, -0.8}, {-0.5, 0.8660254037844386}, {-0.5, -0.8660254037844386}, {2.0, 0.0}, {0.5, 0.0}
    };
    double max_err = 0.0;
    for (auto ex : exact) {
        double bd = 1e9;
        for (auto z : eigs) bd = std::min(bd, std::abs(z - ex));
        max_err = std::max(max_err, bd);
    }
    std::cout << "Francis QR Eigensolver Self-Test max err = " << max_err << "\n";
    return max_err < 1e-12;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  search_4body --verify-eig\n"
                  << "  search_4body --hiphop <mb> <rb0> <va0> <vb0> <vza0> <Th0> <a0> <dvza> <nsteps> [--json]\n"
                  << "  search_4body --hiphop-mb <mb0> <rb0> <va0> <vb0> <vza> <Th0> <a0> <dmb> <nsteps> [--json]\n"
                  << "  search_4body --axial-3plus1 <m4> <vz0> <dvz> <nsteps> [--json]\n"
                  << "  search_4body --hybrid-2d <m1> <m2> <m3> <x2> <vy1> <vy2> <Th> <x4> <vy4> <m4_0> <dm4> <nsteps> <label> [--json]\n"
                  << "  search_4body --hybrid-3d <m1> <m2> <m3> <x2> <vy1> <vy2> <vz2> <Th> <alpha> <x4> <vy4> <vz4> <m4_0> <dm4> <nsteps> <label> [--json]\n";
        return 1;
    }
    std::string mode = argv[1];
    std::cout << std::setprecision(12);
    const double uz[3] = {0, 0, 1};

    if (mode == "--verify-eig") {
        return verify_eig() ? 0 : 2;
    } else if (mode == "--hiphop" || mode == "--hiphop-mb") {
        bool vary_mb = (mode == "--hiphop-mb");
        double mb0 = std::atof(argv[2]);
        double p[5] = {std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]),
                       std::atof(argv[7]), std::atof(argv[8])};
        double vza0 = std::atof(argv[6]);
        double step = std::atof(argv[9]);
        int nsteps = std::atoi(argv[10]);
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        for (int s = 0; s <= nsteps; ++s) {
            double mb = vary_mb ? (mb0 + s * step) : mb0;
            double vza = vary_mb ? vza0 : (vza0 + s * step);
            if (mb <= 0.002 || vza < 1e-4) break;
            double ferr = 0;
            if (!refine_hiphop(p, mb, vza, ferr)) break;
            Masses4B m = {1.0, mb, 1.0, mb};
            State4B y0 = make_hiphop(mb, p[0], p[1], p[2], vza);
            double T = 2.0 * p[3], two_a = 2.0 * p[4];
            State4B yT = y0; integrate4b(yT, T, m, 1e-13, 500000);
            double ret = state_dist4b(rotate_state4b(yT, uz, -two_a), y0);
            auto gm = measure_geom4b(y0, T, m);
            auto st = analyze_stab_4b(y0, T, uz, two_a, m, 3, emit_j);
            std::cerr << std::setprecision(10) << "hh mb=" << mb << " vza=" << vza << " rb=" << p[0]
                      << " va=" << p[1] << " vb=" << p[2] << " T=" << T << " 2a=" << two_a
                      << " thick=" << gm.thick << " rmin=" << gm.rmin << " rad_err=" << st.max_rad_err
                      << " triv=" << st.trivial_err << " stab=" << st.stable << " marg=" << st.margin
                      << " surv=" << st.surv100 << " ret=" << ret << "\n";
            if (emit_j && st.stable && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7) {
                emit_4b_json(3, "4b_3d_hiphop", m, y0, T, two_a, ret, st, gm);
            } else if (!emit_j) {
                std::cout << mb << " " << vza << " " << p[0] << " " << p[1] << " " << p[2] << " "
                          << T << " " << two_a << " " << gm.thick << " " << gm.rmin << " "
                          << (st.stable?1:0) << " " << st.margin << " " << st.max_rad_err << "\n" << std::flush;
            }
        }
    } else if (mode == "--axial-3plus1") {
        double m4 = std::atof(argv[2]), vz0 = std::atof(argv[3]), dvz = std::atof(argv[4]);
        int nsteps = std::atoi(argv[5]);
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        // Circular equilibrium tangential speed of 3 equal masses on unit circle around central mass m4:
        // Centripetal accel = vt^2 / 1 = m4 + 1/sqrt(3)
        double vt0 = std::sqrt(m4 + 1.0 / std::sqrt(3.0));
        // Vertical small-oscillation frequency of ring vs central body:
        // d^2(z_ring - z_4)/dt^2 = -(m4 + 3) * (z_ring - z_4) / 1^3 => omega_z = sqrt(m4 + 3), Th = pi / omega_z
        double omega_z = std::sqrt(m4 + 3.0);
        double Th0 = M_PI / omega_z;
        double a0 = vt0 * Th0;
        double p[3] = {vt0, Th0, a0};
        Masses4B m = {1.0, 1.0, 1.0, m4};
        for (int s = 0; s <= nsteps; ++s) {
            double vz = vz0 + s * dvz;
            if (vz < 1e-4) break;
            double ferr = 0;
            if (!refine_axial_3plus1(p, m4, vz, ferr)) break;
            State4B y0 = make_axial_3plus1(m4, p[0], vz);
            double T = 2.0 * p[1], two_a = 2.0 * p[2];
            State4B yT = y0; integrate4b(yT, T, m, 1e-13, 500000);
            double ret = state_dist4b(rotate_state4b(yT, uz, -two_a), y0);
            auto gm = measure_geom4b(y0, T, m);
            auto st = analyze_stab_4b(y0, T, uz, two_a, m, 3, emit_j);
            std::cerr << std::setprecision(10) << "axial m4=" << m4 << " vz=" << vz
                      << " vt=" << p[0] << " T=" << T << " 2a=" << two_a
                      << " thick=" << gm.thick << " rmin=" << gm.rmin << " rad_err=" << st.max_rad_err
                      << " triv=" << st.trivial_err << " stab=" << st.stable << " marg=" << st.margin
                      << " surv=" << st.surv100 << " ret=" << ret << "\n";
            if (emit_j && st.stable && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7) {
                emit_4b_json(3, "4b_3d_central_axial", m, y0, T, two_a, ret, st, gm);
            } else if (!emit_j) {
                std::cout << m4 << " " << vz << " " << p[0] << " " << T << " " << two_a << " "
                          << gm.thick << " " << gm.rmin << " " << (st.stable?1:0) << " "
                          << st.margin << " " << st.max_rad_err << "\n" << std::flush;
            }
        }
    } else if (mode == "--hybrid-2d") {
        Masses4B m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[11])};
        double p[6] = {std::atof(argv[5]), std::atof(argv[9]), std::atof(argv[6]),
                       std::atof(argv[7]), std::atof(argv[10]), std::atof(argv[8])};
        double dm4 = std::atof(argv[12]);
        int nsteps = std::atoi(argv[13]);
        std::string label = (argc > 14) ? argv[14] : "4b_2d_3plus1";
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        double m4_0 = m[3];
        for (int s = 0; s <= nsteps; ++s) {
            m[3] = m4_0 + s * dm4;
            if (m[3] <= 1e-6) break;
            double ferr = 0;
            if (!refine_perp4b_2d(p, m, ferr)) break;
            State4B y0 = make_perp_4b(p[0], p[1], p[2], p[3], p[4], 0.0, 0.0, m);
            double T = 2.0 * p[5];
            State4B yT = y0; integrate4b(yT, T, m, 1e-13, 600000);
            double ret = state_dist4b(yT, y0);
            auto gm = measure_geom4b(y0, T, m);
            auto st = analyze_stab_4b(y0, T, uz, 0.0, m, 2, emit_j);
            std::cerr << std::setprecision(10) << label << " m4=" << m[3] << " x4=" << p[1]
                      << " vy4=" << p[4] << " T=" << T << " rmin=" << gm.rmin << " swaps=" << gm.swaps
                      << " rad_err=" << st.max_rad_err << " triv=" << st.trivial_err
                      << " stab=" << st.stable << " marg=" << st.margin << " surv=" << st.surv100
                      << " ret=" << ret << "\n";
            if (emit_j && st.stable && st.surv100 >= 100 && ret < 1e-7) {
                emit_4b_json(2, label, m, y0, T, 0.0, ret, st, gm);
            } else if (!emit_j) {
                std::cout << m[3] << " " << p[1] << " " << p[4] << " " << T << " " << gm.rmin
                          << " " << (st.stable?1:0) << " " << st.margin << " " << st.max_rad_err << "\n" << std::flush;
            }
        }
    } else if (mode == "--hybrid-3d") {
        // Args: m1 m2 m3 x2 vy1 vy2 vz2 Th alpha x4 vy4 vz4 m4_0 dm4 nsteps label [--json]
        Masses4B m = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[14])};
        double p[8] = {std::atof(argv[5]), std::atof(argv[11]), std::atof(argv[6]),
                       std::atof(argv[7]), std::atof(argv[12]), std::atof(argv[8]),
                       std::atof(argv[9]), std::atof(argv[10])};
        double vz4 = std::atof(argv[13]);
        double dm4 = std::atof(argv[15]);
        int nsteps = std::atoi(argv[16]);
        std::string label = (argc > 17) ? argv[17] : "4b_3d_3plus1";
        bool emit_j = false; for (int i=1;i<argc;++i) if (std::string(argv[i]) == "--json") emit_j = true;
        double m4_0 = m[3];
        for (int s = 0; s <= nsteps; ++s) {
            m[3] = m4_0 + s * dm4;
            if (m[3] <= 1e-6) break;
            double ferr = 0;
            if (!refine_perp4b_3d(p, vz4, m, ferr)) break;
            State4B y0 = make_perp_4b(p[0], p[1], p[2], p[3], p[4], p[5], vz4, m);
            double T = 2.0 * p[6], two_a = 2.0 * p[7];
            auto Lv = angmom4b(y0, m);
            double Lm = std::sqrt(Lv[0]*Lv[0]+Lv[1]*Lv[1]+Lv[2]*Lv[2]);
            double u[3] = {Lv[0]/Lm, Lv[1]/Lm, Lv[2]/Lm};
            State4B yT = y0; integrate4b(yT, T, m, 1e-13, 800000);
            double ret = state_dist4b(rotate_state4b(yT, u, -two_a), y0);
            auto gm = measure_geom4b(y0, T, m);
            auto st = analyze_stab_4b(y0, T, u, two_a, m, 3, emit_j);
            std::cerr << std::setprecision(8) << label << " m4=" << m[3] << " x4=" << p[1]
                      << " vy4=" << p[4] << " vz4=" << vz4 << " T=" << T << " 2a=" << two_a
                      << " thick=" << gm.thick << " rmin=" << gm.rmin << " swaps=" << gm.swaps
                      << " rad_err=" << st.max_rad_err << " triv=" << st.trivial_err
                      << " stab=" << st.stable << " marg=" << st.margin << " surv=" << st.surv100
                      << " ret=" << ret << " ks=[";
            for (size_t i=0;i<st.k_indices.size();++i) std::cerr << st.k_indices[i] << (i+1<st.k_indices.size()?",":"");
            std::cerr << "]\n";
            if (emit_j && st.stable && st.surv100 >= 100 && gm.thick >= 0.05 && ret < 1e-7) {
                emit_4b_json(3, label, m, y0, T, two_a, ret, st, gm);
            } else if (!emit_j) {
                std::cout << m[3] << " " << p[1] << " " << p[4] << " " << T << " " << two_a << " "
                          << gm.thick << " " << gm.rmin << " " << (st.stable?1:0) << " "
                          << st.margin << " " << st.max_rad_err << "\n" << std::flush;
            }
        }
    }
    return 0;
}
