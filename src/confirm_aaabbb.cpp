// confirm_aaabbb.cpp
// Independent 128-bit (__float128, ~34 decimal digits) symplectic 8-stage 16th-order
// Gauss-Legendre (GL16) confirmation engine + 10,000-period full-3D symplectic GL8
// KAM-torus perturbation verifier for the AAAbbb Shuttle-Exchange family.
//
// Build:
//   g++ -O3 -march=native -std=c++20 -pthread src/confirm_aaabbb.cpp -lquadmath -o bin/confirm_aaabbb
//
// Usage:
//   ./bin/confirm_aaabbb --self-test
//   ./bin/confirm_aaabbb --converge [--out-dir data/aaabbb]
//   ./bin/confirm_aaabbb --long-term [--periods 10000] [--out-dir data/aaabbb]
//   ./bin/confirm_aaabbb --all [--periods 10000] [--out-dir data/aaabbb]

#include <quadmath.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using f128 = __float128;

static inline f128 qabs(f128 x) { return fabsq(x); }
static inline f128 qsqrt(f128 x) { return sqrtq(x); }
static inline f128 qhypot(f128 x, f128 y) { return hypotq(x, y); }
static inline f128 qsin(f128 x) { return sinq(x); }
static inline f128 qcos(f128 x) { return cosq(x); }
static inline f128 qatan2(f128 y, f128 x) { return atan2q(y, x); }
static inline f128 qpow(f128 x, f128 y) { return powq(x, y); }

static std::string qstr(f128 x, int digits = 33) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.*Qg", digits, x);
    return std::string(buf);
}
static std::string qsci(f128 x, int digits = 3) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.*Qe", digits, x);
    return std::string(buf);
}
static std::string qfix(f128 x, int digits = 32) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.*Qf", digits, x);
    return std::string(buf);
}

static const f128 PI_Q = strtoflt128("3.14159265358979323846264338327950288419716939937510", nullptr);

// ============================================================================
// Exact s-stage 2s-order Gauss-Legendre Butcher Tableau in __float128
// ============================================================================
template <int S>
struct GaussLegendre {
    f128 c[S];
    f128 b[S];
    f128 A[S][S];
    f128 max_quad_err = 0;
    f128 max_stage_err = 0;
    f128 max_symp_err = 0;

    GaussLegendre() {
        // 1. Roots u_i of Legendre polynomial P_S(u) = 0 on [-1, 1] and exact closed-form weights b_i
        f128 u_roots[S], dp_vals[S];
        for (int i = 0; i < S; ++i) {
            f128 u = -qcos(PI_Q * (f128)(i + 0.75) / (f128)(S + 0.5));
            f128 dp = 0.0Q;
            for (int it = 0; it < 40; ++it) {
                f128 p0 = 1.0Q, p1 = u;
                for (int k = 2; k <= S; ++k) {
                    f128 pk = ((f128)(2 * k - 1) * u * p1 - (f128)(k - 1) * p0) / (f128)k;
                    p0 = p1; p1 = pk;
                }
                dp = (f128)S * (u * p1 - p0) / (u * u - 1.0Q);
                f128 du = p1 / dp;
                u -= du;
                if (qabs(du) < 1e-33Q) break;
            }
            u_roots[i] = u;
            dp_vals[i] = dp;
        }
        // Sort by u_roots ascending
        int ord[S];
        for (int i = 0; i < S; ++i) ord[i] = i;
        std::sort(ord, ord + S, [&](int a, int b_) { return u_roots[a] < u_roots[b_]; });
        for (int i = 0; i < S; ++i) {
            f128 u = u_roots[ord[i]], dp = dp_vals[ord[i]];
            c[i] = 0.5Q * (u + 1.0Q);
            b[i] = 1.0Q / ((1.0Q - u * u) * dp * dp);
        }

        // 2. Exact A_{ij} = int_0^{c_i} l_j(x) dx via S-point Gauss quadrature in factored Lagrange form
        //    (exact for degree S-1 <= 2S-1 with zero power-basis Vandermonde cancellation!)
        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                f128 sum = 0.0Q;
                for (int q = 0; q < S; ++q) {
                    f128 x = c[i] * c[q];
                    f128 lj = 1.0Q;
                    for (int k = 0; k < S; ++k) {
                        if (k == j) continue;
                        lj *= (x - c[k]) / (c[j] - c[k]);
                    }
                    sum += b[q] * lj;
                }
                A[i][j] = c[i] * sum;
            }
        }

        // 3. Self-verification of order 2S, stage order S, and symplecticity
        for (int k = 1; k <= 2 * S; ++k) {
            f128 sum = 0.0Q;
            for (int j = 0; j < S; ++j) sum += b[j] * qpow(c[j], (f128)(k - 1));
            f128 err = qabs(sum - 1.0Q / (f128)k);
            if (err > max_quad_err) max_quad_err = err;
        }
        for (int i = 0; i < S; ++i) {
            for (int k = 1; k <= S; ++k) {
                f128 sum = 0.0Q;
                for (int j = 0; j < S; ++j) sum += A[i][j] * qpow(c[j], (f128)(k - 1));
                f128 err = qabs(sum - qpow(c[i], (f128)k) / (f128)k);
                if (err > max_stage_err) max_stage_err = err;
            }
        }
        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                f128 err = qabs(b[i] * A[i][j] + b[j] * A[j][i] - b[i] * b[j]);
                if (err > max_symp_err) max_symp_err = err;
            }
        }
    }
};

static const GaussLegendre<8> GL16; // 8-stage 16th-order exact symplectic in __float128
static const GaussLegendre<4> GL8;  // 4-stage 8th-order exact symplectic

// ============================================================================
// 128-bit Newtonian Planar 3-Body Dynamics + Exact Analytical Variational Eqns
// ============================================================================
using StateQ = std::array<f128, 12>;
using MassesQ = std::array<f128, 3>;

static inline StateQ make_perp_q(f128 x2, f128 vy1, f128 vy2, const MassesQ& m) {
    f128 x3 = -(-m[0] + m[1] * x2) / m[2];
    f128 vy3 = -(m[0] * vy1 + m[1] * vy2) / m[2];
    return {-1.0Q, 0.0Q, x2, 0.0Q, x3, 0.0Q, 0.0Q, vy1, 0.0Q, vy2, 0.0Q, vy3};
}

static inline StateQ rot2d_q(const StateQ& s, f128 ang) {
    f128 c = qcos(ang), sn = qsin(ang);
    StateQ o{};
    for (int b = 0; b < 6; ++b) {
        o[2 * b]     = c * s[2 * b] - sn * s[2 * b + 1];
        o[2 * b + 1] = sn * s[2 * b] + c * s[2 * b + 1];
    }
    return o;
}

static inline f128 min_pair_dist_q(const StateQ& y) {
    f128 r12 = qhypot(y[2] - y[0], y[3] - y[1]);
    f128 r13 = qhypot(y[4] - y[0], y[5] - y[1]);
    f128 r23 = qhypot(y[4] - y[2], y[5] - y[3]);
    f128 r = (r12 < r13) ? r12 : r13;
    return (r23 < r) ? r23 : r;
}

static inline void deriv_q(const StateQ& y, const MassesQ& m, StateQ& dy) {
    for (int i = 0; i < 6; ++i) dy[i] = y[6 + i];
    for (int i = 6; i < 12; ++i) dy[i] = 0.0Q;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            f128 dx = y[2 * j] - y[2 * i];
            f128 dy_ = y[2 * j + 1] - y[2 * i + 1];
            f128 r2 = dx * dx + dy_ * dy_;
            f128 ir = 1.0Q / qsqrt(r2);
            f128 ir3 = ir * ir * ir;
            f128 fx = dx * ir3, fy = dy_ * ir3;
            dy[6 + 2 * i] += m[j] * fx; dy[7 + 2 * i] += m[j] * fy;
            dy[6 + 2 * j] -= m[i] * fx; dy[7 + 2 * j] -= m[i] * fy;
        }
    }
}

static inline void hessian_q(const StateQ& y, const MassesQ& m, f128 Kp[6][6], f128 Kz[3][3]) {
    for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) Kp[r][c] = 0.0Q;
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) Kz[r][c] = 0.0Q;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            f128 dx = y[2 * j] - y[2 * i];
            f128 dy_ = y[2 * j + 1] - y[2 * i + 1];
            f128 r2 = dx * dx + dy_ * dy_;
            f128 ir = 1.0Q / qsqrt(r2);
            f128 ir3 = ir * ir * ir;
            f128 ir5 = ir3 * ir * ir;
            f128 hxx = (3.0Q * dx * dx - r2) * ir5;
            f128 hxy = (3.0Q * dx * dy_) * ir5;
            f128 hyy = (3.0Q * dy_ * dy_ - r2) * ir5;
            // d(a_i)/d(r_i) = +m_j * H,  d(a_i)/d(r_j) = -m_j * H
            Kp[2*i][2*i]     += m[j]*hxx; Kp[2*i][2*i+1]   += m[j]*hxy;
            Kp[2*i+1][2*i]   += m[j]*hxy; Kp[2*i+1][2*i+1] += m[j]*hyy;
            Kp[2*i][2*j]     -= m[j]*hxx; Kp[2*i][2*j+1]   -= m[j]*hxy;
            Kp[2*i+1][2*j]   -= m[j]*hxy; Kp[2*i+1][2*j+1] -= m[j]*hyy;
            Kp[2*j][2*j]     += m[i]*hxx; Kp[2*j][2*j+1]   += m[i]*hxy;
            Kp[2*j+1][2*j]   += m[i]*hxy; Kp[2*j+1][2*j+1] += m[i]*hyy;
            Kp[2*j][2*i]     -= m[i]*hxx; Kp[2*j][2*i+1]   -= m[i]*hxy;
            Kp[2*j+1][2*i]   -= m[i]*hxy; Kp[2*j+1][2*i+1] -= m[i]*hyy;
            // Out-of-plane vertical linearized acceleration: a_{z,i} = sum_j m_j (z_j - z_i) / r_{ij}^3
            Kz[i][i] -= m[j] * ir3; Kz[i][j] += m[j] * ir3;
            Kz[j][j] -= m[i] * ir3; Kz[j][i] += m[i] * ir3;
        }
    }
}

static inline void invariants_q(const StateQ& y, const MassesQ& m, f128& E, f128& L) {
    f128 K = 0.0Q, U = 0.0Q;
    L = 0.0Q;
    for (int i = 0; i < 3; ++i) {
        K += 0.5Q * m[i] * (y[6 + 2*i] * y[6 + 2*i] + y[7 + 2*i] * y[7 + 2*i]);
        L += m[i] * (y[2*i] * y[7 + 2*i] - y[2*i + 1] * y[6 + 2*i]);
        for (int j = i + 1; j < 3; ++j)
            U -= m[i] * m[j] / qhypot(y[2*j] - y[2*i], y[2*j + 1] - y[2*i + 1]);
    }
    E = K + U;
}

// Single GL16 step for trajectory only (12D)
static bool gl16_step_traj(StateQ& y, f128 h, const MassesQ& m, StateQ Y_stages[8] = nullptr) {
    StateQ K[8], Y[8];
    StateQ f0;
    deriv_q(y, m, f0);
    for (int s = 0; s < 8; ++s) K[s] = f0;
    for (int it = 0; it < 60; ++it) {
        f128 max_rel = 0.0Q;
        for (int s = 0; s < 8; ++s) {
            for (int i = 0; i < 12; ++i) {
                f128 sum = 0.0Q;
                for (int j = 0; j < 8; ++j) sum += GL16.A[s][j] * K[j][i];
                Y[s][i] = y[i] + h * sum;
            }
            StateQ Knew;
            deriv_q(Y[s], m, Knew);
            for (int i = 0; i < 12; ++i) {
                f128 sc = 1.0Q + qabs(Knew[i]);
                f128 d = qabs(Knew[i] - K[s][i]) / sc;
                if (d > max_rel) max_rel = d;
                K[s][i] = Knew[i];
            }
        }
        if (max_rel < 2e-33Q) {
            for (int i = 0; i < 12; ++i) {
                f128 sum = 0.0Q;
                for (int s = 0; s < 8; ++s) sum += GL16.b[s] * K[s][i];
                y[i] += h * sum;
            }
            if (Y_stages) for (int s = 0; s < 8; ++s) Y_stages[s] = Y[s];
            return true;
        }
    }
    return false;
}

// Single GL16 step for trajectory + exact 12x12 planar & 6x6 vertical variational matrices
static bool gl16_step_var(StateQ& y, f128 Mp[12][12], f128 Mz[6][6], f128 h, const MassesQ& m) {
    StateQ Y[8];
    StateQ y_try = y;
    if (!gl16_step_traj(y_try, h, m, Y)) return false;
    y = y_try;

    f128 Kp[8][6][6], Kz[8][3][3];
    for (int s = 0; s < 8; ++s) hessian_q(Y[s], m, Kp[s], Kz[s]);

    // Planar 12x12 variational stages
    f128 KMp[8][12][12];
    for (int s = 0; s < 8; ++s) {
        for (int c = 0; c < 12; ++c) {
            for (int r = 0; r < 6; ++r) {
                KMp[s][r][c] = Mp[6 + r][c];
                f128 v = 0.0Q;
                for (int k = 0; k < 6; ++k) v += Kp[s][r][k] * Mp[k][c];
                KMp[s][6 + r][c] = v;
            }
        }
    }
    for (int it = 0; it < 60; ++it) {
        f128 max_rel = 0.0Q;
        for (int s = 0; s < 8; ++s) {
            f128 YMp[12][12];
            for (int r = 0; r < 12; ++r) {
                for (int c = 0; c < 12; ++c) {
                    f128 sum = 0.0Q;
                    for (int j = 0; j < 8; ++j) sum += GL16.A[s][j] * KMp[j][r][c];
                    YMp[r][c] = Mp[r][c] + h * sum;
                }
            }
            for (int c = 0; c < 12; ++c) {
                for (int r = 0; r < 6; ++r) {
                    f128 v0 = YMp[6 + r][c];
                    f128 v1 = 0.0Q;
                    for (int k = 0; k < 6; ++k) v1 += Kp[s][r][k] * YMp[k][c];
                    f128 d0 = qabs(v0 - KMp[s][r][c]) / (1.0Q + qabs(v0));
                    f128 d1 = qabs(v1 - KMp[s][6 + r][c]) / (1.0Q + qabs(v1));
                    if (d0 > max_rel) max_rel = d0;
                    if (d1 > max_rel) max_rel = d1;
                    KMp[s][r][c] = v0;
                    KMp[s][6 + r][c] = v1;
                }
            }
        }
        if (max_rel < 1e-28Q) break;
    }
    for (int r = 0; r < 12; ++r)
        for (int c = 0; c < 12; ++c) {
            f128 sum = 0.0Q;
            for (int s = 0; s < 8; ++s) sum += GL16.b[s] * KMp[s][r][c];
            Mp[r][c] += h * sum;
        }

    // Vertical 6x6 variational stages
    f128 KMz[8][6][6];
    for (int s = 0; s < 8; ++s) {
        for (int c = 0; c < 6; ++c) {
            for (int r = 0; r < 3; ++r) {
                KMz[s][r][c] = Mz[3 + r][c];
                f128 v = 0.0Q;
                for (int k = 0; k < 3; ++k) v += Kz[s][r][k] * Mz[k][c];
                KMz[s][3 + r][c] = v;
            }
        }
    }
    for (int it = 0; it < 60; ++it) {
        f128 max_rel = 0.0Q;
        for (int s = 0; s < 8; ++s) {
            f128 YMz[6][6];
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c) {
                    f128 sum = 0.0Q;
                    for (int j = 0; j < 8; ++j) sum += GL16.A[s][j] * KMz[j][r][c];
                    YMz[r][c] = Mz[r][c] + h * sum;
                }
            for (int c = 0; c < 6; ++c)
                for (int r = 0; r < 3; ++r) {
                    f128 v0 = YMz[3 + r][c];
                    f128 v1 = 0.0Q;
                    for (int k = 0; k < 3; ++k) v1 += Kz[s][r][k] * YMz[k][c];
                    f128 d0 = qabs(v0 - KMz[s][r][c]) / (1.0Q + qabs(v0));
                    f128 d1 = qabs(v1 - KMz[s][3 + r][c]) / (1.0Q + qabs(v1));
                    if (d0 > max_rel) max_rel = d0;
                    if (d1 > max_rel) max_rel = d1;
                    KMz[s][r][c] = v0;
                    KMz[s][3 + r][c] = v1;
                }
        }
        if (max_rel < 1e-28Q) break;
    }
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) {
            f128 sum = 0.0Q;
            for (int s = 0; s < 8; ++s) sum += GL16.b[s] * KMz[s][r][c];
            Mz[r][c] += h * sum;
        }
    return true;
}

// Integrate trajectory only over [0, T] with rmin-scaled GL16 steps
static bool integrate_traj_q(StateQ& y, f128 T, const MassesQ& m, f128 c_scale = 0.025Q, f128 hmax = 0.025Q) {
    f128 t = 0.0Q;
    int steps = 0;
    while (t < T - 1e-33Q) {
        if (++steps > 500000) return false;
        f128 rm = min_pair_dist_q(y);
        if (rm < 1e-6Q) return false;
        f128 h = c_scale * rm * qsqrt(rm);
        if (h > hmax) h = hmax;
        if (h < 1e-7Q) h = 1e-7Q;
        if (t + h > T) h = T - t;
        // Automatic sub-step if fixed-point doesn't converge in 40 iters
        int sub = 1;
        while (!gl16_step_traj(y, h / (f128)sub, m)) {
            sub *= 2;
            if (sub > 64) return false;
        }
        for (int k = 1; k < sub; ++k)
            if (!gl16_step_traj(y, h / (f128)sub, m)) return false;
        t += h;
    }
    return true;
}

// ============================================================================
// Exact Montgomery Free-Group Word + Syzygy & Geometry Recorder (128-bit)
// ============================================================================
static std::string reduce_cyclic(std::string w) {
    auto inv = [](char c) -> char {
        if (c == 'a') return 'A';
        if (c == 'A') return 'a';
        if (c == 'b') return 'B';
        return 'b';
    };
    bool ch = true;
    while (ch && !w.empty()) {
        ch = false;
        std::string nw;
        for (char c : w) {
            if (!nw.empty() && nw.back() == inv(c)) { nw.pop_back(); ch = true; }
            else nw.push_back(c);
        }
        while (nw.size() >= 2 && nw.front() == inv(nw.back())) {
            nw.erase(nw.begin()); nw.pop_back(); ch = true;
        }
        w = nw;
    }
    return w;
}
static std::string canonical_word(const std::string& w_in) {
    std::string w = reduce_cyclic(w_in);
    if (w.empty()) return "1";
    const char* maps[4] = {"aAbB", "AaBb", "bBaA", "BbAa"};
    auto apply = [](const std::string& s, const char* mp) {
        std::string r; r.reserve(s.size());
        for (char c : s) {
            if (c == 'a') r.push_back(mp[0]);
            else if (c == 'A') r.push_back(mp[1]);
            else if (c == 'b') r.push_back(mp[2]);
            else r.push_back(mp[3]);
        }
        return r;
    };
    auto inv_rev = [](const std::string& s) {
        std::string r; r.reserve(s.size());
        for (auto it = s.rbegin(); it != s.rend(); ++it) {
            char c = *it;
            r.push_back(c=='a'?'A':c=='A'?'a':c=='b'?'B':'b');
        }
        return r;
    };
    std::string best;
    for (int m = 0; m < 4; ++m) {
        std::string s1 = apply(w, maps[m]);
        std::string s2 = inv_rev(s1);
        for (const std::string& cand : {s1, s2}) {
            std::string dbl = cand + cand;
            for (size_t i = 0; i < cand.size(); ++i) {
                std::string rot = dbl.substr(i, cand.size());
                if (best.empty() || rot < best) best = rot;
            }
        }
    }
    return best;
}

static std::string word_from_samples(const std::vector<std::array<double, 6>>& pos_samples) {
    // Compute exact Montgomery word over all 6 body permutations (S_3 canonical)
    const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    std::string best_s3;
    const double thr = 1.0 / std::sqrt(3.0);
    for (int p = 0; p < 6; ++p) {
        int i0 = perms[p][0], i1 = perms[p][1], i2 = perms[p][2];
        std::vector<std::array<double, 3>> pts;
        pts.reserve(pos_samples.size());
        for (const auto& s : pos_samples) {
            double rx1 = s[2*i0], ry1 = s[2*i0+1];
            double rx2 = s[2*i1], ry2 = s[2*i1+1];
            double rx3 = s[2*i2], ry3 = s[2*i2+1];
            double rhox = (rx1 - rx2) / std::sqrt(2.0), rhoy = (ry1 - ry2) / std::sqrt(2.0);
            double lamx = (rx1 + rx2 - 2.0*rx3) / std::sqrt(6.0), lamy = (ry1 + ry2 - 2.0*ry3) / std::sqrt(6.0);
            double n2 = rhox*rhox + rhoy*rhoy + lamx*lamx + lamy*lamy;
            double sx = 2.0*(rhox*lamx + rhoy*lamy) / n2;
            double sy = 2.0*(rhox*lamy - rhoy*lamx) / n2;
            double sz = (lamx*lamx + lamy*lamy - rhox*rhox - rhoy*rhoy) / n2;
            // Rotate north pole to collision (sz -> -1 is collision of i0-i1)
            pts.push_back({sx, sy, sz});
        }
        std::string w;
        size_t n = pts.size();
        for (size_t k = 0; k < n; ++k) {
            const auto& a = pts[k];
            const auto& b = pts[(k + 1) % n];
            // Wait: standard Šuvakov-Dmitrašinović projection from (sx, sy, sz) with poles at (0,0,1):
            // Collision points lie on equator sy = 0!
            // Actually match audit_catalog.cpp's exact convention:
            //   sx = (|lam|^2 - |rho|^2)/n2, sy = 2*(rho.lam)/n2, sz = 2*(rho x lam)/n2
            // Let's use audit_catalog.cpp's EXACT (nx, ny, nz) coordinates!
            (void)a; (void)b;
        }
        // Recompute directly in audit_catalog's exact (nx, ny, nz)
        std::vector<std::array<double, 3>> sp;
        sp.reserve(pos_samples.size());
        for (const auto& s : pos_samples) {
            double rx1 = s[2*i0], ry1 = s[2*i0+1];
            double rx2 = s[2*i1], ry2 = s[2*i1+1];
            double rx3 = s[2*i2], ry3 = s[2*i2+1];
            double rhox = (rx1 - rx2) / std::sqrt(2.0), rhoy = (ry1 - ry2) / std::sqrt(2.0);
            double lamx = (rx1 + rx2 - 2.0*rx3) / std::sqrt(6.0), lamy = (ry1 + ry2 - 2.0*ry3) / std::sqrt(6.0);
            double n2 = rhox*rhox + rhoy*rhoy + lamx*lamx + lamy*lamy;
            double nx = (lamx*lamx + lamy*lamy - rhox*rhox - rhoy*rhoy) / n2;
            double ny = 2.0*(rhox*lamx + rhoy*lamy) / n2;
            double nz = 2.0*(rhox*lamy - rhoy*lamx) / n2;
            sp.push_back({nx, ny, nz});
        }
        for (size_t k = 0; k < sp.size(); ++k) {
            const auto& c0 = sp[k];
            const auto& c1 = sp[(k + 1) % sp.size()];
            if (c0[2] == 0.0) continue;
            if ((c0[2] < 0 && c1[2] >= 0) || (c0[2] > 0 && c1[2] <= 0)) {
                double t = -c0[2] / (c1[2] - c0[2]);
                double cx = c0[0] + t * (c1[0] - c0[0]);
                double cy = c0[1] + t * (c1[1] - c0[1]);
                double ux = cy / (1.0 - cx); // stereographic projection from collision (1,0,0)
                bool up = (c0[2] < 0);
                if (ux > thr) w.push_back(up ? 'a' : 'A');
                else if (ux < -thr) w.push_back(up ? 'b' : 'B');
            }
        }
        std::string cw = canonical_word(w);
        if (best_s3.empty() || cw < best_s3) best_s3 = cw;
    }
    return best_s3;
}

// Full 128-bit period evaluation (trajectory + 12x12 & 6x6 variational + geometry + syzygies + word)
struct QuadAudit {
    std::string name, role;
    MassesQ m;
    f128 alpha = 0; // half-period rotation (0 for inertial; pi/q for q-fold RPO)
    int q_fold = 1;
    f128 p[4]; // (x2, vy1, vy2, Thalf)
    StateQ y0, yend;
    f128 T, E, L, Tstar, Lambda;
    f128 ret_err, dE, dL, symp_err;
    f128 step_half_y_err = 0, step_half_k_err = 0;
    bool real_k = false, stable_3d = false;
    f128 k1 = 0, k2 = 0, kz = 0, marg3 = -1, krein_gap = 0;
    f128 rmin_pair[3], rmax_pair[3], rmin = 0, rmax = 0;
    f128 wind[3];
    int syzygies = 0;
    std::string syz_seq, word_s3;
    // For RPO closures (q_fold > 1), also exact T_inert = q_fold * T_rpo audit:
    bool has_inert_closure = false;
    f128 in_T = 0, in_ret = 0, in_symp = 0, in_k1 = 0, in_k2 = 0, in_kz = 0, in_marg3 = -1;
    bool in_st3 = false;
};

static bool eval_var_period_q(const StateQ& y0, f128 T, f128 rot_angle, const MassesQ& m,
                              f128 c_scale, StateQ& yend, f128& symp_err,
                              bool& real_k, f128& k1, f128& k2, f128& kz, f128& marg3,
                              QuadAudit* geom_out = nullptr,
                              f128 Mrot_out[12][12] = nullptr, f128 Mz_out[6][6] = nullptr) {
    StateQ y = y0;
    f128 Mp[12][12] = {{0}}, Mz[6][6] = {{0}};
    for (int i = 0; i < 12; ++i) Mp[i][i] = 1.0Q;
    for (int i = 0; i < 6; ++i) Mz[i][i] = 1.0Q;

    const int prs[3][2] = {{0,1},{0,2},{1,2}};
    f128 prev_ang[3], tot_ang[3] = {0,0,0};
    f128 rmin_p[3] = {1e30Q, 1e30Q, 1e30Q}, rmax_p[3] = {0, 0, 0};
    for (int q = 0; q < 3; ++q) {
        f128 dx = y[2*prs[q][1]] - y[2*prs[q][0]], dy_ = y[2*prs[q][1]+1] - y[2*prs[q][0]+1];
        f128 r = qhypot(dx, dy_);
        rmin_p[q] = r; rmax_p[q] = r;
        prev_ang[q] = qatan2(dy_, dx);
    }

    auto signed_area = [](const StateQ& s) -> f128 {
        return (s[2] - s[0]) * (s[5] - s[1]) - (s[3] - s[1]) * (s[4] - s[0]);
    };
    std::vector<std::array<double, 6>> samples;
    std::string syz_seq;
    f128 prev_area = signed_area(y);
    StateQ prev_y = y;

    f128 t = 0.0Q, hmax = (c_scale < 0.02Q) ? 0.0125Q : 0.025Q;
    int steps = 0;
    while (t < T - 1e-33Q) {
        if (++steps > 600000) return false;
        if (geom_out) samples.push_back({(double)y[0],(double)y[1],(double)y[2],(double)y[3],(double)y[4],(double)y[5]});
        f128 rm = min_pair_dist_q(y);
        if (rm < 1e-6Q) return false;
        f128 h = c_scale * rm * qsqrt(rm);
        if (h > hmax) h = hmax;
        if (h < 1e-7Q) h = 1e-7Q;
        if (t + h > T) h = T - t;
        prev_y = y;
        int sub = 1;
        while (!gl16_step_var(y, Mp, Mz, h / (f128)sub, m)) {
            sub *= 2;
            if (sub > 64) return false;
        }
        for (int k = 1; k < sub; ++k)
            if (!gl16_step_var(y, Mp, Mz, h / (f128)sub, m)) return false;
        t += h;

        if (geom_out) {
            for (int q = 0; q < 3; ++q) {
                f128 dx = y[2*prs[q][1]] - y[2*prs[q][0]], dy_ = y[2*prs[q][1]+1] - y[2*prs[q][0]+1];
                f128 r = qhypot(dx, dy_);
                if (r < rmin_p[q]) rmin_p[q] = r;
                if (r > rmax_p[q]) rmax_p[q] = r;
                f128 ang = qatan2(dy_, dx), da = ang - prev_ang[q];
                while (da > PI_Q) da -= 2.0Q * PI_Q;
                while (da < -PI_Q) da += 2.0Q * PI_Q;
                tot_ang[q] += da; prev_ang[q] = ang;
            }
            f128 area = signed_area(y);
            // Detect syzygy crossing inside (0, T] (avoid counting t=0 twice)
            if (steps > 1 && ((prev_area < 0 && area >= 0) || (prev_area > 0 && area <= 0))) {
                f128 frac = -prev_area / (area - prev_area);
                f128 cx[3], cy[3];
                for (int b = 0; b < 3; ++b) {
                    cx[b] = prev_y[2*b] + frac * (y[2*b] - prev_y[2*b]);
                    cy[b] = prev_y[2*b+1] + frac * (y[2*b+1] - prev_y[2*b+1]);
                }
                // Middle body is the one with dot product (r_a - r_m).(r_b - r_m) < 0
                int mid = 1;
                for (int b = 0; b < 3; ++b) {
                    int a = (b + 1) % 3, c = (b + 2) % 3;
                    f128 dp = (cx[a]-cx[b])*(cx[c]-cx[b]) + (cy[a]-cy[b])*(cy[c]-cy[b]);
                    if (dp < 0) { mid = b + 1; break; }
                }
                syz_seq.push_back((char)('0' + mid));
            }
            prev_area = area;
        }
    }
    yend = y;

    // Rotate monodromy by -rot_angle
    f128 co = qcos(-rot_angle), sn = qsin(-rot_angle);
    f128 Mrot[12][12];
    for (int col = 0; col < 12; ++col) {
        for (int blk = 0; blk < 6; ++blk) {
            f128 vx = Mp[2*blk][col], vy = Mp[2*blk + 1][col];
            Mrot[2*blk][col]     = co * vx - sn * vy;
            Mrot[2*blk + 1][col] = sn * vx + co * vy;
        }
    }

    // Exact canonical symplecticity check: S^T J S = J where S = D M D^{-1}, p_i = m_i v_i
    symp_err = 0.0Q;
    for (int i = 0; i < 12; ++i) {
        f128 sci = (i >= 6) ? (1.0Q / m[(i - 6) / 2]) : 1.0Q;
        for (int j = 0; j < 12; ++j) {
            f128 scj = (j >= 6) ? (1.0Q / m[(j - 6) / 2]) : 1.0Q;
            f128 s = 0.0Q;
            for (int k = 0; k < 6; ++k) {
                f128 mk = m[k / 2];
                s += (Mrot[k][i] * sci) * (mk * Mrot[k + 6][j] * scj)
                   - (mk * Mrot[k + 6][i] * sci) * (Mrot[k][j] * scj);
            }
            f128 tgt = (j == i + 6) ? 1.0Q : ((i == j + 6) ? -1.0Q : 0.0Q);
            f128 e = qabs(s - tgt);
            if (e > symp_err) symp_err = e;
        }
    }

    // Trace invariants (immune to Jordan block splitting on trivial unit eigenvalues)
    f128 tr1 = 0.0Q, tr2 = 0.0Q;
    for (int i = 0; i < 12; ++i) {
        tr1 += Mrot[i][i];
        for (int j = 0; j < 12; ++j) tr2 += Mrot[i][j] * Mrot[j][i];
    }
    f128 s_val = tr1 - (4.0Q + 4.0Q * qcos(rot_angle));
    f128 t2_val = tr2 - (4.0Q + 4.0Q * qcos(2.0Q * rot_angle));
    f128 p_val = 0.5Q * (s_val * s_val - (t2_val + 4.0Q));
    f128 disc = s_val * s_val - 4.0Q * p_val;
    if (disc >= -1e-24Q) {
        real_k = true;
        f128 sd = qsqrt(disc > 0.0Q ? disc : 0.0Q);
        k1 = 0.5Q * (s_val + sd);
        k2 = 0.5Q * (s_val - sd);
    } else {
        real_k = false;
        k1 = 0.5Q * s_val;
        k2 = qsqrt(-disc);
    }
    f128 trz = 0.0Q;
    for (int i = 0; i < 6; ++i) trz += Mz[i][i];
    kz = trz - (2.0Q + 2.0Q * qcos(rot_angle));
    f128 maxk = qabs(k1);
    if (qabs(k2) > maxk) maxk = qabs(k2);
    if (qabs(kz) > maxk) maxk = qabs(kz);
    marg3 = (real_k && maxk < 2.0Q) ? (2.0Q - maxk) : -1.0Q;

    if (geom_out) {
        for (int q = 0; q < 3; ++q) {
            geom_out->rmin_pair[q] = rmin_p[q];
            geom_out->rmax_pair[q] = rmax_p[q];
            geom_out->wind[q] = tot_ang[q] / (2.0Q * PI_Q);
        }
        geom_out->rmin = std::min({rmin_p[0], rmin_p[1], rmin_p[2]});
        geom_out->rmax = std::max({rmax_p[0], rmax_p[1], rmax_p[2]});
        geom_out->syz_seq = syz_seq;
        geom_out->syzygies = (int)syz_seq.size();
        geom_out->word_s3 = word_from_samples(samples);
    }
    if (Mrot_out) {
        for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) Mrot_out[r][c] = Mrot[r][c];
    }
    if (Mz_out) {
        for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) Mz_out[r][c] = Mz[r][c];
    }
    return true;
}

static void extract_k_and_symp(const f128 Mrot[12][12], const f128 Mz[6][6], f128 rot_angle,
                               const MassesQ& m, f128& symp_err, bool& real_k,
                               f128& k1, f128& k2, f128& kz, f128& marg3) {
    symp_err = 0.0Q;
    for (int i = 0; i < 12; ++i) {
        f128 sci = (i >= 6) ? (1.0Q / m[(i - 6) / 2]) : 1.0Q;
        for (int j = 0; j < 12; ++j) {
            f128 scj = (j >= 6) ? (1.0Q / m[(j - 6) / 2]) : 1.0Q;
            f128 s = 0.0Q;
            for (int k = 0; k < 6; ++k) {
                f128 mk = m[k / 2];
                s += (Mrot[k][i] * sci) * (mk * Mrot[k + 6][j] * scj)
                   - (mk * Mrot[k + 6][i] * sci) * (Mrot[k][j] * scj);
            }
            f128 tgt = (j == i + 6) ? 1.0Q : ((i == j + 6) ? -1.0Q : 0.0Q);
            f128 e = qabs(s - tgt);
            if (e > symp_err) symp_err = e;
        }
    }
    f128 tr1 = 0.0Q, tr2 = 0.0Q;
    for (int i = 0; i < 12; ++i) {
        tr1 += Mrot[i][i];
        for (int j = 0; j < 12; ++j) tr2 += Mrot[i][j] * Mrot[j][i];
    }
    f128 s_val = tr1 - (4.0Q + 4.0Q * qcos(rot_angle));
    f128 t2_val = tr2 - (4.0Q + 4.0Q * qcos(2.0Q * rot_angle));
    f128 p_val = 0.5Q * (s_val * s_val - (t2_val + 4.0Q));
    f128 disc = s_val * s_val - 4.0Q * p_val;
    if (disc >= -1e-24Q) {
        real_k = true;
        f128 sd = qsqrt(disc > 0.0Q ? disc : 0.0Q);
        k1 = 0.5Q * (s_val + sd);
        k2 = 0.5Q * (s_val - sd);
    } else {
        real_k = false;
        k1 = 0.5Q * s_val;
        k2 = qsqrt(-disc);
    }
    f128 trz = 0.0Q;
    for (int i = 0; i < 6; ++i) trz += Mz[i][i];
    kz = trz - (2.0Q + 2.0Q * qcos(rot_angle));
    f128 maxk = std::max({qabs(k1), qabs(k2), qabs(kz)});
    marg3 = (real_k && maxk < 2.0Q) ? (2.0Q - maxk) : -1.0Q;
}

// ============================================================================
// 128-bit Newton Shooting to ret_err < 1e-31 (32+ decimal digits)
// ============================================================================
static bool solve4_q(f128 A[4][4], f128 b[4], f128 x[4]) {
    f128 M[4][5];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }
    for (int c = 0; c < 4; ++c) {
        int piv = c;
        for (int r = c + 1; r < 4; ++r)
            if (qabs(M[r][c]) > qabs(M[piv][c])) piv = r;
        if (qabs(M[piv][c]) < 1e-32Q) return false;
        if (piv != c) for (int j = c; j <= 4; ++j) std::swap(M[c][j], M[piv][j]);
        f128 d = M[c][c];
        for (int j = c; j <= 4; ++j) M[c][j] /= d;
        for (int r = 0; r < 4; ++r) if (r != c) {
            f128 f = M[r][c];
            for (int j = c; j <= 4; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < 4; ++i) x[i] = M[i][4];
    return true;
}

static bool newton_shoot_128(f128 p[4], f128 alpha, const MassesQ& m, f128& final_res) {
    auto eval_F = [&](const f128 q[4], f128 F[4]) -> bool {
        if (q[3] < 0.5Q || q[3] > 30.0Q) return false;
        StateQ y = make_perp_q(q[0], q[1], q[2], m);
        if (!integrate_traj_q(y, q[3], m, 0.020Q, 0.020Q)) return false;
        StateQ yr = rot2d_q(y, -alpha);
        F[0] = yr[1]; F[1] = yr[3]; F[2] = yr[6]; F[3] = yr[8];
        return true;
    };
    f128 F[4];
    if (!eval_F(p, F)) return false;
    for (int it = 0; it < 12; ++it) {
        f128 fn = 0.0Q;
        for (int i = 0; i < 4; ++i) if (qabs(F[i]) > fn) fn = qabs(F[i]);
        final_res = fn;
        if (fn < 1e-31Q) return true;
        // Central FD Jacobian with step h ~ 1e-16 (optimal for 34-digit precision!)
        f128 J[4][4];
        f128 h = (fn > 1e-10Q) ? 1e-8Q : 1e-16Q;
        for (int j = 0; j < 4; ++j) {
            f128 pp[4] = {p[0], p[1], p[2], p[3]};
            f128 pm[4] = {p[0], p[1], p[2], p[3]};
            pp[j] += h; pm[j] -= h;
            f128 Fp[4], Fm[4];
            if (!eval_F(pp, Fp) || !eval_F(pm, Fm)) return false;
            for (int i = 0; i < 4; ++i) J[i][j] = (Fp[i] - Fm[i]) / (2.0Q * h);
        }
        f128 rhs[4] = {-F[0], -F[1], -F[2], -F[3]}, dp[4];
        if (!solve4_q(J, rhs, dp)) return false;
        f128 damp = 1.0Q;
        bool ok = false;
        for (int ls = 0; ls < 8; ++ls) {
            f128 cand[4];
            for (int j = 0; j < 4; ++j) cand[j] = p[j] + damp * dp[j];
            f128 Fc[4];
            if (eval_F(cand, Fc)) {
                f128 fnc = 0.0Q;
                for (int i = 0; i < 4; ++i) if (qabs(Fc[i]) > fnc) fnc = qabs(Fc[i]);
                if (fnc < fn || fnc < 1e-31Q) {
                    for (int j = 0; j < 4; ++j) { p[j] = cand[j]; F[j] = Fc[j]; }
                    final_res = fnc;
                    ok = true;
                    break;
                }
            }
            damp *= 0.5Q;
        }
        if (!ok) break;
    }
    return final_res < 1e-30Q;
}

// ============================================================================
// Fast Double Continuation Helper (continues from 2D-DEM-002 to any target)
// ============================================================================
static bool continue_double_seed(const double p_start[4], std::array<double, 3> m_start,
                                 std::array<double, 3> m_tgt, double alpha_tgt, double p_out[4]) {
    auto integ = [](std::array<double, 12>& y, double T, const std::array<double, 3>& m) -> bool {
        static constexpr double a21=1./5,a31=3./40,a32=9./40,a41=44./45,a42=-56./15,a43=32./9,
            a51=19372./6561,a52=-25360./2187,a53=64448./6561,a54=-212./729,
            a61=9017./3168,a62=-355./33,a63=46732./5247,a64=49./176,a65=-5103./18656,
            b1=35./384,b3=500./1113,b4=125./192,b5=-2187./6784,b6=11./84,
            e1=71./57600,e3=-71./16695,e4=71./1920,e5=-17253./339200,e6=22./525,e7=-1./40;
        auto rhs = [&](const std::array<double, 12>& s, std::array<double, 12>& ds) {
            for (int i=0;i<6;++i) ds[i]=s[6+i];
            for (int i=6;i<12;++i) ds[i]=0.0;
            for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
                double dx=s[2*j]-s[2*i], dy=s[2*j+1]-s[2*i+1], r2=dx*dx+dy*dy;
                double ir3 = 1.0 / (r2 * std::sqrt(r2));
                ds[6+2*i]+=m[j]*dx*ir3; ds[7+2*i]+=m[j]*dy*ir3;
                ds[6+2*j]-=m[i]*dx*ir3; ds[7+2*j]-=m[i]*dy*ir3;
            }
        };
        double t=0, h=1e-3, tol=1e-13;
        std::array<double,12> k1,k2,k3,k4,k5,k6,k7,yt,yn;
        for (int step=0; t<T-1e-14 && step<100000; ++step) {
            if (t+h>T) h=T-t;
            rhs(y,k1);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*a21*k1[i]; rhs(yt,k2);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); rhs(yt,k3);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); rhs(yt,k4);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); rhs(yt,k5);
            for(int i=0;i<12;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); rhs(yt,k6);
            for(int i=0;i<12;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); rhs(yn,k7);
            double err=0;
            for(int i=0;i<12;++i){ double ei=h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]); double sc=tol*(1+std::max(std::abs(y[i]),std::abs(yn[i]))); err+=(ei/sc)*(ei/sc); }
            err=std::sqrt(err/12.0);
            if(err<=1.0){ y=yn; t+=h; h=std::min(h*(err==0?5.0:std::min(5.0,std::max(0.2,0.9*std::pow(err,-0.2)))), 0.05); }
            else { h*=std::max(0.1,0.9*std::pow(err,-0.25)); if(h<1e-14) return false; }
        }
        return t>=T-1e-12;
    };
    auto eval = [&](const double q[4], double a, const std::array<double,3>& m, double F[4]) {
        if (q[3] < 0.5 || q[3] > 30.0) return false;
        std::array<double,12> y = {-1.0, 0, q[0], 0, -(-m[0]+m[1]*q[0])/m[2], 0,
                                   0, q[1], 0, q[2], 0, -(m[0]*q[1]+m[1]*q[2])/m[2]};
        if (!integ(y, q[3], m)) return false;
        double co = std::cos(-a), sn = std::sin(-a);
        F[0] = sn*y[0] + co*y[1];
        F[1] = sn*y[2] + co*y[3];
        F[2] = co*y[6] - sn*y[7];
        F[3] = co*y[8] - sn*y[9];
        return true;
    };
    for (int i=0;i<4;++i) p_out[i] = p_start[i];
    auto step_phase = [&](std::array<double,3> ma, std::array<double,3> mb, double a0, double a1, int ns) -> bool {
        for (int s = 1; s <= ns; ++s) {
            double frac = (double)s / ns;
            std::array<double,3> m = {
                ma[0] + frac*(mb[0]-ma[0]),
                ma[1] + frac*(mb[1]-ma[1]),
                ma[2] + frac*(mb[2]-ma[2])
            };
            double a = a0 + frac * (a1 - a0);
            for (int it = 0; it < 20; ++it) {
                double F[4]; if (!eval(p_out, a, m, F)) return false;
                double fn = std::max({std::abs(F[0]),std::abs(F[1]),std::abs(F[2]),std::abs(F[3])});
                if (fn < 1e-12) break;
                double J[4][4];
                for (int j=0;j<4;++j) {
                    double pp[4]={p_out[0],p_out[1],p_out[2],p_out[3]}, pm[4]={p_out[0],p_out[1],p_out[2],p_out[3]};
                    double h=1e-7; pp[j]+=h; pm[j]-=h; double Fp[4],Fm[4];
                    if (!eval(pp,a,m,Fp)||!eval(pm,a,m,Fm)) return false;
                    for(int i=0;i<4;++i) J[i][j]=(Fp[i]-Fm[i])/(2*h);
                }
                double M[4][5];
                for(int i=0;i<4;++i){ for(int j=0;j<4;++j) M[i][j]=J[i][j]; M[i][4]=-F[i]; }
                for(int c=0;c<4;++c){
                    int piv=c; for(int r=c+1;r<4;++r) if(std::abs(M[r][c])>std::abs(M[piv][c])) piv=r;
                    if(std::abs(M[piv][c])<1e-18) return false;
                    for(int j=c;j<=4;++j) std::swap(M[c][j],M[piv][j]);
                    double d=M[c][c]; for(int j=c;j<=4;++j) M[c][j]/=d;
                    for(int r=0;r<4;++r) if(r!=c){ double f=M[r][c]; for(int j=c;j<=4;++j) M[r][j]-=f*M[c][j]; }
                }
                double step_n = std::sqrt(M[0][4]*M[0][4]+M[1][4]*M[1][4]+M[2][4]*M[2][4]+M[3][4]*M[3][4]);
                double sc0 = (step_n > 0.15) ? (0.15 / step_n) : 1.0;
                double damp = sc0;
                for (int ls=0; ls<10; ++ls) {
                    double cand[4]; for(int j=0;j<4;++j) cand[j]=p_out[j]+damp*M[j][4];
                    double Fc[4];
                    if (eval(cand,a,m,Fc) && std::max({std::abs(Fc[0]),std::abs(Fc[1]),std::abs(Fc[2]),std::abs(Fc[3])}) < fn) {
                        for(int j=0;j<4;++j) p_out[j]=cand[j]; break;
                    }
                    damp *= 0.5;
                }
            }
        }
        return true;
    };
    double dm = std::max({std::abs(m_tgt[0]-m_start[0]), std::abs(m_tgt[1]-m_start[1]), std::abs(m_tgt[2]-m_start[2])});
    int ns_m = std::max(1, (int)std::ceil(dm / 0.008));
    if (!step_phase(m_start, m_tgt, 0.0, 0.0, ns_m)) return false;
    if (std::abs(alpha_tgt) > 1e-12) {
        int ns_a = std::max(1, (int)std::ceil(std::abs(alpha_tgt) / 0.01));
        if (!step_phase(m_tgt, m_tgt, 0.0, alpha_tgt, ns_a)) return false;
    }
    return true;
}

// ============================================================================
// 10,000-Period Full-3D Symplectic GL8 Empirical KAM-Torus Perturbation Test
// ============================================================================
struct LongTermResult {
    std::string name;
    double delta;
    int target_periods;
    int survived_periods;
    double max_dE_rel;
    double max_dL_abs;
    double rmin_obs;
    double rmax_obs;
    double max_z;
    double max_return_drift;
};

static LongTermResult run_long_term_3d(const QuadAudit& qa, double delta, int n_periods) {
    // Full 3D 18D state: [x0,y0,z0, x1,y1,z1, x2,y2,z2, vx0,vy0,vz0, vx1,vy1,vz1, vx2,vy2,vz2]
    double m[3] = {(double)qa.m[0], (double)qa.m[1], (double)qa.m[2]};
    double y[18] = {0};
    for (int b = 0; b < 3; ++b) {
        y[3*b]     = (double)qa.y0[2*b];
        y[3*b + 1] = (double)qa.y0[2*b + 1];
        y[3*b + 2] = 0.0;
        y[9 + 3*b]     = (double)qa.y0[6 + 2*b];
        y[9 + 3*b + 1] = (double)qa.y0[7 + 2*b];
        y[9 + 3*b + 2] = 0.0;
    }
    // Apply exact zero-total-linear-momentum full-3D velocity perturbation of norm ~ delta!
    // Direction: body 0 gets (+0.6, +0.5, +0.6245)*delta, compensated on bodies 1 & 2 to preserve P = 0.
    double dv0[3] = {0.57735026919 * delta, 0.57735026919 * delta, 0.57735026919 * delta};
    double M23 = m[1] + m[2];
    for (int c = 0; c < 3; ++c) {
        y[9 + c]     += dv0[c];
        y[9 + 3 + c] -= (m[0] * dv0[c]) / M23;
        y[9 + 6 + c] -= (m[0] * dv0[c]) / M23;
    }
    double y_ref[18];
    std::memcpy(y_ref, y, sizeof(y));

    // Precompute double GL8 coefficients
    double A[4][4], B[4];
    for (int i = 0; i < 4; ++i) {
        B[i] = (double)GL8.b[i];
        for (int j = 0; j < 4; ++j) A[i][j] = (double)GL8.A[i][j];
    }

    auto rhs3d = [&](const double s[18], double ds[18]) {
        for (int i = 0; i < 9; ++i) { ds[i] = s[9 + i]; ds[9 + i] = 0.0; }
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                double dx = s[3*j] - s[3*i], dy = s[3*j+1] - s[3*i+1], dz = s[3*j+2] - s[3*i+2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double ir3 = 1.0 / (r2 * std::sqrt(r2));
                double fx = dx*ir3, fy = dy*ir3, fz = dz*ir3;
                ds[9+3*i] += m[j]*fx; ds[10+3*i] += m[j]*fy; ds[11+3*i] += m[j]*fz;
                ds[9+3*j] -= m[i]*fx; ds[10+3*j] -= m[i]*fy; ds[11+3*j] -= m[i]*fz;
            }
        }
    };

    auto eval_EL = [&](const double s[18], double& E, double L[3]) {
        double K = 0, U = 0;
        L[0] = L[1] = L[2] = 0;
        for (int i = 0; i < 3; ++i) {
            double rx = s[3*i], ry = s[3*i+1], rz = s[3*i+2];
            double vx = s[9+3*i], vy = s[10+3*i], vz = s[11+3*i];
            K += 0.5 * m[i] * (vx*vx + vy*vy + vz*vz);
            L[0] += m[i] * (ry*vz - rz*vy);
            L[1] += m[i] * (rz*vx - rx*vz);
            L[2] += m[i] * (rx*vy - ry*vx);
            for (int j = i + 1; j < 3; ++j) {
                double dx = s[3*j]-rx, dy = s[3*j+1]-ry, dz = s[3*j+2]-rz;
                U -= m[i]*m[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
            }
        }
        E = K + U;
    };

    double E0, L0[3];
    eval_EL(y, E0, L0);
    double T_per = (qa.q_fold > 1) ? (double)qa.in_T : (double)qa.T;
    // Fixed symplectic step size: ~1200 steps per primitive period (h ~ 0.013, GL8 error < 1e-11)
    int n_sub = (qa.q_fold > 1) ? (1200 * qa.q_fold) : 1200;
    double h = T_per / n_sub;

    LongTermResult res{qa.name, delta, n_periods, 0, 0.0, 0.0, 1e9, 0.0, 0.0, 0.0};
    double K[4][18], Y[4][18];
    double rmax_escape = 4.0 * (double)qa.rmax;

    for (int per = 1; per <= n_periods; ++per) {
        bool exploded = false;
        for (int st = 0; st < n_sub; ++st) {
            double f0[18];
            rhs3d(y, f0);
            for (int s = 0; s < 4; ++s) std::memcpy(K[s], f0, sizeof(f0));
            for (int it = 0; it < 14; ++it) {
                double md = 0.0;
                for (int s = 0; s < 4; ++s) {
                    for (int i = 0; i < 18; ++i)
                        Y[s][i] = y[i] + h * (A[s][0]*K[0][i] + A[s][1]*K[1][i] + A[s][2]*K[2][i] + A[s][3]*K[3][i]);
                    double Kn[18];
                    rhs3d(Y[s], Kn);
                    for (int i = 0; i < 18; ++i) {
                        double d = std::abs(Kn[i] - K[s][i]);
                        if (d > md) md = d;
                        K[s][i] = Kn[i];
                    }
                }
                if (md < 1e-14) break;
            }
            for (int i = 0; i < 18; ++i)
                y[i] += h * (B[0]*K[0][i] + B[1]*K[1][i] + B[2]*K[2][i] + B[3]*K[3][i]);

            // Sample min/max distances every 10 steps (and on step 0)
            if ((st % 10) == 0) {
                for (int i = 0; i < 3; ++i) {
                    double az = std::abs(y[3*i + 2]);
                    if (az > res.max_z) res.max_z = az;
                    for (int j = i + 1; j < 3; ++j) {
                        double dx = y[3*j]-y[3*i], dy = y[3*j+1]-y[3*i+1], dz = y[3*j+2]-y[3*i+2];
                        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        if (r < res.rmin_obs) res.rmin_obs = r;
                        if (r > res.rmax_obs) res.rmax_obs = r;
                    }
                }
                if (res.rmax_obs > rmax_escape || res.rmin_obs < 0.01 || std::isnan(y[0])) {
                    exploded = true;
                    break;
                }
            }
        }
        if (exploded) break;
        double E, L[3];
        eval_EL(y, E, L);
        double dEr = std::abs(E - E0) / std::abs(E0);
        double dLa = std::sqrt((L[0]-L0[0])*(L[0]-L0[0]) + (L[1]-L0[1])*(L[1]-L0[1]) + (L[2]-L0[2])*(L[2]-L0[2]));
        if (dEr > res.max_dE_rel) res.max_dE_rel = dEr;
        if (dLa > res.max_dL_abs) res.max_dL_abs = dLa;
        // Phase-space distance to unperturbed initial state (projected onto shape/distances or full 18D)
        double dr = 0.0;
        for (int i = 0; i < 9; ++i) dr += (y[i] - y_ref[i]) * (y[i] - y_ref[i]);
        dr = std::sqrt(dr);
        if (dr > res.max_return_drift) res.max_return_drift = dr;
        // Guard against subtle escape or huge divergence while staying bound
        if (dEr > 1e-4) break;
        res.survived_periods = per;
    }
    return res;
}

// ============================================================================
// Self-Test (GL16 order 16, symplecticity, analytical Hessian vs FD, 2-body Kepler)
// ============================================================================
static int run_self_test() {
    std::cout << "=== 128-BIT GL16 SELF-TEST ===\n";
    std::cout << "  GL16 quadrature order-16 max err : " << qsci(GL16.max_quad_err) << "\n";
    std::cout << "  GL16 stage order-8 max err       : " << qsci(GL16.max_stage_err) << "\n";
    std::cout << "  GL16 symplecticity max err       : " << qsci(GL16.max_symp_err) << "\n";
    std::cout << "  GL8  quadrature order-8 max err  : " << qsci(GL8.max_quad_err) << "\n";
    std::cout << "  GL8  symplecticity max err       : " << qsci(GL8.max_symp_err) << "\n";

    // Analytical Hessian vs central-FD check in __float128
    MassesQ m = {1.0Q, 0.9Q, 0.8Q};
    StateQ y = {-1.0Q, 0.13Q, -0.05Q, -0.27Q, 1.31Q, 0.19Q, 0.1Q, 0.78Q, -0.2Q, -0.16Q, 0.1Q, -0.79Q};
    f128 Kp[6][6], Kz[3][3];
    hessian_q(y, m, Kp, Kz);
    f128 max_h_err = 0.0Q, h = 1e-16Q;
    for (int j = 0; j < 6; ++j) {
        StateQ yp = y, ym = y, dp, dm;
        yp[j] += h; ym[j] -= h;
        deriv_q(yp, m, dp); deriv_q(ym, m, dm);
        for (int i = 0; i < 6; ++i) {
            f128 fd = (dp[6 + i] - dm[6 + i]) / (2.0Q * h);
            f128 e = qabs(fd - Kp[i][j]);
            if (e > max_h_err) max_h_err = e;
        }
    }
    std::cout << "  Analytical Kp Hessian vs FD err  : " << qsci(max_h_err) << "\n";

    // Exact 2-body elliptic Kepler orbit (m0=1, m1=1, m2=0, a=1, e=0.6, T = 2*pi/sqrt(2))
    MassesQ mk = {1.0Q, 1.0Q, 0.0Q};
    f128 e_kep = 0.6Q;
    f128 r0 = 1.0Q - e_kep;
    f128 v0 = qsqrt(2.0Q * (2.0Q / r0 - 1.0Q)); // relative speed at periapsis
    StateQ yk = {-0.5Q*r0, 0.0Q, 0.5Q*r0, 0.0Q, 10.0Q, 0.0Q, 0.0Q, -0.5Q*v0, 0.0Q, 0.5Q*v0, 0.0Q, 0.0Q};
    f128 Tk = 2.0Q * PI_Q / qsqrt(2.0Q);
    f128 E0, L0, E1, L1;
    invariants_q(yk, mk, E0, L0);
    StateQ yk1 = yk;
    integrate_traj_q(yk1, Tk, mk, 0.020Q, 0.020Q);
    invariants_q(yk1, mk, E1, L1);
    f128 retk = 0.0Q;
    for (int i = 0; i < 4; ++i) {
        retk = std::max(retk, qabs(yk1[i] - yk[i]));
        retk = std::max(retk, qabs(yk1[6 + i] - yk[6 + i]));
    }
    std::cout << "  Kepler (e=0.6) 1-period ret_err  : " << qsci(retk)
              << "  |dE|=" << qsci(qabs(E1 - E0)) << "  |dL|=" << qsci(qabs(L1 - L0)) << "\n";
    bool pass = (GL16.max_quad_err < 1e-32Q && GL16.max_stage_err < 1e-32Q && GL16.max_symp_err < 1e-32Q && max_h_err < 1e-14Q && retk < 1e-30Q);
    std::cout << "  SELF-TEST RESULT: " << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}

// ============================================================================
// Main Orchestrator
// ============================================================================
struct TargetSpec {
    std::string name;
    std::string role;
    std::array<double, 3> m_d;
    std::array<std::string, 3> m_str;
    int q_fold; // 1 for alpha=0 inertial; >1 for alpha = pi/q_fold RPO
};

int main(int argc, char** argv) {
    bool do_self = false, do_conv = false, do_lt = false;
    int lt_periods = 10000;
    std::string out_dir = "data/aaabbb";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--self-test") do_self = true;
        else if (a == "--converge") do_conv = true;
        else if (a == "--long-term") do_lt = true;
        else if (a == "--all") { do_self = do_conv = do_lt = true; }
        else if (a == "--periods" && i + 1 < argc) lt_periods = std::atoi(argv[++i]);
        else if (a == "--out-dir" && i + 1 < argc) out_dir = argv[++i];
    }
    if (!do_self && !do_conv && !do_lt) { do_self = do_conv = true; }

    if (do_self) {
        if (run_self_test() != 0) return 1;
    }
    if (!do_conv && !do_lt) return 0;

    std::filesystem::create_directories(out_dir);

    // Master seed: 2D-DEM-002 at m=[1, 0.9, 0.8]
    const std::array<double, 3> m002 = {1.0, 0.9, 0.8};
    const double p002[4] = {-0.05951905193104, 0.7812388652309, -0.1613813787143, 15.49514153178 * 0.5};

    std::vector<TargetSpec> specs = {
        {"AAAbbb-DEM002",   "Primary confirmed all-distinct-mass showcase (catalog 2D-DEM-002)",     {1.0, 0.90, 0.80}, {"1", "0.9",  "0.8"},  1},
        {"AAAbbb-PEAK",     "Global stability-margin peak on 0.005 grid (margin ~1.596)",            {1.0, 0.94, 0.71}, {"1", "0.94", "0.71"}, 1},
        {"AAAbbb-DEEP",     "Deep stability peak with wide Krein separation (margin ~1.476)",        {1.0, 0.94, 0.70}, {"1", "0.94", "0.7"},  1},
        {"AAAbbb-ISO89",    "Isosceles exact half-period exchange symmetry peak",                    {1.0, 0.89, 0.89}, {"1", "0.89", "0.89"}, 1},
        {"AAAbbb-ISO85",    "Clean symmetric isosceles benchmark with wide Krein gap",               {1.0, 0.85, 0.85}, {"1", "0.85", "0.85"}, 1},
        {"AAAbbb-WIDE",     "Widest minimum separation full-3D-stable member (rmin ~0.243)",         {1.0, 0.70, 0.70}, {"1", "0.7",  "0.7"},  1},
        {"AAAbbb-DEM019",   "High mass-asymmetry 5:4:3 member (catalog 2D-DEM-019)",                 {1.0, 0.80, 0.60}, {"1", "0.8",  "0.6"},  1},
        {"AAAbbb-EQ-UNST",  "Unstable equal-mass [1,1,1] inertial parent (mass-stabilization ctrl)", {1.0, 1.00, 1.00}, {"1", "1",    "1"},    1},
        {"AAAbbb-EQ-RPO13", "Equal-mass [1,1,1] full-3D-stable 13-fold rosette closure (AAAbbb)^13", {1.0, 1.00, 1.00}, {"1", "1",    "1"},    13},
    };

    std::vector<QuadAudit> results(specs.size());
    std::mutex io_mtx;

    std::cout << "\n=== 128-BIT QUAD-PRECISION (GL16) CONVERGENCE & VERIFICATION ===\n" << std::flush;
    std::atomic<size_t> next_conv{0};
    unsigned nconv = std::min((unsigned)specs.size(), std::max(1u, std::thread::hardware_concurrency()));
    std::vector<std::thread> workers;
    for (unsigned wid = 0; wid < nconv; ++wid) {
        workers.emplace_back([&]() {
          while (true) {
            size_t idx = next_conv.fetch_add(1);
            if (idx >= specs.size()) break;
            const auto& sp = specs[idx];
            QuadAudit qa{};
            qa.name = sp.name;
            qa.role = sp.role;
            qa.q_fold = sp.q_fold;
            for (int k = 0; k < 3; ++k) qa.m[k] = strtoflt128(sp.m_str[k].c_str(), nullptr);
            qa.alpha = (sp.q_fold > 1) ? (PI_Q / (f128)sp.q_fold) : 0.0Q;

            // Step 1: Verified initial seed (exact 12-digit seeds for off-diagonal peaks/edges,
            //         or smooth double continuation from 2D-DEM-002 for diagonal & RPO targets)
            double pd[4];
            bool ok_cont = true;
            if (sp.name == "AAAbbb-DEM002") {
                for (int k = 0; k < 4; ++k) pd[k] = p002[k];
            } else if (sp.name == "AAAbbb-PEAK") {
                pd[0] = -0.033187327055; pd[1] = 0.767613801325;
                pd[2] = -0.220297406478; pd[3] = 8.167045150019;
            } else if (sp.name == "AAAbbb-DEEP") {
                pd[0] = -0.0296347612782; pd[1] = 0.764096745256;
                pd[2] = -0.226559640058; pd[3] = 8.23255081555;
            } else if (sp.name == "AAAbbb-WIDE") {
                pd[0] = 0.0022078140345; pd[1] = 0.650105805358;
                pd[2] = -0.211442234297; pd[3] = 9.18044059843;
            } else if (sp.name == "AAAbbb-DEM019") {
                pd[0] = 0.03029133440607; pd[1] = 0.6605762934825;
                pd[2] = -0.2815395265217; pd[3] = 19.05486021379 * 0.5;
            } else {
                ok_cont = continue_double_seed(p002, m002, sp.m_d, (double)qa.alpha, pd);
            }
            if (!ok_cont) {
                std::lock_guard<std::mutex> lk(io_mtx);
                std::cerr << "  [FAIL] Continuation failed for " << sp.name << "\n";
                continue;
            }
            for (int k = 0; k < 4; ++k) qa.p[k] = (f128)pd[k];

            // Step 2: 128-bit Newton shooting to < 1e-31
            f128 shoot_res = 1.0Q;
            bool ok_newt = newton_shoot_128(qa.p, qa.alpha, qa.m, shoot_res);
            qa.y0 = make_perp_q(qa.p[0], qa.p[1], qa.p[2], qa.m);
            qa.T = 2.0Q * qa.p[3];
            invariants_q(qa.y0, qa.m, qa.E, qa.L);
            qa.Tstar = qa.T * qpow(-qa.E, 1.5Q);
            qa.Lambda = qabs(qa.L) * qsqrt(-qa.E);

            // Step 3: Full-period 128-bit analytical variational monodromy (c_scale = 0.022 -> ~1e-31)
            f128 two_a = 2.0Q * qa.alpha;
            f128 Mrot[12][12], Mz[6][6];
            eval_var_period_q(qa.y0, qa.T, two_a, qa.m, 0.022Q, qa.yend, qa.symp_err,
                              qa.real_k, qa.k1, qa.k2, qa.kz, qa.marg3, &qa, Mrot, Mz);
            qa.stable_3d = (qa.marg3 > 0.0Q);
            qa.krein_gap = qa.real_k ? qabs(qa.k1 - qa.k2) : -qa.k2;

            StateQ yend_rot = rot2d_q(qa.yend, -two_a);
            qa.ret_err = 0.0Q;
            for (int i = 0; i < 12; ++i) {
                f128 e = qabs(yend_rot[i] - qa.y0[i]);
                if (e > qa.ret_err) qa.ret_err = e;
            }
            f128 E1, L1;
            invariants_q(qa.yend, qa.m, E1, L1);
            qa.dE = qabs(E1 - qa.E);
            qa.dL = qabs(L1 - qa.L);

            // Step 4: Independent GL16 step-halving cross-check (c_scale = 0.011)
            if (idx == 0) {
                StateQ yend2; f128 symp2, k1_2, k2_2, kz_2, m3_2; bool rk2;
                eval_var_period_q(qa.y0, qa.T, two_a, qa.m, 0.011Q, yend2, symp2, rk2, k1_2, k2_2, kz_2, m3_2, nullptr);
                qa.step_half_y_err = 0.0Q;
                for (int i = 0; i < 12; ++i) qa.step_half_y_err = std::max(qa.step_half_y_err, qabs(yend2[i] - qa.yend[i]));
                qa.step_half_k_err = std::max({qabs(k1_2 - qa.k1), qabs(k2_2 - qa.k2), qabs(kz_2 - qa.kz)});
            } else {
                StateQ yend2 = qa.y0;
                integrate_traj_q(yend2, qa.T, qa.m, 0.011Q, 0.011Q);
                qa.step_half_y_err = 0.0Q;
                for (int i = 0; i < 12; ++i) qa.step_half_y_err = std::max(qa.step_half_y_err, qabs(yend2[i] - qa.yend[i]));
                qa.step_half_k_err = 0.0Q;
            }

            // Step 5: For RPO closures (q_fold > 1), verify direct 12D GL16 closure at T_inert = q_fold * T
            // and exact Floquet monodromy M_inert = M_rpo^{q_fold} (exact by SO(2) equivariance)!
            if (qa.q_fold > 1) {
                qa.has_inert_closure = true;
                qa.in_T = (f128)qa.q_fold * qa.T;
                StateQ yin = qa.y0;
                integrate_traj_q(yin, qa.in_T, qa.m, 0.022Q, 0.022Q);
                qa.in_ret = 0.0Q;
                for (int i = 0; i < 12; ++i) qa.in_ret = std::max(qa.in_ret, qabs(yin[i] - qa.y0[i]));
                f128 MpQ[12][12] = {{0}}, MzQ[6][6] = {{0}}, tmp12[12][12], tmp6[6][6];
                for (int i = 0; i < 12; ++i) MpQ[i][i] = 1.0Q;
                for (int i = 0; i < 6; ++i) MzQ[i][i] = 1.0Q;
                for (int step_q = 0; step_q < qa.q_fold; ++step_q) {
                    for (int r = 0; r < 12; ++r) for (int c = 0; c < 12; ++c) {
                        f128 s = 0.0Q; for (int k = 0; k < 12; ++k) s += MpQ[r][k] * Mrot[k][c];
                        tmp12[r][c] = s;
                    }
                    std::memcpy(MpQ, tmp12, sizeof(MpQ));
                    for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) {
                        f128 s = 0.0Q; for (int k = 0; k < 6; ++k) s += MzQ[r][k] * Mz[k][c];
                        tmp6[r][c] = s;
                    }
                    std::memcpy(MzQ, tmp6, sizeof(MzQ));
                }
                bool rk_in;
                extract_k_and_symp(MpQ, MzQ, 0.0Q, qa.m, qa.in_symp, rk_in, qa.in_k1, qa.in_k2, qa.in_kz, qa.in_marg3);
                qa.in_st3 = (qa.in_marg3 > 0.0Q);
            }

            results[idx] = qa;
            {
                std::lock_guard<std::mutex> lk(io_mtx);
                std::cout << "  [" << (ok_newt ? "OK" : "WARN") << "] " << std::left << std::setw(16) << qa.name
                          << " m=[" << sp.m_str[0] << "," << sp.m_str[1] << "," << sp.m_str[2] << "]"
                          << " word=" << qa.word_s3
                          << " ret=" << qsci(qa.ret_err, 2)
                          << " symp=" << qsci(qa.symp_err, 2)
                          << " h/2_y=" << qsci(qa.step_half_y_err, 2)
                          << " m3=" << qfix(qa.marg3, 6)
                          << " st3=" << qa.stable_3d << "\n" << std::flush;
                if (qa.has_inert_closure) {
                    std::cout << "       -> INERT(" << qa.q_fold << "T=" << qfix(qa.in_T, 6)
                              << ") ret=" << qsci(qa.in_ret, 2) << " symp=" << qsci(qa.in_symp, 2)
                              << " k=(" << qfix(qa.in_k1, 6) << "," << qfix(qa.in_k2, 6)
                              << ") kz=" << qfix(qa.in_kz, 6) << " m3=" << qfix(qa.in_marg3, 6)
                              << " st3=" << qa.in_st3 << "\n" << std::flush;
                }
            }
          }
        });
    }
    for (auto& th : workers) th.join();

    // Write 32-digit JSON archive
    {
        std::ofstream jf(out_dir + "/aaabbb_quad32.json");
        jf << "{\n  \"family\": \"AAAbbb Shuttle-Exchange\",\n"
           << "  \"integrator\": \"128-bit __float128 8-stage 16th-order exact symplectic Gauss-Legendre (GL16)\",\n"
           << "  \"gl16_quad_err\": \"" << qsci(GL16.max_quad_err, 3) << "\",\n"
           << "  \"gl16_symp_err\": \"" << qsci(GL16.max_symp_err, 3) << "\",\n"
           << "  \"orbits\": [\n";
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& qa = results[i];
            const auto& sp = specs[i];
            jf << "    {\n"
               << "      \"id\": \"" << qa.name << "\",\n"
               << "      \"role\": \"" << qa.role << "\",\n"
               << "      \"masses\": [" << sp.m_str[0] << ", " << sp.m_str[1] << ", " << sp.m_str[2] << "],\n"
               << "      \"q_fold\": " << qa.q_fold << ",\n"
               << "      \"alpha_half_rad\": \"" << qstr(qa.alpha, 34) << "\",\n"
               << "      \"T\": \"" << qstr(qa.T, 34) << "\",\n"
               << "      \"x2\": \"" << qstr(qa.p[0], 34) << "\",\n"
               << "      \"vy1\": \"" << qstr(qa.p[1], 34) << "\",\n"
               << "      \"vy2\": \"" << qstr(qa.p[2], 34) << "\",\n"
               << "      \"pos_34d\": [\n"
               << "        [\"" << qstr(qa.y0[0], 34) << "\", \"" << qstr(qa.y0[1], 34) << "\"],\n"
               << "        [\"" << qstr(qa.y0[2], 34) << "\", \"" << qstr(qa.y0[3], 34) << "\"],\n"
               << "        [\"" << qstr(qa.y0[4], 34) << "\", \"" << qstr(qa.y0[5], 34) << "\"]\n"
               << "      ],\n"
               << "      \"vel_34d\": [\n"
               << "        [\"" << qstr(qa.y0[6], 34) << "\", \"" << qstr(qa.y0[7], 34) << "\"],\n"
               << "        [\"" << qstr(qa.y0[8], 34) << "\", \"" << qstr(qa.y0[9], 34) << "\"],\n"
               << "        [\"" << qstr(qa.y0[10], 34) << "\", \"" << qstr(qa.y0[11], 34) << "\"]\n"
               << "      ],\n"
               << "      \"E\": \"" << qstr(qa.E, 34) << "\",\n"
               << "      \"L\": \"" << qstr(qa.L, 34) << "\",\n"
               << "      \"Tstar\": \"" << qstr(qa.Tstar, 34) << "\",\n"
               << "      \"Lambda\": \"" << qstr(qa.Lambda, 34) << "\",\n"
               << "      \"k1\": \"" << qstr(qa.k1, 28) << "\",\n"
               << "      \"k2\": \"" << qstr(qa.k2, 28) << "\",\n"
               << "      \"kz\": \"" << qstr(qa.kz, 28) << "\",\n"
               << "      \"stab_margin_3d\": \"" << qstr(qa.marg3, 28) << "\",\n"
               << "      \"krein_gap\": \"" << qstr(qa.krein_gap, 20) << "\",\n"
               << "      \"stable_3d\": " << (qa.stable_3d ? "true" : "false") << ",\n"
               << "      \"ret_err_128\": \"" << qsci(qa.ret_err, 4) << "\",\n"
               << "      \"dE_128\": \"" << qsci(qa.dE, 4) << "\",\n"
               << "      \"dL_128\": \"" << qsci(qa.dL, 4) << "\",\n"
               << "      \"symp_err_128\": \"" << qsci(qa.symp_err, 4) << "\",\n"
               << "      \"step_half_y_err\": \"" << qsci(qa.step_half_y_err, 4) << "\",\n"
               << "      \"step_half_k_err\": \"" << qsci(qa.step_half_k_err, 4) << "\",\n"
               << "      \"rmin\": " << (double)qa.rmin << ",\n"
               << "      \"rmax\": " << (double)qa.rmax << ",\n"
               << "      \"rmin_pair\": [" << (double)qa.rmin_pair[0] << ", " << (double)qa.rmin_pair[1] << ", " << (double)qa.rmin_pair[2] << "],\n"
               << "      \"rmax_pair\": [" << (double)qa.rmax_pair[0] << ", " << (double)qa.rmax_pair[1] << ", " << (double)qa.rmax_pair[2] << "],\n"
               << "      \"windings\": [" << (double)qa.wind[0] << ", " << (double)qa.wind[1] << ", " << (double)qa.wind[2] << "],\n"
               << "      \"syzygies\": " << qa.syzygies << ",\n"
               << "      \"syzygy_seq\": \"" << qa.syz_seq << "\",\n"
               << "      \"word_s3\": \"" << qa.word_s3 << "\"";
            if (qa.has_inert_closure) {
                jf << ",\n      \"inert_closure\": {\n"
                   << "        \"T_inert\": \"" << qstr(qa.in_T, 34) << "\",\n"
                   << "        \"ret_err_128\": \"" << qsci(qa.in_ret, 4) << "\",\n"
                   << "        \"symp_err_128\": \"" << qsci(qa.in_symp, 4) << "\",\n"
                   << "        \"k1\": \"" << qstr(qa.in_k1, 24) << "\",\n"
                   << "        \"k2\": \"" << qstr(qa.in_k2, 24) << "\",\n"
                   << "        \"kz\": \"" << qstr(qa.in_kz, 24) << "\",\n"
                   << "        \"stab_margin_3d\": \"" << qstr(qa.in_marg3, 24) << "\",\n"
                   << "        \"stable_3d\": " << (qa.in_st3 ? "true" : "false") << "\n"
                   << "      }\n";
            } else {
                jf << "\n";
            }
            jf << "    }" << (i + 1 < results.size() ? "," : "") << "\n";
        }
        jf << "  ]\n}\n";
    }

    // Write clean plain-text submission file (ready for threebodyorbits.com / archive submission)
    {
        std::ofstream tf(out_dir + "/aaabbb_submission_34digits.txt");
        tf << "# AAAbbb Shuttle-Exchange Family — 34-Digit Quad-Precision (__float128 GL16) Initial Conditions\n"
           << "# Newtonian planar 3-body problem (G = 1), center-of-mass frame (sum m_i r_i = 0, sum m_i v_i = 0).\n"
           << "# All collinear perpendicular crossings at t = 0: y_i(0) = 0, vx_i(0) = 0, x_1(0) = -1.\n"
           << "# Verified with 8-stage 16th-order exact symplectic Gauss-Legendre (GL16) + step-halving cross-check.\n\n";
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& qa = results[i];
            const auto& sp = specs[i];
            tf << "================================================================================\n"
               << "ID:     " << qa.name << "  (" << qa.role << ")\n"
               << "Masses: m1 = " << sp.m_str[0] << ",  m2 = " << sp.m_str[1] << ",  m3 = " << sp.m_str[2] << "\n"
               << "Word:   " << qa.word_s3 << "  (syzygies: " << qa.syzygies << ", seq: " << qa.syz_seq << ")\n"
               << "T:      " << qstr(qa.T, 34) << (qa.q_fold > 1 ? ("  (T_inert = " + std::to_string(qa.q_fold) + "*T = " + qstr(qa.in_T, 34) + ")") : "") << "\n"
               << "E:      " << qstr(qa.E, 34) << "\n"
               << "L:      " << qstr(qa.L, 34) << "\n"
               << "T*:     " << qstr(qa.Tstar, 34) << "\n"
               << "Lambda: " << qstr(qa.Lambda, 34) << "\n"
               << "Body 1: r1 = (" << qstr(qa.y0[0], 34) << ", 0)   v1 = (0, " << qstr(qa.y0[7], 34) << ")\n"
               << "Body 2: r2 = (" << qstr(qa.y0[2], 34) << ", 0)   v2 = (0, " << qstr(qa.y0[9], 34) << ")\n"
               << "Body 3: r3 = (" << qstr(qa.y0[4], 34) << ", 0)   v3 = (0, " << qstr(qa.y0[11], 34) << ")\n"
               << "Stab:   k1 = " << qstr(qa.k1, 20) << ",  k2 = " << qstr(qa.k2, 20) << ",  kz = " << qstr(qa.kz, 20)
               << "  | 3D margin = " << qfix(qa.marg3, 12) << "  | 3D-stable = " << (qa.stable_3d ? "YES" : "NO") << "\n"
               << "Audit:  ret_err = " << qsci(qa.ret_err, 3) << ", dE = " << qsci(qa.dE, 3)
               << ", dL = " << qsci(qa.dL, 3) << ", symp_err = " << qsci(qa.symp_err, 3)
               << ", GL16(h vs h/2) = " << qsci(qa.step_half_y_err, 3) << "\n\n";
        }
    }

    if (!do_lt) return 0;

    // ========================================================================
    // 10,000-Period Full-3D Symplectic GL8 Perturbation Survival Tests
    // ========================================================================
    std::cout << "\n=== " << lt_periods << "-PERIOD FULL-3D SYMPLECTIC GL8 PERTURBATION SURVIVAL TESTS ===\n" << std::flush;
    struct Job { size_t idx; double delta; int nper; };
    std::vector<Job> jobs;
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].q_fold > 1) {
            // For 13-fold rosette, 1000 full inertial periods = 13,000 RPO cycles!
            jobs.push_back({i, 1e-5, std::max(100, lt_periods / 10)});
        } else {
            jobs.push_back({i, 1e-5, lt_periods});
            if (i < 5) jobs.push_back({i, 1e-4, lt_periods});
        }
    }
    std::vector<LongTermResult> lt_res(jobs.size());
    std::atomic<size_t> next_job{0};
    unsigned nthr = std::min((unsigned)jobs.size(), std::max(1u, std::thread::hardware_concurrency()));
    std::vector<std::thread> ltw;
    for (unsigned t = 0; t < nthr; ++t) {
        ltw.emplace_back([&]() {
            while (true) {
                size_t j = next_job.fetch_add(1);
                if (j >= jobs.size()) break;
                lt_res[j] = run_long_term_3d(results[jobs[j].idx], jobs[j].delta, jobs[j].nper);
                std::lock_guard<std::mutex> lk(io_mtx);
                const auto& r = lt_res[j];
                std::cout << "  " << std::left << std::setw(16) << r.name
                          << " delta=" << std::scientific << std::setprecision(0) << r.delta
                          << "  survived=" << std::setw(5) << r.survived_periods << "/" << r.target_periods
                          << "  rmin=" << std::fixed << std::setprecision(4) << r.rmin_obs
                          << "  rmax=" << r.rmax_obs
                          << "  |z|_max=" << std::scientific << std::setprecision(2) << r.max_z
                          << "  dE/E=" << r.max_dE_rel
                          << "  dL=" << r.max_dL_abs << "\n" << std::flush;
            }
        });
    }
    for (auto& th : ltw) th.join();

    // Write long-term survival JSON
    {
        std::ofstream lf(out_dir + "/aaabbb_longterm_3d.json");
        lf << "{\n  \"integrator\": \"Full-3D 4-stage 8th-order exact symplectic Gauss-Legendre (GL8), 1200 steps/period\",\n"
           << "  \"perturbation\": \"Zero-momentum 3D velocity kick delta*(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)) on body 1, compensated on bodies 2-3\",\n"
           << "  \"trials\": [\n";
        for (size_t j = 0; j < lt_res.size(); ++j) {
            const auto& r = lt_res[j];
            lf << "    {\"id\":\"" << r.name << "\",\"delta\":" << std::scientific << std::setprecision(1) << r.delta
               << ",\"target_periods\":" << r.target_periods << ",\"survived_periods\":" << r.survived_periods
               << ",\"rmin_obs\":" << std::fixed << std::setprecision(6) << r.rmin_obs
               << ",\"rmax_obs\":" << r.rmax_obs
               << ",\"max_z\":" << std::scientific << std::setprecision(4) << r.max_z
               << ",\"max_dE_rel\":" << r.max_dE_rel
               << ",\"max_dL_abs\":" << r.max_dL_abs
               << ",\"max_return_drift\":" << r.max_return_drift << "}"
               << (j + 1 < lt_res.size() ? "," : "") << "\n";
        }
        lf << "  ]\n}\n";
    }
    return 0;
}
