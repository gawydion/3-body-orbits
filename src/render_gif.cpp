#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

static inline void deriv_newt(const double* s, const double* m, int nb, int dim, double* ds) {
    int half = nb * dim;
    for (int i = 0; i < half; ++i) ds[i] = s[half + i];
    for (int i = half; i < 2 * half; ++i) ds[i] = 0.0;
    for (int i = 0; i < nb; ++i) {
        for (int j = i + 1; j < nb; ++j) {
            double dx = s[dim*j] - s[dim*i];
            double dy = s[dim*j+1] - s[dim*i+1];
            double dz = (dim == 3) ? (s[dim*j+2] - s[dim*i+2]) : 0.0;
            double r2 = dx*dx + dy*dy + dz*dz;
            double ir = 1.0 / std::sqrt(r2), ir3 = ir*ir*ir;
            ds[half + dim*i]   += m[j]*dx*ir3;
            ds[half + dim*i+1] += m[j]*dy*ir3;
            ds[half + dim*j]   -= m[i]*dx*ir3;
            ds[half + dim*j+1] -= m[i]*dy*ir3;
            if (dim == 3) {
                ds[half + dim*i+2] += m[j]*dz*ir3;
                ds[half + dim*j+2] -= m[i]*dz*ir3;
            }
        }
    }
}

// Exact 1PN ADM Hamiltonian equations of motion for 3 bodies in 2D or 3D
static inline void deriv_1pn(const double* s, const double* m, int dim, double beta, double* ds) {
    // Expand to 3D canonical state (18 vars)
    double r[3][3] = {{0}}, p[3][3] = {{0}};
    int half = 3 * dim;
    for (int a = 0; a < 3; ++a) {
        r[a][0] = s[dim*a]; r[a][1] = s[dim*a+1]; r[a][2] = (dim == 3) ? s[dim*a+2] : 0.0;
        p[a][0] = s[half + dim*a]; p[a][1] = s[half + dim*a+1]; p[a][2] = (dim == 3) ? s[half + dim*a+2] : 0.0;
    }
    double dr[3][3] = {{0}}, dp[3][3] = {{0}};
    double phi[3] = {0, 0, 0}, p2[3] = {0, 0, 0};
    for (int a = 0; a < 3; ++a) {
        p2[a] = p[a][0]*p[a][0] + p[a][1]*p[a][1] + p[a][2]*p[a][2];
        double f0 = 1.0 / m[a] - beta * 0.5 * p2[a] / (m[a]*m[a]*m[a]);
        for (int k = 0; k < 3; ++k) dr[a][k] = f0 * p[a][k];
    }
    double dvec[3][3][3], rdist[3][3], invr[3][3], invr2[3][3], invr3[3][3];
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            double dx = r[a][0]-r[b][0], dy = r[a][1]-r[b][1], dz = r[a][2]-r[b][2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            double ir = 1.0 / dist, ir2 = ir*ir, ir3 = ir2*ir;
            dvec[a][b][0]=dx; dvec[a][b][1]=dy; dvec[a][b][2]=dz;
            rdist[a][b]=dist; invr[a][b]=ir; invr2[a][b]=ir2; invr3[a][b]=ir3;
            phi[a] += m[b]*ir; phi[b] += m[a]*ir;
            for (int k = 0; k < 3; ++k) {
                double f = -m[a]*m[b]*ir3*dvec[a][b][k];
                dp[a][k] += f; dp[b][k] -= f;
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = a + 1; b < 3; ++b) {
            double ir = invr[a][b], ir2 = invr2[a][b], ir3 = invr3[a][b];
            double nx = dvec[a][b][0]*ir, ny = dvec[a][b][1]*ir, nz = dvec[a][b][2]*ir;
            double pa_pb = p[a][0]*p[b][0] + p[a][1]*p[b][1] + p[a][2]*p[b][2];
            double npa = nx*p[a][0] + ny*p[a][1] + nz*p[a][2];
            double npb = nx*p[b][0] + ny*p[b][1] + nz*p[b][2];
            double K_ab = 3.0 * ((m[b]/m[a])*p2[a] + (m[a]/m[b])*p2[b]) - 7.0*pa_pb - npa*npb;
            double nv[3] = {nx, ny, nz};
            for (int k = 0; k < 3; ++k) {
                dr[a][k] += beta * (-0.5*ir) * (6.0*(m[b]/m[a])*p[a][k] - 7.0*p[b][k] - npb*nv[k]);
                dr[b][k] += beta * (-0.5*ir) * (6.0*(m[a]/m[b])*p[b][k] - 7.0*p[a][k] - npa*nv[k]);
                double dK_dra = -ir * ( (p[a][k] - npa*nv[k])*npb + (p[b][k] - npb*nv[k])*npa );
                double dH1_kin = 0.5*K_ab*ir2*nv[k] - 0.5*ir*dK_dra;
                double dH1_pot = -m[a]*m[b]*(phi[a] + phi[b])*ir3*dvec[a][b][k];
                double dH1_tot = beta * (dH1_kin + dH1_pot);
                dp[a][k] -= dH1_tot;
                dp[b][k] += dH1_tot;
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int k = 0; k < dim; ++k) {
            ds[dim*a + k] = dr[a][k];
            ds[half + dim*a + k] = dp[a][k];
        }
    }
}

static void advance(double* y, double dt, const double* m, int nb, int dim, double beta, double tol = 1e-11) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    int n = 2 * nb * dim;
    double k1[24], k2[24], k3[24], k4[24], k5[24], k6[24], k7[24], yt[24], yn[24];
    auto eval = [&](const double* s, double* ds) {
        if (beta > 0.0 && nb == 3) deriv_1pn(s, m, dim, beta, ds);
        else deriv_newt(s, m, nb, dim, ds);
    };
    double t = 0.0, h = std::min(0.005, dt);
    int steps = 0;
    while (t < dt - 1e-14 && steps < 40000) {
        if (t + h > dt) h = dt - t;
        eval(y, k1);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*a21*k1[i]; eval(yt, k2);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); eval(yt, k3);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); eval(yt, k4);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); eval(yt, k5);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); eval(yt, k6);
        for (int i=0;i<n;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); eval(yn, k7);
        double err = 0.0;
        for (int i=0;i<n;++i) {
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = tol*(1.0 + std::max(std::abs(y[i]), std::abs(yn[i])));
            double r = ei/sc; err += r*r;
        }
        err = std::sqrt(err / n);
        if (err <= 1.0) {
            for (int i=0;i<n;++i) y[i] = yn[i];
            t += h; ++steps;
            h = std::min(h * ((err==0)?5.0:std::min(5.0, std::max(0.2, 0.9*std::pow(err,-0.2)))), 0.05);
        } else {
            h *= std::max(0.1, 0.9*std::pow(err,-0.25));
            if (h < 1e-12) break;
        }
    }
}

struct Canvas {
    int W, H;
    std::vector<float> buf;

    Canvas(int w, int h) : W(w), H(h), buf(w * h * 3, 0.0f) {}

    void clear(float r, float g, float b) {
        for (int i = 0; i < W * H; ++i) {
            buf[3*i] = r; buf[3*i+1] = g; buf[3*i+2] = b;
        }
    }

    inline void blend(int x, int y, float r, float g, float b, float a) {
        if (x < 0 || x >= W || y < 0 || y >= H || a <= 0.0f) return;
        a = std::min(1.0f, a);
        int idx = 3 * (y * W + x);
        buf[idx]   = buf[idx]   * (1.0f - a) + r * a;
        buf[idx+1] = buf[idx+1] * (1.0f - a) + g * a;
        buf[idx+2] = buf[idx+2] * (1.0f - a) + b * a;
    }

    void draw_line(float x0, float y0, float x1, float y1, float r, float g, float b, float a, float thick = 1.2f) {
        float dx = x1 - x0, dy = y1 - y0;
        float len = std::hypot(dx, dy);
        if (len < 1e-4f) return;
        int steps = std::max(1, (int)std::ceil(len * 1.5f));
        float rad = thick * 0.5f;
        int ir = (int)std::ceil(rad + 1.0f);
        for (int s = 0; s <= steps; ++s) {
            float t = (float)s / steps;
            float px = x0 + dx * t, py = y0 + dy * t;
            int ix0 = (int)std::floor(px) - ir, ix1 = (int)std::ceil(px) + ir;
            int iy0 = (int)std::floor(py) - ir, iy1 = (int)std::ceil(py) + ir;
            for (int iy = iy0; iy <= iy1; ++iy) {
                for (int ix = ix0; ix <= ix1; ++ix) {
                    float d = std::hypot(ix - px, iy - py);
                    if (d < rad + 1.0f) {
                        float cov = std::clamp(rad + 0.8f - d, 0.0f, 1.0f);
                        blend(ix, iy, r, g, b, a * cov * 0.65f);
                    }
                }
            }
        }
    }

    void draw_body(float cx, float cy, float rad, float r, float g, float b) {
        int glow_r = (int)std::ceil(rad * 2.8f);
        for (int dy = -glow_r; dy <= glow_r; ++dy) {
            for (int dx = -glow_r; dx <= glow_r; ++dx) {
                float d = std::hypot((float)dx, (float)dy);
                if (d < rad * 2.8f) {
                    float t = d / (rad * 2.8f);
                    float ga = (1.0f - t) * (1.0f - t) * 0.35f;
                    blend((int)std::round(cx) + dx, (int)std::round(cy) + dy, r, g, b, ga);
                }
            }
        }
        int ir = (int)std::ceil(rad + 2.0f);
        for (int dy = -ir; dy <= ir; ++dy) {
            for (int dx = -ir; dx <= ir; ++dx) {
                float px = std::round(cx) + dx, py = std::round(cy) + dy;
                float d = std::hypot(px - cx, py - cy);
                if (d <= rad + 1.2f) {
                    float cov = std::clamp(rad + 0.7f - d, 0.0f, 1.0f);
                    float hl = std::max(0.0f, 1.0f - std::hypot(dx + rad*0.28f, dy + rad*0.28f) / (rad * 0.85f));
                    float rr = std::min(1.0f, r + hl * 0.45f);
                    float gg = std::min(1.0f, g + hl * 0.45f);
                    float bb = std::min(1.0f, b + hl * 0.45f);
                    blend((int)px, (int)py, rr, gg, bb, cov);
                    if (d > rad - 1.1f && d < rad + 0.8f) {
                        float rim = 1.0f - std::abs(d - (rad - 0.15f)) / 0.95f;
                        blend((int)px, (int)py, 1.0f, 1.0f, 1.0f, std::max(0.0f, rim) * 0.65f);
                    }
                }
            }
        }
    }

    void write_rgb24(FILE* fp) const {
        std::vector<uint8_t> out(W * H * 3);
        for (int i = 0; i < W * H * 3; ++i) {
            out[i] = (uint8_t)std::clamp((int)std::round(buf[i] * 255.0f), 0, 255);
        }
        fwrite(out.data(), 1, out.size(), fp);
    }
};

int main(int argc, char** argv) {
    // Format:
    // nb dim beta W H frames T_ghost T_anim yaw0 yaw_span pitch m[0..nb-1] y0[0..2*nb*dim-1]
    int nb, dim, W, H, frames;
    double beta, T_ghost, T_anim, yaw0, yaw_span, pitch;
    if (!(std::cin >> nb >> dim >> beta >> W >> H >> frames >> T_ghost >> T_anim >> yaw0 >> yaw_span >> pitch))
        return 1;
    if (nb < 3 || nb > 4 || (dim != 2 && dim != 3)) return 2;
    double m[4] = {0}, y0[24] = {0};
    for (int i = 0; i < nb; ++i) std::cin >> m[i];
    int n = 2 * nb * dim;
    for (int i = 0; i < n; ++i) std::cin >> y0[i];

    int N_ghost = (T_ghost > 120.0) ? 3600 : ((T_ghost > 50.0) ? 2400 : 800);
    std::vector<std::array<double, 12>> ghost(N_ghost + 1);
    double yt[24];
    for (int i = 0; i < n; ++i) yt[i] = y0[i];
    double maxR = 0.5;
    for (int k = 0; k <= N_ghost; ++k) {
        std::array<double, 12> pt = {0};
        for (int b = 0; b < nb; ++b) {
            pt[3*b]   = yt[dim*b];
            pt[3*b+1] = yt[dim*b+1];
            pt[3*b+2] = (dim == 3) ? yt[dim*b+2] : 0.0;
            double r = std::sqrt(pt[3*b]*pt[3*b] + pt[3*b+1]*pt[3*b+1] + pt[3*b+2]*pt[3*b+2]);
            if (r > maxR) maxR = r;
        }
        ghost[k] = pt;
        if (k < N_ghost) advance(yt, T_ghost / N_ghost, m, nb, dim, beta, 1e-11);
    }

    int sub = 4;
    int N_anim = frames * sub;
    std::vector<std::array<double, 12>> anim(N_anim + 1);
    for (int i = 0; i < n; ++i) yt[i] = y0[i];
    for (int k = 0; k <= N_anim; ++k) {
        std::array<double, 12> pt = {0};
        for (int b = 0; b < nb; ++b) {
            pt[3*b]   = yt[dim*b];
            pt[3*b+1] = yt[dim*b+1];
            pt[3*b+2] = (dim == 3) ? yt[dim*b+2] : 0.0;
        }
        anim[k] = pt;
        if (k < N_anim) advance(yt, T_anim / N_anim, m, nb, dim, beta, 1e-11);
    }

    double scale = std::min(W, H) * 0.42 / maxR;
    const float cols[4][3] = {
        {0.231f, 0.510f, 0.965f}, // blue    #3b82f6
        {0.961f, 0.620f, 0.043f}, // amber   #f59e0b
        {0.925f, 0.282f, 0.600f}, // pink    #ec4899
        {0.063f, 0.725f, 0.506f}  // emerald #10b981
    };

    Canvas cv(W, H);
    int trail_len = std::min(N_anim - 1, (int)(N_anim * 0.42));

    for (int f = 0; f < frames; ++f) {
        cv.clear(0.043f, 0.055f, 0.078f);

        double yaw = yaw0 + yaw_span * ((double)f / frames);
        double cy = std::cos(yaw), sy = std::sin(yaw);
        double cp = std::cos(pitch), sp = std::sin(pitch);

        auto proj = [&](double x, double y, double z, float& sx, float& sy_out, float& sz) {
            if (dim == 2) {
                sx = (float)(W * 0.5 + x * scale);
                sy_out = (float)(H * 0.5 - y * scale);
                sz = 0.0f;
            } else {
                double rx = x * cy - y * sy;
                double ry = x * sy + y * cy;
                double px = rx;
                double py = ry * sp + z * cp;
                double pz = ry * cp - z * sp;
                sx = (float)(W * 0.5 + px * scale);
                sy_out = (float)(H * 0.5 - py * scale);
                sz = (float)pz;
            }
        };

        float ox, oy, oz;
        proj(0, 0, 0, ox, oy, oz);
        if (dim == 2) {
            cv.draw_line(ox - 10, oy, ox + 10, oy, 0.4f, 0.48f, 0.6f, 0.28f, 1.0f);
            cv.draw_line(ox, oy - 10, ox, oy + 10, 0.4f, 0.48f, 0.6f, 0.28f, 1.0f);
        } else {
            double ax = maxR * 0.35;
            float axx, axy, axz;
            proj(ax, 0, 0, axx, axy, axz); cv.draw_line(ox, oy, axx, axy, 0.5f, 0.55f, 0.65f, 0.25f, 1.0f);
            proj(0, ax, 0, axx, axy, axz); cv.draw_line(ox, oy, axx, axy, 0.5f, 0.55f, 0.65f, 0.25f, 1.0f);
            proj(0, 0, ax, axx, axy, axz); cv.draw_line(ox, oy, axx, axy, 0.3f, 0.8f, 0.9f, 0.35f, 1.2f);
        }

        float ghost_alpha = (T_ghost > 120.0) ? 0.11f : ((T_ghost > 50.0) ? 0.14f : 0.22f);
        for (int b = 0; b < nb; ++b) {
            float px0, py0, pz0;
            proj(ghost[0][3*b], ghost[0][3*b+1], ghost[0][3*b+2], px0, py0, pz0);
            for (int k = 1; k <= N_ghost; ++k) {
                float px1, py1, pz1;
                proj(ghost[k][3*b], ghost[k][3*b+1], ghost[k][3*b+2], px1, py1, pz1);
                cv.draw_line(px0, py0, px1, py1, cols[b][0], cols[b][1], cols[b][2], ghost_alpha, 1.2f);
                px0 = px1; py0 = py1;
            }
        }

        int cur_idx = f * sub;
        double jump_thresh = std::max(0.25, 0.08 * maxR * maxR);
        for (int b = 0; b < nb; ++b) {
            for (int step = 0; step < trail_len; ++step) {
                int i0 = (cur_idx - trail_len + step + N_anim) % N_anim;
                int i1 = (i0 + 1) % N_anim;
                double dx = anim[i1][3*b] - anim[i0][3*b];
                double dy = anim[i1][3*b+1] - anim[i0][3*b+1];
                double dz = anim[i1][3*b+2] - anim[i0][3*b+2];
                if (dx*dx + dy*dy + dz*dz > jump_thresh) continue;
                float frac = (float)(step + 1) / trail_len;
                float alpha = frac * frac * 0.92f;
                float x0, y0p, z0, x1, y1p, z1;
                proj(anim[i0][3*b], anim[i0][3*b+1], anim[i0][3*b+2], x0, y0p, z0);
                proj(anim[i1][3*b], anim[i1][3*b+1], anim[i1][3*b+2], x1, y1p, z1);
                cv.draw_line(x0, y0p, x1, y1p, cols[b][0], cols[b][1], cols[b][2], alpha, 2.2f);
            }
        }

        struct BDraw { int b; float x, y, z, r; };
        std::vector<BDraw> bds(nb);
        for (int b = 0; b < nb; ++b) {
            float sx, sy_out, sz;
            proj(anim[cur_idx][3*b], anim[cur_idx][3*b+1], anim[cur_idx][3*b+2], sx, sy_out, sz);
            float rad = 4.6f + 3.4f * (float)std::cbrt(std::max(0.02, m[b]));
            bds[b] = {b, sx, sy_out, sz, rad};
        }
        std::sort(bds.begin(), bds.end(), [](const BDraw& a, const BDraw& b) { return a.z < b.z; });
        for (const auto& bd : bds) {
            cv.draw_body(bd.x, bd.y, bd.r, cols[bd.b][0], cols[bd.b][1], cols[bd.b][2]);
        }

        cv.write_rgb24(stdout);
    }
    return 0;
}
