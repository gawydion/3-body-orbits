#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

struct OrbitInput {
    int dim;
    double m[3];
    double y0[18];
    double T_ghost; // full period for ghost track
    double T_anim;  // period for 1 seamless animation loop
    int frames;
    double yaw0, yaw_span, pitch;
};

static inline void deriv(const double* s, const double m[3], int dim, double* ds) {
    int half = (dim == 3) ? 9 : 6;
    for (int i = 0; i < half; ++i) ds[i] = s[half + i];
    for (int i = half; i < 2 * half; ++i) ds[i] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
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

static void advance(double* y, double dt, const double m[3], int dim, double tol = 1e-11) {
    static constexpr double a21=1.0/5.0, a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        b1=35.0/384.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
    int n = (dim == 3) ? 18 : 12;
    double k1[18], k2[18], k3[18], k4[18], k5[18], k6[18], k7[18], yt[18], yn[18];
    double t = 0.0, h = std::min(0.005, dt);
    int steps = 0;
    while (t < dt - 1e-14 && steps < 20000) {
        if (t + h > dt) h = dt - t;
        deriv(y, m, dim, k1);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*a21*k1[i]; deriv(yt, m, dim, k2);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a31*k1[i]+a32*k2[i]); deriv(yt, m, dim, k3);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(yt, m, dim, k4);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(yt, m, dim, k5);
        for (int i=0;i<n;++i) yt[i]=y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(yt, m, dim, k6);
        for (int i=0;i<n;++i) yn[i]=y[i]+h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]); deriv(yn, m, dim, k7);
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
    std::vector<float> buf; // RGB float [0,1]

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
        // Soft glow halo
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
        // Solid body with crisp AA edge + white rim
        int ir = (int)std::ceil(rad + 2.0f);
        for (int dy = -ir; dy <= ir; ++dy) {
            for (int dx = -ir; dx <= ir; ++dx) {
                float px = std::round(cx) + dx, py = std::round(cy) + dy;
                float d = std::hypot(px - cx, py - cy);
                if (d <= rad + 1.2f) {
                    float cov = std::clamp(rad + 0.7f - d, 0.0f, 1.0f);
                    // Subtle specular highlight
                    float hl = std::max(0.0f, 1.0f - std::hypot(dx + rad*0.28f, dy + rad*0.28f) / (rad * 0.85f));
                    float rr = std::min(1.0f, r + hl * 0.45f);
                    float gg = std::min(1.0f, g + hl * 0.45f);
                    float bb = std::min(1.0f, b + hl * 0.45f);
                    blend((int)px, (int)py, rr, gg, bb, cov);
                    // Thin crisp outer rim
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
    // Reads from stdin:
    // dim W H frames T_ghost T_anim yaw0 yaw_span pitch m1 m2 m3 y0[0..2*3*dim-1]
    int dim, W, H, frames;
    double T_ghost, T_anim, yaw0, yaw_span, pitch;
    double m[3], y0[18] = {0};
    if (!(std::cin >> dim >> W >> H >> frames >> T_ghost >> T_anim >> yaw0 >> yaw_span >> pitch
          >> m[0] >> m[1] >> m[2])) return 1;
    int n = (dim == 3) ? 18 : 12;
    for (int i = 0; i < n; ++i) std::cin >> y0[i];

    // Precompute ghost track (high resolution)
    int N_ghost = (T_ghost > 50.0) ? 2400 : 800;
    std::vector<std::array<double, 9>> ghost(N_ghost + 1);
    double yt[18];
    for (int i = 0; i < n; ++i) yt[i] = y0[i];
    double maxR = 0.5;
    for (int k = 0; k <= N_ghost; ++k) {
        std::array<double, 9> pt = {0};
        for (int b = 0; b < 3; ++b) {
            pt[3*b]   = yt[dim*b];
            pt[3*b+1] = yt[dim*b+1];
            pt[3*b+2] = (dim == 3) ? yt[dim*b+2] : 0.0;
            double r = std::sqrt(pt[3*b]*pt[3*b] + pt[3*b+1]*pt[3*b+1] + pt[3*b+2]*pt[3*b+2]);
            if (r > maxR) maxR = r;
        }
        ghost[k] = pt;
        if (k < N_ghost) advance(yt, T_ghost / N_ghost, m, dim, 1e-11);
    }

    // Precompute dense animation samples (substeps = 4 per frame for silky trails)
    int sub = 4;
    int N_anim = frames * sub;
    std::vector<std::array<double, 9>> anim(N_anim + 1);
    for (int i = 0; i < n; ++i) yt[i] = y0[i];
    for (int k = 0; k <= N_anim; ++k) {
        std::array<double, 9> pt = {0};
        for (int b = 0; b < 3; ++b) {
            pt[3*b]   = yt[dim*b];
            pt[3*b+1] = yt[dim*b+1];
            pt[3*b+2] = (dim == 3) ? yt[dim*b+2] : 0.0;
        }
        anim[k] = pt;
        if (k < N_anim) advance(yt, T_anim / N_anim, m, dim, 1e-11);
    }

    double scale = std::min(W, H) * 0.42 / maxR;
    const float cols[3][3] = {
        {0.231f, 0.510f, 0.965f}, // blue #3b82f6
        {0.961f, 0.620f, 0.043f}, // amber #f59e0b
        {0.925f, 0.282f, 0.600f}  // pink #ec4899
    };

    Canvas cv(W, H);
    int trail_len = std::min(N_anim - 1, (int)(N_anim * 0.42));

    for (int f = 0; f < frames; ++f) {
        // Dark slate background #0b0e14
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

        // Subtle origin crosshair (or 3D triad)
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

        // Faint ghost track
        float ghost_alpha = (T_ghost > 50.0) ? 0.14f : 0.22f;
        for (int b = 0; b < 3; ++b) {
            float px0, py0, pz0;
            proj(ghost[0][3*b], ghost[0][3*b+1], ghost[0][3*b+2], px0, py0, pz0);
            for (int k = 1; k <= N_ghost; ++k) {
                float px1, py1, pz1;
                proj(ghost[k][3*b], ghost[k][3*b+1], ghost[k][3*b+2], px1, py1, pz1);
                cv.draw_line(px0, py0, px1, py1, cols[b][0], cols[b][1], cols[b][2], ghost_alpha, 1.2f);
                px0 = px1; py0 = py1;
            }
        }

        // Live fading trail (wrapping seamlessly for periodic orbits!)
        int cur_idx = f * sub;
        for (int b = 0; b < 3; ++b) {
            for (int step = 0; step < trail_len; ++step) {
                int i0 = (cur_idx - trail_len + step + N_anim) % N_anim;
                int i1 = (i0 + 1) % N_anim;
                // Avoid wrap jump if T_anim != T_ghost (e.g. non-closed RPO)
                double dx = anim[i1][3*b] - anim[i0][3*b];
                double dy = anim[i1][3*b+1] - anim[i0][3*b+1];
                double dz = anim[i1][3*b+2] - anim[i0][3*b+2];
                if (dx*dx + dy*dy + dz*dz > 0.25) continue;
                float frac = (float)(step + 1) / trail_len;
                float alpha = frac * frac * 0.92f;
                float x0, y0p, z0, x1, y1p, z1;
                proj(anim[i0][3*b], anim[i0][3*b+1], anim[i0][3*b+2], x0, y0p, z0);
                proj(anim[i1][3*b], anim[i1][3*b+1], anim[i1][3*b+2], x1, y1p, z1);
                cv.draw_line(x0, y0p, x1, y1p, cols[b][0], cols[b][1], cols[b][2], alpha, 2.2f);
            }
        }

        // Bodies (depth-sorted in 3D)
        struct BDraw { int b; float x, y, z, r; };
        std::array<BDraw, 3> bds;
        for (int b = 0; b < 3; ++b) {
            float sx, sy_out, sz;
            proj(anim[cur_idx][3*b], anim[cur_idx][3*b+1], anim[cur_idx][3*b+2], sx, sy_out, sz);
            float rad = 5.0f + 3.4f * (float)std::cbrt(m[b]);
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
