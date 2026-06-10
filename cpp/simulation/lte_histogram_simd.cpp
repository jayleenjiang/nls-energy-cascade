// lte_histogram_simd.cpp
//
// SIMD-over-TRAJECTORIES rewrite of lte_histogram.cpp.
//
// Key idea (answers Prof. Li's Q1): the previous version vectorized Eigen
// array ops over the n=25 MODES, but (a) heap-allocated ~10 length-n vectors
// EVERY time step (2e8 steps -> billions of allocs), (b) sin/cos went through
// scalar libm, (c) n=25 is too short for good SIMD. The CPU was NOT doing
// 8/16 lanes per cycle.
//
// Here we instead run W independent chains in lockstep (structure-of-arrays):
// state is I[j*W + lane], phi[j*W + lane]. The inner loop over `lane` is
// marked `#pragma omp simd`, so each arithmetic op processes W trajectories
// at once (W=16 floats -> one AVX-512 op, or two AVX2 ops). No per-step
// allocation. fast_sin/fast_cos replace libm in the hot path (Q2).
//
// Q3: 1D marginal histograms P(I_a), P(I_b), P(theta) are accumulated in C++
//     alongside the 3D histogram and written to a .marg file (finer bins).
// Q4: defaults give ~1.0e9 samples (N_GROUPS*W*samples_per_traj).
//
// Bug fixes carried over (unchanged physics vs lte_histogram.cpp):
//   (1) boundary phi-drift sign (+=),  (2) boundary noise sqrt(2) factor,
//   Case-0 closed-chain Gibbs bath.
//
// Float safety: phi is WRAPPED to [-pi,pi) every step. The drift depends on
// phi only through cos/sin of 2*(phi_j - phi_{j±1}), which are 2*pi-periodic,
// so per-mode wrapping leaves the dynamics exactly invariant while keeping phi
// O(1) (the old code let phi grow ~1e6 over T=2e5 -> float would lose phase).
//
// Build (Linux):
//   g++ -O3 -march=native -ffast-math -fopenmp lte_histogram_simd.cpp -o lte_histogram_simd
// Build (mac, Apple Silicon, clang):
//   clang++ -O3 -mcpu=native -ffast-math -std=c++17 \
//       -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
//       -L/opt/homebrew/opt/libomp/lib -lomp \
//       lte_histogram_simd.cpp -o lte_histogram_simd
//
// Usage (same as before, sites at the end):
//   ./lte_histogram_simd <T1> <Tn> <n> <out_prefix> [site1 site2 ...]
//   ./lte_histogram_simd 10 2 25 histo_N25 6 12 18

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <chrono>
#include <string>
#include <omp.h>

// ===================== simulation parameters =====================
static const double gamma_val = 0.1;
static const double T_final    = 200000.0;   // long run
static const double T_burnin   = 2000.0;
static const double dt         = 0.001;
static const int    dump_every_steps = 100;  // keep samples ~decorrelated
static const int    N_GROUPS   = 32;          // groups of W chains  (32*16=512 chains)
static const int    N_thread   = 8;

// ===================== SIMD width ================================
// W=16 floats -> AVX-512 single op (or 2x AVX2). This is the "16 entries
// in 1 cycle" the question asks about; the inner lane loop vectorizes.
static const int    W = 16;

// ===================== 3D histogram params =======================
static const int    NB    = 80;     // bins per axis for the joint 3D hist
static const float  I_LO  = 0.0f;
static const float  I_HI  = 4.0f;
static const float  TH_LO = -3.14159265358979323846f;
static const float  TH_HI =  3.14159265358979323846f;

// ===================== 1D marginal params ========================
// finer 1D bins (the 3D is memory-limited to NB^3; marginals can afford more)
static const int    NB1   = 400;

static const float  PI_f  = 3.14159265358979323846f;
static const double PI    = 3.14159265358979323846;

// ============================================================
// fast, branchless transcendentals (vectorize cleanly under omp simd)
// input to fast_sin must be in [-pi, pi]
// ============================================================
static inline float wrap_pi(float x) {
    const float inv2pi = 0.159154943091895f, twopi = 6.28318530717959f;
    return x - twopi * std::rint(x * inv2pi);
}
static inline float fast_sin(float x) {            // Bhaskara-style, x in [-pi,pi]
    const float B = 1.27323954f, C = -0.40528473f, P = 0.225f;
    float y = B * x + C * x * std::fabs(x);
    return P * (y * std::fabs(y) - y) + y;
}
static inline float fast_cos(float x) {
    return fast_sin(wrap_pi(x + 1.5707963268f));
}
// fast natural log (fastapprox-style), x>0; ~1e-3 accuracy is plenty for noise
static inline float fast_log(float x) {
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFFu) | 0x3f000000u };
    float y = (float)vx.i * 1.1920928955078125e-7f;
    return 0.69314718f * (y - 124.22551499f
           - 1.498030302f * mx.f - 1.72587999f / (0.3520887068f + mx.f));
}

// ============================================================
// xoshiro256++  (fast scalar PRNG)
// ============================================================
struct Xoshiro {
    uint64_t s[4];
    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    void seed(uint64_t v) { uint64_t x = v; for (int i = 0; i < 4; ++i) s[i] = splitmix64(x); }
    static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    inline uint64_t next() {
        uint64_t r = rotl(s[0] + s[3], 23) + s[0];
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return r;
    }
    inline float uniform() { return (float)((next() >> 11) * (1.0 / 9007199254740992.0)); }
};

// fill `out[W]` with W i.i.d. standard normals (Box-Muller, fast_log + fast_sin)
static inline void fill_gauss(float* out, Xoshiro& rng) {
    for (int l = 0; l < W; l += 2) {
        float u1 = rng.uniform(); if (u1 < 1e-30f) u1 = 1e-30f;
        float u2 = rng.uniform();
        float r  = std::sqrt(-2.0f * fast_log(u1));
        float a  = wrap_pi(6.28318530717959f * u2 - 3.14159265358979f);
        out[l]   = r * fast_sin(a);
        if (l + 1 < W) out[l + 1] = r * fast_cos(a);
    }
}

// ============================================================
// 3D histogram (I_a, I_b, theta) for one site pair
// ============================================================
struct Hist3D {
    std::vector<uint64_t> bins;
    uint64_t overflow = 0, total = 0;
    Hist3D() : bins((size_t)NB * NB * NB, 0) {}
    inline void add(float Ia, float Ib, float th) {
        ++total;
        if (Ia < I_LO || Ia >= I_HI || Ib < I_LO || Ib >= I_HI) { ++overflow; return; }
        int ia = int((Ia - I_LO) / (I_HI - I_LO) * NB);
        int ib = int((Ib - I_LO) / (I_HI - I_LO) * NB);
        int it = int((th - TH_LO) / (TH_HI - TH_LO) * NB);
        if (ia < 0) ia = 0; if (ia >= NB) ia = NB - 1;
        if (ib < 0) ib = 0; if (ib >= NB) ib = NB - 1;
        if (it < 0) it = 0; if (it >= NB) it = NB - 1;
        ++bins[((size_t)ia * NB + ib) * NB + it];
    }
    void merge(const Hist3D& o) {
        for (size_t k = 0; k < bins.size(); ++k) bins[k] += o.bins[k];
        overflow += o.overflow; total += o.total;
    }
};

// ============================================================
// 1D marginals P(I_a), P(I_b), P(theta) for one site pair
// (finer NB1 bins; I-overflow counted separately)
// ============================================================
struct Marg1D {
    std::vector<uint64_t> Ia, Ib, th;
    uint64_t Ia_of = 0, Ib_of = 0, total = 0;
    Marg1D() : Ia(NB1, 0), Ib(NB1, 0), th(NB1, 0) {}
    inline void add(float a, float b, float t) {
        ++total;
        // I marginals share [I_LO, I_HI]
        if (a >= I_LO && a < I_HI) { int k = int((a - I_LO) / (I_HI - I_LO) * NB1);
            if (k < 0) k = 0; if (k >= NB1) k = NB1 - 1; ++Ia[k]; } else ++Ia_of;
        if (b >= I_LO && b < I_HI) { int k = int((b - I_LO) / (I_HI - I_LO) * NB1);
            if (k < 0) k = 0; if (k >= NB1) k = NB1 - 1; ++Ib[k]; } else ++Ib_of;
        // theta always in [-pi,pi)
        int kt = int((t - TH_LO) / (TH_HI - TH_LO) * NB1);
        if (kt < 0) kt = 0; if (kt >= NB1) kt = NB1 - 1; ++th[kt];
    }
    void merge(const Marg1D& o) {
        for (int k = 0; k < NB1; ++k) { Ia[k] += o.Ia[k]; Ib[k] += o.Ib[k]; th[k] += o.th[k]; }
        Ia_of += o.Ia_of; Ib_of += o.Ib_of; total += o.total;
    }
};

// ============================================================
// One GROUP = W independent chains evolved in lockstep (SoA + omp simd).
// Accumulates into thread-local hists/margs/energy.
// ============================================================
static void run_group(int n, double T1d, double Tnd,
                      const std::vector<int>& sites, uint64_t seed_val,
                      std::vector<Hist3D>& hists, std::vector<Marg1D>& margs,
                      std::vector<double>& energy_accum, uint64_t& sample_count)
{
    const float  g    = (float)gamma_val;
    const float  T1   = (float)T1d, Tn = (float)Tnd;
    const float  dtf  = (float)dt, sdt = std::sqrt((float)dt);
    const int    NS   = (int)sites.size();

    // state: I[j*W+l], phi[j*W+l]   (allocated ONCE, not per step)
    std::vector<float> I((size_t)n * W), phi((size_t)n * W);
    std::vector<float> dI((size_t)n * W), dph((size_t)n * W);
    for (int j = 0; j < n; ++j)
        for (int l = 0; l < W; ++l) { I[j*W+l] = (j==0 ? 1.0f : 0.1f); phi[j*W+l] = 0.0f; }

    float M[W];
    float nI0[W], nIn[W], nph0[W], nphn[W];
    Xoshiro rng; rng.seed(seed_val);

    long step = 0; double t = 0.0; bool measuring = false;
    while (t < T_final) {
        // --- per-lane total mass M = sum_j I[j] ---
        #pragma omp simd
        for (int l = 0; l < W; ++l) M[l] = 0.0f;
        for (int j = 0; j < n; ++j) {
            #pragma omp simd
            for (int l = 0; l < W; ++l) M[l] += I[j*W+l];
        }

        // --- bulk drift for every mode (inner loop over lanes = SIMD) ---
        for (int j = 0; j < n; ++j) {
            const bool hasL = (j > 0), hasR = (j < n - 1);
            const int jm = (j-1)*W, jp = (j+1)*W, jj = j*W;
            #pragma omp simd
            for (int l = 0; l < W; ++l) {
                float Ij = I[jj+l], pj = phi[jj+l];
                float Iprev = hasL ? I[jm+l] : 0.0f;
                float Inext = hasR ? I[jp+l] : 0.0f;
                float pprev = hasL ? phi[jm+l] : 0.0f;
                float pnext = hasR ? phi[jp+l] : 0.0f;
                float dprev = wrap_pi(2.0f * (pj - pprev));
                float dnext = wrap_pi(2.0f * (pj - pnext));
                float sprev = fast_sin(dprev), cprev = fast_cos(dprev);
                float snext = fast_sin(dnext), cnext = fast_cos(dnext);
                dI[jj+l]  = 4.0f * Ij * (Iprev * sprev + Inext * snext);
                dph[jj+l] = 2.0f * M[l] - Ij + 2.0f * Iprev * cprev + 2.0f * Inext * cnext;
            }
        }

        // --- Case-0 boundary bath drift (j=0 and j=n-1) ---
        {
            const int j0 = 0, j1 = 1*W;
            #pragma omp simd
            for (int l = 0; l < W; ++l) {
                float I0 = I[j0+l], I1 = I[j1+l];
                float dnext = wrap_pi(2.0f * (phi[j0+l] - phi[j1+l]));
                dI[j0+l]  += 2.0f * g * (2.0f * T1 -
                              (2.0f * M[l] * I0 - I0 * I0 + 2.0f * I1 * I0 * fast_cos(dnext)));
                dph[j0+l] += g * (2.0f * I1 * fast_sin(dnext));
            }
            const int jn = (n-1)*W, jnm = (n-2)*W;
            #pragma omp simd
            for (int l = 0; l < W; ++l) {
                float In = I[jn+l], Inm = I[jnm+l];
                float dprev = wrap_pi(2.0f * (phi[jn+l] - phi[jnm+l]));
                dI[jn+l]  += 2.0f * g * (2.0f * Tn -
                              (2.0f * M[l] * In - In * In + 2.0f * Inm * In * fast_cos(dprev)));
                dph[jn+l] += g * (2.0f * Inm * fast_sin(dprev));
            }
        }

        // --- boundary noise (sqrt(2) factor; sigma from CURRENT I) ---
        fill_gauss(nI0, rng); fill_gauss(nIn, rng);
        fill_gauss(nph0, rng); fill_gauss(nphn, rng);

        // save start-of-step boundary I for EM noise sigma (Ito: evaluate at step start)
        float I0s[W], Ins[W];
        #pragma omp simd
        for (int l = 0; l < W; ++l) { I0s[l] = I[l]; Ins[l] = I[(n-1)*W+l]; }

        // --- update all modes (in place; dI/dph already computed) ---
        for (int j = 0; j < n; ++j) {
            const int jj = j*W;
            #pragma omp simd
            for (int l = 0; l < W; ++l) {
                float Inew = I[jj+l] + dI[jj+l] * dtf;
                phi[jj+l]  = wrap_pi(phi[jj+l] + dph[jj+l] * dtf);
                I[jj+l]    = Inew < 1e-10f ? 1e-10f : Inew;
            }
        }
        // add boundary noise (sigma uses start-of-step I -> exact Euler-Maruyama)
        {
            const int j0 = 0, jn = (n-1)*W;
            #pragma omp simd
            for (int l = 0; l < W; ++l) {
                float I0 = I0s[l] > 1e-10f ? I0s[l] : 1e-10f;
                float In = Ins[l] > 1e-10f ? Ins[l] : 1e-10f;
                float a0 = 2.0f * std::sqrt(2.0f * g * T1 * I0) * sdt * nI0[l];
                float an = 2.0f * std::sqrt(2.0f * g * Tn * In) * sdt * nIn[l];
                float p0 =        std::sqrt(2.0f * g * T1 / I0) * sdt * nph0[l];
                float pn =        std::sqrt(2.0f * g * Tn / In) * sdt * nphn[l];
                float v0 = I[j0+l] + a0; I[j0+l] = v0 < 1e-10f ? 1e-10f : v0;
                float vn = I[jn+l] + an; I[jn+l] = vn < 1e-10f ? 1e-10f : vn;
                phi[j0+l] = wrap_pi(phi[j0+l] + p0);
                phi[jn+l] = wrap_pi(phi[jn+l] + pn);
            }
        }

        t += dt; ++step;
        if (!measuring && t >= T_burnin) measuring = true;

        if (measuring && (step % dump_every_steps == 0)) {
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int l = 0; l < W; ++l) s += I[j*W+l];
                energy_accum[j] += s;
            }
            sample_count += W;
            for (int sIdx = 0; sIdx < NS; ++sIdx) {
                int j = sites[sIdx];
                for (int l = 0; l < W; ++l) {
                    float Ia = I[j*W+l], Ib = I[(j+1)*W+l];
                    float th = wrap_pi(2.0f * (phi[(j+1)*W+l] - phi[j*W+l]));
                    hists[sIdx].add(Ia, Ib, th);
                    margs[sIdx].add(Ia, Ib, th);
                }
            }
        }
    }
}

// ============================================================
static void write_hist(const std::string& fn, int j, const Hist3D& h,
                       double T1, double Tn, int n) {
    std::ofstream f(fn);
    f << "# LTE 3D histogram (I_" << j << ", I_" << j+1 << ", theta_" << j << ")\n";
    f << "# T1=" << T1 << " Tn=" << Tn << " n=" << n << "\n";
    f << "NB " << NB << "\nI_LO " << I_LO << "\nI_HI " << I_HI
      << "\nTH_LO " << TH_LO << "\nTH_HI " << TH_HI << "\n";
    f << "TOTAL " << h.total << "\nOVERFLOW " << h.overflow << "\n";
    f << "# format: ia ib it count (nonzero only)\n";
    for (int ia = 0; ia < NB; ++ia)
      for (int ib = 0; ib < NB; ++ib)
        for (int it = 0; it < NB; ++it) {
            uint64_t c = h.bins[((size_t)ia*NB+ib)*NB+it];
            if (c) f << ia << " " << ib << " " << it << " " << c << "\n";
        }
}

static void write_marg(const std::string& fn, int j, const Marg1D& m,
                       double T1, double Tn, int n) {
    std::ofstream f(fn);
    f << "# LTE 1D marginals for site pair (" << j << "," << j+1 << ")\n";
    f << "# T1=" << T1 << " Tn=" << Tn << " n=" << n << "\n";
    f << "NB1 " << NB1 << "\nI_LO " << I_LO << "\nI_HI " << I_HI
      << "\nTH_LO " << TH_LO << "\nTH_HI " << TH_HI << "\n";
    f << "TOTAL " << m.total << " IA_OVERFLOW " << m.Ia_of
      << " IB_OVERFLOW " << m.Ib_of << "\n";
    f << "# format: k  count_Ia  count_Ib  count_theta\n";
    for (int k = 0; k < NB1; ++k)
        f << k << " " << m.Ia[k] << " " << m.Ib[k] << " " << m.th[k] << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <T1> <Tn> <n> <out_prefix> [site1 site2 ...]\n"
                  << "Example: " << argv[0] << " 10 2 25 histo_N25 6 12 18\n";
        return 1;
    }
    double T1 = std::atof(argv[1]), Tn = std::atof(argv[2]);
    int    n  = std::atoi(argv[3]);
    std::string prefix = argv[4];
    std::vector<int> sites;
    for (int a = 5; a < argc; ++a) sites.push_back(std::atoi(argv[a]));
    if (sites.empty()) sites.push_back(n / 2);
    for (int j : sites) if (j < 0 || j + 1 >= n) {
        std::cerr << "Bad site " << j << " (need 0<=j, j+1<n)\n"; return 1; }

    long spt = long((T_final - T_burnin) / (dump_every_steps * dt));   // per chain
    double total_expected = (double)N_GROUPS * W * spt;
    std::cout << "LTE SIMD histogram (Case 0 closed chain)\n"
              << "T1=" << T1 << " Tn=" << Tn << " n=" << n
              << "  W=" << W << " lanes  N_GROUPS=" << N_GROUPS
              << " -> " << N_GROUPS * W << " chains\n"
              << "samples/chain=" << spt
              << "  total expected=" << total_expected
              << " (~" << total_expected / 1e9 << "e9)\nsites:";
    for (int j : sites) std::cout << " (" << j << "," << j+1 << ")";
    std::cout << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    const int NS = (int)sites.size();
    std::vector<Hist3D> g_h(NS);
    std::vector<Marg1D> g_m(NS);
    std::vector<double> g_E(n, 0.0);
    uint64_t g_N = 0;

    #pragma omp parallel num_threads(N_thread)
    {
        std::vector<Hist3D> lh(NS);
        std::vector<Marg1D> lm(NS);
        std::vector<double> lE(n, 0.0);
        uint64_t lN = 0;

        #pragma omp for schedule(dynamic)
        for (int grp = 0; grp < N_GROUPS; ++grp) {
            uint64_t seed = 0xC0FFEEULL
                + 0x9E3779B97F4A7C15ULL * (uint64_t)grp;
            run_group(n, T1, Tn, sites, seed, lh, lm, lE, lN);
            #pragma omp critical
            std::cout << "  group " << grp << " done\n";
        }
        #pragma omp critical
        {
            for (int s = 0; s < NS; ++s) { g_h[s].merge(lh[s]); g_m[s].merge(lm[s]); }
            for (int j = 0; j < n; ++j) g_E[j] += lE[j];
            g_N += lN;
        }
    }

    for (int s = 0; s < NS; ++s) {
        std::ostringstream a, b;
        a << prefix << "_j" << sites[s] << ".hist";
        b << prefix << "_j" << sites[s] << ".marg";
        write_hist(a.str(), sites[s], g_h[s], T1, Tn, n);
        write_marg(b.str(), sites[s], g_m[s], T1, Tn, n);
        std::cout << "Wrote " << a.str() << " and " << b.str()
                  << "  (" << g_h[s].total << " samples, "
                  << g_h[s].overflow << " I-overflow)\n";
    }
    {
        std::ofstream f((prefix + "_profile.txt").c_str());
        f << "# mode_index   <I_j>\n# total_samples=" << g_N << "\n";
        for (int j = 0; j < n; ++j) f << j << " " << g_E[j] / (double)g_N << "\n";
        std::cout << "Wrote " << prefix << "_profile.txt\n";
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
    std::cout << "Total samples: " << g_N << " (~" << (double)g_N/1e9 << "e9)\n"
              << "Done in " << sec << "s.\n";
    return 0;
}
