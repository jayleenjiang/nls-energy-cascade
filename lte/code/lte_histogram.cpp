// lte_histogram.cpp

//   * per-thread private histograms, reduced at the end (no locking)
//   * xoshiro256++ RNG + Box-Muller Gaussian (cached, no std overhead)
//
// Build:
//   clang++ -O3 -mcpu=native -ffast-math -std=c++17 \
//       -I/opt/homebrew/include/eigen3 \
//       -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
//       -L/opt/homebrew/opt/libomp/lib -lomp \
//       lte_histogram.cpp -o lte_histogram
//
// Usage:
//   ./lte_histogram <T1> <Tn> <n> <out_prefix> [site1 site2 ...]
// Example (middle + two others, n=25):
//   ./lte_histogram 10 2 25 histo_N25 6 12 18
//
// Output: one file <out_prefix>_jSITE.hist per requested site, plus
//         <out_prefix>_profile.txt with <I_j> for all modes.

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <chrono>
#include <string>
#include <Eigen/Dense>
#include <omp.h>

// ===================== simulation parameters =====================
const double gamma_val = 0.1;
const double T_final   = 200000.0;   // long run for ~1e9 samples
const double T_burnin  = 2000.0;
const double dt        = 0.001;
const int    dump_every_steps = 100; // keep samples ~decorrelated
const int    N_traj    = 256;        // many trajectories
const int    N_thread  = 8;          // set to your core count

// ===================== histogram parameters ======================
const int    NB      = 80;           // bins per axis (I_a, I_b, theta)
const double I_LO    = 0.0;
const double I_HI    = 4.0;           // I range; samples above go to overflow
const double TH_LO   = -M_PI;
const double TH_HI   =  M_PI;

const double PI = 3.14159265358979323846;
using Vector = Eigen::VectorXd;

// ============================================================
// xoshiro256++  --  fast PRNG
// ============================================================
struct Xoshiro {
    uint64_t s[4];

    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    void seed(uint64_t seed_val) {
        uint64_t x = seed_val;
        for (int i = 0; i < 4; ++i) s[i] = splitmix64(x);
    }
    static inline uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    inline uint64_t next() {
        uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1];
        s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }
    // uniform double in [0,1)
    inline double uniform() {
        return (next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

// ============================================================
// Gaussian sampler: Box-Muller with one cached value.
// (Simple, fast enough; avoids std::normal_distribution overhead.)
// With -O3 -ffast-math the compiler fuses the adjacent sin/cos
// of the same argument into a single sincos.
// ============================================================
struct Gaussian {
    Xoshiro rng;
    bool have_cached = false;
    double cached = 0.0;

    void seed(uint64_t s) { rng.seed(s); have_cached = false; }

    inline double sample() {
        if (have_cached) { have_cached = false; return cached; }
        // draw two uniforms, avoid log(0)
        double u1 = rng.uniform();
        double u2 = rng.uniform();
        if (u1 < 1e-300) u1 = 1e-300;
        double r   = std::sqrt(-2.0 * std::log(u1));
        double ang = 2.0 * PI * u2;
        double s = std::sin(ang);
        double c = std::cos(ang);
        cached = r * s;
        have_cached = true;
        return r * c;
    }
};

// ============================================================
// 3D histogram for one site pair (I_a, I_b, theta)
// ============================================================
struct Hist3D {
    std::vector<uint64_t> bins;   // NB*NB*NB
    uint64_t overflow = 0;        // samples outside [I_LO,I_HI] in I
    uint64_t total    = 0;

    Hist3D() : bins(static_cast<size_t>(NB)*NB*NB, 0) {}

    inline void add(double Ia, double Ib, double th) {
        ++total;
        if (Ia < I_LO || Ia >= I_HI || Ib < I_LO || Ib >= I_HI) {
            ++overflow; return;
        }
        // theta should already be in [-pi, pi)
        int ia = int((Ia - I_LO) / (I_HI - I_LO) * NB);
        int ib = int((Ib - I_LO) / (I_HI - I_LO) * NB);
        int it = int((th - TH_LO) / (TH_HI - TH_LO) * NB);
        if (ia < 0) ia = 0; if (ia >= NB) ia = NB-1;
        if (ib < 0) ib = 0; if (ib >= NB) ib = NB-1;
        if (it < 0) it = 0; if (it >= NB) it = NB-1;
        ++bins[(size_t(ia)*NB + ib)*NB + it];
    }

    void merge(const Hist3D& o) {
        for (size_t k = 0; k < bins.size(); ++k) bins[k] += o.bins[k];
        overflow += o.overflow;
        total    += o.total;
    }
};

// ============================================================
// One trajectory: evolve SDE, accumulate into per-thread histograms.
// site_list gives the j values; we histogram (I_j, I_{j+1}, theta_j).
// ============================================================
void run_one_trajectory(int n, double T1, double Tn,
                         const std::vector<int>& site_list,
                         uint64_t seed_val,
                         std::vector<Hist3D>& hists,   // one per site
                         Vector& energy_accum,         // sum of I over samples
                         uint64_t& sample_count)
{
    Gaussian g; g.seed(seed_val);

    Vector I   = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);

    double current_time = 0.0;
    long step = 0;
    bool measuring = false;

    const double sqrt_dt = std::sqrt(dt);

    while (current_time < T_final) {
        Vector I_padded = Vector::Zero(n + 2);
        I_padded.segment(1, n) = I;
        Vector I_prev = I_padded.head(n);
        Vector I_next = I_padded.tail(n);

        Vector phi_padded = Vector::Zero(n + 2);
        phi_padded.segment(1, n) = phi;
        Vector phi_prev = phi_padded.head(n);
        Vector phi_next = phi_padded.tail(n);

        Vector d_phi_prev = 2.0 * (phi - phi_prev);
        Vector d_phi_next = 2.0 * (phi - phi_next);

        Vector drift_I = 4.0 * I.array() *
            (I_prev.array() * d_phi_prev.array().sin() +
             I_next.array() * d_phi_next.array().sin());

        double total_mass_I = I.sum();

        Vector drift_phi = (Vector::Constant(n, 2.0 * total_mass_I) - I).array()
            + 2.0 * I_prev.array() * d_phi_prev.array().cos()
            + 2.0 * I_next.array() * d_phi_next.array().cos();

        // boundary terms (Case 0, bug-fixed)
        drift_I(0)   += 2.0 * gamma_val * (2.0 * T1 -
            (2.0 * total_mass_I * I(0) - std::pow(I(0), 2)
             + 2.0 * I_next(0) * I(0) * std::cos(d_phi_next(0))));
        drift_phi(0) += gamma_val * (2.0 * I_next(0) * std::sin(d_phi_next(0)));

        drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn -
            (2.0 * total_mass_I * I(n-1) - std::pow(I(n-1), 2)
             + 2.0 * I_prev(n-1) * I(n-1) * std::cos(d_phi_prev(n-1))));
        drift_phi(n-1) += gamma_val * (2.0 * I_prev(n-1) * std::sin(d_phi_prev(n-1)));

        // noise (bug-fixed sqrt(2))
        double nI0  = 2.0 * std::sqrt(2.0 * gamma_val * T1 * I(0))     * g.sample();
        double nIn  = 2.0 * std::sqrt(2.0 * gamma_val * Tn * I(n-1))   * g.sample();
        double nph0 =       std::sqrt(2.0 * gamma_val * T1 / I(0))     * g.sample();
        double nphn =       std::sqrt(2.0 * gamma_val * Tn / I(n-1))   * g.sample();

        I   += drift_I   * dt;
        phi += drift_phi * dt;
        I(0)   += nI0  * sqrt_dt;
        I(n-1) += nIn  * sqrt_dt;
        phi(0)   += nph0 * sqrt_dt;
        phi(n-1) += nphn * sqrt_dt;

        I = I.cwiseMax(1e-10);

        current_time += dt;
        ++step;

        if (!measuring && current_time >= T_burnin) measuring = true;

        if (measuring && (step % dump_every_steps == 0)) {
            energy_accum += I;
            ++sample_count;
            for (size_t s = 0; s < site_list.size(); ++s) {
                int j = site_list[s];
                double Ia = I(j);
                double Ib = I(j+1);
                double th = 2.0 * (phi(j+1) - phi(j));
                // wrap theta to [-pi, pi)
                th = std::fmod(th + PI, 2.0*PI);
                if (th < 0) th += 2.0*PI;
                th -= PI;
                hists[s].add(Ia, Ib, th);
            }
        }
    }
}

// ============================================================
void write_histogram(const std::string& fname, int j,
                     const Hist3D& h, double T1, double Tn, int n)
{
    std::ofstream f(fname);
    f << "# LTE 3D histogram for site pair (" << j << "," << j+1 << ")\n";
    f << "# (I_" << j << ", I_" << j+1 << ", theta_" << j << ")\n";
    f << "# T1=" << T1 << " Tn=" << Tn << " n=" << n << "\n";
    f << "# NB=" << NB << "\n";
    f << "# I_LO=" << I_LO << " I_HI=" << I_HI << "\n";
    f << "# TH_LO=" << TH_LO << " TH_HI=" << TH_HI << "\n";
    f << "# total_samples=" << h.total << " overflow=" << h.overflow << "\n";
    f << "# format: ia ib it count   (only nonzero bins)\n";
    f << "NB " << NB << "\n";
    f << "I_LO " << I_LO << "\n";
    f << "I_HI " << I_HI << "\n";
    f << "TH_LO " << TH_LO << "\n";
    f << "TH_HI " << TH_HI << "\n";
    f << "TOTAL " << h.total << "\n";
    f << "OVERFLOW " << h.overflow << "\n";
    for (int ia = 0; ia < NB; ++ia)
      for (int ib = 0; ib < NB; ++ib)
        for (int it = 0; it < NB; ++it) {
            uint64_t c = h.bins[(size_t(ia)*NB + ib)*NB + it];
            if (c) f << ia << " " << ib << " " << it << " " << c << "\n";
        }
    f.close();
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <T1> <Tn> <n> <out_prefix> [site1 site2 ...]\n"
                  << "Example: " << argv[0] << " 10 2 25 histo_N25 6 12 18\n";
        return 1;
    }
    double T1 = std::atof(argv[1]);
    double Tn = std::atof(argv[2]);
    int    n  = std::atoi(argv[3]);
    std::string prefix = argv[4];

    std::vector<int> site_list;
    for (int a = 5; a < argc; ++a) site_list.push_back(std::atoi(argv[a]));
    if (site_list.empty()) site_list.push_back(n/2);  // default: middle

    for (int j : site_list) {
        if (j < 0 || j+1 >= n) {
            std::cerr << "Bad site " << j << " (need 0 <= j, j+1 < n)\n";
            return 1;
        }
    }

    long samples_per_traj = long((T_final - T_burnin) / (dump_every_steps * dt));
    long total_expected   = (long)N_traj * samples_per_traj;

    std::cout << "LTE histogram accumulator (Case 0 closed chain)\n"
              << "T1=" << T1 << " Tn=" << Tn << " n=" << n << "\n"
              << "sites: ";
    for (int j : site_list) std::cout << "(" << j << "," << j+1 << ") ";
    std::cout << "\n"
              << "T_final=" << T_final << " T_burnin=" << T_burnin
              << " dt=" << dt << " dump_every=" << dump_every_steps << "\n"
              << "N_traj=" << N_traj << " NB=" << NB << "\n"
              << "samples/traj=" << samples_per_traj
              << "  total expected=" << total_expected
              << " (~" << (double)total_expected/1e9 << "e9)\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    int n_sites = site_list.size();
    // global (merged) histograms + energy profile
    std::vector<Hist3D> global_hists(n_sites);
    Vector global_energy = Vector::Zero(n);
    uint64_t global_samples = 0;

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();

        // per-thread private histograms + energy accumulator
        std::vector<Hist3D> local_hists(n_sites);
        Vector local_energy = Vector::Zero(n);
        uint64_t local_samples = 0;

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N_traj; ++t) {
            // unique seed per trajectory
            uint64_t seed_val = 0x12345ULL
                + 0x9E3779B97F4A7C15ULL * (uint64_t)t
                + 0xD1B54A32D192ED03ULL * (uint64_t)rank;
            run_one_trajectory(n, T1, Tn, site_list, seed_val,
                               local_hists, local_energy, local_samples);
            if (rank == 0) {
                #pragma omp critical
                std::cout << "  trajectory " << t << " done\n";
            }
        }

        // reduce into global
        #pragma omp critical
        {
            for (int s = 0; s < n_sites; ++s)
                global_hists[s].merge(local_hists[s]);
            global_energy += local_energy;
            global_samples += local_samples;
        }
    }

    // ---- write outputs ----
    for (int s = 0; s < n_sites; ++s) {
        std::ostringstream fn;
        fn << prefix << "_j" << site_list[s] << ".hist";
        write_histogram(fn.str(), site_list[s], global_hists[s], T1, Tn, n);
        std::cout << "Wrote " << fn.str()
                  << "  (" << global_hists[s].total << " samples, "
                  << global_hists[s].overflow << " overflow)\n";
    }

    // energy profile
    {
        std::ostringstream fn;
        fn << prefix << "_profile.txt";
        std::ofstream f(fn.str());
        f << "# mode_index   <I_j>\n";
        f << "# total_samples=" << global_samples << "\n";
        for (int j = 0; j < n; ++j)
            f << j << " " << global_energy(j) / (double)global_samples << "\n";
        f.close();
        std::cout << "Wrote " << fn.str() << "\n";
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
    std::cout << "Total samples: " << global_samples
              << " (~" << (double)global_samples/1e9 << "e9)\n";
    std::cout << "Done in " << sec << "s.\n";
    return 0;
}
