// sec4_3_return.cpp  --  Section 4.3 (Fig. 19): relaxation of the three-mode chain from a fixed initial state.
//
// WHAT IT COMPUTES
//   Many trajectories of the three-mode chain (n = 3) all started at the same point
//     x0 = (I_1, I_2, I_3, theta_1, theta_3),   theta_1 = 2 (phi_1 - phi_2),  theta_3 = 2 (phi_3 - phi_2),
//   and, at the times t = 0, 0.01, 0.02, ..., the averages over trajectories (i.e. E_x0[g(X_t)]) of
//     * 12 indicator functions: 1 if X_t lies in a ball around one of 6 centres, for 2 radii (0.5 and 1);
//       distance = Euclidean in (I_1, I_2, I_3, theta_1, theta_3), angle differences taken in [-pi, pi);
//     * 6 observables: I_1, I_2, I_3, cos theta_1, cos theta_3, I_1 - I_3.
//   Fig. 19 shows how E_x0[1_ball(X_t)] for the ball around x0 decays and oscillates in t for several gamma.
//   The trajectories are split into 32 groups (batch b goes to group b mod 32), which gives error bars.
//
// MODEL AND SCHEME
//   As in sec4_3_spectrum.cpp: Cartesian amplitudes c_j, implicit midpoint for the Hamiltonian part (fixed frame,
//   tolerance 2e-13, at most 20 iterations), Euler--Maruyama for the canonical baths at modes 1 and 3 (order
//   alternating); batches of 16 trajectories (Vec16), one batch per OpenMP thread at a time.
//   Initial point: phi_2 = 0, phi_1 = theta_1/2, phi_3 = theta_3/2, |c_j|^2 = I_j.
//
// OUTPUT  (default name return_x<I1>_<I2>_<I3>_<theta1>_<theta3>_g<gamma>.json)
//   One JSON file: parameters, x0, the centres and radii of the balls, and "EX": for every group g a list of
//   nlag x 18 sums (time index outer; 12 indicators, then the 6 observables); divide by the number of trajectories
//   of a group (batches * 16 / 32) to get the averages.  Also the number of midpoint failures and of batches that
//   became non-finite.
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 -Ieigen sec4_3_return.cpp -o sec4_3_return
//   ./sec4_3_return --gamma 0.1 --seed 101
//   ./sec4_3_return --gamma 0.01 --seed 107 --x0 2,1,1,0,0
// Options and their defaults (the paper settings):
//   --T1 2  --T3 8  --gamma 0.1  --dt 1e-3  --snapshot 0.01  --nlag 501 (times 0 .. 5)  --batches 4096 (65536 trajectories)
//   --seed 1  --threads 2  --x0 1,1,1,0,0  --out (see above)   --dry-run 1   print the parameters and stop.
//
// The arithmetic is exactly that of the program used for the paper (local_test_v2.cpp, mode "ex"): with the same
// compiler, flags and 1 thread the output is identical bit for bit.  (With more threads, the batches are added to
// their group in the order they finish, so the sums can differ in the last digit.)
//
// LAYOUT OF THIS FILE
//   1. constants and types
//   2. random numbers            Rng, gaussian_pair
//   3. chain state and forces    State, compute_force
//   4. time step                 midpoint_step, bath_step, wrap_angle, finite_state
//   5. measured functions        ball centres and radii, reduced, indicators, set_point
//   6. parameters                Params, read_options, print_parameters
//   7. simulation and output     run
//   8. main

#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {   // internal linkage, as in the original program (this also keeps gcc's rounding identical)

// ============================================================================ 1. constants and types
constexpr int N = 3;                          // number of modes
constexpr int LANES = 16;                     // trajectories per batch
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_ITERATIONS = 20;            // fixed-point iterations of the midpoint rule
constexpr double TOLERANCE = 2.0e-13;         // stop when max |change| < TOLERANCE
constexpr int NGROUPS = 32;                   // groups of batches (for error bars)

using Vec16 = Eigen::Array<double, LANES, 1>;                      // one value per trajectory of a batch
using Modes = std::vector<Vec16, Eigen::aligned_allocator<Vec16>>;  // one Vec16 per mode

// ============================================================================ 2. random numbers
// Rng: Xoshiro256++, one per trajectory, seeded from the base seed and the trajectory number only.
struct Rng {
    std::array<std::uint64_t, 4> s{};
    static std::uint64_t splitmix64(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    void seed(std::uint64_t base_seed, std::uint64_t trajectory) {
        std::uint64_t x = base_seed ^ (0xD2B74407B1CE6E93ULL * (trajectory + 1ULL));
        std::uint64_t value = splitmix64(x);
        for (auto& v : s) v = splitmix64(value);
    }
    std::uint64_t next() {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        return result;
    }
    double uniform() { return ((next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }   // in (0, 1)
};

// gaussian_pair: two independent N(0,1) numbers per trajectory (Box--Muller)
void gaussian_pair(std::array<Rng, LANES>& rng, Vec16& first, Vec16& second) {
    for (int lane = 0; lane < LANES; ++lane) {
        const double u1 = rng[lane].uniform();
        const double u2 = rng[lane].uniform();
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double angle = TWO_PI * u2;
        first(lane) = radius * std::cos(angle);
        second(lane) = radius * std::sin(angle);
    }
}

// ============================================================================ 3. chain state and forces
// State of one batch: x, y = real and imaginary parts of c_j; compute_force fills I = |c|^2, (sq_re, sq_im) = c^2,
// (F_re, F_im) = F and M = I_1 + I_2 + I_3.
struct State {
    Modes x{N}, y{N}, I{N}, sq_re{N}, sq_im{N}, F_re{N}, F_im{N};
    Vec16 M = Vec16::Zero();
};

// compute_force: F_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2) (missing neighbours count 0);
// with s = c_{j-1}^2 + c_{j+1}^2:  conj(c) s = (x s_re + y s_im) + i (x s_im - y s_re).
void compute_force(State& st) {
    st.M.setZero();
    for (int j = 0; j < N; ++j) {
        st.I[j] = st.x[j].square() + st.y[j].square();
        st.sq_re[j] = st.x[j].square() - st.y[j].square();
        st.sq_im[j] = 2.0 * st.x[j] * st.y[j];
        st.M += st.I[j];
    }
    for (int j = 0; j < N; ++j) {
        Vec16 nb_re = Vec16::Zero();
        Vec16 nb_im = Vec16::Zero();
        if (j > 0) { nb_re += st.sq_re[j - 1]; nb_im += st.sq_im[j - 1]; }
        if (j + 1 < N) { nb_re += st.sq_re[j + 1]; nb_im += st.sq_im[j + 1]; }
        const Vec16 onsite = 2.0 * st.M - st.I[j];
        st.F_re[j] = onsite * st.x[j] + 2.0 * (nb_re * st.x[j] + nb_im * st.y[j]);
        st.F_im[j] = onsite * st.y[j] + 2.0 * (nb_im * st.x[j] - nb_re * st.y[j]);
    }
}

// ============================================================================ 4. time step
// MidpointWork: scratch space (the midpoint and the current guess for c').
struct MidpointWork {
    State mid;
    Modes cx{N}, cy{N};
};

// midpoint_step: c' = c + dt i F((c + c')/2) (fixed frame), i.e. x' = x - dt F_im(mid), y' = y + dt F_re(mid),
// by fixed-point iteration from the explicit Euler step until no coordinate changes by TOLERANCE or more.
// At the end st holds c' and its forces.  Returns false if it did not converge in MAX_ITERATIONS iterations.
bool midpoint_step(State& st, MidpointWork& w, double dt, int& iterations) {
    compute_force(st);
    for (int j = 0; j < N; ++j) {
        w.cx[j] = st.x[j] - dt * st.F_im[j];
        w.cy[j] = st.y[j] + dt * st.F_re[j];
    }
    bool converged = false;
    for (iterations = 1; iterations <= MAX_ITERATIONS; ++iterations) {
        for (int j = 0; j < N; ++j) {
            w.mid.x[j] = 0.5 * (st.x[j] + w.cx[j]);
            w.mid.y[j] = 0.5 * (st.y[j] + w.cy[j]);
        }
        compute_force(w.mid);
        double change = 0.0;
        for (int j = 0; j < N; ++j) {
            const Vec16 new_x = st.x[j] - dt * w.mid.F_im[j];
            const Vec16 new_y = st.y[j] + dt * w.mid.F_re[j];
            change = std::max(change, (new_x - w.cx[j]).abs().maxCoeff());
            change = std::max(change, (new_y - w.cy[j]).abs().maxCoeff());
            w.cx[j] = new_x; w.cy[j] = new_y;
        }
        if (change < TOLERANCE) { converged = true; break; }
    }
    for (int j = 0; j < N; ++j) { st.x[j] = w.cx[j]; st.y[j] = w.cy[j]; }
    compute_force(st);
    return converged;
}

// bath_step: Euler--Maruyama step of the canonical bath at mode r,
//   c_r += -gamma F_r dt + sqrt(2 gamma T dt) (N_x + i N_y),  then recompute the forces.
void bath_step(State& st, int r, double T, double gamma, const Vec16& Nx, const Vec16& Ny, double dt) {
    const double noise = std::sqrt(2.0 * gamma * T * dt);
    st.x[r] += -gamma * st.F_re[r] * dt + noise * Nx;
    st.y[r] += -gamma * st.F_im[r] * dt + noise * Ny;
    compute_force(st);
}

// wrap_angle: the angle moved into [-pi, pi)
double wrap_angle(double value) {
    value = std::fmod(value + PI, TWO_PI);
    if (value < 0.0) value += TWO_PI;
    return value - PI;
}
// finite_state: all coordinates of all 16 trajectories are finite
bool finite_state(const State& st) {
    for (int j = 0; j < N; ++j)
        if (!st.x[j].isFinite().all() || !st.y[j].isFinite().all()) return false;
    return true;
}

// ============================================================================ 5. measured functions
constexpr int NC = 6, NR = 2, NF = NC * NR, NOBSX = 6;   // centres, radii, indicator functions, observables
// ball centres (I1, I2, I3, theta1, theta3) and radii
const double CEN[NC][5] = {{1,1,1,0,0},{1,1,1,PI,PI},{2,1,1,0,0},{1,1,2,0,0},{1,1,1,PI/2,-PI/2},{1,2,1,PI,PI}};
const double RAD[NR] = {0.5, 1.0};

// reduced: the reduced coordinates o = (I1, I2, I3, theta1, theta3) of trajectory l
inline void reduced(const State& st, int l, double* o) {
    const double p1 = std::atan2(st.y[0](l), st.x[0](l)), p2 = std::atan2(st.y[1](l), st.x[1](l)),
                 p3 = std::atan2(st.y[2](l), st.x[2](l));
    o[0] = st.I[0](l); o[1] = st.I[1](l); o[2] = st.I[2](l);
    o[3] = wrap_angle(2.0 * (p1 - p2)); o[4] = wrap_angle(2.0 * (p3 - p2));
}
// indicators: g[c * NR + r] = 1 if o lies in the ball of radius RAD[r] around centre c
inline void indicators(const double* o, unsigned char* g) {
    for (int c = 0; c < NC; ++c) {
        double d2 = 0;
        for (int j = 0; j < 3; ++j) d2 += (o[j] - CEN[c][j]) * (o[j] - CEN[c][j]);
        for (int j = 3; j < 5; ++j) { const double a = wrap_angle(o[j] - CEN[c][j]); d2 += a * a; }
        for (int r = 0; r < NR; ++r) g[c * NR + r] = d2 < RAD[r] * RAD[r];
    }
}
// set_point: every trajectory of the batch at x0 (phi_2 = 0, phi_1 = theta1/2, phi_3 = theta3/2, |c_j|^2 = I_j)
void set_point(State& st, const double* x0) {
    const double ph[3] = {0.5 * x0[3], 0.0, 0.5 * x0[4]};
    for (int j = 0; j < N; ++j) { st.x[j].setConstant(std::sqrt(x0[j]) * std::cos(ph[j])); st.y[j].setConstant(std::sqrt(x0[j]) * std::sin(ph[j])); }
}

// ============================================================================ 6. parameters
struct Params {
    double T1 = 2, T3 = 8, gamma = 0.1, dt = 1e-3, snapshot = 0.01;
    int nlag = 501, batches = 4096, threads = 2;
    std::uint64_t seed = 1;
    std::string x0_text = "1,1,1,0,0";            // initial point as given
    double x0[5] = {1, 1, 1, 0, 0};
    std::string out;                              // output file; empty: return_x<..>_g<gamma>.json
    std::string gamma_text = "0.1";               // gamma as given (used in the default file name)
    bool dry_run = false;
};

std::string fmt(double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return b; }

Params read_options(int argc, char** argv) {
    Params p;
    for (int i = 1; i < argc; i += 2) {
        const std::string k = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        const std::string v = argv[i + 1];
        if (k == "--T1") p.T1 = std::stod(v); else if (k == "--T3") p.T3 = std::stod(v);
        else if (k == "--gamma") { p.gamma = std::stod(v); p.gamma_text = v; }
        else if (k == "--dt") p.dt = std::stod(v); else if (k == "--snapshot") p.snapshot = std::stod(v);
        else if (k == "--nlag") p.nlag = std::stoi(v); else if (k == "--batches") p.batches = std::stoi(v);
        else if (k == "--threads") p.threads = std::stoi(v); else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--x0") p.x0_text = v;
        else if (k == "--out") p.out = v;
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    std::stringstream ss(p.x0_text); std::string t; int k = 0;
    while (std::getline(ss, t, ',')) { if (k < 5) p.x0[k] = std::stod(t); ++k; }
    if (k != 5) throw std::runtime_error("--x0 needs 5 values I1,I2,I3,theta1,theta3");
    if (p.out.empty()) { std::string tag = p.x0_text; for (auto& ch : tag) if (ch == ',') ch = '_'; p.out = "return_x" + tag + "_g" + p.gamma_text + ".json"; }
    return p;
}

// print_parameters: the parameters in the command format of the original program
//   local_test_v2 ex OUT T1 T3 gamma dt burn measure snapshot nlag batches seed threads I1 I2 I3 theta1 theta3
// (burn and measure are not used in mode ex and are written as 0).
void print_parameters(const Params& p) {
    std::cout << "command: ex " << p.out << " " << fmt(p.T1) << " " << fmt(p.T3) << " " << p.gamma_text << " " << fmt(p.dt) << " 0 0 "
              << fmt(p.snapshot) << " " << p.nlag << " " << p.batches << " " << p.seed << " " << p.threads;
    std::stringstream ss(p.x0_text); std::string t;
    while (std::getline(ss, t, ',')) std::cout << " " << t;
    std::cout << std::endl;
}

// ============================================================================ 7. simulation and output
int run(const Params& p) {
    const int NL = p.nlag, batches = p.batches;
    const auto snap_steps = (std::int64_t)std::llround(p.snapshot / p.dt);   // steps between two recorded times
    // EX[g]: sums of group g, index (time * 18 + function)
    std::vector<std::vector<double>> EX(NGROUPS, std::vector<double>((size_t)(NF + NOBSX) * NL, 0));
    std::vector<std::uint64_t> fails(batches, 0), nonfinite(batches, 0);
#pragma omp parallel for num_threads(p.threads) schedule(dynamic, 1)
    for (int b = 0; b < batches; ++b) {
        // ---- 7a. one batch: all 16 trajectories at x0, their generators, the sums of this batch
        State st; MidpointWork ws;
        set_point(st, p.x0);
        compute_force(st);
        std::array<Rng, LANES> rng;
        for (int l = 0; l < LANES; ++l) rng[l].seed(p.seed, (std::uint64_t)(b * LANES + l));
        std::vector<double> EXb((size_t)(NF + NOBSX) * NL, 0.0);
        Vec16 lx, ly, rx, ry; std::int64_t k = 0;          // k = index of the next recorded time
        // record: add the indicators and observables of the 16 trajectories at time index kk
        auto record = [&](std::int64_t kk) {
            for (int l = 0; l < LANES; ++l) {
                double o[5]; unsigned char g[NF]; reduced(st, l, o); indicators(o, g);
                double* e = &EXb[(size_t)kk * (NF + NOBSX)];
                for (int f = 0; f < NF; ++f) e[f] += g[f];
                e[NF] += o[0]; e[NF + 1] += o[1]; e[NF + 2] += o[2]; e[NF + 3] += std::cos(o[3]); e[NF + 4] += std::cos(o[4]);
                e[NF + 5] += o[0] - o[2];
            }
        };
        record(k++);                                         // t = 0
        const std::int64_t total = (NL - 1) * snap_steps;
        for (std::int64_t step = 1; step <= total; ++step) {
            // ---- 7b. one time step: noise (x of both baths, then y), midpoint step, two bath steps (alternating)
            gaussian_pair(rng, lx, rx);
            gaussian_pair(rng, ly, ry);
            int it = 0;
            if (!midpoint_step(st, ws, p.dt, it)) ++fails[b];
            if (step % 2 == 0) { bath_step(st, 0, p.T1, p.gamma, lx, ly, p.dt); bath_step(st, N - 1, p.T3, p.gamma, rx, ry, p.dt); }
            else { bath_step(st, N - 1, p.T3, p.gamma, rx, ry, p.dt); bath_step(st, 0, p.T1, p.gamma, lx, ly, p.dt); }
            if (!finite_state(st)) { nonfinite[b] = 1; break; }
            if (step % snap_steps == 0) record(k++);
        }
        // ---- 7c. add the batch to its group
#pragma omp critical
        { auto& G = EX[b % NGROUPS]; for (size_t i = 0; i < EXb.size(); ++i) G[i] += EXb[i]; }
    }
    // ---- 7d. the JSON file (17 significant digits)
    std::ofstream f(p.out);
    f << std::setprecision(17) << "{\"mode\":\"ex\",\"T1\":" << p.T1 << ",\"T3\":" << p.T3 << ",\"gamma\":" << p.gamma << ",\"dt\":" << p.dt
      << ",\"snapshot\":" << p.snapshot << ",\"nlag\":" << NL << ",\"batches\":" << batches << ",\"lanes\":" << LANES << ",\"x0\":[";
    for (int j = 0; j < 5; ++j) f << p.x0[j] << (j < 4 ? "," : "]");
    f << ",\"centers\":["; for (int c = 0; c < NC; ++c) { f << "["; for (int j = 0; j < 5; ++j) f << CEN[c][j] << (j < 4 ? "," : "]"); f << (c < NC - 1 ? "," : "]"); }
    f << ",\"radii\":[" << RAD[0] << "," << RAD[1] << "]";
    f << ",\"EX\":[";
    for (int g = 0; g < NGROUPS; ++g) {
        f << "[";
        for (size_t i = 0; i < EX[g].size(); ++i) f << EX[g][i] << (i + 1 < EX[g].size() ? "," : "");
        f << "]" << (g < NGROUPS - 1 ? "," : "");
    }
    f << "]" << ",\"ngroups\":" << NGROUPS;
    std::uint64_t F = 0, NFin = 0; for (int b = 0; b < batches; ++b) { F += fails[b]; NFin += nonfinite[b]; }
    f << ",\"midpoint_failures\":" << F << ",\"nonfinite_batches\":" << NFin << "}\n";
    std::cout << "done " << p.out << " failures=" << F << " nonfinite=" << NFin << std::endl;
    return NFin ? 2 : 0;
}

}  // namespace

// ============================================================================ 8. main
int main(int argc, char** argv) {
    try {
        const Params p = read_options(argc, argv);
        print_parameters(p);
        return p.dry_run ? 0 : run(p);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
