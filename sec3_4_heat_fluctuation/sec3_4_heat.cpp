// sec3_4_heat.cpp  --  Section 3.4: heat exchanged with the baths (fluctuation relation).
//
// WHAT IT COMPUTES
//   Many independent trajectories of the open chain c_1, ..., c_n (c_j = x_j + i y_j complex) with baths at
//   modes 1 and n (temperatures T1, Tn).  After a burn-in the measurement time is cut into consecutive blocks
//   of length tau0 (--block-time, 5 in the paper).  For every block of every trajectory it records
//     Q_left, Q_right   heat delivered by the left and the right bath to the system during the block,
//     dE                change of the energy E = H/2 during the block,
//     Sigma_medium      entropy produced in the baths, -Q_left/T1 - Q_right/Tn,
//     J                 time-averaged action current at one bond (default the bond (n/2 - 1, n/2)).
//   The analysis adds consecutive blocks to get the heat Q_tau over windows tau = 5, 10, 20, ..., 640 and fits
//   the slope of log p(Q) / p(-Q).
//
// MODEL
//   E = H/2,  H = M^2 - 1/2 sum_j I_j^2 + 2 sum_{j>=2} I_{j-1} I_j cos(theta_j),  I_j = |c_j|^2,  M = sum_j I_j.
//   F_j = dH/d conj(c_j) = dE/dx_j + i dE/dy_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2).
//   Interior modes:  dc_j = i F_j dt.   Bath modes r = 1, n (canonical bath):
//     dc_r = (i - gamma) F_r dt + sqrt(2 gamma T_r) (dW_x + i dW_y),
//   whose invariant measure in equilibrium (T1 = Tn = T) is the Gibbs density exp(-E/T) = exp(-H/(2T)).
//   The heat from bath r is the work it does on E; positive heat means energy delivered to the system.
//
// NUMERICAL SCHEME (one step of size dt)
//   1. Hamiltonian part: implicit midpoint rule  c' = c + dt i F((c + c')/2)  by fixed-point iteration
//      (tolerance 2e-13, at most 20 iterations), in the FIXED frame (the full F, including 2M c_j): the heat
//      is measured as a change of E, so the frame is not changed here.
//   2. Bath part: one Euler--Maruyama step at mode 1 and one at mode n (order alternating).  The heat of a
//      bath in this step is the exact change of E across its substep.  This converges to the Stratonovich
//      heat as dt -> 0, has no 1/I singularity, and makes Q_left + Q_right - dE equal to the energy error of
//      the midpoint steps (written in the blocks file as energy_balance_error).
//   Every mode starts at I_j = sqrt((T1 + Tn) / (2n)), all phases 0.  gamma = 0.1.
//
// PARALLEL STRUCTURE
//   Trajectories are grouped in batches of 16 (Vec16 holds one value per trajectory of a batch).  Each batch
//   is simulated from start to end by one OpenMP thread (batches are given to the threads in a fixed order).
//
// OUTPUT  (prefix OUT, default heat_T<T1>_T<Tn>_n<n>)
//   OUT_blocks.csv    one row per (trajectory, block): stream_id, block_id, q_left, q_right, delta_energy,
//                     entropy_medium, entropy_rate (= entropy_medium / tau0), action_current,
//                     energy_balance_error (= q_left + q_right - delta_energy)
//   OUT_profile.csv   time-averaged action of every mode during the measurement
//   OUT_burnin.csv    E of every trajectory at 10 equally spaced times of the burn-in (to check stationarity)
//   OUT_summary.csv   parameters and means over all blocks (heat rates, entropy rate, current, energy balance),
//                     midpoint failures and iterations, mean energy error rate of the midpoint steps, run time
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 -Ieigen sec3_4_heat.cpp -o sec3_4_heat
//   ./sec3_4_heat --T1 8 --Tn 4 --dt 1.25e-4 --seed 2026091104
// Options and their defaults (the paper settings):
//   --T1 6.5  --Tn 5.5  --n 10  --batches 8 (128 trajectories)  --burnin 500  --block-time 5
//   --blocks 31252 (blocks per trajectory)  --dt 5e-4  --seed 1  --threads 2  --bond n/2  --out (see above)
//   --dry-run 1   print the parameters and stop.
//   burnin/dt and block-time/dt must be integers.
//
// The arithmetic is exactly that of the program used for the paper (NLS_entropy_ft.cpp, command "sample"):
// with the same compiler and flags the output files are identical bit for bit (with 1 thread; with more
// threads the means can differ in the last digit because the order of the sums changes).
//
// LAYOUT OF THIS FILE
//   1. constants and types
//   2. random numbers                 Rng, gaussian_pair
//   3. chain state, forces, energy    Chain, compute_force, energy, action_current
//   4. Hamiltonian part               MidpointWork, midpoint_step
//   5. bath part                      bath_step (returns the heat)
//   6. parameters                     Params, step_count, read_options, print_parameters
//   7. simulation of one batch        BlockRecord, BatchResult, burnin_monitor_steps, run_batch
//   8. all batches and output         run
//   9. main

#include <Eigen/Dense>   // Eigen::Array: fixed-size arrays with element-wise arithmetic
#include <omp.h>         // OpenMP: one batch per thread at a time

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================================ 1. constants and types
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr int LANES = 16;                     // trajectories per batch
constexpr int MAX_ITERATIONS = 20;            // fixed-point iterations of the midpoint rule
constexpr double TOLERANCE = 2.0e-13;         // stop when max |change| < TOLERANCE
constexpr double GAMMA = 0.1;                 // bath coupling
constexpr int BURNIN_MONITORS = 10;           // number of times E is written during the burn-in
constexpr const char* MODEL_VERSION = "gibbs-cartesian-entropy-ft-v1";   // written in the summary

// Vec16: 16 doubles, one per trajectory of a batch.  Arithmetic on Vec16 is element-wise.
using Vec16 = Eigen::Array<double, LANES, 1>;
// Modes: one Vec16 per mode j = 0..n-1 (zero-based here; mode j in the code is mode j+1 in the paper).
using Modes = std::vector<Vec16, Eigen::aligned_allocator<Vec16>>;

// ============================================================================ 2. random numbers
// Rng: the Xoshiro256++ generator, one per trajectory, seeded from the base seed S and the trajectory
// number k only (so the result does not depend on the number of threads).
struct Rng {
    std::array<std::uint64_t, 4> s{};   // the 256-bit state

    // splitmix64: a simple generator used only to turn a seed into well-mixed 64-bit numbers.
    static std::uint64_t splitmix64(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    // seed trajectory k: mix S with k, then fill the four state words from splitmix64
    void seed(std::uint64_t base_seed, std::uint64_t trajectory) {
        std::uint64_t x = base_seed ^ (0xD2B74407B1CE6E93ULL * (trajectory + 1ULL));
        std::uint64_t value = splitmix64(x);
        for (auto& v : s) v = splitmix64(value);
    }
    // next 64-bit random number (the Xoshiro256++ update)
    std::uint64_t next() {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        return result;
    }
    // uniform number in the open interval (0, 1) (0 never occurs, so log(u) is finite)
    double uniform() { return ((next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
};

// gaussian_pair: two independent N(0,1) numbers for each of the 16 trajectories (Box--Muller:
// R = sqrt(-2 log u1), angle 2 pi u2, the numbers R cos(angle) and R sin(angle)).
void gaussian_pair(std::array<Rng, LANES>& rng, Vec16& first, Vec16& second) {
    for (int lane = 0; lane < LANES; ++lane) {
        const double u1 = rng[lane].uniform(), u2 = rng[lane].uniform();
        const double radius = std::sqrt(-2.0 * std::log(u1)), angle = 2.0 * PI * u2;
        first(lane) = radius * std::cos(angle);
        second(lane) = radius * std::sin(angle);
    }
}

// ============================================================================ 3. chain state, forces, energy
// Chain: the state of one batch.  x[j], y[j] = real and imaginary part of c_j; compute_force fills
// I[j] = |c_j|^2, (sq_re, sq_im) = c_j^2, (F_re, F_im) = F_j and M = sum_j I_j.
struct Chain {
    int n;
    Modes x, y, I, sq_re, sq_im, F_re, F_im;
    Vec16 M = Vec16::Zero();
    explicit Chain(int n_) : n(n_), x(n_), y(n_), I(n_), sq_re(n_), sq_im(n_), F_re(n_), F_im(n_) {}
};

// compute_force: F_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2)  (missing neighbours count 0).
// With s = c_{j-1}^2 + c_{j+1}^2 = s_re + i s_im:  conj(c) s = (x s_re + y s_im) + i (x s_im - y s_re).
void compute_force(Chain& s) {
    s.M.setZero();
    for (int j = 0; j < s.n; ++j) {
        s.I[j] = s.x[j].square() + s.y[j].square();
        s.sq_re[j] = s.x[j].square() - s.y[j].square();
        s.sq_im[j] = 2.0 * s.x[j] * s.y[j];
        s.M += s.I[j];
    }
    for (int j = 0; j < s.n; ++j) {
        Vec16 nb_re = Vec16::Zero(), nb_im = Vec16::Zero();          // s = c_{j-1}^2 + c_{j+1}^2
        if (j > 0) { nb_re += s.sq_re[j - 1]; nb_im += s.sq_im[j - 1]; }
        if (j + 1 < s.n) { nb_re += s.sq_re[j + 1]; nb_im += s.sq_im[j + 1]; }
        const Vec16 onsite = Vec16(2.0 * s.M - s.I[j]);
        s.F_re[j] = onsite * s.x[j] + 2.0 * (nb_re * s.x[j] + nb_im * s.y[j]);
        s.F_im[j] = onsite * s.y[j] + 2.0 * (nb_im * s.x[j] - nb_re * s.y[j]);
    }
}

// energy: E = H/2 = M^2/2 - 1/4 sum_j I_j^2 + sum_{j>=2} Re(c_j^2 conj(c_{j-1})^2)
// (uses the arrays filled by compute_force; Re(c_j^2 conj(c_{j-1}^2)) = I_{j-1} I_j cos theta_j).
Vec16 energy(const Chain& s) {
    Vec16 sum_sq = Vec16::Zero(), coupling = Vec16::Zero();
    for (int j = 0; j < s.n; ++j) {
        sum_sq += s.I[j].square();
        if (j > 0) coupling += s.sq_re[j] * s.sq_re[j - 1] + s.sq_im[j] * s.sq_im[j - 1];
    }
    return 0.5 * s.M.square() - 0.25 * sum_sq + coupling;
}

// action_current: J at the bond (bond-1, bond) (zero-based modes),
//   J = 4 I_{j-1} I_j sin theta_j = 4 Im(c_j^2 conj(c_{j-1})^2).
// (Declared inline as in the original program: with gcc this keeps the rounding of a*b - c*d identical.)
inline Vec16 action_current(const Chain& s, int bond) {
    return 4.0 * (s.sq_im[bond] * s.sq_re[bond - 1] - s.sq_re[bond] * s.sq_im[bond - 1]);
}

// ============================================================================ 4. Hamiltonian part
// MidpointWork: scratch space for one batch.
struct MidpointWork {
    Chain mid;        // the midpoint (c + c')/2, with its forces
    Modes cx, cy;     // current guess for c'
    explicit MidpointWork(int n) : mid(n), cx(n), cy(n) {}
};

// midpoint_step: implicit midpoint step c' = c + dt i F((c + c')/2) with the full F (fixed frame).
// In real form x' = x - dt F_im(mid), y' = y + dt F_re(mid).  Start from the explicit Euler step and iterate
// until no coordinate of any of the 16 trajectories changes by TOLERANCE or more.  At the end s holds c' and
// its forces.  Returns false if the iteration did not converge in MAX_ITERATIONS iterations; `iterations`
// returns the number of iterations used (at most MAX_ITERATIONS).
bool midpoint_step(Chain& s, MidpointWork& w, double dt, int& iterations) {
    compute_force(s);
    for (int j = 0; j < s.n; ++j) {                   // first guess: explicit Euler step
        w.cx[j] = s.x[j] - dt * s.F_im[j];
        w.cy[j] = s.y[j] + dt * s.F_re[j];
    }
    bool converged = false;
    int it;
    for (it = 1; it <= MAX_ITERATIONS; ++it) {
        for (int j = 0; j < s.n; ++j) {               // midpoint and its forces
            w.mid.x[j] = 0.5 * (s.x[j] + w.cx[j]);
            w.mid.y[j] = 0.5 * (s.y[j] + w.cy[j]);
        }
        compute_force(w.mid);
        double change = 0.0;
        for (int j = 0; j < s.n; ++j) {               // new guess c' = c + dt i F(midpoint)
            const Vec16 new_x = s.x[j] - dt * w.mid.F_im[j];
            const Vec16 new_y = s.y[j] + dt * w.mid.F_re[j];
            change = std::max(change, (new_x - w.cx[j]).abs().maxCoeff());
            change = std::max(change, (new_y - w.cy[j]).abs().maxCoeff());
            w.cx[j] = new_x; w.cy[j] = new_y;
        }
        if (change < TOLERANCE) { converged = true; break; }
    }
    iterations = std::min(it, MAX_ITERATIONS);
    for (int j = 0; j < s.n; ++j) { s.x[j] = w.cx[j]; s.y[j] = w.cy[j]; }
    compute_force(s);
    return converged;
}

// ============================================================================ 5. bath part
// bath_step: Euler--Maruyama step of the canonical bath at mode r,
//     c_r += -gamma F_r dt + sqrt(2 gamma T) sqrt(dt) (N_x + i N_y),
// using the forces of the current state; recomputes the forces and returns the heat of this substep,
// i.e. the change of E (one value per trajectory).
Vec16 bath_step(Chain& s, int r, double T, const Vec16& Nx, const Vec16& Ny, double dt, double sqrt_dt) {
    const Vec16 E_before = energy(s);
    const Vec16 dx = -GAMMA * s.F_re[r], dy = -GAMMA * s.F_im[r];
    const double noise = std::sqrt(2.0 * GAMMA * T);
    s.x[r] += dx * dt + noise * sqrt_dt * Nx;
    s.y[r] += dy * dt + noise * sqrt_dt * Ny;
    compute_force(s);
    return energy(s) - E_before;
}

// ============================================================================ 6. parameters
// Params: the command-line options.  The initial values are the paper settings.
struct Params {
    double T1 = 6.5, Tn = 5.5;                    // bath temperatures at modes 1 and n
    int n = 10;                                   // number of modes
    int batches = 8;                              // batches of 16 trajectories
    double burnin = 500;                          // burn-in time
    double block_time = 5;                        // length tau0 of one block
    int blocks = 31252;                           // blocks per trajectory
    double dt = 5e-4;                             // time step
    std::uint64_t seed = 1;                       // base seed
    int threads = 2;                              // OpenMP threads
    int bond = -1;                                // bond of the current (right mode, zero-based); < 0: n/2
    bool bond_given = false;                      // --bond was given (printed in the parameter line only then)
    std::string out;                              // output prefix; empty: heat_T<T1>_T<Tn>_n<n>
    bool dry_run = false;
    std::int64_t burn_steps = 0, block_steps = 0; // burnin/dt and block_time/dt
};

std::string fmt(double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return b; }

// step_count: time/dt as an integer; an error if time/dt is not an integer (up to rounding).
std::int64_t step_count(double time, double dt, const std::string& label, bool allow_zero) {
    const double raw = time / dt;
    const auto steps = static_cast<std::int64_t>(std::llround(raw));
    if ((!allow_zero && steps <= 0) || (allow_zero && steps < 0)) throw std::invalid_argument(label + " has invalid step count");
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(time));
    if (std::abs(steps * dt - time) > tolerance) {
        std::ostringstream message;
        message << label << "/dt must be an integer; got " << std::setprecision(17) << raw;
        throw std::invalid_argument(message.str());
    }
    return steps;
}

// read_options: read "--name value" pairs, fill in the defaults that depend on other options, check values.
Params read_options(int argc, char** argv) {
    Params p;
    for (int i = 1; i < argc; i += 2) {
        const std::string k = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        const std::string v = argv[i + 1];
        if (k == "--T1") p.T1 = std::stod(v);
        else if (k == "--Tn") p.Tn = std::stod(v);
        else if (k == "--n") p.n = std::stoi(v);
        else if (k == "--batches") p.batches = std::stoi(v);
        else if (k == "--burnin") p.burnin = std::stod(v);
        else if (k == "--block-time") p.block_time = std::stod(v);
        else if (k == "--blocks") p.blocks = std::stoi(v);
        else if (k == "--dt") p.dt = std::stod(v);
        else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--threads") p.threads = std::stoi(v);
        else if (k == "--bond") { p.bond = std::stoi(v); p.bond_given = true; }
        else if (k == "--out") p.out = v;
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    if (p.bond < 0) p.bond = p.n / 2;
    if (p.out.empty()) p.out = "heat_T" + fmt(p.T1) + "_T" + fmt(p.Tn) + "_n" + std::to_string(p.n);
    if (!(p.T1 > 0.0) || !(p.Tn > 0.0)) throw std::invalid_argument("T1 and Tn must be positive");
    if (p.n < 2) throw std::invalid_argument("n must be at least 2");
    if (p.batches <= 0 || p.blocks <= 0 || p.threads <= 0) throw std::invalid_argument("batches, blocks and threads must be positive");
    if (!(p.dt > 0.0) || !(p.block_time > 0.0) || p.burnin < 0.0) throw std::invalid_argument("dt and block-time must be positive and burnin nonnegative");
    if (p.bond < 1 || p.bond >= p.n) throw std::invalid_argument("bond must satisfy 1 <= bond < n");
    p.burn_steps = step_count(p.burnin, p.dt, "burnin", true);
    p.block_steps = step_count(p.block_time, p.dt, "block_time", false);
    return p;
}

// print_parameters: the parameters in the command format of the original program
//   NLS_entropy_ft sample T1 Tn n batches burnin block_time blocks dt seed threads prefix [bond]
// (the check scripts give this line to the original program and compare the outputs; as in the original,
// the bond is written only if it was given, otherwise the original takes n/2 as well).
void print_parameters(const Params& p) {
    std::cout << "command: sample " << fmt(p.T1) << " " << fmt(p.Tn) << " " << p.n << " " << p.batches << " " << fmt(p.burnin) << " "
              << fmt(p.block_time) << " " << p.blocks << " " << fmt(p.dt) << " " << p.seed << " " << p.threads << " " << p.out;
    if (p.bond_given) std::cout << " " << p.bond;
    std::cout << std::endl;
}

// ============================================================================ 7. simulation of one batch
// BlockRecord: what is recorded for one block of one trajectory.
struct BlockRecord {
    double q_left = 0.0, q_right = 0.0;           // heat from the two baths during the block
    double delta_energy = 0.0;                    // E(end of block) - E(start of block)
    double entropy_medium = 0.0;                  // -q_left/T1 - q_right/Tn
    double action_current = 0.0;                  // time average of J at the chosen bond
};

// BatchResult: everything one batch returns.
struct BatchResult {
    std::vector<BlockRecord> blocks;              // index block * 16 + lane
    std::vector<double> burnin_energy;            // index monitor * 16 + lane
    Modes action_integral;                        // integral of I_j over the measurement, per mode
    std::uint64_t failures = 0, iteration_sum = 0;
    Vec16 hamiltonian_energy_error = Vec16::Zero();   // sum of E(after) - E(before) over the midpoint steps
};

// burnin_monitor_steps: the steps after which E is written during the burn-in, at 1/10, 2/10, ..., 1 of it.
std::vector<std::int64_t> burnin_monitor_steps(std::int64_t burn_steps) {
    std::vector<std::int64_t> result;
    if (burn_steps <= 0) return result;
    for (int k = 1; k <= BURNIN_MONITORS; ++k) {
        const auto step = std::max<std::int64_t>(1, (burn_steps * k) / BURNIN_MONITORS);
        if (result.empty() || result.back() != step) result.push_back(step);
    }
    return result;
}

// run_batch: simulate the 16 trajectories of batch b from the initial state to the end of the last block.
BatchResult run_batch(const Params& p, int b) {
    // ---- 7a. set-up: records, initial state c_j = I0^(1/2), I0 = sqrt((T1+Tn)/(2n)), generators
    BatchResult result;
    result.blocks.assign(static_cast<std::size_t>(p.blocks) * LANES, {});
    const auto monitors = burnin_monitor_steps(p.burn_steps);
    result.burnin_energy.assign(monitors.size() * LANES, 0.0);
    result.action_integral.resize(p.n);
    for (auto& v : result.action_integral) v.setZero();

    Chain s(p.n);
    MidpointWork w(p.n);
    const double amplitude = std::sqrt(std::sqrt(0.5 * (p.T1 + p.Tn) / static_cast<double>(p.n)));
    for (int j = 0; j < p.n; ++j) { s.x[j].setConstant(amplitude); s.y[j].setZero(); }
    std::array<Rng, LANES> rng;
    for (int l = 0; l < LANES; ++l) rng[l].seed(p.seed, static_cast<std::uint64_t>(b) * LANES + l);

    const double sqrt_dt = std::sqrt(p.dt);
    const auto total_steps = p.burn_steps + p.block_steps * static_cast<std::int64_t>(p.blocks);

    // running values of the current block
    Vec16 q_left = Vec16::Zero(), q_right = Vec16::Zero(), current_integral = Vec16::Zero(), block_start_energy = Vec16::Zero();
    std::size_t monitor = 0;
    int block = 0;

    for (std::int64_t step = 0; step < total_steps; ++step) {
        // ---- 7b. values at the start of the step: E (start of the first block) and J at the chosen bond
        compute_force(s);
        const Vec16 E_now = energy(s);
        const Vec16 J_now = action_current(s, p.bond);
        if (step == p.burn_steps) block_start_energy = E_now;

        // ---- 7c. noise for this step: (left x, right x), then (left y, right y)
        Vec16 N1x, Nnx, N1y, Nny;
        gaussian_pair(rng, N1x, Nnx);
        gaussian_pair(rng, N1y, Nny);

        // ---- 7d. during the measurement: integrals of J and of every I_j (left-point rule)
        if (step >= p.burn_steps) {
            current_integral += J_now * p.dt;
            for (int j = 0; j < p.n; ++j) result.action_integral[j] += s.I[j] * p.dt;
        }

        // ---- 7e. Hamiltonian step, with its energy error (E after - E before; zero for the exact flow)
        const Vec16 E_before = energy(s);
        int it = 0;
        if (!midpoint_step(s, w, p.dt, it)) ++result.failures;
        result.iteration_sum += static_cast<std::uint64_t>(it);
        result.hamiltonian_energy_error += energy(s) - E_before;

        // ---- 7f. the two bath steps (left first on even steps); add their heat during the measurement
        Vec16 dq_left = Vec16::Zero(), dq_right = Vec16::Zero();
        if (step % 2 == 0) {
            dq_left = bath_step(s, 0, p.T1, N1x, N1y, p.dt, sqrt_dt);
            dq_right = bath_step(s, p.n - 1, p.Tn, Nnx, Nny, p.dt, sqrt_dt);
        } else {
            dq_right = bath_step(s, p.n - 1, p.Tn, Nnx, Nny, p.dt, sqrt_dt);
            dq_left = bath_step(s, 0, p.T1, N1x, N1y, p.dt, sqrt_dt);
        }
        if (step >= p.burn_steps) { q_left += dq_left; q_right += dq_right; }

        // ---- 7g. burn-in monitor: store E at the monitor steps
        if (monitor < monitors.size() && step + 1 == monitors[monitor]) {
            compute_force(s);
            const Vec16 E_monitor = energy(s);
            for (int l = 0; l < LANES; ++l) result.burnin_energy[monitor * LANES + l] = E_monitor(l);
            ++monitor;
        }

        // ---- 7h. end of a block: store the block record of every trajectory and start a new block
        if (step >= p.burn_steps && (step - p.burn_steps + 1) % p.block_steps == 0) {
            compute_force(s);
            const Vec16 block_end_energy = energy(s);
            for (int l = 0; l < LANES; ++l) {
                BlockRecord& rec = result.blocks[static_cast<std::size_t>(block) * LANES + l];
                rec.q_left = q_left(l);
                rec.q_right = q_right(l);
                rec.delta_energy = block_end_energy(l) - block_start_energy(l);
                rec.entropy_medium = -q_left(l) / p.T1 - q_right(l) / p.Tn;
                rec.action_current = current_integral(l) / p.block_time;
            }
            q_left.setZero(); q_right.setZero(); current_integral.setZero();
            block_start_energy = block_end_energy;
            ++block;
        }
    }
    return result;
}

// ============================================================================ 8. all batches and output
double mean_of(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size()); }
double rms_of(const std::vector<double>& v) { double t = 0.0; for (double x : v) t += x * x; return std::sqrt(t / static_cast<double>(v.size())); }

// open_output: a text file written with 17 significant digits (exact round trip of doubles)
std::ofstream open_output(const std::string& path) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot open output file: " + path);
    stream << std::setprecision(17);
    return stream;
}

void run(const Params& p) {
    const std::filesystem::path prefix_path(p.out);
    if (prefix_path.has_parent_path()) std::filesystem::create_directories(prefix_path.parent_path());
    const int L = p.batches * LANES;                                  // number of trajectories
    const std::size_t total_blocks = static_cast<std::size_t>(L) * p.blocks;
    std::cout << "sec3_4_heat n=" << p.n << " T1=" << p.T1 << " Tn=" << p.Tn << " dt=" << p.dt << " burnin=" << p.burnin
              << " block_time=" << p.block_time << " blocks_per_stream=" << p.blocks << " streams=" << L << " bond=" << p.bond
              << " seed=" << p.seed << " threads=" << p.threads << "\n";

    // ---- 8a. run the batches; copy each batch's records into stream order (stream = 16 b + lane)
    const auto start = std::chrono::steady_clock::now();
    const auto monitors = burnin_monitor_steps(p.burn_steps);
    std::vector<BlockRecord> all_blocks(total_blocks);
    std::vector<double> all_burnin(static_cast<std::size_t>(L) * monitors.size(), 0.0);
    Modes action_integral(p.n);
    for (auto& v : action_integral) v.setZero();
    std::uint64_t failures = 0, iteration_sum = 0;
    Vec16 hamiltonian_energy_error = Vec16::Zero();
    std::atomic<int> completed{0};

#pragma omp parallel for schedule(static) num_threads(p.threads)
    for (int b = 0; b < p.batches; ++b) {
        auto r = run_batch(p, b);
        for (int k = 0; k < p.blocks; ++k)
            for (int l = 0; l < LANES; ++l)
                all_blocks[static_cast<std::size_t>(b * LANES + l) * p.blocks + k] = r.blocks[static_cast<std::size_t>(k) * LANES + l];
        for (std::size_t m = 0; m < monitors.size(); ++m)
            for (int l = 0; l < LANES; ++l)
                all_burnin[static_cast<std::size_t>(b * LANES + l) * monitors.size() + m] = r.burnin_energy[m * LANES + l];
#pragma omp critical
        {
            for (int j = 0; j < p.n; ++j) action_integral[j] += r.action_integral[j];
            failures += r.failures;
            iteration_sum += r.iteration_sum;
            hamiltonian_energy_error += r.hamiltonian_energy_error;
        }
        const int done = ++completed;
        if (done == p.batches || done % std::max(1, p.batches / 10) == 0) {
#pragma omp critical
            std::cout << "completed " << done << "/" << p.batches << " batches\n";
        }
    }
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    // ---- 8b. means over all blocks
    const double midpoint_calls = static_cast<double>(p.batches) * (p.burn_steps + p.block_steps * static_cast<std::int64_t>(p.blocks));
    const double failure_rate = static_cast<double>(failures) / midpoint_calls;
    const double mean_iterations = static_cast<double>(iteration_sum) / midpoint_calls;
    const double energy_error_rate = hamiltonian_energy_error.sum() / (static_cast<double>(L) * (p.burnin + p.block_time * p.blocks));
    std::vector<double> q_left_rate, q_right_rate, energy_drift_rate, entropy_rate, current, balance_rate;
    for (const auto& r : all_blocks) {
        q_left_rate.push_back(r.q_left / p.block_time);
        q_right_rate.push_back(r.q_right / p.block_time);
        energy_drift_rate.push_back(r.delta_energy / p.block_time);
        entropy_rate.push_back(r.entropy_medium / p.block_time);
        current.push_back(r.action_current);
        balance_rate.push_back((r.q_left + r.q_right - r.delta_energy) / p.block_time);
    }

    // ---- 8c. OUT_blocks.csv
    {
        auto o = open_output(p.out + "_blocks.csv");
        o << "stream_id,block_id,q_left,q_right,delta_energy,entropy_medium,entropy_rate,action_current,energy_balance_error\n";
        for (int stream = 0; stream < L; ++stream)
            for (int k = 0; k < p.blocks; ++k) {
                const auto& r = all_blocks[static_cast<std::size_t>(stream) * p.blocks + k];
                o << stream << ',' << k << ',' << r.q_left << ',' << r.q_right << ',' << r.delta_energy << ',' << r.entropy_medium << ','
                  << r.entropy_medium / p.block_time << ',' << r.action_current << ',' << r.q_left + r.q_right - r.delta_energy << '\n';
            }
    }
    // ---- 8d. OUT_profile.csv: time and trajectory average of I_j during the measurement
    {
        auto o = open_output(p.out + "_profile.csv");
        o << "site,mean_action\n";
        const double denominator = static_cast<double>(L) * p.blocks * p.block_time;
        for (int j = 0; j < p.n; ++j) o << j + 1 << ',' << action_integral[j].sum() / denominator << '\n';
    }
    // ---- 8e. OUT_burnin.csv: E during the burn-in
    {
        auto o = open_output(p.out + "_burnin.csv");
        o << "stream_id,time,energy\n";
        for (int stream = 0; stream < L; ++stream)
            for (std::size_t m = 0; m < monitors.size(); ++m)
                o << stream << ',' << monitors[m] * p.dt << ',' << all_burnin[static_cast<std::size_t>(stream) * monitors.size() + m] << '\n';
    }
    // ---- 8f. OUT_summary.csv
    {
        auto o = open_output(p.out + "_summary.csv");
        o << "model_version,coordinates,energy_convention,n,T1,Tn,gamma,dt,burnin,block_time,blocks_per_stream,batches,lanes,n_streams,"
          << "n_blocks,bond,seed,threads,mean_q_left_rate,mean_q_right_rate,mean_energy_drift_rate,mean_entropy_rate,mean_action_current,"
          << "mean_energy_balance_error_rate,rms_energy_balance_error_rate,midpoint_failure_count,midpoint_failure_rate,"
          << "mean_midpoint_iterations,mean_hamiltonian_energy_error_rate,elapsed_seconds\n";
        o << MODEL_VERSION << ",cartesian,E=H/2," << p.n << ',' << p.T1 << ',' << p.Tn << ',' << GAMMA << ',' << p.dt << ',' << p.burnin << ','
          << p.block_time << ',' << p.blocks << ',' << p.batches << ',' << LANES << ',' << L << ',' << total_blocks << ',' << p.bond << ','
          << p.seed << ',' << p.threads << ',' << mean_of(q_left_rate) << ',' << mean_of(q_right_rate) << ',' << mean_of(energy_drift_rate) << ','
          << mean_of(entropy_rate) << ',' << mean_of(current) << ',' << mean_of(balance_rate) << ',' << rms_of(balance_rate) << ','
          << failures << ',' << failure_rate << ',' << mean_iterations << ',' << energy_error_rate << ',' << elapsed << '\n';
    }
    std::cout << "mean q_left rate=" << mean_of(q_left_rate) << " mean q_right rate=" << mean_of(q_right_rate)
              << " mean entropy rate=" << mean_of(entropy_rate) << " mean action current=" << mean_of(current)
              << " balance RMS rate=" << rms_of(balance_rate) << " midpoint failures=" << failures
              << " mean midpoint iterations=" << mean_iterations << " elapsed=" << elapsed << " s\n";
}

// ============================================================================ 9. main
int main(int argc, char** argv) {
    try {
        const Params p = read_options(argc, argv);
        print_parameters(p);
        if (!p.dry_run) run(p);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
