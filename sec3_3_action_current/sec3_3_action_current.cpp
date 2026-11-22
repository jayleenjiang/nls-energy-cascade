// sec3_3_action_current.cpp  
// WHAT IT COMPUTES
//   Many independent trajectories of the chain c_1, ..., c_n (c_j = x_j + i y_j complex), with baths at
//   modes 1 and n at temperatures T1 and Tn.  After a burn-in it records, every 0.1 time units,
//     the actions          I_j = |c_j|^2,                                  j = 1..n
//     the action currents  J_j = 4 I_{j-1} I_j sin(theta_j) = 4 Im( c_j^2 conj(c_{j-1})^2 ),   j = 2..n
//   and writes their time averages (averaged over trajectories, with standard errors).
//
// MODEL
//   H = M^2 - 1/2 sum_j I_j^2 + 2 sum_{j>=2} I_{j-1} I_j cos(theta_j),   M = I_1 + ... + I_n.
//   Interior modes:  dc_j/dt = i F_j,   F_j = dH/d conj(c_j) = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2).
//   Bath modes r = 1, n:   dc_r = i F_r dt - gamma D_r dt + sqrt(2 gamma T_r) (dW_x + i dW_y),
//   with the drift D_r chosen by --bc:
//     1  : D = F                                  (canonical bath)
//     2  : D = F_rad = Re(conj(c) F)/|c|^2 * c    (radial part of F only)
//     3  : D = |c|^2 c                            (local action drift)
//     3b : D = |c|^2 c + F_tan,  F_tan = Im(conj(c) F)/|c|^2 * i c
//
// NUMERICAL SCHEME (one step of size dt)
//   1. Hamiltonian part: implicit midpoint rule  c' = c + dt i F((c + c')/2),  solved by fixed-point
//      iteration (stop when the change is < 2e-13, at most 20 iterations).  It is done in the co-rotating
//      frame: the term 2M c_j of F_j only rotates all phases together and does not change I_j or theta_j,
//      so it is left out here.  M is still conserved exactly by the midpoint rule.
//   2. Bath part: one Euler--Maruyama step at mode 1 and one at mode n, using the full F.  The order of the
//      two baths alternates from step to step.
//   Every mode starts at I_j = sqrt((T1 + Tn) / (2n)), all phases 0.
//
// OUTPUT  (prefix OUT, default BC<bc>_T<T1>_T<Tn>_n<n>)
//   OUT_summary.txt          parameters; mean_I (j = 1..n), se_I, mean_current (j = 2..n), se_current
//                            (standard errors over the independent trajectories); the profile <I_j> at
//                            4 checkpoints during the run (to see that it has stopped changing).
//   OUT_current_streams.txt  one row per trajectory: time-averaged M, then time-averaged J_2 .. J_n.
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 -Ieigen sec3_3_action_current.cpp -o sec3_3_action_current
//   ./sec3_3_action_current --n 100 --bc 3b --seed 202611090123
//   ./sec3_3_action_current --n 15 --T1 4.875 --Tn 5.125 --seed 202611090000     (linear response, dT = 0.25)
//
// LAYOUT OF THIS FILE
//   1. constants and types
//   2. random numbers                 Rng, gaussian_pair
//   3. chain state and forces         Chain, compute_force
//   4. Hamiltonian part               MidpointWork, midpoint_step
//   5. bath part                      Coupling, bath_step
//   6. parameters                     Params, read_options, print_parameters
//   7. simulation                     run: set-up, one time step, recording, burn-in, measurement, output
//   8. main

#include <Eigen/Dense>   // Eigen::Array: fixed-size arrays with element-wise arithmetic
#include <omp.h>         // OpenMP: one batch per thread at a time

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================================ 1. constants and types
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr int LANES = 16;                     // trajectories per batch
constexpr int MAX_ITERATIONS = 20;            // fixed-point iterations of the midpoint rule
constexpr double TOLERANCE = 2.0e-13;         // stop when max |change| < TOLERANCE
constexpr double SAMPLE_INTERVAL = 0.1;       // time between two recorded samples
constexpr int CHECKPOINTS = 4;                // profile printed after each quarter of the measurement

// Vec16: 16 doubles, one per trajectory of a batch.  Arithmetic on Vec16 is element-wise, e.g.
// a * b multiplies the 16 pairs (a(l), b(l)), and a.square() squares each entry.
using Vec16 = Eigen::Array<double, LANES, 1>;
// Modes: one Vec16 per mode j = 0..n-1 (zero-based here; mode j in the code is mode j+1 in the paper).
// The aligned allocator is required by Eigen for vectors of fixed-size arrays.
using Modes = std::vector<Vec16, Eigen::aligned_allocator<Vec16>>;

// ============================================================================ 2. random numbers
// Rng: the Xoshiro256++ generator (a standard 64-bit generator with a 256-bit state), one per trajectory.
// Each trajectory gets its own generator, seeded from the base seed S and the trajectory number k only.
// So trajectory k is the same whatever the number of threads or the order in which batches are run.
struct Rng {
    std::array<std::uint64_t, 4> s{};   // the 256-bit state

    // splitmix64: a simple generator used only to turn a seed into well-mixed 64-bit numbers.
    static std::uint64_t splitmix64(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    // rotate the bits of x left by k places
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

    // Seed trajectory k: mix S with k (multiplication by a large odd constant), then fill the four
    // state words from splitmix64.
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
    // uniform number in the open interval (0, 1): the top 53 bits, shifted by half a unit so that 0 never
    // occurs (log(0) in Box--Muller would be infinite)
    double uniform() { return ((next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
};

// gaussian_pair: two independent N(0,1) numbers for each of the 16 trajectories, by the Box--Muller
// transform: with u1, u2 uniform on (0,1), R = sqrt(-2 log u1) and angle 2 pi u2, the numbers
// R cos(angle) and R sin(angle) are independent standard normals.  Used for the noise in x and in y.
void gaussian_pair(std::array<Rng, LANES>& rng, Vec16& first, Vec16& second) {
    for (int lane = 0; lane < LANES; ++lane) {
        const double u1 = rng[lane].uniform(), u2 = rng[lane].uniform();
        const double radius = std::sqrt(-2.0 * std::log(u1)), angle = 2.0 * PI * u2;
        first(lane) = radius * std::cos(angle);
        second(lane) = radius * std::sin(angle);
    }
}

// ============================================================================ 3. chain state and forces
// Chain: the state of one batch (16 trajectories of the whole chain).
//   x[j], y[j]          real and imaginary part of c_j  (these are the variables that are integrated)
// and quantities computed from them by compute_force, kept because several parts of the step use them:
//   I[j]                I_j = |c_j|^2 = x^2 + y^2
//   sq_re[j], sq_im[j]  real and imaginary part of c_j^2 = (x^2 - y^2) + i (2 x y)
//   F_re[j], F_im[j]    real and imaginary part of the force F_j
//   M                   M = I_1 + ... + I_n
struct Chain {
    int n;
    Modes x, y, I, sq_re, sq_im, F_re, F_im;
    Vec16 M = Vec16::Zero();
    explicit Chain(int n_) : n(n_), x(n_), y(n_), I(n_), sq_re(n_), sq_im(n_), F_re(n_), F_im(n_) {}
};

// compute_force: fill I, c^2, M and the force from x, y.
//   F_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2),
// where a neighbour that does not exist (j-1 = 0 or j+1 = n+1) contributes 0.
// Writing s = c_{j-1}^2 + c_{j+1}^2 = s_re + i s_im and c = x + i y,
//   conj(c) s = (x s_re + y s_im) + i (x s_im - y s_re),
// which gives the two lines for F_re and F_im below.
// include_rotation = false leaves out the term 2M c_j (used in the co-rotating Hamiltonian step).
void compute_force(Chain& s, bool include_rotation = true) {
    // pass 1: I_j, c_j^2 and M (needed by every F_j)
    s.M.setZero();
    for (int j = 0; j < s.n; ++j) {
        s.I[j] = s.x[j].square() + s.y[j].square();
        s.sq_re[j] = s.x[j].square() - s.y[j].square();
        s.sq_im[j] = 2.0 * s.x[j] * s.y[j];
        s.M += s.I[j];
    }
    // pass 2: F_j
    for (int j = 0; j < s.n; ++j) {
        Vec16 nb_re = Vec16::Zero(), nb_im = Vec16::Zero();          // s = c_{j-1}^2 + c_{j+1}^2
        if (j > 0) { nb_re += s.sq_re[j - 1]; nb_im += s.sq_im[j - 1]; }
        if (j + 1 < s.n) { nb_re += s.sq_re[j + 1]; nb_im += s.sq_im[j + 1]; }
        // on-site factor multiplying c_j: (2M - I_j), or -I_j without the rotation term
        const Vec16 onsite = include_rotation ? Vec16(2.0 * s.M - s.I[j]) : Vec16(-s.I[j]);
        s.F_re[j] = onsite * s.x[j] + 2.0 * (nb_re * s.x[j] + nb_im * s.y[j]);
        s.F_im[j] = onsite * s.y[j] + 2.0 * (nb_im * s.x[j] - nb_re * s.y[j]);
    }
}

// ============================================================================ 4. Hamiltonian part
// The Hamiltonian flow dc/dt = i F(c) is advanced by the implicit midpoint rule
//     c' = c + dt i F((c + c')/2).
// This is symplectic and time-reversible and conserves M exactly (M is a quadratic invariant).
// In real form, i F = -F_im + i F_re, so  x' = x - dt F_im(mid),  y' = y + dt F_re(mid).
// The equation for c' is solved by fixed-point iteration: start from the explicit Euler step
// c' = c + dt i F(c), then repeat  c' <- c + dt i F((c + c')/2)  until c' changes by less than TOLERANCE.

// MidpointWork: scratch space for one thread (allocated once per thread, reused every step).
struct MidpointWork {
    Chain mid;        // the midpoint (c + c')/2, with its forces
    Modes cx, cy;     // current guess for c' = cx + i cy
    explicit MidpointWork(int n) : mid(n), cx(n), cy(n) {}
};

// midpoint_step: one implicit midpoint step of size dt for the batch s, in the co-rotating frame (forces
// without 2M c_j).  At the end s holds c' and the FULL forces of c' (the bath step needs F with the
// 2M c_j term).  Returns false if the iteration did not converge in MAX_ITERATIONS iterations (the last
// iterate is kept; this is counted as midpoint_failures in the summary and was 0 in all paper runs).
// `iterations` returns the number of iterations used (MAX_ITERATIONS + 1 if it did not converge).
bool midpoint_step(Chain& s, MidpointWork& w, double dt, int& iterations) {
    // first guess: explicit Euler step with the forces of the current state
    compute_force(s, false);
    for (int j = 0; j < s.n; ++j) {
        w.cx[j] = s.x[j] - dt * s.F_im[j];
        w.cy[j] = s.y[j] + dt * s.F_re[j];
    }
    bool converged = false;
    for (iterations = 1; iterations <= MAX_ITERATIONS; ++iterations) {
        // midpoint (c + c')/2 and its forces
        for (int j = 0; j < s.n; ++j) {
            w.mid.x[j] = 0.5 * (s.x[j] + w.cx[j]);
            w.mid.y[j] = 0.5 * (s.y[j] + w.cy[j]);
        }
        compute_force(w.mid, false);
        // new guess c' = c + dt i F(midpoint); `change` is the largest change of any coordinate in any
        // of the 16 trajectories (all 16 must converge)
        double change = 0.0;
        for (int j = 0; j < s.n; ++j) {
            const Vec16 new_x = s.x[j] - dt * w.mid.F_im[j];
            const Vec16 new_y = s.y[j] + dt * w.mid.F_re[j];
            change = std::max(change, (new_x - w.cx[j]).abs().maxCoeff());
            change = std::max(change, (new_y - w.cy[j]).abs().maxCoeff());
            w.cx[j] = new_x; w.cy[j] = new_y;
        }
        if (change < TOLERANCE) { converged = true; break; }
    }
    // accept c' and compute its full forces (with the 2M c_j term) for the bath step
    for (int j = 0; j < s.n; ++j) { s.x[j] = w.cx[j]; s.y[j] = w.cy[j]; }
    compute_force(s);
    return converged;
}

// ============================================================================ 5. bath part
// The four boundary couplings of Sec. 3.1 (option --bc 1, 2, 3, 3b).
enum class Coupling { BC1, BC2, BC3, BC3b };

// parse_coupling: turn the text of --bc into a Coupling (error for anything else).
Coupling parse_coupling(const std::string& s) {
    if (s == "1") return Coupling::BC1;
    if (s == "2") return Coupling::BC2;
    if (s == "3") return Coupling::BC3;
    if (s == "3b") return Coupling::BC3b;
    throw std::runtime_error("--bc must be 1, 2, 3 or 3b");
}
// coupling_name: the name written in the summary file.
const char* coupling_name(Coupling b) {
    switch (b) { case Coupling::BC1: return "BC1"; case Coupling::BC2: return "BC2"; case Coupling::BC3: return "BC3"; default: return "BC3b"; }
}

// bath_step: one Euler--Maruyama step of the bath at mode r (r = 0 or n-1 here):
//     c_r += -gamma D(c_r) dt + sqrt(2 gamma T) sqrt(dt) (N_x + i N_y),
// with N_x, N_y independent N(0,1) (one pair per trajectory, given by the caller).  The drift D depends
// on the coupling (see the header).  It uses the forces of the current state (so it must be called after
// compute_force) and recomputes them afterwards, because the second bath step of the same time step
// needs the forces of the updated state.
void bath_step(Chain& s, int r, double T, double gamma, Coupling bc, const Vec16& Nx, const Vec16& Ny, double dt, double sqrt_dt) {
    const Vec16& x = s.x[r];
    const Vec16& y = s.y[r];
    const Vec16& Fx = s.F_re[r];
    const Vec16& Fy = s.F_im[r];
    Vec16 dx, dy;                                     // the drift -gamma D, real and imaginary part
    switch (bc) {
    case Coupling::BC1:                               // D = F
        dx = -gamma * Fx; dy = -gamma * Fy; break;
    case Coupling::BC2: {                             // D = Re(conj(c) F)/|c|^2 c   (Re(conj(c) F) = x Fx + y Fy)
        const Vec16 I = x.square() + y.square();
        const Vec16 radial = (x * Fx + y * Fy) / I.max(1e-300);   // max(.,1e-300) avoids 0/0 at c = 0
        dx = -gamma * radial * x; dy = -gamma * radial * y; break; }
    case Coupling::BC3: {                             // D = |c|^2 c
        const Vec16 I = x.square() + y.square();
        dx = -gamma * I * x; dy = -gamma * I * y; break; }
    case Coupling::BC3b: {                            // D = |c|^2 c + Im(conj(c) F)/|c|^2 i c   (Im(conj(c) F) = x Fy - y Fx,
        const Vec16 I = x.square() + y.square();      //                                          i c = -y + i x)
        const Vec16 tangential = (x * Fy - y * Fx) / I.max(1e-300);
        dx = -gamma * (I * x - tangential * y); dy = -gamma * (I * y + tangential * x); break; }
    }
    // Euler--Maruyama update: drift times dt plus noise
    const double noise = std::sqrt(2.0 * gamma * T);
    s.x[r] += dx * dt + noise * sqrt_dt * Nx;
    s.y[r] += dy * dt + noise * sqrt_dt * Ny;
    compute_force(s);
}

// ============================================================================ 6. parameters
// Params: everything that can be set on the command line.  The initial values are the paper settings,
// so an option that is not given takes its paper value.
struct Params {
    int n = 50;                                   // number of modes
    double T1 = 4, Tn = 6;                        // bath temperatures at modes 1 and n
    double gamma = 0.1, dt = 1e-3;                // coupling strength, time step
    double burn = -1;                             // burn-in time; < 0 means max(2000, 3.2 n^2) (set in read_options)
    double measure = 20000;                       // measurement time (after the burn-in)
    std::string bc = "1";                         // boundary coupling: 1, 2, 3 or 3b
    std::string out;                              // output prefix; empty means BC<bc>_T<T1>_T<Tn>_n<n>
    int batches = 10;                             // number of batches of 16 trajectories
    int threads = omp_get_max_threads();          // OpenMP threads (all cores by default)
    std::uint64_t seed = 1;                       // base seed
    bool dry_run = false;                         // --dry-run 1: only print the parameters
};

// fmt: a number in the shortest form with up to 10 significant digits (4 -> "4", 0.001 -> "0.001").
std::string fmt(double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return b; }

// read_options: read "--name value" pairs into Params, then fill in the two defaults that depend on
// other parameters (burn-in time and output name) and check --bc.
Params read_options(int argc, char** argv) {
    Params p;
    for (int i = 1; i < argc; i += 2) {
        const std::string k = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        const std::string v = argv[i + 1];
        if (k == "--n") p.n = std::stoi(v);
        else if (k == "--T1") p.T1 = std::stod(v);
        else if (k == "--Tn") p.Tn = std::stod(v);
        else if (k == "--bc") p.bc = v;
        else if (k == "--gamma") p.gamma = std::stod(v);
        else if (k == "--dt") p.dt = std::stod(v);
        else if (k == "--burn") p.burn = std::stod(v);
        else if (k == "--measure") p.measure = std::stod(v);
        else if (k == "--batches") p.batches = std::stoi(v);
        else if (k == "--threads") p.threads = std::stoi(v);
        else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--out") p.out = v;
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    if (p.n < 2) throw std::runtime_error("n must be at least 2");
    if (p.burn < 0) p.burn = std::max(2000.0, 16.0 * p.n * p.n / 5.0);      // 16 n^2 / 5 = 3.2 n^2
    if (p.out.empty()) p.out = "BC" + p.bc + "_T" + fmt(p.T1) + "_T" + fmt(p.Tn) + "_n" + std::to_string(p.n);
    parse_coupling(p.bc);                                                    // stops here if --bc is wrong
    return p;
}

// print_parameters: print all parameters as one line, in the option format of the original program
// nls_cartesian v3 (the check scripts give this line to the original program and compare the outputs).
// --closed 0, --interval 0.1, --checkpoints 4 and --stream-currents 1 are fixed in this program.
void print_parameters(const Params& p) {
    std::cout << "command: --n " << p.n << " --T1 " << fmt(p.T1) << " --Tn " << fmt(p.Tn) << " --bc " << p.bc << " --closed 0"
              << " --gamma " << fmt(p.gamma) << " --dt " << fmt(p.dt) << " --burn " << fmt(p.burn) << " --measure " << fmt(p.measure)
              << " --interval " << fmt(SAMPLE_INTERVAL) << " --batches " << p.batches << " --threads " << p.threads
              << " --checkpoints " << CHECKPOINTS << " --seed " << p.seed << " --stream-currents 1 --out " << p.out << std::endl;
}

// ============================================================================ 7. simulation
// ThreadSums: running sums kept by each thread over all samples it recorded.  They are added up over the
// threads at the end.  (Keeping these totals per thread, rather than computing them from the
// per-trajectory sums, reproduces the rounding of the original program exactly.)
struct ThreadSums {
    std::vector<double> I, J;          // sum over samples of I_j and of J_j
    std::uint64_t samples = 0;         // number of recorded samples (one sample = one trajectory at one time)
    std::uint64_t failures = 0;        // midpoint steps that did not converge
    double iterations = 0;             // total number of midpoint iterations (for the mean in the summary)
    bool nonfinite = false;            // true if any recorded I_j was NaN or infinite
};

// run: the whole simulation for parameters p.
int run(const Params& p) {
    // ---- 7a. step counts
    const int n = p.n, L = p.batches * LANES;                         // L = number of trajectories
    const Coupling bc = parse_coupling(p.bc);
    const long long burn_steps = std::llround(p.burn / p.dt);        // time steps of the burn-in
    const long long stride = std::max(1LL, std::llround(SAMPLE_INTERVAL / p.dt));   // time steps between samples
    const long long samples_per_traj = std::llround(p.measure / SAMPLE_INTERVAL);  // samples per trajectory
    const long long chunk = (samples_per_traj + CHECKPOINTS - 1) / CHECKPOINTS;    // samples between checkpoints
    const double sqrt_dt = std::sqrt(p.dt);

    std::cout << "sec3_3_action_current n=" << n << " " << coupling_name(bc) << " T=(" << p.T1 << "," << p.Tn << ") gamma=" << p.gamma
              << " dt=" << p.dt << " burn=" << p.burn << " measure=" << p.measure << " trajectories=" << L << std::endl;

    // ---- 7b. initial state and random generators
    // one Chain and 16 generators per batch; every mode starts at c_j = I0^(1/2) with I0 = sqrt((T1+Tn)/(2n)),
    // i.e. x_j = I0^(1/2), y_j = 0 (all phases 0).  step_count[b] counts the time steps done by batch b.
    std::vector<Chain> chain;
    std::vector<std::array<Rng, LANES>> rng(p.batches);
    std::vector<long long> step_count(p.batches, 0);
    const double amplitude = std::sqrt(std::sqrt(0.5 * (p.T1 + p.Tn) / n));
    for (int b = 0; b < p.batches; ++b) {
        chain.emplace_back(n);
        for (int j = 0; j < n; ++j) { chain[b].x[j].setConstant(amplitude); chain[b].y[j].setZero(); }
        for (int l = 0; l < LANES; ++l) rng[b][l].seed(p.seed, static_cast<std::uint64_t>(b) * LANES + l);
    }

    // ---- 7c. sums for the averages
    // sums[t]: totals of thread t (see ThreadSums).
    // traj_I[j][k], traj_J[j][k]: sums of I_j and J_j over the samples of trajectory k = 16 b + lane; used
    // for the per-trajectory averages (standard errors and the _current_streams file).  No two threads
    // write the same entry, because a batch is advanced by one thread at a time.
    std::vector<ThreadSums> sums(p.threads);
    for (auto& t : sums) { t.I.assign(n, 0.0); t.J.assign(n, 0.0); }
    std::vector<std::vector<double>> traj_I(n, std::vector<double>(L, 0.0)), traj_J(n, std::vector<double>(L, 0.0));

    // ---- 7d. one time step, repeated: advance batch b by `steps` time steps
    // Each time step: midpoint step (Hamiltonian part), then the two bath steps.  The noise for both baths
    // is drawn first (mode 1, then mode n, always in this order), so the random numbers do not depend on
    // the bath order.  On even steps the bath at mode 1 goes first, on odd steps the bath at mode n; the
    // alternation makes the splitting symmetric on average between the two ends.
    auto advance = [&](int b, long long steps, ThreadSums& t, MidpointWork& w) {
        Chain& s = chain[b];
        Vec16 N1x, N1y, Nnx, Nny;                     // noise for mode 1 and mode n
        for (long long k = 0; k < steps; ++k) {
            int it = 0;
            if (!midpoint_step(s, w, p.dt, it)) ++t.failures;
            t.iterations += it;
            gaussian_pair(rng[b], N1x, N1y);
            gaussian_pair(rng[b], Nnx, Nny);
            if (step_count[b] % 2 == 0) {
                bath_step(s, 0, p.T1, p.gamma, bc, N1x, N1y, p.dt, sqrt_dt);
                bath_step(s, n - 1, p.Tn, p.gamma, bc, Nnx, Nny, p.dt, sqrt_dt);
            } else {
                bath_step(s, n - 1, p.Tn, p.gamma, bc, Nnx, Nny, p.dt, sqrt_dt);
                bath_step(s, 0, p.T1, p.gamma, bc, N1x, N1y, p.dt, sqrt_dt);
            }
            ++step_count[b];
        }
    };

    // ---- 7e. recording: add the current I_j and J_j of the 16 trajectories of batch b to the sums
    // (I and c^2 are up to date, because bath_step ends with compute_force).
    //   J_j = 4 Im(c_j^2 conj(c_{j-1}^2)) = 4 (Im c_j^2 Re c_{j-1}^2 - Re c_j^2 Im c_{j-1}^2)
    auto record = [&](int b, ThreadSums& t) {
        const Chain& s = chain[b];
        for (int l = 0; l < LANES; ++l) {
            const int k = b * LANES + l;              // trajectory number
            for (int j = 0; j < n; ++j) {
                const double I = s.I[j](l);
                if (!std::isfinite(I)) t.nonfinite = true;
                t.I[j] += I; traj_I[j][k] += I;
            }
            for (int j = 1; j < n; ++j) {             // bond between modes j-1 and j (zero-based)
                const double J = 4.0 * (s.sq_im[j](l) * s.sq_re[j - 1](l) - s.sq_re[j](l) * s.sq_im[j - 1](l));
                t.J[j] += J; traj_J[j][k] += J;
            }
            ++t.samples;
        }
    };
    // total number of samples recorded so far by all threads
    auto total_samples = [&]() { std::uint64_t N = 0; for (auto& t : sums) N += t.samples; return N; };

    const auto t0 = std::chrono::steady_clock::now();
    auto seconds = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

    // ---- 7f. burn-in: advance every batch by burn_steps without recording
    // Each thread gets its own MidpointWork; schedule(dynamic) hands the next batch to whichever thread is free.
#pragma omp parallel num_threads(p.threads)
    {
        MidpointWork w(n);
        ThreadSums& t = sums[omp_get_thread_num()];
#pragma omp for schedule(dynamic)
        for (int b = 0; b < p.batches; ++b) advance(b, burn_steps, t, w);
    }

    // ---- 7g. measurement: samples_per_traj times (advance by `stride` steps, then record), done in
    // CHECKPOINTS pieces.  After each piece the profile <I_j> so far is stored and printed, so one can see
    // whether it is still changing (written to the summary as checkpoint_k_mean_I).
    std::vector<std::vector<double>> checkpoint_profile;
    for (long long done = 0; done < samples_per_traj && (int)checkpoint_profile.size() < CHECKPOINTS;) {
        const long long todo = std::min(chunk, samples_per_traj - done);
#pragma omp parallel num_threads(p.threads)
        {
            MidpointWork w(n);
            ThreadSums& t = sums[omp_get_thread_num()];
#pragma omp for schedule(dynamic)
            for (int b = 0; b < p.batches; ++b)
                for (long long q = 0; q < todo; ++q) { advance(b, stride, t, w); record(b, t); }
        }
        done += todo;
        const std::uint64_t N = total_samples();
        std::vector<double> profile(n, 0.0);         // <I_j> over all samples so far
        for (int j = 0; j < n; ++j) { for (auto& t : sums) profile[j] += t.I[j]; profile[j] /= N; }
        checkpoint_profile.push_back(profile);
        std::cout << "checkpoint " << checkpoint_profile.size() << "/" << CHECKPOINTS << " samples=" << N << " elapsed=" << seconds() << "s  <I>=";
        for (int j = 0; j < n; ++j) std::cout << (j ? "," : "") << std::setprecision(5) << profile[j];
        std::cout << std::endl;
    }

    // ---- 7h. averages and standard errors
    // mean_I_j, mean_current_j: average over all samples of all trajectories.
    // Standard error: each trajectory gives one time average; the L trajectories are independent, so the
    // standard error of the mean is  sqrt( sample variance of the L time averages / L ).
    const std::uint64_t N = total_samples(), per_traj = N / L;
    std::uint64_t failures = 0; double iterations = 0; bool nonfinite = false;
    for (auto& t : sums) { failures += t.failures; iterations += t.iterations; nonfinite |= t.nonfinite; }
    for (int j = 0; j < n; ++j) for (auto& v : traj_I[j]) v /= per_traj;   // sums -> time averages per trajectory
    for (int j = 1; j < n; ++j) for (auto& v : traj_J[j]) v /= per_traj;
    auto standard_error = [&](const std::vector<double>& v) {
        double m = 0; for (double a : v) m += a; m /= L;
        double var = 0; for (double a : v) var += (a - m) * (a - m);
        return std::sqrt(var / (L - 1) / L);
    };

    // ---- 7i. OUT_summary.txt: one "name value(s)" line per item (the format read by the analysis scripts)
    std::ofstream f(p.out + "_summary.txt");
    f << std::setprecision(10);
    f << "program nls_cartesian\nn " << n << "\nbc " << coupling_name(bc) << "\nclosed 0\ncorotating 1"
      << "\nT1 " << p.T1 << "\nTn " << p.Tn << "\nT" << n << " " << p.Tn << "\ngamma " << p.gamma << "\ndt " << p.dt
      << "\nburn " << p.burn << "\nmeasure " << p.measure << "\ninterval " << SAMPLE_INTERVAL << "\nbatches " << p.batches
      << "\nstreams " << L << "\nseed " << p.seed << "\nthreads " << p.threads << "\nsamples " << N
      << "\nsnapshots_per_stream " << samples_per_traj << "\nmidpoint_failures " << failures
      << "\nmean_midpoint_iterations " << iterations / (static_cast<double>(p.batches) * (burn_steps + samples_per_traj * stride))
      << "\nnonfinite " << (nonfinite ? 1 : 0) << "\nelapsed_seconds " << seconds() << "\n";
    f << "mean_I";                                    // <I_j>, j = 1..n
    for (int j = 0; j < n; ++j) { double v = 0; for (auto& t : sums) v += t.I[j]; f << " " << v / N; }
    f << "\nse_I";
    for (int j = 0; j < n; ++j) f << " " << standard_error(traj_I[j]);
    f << "\nmean_current";                            // <J_j>, j = 2..n
    for (int j = 1; j < n; ++j) { double v = 0; for (auto& t : sums) v += t.J[j]; f << " " << v / N; }
    f << "\nse_current";
    for (int j = 1; j < n; ++j) f << " " << standard_error(traj_J[j]);
    f << "\n";
    for (std::size_t k = 0; k < checkpoint_profile.size(); ++k) {
        f << "checkpoint_" << k + 1 << "_mean_I";
        for (double v : checkpoint_profile[k]) f << " " << v;
        f << "\n";
    }

    // ---- 7j. OUT_current_streams.txt: per trajectory, time-averaged M and J_2 .. J_n
    // (used in Sec. 3.3 to check that the current is the same at every bond and to compute error bars)
    std::ofstream g(p.out + "_current_streams.txt");
    g << std::setprecision(10);
    g << "# per stream (row = batch*16 + lane): first column the time-averaged total action M, then the time-averaged bond current J_j = 4 I_{j-1} I_j sin(theta_j) for j = 2.." << n << "\n";
    for (int k = 0; k < L; ++k) {
        double M = 0; for (int j = 0; j < n; ++j) M += traj_I[j][k];
        g << M;
        for (int j = 1; j < n; ++j) g << " " << traj_J[j][k];
        g << "\n";
    }
    std::cout << "samples=" << N << " seconds=" << seconds() << " midpoint_failures=" << failures << " nonfinite=" << nonfinite << std::endl;
    return 0;
}

// ============================================================================ 8. main
// Read the options, print the full parameter line, then run (unless --dry-run 1).  Any error (unknown
// option, bad value) is printed and the program exits with status 1.
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
