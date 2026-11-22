// sec4_1_density.cpp  

// WHAT IT COMPUTES
//   Long trajectories of the short open chain (n = 3 or 4) with baths at modes 1 and n, in the reduced
//   coordinates of the density,
//     x = (I_1, ..., I_n, theta_2, ..., theta_n),   I_j = |c_j|^2,   theta_j = arg( c_j^2 conj(c_{j-1})^2 ),
//   a point of dimension D = 2n - 1.  Two uses (the mode is chosen by --queries):
//   * samples (default): every 0.005 time units the state x of every trajectory is stored; the file is
//     used as query points (where the density is estimated) and as independent validation data;
//   * KDE labels (--queries Q.bin): the states are not stored; instead every sample adds its Gaussian kernel
//     to a kernel density estimate at the fixed query points of Q.bin ("streamed KDE").  With 10^9 samples
//     this gives the density values used to train the network.
//     Kernel at a query q:  exp(-|x - q|^2 / 2), with the distance measured in units of hI for the actions
//     and hT for the angles (angle differences taken in (-pi, pi]); only pairs with distance <= rcut count.
//     Mirror images I_j -> -I_j of a sample are added near I_j = 0 (reflection boundary correction).
//     Estimate:  rho(q) = sum over samples and images / ( N (2 pi)^(D/2) hI^n hT^(n-1) ).
//
// MODEL, NUMERICAL SCHEME, PARALLEL STRUCTURE
//   As in sec3_1_action_profile.cpp with the canonical bath (BC1): implicit midpoint for the Hamiltonian
//   part in the co-rotating frame, Euler--Maruyama for the baths (order alternating), batches of 16
//   trajectories advanced together, batches distributed over OpenMP threads.
//   H = M^2 - 1/2 sum_j I_j^2 + 2 sum_{j>=2} I_{j-1} I_j cos(theta_j),   M = I_1 + ... + I_n,
//   F_j = dH/d conj(c_j) = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2),
//   bath modes: dc_r = i F_r dt - gamma F_r dt + sqrt(2 gamma T_r) (dW_x + i dW_y).
//
// OUTPUT  (prefix OUT, default n<n>_T<T1>_T<Tn>)
//   OUT_summary.txt   parameters, counters, mean_I, se_I, mean_current, checkpoint profiles
//                     (KDE mode: also the kernel parameters and the mean number of candidate queries per sample)
//   OUT.bin           samples mode: the states as float32.  Header: 8 bytes "NLS<D>D01\0" (e.g. "NLS5D01\0"),
//                     uint64 number of states, uint32 D, uint32 4 (bytes per value); then the states, one after
//                     the other.  With 2 threads the batches are written in the order the threads finish them.
//   OUT_kde.bin       KDE mode: 8 bytes "NLS<D>KDE1", uint64 number of queries, uint64 number of samples N,
//                     then the estimate rho(q) (float64) for each query in the order of Q.bin.
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 -Ieigen short_chain_density.cpp -o sec4_1_density
//   ./sec4_1_density --n 3 --seed 2026101032                              (samples, NESS (2,8))
//   ./sec4_1_density --n 3 --seed 12 --queries n3_queries_neq.bin         (KDE labels at the query points)
// Options and their defaults (the paper settings):
//   --n 3  --T1 2  --Tn 8  --gamma 0.1  --dt 1e-3  --burn 20  --seed 1  --out (see above)
//   samples mode:  --measure 156.25  --batches 2  --threads 2  --checkpoints 1
//   KDE mode:      --measure 19531.25  --batches 16  --threads (all cores)  --checkpoints 20
//                  --hI 0.2  --hT 0.2  --rcut 4  --reflect 1
//   --dry-run 1   print the parameters and stop.
//

// LAYOUT OF THIS FILE
//   1. constants and types
//   2. random numbers                 Rng, gaussian_pair
//   3. chain state and forces         Chain, compute_force
//   4. Hamiltonian part               MidpointWork, midpoint_step
//   5. bath part                      Coupling, bath_step
//   6. parameters                     Params, read_options, print_parameters
//   7. state files                    header_tag, read_states, write_states
//   8. streamed KDE                   QueryGrid, add_kernels, add_sample, write_kde
//   9. simulation                     run
//  10. main

#include <Eigen/Dense>   // Eigen::Array: fixed-size arrays with element-wise arithmetic
#include <omp.h>         // OpenMP: one batch per thread at a time

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================================ 1. constants and types
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int LANES = 16;                     // trajectories per batch
constexpr int MAX_ITERATIONS = 20;            // fixed-point iterations of the midpoint rule
constexpr double TOLERANCE = 2.0e-13;         // stop when max |change| < TOLERANCE
constexpr double SAMPLE_INTERVAL = 0.005;     // time between two recorded samples

// Vec16: 16 doubles, one per trajectory of a batch.  Arithmetic on Vec16 is element-wise, e.g.
// a * b multiplies the 16 pairs (a(l), b(l)), and a.square() squares each entry.
using Vec16 = Eigen::Array<double, LANES, 1>;
// Modes: one Vec16 per mode j = 0..n-1 (zero-based here; mode j in the code is mode j+1 in the paper).
// The aligned allocator is required by Eigen for vectors of fixed-size arrays.
using Modes = std::vector<Vec16, Eigen::aligned_allocator<Vec16>>;

// wrap_pi: the angle x moved into [-pi, pi] by a multiple of 2 pi
double wrap_pi(double x) { return x - std::round(x / TWO_PI) * TWO_PI; }

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
// The four boundary couplings (option --bc 1, 2, 3, 3b).
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
// Params: everything that can be set on the command line.  Values < 0 (or empty) are filled in by
// read_options with the paper setting of the chosen mode (samples or KDE).
struct Params {
    int n = 3;                                    // number of modes
    double T1 = 2, Tn = 8;                        // bath temperatures at modes 1 and n
    double gamma = 0.1, dt = 1e-3;                // coupling strength, time step
    double burn = 20;                             // burn-in time
    double measure = -1;                          // measurement time: 156.25 (samples) or 19531.25 (KDE)
    int batches = -1;                             // batches of 16 trajectories: 2 (samples) or 16 (KDE)
    int threads = -1;                             // OpenMP threads: 2 (samples) or all cores (KDE)
    int checkpoints = -1;                         // pieces of the measurement: 1 (samples) or 20 (KDE)
    std::string queries;                          // query file; non-empty = KDE mode
    double hI = 0.2, hT = 0.2, rcut = 4;          // kernel widths (actions, angles) and cutoff (in widths)
    int reflect = 1;                              // 1: add mirror images I_j -> -I_j
    std::string out;                              // output prefix; empty means n<n>_T<T1>_T<Tn>
    std::uint64_t seed = 1;                       // base seed
    bool dry_run = false;                         // --dry-run 1: only print the parameters
    bool kde() const { return !queries.empty(); }
};

// fmt: a number in the shortest form with up to 10 significant digits (4 -> "4", 0.001 -> "0.001").
std::string fmt(double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return b; }

// read_options: read "--name value" pairs into Params, then fill in the defaults of the mode.
Params read_options(int argc, char** argv) {
    Params p;
    for (int i = 1; i < argc; i += 2) {
        const std::string k = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        const std::string v = argv[i + 1];
        if (k == "--n") p.n = std::stoi(v);
        else if (k == "--T1") p.T1 = std::stod(v);
        else if (k == "--Tn") p.Tn = std::stod(v);
        else if (k == "--gamma") p.gamma = std::stod(v);
        else if (k == "--dt") p.dt = std::stod(v);
        else if (k == "--burn") p.burn = std::stod(v);
        else if (k == "--measure") p.measure = std::stod(v);
        else if (k == "--batches") p.batches = std::stoi(v);
        else if (k == "--threads") p.threads = std::stoi(v);
        else if (k == "--checkpoints") p.checkpoints = std::stoi(v);
        else if (k == "--queries") p.queries = v;
        else if (k == "--hI") p.hI = std::stod(v);
        else if (k == "--hT") p.hT = std::stod(v);
        else if (k == "--rcut") p.rcut = std::stod(v);
        else if (k == "--reflect") p.reflect = std::stoi(v);
        else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--out") p.out = v;
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    if (p.n < 2) throw std::runtime_error("n must be at least 2");
    if (p.kde() && 2 * p.n - 1 > 16) throw std::runtime_error("streamed KDE supports n <= 8 (D <= 16)");
    if (p.measure < 0) p.measure = p.kde() ? 19531.25 : 156.25;
    if (p.batches < 0) p.batches = p.kde() ? 16 : 2;
    if (p.threads < 0) p.threads = p.kde() ? omp_get_max_threads() : 2;
    if (p.checkpoints < 0) p.checkpoints = p.kde() ? 20 : 1;
    if (p.out.empty()) p.out = "n" + std::to_string(p.n) + "_T" + fmt(p.T1) + "_T" + fmt(p.Tn);
    return p;
}

// print_parameters: print all parameters as one line, in the option format of the original program
// nls_cartesian v2 (the check scripts give this line to the original program and compare the outputs).
void print_parameters(const Params& p) {
    std::cout << "command: --n " << p.n << " --bc 1 --T1 " << fmt(p.T1) << " --Tn " << fmt(p.Tn) << " --gamma " << fmt(p.gamma)
              << " --dt " << fmt(p.dt) << " --burn " << fmt(p.burn) << " --interval " << fmt(SAMPLE_INTERVAL)
              << " --measure " << fmt(p.measure) << " --batches " << p.batches << " --threads " << p.threads
              << " --checkpoints " << p.checkpoints;
    if (p.kde()) std::cout << " --queries " << p.queries << " --hI " << fmt(p.hI) << " --hT " << fmt(p.hT) << " --rcut " << fmt(p.rcut) << " --reflect " << p.reflect;
    else std::cout << " --dump 1";
    std::cout << " --seed " << p.seed << " --out " << p.out << std::endl;
}

// ============================================================================ 7. state files
// header_tag: the 8-byte file tag "NLS<D>D01\0" (states) or "NLS<D>KDE1" (KDE values).
std::string header_tag(int D, const char* kind) {
    std::string h = "NLS" + std::to_string(D) + kind;
    if (std::strcmp(kind, "D01") == 0) h.push_back('\0');
    if (h.size() != 8) throw std::runtime_error("header length");
    return h;
}

// read_states: read a state file (format above) and return the states as doubles, D values per state.
std::vector<double> read_states(const std::string& path, int D) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open: " + path);
    char magic[8]; in.read(magic, 8);
    const std::string want = header_tag(D, "D01");
    if (std::memcmp(magic, want.data(), 8) != 0) throw std::runtime_error("bad header (expected " + want.substr(0, 7) + "): " + path);
    std::uint64_t count; std::uint32_t dim, bytes;
    in.read(reinterpret_cast<char*>(&count), 8); in.read(reinterpret_cast<char*>(&dim), 4); in.read(reinterpret_cast<char*>(&bytes), 4);
    if (dim != static_cast<std::uint32_t>(D) || bytes != 4) throw std::runtime_error("unsupported layout: " + path);
    std::vector<float> f(count * D);
    in.read(reinterpret_cast<char*>(f.data()), f.size() * 4);
    if (!in) throw std::runtime_error("truncated: " + path);
    return std::vector<double>(f.begin(), f.end());
}

// write_states: write the float32 states x (D values per state) with the header above.
void write_states(const std::string& path, const std::vector<float>& x, int D) {
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("cannot write " + path);
    const std::string h = header_tag(D, "D01");
    o.write(h.data(), 8);
    const std::uint64_t count = x.size() / D; const std::uint32_t dim = D, bytes = 4;
    o.write(reinterpret_cast<const char*>(&count), 8); o.write(reinterpret_cast<const char*>(&dim), 4); o.write(reinterpret_cast<const char*>(&bytes), 4);
    o.write(reinterpret_cast<const char*>(x.data()), x.size() * 4);
}

// ============================================================================ 8. streamed KDE
// QueryGrid: the query points sorted into a grid of cells, so that for a sample only the queries in the
// neighbouring cells need to be checked.
//   Action axes: cells of width rcut*hI starting at 0.  Angle axes: [-pi, pi) cut into
//   floor(2 pi / (rcut*hT)) equal cells (periodic).  A query within the cutoff of a sample is then always in
//   one of the 3 neighbouring cells on every axis.
//   q holds the queries sorted by cell, `original` their positions in the query file, and the queries of
//   cell k are q[start[k] .. start[k+1]-1] (cells numbered row-major over the D axes).
struct QueryGrid {
    int D = 0, NA = 0, nq = 0;    // dimension, number of action axes (= n), number of queries
    std::vector<double> q;        // nq x D, sorted by cell
    std::vector<int> original;    // query i of q is query original[i] of the file
    std::vector<int> ncell;       // cells per axis
    std::vector<double> width;    // cell width per axis
    std::vector<int> start;       // start of every cell in q (ncells + 1 entries)
    long long ncells = 1;

    // cell index of the value v on axis d (angles periodic)
    int cell_of(int d, double v) const {
        if (d < NA) return static_cast<int>(std::floor(v / width[d]));
        int k = static_cast<int>(std::floor((v + PI) / width[d]));
        return ((k % ncell[d]) + ncell[d]) % ncell[d];
    }
    // build the grid from the raw queries (D values each, in file order)
    void build(const std::vector<double>& raw, const Params& p, int D_) {
        D = D_; NA = p.n; nq = static_cast<int>(raw.size() / D);
        ncell.assign(D, 1); width.assign(D, 1.0);
        double imax = 0;                                       // largest action of any query
        for (int i = 0; i < nq; ++i) for (int j = 0; j < NA; ++j) imax = std::max(imax, raw[i * D + j]);
        for (int d = 0; d < NA; ++d) { width[d] = p.rcut * p.hI; ncell[d] = static_cast<int>(std::floor(imax / width[d])) + 1; }
        for (int d = NA; d < D; ++d) { ncell[d] = std::max(1, static_cast<int>(std::floor(TWO_PI / (p.rcut * p.hT)))); width[d] = TWO_PI / ncell[d]; }
        ncells = 1;
        for (int d = 0; d < D; ++d) ncells *= ncell[d];
        // sort the queries by cell number (stable, so queries of the same cell keep their file order)
        std::vector<long long> key(nq);
        for (int i = 0; i < nq; ++i) { long long k = 0; for (int d = 0; d < D; ++d) k = k * ncell[d] + cell_of(d, raw[i * D + d]); key[i] = k; }
        original.resize(nq);
        for (int i = 0; i < nq; ++i) original[i] = i;
        std::stable_sort(original.begin(), original.end(), [&](int a, int b) { return key[a] < key[b]; });
        q.resize(raw.size());
        for (int i = 0; i < nq; ++i) for (int d = 0; d < D; ++d) q[i * D + d] = raw[original[i] * D + d];
        // start[k]: number of queries in cells before k (counting sort offsets)
        start.assign(ncells + 1, 0);
        for (int i = 0; i < nq; ++i) start[key[i] + 1]++;
        for (long long k = 0; k < ncells; ++k) start[k + 1] += start[k];
    }
};

// add_kernels: add exp(-d^2/2) of the point x to acc[i] for every query i within the cutoff (d <= rcut).
// The cells to visit are the 3 neighbours on every axis (fewer at the edges of the action range, and no
// repeats on angle axes with fewer than 3 cells).  For each axis the distance from x to the cell (0 if x is
// inside) is stored, and a combination of cells is skipped as soon as the partial sum of these distances
// exceeds the cutoff.  This gives the same sums as checking every query.  `candidates` counts the queries
// whose distance was computed.
void add_kernels(const double* x, const QueryGrid& g, const Params& p, double* acc, std::uint64_t& candidates) {
    const int D = g.D, NA = g.NA;
    const double invI = 1.0 / p.hI, invT = 1.0 / p.hT, r2 = p.rcut * p.rcut;
    int nbn[16]; int nb[16][3]; double bd[16][3];      // per axis: number of cells to visit, the cells, squared distances
    for (int d = 0; d < D; ++d) {
        const int base = g.cell_of(d, x[d]);
        nbn[d] = 0;
        for (int o = -1; o <= 1; ++o) {
            int k = base + o;
            if (d < NA) { if (k < 0 || k >= g.ncell[d]) continue; }
            else {
                k = ((k % g.ncell[d]) + g.ncell[d]) % g.ncell[d];
                bool dup = false;
                for (int t = 0; t < nbn[d]; ++t) dup |= (nb[d][t] == k);
                if (dup) continue;
            }
            double dist;
            if (d < NA) { const double lo = k * g.width[d], hi = lo + g.width[d]; dist = std::max(0.0, std::max(lo - x[d], x[d] - hi)) * invI; }
            else { const double centre = -PI + (k + 0.5) * g.width[d]; dist = std::max(0.0, std::fabs(wrap_pi(x[d] - centre)) - 0.5 * g.width[d]) * invT; }
            nb[d][nbn[d]] = k; bd[d][nbn[d]] = dist * dist; ++nbn[d];
        }
        if (nbn[d] == 0) return;                          // x is too far outside the query range
    }
    // depth-first walk over the combinations of cells: idx[d] = which of the nbn[d] cells on axis d,
    // part[d] = squared distance over axes 0..d-1, flat[d] = cell number over axes 0..d-1
    int idx[16]; double part[17]; long long flat[17];
    for (int d = 0; d < D; ++d) idx[d] = 0;
    part[0] = 0.0; flat[0] = 0;
    int d = 0;
    while (true) {
        if (d < 0) break;
        if (idx[d] >= nbn[d]) { idx[d] = 0; --d; if (d >= 0) ++idx[d]; continue; }
        const double partial = part[d] + bd[d][idx[d]];
        if (partial > r2) { ++idx[d]; continue; }         // this cell is beyond the cutoff
        part[d + 1] = partial; flat[d + 1] = flat[d] * g.ncell[d] + nb[d][idx[d]];
        if (d == D - 1) {                                 // a complete cell: check its queries
            const long long cell = flat[D];
            for (int i = g.start[cell]; i < g.start[cell + 1]; ++i) {
                const double* qq = &g.q[static_cast<std::size_t>(i) * D];
                double d2 = 0;
                for (int j = 0; j < NA; ++j) { const double u = (x[j] - qq[j]) * invI; d2 += u * u; }
                if (d2 > r2) continue;
                for (int j = NA; j < D; ++j) { const double u = wrap_pi(x[j] - qq[j]) * invT; d2 += u * u; }
                ++candidates;
                if (d2 <= r2) acc[i] += std::exp(-0.5 * d2);
            }
            ++idx[d];
        } else { ++d; idx[d] = 0; }
    }
}

// add_sample: add one sample x and, with --reflect 1, its mirror images (I_j -> -I_j for every subset of the
// actions that are closer to 0 than the kernel reach rcut*hI; only those images can reach a query).
void add_sample(const double* x, const QueryGrid& g, const Params& p, double* acc, std::uint64_t& candidates) {
    const int D = g.D, NA = g.NA;
    double y[16];
    if (!p.reflect) { add_kernels(x, g, p, acc, candidates); return; }
    const double reach = p.rcut * p.hI;
    int flip[16], nf = 0;
    for (int j = 0; j < NA; ++j) if (x[j] < reach) flip[nf++] = j;
    for (int m = 0; m < (1 << nf); ++m) {                // m runs over the subsets of the actions in `flip`
        for (int e = 0; e < D; ++e) y[e] = x[e];
        for (int t = 0; t < nf; ++t) if (m & (1 << t)) y[flip[t]] = -y[flip[t]];
        add_kernels(y, g, p, acc, candidates);
    }
}

// write_kde: normalise the kernel sums of all threads and write OUT_kde.bin (queries in file order).
void write_kde(const Params& p, const QueryGrid& g, const std::vector<std::vector<double>>& acc, std::uint64_t N) {
    const int D = g.D;
    std::vector<double> dens(g.nq, 0.0);
    const double norm = 1.0 / (static_cast<double>(N) * std::pow(2.0 * PI, 0.5 * D) * std::pow(p.hI, g.NA) * std::pow(p.hT, D - g.NA));
    for (int i = 0; i < g.nq; ++i) { double s = 0; for (const auto& a : acc) s += a[i]; dens[g.original[i]] = s * norm; }
    std::ofstream o(p.out + "_kde.bin", std::ios::binary);
    const std::string h = header_tag(D, "KDE1");
    o.write(h.data(), 8);
    const std::uint64_t Q = g.nq;
    o.write(reinterpret_cast<const char*>(&Q), 8); o.write(reinterpret_cast<const char*>(&N), 8);
    o.write(reinterpret_cast<const char*>(dens.data()), dens.size() * 8);
}

// ============================================================================ 9. simulation
// ThreadSums: running sums and outputs kept by each thread; added up (or concatenated) over the threads.
// (Keeping the totals per thread reproduces the rounding of the original program exactly.)
struct ThreadSums {
    std::vector<double> I, J;          // sum over samples of I_j and of J_j
    std::vector<double> kde;           // KDE mode: kernel sums at the (sorted) queries
    std::vector<float> states;         // samples mode: the recorded states
    std::uint64_t samples = 0;         // number of recorded samples (one sample = one trajectory at one time)
    std::uint64_t failures = 0;        // midpoint steps that did not converge
    std::uint64_t candidates = 0;      // KDE mode: queries whose distance was computed
    double iterations = 0;             // total number of midpoint iterations (for the mean in the summary)
    bool nonfinite = false;            // true if any recorded I_j was NaN or infinite
};

// run: the whole simulation for parameters p.
int run(const Params& p) {
    // ---- 9a. step counts and query grid
    const int n = p.n, D = 2 * n - 1, L = p.batches * LANES;          // L = number of trajectories
    const Coupling bc = Coupling::BC1;
    QueryGrid grid;
    if (p.kde()) grid.build(read_states(p.queries, D), p, D);
    const long long burn_steps = std::llround(p.burn / p.dt);        // time steps of the burn-in
    const long long stride = std::max(1LL, std::llround(SAMPLE_INTERVAL / p.dt));   // time steps between samples
    const long long samples_per_traj = std::llround(p.measure / SAMPLE_INTERVAL);  // samples per trajectory
    const int ncp = std::max(1, p.checkpoints);
    const long long chunk = (samples_per_traj + ncp - 1) / ncp;                    // samples between checkpoints
    const double sqrt_dt = std::sqrt(p.dt);

    std::cout << "sec4_1_density n=" << n << " T=(" << p.T1 << "," << p.Tn << ") gamma=" << p.gamma << " dt=" << p.dt << " burn=" << p.burn
              << " measure=" << p.measure << " trajectories=" << L << (p.kde() ? " kde queries=" + std::to_string(grid.nq) : " samples") << std::endl;

    // ---- 9b. initial state and random generators
    // one Chain and 16 generators per batch; every mode starts at c_j = I0^(1/2) with I0 = sqrt((T1+Tn)/(2n)).
    std::vector<Chain> chain;
    std::vector<std::array<Rng, LANES>> rng(p.batches);
    std::vector<long long> step_count(p.batches, 0);
    const double amplitude = std::sqrt(std::sqrt(0.5 * (p.T1 + p.Tn) / n));
    for (int b = 0; b < p.batches; ++b) {
        chain.emplace_back(n);
        for (int j = 0; j < n; ++j) { chain[b].x[j].setConstant(amplitude); chain[b].y[j].setZero(); }
        for (int l = 0; l < LANES; ++l) rng[b][l].seed(p.seed, static_cast<std::uint64_t>(b) * LANES + l);
    }

    // ---- 9c. sums: sums[t] per thread; traj_I[j][k] per trajectory k = 16 b + lane (for the standard errors)
    std::vector<ThreadSums> sums(p.threads);
    for (auto& t : sums) { t.I.assign(n, 0.0); t.J.assign(n, 0.0); if (p.kde()) t.kde.assign(grid.nq, 0.0); }
    std::vector<std::vector<double>> traj_I(n, std::vector<double>(L, 0.0));

    // ---- 9d. one time step, repeated: midpoint step, then the two bath steps (order alternating; the
    // noise is always drawn for mode 1 first, then for mode n)
    auto advance = [&](int b, long long steps, ThreadSums& t, MidpointWork& w) {
        Chain& s = chain[b];
        Vec16 N1x, N1y, Nnx, Nny;
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

    // ---- 9e. recording: for each trajectory of batch b form x = (I_1..I_n, theta_2..theta_n), add I_j and
    // J_j to the sums, and either store x (samples mode) or add its kernels (KDE mode).
    //   im = Im( c_j^2 conj(c_{j-1})^2 ),  re = Re( ... );  theta_j = atan2(im, re),  J_j = 4 im.
    auto record = [&](int b, ThreadSums& t) {
        const Chain& s = chain[b];
        double x[16];
        for (int l = 0; l < LANES; ++l) {
            const int k = b * LANES + l;              // trajectory number
            for (int j = 0; j < n; ++j) {
                x[j] = s.I[j](l);
                if (!std::isfinite(x[j])) t.nonfinite = true;
                t.I[j] += x[j]; traj_I[j][k] += x[j];
            }
            for (int j = 1; j < n; ++j) {
                const double re = s.sq_re[j](l) * s.sq_re[j - 1](l) + s.sq_im[j](l) * s.sq_im[j - 1](l);
                const double im = s.sq_im[j](l) * s.sq_re[j - 1](l) - s.sq_re[j](l) * s.sq_im[j - 1](l);
                t.J[j] += 4.0 * im;
                x[n + j - 1] = std::atan2(im, re);
            }
            ++t.samples;
            if (p.kde()) add_sample(x, grid, p, t.kde.data(), t.candidates);
            else for (int e = 0; e < D; ++e) t.states.push_back(static_cast<float>(x[e]));
        }
    };
    auto total_samples = [&]() { std::uint64_t N = 0; for (auto& t : sums) N += t.samples; return N; };

    const auto t0 = std::chrono::steady_clock::now();
    auto seconds = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

    // ---- 9f. burn-in
#pragma omp parallel num_threads(p.threads)
    {
        MidpointWork w(n);
        ThreadSums& t = sums[omp_get_thread_num()];
#pragma omp for schedule(dynamic)
        for (int b = 0; b < p.batches; ++b) advance(b, burn_steps, t, w);
    }

    // ---- 9g. measurement in ncp pieces; in KDE mode the estimate so far is written after every piece, so a
    // run that is stopped early still leaves a usable file
    std::vector<std::vector<double>> checkpoint_profile;
    for (long long done = 0; done < samples_per_traj && (int)checkpoint_profile.size() < ncp;) {
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
        std::vector<double> profile(n, 0.0);
        for (int j = 0; j < n; ++j) { for (auto& t : sums) profile[j] += t.I[j]; profile[j] /= N; }
        checkpoint_profile.push_back(profile);
        if (p.kde()) { std::vector<std::vector<double>> parts; for (auto& t : sums) parts.push_back(t.kde); write_kde(p, grid, parts, N); }
        std::cout << "checkpoint " << checkpoint_profile.size() << "/" << ncp << " samples=" << N << " elapsed=" << seconds() << "s  <I>=";
        for (int j = 0; j < n; ++j) std::cout << (j ? "," : "") << std::setprecision(5) << profile[j];
        std::cout << std::endl;
    }

    // ---- 9h. OUT_summary.txt
    const std::uint64_t N = total_samples(), per_traj = N / L;
    std::uint64_t failures = 0, candidates = 0; double iterations = 0; bool nonfinite = false;
    for (auto& t : sums) { failures += t.failures; candidates += t.candidates; iterations += t.iterations; nonfinite |= t.nonfinite; }
    std::ofstream f(p.out + "_summary.txt");
    f << std::setprecision(10);
    f << "program nls_cartesian\nn " << n << "\nbc " << coupling_name(bc) << "\nclosed 0\ncorotating 1"
      << "\nT1 " << p.T1 << "\nTn " << p.Tn << "\nT" << n << " " << p.Tn << "\ngamma " << p.gamma << "\ndt " << p.dt
      << "\nburn " << p.burn << "\nmeasure " << p.measure << "\ninterval " << SAMPLE_INTERVAL << "\nbatches " << p.batches
      << "\nstreams " << L << "\nseed " << p.seed << "\nthreads " << p.threads << "\nsamples " << N
      << "\nsnapshots_per_stream " << samples_per_traj << "\nmidpoint_failures " << failures
      << "\nmean_midpoint_iterations " << iterations / (static_cast<double>(p.batches) * (burn_steps + samples_per_traj * stride))
      << "\nnonfinite " << (nonfinite ? 1 : 0) << "\nelapsed_seconds " << seconds() << "\n";
    if (p.kde()) f << "queries " << p.queries << "\nhI " << p.hI << "\nhT " << p.hT << "\nrcut " << p.rcut << "\nreflect " << p.reflect
                   << "\nmean_kernel_candidates_per_sample " << static_cast<double>(candidates) / N << "\n";
    f << "mean_I";
    for (int j = 0; j < n; ++j) { double v = 0; for (auto& t : sums) v += t.I[j]; f << " " << v / N; }
    f << "\nse_I";                                    // standard error over the L independent trajectories
    for (int j = 0; j < n; ++j) {
        std::vector<double> a(traj_I[j]);
        double m = 0; for (auto& v : a) { v /= per_traj; m += v; }
        m /= L;
        double var = 0; for (double v : a) var += (v - m) * (v - m);
        f << " " << std::sqrt(var / (L - 1) / L);
    }
    f << "\nmean_current";
    for (int j = 1; j < n; ++j) { double v = 0; for (auto& t : sums) v += t.J[j]; f << " " << v / N; }
    f << "\n";
    for (std::size_t k = 0; k < checkpoint_profile.size(); ++k) {
        f << "checkpoint_" << k + 1 << "_mean_I";
        for (double v : checkpoint_profile[k]) f << " " << v;
        f << "\n";
    }
    f.close();

    // ---- 9i. samples mode: the states of all threads, thread 0 first
    if (!p.kde()) {
        std::vector<float> all;
        for (auto& t : sums) all.insert(all.end(), t.states.begin(), t.states.end());
        write_states(p.out + ".bin", all, D);
    }
    std::cout << "samples=" << N << " seconds=" << seconds() << " midpoint_failures=" << failures << " nonfinite=" << nonfinite << std::endl;
    return 0;
}

// ============================================================================ 10. main
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
