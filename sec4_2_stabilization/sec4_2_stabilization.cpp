// sec4_2_stabilization.cpp  --  Section 4.2: mechanism of stabilization of the short chain.
//
// WHAT IT COMPUTES
//   Independent trajectories ("streams") of the open chain with n modes and baths at modes 1 and n, started
//   from a chosen state: e.g. tiny end actions, a large total action M0 shared by the middle modes, and chosen
//   initial bond angles.  It records the time averages listed under OUTPUT and, every --interval time units, the
//   state (I_1..I_n, theta_2..theta_n) of every stream, from which the analysis gets M(t) and the angles.
//   A mode can also be "held" (--hold): it starts at exactly I = 0 and then stays at 0, and only its phase moves;
//   this gives the reduced boundary dynamics of the paper.
//
// MODEL (written with c_j = sqrt(I_j) exp(i phi_j))
//   H = M^2 - 1/2 sum_j I_j^2 + 2 sum_{j>=2} I_{j-1} I_j cos(theta_j),  theta_j = 2 (phi_j - phi_{j-1}),  M = sum_j I_j.
//   Hamiltonian part in action--angle variables (fixed frame):
//     dI_j/dt   = 4 I_j ( I_{j-1} sin(theta_j) - I_{j+1} sin(theta_{j+1}) ),
//     dphi_j/dt = (2M - I_j) + 2 I_{j-1} cos(theta_j) + 2 I_{j+1} cos(theta_{j+1})     (missing neighbours count 0).
//   Bath at modes r = 1, n (canonical bath, in Cartesian form):
//     dc_r = -gamma F_r dt + sqrt(2 gamma T_r) (dW_x + i dW_y),   F_r = (2M - I_r) c_r + 2 conj(c_r) (c_{r-1}^2 + c_{r+1}^2).
//
// NUMERICAL SCHEME (one step of size dt)
//   1. Hamiltonian part: classical fourth-order Runge--Kutta (RK4) step in (I, phi).  An action that is exactly 0
//      stays exactly 0, because dI_j/dt is proportional to I_j.  (The midpoint scheme of Sec. 3 works in Cartesian
//      coordinates, where a mode at I = 0 has no phase; for the held modes the action--angle form is needed.)
//   2. Bath part: convert to c = sqrt(I) e^{i phi}, one Euler--Maruyama step at mode 1 and one at mode n (order
//      alternating every step), and convert the two end modes back to (I, phi).  Only the end modes are converted
//      back: the bath does not change the others, and this keeps the phase of a mode at I = 0.
//
// RANDOM NUMBERS
//   Stream s uses its own std::mt19937_64 generator, seeded with seed * 1000003 + 7919 s + 17, and
//   std::normal_distribution.  The initial phases are always drawn (uniform on [0, 2 pi)); with --phases or --start
//   they are then replaced, so the noise of stream s is the same for every choice of initial phases.
//   Note: std::normal_distribution is implemented differently by different C++ libraries (libstdc++ of gcc,
//   libc++ of clang on a Mac), so results are reproduced only with the same library.
//
// OUTPUT  (prefix OUT, default sec42_n<n>[_M<M0>]_<start>)
//   OUT_summary.txt  parameters, then for every mode j the mean over streams and its standard error of
//                    meanI_j   time average of I_j,
//                    rate_j    time average of the Hamiltonian d log I_j / dt = 4 (I_{j-1} sin theta_j - I_{j+1} sin theta_{j+1}),
//                    Sleft_j   time average of -I_{j-1} sin theta_j,   Sright_j  time average of I_{j+1} sin theta_{j+1};
//                    and the number of steps with a non-finite state (nonfinite_steps).
//   OUT_dump.f32     (if --interval > 0) float32 states, layout (streams/16, snapshots, 16, 2n-1): columns
//                    I_1..I_n, then theta_2..theta_n in (-pi, pi]; snapshot q is the state at time (q+1) * interval.
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 sec4_2_stabilization.cpp -o sec4_2_stabilization
//   ./sec4_2_stabilization --n 4 --M0 100 --start anti --measure 800
//   ./sec4_2_stabilization --n 3 --M0 100 --start worst --gamma 1e-3 --measure ... --streams 16384   (see PAPER_RUNS.md)
// Options and their defaults (the paper settings):
//   --n (required)  --measure (required)  --T1 2  --Tn 8  --gamma 0.1  --dt 1e-4  --burn 0  --interval 0.1
//   --streams 128  --seed 701  --threads (all cores)  --out (see above)
//   initial state:
//   --M0 X                 I_1 = I_n = 0.01, the middle modes share X equally
//   --init a1,...,an       the initial actions (instead of --M0); default for every mode (T1+Tn)/(2n)
//   --start worst|anti|random   initial angles: worst = every bond angle 4 pi/3; anti = end bonds 4 pi/3, interior
//                          bonds pi; random = random phases (default)
//   --phases p1,...,pn     the initial phases (instead of --start)
//   --hold z1,z2,...       middle modes (one-based) that start and stay at I = 0
//   --dry-run 1            print the parameters and stop.
//   Angle convention of --start: the end bonds are measured outward, theta_1 = 2 (phi_1 - phi_2) and
//   theta_n = 2 (phi_n - phi_{n-1}); interior bonds left to right, theta_j = 2 (phi_j - phi_{j-1}).
//   --streams must be a multiple of 16 when --interval > 0 (the dump is written in blocks of 16 streams).
//
// The arithmetic is exactly that of the program used for the paper (rk4a_chain_v2.cpp): with the same compiler,
// flags and C++ library the output files are identical bit for bit, for any number of threads.
//
// LAYOUT OF THIS FILE
//   1. constants and types
//   2. model                       Model, cartesian_force, bath_step
//   3. Hamiltonian part            aa_rhs, rk4_step
//   4. one full step               to_cartesian, step
//   5. parameters                  Params, read_options, initial-state options, print_parameters
//   6. simulation                  run: streams, time averages, snapshots, output
//   7. main

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================================ 1. constants and types
using cd = std::complex<double>;
const double PI = 3.14159265358979323846;

// ============================================================================ 2. model
// Model: the chain length, bath temperatures and coupling, used by every step.
struct Model {
    int n;
    double T1, Tn, g;     // temperatures at modes 1 and n, bath coupling gamma
};

// cartesian_force: F_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2) for all j (missing neighbours 0).
void cartesian_force(const Model& m, const std::vector<cd>& c, std::vector<cd>& F) {
    double M = 0; for (auto z : c) M += std::norm(z);
    for (int j = 0; j < m.n; ++j) {
        cd N = 0;                                  // c_{j-1}^2 + c_{j+1}^2
        if (j > 0) N += c[j-1] * c[j-1];
        if (j < m.n - 1) N += c[j+1] * c[j+1];
        F[j] = (2 * M - std::norm(c[j])) * c[j] + 2.0 * N * std::conj(c[j]);
    }
}

// bath_step: Euler--Maruyama step of the bath at mode r with temperature T:
//   c_r += -gamma F_r h + sqrt(2 gamma T) (w1 + i w2),   (w1, w2) = sqrt(h) * two N(0,1) numbers.
// F is scratch space for the forces (computed here from the current c).
void bath_step(const Model& m, std::vector<cd>& c, int r, double T, double h, double w1, double w2, std::vector<cd>& F) {
    cartesian_force(m, c, F);
    c[r] += -m.g * F[r] * h + std::sqrt(2 * m.g * T) * cd(w1, w2);
}

// ============================================================================ 3. Hamiltonian part
// aa_rhs: right-hand side of the Hamiltonian flow in action--angle variables (see MODEL),
//   dI_j = 4 I_j ( I_{j-1} sin d_left - I_{j+1} sin d_right ),
//   dp_j = (2M - I_j) + 2 I_{j-1} cos d_left + 2 I_{j+1} cos d_right,
// with d_left = 2 (p_j - p_{j-1}) = theta_j and d_right = 2 (p_{j+1} - p_j) = theta_{j+1}.
void aa_rhs(const Model& m, const std::vector<double>& I, const std::vector<double>& p, std::vector<double>& dI, std::vector<double>& dp) {
    double M = 0; for (auto v : I) M += v;
    for (int j = 0; j < m.n; ++j) {
        double tI = 0, tp = 0;
        if (j > 0)       { double d = 2*(p[j]-p[j-1]); tI += I[j-1]*std::sin(d); tp += 2*I[j-1]*std::cos(d); }
        if (j < m.n - 1) { double d = 2*(p[j+1]-p[j]); tI -= I[j+1]*std::sin(d); tp += 2*I[j+1]*std::cos(d); }
        dI[j] = 4*I[j]*tI;
        dp[j] = (2*M - I[j]) + tp;
    }
}

// rk4_step: one classical RK4 step of size h for (I, p):
//   k1 = f(y), k2 = f(y + h/2 k1), k3 = f(y + h/2 k2), k4 = f(y + h k3),  y += h/6 (k1 + 2 k2 + 2 k3 + k4).
void rk4_step(const Model& m, std::vector<double>& I, std::vector<double>& p, double h) {
    const int n = m.n;
    std::vector<double> a1(n), b1(n), a2(n), b2(n), a3(n), b3(n), a4(n), b4(n), It(n), pt(n);
    aa_rhs(m, I, p, a1, b1);
    for (int j = 0; j < n; ++j) { It[j] = I[j] + 0.5*h*a1[j]; pt[j] = p[j] + 0.5*h*b1[j]; }
    aa_rhs(m, It, pt, a2, b2);
    for (int j = 0; j < n; ++j) { It[j] = I[j] + 0.5*h*a2[j]; pt[j] = p[j] + 0.5*h*b2[j]; }
    aa_rhs(m, It, pt, a3, b3);
    for (int j = 0; j < n; ++j) { It[j] = I[j] + h*a3[j]; pt[j] = p[j] + h*b3[j]; }
    aa_rhs(m, It, pt, a4, b4);
    for (int j = 0; j < n; ++j) {
        I[j] += h/6*(a1[j] + 2*a2[j] + 2*a3[j] + a4[j]);
        p[j] += h/6*(b1[j] + 2*b2[j] + 2*b3[j] + b4[j]);
    }
}

// ============================================================================ 4. one full step
// to_cartesian: c_j = sqrt(I_j) e^{i p_j}.
void to_cartesian(const Model& m, const std::vector<double>& I, const std::vector<double>& p, std::vector<cd>& c) {
    for (int j = 0; j < m.n; ++j) c[j] = std::polar(std::sqrt(I[j]), p[j]);
}

// step: one time step number k of size h.  dW[0..3] = sqrt(h) * N(0,1) for (mode 1 x, mode 1 y, mode n x, mode n y).
//   1. RK4 step of the Hamiltonian part;
//   2. bath steps in Cartesian form: mode 1 then mode n on even k, the other order on odd k;
//   3. back to (I, p) for the two end modes.  The new phase is taken as the old phase plus the change of arg c
//      reduced to [-pi, pi], so the phase does not jump by 2 pi.
void step(const Model& m, std::vector<double>& I, std::vector<double>& p, double h, const double* dW, long k,
          std::vector<cd>& c, std::vector<cd>& F) {
    const int n = m.n;
    rk4_step(m, I, p, h);
    std::vector<double> pb = p;                    // phases after the Hamiltonian part
    to_cartesian(m, I, p, c);
    if (k % 2 == 0) { bath_step(m, c, 0, m.T1, h, dW[0], dW[1], F); bath_step(m, c, n-1, m.Tn, h, dW[2], dW[3], F); }
    else            { bath_step(m, c, n-1, m.Tn, h, dW[2], dW[3], F); bath_step(m, c, 0, m.T1, h, dW[0], dW[1], F); }
    const int ends[2] = {0, n - 1};
    for (int e = 0; e < 2; ++e) {
        const int j = ends[e];
        I[j] = std::norm(c[j]);
        p[j] = pb[j] + std::remainder(std::arg(c[j]) - pb[j], 2*PI);
    }
}

// ============================================================================ 5. parameters
struct Params {
    int n = 0;                                    // number of modes (required)
    double T1 = 2, Tn = 8, gamma = 0.1;           // bath temperatures and coupling
    double dt = 1e-4, burn = 0, measure = -1;     // time step, burn-in time, measurement time (required)
    double interval = 0.1;                        // time between two snapshots of the dump (0: no dump)
    int streams = 128;                            // number of independent trajectories
    int threads = omp_get_max_threads();          // OpenMP threads
    unsigned long long seed = 701;                // base seed
    std::vector<int> hold;                        // held middle modes (one-based)
    std::vector<double> init, phases;             // initial actions and phases (empty: default / random)
    std::string init_text, phases_text;           // the same lists as text (printed by --dry-run)
    std::string M0, start;                        // --M0 and --start as given (they fill init and phases)
    std::string out;                              // output prefix
    bool dry_run = false;
};

// list readers for "a,b,c"
std::vector<int> int_list(const std::string& s) { std::vector<int> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) v.push_back(std::stoi(t)); return v; }
std::vector<double> double_list(const std::string& s) { std::vector<double> v; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) v.push_back(std::stod(t)); return v; }
std::string fmt(double v, const char* f = "%.10g") { char b[64]; std::snprintf(b, sizeof b, f, v); return b; }

// read_options: read "--name value" pairs, turn --M0 / --start into initial actions / phases, fill defaults.
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
        else if (k == "--interval") p.interval = std::stod(v);
        else if (k == "--streams") p.streams = std::stoi(v);
        else if (k == "--threads") p.threads = std::stoi(v);
        else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--hold") p.hold = int_list(v);
        else if (k == "--init") { p.init = double_list(v); p.init_text = v; }
        else if (k == "--phases") { p.phases = double_list(v); p.phases_text = v; }
        else if (k == "--M0") p.M0 = v;
        else if (k == "--start") p.start = v;
        else if (k == "--out") p.out = v;
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    if (p.n < 3) throw std::runtime_error("--n is required (n >= 3)");
    if (!(p.measure > 0)) throw std::runtime_error("--measure is required");
    const int n = p.n;
    // --M0 X: actions 0.01, X/(n-2), ..., X/(n-2), 0.01.  X/(n-2) is rounded to 10 significant digits, as in the
    // task tables of the paper runs.
    if (!p.M0.empty()) {
        if (!p.init.empty()) throw std::runtime_error("give --M0 or --init, not both");
        p.init_text = "0.01";
        for (int j = 0; j < n - 2; ++j) p.init_text += "," + fmt(std::stod(p.M0) / (n - 2));
        p.init_text += ",0.01";
        p.init = double_list(p.init_text);
    }
    // --start worst|anti: bond angles theta (end bonds 4 pi/3; interior bonds 4 pi/3 or pi) turned into phases
    // with phi_2 = 0:  phi_1 = theta_1 / 2 (end bond measured outward), phi_j = phi_{j-1} + theta_j / 2 for j >= 3.
    // The phases are rounded to 12 decimals, as in the task tables of the paper runs.
    if (!p.start.empty() && p.start != "random") {
        if (!p.phases.empty()) throw std::runtime_error("give --start or --phases, not both");
        const double pi = 3.141592653589793238462643383279502884;
        const double interior = p.start == "worst" ? 4 * pi / 3 : (p.start == "anti" ? pi : -1);
        if (interior < 0) throw std::runtime_error("--start must be worst, anti or random");
        std::vector<double> theta; theta.push_back(4 * pi / 3);
        for (int j = 0; j < n - 3; ++j) theta.push_back(interior);
        theta.push_back(4 * pi / 3);
        std::vector<double> phi(n, 0.0); phi[0] = theta[0] / 2;
        for (int j = 2; j < n; ++j) phi[j] = phi[j - 1] + theta[j - 1] / 2;
        for (int j = 0; j < n; ++j) p.phases_text += (j ? "," : "") + fmt(phi[j], "%.12f");
        p.phases = double_list(p.phases_text);
    }
    if (!p.init.empty() && (int)p.init.size() != n) throw std::runtime_error("--init needs n values");
    if (!p.phases.empty() && (int)p.phases.size() != n) throw std::runtime_error("--phases needs n values");
    for (int z : p.hold) if (z <= 1 || z >= n) throw std::runtime_error("only middle modes can be held");
    if (p.interval > 0 && p.streams % 16 != 0) throw std::runtime_error("--streams must be a multiple of 16 when --interval > 0");
    if (p.out.empty()) p.out = "sec42_n" + std::to_string(n) + (p.M0.empty() ? "" : "_M" + p.M0) + "_" + (p.start.empty() ? "run" : p.start);
    return p;
}

// print_parameters: all parameters in the option format of the original program rk4a_chain_v2 (used by the checks)
void print_parameters(const Params& p) {
    std::cout << "command: --n " << p.n << " --T1 " << fmt(p.T1) << " --Tn " << fmt(p.Tn) << " --gamma " << fmt(p.gamma) << " --dt " << fmt(p.dt)
              << " --burn " << fmt(p.burn) << " --measure " << fmt(p.measure) << " --interval " << fmt(p.interval) << " --streams " << p.streams
              << " --threads " << p.threads << " --seed " << p.seed;
    if (!p.hold.empty()) { std::cout << " --hold "; for (std::size_t i = 0; i < p.hold.size(); ++i) std::cout << (i ? "," : "") << p.hold[i]; }
    if (!p.init.empty()) std::cout << " --init " << p.init_text;
    if (!p.phases.empty()) std::cout << " --phases " << p.phases_text;
    std::cout << " --out " << p.out << std::endl;
}

// ============================================================================ 6. simulation
int run(const Params& P) {
    const Model m{P.n, P.T1, P.Tn, P.gamma};
    const int n = P.n, streams = P.streams, D = 2 * n - 1;
    std::vector<char> held(n, 0); for (int z : P.hold) held[z - 1] = 1;
    // ---- 6a. step counts and storage
    const long nburn = std::lround(P.burn / P.dt), nmeas = std::lround(P.measure / P.dt);
    const long every = P.interval > 0 ? std::lround(P.interval / P.dt) : 0;   // steps between snapshots
    const long nsnap = every > 0 ? nmeas / every : 0;                          // snapshots per stream
    std::vector<float> dump((size_t)nsnap * streams * D);                      // layout (snapshot, stream, D)
    std::vector<double> sumI((size_t)streams * n, 0), sumRate((size_t)streams * n, 0), sumS((size_t)streams * 2 * n, 0);
    std::vector<long> nonfinite(streams, 0);
    auto t0 = std::chrono::steady_clock::now();
    std::cout << "sec4_2_stabilization n=" << n << " T=(" << P.T1 << "," << P.Tn << ") gamma=" << P.gamma << " dt=" << P.dt
              << " measure=" << P.measure << " streams=" << streams << " snapshots=" << nsnap << std::endl;

    // ---- 6b. the streams (each one by one thread; schedule(dynamic, 1) gives the next stream to a free thread)
    omp_set_num_threads(P.threads);
#pragma omp parallel for schedule(dynamic, 1)
    for (int s = 0; s < streams; ++s) {
        std::mt19937_64 rng(P.seed * 1000003ULL + (unsigned long long)s * 7919ULL + 17ULL);
        std::normal_distribution<double> N01(0.0, 1.0);
        std::uniform_real_distribution<double> U(0.0, 2 * PI);
        // initial state: held modes at I = 0, otherwise --init or (T1+Tn)/(2n); phases random unless given
        std::vector<double> I(n), p(n), sn(n);
        for (int j = 0; j < n; ++j) {
            I[j] = held[j] ? 0.0 : (P.init.empty() ? 0.5 * (P.T1 + P.Tn) / n : P.init[j]);
            double u = U(rng);
            p[j] = P.phases.empty() ? u : P.phases[j];
        }
        std::vector<cd> c(n), F(n);
        long k = 0; double dW[4]; const double sq = std::sqrt(P.dt);
        auto one = [&]() { for (auto& v : dW) v = sq * N01(rng); step(m, I, p, P.dt, dW, k, c, F); ++k; };
        for (long it = 0; it < nburn; ++it) one();
        for (long it = 0; it < nmeas; ++it) {
            // ---- 6c. time averages (left-point rule): I_j, the rate 4 (I_{j-1} sin theta_j - I_{j+1} sin theta_{j+1})
            // and its two parts -I_{j-1} sin theta_j (Sleft) and I_{j+1} sin theta_{j+1} (Sright)
            for (int j = 1; j < n; ++j) sn[j] = std::sin(2 * (p[j] - p[j - 1]));
            for (int j = 0; j < n; ++j) {
                sumI[(size_t)s * n + j] += I[j];
                double Il = j > 0 ? I[j - 1] : 0, Ir = j + 1 < n ? I[j + 1] : 0;
                double sl = j > 0 ? sn[j] : 0, sr = j + 1 < n ? sn[j + 1] : 0;
                sumRate[(size_t)s * n + j] += 4 * (Il * sl - Ir * sr);
                sumS[(size_t)s * 2 * n + 2 * j] += -Il * sl;
                sumS[(size_t)s * 2 * n + 2 * j + 1] += Ir * sr;
            }
            // ---- 6d. snapshot: the state before the step (it + 1) when (it + 1) is a multiple of `every`
            if (every > 0 && (it + 1) % every == 0) {
                long q = (it + 1) / every - 1;
                if (q < nsnap) {
                    float* row = &dump[((size_t)q * streams + s) * D];
                    for (int j = 0; j < n; ++j) row[j] = (float)I[j];
                    for (int j = 1; j < n; ++j) row[n + j - 1] = (float)std::remainder(2 * (p[j] - p[j - 1]), 2 * PI);
                }
            }
            one();
            for (int j = 0; j < n; ++j) if (!std::isfinite(I[j]) || !std::isfinite(p[j])) { nonfinite[s] += 1; break; }
        }
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // ---- 6e. OUT_summary.txt (keys as in the original program, so the analysis scripts read it unchanged)
    std::ofstream f(P.out + "_summary.txt");
    f << "program rk4a_chain_v2 (schemes4 RK4_A + Cartesian bath; --phases)\nn " << n << "\nT1 " << P.T1 << "\nTn " << P.Tn << "\ngamma " << P.gamma
      << "\ndt " << P.dt << "\nburn " << P.burn << "\nmeasure " << P.measure << "\nstreams " << streams << "\nseed " << P.seed << "\nhold";
    for (int z : P.hold) f << " " << z;
    f << "\ninit"; for (double v : P.init) f << " " << v;
    f << "\nphases"; for (double v : P.phases) f << " " << v;
    long total_nonfinite = 0; for (long v : nonfinite) total_nonfinite += v;
    f << "\nsubstep_events 0\nfloor 0\nreflections_at_floor 0\nnonfinite_steps " << total_nonfinite << "\nelapsed_seconds " << el << "\n";
    // mean over streams of the time averages, and its standard error sqrt(sample variance / streams)
    auto stat = [&](std::vector<double>& a, int stride, int off, const std::string& name) {
        double mean = 0, m2 = 0;
        for (int s = 0; s < streams; ++s) { double v = a[(size_t)s * stride + off] / nmeas; mean += v; m2 += v * v; }
        mean /= streams; double var = (m2 / streams - mean * mean) * streams / (streams - 1);
        f << name << " " << mean << " " << std::sqrt(var / streams) << "\n";
    };
    for (int j = 0; j < n; ++j) {
        const std::string J = std::to_string(j + 1);
        stat(sumI, n, j, "meanI_" + J); stat(sumRate, n, j, "rate_" + J);
        stat(sumS, 2 * n, 2 * j, "Sleft_" + J); stat(sumS, 2 * n, 2 * j + 1, "Sright_" + J);
    }
    f.close();

    // ---- 6f. OUT_dump.f32: reorder (snapshot, stream) into (stream/16, snapshot, stream%16)
    if (nsnap > 0) {
        std::vector<float> outb((size_t)nsnap * streams * D);
        for (long q = 0; q < nsnap; ++q)
            for (int s = 0; s < streams; ++s) {
                int b = s / 16, lane = s % 16;
                std::memcpy(&outb[(((size_t)b * nsnap + q) * 16 + lane) * D], &dump[((size_t)q * streams + s) * D], D * sizeof(float));
            }
        std::ofstream g(P.out + "_dump.f32", std::ios::binary);
        g.write((const char*)outb.data(), outb.size() * sizeof(float));
    }
    std::printf("done %s seconds=%.1f nonfinite=%ld\n", P.out.c_str(), el, total_nonfinite);
    return 0;
}

// ============================================================================ 7. main
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
