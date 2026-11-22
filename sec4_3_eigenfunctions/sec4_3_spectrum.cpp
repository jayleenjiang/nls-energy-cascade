// sec4_3_spectrum.cpp  --  Section 4.3: slow spectrum of the generator of the three-mode chain (gEDMD).
//
// WHAT IT COMPUTES
//   Long trajectories of the three-mode chain (n = 3) with baths at modes 1 and 3, sampled every 0.1 time
//   units in the reduced coordinates
//     x = (I_1, I_2, I_3, theta_1, theta_3),   theta_1 = 2 (phi_1 - phi_2),  theta_3 = 2 (phi_3 - phi_2).
//   Every sample is used at once and then dropped ("streaming"): for a dictionary of K functions phi_k(x) it
//   adds up
//     A = sum over samples of phi(x) phi(x)^T,      B = sum over samples of phi(x) (L phi)(x)^T,
//   where L is the generator of the SDE.  The analysis (analyze_step2.py) solves the generalized eigenvalue
//   problem B^T v = lambda A v (generator EDMD): the eigenvalues lambda approximate the slow spectrum of L, and
//   v the coefficients of the eigenfunctions.  The samples are split into 4 groups (batch b goes to group b mod 4)
//   so the analysis can compare eigenvalues and eigenfunctions between groups (out-of-sample check).
//   The program also accumulates the autocorrelation functions of 6 observables up to 50 sampling intervals.
//
// MODEL (n = 3; the same chain as in Sec. 3, canonical baths)
//   H = M^2 - 1/2 sum_j I_j^2 + 2 I_1 I_2 cos(theta_1) + 2 I_2 I_3 cos(theta_3),   M = I_1 + I_2 + I_3.
//   Cartesian form used for the simulation: dc_j = i F_j dt (j = 2),  dc_r = (i - gamma) F_r dt
//   + sqrt(2 gamma T_r) (dW_x + i dW_y) (r = 1, 3),  F_j = (2M - I_j) c_j + 2 conj(c_j) (c_{j-1}^2 + c_{j+1}^2).
//   Generator in x (paper (2.18)):  L g = f . grad g + sum_i a_i d_i^2 g  with the drift f given in eval_dict and
//     a_{I1} = 4 gamma T1 I_1,  a_{I3} = 4 gamma T3 I_3,  a_{theta1} = 4 gamma T1 / I_1,  a_{theta3} = 4 gamma T3 / I_3.
//
// DICTIONARY
//   phi = prod_j (I_j / s_j)^{p_j} (I_j / (s_j + I_j))^{q_j} * cos or sin (m1 theta_1 + m3 theta_3),
//   with p_1 + p_2 + p_3 <= d, 0 <= m1 <= M, -M <= m3 <= M (and the constant mode), and
//   q = (|m1|, |m1 + m3|, |m3|), which makes every function smooth in the Cartesian amplitudes c_j.
//   s_j are scales of the actions (about their mean values).  With d = 5 and M = 3, K = 2744.  Dictionaries with
//   smaller d, M are subsets of this one, so only the largest is accumulated.
//
// NUMERICAL SCHEME (one step of size dt)
//   Implicit midpoint rule for the Hamiltonian part (fixed frame, fixed-point iteration, tolerance 2e-13, at most
//   20 iterations), then one Euler--Maruyama step at mode 1 and one at mode 3 (order alternating).  Batches of 16
//   trajectories (Vec16 holds one value per trajectory); each batch is run by one OpenMP thread.
//
// OUTPUT  (prefix OUT, default spectrum_T<T1>_T<T3>_dt<dt>)
//   OUT.bin    "NLS3GED1", int64 K, int64 number of groups, int64 number of lags; then for every group: number of
//              samples, A (K x K), B (K x K, column-major, B(k,l) = sum phi_k L phi_l), sum of L phi, sum of
//              (L phi)^2, sum of x, sum of the observables, lagged products of the observables, counts per lag.
//   OUT.json   parameters, run statistics and the list of dictionary terms (p1,p2,p3,q1,q2,q3,m1,m3,sine).
//
// USAGE  (compile and run inside a batch job, not on a login node)
//   g++ -O3 -march=native -fopenmp -std=c++17 -Ieigen sec4_3_spectrum.cpp -o sec4_3_spectrum
//   ./sec4_3_spectrum --seed 5101                                   (NESS (2,8), dt 1e-3)
//   ./sec4_3_spectrum --T1 5 --T3 5 --dt 2.5e-4 --seed 5104         (equilibrium, smaller dt)
// Options and their defaults (the paper settings):
//   --T1 2  --T3 8  --gamma 0.1  --dt 1e-3  --burn 50  --measure 10000  --snapshot 0.1  --batches 16 (256 trajectories)
//   --groups 4  --d 5  --M 3  --seed 1  --threads (all cores)  --out (see above)
//   --s1 --s2 --s3   action scales; default 1.32, 1.71, 1.65 for (T1,T3) = (2,8) and 1.48, 1.68, 1.48 for (5,5)
//                    (about the mean actions); they must be given for other temperatures.
//   --dry-run 1      print the parameters and stop.
//
// The arithmetic is exactly that of the program used for the paper (NLS3_gedmd_stream.cpp): with the same compiler,
// flags and 1 thread the output files are identical bit for bit.  (With more threads, the groups are summed in
// the order the batches finish, so the sums can differ in the last digit.)
//
// LAYOUT OF THIS FILE
//   1. constants and types
//   2. random numbers            Rng, gaussian_pair
//   3. chain state and forces    State, compute_force
//   4. time step                 midpoint_step, bath_step, wrap_angle, finite_state
//   5. dictionary and generator  Term, make_dictionary, Dict, eval_dict
//   6. accumulators              observables, Acc, Feeder
//   7. parameters and output     Params, read_options, print_parameters, write_output
//   8. simulation                run
//   9. main

#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
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

namespace {   // everything has internal linkage, as in the original program (this also keeps gcc's rounding identical)

// ============================================================================ 1. constants and types
constexpr int N = 3;                          // number of modes
constexpr int LANES = 16;                     // trajectories per batch
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_ITERATIONS = 20;            // fixed-point iterations of the midpoint rule
constexpr double TOLERANCE = 2.0e-13;         // stop when max |change| < TOLERANCE
constexpr int NOBS_ACF = 6;                   // observables for the autocorrelation
constexpr int NLAG = 51;                      // lags 0..50 sampling intervals

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
// State of one batch: x, y = real and imaginary parts of c_j; the rest is filled by compute_force:
// I = |c|^2, (sq_re, sq_im) = c^2, (F_re, F_im) = F, M = I_1 + I_2 + I_3.
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
//   c_r += -gamma F_r dt + sqrt(2 gamma T dt) (N_x + i N_y),
// then recompute the forces (the second bath step of the time step needs the updated ones).
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

// ============================================================================ 5. dictionary and generator
// Term: one dictionary function (see DICTIONARY in the header): powers p, q, angle numbers m1, m3, cos or sin.
struct Term { int p[3]; int q[3]; int m1, m3; int sine; };

// make_dictionary: all terms for the action degree d and angle order M, in the order of gedmd_n3.py:
// angle modes outer (constant first, then (m1, m3) with cos then sin), action powers inner.
std::vector<Term> make_dictionary(int d, int M) {
    std::vector<std::array<int, 3>> acts;                     // (p1, p2, p3) with p1 + p2 + p3 <= d
    for (int a = 0; a <= d; ++a) for (int b = 0; b <= d; ++b) for (int c = 0; c <= d; ++c)
        if (a + b + c <= d) acts.push_back({a, b, c});
    struct Mode { int m1, m3, sine; };
    std::vector<Mode> modes{{0, 0, 0}};                       // the constant angle mode
    for (int m1 = 0; m1 <= M; ++m1)
        for (int m3 = -M; m3 <= M; ++m3) {
            if (m1 == 0 && m3 <= 0) continue;                 // (m1, m3) and (-m1, -m3) give the same functions
            modes.push_back({m1, m3, 0});
            modes.push_back({m1, m3, 1});
        }
    std::vector<Term> t;
    for (auto& md : modes)
        for (auto& p : acts) {
            Term x;
            x.p[0] = p[0]; x.p[1] = p[1]; x.p[2] = p[2];
            x.q[0] = std::abs(md.m1); x.q[1] = std::abs(md.m1 + md.m3); x.q[2] = std::abs(md.m3);
            x.m1 = md.m1; x.m3 = md.m3; x.sine = md.sine;
            t.push_back(x);
        }
    return t;
}

// Dict: the dictionary, its size K, the largest powers (for the power tables) and the action scales s_j.
struct Dict {
    std::vector<Term> terms;
    int K = 0, dmax = 0, qmax = 0;
    double s[3];
};

// eval_dict: phi_k(x) and (L phi_k)(x) for all k at one state x = (I1, I2, I3, theta1, theta3).
//   Drift of the SDE in x (Ito form), with dH_j = dH/dI_j:
//     f_I1 = 4 I1 I2 sin th1 + 2 gamma (2 T1 - I1 dH_1),   f_I2 = -4 I2 (I1 sin th1 + I3 sin th3),
//     f_I3 = 4 I2 I3 sin th3 + 2 gamma (2 T3 - I3 dH_3),
//     f_th1 = 2 (dH_1 - dH_2) + 4 gamma I2 sin th1,        f_th3 = 2 (dH_3 - dH_2) + 4 gamma I2 sin th3.
//   Each term is phi = U(I) * T(angles) with U = u_1(I1) u_2(I2) u_3(I3), u_j = z_j^p h_j^q, z_j = I_j/s_j,
//   h_j = I_j/(s_j + I_j); the derivatives of u_j come from power tables computed once per state.
//   L phi = f . grad phi + a_I1 d^2phi/dI1^2 + a_I3 d^2phi/dI3^2 - (a_th1 m1^2 + a_th3 m3^2) phi.
void eval_dict(const Dict& D, const double* x, double T1, double T3, double gam, double* phi, double* Lphi) {
    const double I[3] = {x[0], x[1], x[2]};
    const double t1 = x[3], t3 = x[4];
    const double c1 = std::cos(t1), s1 = std::sin(t1), c3 = std::cos(t3), s3 = std::sin(t3);
    const double M = I[0] + I[1] + I[2];
    const double dH1 = 2 * M - I[0] + 2 * I[1] * c1;
    const double dH2 = 2 * M - I[1] + 2 * I[0] * c1 + 2 * I[2] * c3;
    const double dH3 = 2 * M - I[2] + 2 * I[1] * c3;
    double f[5];
    f[0] = 4 * I[0] * I[1] * s1 + 2 * gam * (2 * T1 - I[0] * dH1);
    f[1] = -4 * I[1] * (I[0] * s1 + I[2] * s3);
    f[2] = 4 * I[1] * I[2] * s3 + 2 * gam * (2 * T3 - I[2] * dH3);
    f[3] = 2 * (dH1 - dH2) + 4 * gam * I[1] * s1;
    f[4] = 2 * (dH3 - dH2) + 4 * gam * I[1] * s3;
    const double aI1 = 4 * gam * T1 * I[0], aI3 = 4 * gam * T3 * I[2];
    const double at1 = 4 * gam * T1 / I[0], at3 = 4 * gam * T3 / I[2];
    // power tables per action: zp = z^p, zp1 = d(z^p)/dI, zp2 = d^2(z^p)/dI^2, and the same (hq...) for h^q,
    // using dh/dI = s/(s+I)^2 and d^2h/dI^2 = -2 s/(s+I)^3
    double zp[3][16], zp1[3][16], zp2[3][16], hq[3][16], hq1[3][16], hq2[3][16];
    for (int j = 0; j < 3; ++j) {
        const double sj = D.s[j], z = I[j] / sj, h = I[j] / (sj + I[j]);
        const double h1 = sj / ((sj + I[j]) * (sj + I[j])), h2 = -2.0 * sj / ((sj + I[j]) * (sj + I[j]) * (sj + I[j]));
        double zpow[17], hpow[17];
        zpow[0] = 1; hpow[0] = 1;
        for (int k = 1; k <= 16; ++k) { zpow[k] = zpow[k - 1] * z; hpow[k] = hpow[k - 1] * h; }
        for (int p = 0; p <= D.dmax; ++p) {
            zp[j][p] = zpow[p];
            zp1[j][p] = p >= 1 ? p * zpow[p - 1] / sj : 0.0;
            zp2[j][p] = p >= 2 ? p * (p - 1) * zpow[p - 2] / (sj * sj) : 0.0;
        }
        for (int q = 0; q <= D.qmax; ++q) {
            hq[j][q] = hpow[q];
            hq1[j][q] = q >= 1 ? q * hpow[q - 1] * h1 : 0.0;
            hq2[j][q] = (q >= 2 ? q * (q - 1) * hpow[q - 2] * h1 * h1 : 0.0) + (q >= 1 ? q * hpow[q - 1] * h2 : 0.0);
        }
    }
    for (int k = 0; k < D.K; ++k) {
        const Term& t = D.terms[k];
        double u[3], du[3], d2u[3];                        // u_j and its first two derivatives (product rule)
        for (int j = 0; j < 3; ++j) {
            const int p = t.p[j], q = t.q[j];
            u[j] = zp[j][p] * hq[j][q];
            du[j] = zp1[j][p] * hq[j][q] + zp[j][p] * hq1[j][q];
            d2u[j] = zp2[j][p] * hq[j][q] + 2.0 * zp1[j][p] * hq1[j][q] + zp[j][p] * hq2[j][q];
        }
        const double arg = t.m1 * t1 + t.m3 * t3;
        const double ca = std::cos(arg), sa = std::sin(arg);
        const double T = t.sine ? sa : ca, dT = t.sine ? ca : -sa;   // angle factor and its derivative in arg
        const double U = u[0] * u[1] * u[2];
        const double ph = U * T;
        const double gI1 = du[0] * u[1] * u[2] * T, gI2 = u[0] * du[1] * u[2] * T, gI3 = u[0] * u[1] * du[2] * T;
        const double gt1 = U * t.m1 * dT, gt3 = U * t.m3 * dT;
        phi[k] = ph;
        Lphi[k] = f[0] * gI1 + f[1] * gI2 + f[2] * gI3 + f[3] * gt1 + f[4] * gt3
                + aI1 * d2u[0] * u[1] * u[2] * T + aI3 * u[0] * u[1] * d2u[2] * T
                - (at1 * t.m1 * t.m1 + at3 * t.m3 * t.m3) * ph;
    }
}

// ============================================================================ 6. accumulators
// observables for the autocorrelation: I2, I1 + I3, I1 - I3, h1 h2 cos th1, h2 h3 cos th3, cos th1
void observables(const double* x, const double* s, double* o) {
    const double h1 = x[0] / (s[0] + x[0]), h2 = x[1] / (s[1] + x[1]), h3 = x[2] / (s[2] + x[2]);
    o[0] = x[1]; o[1] = x[0] + x[2]; o[2] = x[0] - x[2];
    o[3] = h1 * h2 * std::cos(x[3]); o[4] = h2 * h3 * std::cos(x[4]); o[5] = std::cos(x[3]);
}

// Acc: the sums of one group (or of one batch before it is added to its group).
struct Acc {
    int K = 0;
    double n = 0;                                 // number of samples
    Eigen::MatrixXd A, B;                         // sum phi phi^T (lower triangle only), sum phi (L phi)^T
    Eigen::VectorXd sL, sL2;                      // sum L phi, sum (L phi)^2
    double sx[5] = {0, 0, 0, 0, 0};               // sum x
    double so[NOBS_ACF] = {0};                    // sum of the observables
    double acf[NOBS_ACF][NLAG] = {{0}};           // sum o(t) o(t - lag)
    double acfn[NLAG] = {0};                      // number of pairs per lag
    void init(int K_) {
        K = K_; A = Eigen::MatrixXd::Zero(K, K); B = Eigen::MatrixXd::Zero(K, K);
        sL = Eigen::VectorXd::Zero(K); sL2 = Eigen::VectorXd::Zero(K);
    }
    void add(const Acc& o) {
        n += o.n; A += o.A; B += o.B; sL += o.sL; sL2 += o.sL2;
        for (int i = 0; i < 5; ++i) sx[i] += o.sx[i];
        for (int i = 0; i < NOBS_ACF; ++i) { so[i] += o.so[i]; for (int l = 0; l < NLAG; ++l) acf[i][l] += o.acf[i][l]; }
        for (int l = 0; l < NLAG; ++l) acfn[l] += o.acfn[l];
    }
};

// Feeder: collects phi and L phi of up to 256 samples as columns of two matrices and then adds them to A and B with
// one matrix product each (much faster than one sample at a time).
struct Feeder {
    const Dict* D; double T1, T3, gam; Acc* acc;
    static constexpr int NB = 256;
    Eigen::MatrixXd Phi, LPhi; int fill = 0;
    void init(const Dict* D_, double T1_, double T3_, double g_, Acc* a) {
        D = D_; T1 = T1_; T3 = T3_; gam = g_; acc = a;
        Phi.resize(D->K, NB); LPhi.resize(D->K, NB); fill = 0;
    }
    void push(const double* x) {
        eval_dict(*D, x, T1, T3, gam, Phi.col(fill).data(), LPhi.col(fill).data());
        for (int i = 0; i < 5; ++i) acc->sx[i] += x[i];
        acc->n += 1;
        if (++fill == NB) flush();
    }
    void flush() {
        if (fill == 0) return;
        auto P = Phi.leftCols(fill); auto L = LPhi.leftCols(fill);
        acc->A.selfadjointView<Eigen::Lower>().rankUpdate(P);     // A += P P^T (lower triangle)
        acc->B.noalias() += P * L.transpose();                   // B += P L^T
        acc->sL += L.rowwise().sum();
        acc->sL2 += L.array().square().matrix().rowwise().sum();
        fill = 0;
    }
};

// ============================================================================ 7. parameters and output
struct Params {
    double T1 = 2, T3 = 8, gamma = 0.1, dt = 1e-3, burn = 50, measure = 10000, snapshot = 0.1;
    int batches = 16, groups = 4, threads = omp_get_max_threads(), d = 5, M = 3;
    std::uint64_t seed = 1;
    double s[3] = {-1, -1, -1};                   // action scales (< 0: paper value for (2,8) or (5,5))
    std::string out;                              // output prefix; empty: spectrum_T<T1>_T<T3>_dt<dt>
    bool dry_run = false;
};

std::string fmt(double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return b; }

// read_options: read "--name value" pairs, fill in the scales and the output name, check values.
Params read_options(int argc, char** argv) {
    Params p;
    for (int i = 1; i < argc; i += 2) {
        const std::string k = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        const std::string v = argv[i + 1];
        if (k == "--out") p.out = v;
        else if (k == "--T1") p.T1 = std::stod(v); else if (k == "--T3") p.T3 = std::stod(v);
        else if (k == "--gamma") p.gamma = std::stod(v); else if (k == "--dt") p.dt = std::stod(v);
        else if (k == "--burn") p.burn = std::stod(v); else if (k == "--measure") p.measure = std::stod(v);
        else if (k == "--snapshot") p.snapshot = std::stod(v); else if (k == "--batches") p.batches = std::stoi(v);
        else if (k == "--groups") p.groups = std::stoi(v); else if (k == "--threads") p.threads = std::stoi(v);
        else if (k == "--d") p.d = std::stoi(v); else if (k == "--M") p.M = std::stoi(v);
        else if (k == "--seed") p.seed = std::stoull(v);
        else if (k == "--s1") p.s[0] = std::stod(v); else if (k == "--s2") p.s[1] = std::stod(v); else if (k == "--s3") p.s[2] = std::stod(v);
        else if (k == "--dry-run") p.dry_run = (v == "1");
        else throw std::runtime_error("unknown option " + k);
    }
    const double paper28[3] = {1.32, 1.71, 1.65}, paper55[3] = {1.48, 1.68, 1.48};
    for (int j = 0; j < 3; ++j) if (p.s[j] < 0) {
        if (p.T1 == 2 && p.T3 == 8) p.s[j] = paper28[j];
        else if (p.T1 == 5 && p.T3 == 5) p.s[j] = paper55[j];
        else throw std::runtime_error("give --s1 --s2 --s3 (about the mean actions) for these temperatures");
    }
    if (p.d > 15 || 2 * p.M > 15) throw std::runtime_error("d or M too large");
    if (p.out.empty()) p.out = "spectrum_T" + fmt(p.T1) + "_T" + fmt(p.T3) + "_dt" + fmt(p.dt);
    return p;
}

// print_parameters: all parameters in the option format of the original program NLS3_gedmd_stream (for the checks)
void print_parameters(const Params& p) {
    std::cout << "command: --T1 " << fmt(p.T1) << " --T3 " << fmt(p.T3) << " --gamma " << fmt(p.gamma) << " --dt " << fmt(p.dt)
              << " --burn " << fmt(p.burn) << " --measure " << fmt(p.measure) << " --snapshot " << fmt(p.snapshot)
              << " --batches " << p.batches << " --groups " << p.groups << " --d " << p.d << " --M " << p.M << " --seed " << p.seed
              << " --threads " << p.threads << " --s1 " << fmt(p.s[0]) << " --s2 " << fmt(p.s[1]) << " --s3 " << fmt(p.s[2])
              << " --out " << p.out << std::endl;
}

// write_output: OUT.bin (the sums of every group) and OUT.json (parameters and dictionary), formats in the header.
void write_output(const Params& c, const Dict& D, const std::vector<Acc>& G, double seconds,
                  std::uint64_t fails, std::uint64_t nonfinite, double minact) {
    std::ofstream f(c.out + ".bin", std::ios::binary);
    f.write("NLS3GED1", 8);
    const std::int64_t hdr[3] = {D.K, (std::int64_t)G.size(), NLAG};
    f.write((const char*)hdr, sizeof hdr);
    for (const auto& g : G) {
        f.write((const char*)&g.n, 8);
        Eigen::MatrixXd A = g.A.selfadjointView<Eigen::Lower>();      // full symmetric A from its lower triangle
        f.write((const char*)A.data(), sizeof(double) * D.K * D.K);
        f.write((const char*)g.B.data(), sizeof(double) * D.K * D.K);
        f.write((const char*)g.sL.data(), sizeof(double) * D.K);
        f.write((const char*)g.sL2.data(), sizeof(double) * D.K);
        f.write((const char*)g.sx, sizeof g.sx);
        f.write((const char*)g.so, sizeof g.so);
        f.write((const char*)g.acf, sizeof g.acf);
        f.write((const char*)g.acfn, sizeof g.acfn);
    }
    std::ofstream j(c.out + ".json");
    j << std::setprecision(12) << "{\n \"format\": \"NLS3GED1\", \"K\": " << D.K << ", \"groups\": " << G.size()
      << ", \"nlag\": " << NLAG << ", \"d\": " << c.d << ", \"M\": " << c.M
      << ", \"T1\": " << c.T1 << ", \"T3\": " << c.T3 << ", \"gamma\": " << c.gamma << ", \"dt\": " << c.dt
      << ", \"burn\": " << c.burn << ", \"measure\": " << c.measure << ", \"snapshot\": " << c.snapshot
      << ", \"batches\": " << c.batches << ", \"lanes\": " << c.batches * LANES << ", \"seed\": " << c.seed
      << ", \"from\": \"\", \"scales\": [" << c.s[0] << ", " << c.s[1] << ", " << c.s[2] << "]"
      << ", \"seconds\": " << seconds << ", \"midpoint_failures\": " << fails << ", \"nonfinite_batches\": " << nonfinite
      << ", \"min_action\": " << minact
      << ",\n \"observables\": [\"I2\", \"I1+I3\", \"I1-I3\", \"h1h2cos1\", \"h2h3cos3\", \"cos1\"]"
      << ",\n \"terms\": [";
    for (int k = 0; k < D.K; ++k) {
        const Term& t = D.terms[k];
        j << (k ? "," : "") << "[" << t.p[0] << "," << t.p[1] << "," << t.p[2] << "," << t.q[0] << "," << t.q[1] << ","
          << t.q[2] << "," << t.m1 << "," << t.m3 << "," << t.sine << "]";
    }
    j << "]\n}\n";
}

// ============================================================================ 8. simulation
int run(const Params& c) {
    Eigen::setNbThreads(1);                       // the matrix products run inside the OpenMP threads
    Dict D; D.terms = make_dictionary(c.d, c.M); D.K = (int)D.terms.size();
    D.dmax = c.d; D.qmax = 2 * c.M; for (int j = 0; j < 3; ++j) D.s[j] = c.s[j];
    std::cout << "K=" << D.K << " (d=" << c.d << ", M=" << c.M << ")" << std::endl;
    std::vector<Acc> G(c.groups);
    for (auto& g : G) g.init(D.K);
    const auto t0 = std::chrono::steady_clock::now();
    std::uint64_t fails = 0, nonfinite = 0; double minact = 1e300;

    const auto burn_steps = (std::int64_t)std::llround(c.burn / c.dt);       // steps of the burn-in
    const auto snap_steps = (std::int64_t)std::llround(c.snapshot / c.dt);   // steps between samples
    const auto nsnap = (std::int64_t)std::llround(c.measure / c.snapshot);   // samples per trajectory
    std::atomic<int> done{0};
#pragma omp parallel for num_threads(c.threads) schedule(dynamic, 1)
    for (int b = 0; b < c.batches; ++b) {
        // ---- 8a. one batch: its own sums, state c_j = I0^(1/2) (I0 = sqrt((T1+T3)/6)), generators
        Acc local; local.init(D.K);
        Feeder fd; fd.init(&D, c.T1, c.T3, c.gamma, &local);
        State st; MidpointWork ws;
        const double amplitude = std::sqrt(std::sqrt(0.5 * (c.T1 + c.T3) / N));
        for (int j = 0; j < N; ++j) { st.x[j].setConstant(amplitude); st.y[j].setZero(); }
        std::array<Rng, LANES> rng;
        for (int l = 0; l < LANES; ++l) rng[l].seed(c.seed, (std::uint64_t)(b * LANES + l));
        Vec16 lx, ly, rx, ry;                                      // noise of mode 1 (l) and mode 3 (r)
        std::vector<double> hist((size_t)LANES * NLAG * NOBS_ACF, 0.0);   // last NLAG values of the observables
        std::int64_t k = 0; std::uint64_t myfails = 0; bool bad = false; double mymin = 1e300;
        const std::int64_t total = burn_steps + nsnap * snap_steps;
        for (std::int64_t step = 1; step <= total; ++step) {
            // ---- 8b. one time step: noise (x of both baths, then y), midpoint step, two bath steps (alternating)
            gaussian_pair(rng, lx, rx);
            gaussian_pair(rng, ly, ry);
            int it = 0;
            if (!midpoint_step(st, ws, c.dt, it)) ++myfails;
            if (step % 2 == 0) { bath_step(st, 0, c.T1, c.gamma, lx, ly, c.dt); bath_step(st, N - 1, c.T3, c.gamma, rx, ry, c.dt); }
            else { bath_step(st, N - 1, c.T3, c.gamma, rx, ry, c.dt); bath_step(st, 0, c.T1, c.gamma, lx, ly, c.dt); }
            if (!finite_state(st)) { bad = true; break; }
            // ---- 8c. a sample: x of every trajectory, into the dictionary sums and the autocorrelation sums
            if (step > burn_steps && (step - burn_steps) % snap_steps == 0 && k < nsnap) {
                for (int l = 0; l < LANES; ++l) {
                    const double p1 = std::atan2(st.y[0](l), st.x[0](l));
                    const double p2 = std::atan2(st.y[1](l), st.x[1](l));
                    const double p3 = std::atan2(st.y[2](l), st.x[2](l));
                    double x[5] = {st.I[0](l), st.I[1](l), st.I[2](l), wrap_angle(2.0 * (p1 - p2)), wrap_angle(2.0 * (p3 - p2))};
                    for (int j = 0; j < 3; ++j) mymin = std::min(mymin, x[j]);
                    fd.push(x);
                    double o[NOBS_ACF]; observables(x, c.s, o);
                    double* H = &hist[(size_t)l * NLAG * NOBS_ACF];    // ring buffer of this trajectory
                    const int slot = (int)(k % NLAG);
                    for (int i = 0; i < NOBS_ACF; ++i) H[slot * NOBS_ACF + i] = o[i];
                    const int maxlag = (int)std::min<std::int64_t>(k, NLAG - 1);
                    for (int lag = 0; lag <= maxlag; ++lag) {
                        const int s2 = (int)((k - lag) % NLAG);
                        for (int i = 0; i < NOBS_ACF; ++i) local.acf[i][lag] += o[i] * H[s2 * NOBS_ACF + i];
                        local.acfn[lag] += 1;
                    }
                    for (int i = 0; i < NOBS_ACF; ++i) local.so[i] += o[i];
                }
                ++k;
            }
        }
        fd.flush();
        // ---- 8d. add the batch to its group (b mod groups)
#pragma omp critical
        {
            G[b % c.groups].add(local);
            fails += myfails; nonfinite += bad ? 1 : 0; minact = std::min(minact, mymin);
            std::cout << "batch " << b << " done (" << ++done << "/" << c.batches << ")" << (bad ? " NONFINITE" : "") << std::endl;
        }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    write_output(c, D, G, secs, fails, nonfinite, minact);
    double n = 0; for (auto& g : G) n += g.n;
    std::cout << "samples=" << n << " seconds=" << secs << " midpoint_failures=" << fails
              << " nonfinite_batches=" << nonfinite << " min_action=" << minact << std::endl;
    return nonfinite ? 2 : 0;
}

}  // namespace

// ============================================================================ 9. main
int main(int argc, char** argv) {
    try {
        const Params p = read_options(argc, argv);
        print_parameters(p);
        return p.dry_run ? 0 : run(p);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
