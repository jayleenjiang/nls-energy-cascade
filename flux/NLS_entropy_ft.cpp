// Joint bath-heat, entropy-production, and action-current sampler for the
// Gibbs-preserving boundary-driven NLS cascade chain.
//
// This program is intentionally separate from NLS_flux_canonical.cpp.  It
// represents c_j = x_j + i y_j in Cartesian coordinates, which is equivalent
// to the action-angle SDE away from I_j=0 but avoids the positivity projection
// and angular singularity of an Euler--Maruyama action-angle discretization.
//
// Physical energy convention:
//   E = H/2,
//   H = M^2 - (1/2) sum_j |c_j|^4
//       + 2 sum_j Re(c_j^2 conjugate(c_{j-1}^2)).
//
// With F_j = partial H / partial conjugate(c_j) = grad_{(x_j,y_j)} E,
// the dynamics is
//   dc_j = i F_j dt
// at interior sites and, at boundary r in {1,n},
//   dc_r = (i-gamma) F_r dt
//          + sqrt(2 gamma T_r) (dW_{r,x} + i dW_{r,y}).
//
// The Stratonovich heat delivered to E by bath r is
//   dQ_r = grad_r E o dc_r^bath.
// In Ito form, using Delta_r E = 4 M,
//   dQ_r = [-gamma |F_r|^2 + 4 gamma T_r M] dt
//          + sqrt(2 gamma T_r) F_r . dW_r.
// Positive Q means energy delivered by the bath to the system.  The medium
// entropy increment is dSigma_m = -dQ_1/T_1 - dQ_n/T_n.
//
// Numerically, each step is split into a time-reversible implicit-midpoint
// Hamiltonian update followed by separate left/right bath updates.  The bath
// order alternates every step.  Q_r is the exact change of E across the
// corresponding discrete bath substep, which converges to the Stratonovich
// heat above and makes the block energy-balance residual equal to the
// Hamiltonian integrator's accumulated energy error (up to roundoff).
//
// Build on Apple Silicon:
//   clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
//     -Xpreprocessor -fopenmp \
//     -I/opt/homebrew/include/eigen3 \
//     -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     NLS_entropy_ft.cpp -o entropy_ft
//
// Self-test:
//   ./entropy_ft selftest
//
// Sample:
//   ./entropy_ft sample T1 Tn n batches burnin block_time blocks_per_stream \
//       dt seed threads out_prefix [bond]
//
// One batch contains 16 statistically independent streams.  The optional
// bond is the right endpoint j of the zero-based pair (j-1,j) in the paper's
// one-based convention, so 1 <= bond < n; it defaults to n/2.

#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double GAMMA = 0.1;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int LANES = 16;
constexpr int N_BURNIN_MONITORS = 10;
constexpr int MAX_MIDPOINT_ITERATIONS = 20;
constexpr double MIDPOINT_TOLERANCE = 2.0e-13;
constexpr const char* MODEL_VERSION = "gibbs-cartesian-entropy-ft-v1";
constexpr const char* N3_ENDPOINT_MODEL_VERSION =
    "gibbs-cartesian-entropy-ft-n3-endpoints-v1";

using A16d = Eigen::Array<double, LANES, 1>;
using AlignedVec = std::vector<A16d, Eigen::aligned_allocator<A16d>>;

struct Xoshiro256pp {
    std::array<std::uint64_t, 4> s{};

    static std::uint64_t splitmix64(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    static std::uint64_t rotl(std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    void seed(std::uint64_t value) {
        for (auto& v : s) {
            v = splitmix64(value);
        }
    }

    std::uint64_t next() {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    double uniform_open() {
        constexpr double scale = 1.0 / 9007199254740992.0;
        return ((next() >> 11) + 0.5) * scale;
    }
};

std::uint64_t stream_seed(std::uint64_t base_seed,
                          std::uint64_t stream_id) {
    std::uint64_t x =
        base_seed ^ (0xD2B74407B1CE6E93ULL * (stream_id + 1ULL));
    return Xoshiro256pp::splitmix64(x);
}

void fill_gaussian_pair(std::array<Xoshiro256pp, LANES>& rng,
                        A16d& first,
                        A16d& second) {
    for (int lane = 0; lane < LANES; ++lane) {
        const double u1 = rng[lane].uniform_open();
        const double u2 = rng[lane].uniform_open();
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double angle = TWO_PI * u2;
        first(lane) = radius * std::cos(angle);
        second(lane) = radius * std::sin(angle);
    }
}

std::int64_t checked_step_count(double time,
                                double dt,
                                const std::string& label,
                                bool allow_zero = false) {
    const double raw = time / dt;
    const auto steps = static_cast<std::int64_t>(std::llround(raw));
    if ((!allow_zero && steps <= 0) || (allow_zero && steps < 0)) {
        throw std::invalid_argument(label + " has invalid step count");
    }
    const double reconstructed = steps * dt;
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(time));
    if (std::abs(reconstructed - time) > tolerance) {
        std::ostringstream message;
        message << label << "/dt must be an integer; got "
                << std::setprecision(17) << raw;
        throw std::invalid_argument(message.str());
    }
    return steps;
}

void ensure_parent_directory(const std::string& prefix) {
    const std::filesystem::path path(prefix);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

std::ofstream open_output(const std::string& path) {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error("cannot open output file: " + path);
    }
    stream << std::setprecision(17);
    return stream;
}

struct State {
    int n;
    AlignedVec x;
    AlignedVec y;
    AlignedVec action;
    AlignedVec square_real;
    AlignedVec square_imag;
    AlignedVec force_real;
    AlignedVec force_imag;
    A16d total_action = A16d::Zero();

    explicit State(int n_)
        : n(n_),
          x(n_),
          y(n_),
          action(n_),
          square_real(n_),
          square_imag(n_),
          force_real(n_),
          force_imag(n_) {}
};

struct MidpointWorkspace {
    State midpoint;
    AlignedVec old_x;
    AlignedVec old_y;
    AlignedVec guess_x;
    AlignedVec guess_y;
    AlignedVec candidate_x;
    AlignedVec candidate_y;

    explicit MidpointWorkspace(int n)
        : midpoint(n),
          old_x(n),
          old_y(n),
          guess_x(n),
          guess_y(n),
          candidate_x(n),
          candidate_y(n) {}
};

void initialize_state(State& state, double T1, double Tn) {
    const double initial_action =
        std::sqrt(0.5 * (T1 + Tn) / static_cast<double>(state.n));
    const double initial_amplitude = std::sqrt(initial_action);
    for (int j = 0; j < state.n; ++j) {
        state.x[j].setConstant(initial_amplitude);
        state.y[j].setZero();
    }
}

void compute_force(State& state) {
    state.total_action.setZero();
    for (int j = 0; j < state.n; ++j) {
        state.action[j] = state.x[j].square() + state.y[j].square();
        state.square_real[j] = state.x[j].square() - state.y[j].square();
        state.square_imag[j] = 2.0 * state.x[j] * state.y[j];
        state.total_action += state.action[j];
    }

    for (int j = 0; j < state.n; ++j) {
        A16d neighbor_real = A16d::Zero();
        A16d neighbor_imag = A16d::Zero();
        if (j > 0) {
            neighbor_real += state.square_real[j - 1];
            neighbor_imag += state.square_imag[j - 1];
        }
        if (j + 1 < state.n) {
            neighbor_real += state.square_real[j + 1];
            neighbor_imag += state.square_imag[j + 1];
        }

        const A16d onsite = 2.0 * state.total_action - state.action[j];
        state.force_real[j] =
            onsite * state.x[j] +
            2.0 * (neighbor_real * state.x[j] +
                   neighbor_imag * state.y[j]);
        state.force_imag[j] =
            onsite * state.y[j] +
            2.0 * (neighbor_imag * state.x[j] -
                   neighbor_real * state.y[j]);
    }
}

A16d physical_energy(const State& state) {
    A16d sum_action_squared = A16d::Zero();
    A16d interaction = A16d::Zero();
    for (int j = 0; j < state.n; ++j) {
        sum_action_squared += state.action[j].square();
        if (j > 0) {
            interaction +=
                state.square_real[j] * state.square_real[j - 1] +
                state.square_imag[j] * state.square_imag[j - 1];
        }
    }
    return 0.5 * state.total_action.square() -
           0.25 * sum_action_squared + interaction;
}

A16d action_current(const State& state, int bond) {
    return 4.0 *
           (state.square_imag[bond] * state.square_real[bond - 1] -
            state.square_real[bond] * state.square_imag[bond - 1]);
}

struct MidpointStepResult {
    A16d energy_error = A16d::Zero();
    int iterations = 0;
    bool converged = false;
};

MidpointStepResult hamiltonian_midpoint_step(State& state,
                                             MidpointWorkspace& workspace,
                                             double dt) {
    compute_force(state);
    const A16d energy_before = physical_energy(state);
    for (int j = 0; j < state.n; ++j) {
        workspace.old_x[j] = state.x[j];
        workspace.old_y[j] = state.y[j];
        workspace.guess_x[j] = state.x[j] - dt * state.force_imag[j];
        workspace.guess_y[j] = state.y[j] + dt * state.force_real[j];
    }

    MidpointStepResult result;
    for (int iteration = 1; iteration <= MAX_MIDPOINT_ITERATIONS;
         ++iteration) {
        for (int j = 0; j < state.n; ++j) {
            workspace.midpoint.x[j] =
                0.5 * (workspace.old_x[j] + workspace.guess_x[j]);
            workspace.midpoint.y[j] =
                0.5 * (workspace.old_y[j] + workspace.guess_y[j]);
        }
        compute_force(workspace.midpoint);

        double maximum_change = 0.0;
        for (int j = 0; j < state.n; ++j) {
            workspace.candidate_x[j] =
                workspace.old_x[j] -
                dt * workspace.midpoint.force_imag[j];
            workspace.candidate_y[j] =
                workspace.old_y[j] +
                dt * workspace.midpoint.force_real[j];
            maximum_change = std::max(
                maximum_change,
                (workspace.candidate_x[j] - workspace.guess_x[j])
                    .abs()
                    .maxCoeff());
            maximum_change = std::max(
                maximum_change,
                (workspace.candidate_y[j] - workspace.guess_y[j])
                    .abs()
                    .maxCoeff());
        }
        for (int j = 0; j < state.n; ++j) {
            workspace.guess_x[j] = workspace.candidate_x[j];
            workspace.guess_y[j] = workspace.candidate_y[j];
        }
        result.iterations = iteration;
        if (maximum_change < MIDPOINT_TOLERANCE) {
            result.converged = true;
            break;
        }
    }

    for (int j = 0; j < state.n; ++j) {
        state.x[j] = workspace.guess_x[j];
        state.y[j] = workspace.guess_y[j];
    }
    compute_force(state);
    result.energy_error = physical_energy(state) - energy_before;
    return result;
}

A16d bath_step(State& state,
                int site,
                double temperature,
                const A16d& normal_x,
                const A16d& normal_y,
                double dt,
                double sqrt_dt) {
    const A16d energy_before = physical_energy(state);
    const double noise_scale = std::sqrt(2.0 * GAMMA * temperature);
    state.x[site] +=
        -GAMMA * state.force_real[site] * dt +
        noise_scale * sqrt_dt * normal_x;
    state.y[site] +=
        -GAMMA * state.force_imag[site] * dt +
        noise_scale * sqrt_dt * normal_y;
    compute_force(state);
    return physical_energy(state) - energy_before;
}

struct Parameters {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int batches = 0;
    double burnin = 0.0;
    double block_time = 0.0;
    int blocks_per_stream = 0;
    double dt = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    int bond = 0;
    std::int64_t burn_steps = 0;
    std::int64_t block_steps = 0;
    bool save_n3_endpoints = false;
};

void validate_parameters(const Parameters& p) {
    if (!(p.T1 > 0.0) || !(p.Tn > 0.0)) {
        throw std::invalid_argument("T1 and Tn must be positive");
    }
    if (p.n < 2) {
        throw std::invalid_argument("n must be at least 2");
    }
    if (p.batches <= 0 || p.blocks_per_stream <= 0 || p.threads <= 0) {
        throw std::invalid_argument(
            "batches, blocks_per_stream, and threads must be positive");
    }
    if (!(p.dt > 0.0) || !(p.block_time > 0.0) || p.burnin < 0.0) {
        throw std::invalid_argument(
            "dt and block_time must be positive and burnin nonnegative");
    }
    if (p.bond < 1 || p.bond >= p.n) {
        throw std::invalid_argument("bond must satisfy 1 <= bond < n");
    }
}

Parameters parse_sample(int argc,
                        char* argv[],
                        bool save_n3_endpoints = false) {
    if (argc < 13 || argc > 14) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " sample T1 Tn n batches burnin block_time"
                << " blocks_per_stream dt seed threads out_prefix [bond]";
        throw std::invalid_argument(message.str());
    }
    Parameters p;
    p.T1 = std::stod(argv[2]);
    p.Tn = std::stod(argv[3]);
    p.n = std::stoi(argv[4]);
    p.batches = std::stoi(argv[5]);
    p.burnin = std::stod(argv[6]);
    p.block_time = std::stod(argv[7]);
    p.blocks_per_stream = std::stoi(argv[8]);
    p.dt = std::stod(argv[9]);
    p.seed = static_cast<std::uint64_t>(std::stoull(argv[10]));
    p.threads = std::stoi(argv[11]);
    p.prefix = argv[12];
    p.bond = argc == 13 ? p.n / 2 : std::stoi(argv[13]);
    p.save_n3_endpoints = save_n3_endpoints;
    validate_parameters(p);
    if (p.save_n3_endpoints && p.n != 3) {
        throw std::invalid_argument("sample_n3 requires n=3");
    }
    p.burn_steps = checked_step_count(p.burnin, p.dt, "burnin", true);
    p.block_steps = checked_step_count(p.block_time, p.dt, "block_time");
    return p;
}

struct BlockRecord {
    double q_left = 0.0;
    double q_right = 0.0;
    double delta_energy = 0.0;
    double entropy_medium = 0.0;
    double action_current = 0.0;
    std::array<double, 5> reduced_start{};
    std::array<double, 5> reduced_end{};
};

struct ReducedState3 {
    A16d i1 = A16d::Zero();
    A16d i2 = A16d::Zero();
    A16d i3 = A16d::Zero();
    A16d theta1 = A16d::Zero();
    A16d theta3 = A16d::Zero();
};

double wrap_angle(double value) {
    value = std::fmod(value + PI, TWO_PI);
    if (value < 0.0) {
        value += TWO_PI;
    }
    return value - PI;
}

ReducedState3 reduced_state3(const State& state) {
    if (state.n != 3) {
        throw std::invalid_argument("reduced_state3 requires n=3");
    }
    ReducedState3 reduced;
    reduced.i1 = state.action[0];
    reduced.i2 = state.action[1];
    reduced.i3 = state.action[2];
    for (int lane = 0; lane < LANES; ++lane) {
        const double phi1 = std::atan2(state.y[0](lane), state.x[0](lane));
        const double phi2 = std::atan2(state.y[1](lane), state.x[1](lane));
        const double phi3 = std::atan2(state.y[2](lane), state.x[2](lane));
        reduced.theta1(lane) = wrap_angle(2.0 * (phi1 - phi2));
        reduced.theta3(lane) = wrap_angle(2.0 * (phi3 - phi2));
    }
    return reduced;
}

struct BatchResult {
    std::vector<BlockRecord> blocks;
    std::vector<double> burnin_energy;
    AlignedVec action_integral;
    std::uint64_t midpoint_failure_count = 0;
    std::uint64_t midpoint_iteration_sum = 0;
    A16d hamiltonian_energy_error = A16d::Zero();
};

std::vector<std::int64_t> monitor_steps(std::int64_t burn_steps) {
    std::vector<std::int64_t> result;
    if (burn_steps <= 0) {
        return result;
    }
    result.reserve(N_BURNIN_MONITORS);
    for (int k = 1; k <= N_BURNIN_MONITORS; ++k) {
        const auto step = std::max<std::int64_t>(
            1, (burn_steps * k) / N_BURNIN_MONITORS);
        if (result.empty() || result.back() != step) {
            result.push_back(step);
        }
    }
    return result;
}

BatchResult run_batch(const Parameters& p, int batch_index) {
    BatchResult result;
    result.blocks.assign(
        static_cast<std::size_t>(p.blocks_per_stream) * LANES, {});
    const auto checkpoints = monitor_steps(p.burn_steps);
    result.burnin_energy.assign(checkpoints.size() * LANES, 0.0);
    result.action_integral.resize(p.n);
    for (auto& value : result.action_integral) {
        value.setZero();
    }

    State state(p.n);
    MidpointWorkspace midpoint_workspace(p.n);
    initialize_state(state, p.T1, p.Tn);

    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto id = static_cast<std::uint64_t>(batch_index) * LANES +
                        static_cast<std::uint64_t>(lane);
        rng[lane].seed(stream_seed(p.seed, id));
    }

    const double sqrt_dt = std::sqrt(p.dt);
    const auto measurement_steps =
        p.block_steps * static_cast<std::int64_t>(p.blocks_per_stream);
    const auto total_steps = p.burn_steps + measurement_steps;

    A16d q_left = A16d::Zero();
    A16d q_right = A16d::Zero();
    A16d action_integral = A16d::Zero();
    A16d block_start_energy = A16d::Zero();
    ReducedState3 block_start_reduced;
    std::size_t checkpoint_index = 0;
    int block_index = 0;

    for (std::int64_t step = 0; step < total_steps; ++step) {
        compute_force(state);
        const A16d energy_before = physical_energy(state);
        const A16d current = action_current(state, p.bond);

        if (step == p.burn_steps) {
            block_start_energy = energy_before;
            if (p.save_n3_endpoints) {
                block_start_reduced = reduced_state3(state);
            }
        }

        A16d normal_left_x, normal_right_x;
        A16d normal_left_y, normal_right_y;
        fill_gaussian_pair(rng, normal_left_x, normal_right_x);
        fill_gaussian_pair(rng, normal_left_y, normal_right_y);

        if (step >= p.burn_steps) {
            action_integral += current * p.dt;
            for (int j = 0; j < p.n; ++j) {
                result.action_integral[j] += state.action[j] * p.dt;
            }
        }

        const auto midpoint_result = hamiltonian_midpoint_step(
            state, midpoint_workspace, p.dt);
        result.midpoint_iteration_sum +=
            static_cast<std::uint64_t>(midpoint_result.iterations);
        if (!midpoint_result.converged) {
            ++result.midpoint_failure_count;
        }
        result.hamiltonian_energy_error += midpoint_result.energy_error;

        A16d step_q_left = A16d::Zero();
        A16d step_q_right = A16d::Zero();
        if (step % 2 == 0) {
            step_q_left = bath_step(
                state, 0, p.T1, normal_left_x, normal_left_y, p.dt, sqrt_dt);
            step_q_right = bath_step(
                state, p.n - 1, p.Tn, normal_right_x, normal_right_y,
                p.dt, sqrt_dt);
        } else {
            step_q_right = bath_step(
                state, p.n - 1, p.Tn, normal_right_x, normal_right_y,
                p.dt, sqrt_dt);
            step_q_left = bath_step(
                state, 0, p.T1, normal_left_x, normal_left_y, p.dt, sqrt_dt);
        }
        if (step >= p.burn_steps) {
            q_left += step_q_left;
            q_right += step_q_right;
        }

        if (checkpoint_index < checkpoints.size() &&
            step + 1 == checkpoints[checkpoint_index]) {
            compute_force(state);
            const A16d monitored_energy = physical_energy(state);
            for (int lane = 0; lane < LANES; ++lane) {
                result.burnin_energy[checkpoint_index * LANES + lane] =
                    monitored_energy(lane);
            }
            ++checkpoint_index;
        }

        if (step >= p.burn_steps) {
            const auto measurement_step = step - p.burn_steps + 1;
            if (measurement_step % p.block_steps == 0) {
                compute_force(state);
                const A16d block_end_energy = physical_energy(state);
                ReducedState3 block_end_reduced;
                if (p.save_n3_endpoints) {
                    block_end_reduced = reduced_state3(state);
                }
                for (int lane = 0; lane < LANES; ++lane) {
                    auto& record = result.blocks[
                        static_cast<std::size_t>(block_index) * LANES + lane];
                    record.q_left = q_left(lane);
                    record.q_right = q_right(lane);
                    record.delta_energy =
                        block_end_energy(lane) - block_start_energy(lane);
                    record.entropy_medium =
                        -q_left(lane) / p.T1 - q_right(lane) / p.Tn;
                    record.action_current =
                        action_integral(lane) / p.block_time;
                    if (p.save_n3_endpoints) {
                        record.reduced_start = {
                            block_start_reduced.i1(lane),
                            block_start_reduced.i2(lane),
                            block_start_reduced.i3(lane),
                            block_start_reduced.theta1(lane),
                            block_start_reduced.theta3(lane)};
                        record.reduced_end = {
                            block_end_reduced.i1(lane),
                            block_end_reduced.i2(lane),
                            block_end_reduced.i3(lane),
                            block_end_reduced.theta1(lane),
                            block_end_reduced.theta3(lane)};
                    }
                }
                q_left.setZero();
                q_right.setZero();
                action_integral.setZero();
                block_start_energy = block_end_energy;
                if (p.save_n3_endpoints) {
                    block_start_reduced = block_end_reduced;
                }
                ++block_index;
            }
        }
    }

    return result;
}

double mean_of(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

double rms_of(const std::vector<double>& values) {
    double total = 0.0;
    for (const double value : values) {
        total += value * value;
    }
    return std::sqrt(total / static_cast<double>(values.size()));
}

void run_sample(const Parameters& p) {
    ensure_parent_directory(p.prefix);
    const int stream_count = p.batches * LANES;
    const std::size_t total_blocks =
        static_cast<std::size_t>(stream_count) * p.blocks_per_stream;
    std::cout << "Joint entropy/action sampler\n"
              << "model="
              << (p.save_n3_endpoints ? N3_ENDPOINT_MODEL_VERSION
                                      : MODEL_VERSION)
              << " n=" << p.n
              << " T1=" << p.T1
              << " Tn=" << p.Tn
              << " dt=" << p.dt
              << " burnin=" << p.burnin
              << " block_time=" << p.block_time
              << " blocks_per_stream=" << p.blocks_per_stream
              << " streams=" << stream_count
              << " total_blocks=" << total_blocks
              << " bond=" << p.bond
              << " seed=" << p.seed
              << " threads=" << p.threads << "\n";

    const auto start = std::chrono::steady_clock::now();
    std::vector<BlockRecord> all_blocks(total_blocks);
    std::vector<double> all_burnin(
        static_cast<std::size_t>(stream_count) *
        monitor_steps(p.burn_steps).size(),
        0.0);
    AlignedVec action_integral(p.n);
    for (auto& value : action_integral) {
        value.setZero();
    }
    std::uint64_t midpoint_failure_count = 0;
    std::uint64_t midpoint_iteration_sum = 0;
    A16d hamiltonian_energy_error = A16d::Zero();
    std::atomic<int> completed{0};

    #pragma omp parallel for schedule(static) num_threads(p.threads)
    for (int batch = 0; batch < p.batches; ++batch) {
        auto result = run_batch(p, batch);
        for (int block = 0; block < p.blocks_per_stream; ++block) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int stream = batch * LANES + lane;
                const std::size_t target =
                    static_cast<std::size_t>(stream) * p.blocks_per_stream +
                    block;
                const std::size_t source =
                    static_cast<std::size_t>(block) * LANES + lane;
                all_blocks[target] = result.blocks[source];
            }
        }
        const auto checkpoints = monitor_steps(p.burn_steps);
        for (std::size_t checkpoint = 0; checkpoint < checkpoints.size();
             ++checkpoint) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int stream = batch * LANES + lane;
                all_burnin[static_cast<std::size_t>(stream) *
                               checkpoints.size() + checkpoint] =
                    result.burnin_energy[checkpoint * LANES + lane];
            }
        }
        #pragma omp critical
        {
            for (int j = 0; j < p.n; ++j) {
                action_integral[j] += result.action_integral[j];
            }
            midpoint_failure_count += result.midpoint_failure_count;
            midpoint_iteration_sum += result.midpoint_iteration_sum;
            hamiltonian_energy_error += result.hamiltonian_energy_error;
        }
        const int done = ++completed;
        if (done == p.batches || done % std::max(1, p.batches / 10) == 0) {
            #pragma omp critical
            std::cout << "completed " << done << "/" << p.batches
                      << " batches\n";
        }
    }

    const auto finish = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(finish - start).count();
    const auto steps_per_batch =
        p.burn_steps + p.block_steps *
                           static_cast<std::int64_t>(p.blocks_per_stream);
    const double midpoint_call_count =
        static_cast<double>(p.batches) * steps_per_batch;
    const double midpoint_failure_rate =
        static_cast<double>(midpoint_failure_count) / midpoint_call_count;
    const double mean_midpoint_iterations =
        static_cast<double>(midpoint_iteration_sum) / midpoint_call_count;
    const double total_simulated_time_per_stream =
        p.burnin + p.block_time * p.blocks_per_stream;
    const double mean_hamiltonian_energy_error_rate =
        hamiltonian_energy_error.sum() /
        (static_cast<double>(stream_count) * total_simulated_time_per_stream);

    std::vector<double> q_left_rate;
    std::vector<double> q_right_rate;
    std::vector<double> energy_drift_rate;
    std::vector<double> entropy_rate;
    std::vector<double> action_values;
    std::vector<double> balance_rate;
    q_left_rate.reserve(total_blocks);
    q_right_rate.reserve(total_blocks);
    energy_drift_rate.reserve(total_blocks);
    entropy_rate.reserve(total_blocks);
    action_values.reserve(total_blocks);
    balance_rate.reserve(total_blocks);
    for (const auto& record : all_blocks) {
        q_left_rate.push_back(record.q_left / p.block_time);
        q_right_rate.push_back(record.q_right / p.block_time);
        energy_drift_rate.push_back(record.delta_energy / p.block_time);
        entropy_rate.push_back(record.entropy_medium / p.block_time);
        action_values.push_back(record.action_current);
        balance_rate.push_back(
            (record.q_left + record.q_right - record.delta_energy) /
            p.block_time);
    }

    {
        auto output = open_output(p.prefix + "_blocks.csv");
        output << "stream_id,block_id,q_left,q_right,delta_energy,"
               << "entropy_medium,entropy_rate,action_current,"
               << "energy_balance_error";
        if (p.save_n3_endpoints) {
            output << ",start_I1,start_I2,start_I3,start_theta1,start_theta3,"
                   << "end_I1,end_I2,end_I3,end_theta1,end_theta3";
        }
        output << '\n';
        for (int stream = 0; stream < stream_count; ++stream) {
            for (int block = 0; block < p.blocks_per_stream; ++block) {
                const auto& record = all_blocks[
                    static_cast<std::size_t>(stream) * p.blocks_per_stream +
                    block];
                output << stream << ',' << block << ','
                       << record.q_left << ',' << record.q_right << ','
                       << record.delta_energy << ','
                       << record.entropy_medium << ','
                       << record.entropy_medium / p.block_time << ','
                       << record.action_current << ','
                       << record.q_left + record.q_right -
                              record.delta_energy;
                if (p.save_n3_endpoints) {
                    for (const double value : record.reduced_start) {
                        output << ',' << value;
                    }
                    for (const double value : record.reduced_end) {
                        output << ',' << value;
                    }
                }
                output << '\n';
            }
        }
    }

    {
        auto output = open_output(p.prefix + "_profile.csv");
        output << "site,mean_action\n";
        const double denominator =
            static_cast<double>(stream_count) * p.blocks_per_stream *
            p.block_time;
        for (int j = 0; j < p.n; ++j) {
            output << j + 1 << ','
                   << action_integral[j].sum() / denominator << '\n';
        }
    }

    {
        auto output = open_output(p.prefix + "_burnin.csv");
        output << "stream_id,time,energy\n";
        const auto checkpoints = monitor_steps(p.burn_steps);
        for (int stream = 0; stream < stream_count; ++stream) {
            for (std::size_t checkpoint = 0; checkpoint < checkpoints.size();
                 ++checkpoint) {
                output << stream << ','
                       << checkpoints[checkpoint] * p.dt << ','
                       << all_burnin[static_cast<std::size_t>(stream) *
                                         checkpoints.size() + checkpoint]
                       << '\n';
            }
        }
    }

    {
        auto output = open_output(p.prefix + "_summary.csv");
        output
            << "model_version,coordinates,energy_convention,n,T1,Tn,gamma,dt,"
            << "burnin,block_time,blocks_per_stream,batches,lanes,n_streams,"
            << "n_blocks,bond,seed,threads,mean_q_left_rate,mean_q_right_rate,"
            << "mean_energy_drift_rate,mean_entropy_rate,mean_action_current,"
            << "mean_energy_balance_error_rate,rms_energy_balance_error_rate,"
            << "midpoint_failure_count,midpoint_failure_rate,"
            << "mean_midpoint_iterations,mean_hamiltonian_energy_error_rate,"
            << "elapsed_seconds\n";
        output << (p.save_n3_endpoints ? N3_ENDPOINT_MODEL_VERSION
                                      : MODEL_VERSION)
               << (p.save_n3_endpoints ? ",cartesian+reduced-endpoints,E=H/2,"
                                       : ",cartesian,E=H/2,")
               << p.n << ',' << p.T1 << ',' << p.Tn << ',' << GAMMA << ','
               << p.dt << ',' << p.burnin << ',' << p.block_time << ','
               << p.blocks_per_stream << ',' << p.batches << ',' << LANES
               << ',' << stream_count << ',' << total_blocks << ',' << p.bond
               << ',' << p.seed << ',' << p.threads << ','
               << mean_of(q_left_rate) << ',' << mean_of(q_right_rate) << ','
               << mean_of(energy_drift_rate) << ','
               << mean_of(entropy_rate) << ',' << mean_of(action_values) << ','
               << mean_of(balance_rate) << ',' << rms_of(balance_rate) << ','
               << midpoint_failure_count << ',' << midpoint_failure_rate << ','
               << mean_midpoint_iterations << ','
               << mean_hamiltonian_energy_error_rate << ','
               << elapsed << '\n';
    }

    std::cout << "mean q_left rate=" << mean_of(q_left_rate)
              << " mean q_right rate=" << mean_of(q_right_rate)
              << " mean entropy rate=" << mean_of(entropy_rate)
              << " mean action current=" << mean_of(action_values)
              << " balance RMS rate=" << rms_of(balance_rate)
              << " midpoint failures=" << midpoint_failure_count
              << " mean midpoint iterations=" << mean_midpoint_iterations
              << " elapsed=" << elapsed << " s\n";
}

struct ScalarState {
    std::vector<double> x;
    std::vector<double> y;
};

double scalar_energy(const ScalarState& state) {
    const int n = static_cast<int>(state.x.size());
    std::vector<double> action(n, 0.0);
    std::vector<double> square_real(n, 0.0);
    std::vector<double> square_imag(n, 0.0);
    double total_action = 0.0;
    double action_squared = 0.0;
    for (int j = 0; j < n; ++j) {
        action[j] = state.x[j] * state.x[j] + state.y[j] * state.y[j];
        square_real[j] = state.x[j] * state.x[j] -
                         state.y[j] * state.y[j];
        square_imag[j] = 2.0 * state.x[j] * state.y[j];
        total_action += action[j];
        action_squared += action[j] * action[j];
    }
    double interaction = 0.0;
    for (int j = 1; j < n; ++j) {
        interaction += square_real[j] * square_real[j - 1] +
                       square_imag[j] * square_imag[j - 1];
    }
    return 0.5 * total_action * total_action - 0.25 * action_squared +
           interaction;
}

void scalar_force(const ScalarState& state,
                  std::vector<double>& force_real,
                  std::vector<double>& force_imag) {
    const int n = static_cast<int>(state.x.size());
    std::vector<double> action(n, 0.0);
    std::vector<double> square_real(n, 0.0);
    std::vector<double> square_imag(n, 0.0);
    double total_action = 0.0;
    for (int j = 0; j < n; ++j) {
        action[j] = state.x[j] * state.x[j] + state.y[j] * state.y[j];
        square_real[j] = state.x[j] * state.x[j] -
                         state.y[j] * state.y[j];
        square_imag[j] = 2.0 * state.x[j] * state.y[j];
        total_action += action[j];
    }
    force_real.assign(n, 0.0);
    force_imag.assign(n, 0.0);
    for (int j = 0; j < n; ++j) {
        double neighbor_real = 0.0;
        double neighbor_imag = 0.0;
        if (j > 0) {
            neighbor_real += square_real[j - 1];
            neighbor_imag += square_imag[j - 1];
        }
        if (j + 1 < n) {
            neighbor_real += square_real[j + 1];
            neighbor_imag += square_imag[j + 1];
        }
        const double onsite = 2.0 * total_action - action[j];
        force_real[j] =
            onsite * state.x[j] +
            2.0 * (neighbor_real * state.x[j] +
                   neighbor_imag * state.y[j]);
        force_imag[j] =
            onsite * state.y[j] +
            2.0 * (neighbor_imag * state.x[j] -
                   neighbor_real * state.y[j]);
    }
}

int run_selftest() {
    constexpr int n = 5;
    ScalarState state;
    state.x.resize(n);
    state.y.resize(n);
    std::mt19937_64 generator(20260826ULL);
    std::uniform_real_distribution<double> distribution(-0.9, 0.9);
    for (int j = 0; j < n; ++j) {
        state.x[j] = distribution(generator);
        state.y[j] = distribution(generator);
    }

    std::vector<double> force_real;
    std::vector<double> force_imag;
    scalar_force(state, force_real, force_imag);
    const double epsilon = 1.0e-6;
    double maximum_gradient_error = 0.0;
    for (int j = 0; j < n; ++j) {
        ScalarState plus = state;
        ScalarState minus = state;
        plus.x[j] += epsilon;
        minus.x[j] -= epsilon;
        const double finite_x =
            (scalar_energy(plus) - scalar_energy(minus)) / (2.0 * epsilon);
        plus = state;
        minus = state;
        plus.y[j] += epsilon;
        minus.y[j] -= epsilon;
        const double finite_y =
            (scalar_energy(plus) - scalar_energy(minus)) / (2.0 * epsilon);
        maximum_gradient_error = std::max(
            maximum_gradient_error,
            std::max(std::abs(finite_x - force_real[j]),
                     std::abs(finite_y - force_imag[j])));
    }

    double hamiltonian_energy_derivative = 0.0;
    for (int j = 0; j < n; ++j) {
        hamiltonian_energy_derivative +=
            force_real[j] * (-force_imag[j]) +
            force_imag[j] * force_real[j];
    }

    double maximum_laplacian_error = 0.0;
    double total_action = 0.0;
    for (int j = 0; j < n; ++j) {
        total_action += state.x[j] * state.x[j] +
                        state.y[j] * state.y[j];
    }
    for (const int j : {0, n - 1}) {
        ScalarState plus = state;
        ScalarState minus = state;
        const double center = scalar_energy(state);
        plus.x[j] += epsilon;
        minus.x[j] -= epsilon;
        const double second_x =
            (scalar_energy(plus) - 2.0 * center + scalar_energy(minus)) /
            (epsilon * epsilon);
        plus = state;
        minus = state;
        plus.y[j] += epsilon;
        minus.y[j] -= epsilon;
        const double second_y =
            (scalar_energy(plus) - 2.0 * center + scalar_energy(minus)) /
            (epsilon * epsilon);
        maximum_laplacian_error = std::max(
            maximum_laplacian_error,
            std::abs(second_x + second_y - 4.0 * total_action));
    }

    std::cout << std::setprecision(12)
              << "selftest maximum gradient error="
              << maximum_gradient_error << '\n'
              << "selftest Hamiltonian dE/dt="
              << hamiltonian_energy_derivative << '\n'
              << "selftest maximum boundary Laplacian error="
              << maximum_laplacian_error << '\n';

    const bool passed =
        maximum_gradient_error < 2.0e-8 &&
        std::abs(hamiltonian_energy_derivative) < 1.0e-12 &&
        maximum_laplacian_error < 3.0e-3;
    std::cout << (passed ? "SELFTEST PASS\n" : "SELFTEST FAIL\n");
    return passed ? 0 : 1;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "selftest") {
            return run_selftest();
        }
        if (argc >= 2 && std::string(argv[1]) == "sample") {
            run_sample(parse_sample(argc, argv));
            return 0;
        }
        if (argc >= 2 && std::string(argv[1]) == "sample_n3") {
            run_sample(parse_sample(argc, argv, true));
            return 0;
        }
        std::cerr
            << "Usage:\n"
            << "  " << argv[0] << " selftest\n"
            << "  " << argv[0]
            << " sample T1 Tn n batches burnin block_time blocks_per_stream"
            << " dt seed threads out_prefix [bond]\n"
            << "  " << argv[0]
            << " sample_n3 T1 Tn 3 batches burnin block_time blocks_per_stream"
            << " dt seed threads out_prefix [bond]\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
