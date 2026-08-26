// Transient relaxation and finite-window action-current probes for the
// Gibbs-preserving boundary-driven NLS chain.
//
// This file is intentionally separate from NLS_flux_canonical.cpp so that the
// existing production/current-scaling workflow remains untouched.
//
// Build on Apple Silicon:
//   clang++ -O3 -mcpu=native -std=c++17 \
//     -Xpreprocessor -fopenmp \
//     -I/opt/homebrew/include/eigen3 \
//     -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     NLS_flux_relaxation_tau.cpp -o flux_relax_tau
//
// Modes:
//   transient T1 Tn n batches total_time dt checkpoint_dt seed threads out_prefix [bond]
//   tau       T1 Tn n batches burnin measure dt tau_block seed threads out_prefix [bond]
//
// The measured current is the conserved action current
//   J_j = 4 I_{j-1} I_j sin(2(phi_j - phi_{j-1})).

#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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

namespace {

constexpr double GAMMA = 0.1;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr double ACTION_FLOOR = 1.0e-12;
constexpr int LANES = 16;
constexpr const char* MODEL_VERSION = "gibbs-canonical-relax-tau-v1";

using A16d = Eigen::Array<double, LANES, 1>;
using AlignedVec =
    std::vector<A16d, Eigen::aligned_allocator<A16d>>;

A16d wrap_pi(const A16d& x) {
    return x - (x / TWO_PI).round() * TWO_PI;
}

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

std::uint64_t trajectory_seed(std::uint64_t base_seed,
                              std::uint64_t trajectory_id) {
    std::uint64_t x =
        base_seed ^ (0xD2B74407B1CE6E93ULL * (trajectory_id + 1ULL));
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

std::int64_t checked_step_count(double time, double dt,
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

struct ChainState {
    int n;
    AlignedVec action;
    AlignedVec phase;
    AlignedVec drift_action;
    AlignedVec drift_phase;
    AlignedVec sin_bond;
    AlignedVec cos_bond;

    explicit ChainState(int n_)
        : n(n_),
          action(n_),
          phase(n_),
          drift_action(n_),
          drift_phase(n_),
          sin_bond(n_),
          cos_bond(n_) {}
};

void initialize_state(ChainState& state, double T1, double Tn) {
    const double initial_action =
        std::sqrt(0.5 * (T1 + Tn) / static_cast<double>(state.n));
    for (int j = 0; j < state.n; ++j) {
        state.action[j].setConstant(initial_action);
        state.phase[j].setZero();
    }
}

A16d compute_drift_and_current(ChainState& state, double T1, double Tn,
                               int bond) {
    A16d total_action = A16d::Zero();
    for (int j = 0; j < state.n; ++j) {
        total_action += state.action[j];
    }

    for (int j = 1; j < state.n; ++j) {
        const A16d difference =
            wrap_pi(2.0 * (state.phase[j] - state.phase[j - 1]));
        state.sin_bond[j] = difference.sin();
        state.cos_bond[j] = difference.cos();
    }

    for (int j = 0; j < state.n; ++j) {
        A16d action_left = A16d::Zero();
        A16d action_right = A16d::Zero();
        A16d phase_left = A16d::Zero();
        A16d phase_right = A16d::Zero();

        if (j > 0) {
            action_left = state.action[j - 1] * state.sin_bond[j];
            phase_left = 2.0 * state.action[j - 1] * state.cos_bond[j];
        }
        if (j + 1 < state.n) {
            action_right = -state.action[j + 1] * state.sin_bond[j + 1];
            phase_right = 2.0 * state.action[j + 1] * state.cos_bond[j + 1];
        }

        state.drift_action[j] =
            4.0 * state.action[j] * (action_left + action_right);
        state.drift_phase[j] =
            2.0 * total_action - state.action[j] + phase_left + phase_right;
    }

    const A16d dleft = wrap_pi(2.0 * (state.phase[0] - state.phase[1]));
    state.drift_action[0] += 2.0 * GAMMA *
        (2.0 * T1 -
         (2.0 * total_action * state.action[0] - state.action[0].square() +
          2.0 * state.action[1] * state.action[0] * dleft.cos()));
    state.drift_phase[0] += 2.0 * GAMMA * state.action[1] * dleft.sin();

    const A16d dright =
        wrap_pi(2.0 * (state.phase[state.n - 1] -
                       state.phase[state.n - 2]));
    state.drift_action[state.n - 1] += 2.0 * GAMMA *
        (2.0 * Tn -
         (2.0 * total_action * state.action[state.n - 1] -
          state.action[state.n - 1].square() +
          2.0 * state.action[state.n - 2] *
              state.action[state.n - 1] * dright.cos()));
    state.drift_phase[state.n - 1] +=
        2.0 * GAMMA * state.action[state.n - 2] * dright.sin();

    return 4.0 * state.action[bond - 1] * state.action[bond] *
           state.sin_bond[bond];
}

std::uint64_t euler_maruyama_step(ChainState& state, double T1, double Tn,
                                  double dt, double sqrt_dt,
                                  std::array<Xoshiro256pp, LANES>& rng) {
    A16d normal_action_left, normal_action_right;
    A16d normal_phase_left, normal_phase_right;
    fill_gaussian_pair(rng, normal_action_left, normal_action_right);
    fill_gaussian_pair(rng, normal_phase_left, normal_phase_right);

    const A16d action_left_start = state.action[0].max(ACTION_FLOOR);
    const A16d action_right_start =
        state.action[state.n - 1].max(ACTION_FLOOR);

    for (int j = 0; j < state.n; ++j) {
        state.action[j] += state.drift_action[j] * dt;
        state.phase[j] = wrap_pi(state.phase[j] + state.drift_phase[j] * dt);
    }

    state.action[0] +=
        2.0 * (2.0 * GAMMA * T1 * action_left_start).sqrt() *
        sqrt_dt * normal_action_left;
    state.action[state.n - 1] +=
        2.0 * (2.0 * GAMMA * Tn * action_right_start).sqrt() *
        sqrt_dt * normal_action_right;
    state.phase[0] = wrap_pi(
        state.phase[0] +
        (2.0 * GAMMA * T1 / action_left_start).sqrt() *
            sqrt_dt * normal_phase_left);
    state.phase[state.n - 1] = wrap_pi(
        state.phase[state.n - 1] +
        (2.0 * GAMMA * Tn / action_right_start).sqrt() *
            sqrt_dt * normal_phase_right);

    std::uint64_t projection_count = 0;
    for (int j = 0; j < state.n; ++j) {
        projection_count +=
            static_cast<std::uint64_t>((state.action[j] < ACTION_FLOOR).count());
        state.action[j] = state.action[j].max(ACTION_FLOOR);
    }
    return projection_count;
}

struct CommonParams {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int batches = 0;
    double dt = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    int bond = 0;
};

void validate_common(const CommonParams& p) {
    if (!(p.T1 > 0.0) || !(p.Tn > 0.0)) {
        throw std::invalid_argument("T1 and Tn must be positive");
    }
    if (p.n < 2) {
        throw std::invalid_argument("n must be at least 2");
    }
    if (p.batches <= 0 || p.threads <= 0) {
        throw std::invalid_argument("batches and threads must be positive");
    }
    if (!(p.dt > 0.0)) {
        throw std::invalid_argument("dt must be positive");
    }
    if (p.bond < 1 || p.bond >= p.n) {
        throw std::invalid_argument("bond must satisfy 1 <= bond < n");
    }
}

struct TransientParams {
    CommonParams common;
    double total_time = 0.0;
    double checkpoint_dt = 0.0;
    std::int64_t total_steps = 0;
    std::int64_t checkpoint_steps = 0;
    int checkpoint_count = 0;
};

struct TransientBatchResult {
    std::vector<double> cumulative_current_samples;
    std::vector<double> interval_current_samples;
    std::vector<double> terminal_action_sum;
    std::vector<double> cumulative_action_sum;
    std::uint64_t projection_count = 0;
};

TransientBatchResult run_transient_batch(const TransientParams& p,
                                         int batch_index) {
    const int n = p.common.n;
    const int checkpoints = p.checkpoint_count;
    TransientBatchResult result;
    result.cumulative_current_samples.assign(checkpoints * LANES, 0.0);
    result.interval_current_samples.assign(checkpoints * LANES, 0.0);
    result.terminal_action_sum.assign(checkpoints * n, 0.0);
    result.cumulative_action_sum.assign(checkpoints * n, 0.0);

    ChainState state(n);
    initialize_state(state, p.common.T1, p.common.Tn);

    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto trajectory_id =
            static_cast<std::uint64_t>(batch_index) * LANES + lane;
        rng[lane].seed(trajectory_seed(p.common.seed, trajectory_id));
    }

    AlignedVec cumulative_action_integral(n);
    for (auto& value : cumulative_action_integral) {
        value.setZero();
    }
    A16d cumulative_current_integral = A16d::Zero();
    A16d interval_current_integral = A16d::Zero();
    const double sqrt_dt = std::sqrt(p.common.dt);

    int checkpoint_index = 0;
    for (std::int64_t step = 0; step < p.total_steps; ++step) {
        const A16d current = compute_drift_and_current(
            state, p.common.T1, p.common.Tn, p.common.bond);

        cumulative_current_integral += current * p.common.dt;
        interval_current_integral += current * p.common.dt;
        for (int j = 0; j < n; ++j) {
            cumulative_action_integral[j] += state.action[j] * p.common.dt;
        }

        result.projection_count += euler_maruyama_step(
            state, p.common.T1, p.common.Tn, p.common.dt, sqrt_dt, rng);

        if ((step + 1) % p.checkpoint_steps == 0) {
            const double time = (step + 1) * p.common.dt;
            for (int lane = 0; lane < LANES; ++lane) {
                const int offset = checkpoint_index * LANES + lane;
                result.cumulative_current_samples[offset] =
                    cumulative_current_integral(lane) / time;
                result.interval_current_samples[offset] =
                    interval_current_integral(lane) / p.checkpoint_dt;
            }
            for (int j = 0; j < n; ++j) {
                result.terminal_action_sum[checkpoint_index * n + j] =
                    state.action[j].sum();
                result.cumulative_action_sum[checkpoint_index * n + j] =
                    cumulative_action_integral[j].sum() / time;
            }
            interval_current_integral.setZero();
            ++checkpoint_index;
        }
    }

    return result;
}

struct TauParams {
    CommonParams common;
    double burnin = 0.0;
    double measure = 0.0;
    double tau_block = 0.0;
    std::int64_t burn_steps = 0;
    std::int64_t measure_steps = 0;
    std::int64_t tau_block_steps = 0;
    int block_count = 0;
};

struct TauBatchResult {
    std::vector<double> block_current_samples;
    std::uint64_t projection_count = 0;
};

TauBatchResult run_tau_batch(const TauParams& p, int batch_index) {
    const int n = p.common.n;
    TauBatchResult result;
    result.block_current_samples.assign(p.block_count * LANES, 0.0);

    ChainState state(n);
    initialize_state(state, p.common.T1, p.common.Tn);

    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto trajectory_id =
            static_cast<std::uint64_t>(batch_index) * LANES + lane;
        rng[lane].seed(trajectory_seed(p.common.seed, trajectory_id));
    }

    const double sqrt_dt = std::sqrt(p.common.dt);
    const std::int64_t total_steps = p.burn_steps + p.measure_steps;
    int block_index = 0;
    A16d block_current_integral = A16d::Zero();

    for (std::int64_t step = 0; step < total_steps; ++step) {
        const A16d current = compute_drift_and_current(
            state, p.common.T1, p.common.Tn, p.common.bond);

        if (step >= p.burn_steps) {
            block_current_integral += current * p.common.dt;
            const auto measurement_step = step - p.burn_steps + 1;
            if (measurement_step % p.tau_block_steps == 0) {
                for (int lane = 0; lane < LANES; ++lane) {
                    result.block_current_samples[block_index * LANES + lane] =
                        block_current_integral(lane) / p.tau_block;
                }
                block_current_integral.setZero();
                ++block_index;
            }
        }

        result.projection_count += euler_maruyama_step(
            state, p.common.T1, p.common.Tn, p.common.dt, sqrt_dt, rng);
    }

    return result;
}

double mean_of(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

double sample_sd_of(const std::vector<double>& values, double mean) {
    double sum = 0.0;
    for (const double value : values) {
        const double diff = value - mean;
        sum += diff * diff;
    }
    return std::sqrt(sum / static_cast<double>(values.size() - 1));
}

TransientParams parse_transient(int argc, char* argv[]) {
    if (argc < 12 || argc > 13) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " transient T1 Tn n batches total_time dt checkpoint_dt"
                << " seed threads out_prefix [bond]";
        throw std::invalid_argument(message.str());
    }
    TransientParams p;
    p.common.T1 = std::stod(argv[2]);
    p.common.Tn = std::stod(argv[3]);
    p.common.n = std::stoi(argv[4]);
    p.common.batches = std::stoi(argv[5]);
    p.total_time = std::stod(argv[6]);
    p.common.dt = std::stod(argv[7]);
    p.checkpoint_dt = std::stod(argv[8]);
    p.common.seed = static_cast<std::uint64_t>(std::stoull(argv[9]));
    p.common.threads = std::stoi(argv[10]);
    p.common.prefix = argv[11];
    p.common.bond = argc == 12 ? p.common.n / 2 : std::stoi(argv[12]);
    validate_common(p.common);
    p.total_steps = checked_step_count(p.total_time, p.common.dt, "total_time");
    p.checkpoint_steps =
        checked_step_count(p.checkpoint_dt, p.common.dt, "checkpoint_dt");
    if (p.total_steps % p.checkpoint_steps != 0) {
        throw std::invalid_argument("total_time must be divisible by checkpoint_dt");
    }
    p.checkpoint_count = static_cast<int>(p.total_steps / p.checkpoint_steps);
    return p;
}

TauParams parse_tau(int argc, char* argv[]) {
    if (argc < 13 || argc > 14) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " tau T1 Tn n batches burnin measure dt tau_block"
                << " seed threads out_prefix [bond]";
        throw std::invalid_argument(message.str());
    }
    TauParams p;
    p.common.T1 = std::stod(argv[2]);
    p.common.Tn = std::stod(argv[3]);
    p.common.n = std::stoi(argv[4]);
    p.common.batches = std::stoi(argv[5]);
    p.burnin = std::stod(argv[6]);
    p.measure = std::stod(argv[7]);
    p.common.dt = std::stod(argv[8]);
    p.tau_block = std::stod(argv[9]);
    p.common.seed = static_cast<std::uint64_t>(std::stoull(argv[10]));
    p.common.threads = std::stoi(argv[11]);
    p.common.prefix = argv[12];
    p.common.bond = argc == 13 ? p.common.n / 2 : std::stoi(argv[13]);
    validate_common(p.common);
    p.burn_steps = checked_step_count(p.burnin, p.common.dt, "burnin", true);
    p.measure_steps = checked_step_count(p.measure, p.common.dt, "measure");
    p.tau_block_steps = checked_step_count(p.tau_block, p.common.dt, "tau_block");
    if (p.measure_steps % p.tau_block_steps != 0) {
        throw std::invalid_argument("measure must be divisible by tau_block");
    }
    p.block_count = static_cast<int>(p.measure_steps / p.tau_block_steps);
    return p;
}

void run_transient(const TransientParams& p) {
    ensure_parent_directory(p.common.prefix);
    const int trajectories = p.common.batches * LANES;
    std::cout << "Transient relaxation probe\n"
              << "model=" << MODEL_VERSION
              << " n=" << p.common.n
              << " T1=" << p.common.T1
              << " Tn=" << p.common.Tn
              << " bond=" << p.common.bond
              << " dt=" << p.common.dt
              << " total_time=" << p.total_time
              << " checkpoint_dt=" << p.checkpoint_dt
              << " trajectories=" << trajectories
              << " seed=" << p.common.seed
              << " threads=" << p.common.threads << "\n";

    const auto start = std::chrono::steady_clock::now();
    std::vector<double> cumulative_samples(
        static_cast<std::size_t>(p.checkpoint_count) * trajectories, 0.0);
    std::vector<double> interval_samples(
        static_cast<std::size_t>(p.checkpoint_count) * trajectories, 0.0);
    std::vector<double> terminal_action_sum(
        static_cast<std::size_t>(p.checkpoint_count) * p.common.n, 0.0);
    std::vector<double> cumulative_action_sum(
        static_cast<std::size_t>(p.checkpoint_count) * p.common.n, 0.0);
    std::vector<std::uint64_t> projections(p.common.batches, 0);
    std::atomic<int> completed{0};

    #pragma omp parallel for schedule(static) num_threads(p.common.threads)
    for (int batch = 0; batch < p.common.batches; ++batch) {
        auto result = run_transient_batch(p, batch);
        for (int c = 0; c < p.checkpoint_count; ++c) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int trajectory = batch * LANES + lane;
                cumulative_samples[
                    static_cast<std::size_t>(c) * trajectories + trajectory] =
                    result.cumulative_current_samples[c * LANES + lane];
                interval_samples[
                    static_cast<std::size_t>(c) * trajectories + trajectory] =
                    result.interval_current_samples[c * LANES + lane];
            }
            for (int j = 0; j < p.common.n; ++j) {
                #pragma omp atomic
                terminal_action_sum[
                    static_cast<std::size_t>(c) * p.common.n + j] +=
                    result.terminal_action_sum[c * p.common.n + j];
                #pragma omp atomic
                cumulative_action_sum[
                    static_cast<std::size_t>(c) * p.common.n + j] +=
                    result.cumulative_action_sum[c * p.common.n + j];
            }
        }
        projections[batch] = result.projection_count;
        const int done = ++completed;
        const int stride = std::max(1, p.common.batches / 10);
        if (done == p.common.batches || done % stride == 0) {
            #pragma omp critical
            std::cout << "completed " << done << "/" << p.common.batches
                      << " batches\n";
        }
    }

    const auto stop = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(stop - start).count();

    {
        auto output = open_output(p.common.prefix + "_flux_timeseries.csv");
        output << "n,time,trajectories,mean_cumulative_current,"
               << "sd_cumulative_current,se_cumulative_current,"
               << "mean_last_interval_current,sd_last_interval_current,"
               << "se_last_interval_current\n";
        for (int c = 0; c < p.checkpoint_count; ++c) {
            const double time = (c + 1) * p.checkpoint_dt;
            std::vector<double> cumulative(trajectories), interval(trajectories);
            for (int i = 0; i < trajectories; ++i) {
                cumulative[i] = cumulative_samples[
                    static_cast<std::size_t>(c) * trajectories + i];
                interval[i] = interval_samples[
                    static_cast<std::size_t>(c) * trajectories + i];
            }
            const double cm = mean_of(cumulative);
            const double csd = sample_sd_of(cumulative, cm);
            const double im = mean_of(interval);
            const double isd = sample_sd_of(interval, im);
            output << p.common.n << "," << time << "," << trajectories << ","
                   << cm << "," << csd << ","
                   << csd / std::sqrt(static_cast<double>(trajectories)) << ","
                   << im << "," << isd << ","
                   << isd / std::sqrt(static_cast<double>(trajectories)) << "\n";
        }
    }
    {
        auto output = open_output(p.common.prefix + "_profile_timeseries.csv");
        output << "n,time,mode,mean_terminal_action,mean_cumulative_action\n";
        for (int c = 0; c < p.checkpoint_count; ++c) {
            const double time = (c + 1) * p.checkpoint_dt;
            for (int j = 0; j < p.common.n; ++j) {
                const auto index =
                    static_cast<std::size_t>(c) * p.common.n + j;
                output << p.common.n << "," << time << "," << j << ","
                       << terminal_action_sum[index] /
                              static_cast<double>(trajectories)
                       << ","
                       << cumulative_action_sum[index] /
                              static_cast<double>(trajectories)
                       << "\n";
            }
        }
    }
    {
        const std::uint64_t projection_count =
            std::accumulate(projections.begin(), projections.end(),
                            std::uint64_t{0});
        auto output = open_output(p.common.prefix + "_summary.csv");
        output << "model_version,mode,n,T1,Tn,gamma,dt,total_time,"
               << "checkpoint_dt,batches,lanes,n_trajectories,bond,seed,"
               << "threads,projection_count,elapsed_seconds\n";
        output << MODEL_VERSION << ",transient," << p.common.n << ","
               << p.common.T1 << "," << p.common.Tn << "," << GAMMA << ","
               << p.common.dt << "," << p.total_time << ","
               << p.checkpoint_dt << "," << p.common.batches << ","
               << LANES << "," << trajectories << "," << p.common.bond << ","
               << p.common.seed << "," << p.common.threads << ","
               << projection_count << "," << elapsed << "\n";
    }
}

void run_tau(const TauParams& p) {
    ensure_parent_directory(p.common.prefix);
    const int trajectories = p.common.batches * LANES;
    std::cout << "Finite-tau current-window probe\n"
              << "model=" << MODEL_VERSION
              << " n=" << p.common.n
              << " T1=" << p.common.T1
              << " Tn=" << p.common.Tn
              << " bond=" << p.common.bond
              << " dt=" << p.common.dt
              << " burnin=" << p.burnin
              << " measure=" << p.measure
              << " tau_block=" << p.tau_block
              << " blocks=" << p.block_count
              << " trajectories=" << trajectories
              << " seed=" << p.common.seed
              << " threads=" << p.common.threads << "\n";

    const auto start = std::chrono::steady_clock::now();
    std::vector<double> block_samples(
        static_cast<std::size_t>(p.block_count) * trajectories, 0.0);
    std::vector<std::uint64_t> projections(p.common.batches, 0);
    std::atomic<int> completed{0};

    #pragma omp parallel for schedule(static) num_threads(p.common.threads)
    for (int batch = 0; batch < p.common.batches; ++batch) {
        auto result = run_tau_batch(p, batch);
        for (int block = 0; block < p.block_count; ++block) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int trajectory = batch * LANES + lane;
                block_samples[
                    static_cast<std::size_t>(block) * trajectories +
                    trajectory] =
                    result.block_current_samples[block * LANES + lane];
            }
        }
        projections[batch] = result.projection_count;
        const int done = ++completed;
        const int stride = std::max(1, p.common.batches / 10);
        if (done == p.common.batches || done % stride == 0) {
            #pragma omp critical
            std::cout << "completed " << done << "/" << p.common.batches
                      << " batches\n";
        }
    }

    const auto stop = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(stop - start).count();

    {
        auto output = open_output(p.common.prefix + "_blocks.csv");
        output << "trajectory_id";
        for (int block = 0; block < p.block_count; ++block) {
            output << ",block_" << block << "_current";
        }
        output << "\n";
        for (int i = 0; i < trajectories; ++i) {
            output << i;
            for (int block = 0; block < p.block_count; ++block) {
                output << "," << block_samples[
                    static_cast<std::size_t>(block) * trajectories + i];
            }
            output << "\n";
        }
    }
    {
        auto output = open_output(p.common.prefix + "_prefix_tau_samples.csv");
        output << "trajectory_id";
        for (int k = 1; k <= p.block_count; ++k) {
            output << ",tau_" << (k * p.tau_block) << "_current";
        }
        output << "\n";
        for (int i = 0; i < trajectories; ++i) {
            output << i;
            double cumulative = 0.0;
            for (int block = 0; block < p.block_count; ++block) {
                cumulative += block_samples[
                    static_cast<std::size_t>(block) * trajectories + i];
                output << "," << cumulative / static_cast<double>(block + 1);
            }
            output << "\n";
        }
    }
    {
        const std::uint64_t projection_count =
            std::accumulate(projections.begin(), projections.end(),
                            std::uint64_t{0});
        auto output = open_output(p.common.prefix + "_summary.csv");
        output << "model_version,mode,n,T1,Tn,gamma,dt,burnin,measure,"
               << "tau_block,block_count,batches,lanes,n_trajectories,bond,"
               << "seed,threads,projection_count,elapsed_seconds\n";
        output << MODEL_VERSION << ",tau," << p.common.n << ","
               << p.common.T1 << "," << p.common.Tn << "," << GAMMA << ","
               << p.common.dt << "," << p.burnin << "," << p.measure << ","
               << p.tau_block << "," << p.block_count << ","
               << p.common.batches << "," << LANES << "," << trajectories
               << "," << p.common.bond << "," << p.common.seed << ","
               << p.common.threads << "," << projection_count << ","
               << elapsed << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            throw std::invalid_argument(
                "first argument must be mode: transient or tau");
        }
        const std::string mode = argv[1];
        if (mode == "transient") {
            run_transient(parse_transient(argc, argv));
        } else if (mode == "tau") {
            run_tau(parse_tau(argc, argv));
        } else {
            throw std::invalid_argument("unknown mode: " + mode);
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
