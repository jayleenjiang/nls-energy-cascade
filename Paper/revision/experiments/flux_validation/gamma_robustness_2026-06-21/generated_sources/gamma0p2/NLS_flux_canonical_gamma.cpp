// Canonical Gibbs-preserving boundary-driven NLS chain.
//
// This program measures the finite-time averaged ACTION current
//
//   J_j = 4 I_{j-1} I_j sin(2(phi_j - phi_{j-1}))
//
// across one interior bond.  It intentionally uses double precision and
// Eigen's standard trigonometric functions: publication runs should not mix
// the physical discretization error with an undocumented fast-math error.
//
// Build on Apple Silicon:
//   clang++ -O3 -mcpu=native -std=c++17 \
//     -Xpreprocessor -fopenmp \
//     -I/opt/homebrew/include/eigen3 \
//     -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     NLS_flux_canonical.cpp -o flux_canonical
//
// Usage:
//   ./flux_canonical T1 Tn n batches burnin measure dt seed threads out_prefix [bond]
//
// Example:
//   ./flux_canonical 10 2 20 32 1000 200 0.0005 20260619 8 \
//     results/canonical_n20_dt5e-4
//
// One batch contains LANES=16 independent trajectories.  "bond" is the right
// endpoint j of the bond (j-1,j), with 1 <= j < n; the default is n/2.

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

constexpr double GAMMA = 0.20000000000000001;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr double ACTION_FLOOR = 1.0e-12;
constexpr int LANES = 16;
constexpr int N_BURNIN_MONITORS = 20;
constexpr int N_CURRENT_BLOCKS = 4;
constexpr const char* MODEL_VERSION = "gibbs-canonical-v1-gamma0p2";

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
        // Map to (0,1), avoiding both endpoints in Box--Muller.
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

struct Parameters {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int batches = 0;
    double burnin = 0.0;
    double measure = 0.0;
    double dt = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    int bond = 0;
    std::int64_t burn_steps = 0;
    std::int64_t measure_steps = 0;
};

struct BatchResult {
    AlignedVec action_integral;
    std::array<A16d, N_CURRENT_BLOCKS> current_integral_blocks{};
    std::uint64_t projection_count = 0;
    std::vector<double> burnin_mass;
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
        if (result.empty() || step != result.back()) {
            result.push_back(step);
        }
    }
    return result;
}

BatchResult run_batch(const Parameters& p, int batch_index) {
    BatchResult result;
    result.action_integral.resize(p.n);
    for (auto& value : result.action_integral) {
        value.setZero();
    }
    for (auto& value : result.current_integral_blocks) {
        value.setZero();
    }

    AlignedVec action(p.n), phase(p.n), drift_action(p.n), drift_phase(p.n);
    AlignedVec sin_bond(p.n), cos_bond(p.n);

    const double initial_action =
        std::sqrt(0.5 * (p.T1 + p.Tn) / static_cast<double>(p.n));
    for (int j = 0; j < p.n; ++j) {
        action[j].setConstant(initial_action);
        phase[j].setZero();
    }

    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto trajectory_id =
            static_cast<std::uint64_t>(batch_index) * LANES + lane;
        rng[lane].seed(trajectory_seed(p.seed, trajectory_id));
    }

    const auto checkpoints = monitor_steps(p.burn_steps);
    result.burnin_mass.assign(checkpoints.size(), 0.0);
    std::size_t checkpoint_index = 0;

    const double sqrt_dt = std::sqrt(p.dt);
    const std::int64_t total_steps = p.burn_steps + p.measure_steps;

    for (std::int64_t step = 0; step < total_steps; ++step) {
        A16d total_action = A16d::Zero();
        for (int j = 0; j < p.n; ++j) {
            total_action += action[j];
        }

        for (int j = 1; j < p.n; ++j) {
            const A16d difference =
                wrap_pi(2.0 * (phase[j] - phase[j - 1]));
            sin_bond[j] = difference.sin();
            cos_bond[j] = difference.cos();
        }

        for (int j = 0; j < p.n; ++j) {
            A16d action_left = A16d::Zero();
            A16d action_right = A16d::Zero();
            A16d phase_left = A16d::Zero();
            A16d phase_right = A16d::Zero();

            if (j > 0) {
                action_left = action[j - 1] * sin_bond[j];
                phase_left = 2.0 * action[j - 1] * cos_bond[j];
            }
            if (j + 1 < p.n) {
                action_right = -action[j + 1] * sin_bond[j + 1];
                phase_right = 2.0 * action[j + 1] * cos_bond[j + 1];
            }

            drift_action[j] =
                4.0 * action[j] * (action_left + action_right);
            drift_phase[j] =
                2.0 * total_action - action[j] + phase_left + phase_right;
        }

        // Canonical Gibbs-preserving boundary drift.  The phase convention is
        // dnext=2(phi_boundary-phi_neighbor), hence the plus sign below is
        // equivalent to -2 gamma I_neighbor sin(theta) for
        // theta=2(phi_neighbor-phi_boundary).
        const A16d dleft = wrap_pi(2.0 * (phase[0] - phase[1]));
        drift_action[0] += 2.0 * GAMMA *
            (2.0 * p.T1 -
             (2.0 * total_action * action[0] - action[0].square() +
              2.0 * action[1] * action[0] * dleft.cos()));
        drift_phase[0] +=
            2.0 * GAMMA * action[1] * dleft.sin();

        const A16d dright =
            wrap_pi(2.0 * (phase[p.n - 1] - phase[p.n - 2]));
        drift_action[p.n - 1] += 2.0 * GAMMA *
            (2.0 * p.Tn -
             (2.0 * total_action * action[p.n - 1] -
              action[p.n - 1].square() +
              2.0 * action[p.n - 2] * action[p.n - 1] * dright.cos()));
        drift_phase[p.n - 1] +=
            2.0 * GAMMA * action[p.n - 2] * dright.sin();

        const A16d current =
            4.0 * action[p.bond - 1] * action[p.bond] *
            sin_bond[p.bond];

        if (step >= p.burn_steps) {
            const auto measurement_step = step - p.burn_steps;
            const int block = static_cast<int>(
                measurement_step * N_CURRENT_BLOCKS / p.measure_steps);
            result.current_integral_blocks[block] += current * p.dt;
            for (int j = 0; j < p.n; ++j) {
                result.action_integral[j] += action[j] * p.dt;
            }
        }

        A16d normal_action_left, normal_action_right;
        A16d normal_phase_left, normal_phase_right;
        fill_gaussian_pair(rng, normal_action_left, normal_action_right);
        fill_gaussian_pair(rng, normal_phase_left, normal_phase_right);

        const A16d action_left_start = action[0].max(ACTION_FLOOR);
        const A16d action_right_start =
            action[p.n - 1].max(ACTION_FLOOR);

        for (int j = 0; j < p.n; ++j) {
            action[j] += drift_action[j] * p.dt;
            phase[j] = wrap_pi(phase[j] + drift_phase[j] * p.dt);
        }

        // The sqrt(2) factors are required by the Gibbs-preserving SDE.
        action[0] +=
            2.0 * (2.0 * GAMMA * p.T1 * action_left_start).sqrt() *
            sqrt_dt * normal_action_left;
        action[p.n - 1] +=
            2.0 * (2.0 * GAMMA * p.Tn * action_right_start).sqrt() *
            sqrt_dt * normal_action_right;
        phase[0] = wrap_pi(
            phase[0] +
            (2.0 * GAMMA * p.T1 / action_left_start).sqrt() *
                sqrt_dt * normal_phase_left);
        phase[p.n - 1] = wrap_pi(
            phase[p.n - 1] +
            (2.0 * GAMMA * p.Tn / action_right_start).sqrt() *
                sqrt_dt * normal_phase_right);

        for (int j = 0; j < p.n; ++j) {
            result.projection_count +=
                static_cast<std::uint64_t>((action[j] < ACTION_FLOOR).count());
            action[j] = action[j].max(ACTION_FLOOR);
        }

        if (checkpoint_index < checkpoints.size() &&
            step + 1 == checkpoints[checkpoint_index]) {
            A16d monitored_total = A16d::Zero();
            for (int j = 0; j < p.n; ++j) {
                monitored_total += action[j];
            }
            result.burnin_mass[checkpoint_index] =
                monitored_total.mean();
            ++checkpoint_index;
        }
    }

    return result;
}

std::int64_t checked_step_count(double time, double dt,
                                const std::string& label) {
    const double raw = time / dt;
    const auto steps = static_cast<std::int64_t>(std::llround(raw));
    if (steps <= 0) {
        throw std::invalid_argument(label + " must contain at least one step");
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

Parameters parse_parameters(int argc, char* argv[]) {
    if (argc < 11 || argc > 12) {
        std::ostringstream message;
        message
            << "Usage: " << argv[0]
            << " T1 Tn n batches burnin measure dt seed threads out_prefix [bond]\n"
            << "Example: " << argv[0]
            << " 10 2 20 32 1000 200 0.0005 20260619 8"
            << " results/canonical_n20_dt5e-4\n";
        throw std::invalid_argument(message.str());
    }

    Parameters p;
    p.T1 = std::stod(argv[1]);
    p.Tn = std::stod(argv[2]);
    p.n = std::stoi(argv[3]);
    p.batches = std::stoi(argv[4]);
    p.burnin = std::stod(argv[5]);
    p.measure = std::stod(argv[6]);
    p.dt = std::stod(argv[7]);
    p.seed = static_cast<std::uint64_t>(std::stoull(argv[8]));
    p.threads = std::stoi(argv[9]);
    p.prefix = argv[10];
    p.bond = argc == 12 ? std::stoi(argv[11]) : p.n / 2;

    if (!(p.T1 > 0.0) || !(p.Tn > 0.0)) {
        throw std::invalid_argument("T1 and Tn must be positive");
    }
    if (p.n < 2) {
        throw std::invalid_argument("n must be at least 2");
    }
    if (p.batches <= 0 || p.threads <= 0) {
        throw std::invalid_argument("batches and threads must be positive");
    }
    if (!(p.dt > 0.0) || !(p.burnin > 0.0) || !(p.measure > 0.0)) {
        throw std::invalid_argument("burnin, measure, and dt must be positive");
    }
    if (p.bond < 1 || p.bond >= p.n) {
        throw std::invalid_argument("bond must satisfy 1 <= bond < n");
    }
    p.burn_steps = checked_step_count(p.burnin, p.dt, "burnin");
    p.measure_steps = checked_step_count(p.measure, p.dt, "measure");
    if (p.measure_steps % N_CURRENT_BLOCKS != 0) {
        throw std::invalid_argument(
            "measure/dt must be divisible by four for equal current blocks");
    }
    return p;
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

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Parameters p = parse_parameters(argc, argv);
        ensure_parent_directory(p.prefix);

        const int trajectory_count = p.batches * LANES;
        const double block_duration =
            p.measure / static_cast<double>(N_CURRENT_BLOCKS);
        const auto checkpoints = monitor_steps(p.burn_steps);
        std::vector<double> current_samples(trajectory_count, 0.0);
        std::array<std::vector<double>, N_CURRENT_BLOCKS> current_blocks;
        for (auto& block : current_blocks) {
            block.assign(trajectory_count, 0.0);
        }
        std::vector<double> batch_action(
            static_cast<std::size_t>(p.batches) * p.n, 0.0);
        std::vector<double> batch_burnin(
            static_cast<std::size_t>(p.batches) * checkpoints.size(), 0.0);
        std::vector<std::uint64_t> batch_projections(p.batches, 0);
        std::atomic<int> completed{0};

        std::cout << "Canonical NLS action-current simulation\n"
                  << "model=" << MODEL_VERSION
                  << " n=" << p.n
                  << " T1=" << p.T1
                  << " Tn=" << p.Tn
                  << " bond=" << p.bond
                  << " dt=" << p.dt
                  << " burnin=" << p.burnin
                  << " measure=" << p.measure
                  << " trajectories=" << trajectory_count
                  << " seed=" << p.seed
                  << " threads=" << p.threads << "\n";

        const auto start = std::chrono::steady_clock::now();

        #pragma omp parallel for schedule(static) num_threads(p.threads)
        for (int batch = 0; batch < p.batches; ++batch) {
            BatchResult result = run_batch(p, batch);
            for (int lane = 0; lane < LANES; ++lane) {
                const int trajectory = batch * LANES + lane;
                double integral = 0.0;
                for (int block = 0; block < N_CURRENT_BLOCKS; ++block) {
                    const double block_integral =
                        result.current_integral_blocks[block](lane);
                    current_blocks[block][trajectory] =
                        block_integral / block_duration;
                    integral += block_integral;
                }
                current_samples[trajectory] = integral / p.measure;
            }
            for (int j = 0; j < p.n; ++j) {
                batch_action[static_cast<std::size_t>(batch) * p.n + j] =
                    result.action_integral[j].sum() / p.measure;
            }
            for (std::size_t k = 0; k < checkpoints.size(); ++k) {
                batch_burnin[
                    static_cast<std::size_t>(batch) * checkpoints.size() + k] =
                    result.burnin_mass[k];
            }
            batch_projections[batch] = result.projection_count;

            const int done = ++completed;
            const int stride = std::max(1, p.batches / 10);
            if (done == p.batches || done % stride == 0) {
                #pragma omp critical
                std::cout << "completed " << done << "/" << p.batches
                          << " batches\n";
            }
        }

        const auto stop = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration<double>(stop - start).count();

        const double mean =
            std::accumulate(current_samples.begin(), current_samples.end(), 0.0) /
            trajectory_count;
        double sum_squared = 0.0;
        for (const double value : current_samples) {
            const double difference = value - mean;
            sum_squared += difference * difference;
        }
        const double sample_sd =
            std::sqrt(sum_squared / static_cast<double>(trajectory_count - 1));
        const double standard_error =
            sample_sd / std::sqrt(static_cast<double>(trajectory_count));
        const double half_width = 1.959963984540054 * standard_error;

        std::vector<double> first_half_current(trajectory_count, 0.0);
        std::vector<double> second_half_current(trajectory_count, 0.0);
        std::vector<double> paired_difference(trajectory_count, 0.0);
        for (int i = 0; i < trajectory_count; ++i) {
            first_half_current[i] =
                0.5 * (current_blocks[0][i] + current_blocks[1][i]);
            second_half_current[i] =
                0.5 * (current_blocks[2][i] + current_blocks[3][i]);
            paired_difference[i] =
                second_half_current[i] - first_half_current[i];
        }
        const double first_half_mean =
            std::accumulate(first_half_current.begin(),
                            first_half_current.end(), 0.0) /
            trajectory_count;
        const double second_half_mean =
            std::accumulate(second_half_current.begin(),
                            second_half_current.end(), 0.0) /
            trajectory_count;
        const double paired_difference_mean =
            std::accumulate(paired_difference.begin(),
                            paired_difference.end(), 0.0) /
            trajectory_count;
        double paired_sum_squared = 0.0;
        for (const double value : paired_difference) {
            const double centered = value - paired_difference_mean;
            paired_sum_squared += centered * centered;
        }
        const double paired_difference_sd = std::sqrt(
            paired_sum_squared / static_cast<double>(trajectory_count - 1));
        const double paired_difference_se =
            paired_difference_sd /
            std::sqrt(static_cast<double>(trajectory_count));

        std::vector<double> mean_action(p.n, 0.0);
        for (int batch = 0; batch < p.batches; ++batch) {
            for (int j = 0; j < p.n; ++j) {
                mean_action[j] +=
                    batch_action[static_cast<std::size_t>(batch) * p.n + j];
            }
        }
        for (double& value : mean_action) {
            value /= static_cast<double>(trajectory_count);
        }

        std::vector<double> mean_burnin(checkpoints.size(), 0.0);
        for (int batch = 0; batch < p.batches; ++batch) {
            for (std::size_t k = 0; k < checkpoints.size(); ++k) {
                mean_burnin[k] += batch_burnin[
                    static_cast<std::size_t>(batch) * checkpoints.size() + k];
            }
        }
        for (double& value : mean_burnin) {
            value /= static_cast<double>(p.batches);
        }

        const std::uint64_t projection_count =
            std::accumulate(batch_projections.begin(),
                            batch_projections.end(),
                            std::uint64_t{0});
        const long double update_count =
            static_cast<long double>(trajectory_count) * p.n *
            (p.burn_steps + p.measure_steps);
        const double projection_rate =
            static_cast<double>(projection_count / update_count);
        const double projection_events_per_trajectory_time =
            static_cast<double>(projection_count) /
            (static_cast<double>(trajectory_count) *
             (p.burnin + p.measure));

        {
            auto output = open_output(p.prefix + "_samples.csv");
            output
                << "trajectory_id,time_averaged_action_current,"
                << "first_half_action_current,second_half_action_current,"
                << "block_0_action_current,block_1_action_current,"
                << "block_2_action_current,block_3_action_current\n";
            for (int i = 0; i < trajectory_count; ++i) {
                output << i << "," << current_samples[i] << ","
                       << first_half_current[i] << ","
                       << second_half_current[i];
                for (int block = 0; block < N_CURRENT_BLOCKS; ++block) {
                    output << "," << current_blocks[block][i];
                }
                output << "\n";
            }
        }
        {
            auto output = open_output(p.prefix + "_profile.csv");
            output << "mode,mean_action\n";
            for (int j = 0; j < p.n; ++j) {
                output << j << "," << mean_action[j] << "\n";
            }
        }
        {
            auto output = open_output(p.prefix + "_burnin.csv");
            output << "time,mean_total_action\n";
            for (std::size_t k = 0; k < checkpoints.size(); ++k) {
                output << checkpoints[k] * p.dt << ","
                       << mean_burnin[k] << "\n";
            }
        }
        {
            auto output = open_output(p.prefix + "_summary.csv");
            output
                << "model_version,n,T1,Tn,gamma,dt,burnin,measure,batches,"
                << "lanes,n_trajectories,bond,seed,threads,mean_action_current,"
                << "sample_sd,standard_error,normal95_ci_lower,"
                << "normal95_ci_upper,current_block_duration,"
                << "mean_first_half_current,mean_second_half_current,"
                << "mean_second_minus_first,paired_difference_se,"
                << "action_floor,projection_count,"
                << "projection_rate,projection_events_per_trajectory_time,"
                << "elapsed_seconds\n";
            output
                << MODEL_VERSION << ","
                << p.n << "," << p.T1 << "," << p.Tn << "," << GAMMA << ","
                << p.dt << "," << p.burnin << "," << p.measure << ","
                << p.batches << "," << LANES << "," << trajectory_count << ","
                << p.bond << "," << p.seed << "," << p.threads << ","
                << mean << "," << sample_sd << "," << standard_error << ","
                << mean - half_width << "," << mean + half_width << ","
                << block_duration << "," << first_half_mean << ","
                << second_half_mean << "," << paired_difference_mean << ","
                << paired_difference_se << ","
                << ACTION_FLOOR << "," << projection_count << ","
                << projection_rate << ","
                << projection_events_per_trajectory_time << ","
                << elapsed << "\n";
        }

        std::cout << std::setprecision(10)
                  << "mean action current = " << mean << "\n"
                  << "sample SD = " << sample_sd << "\n"
                  << "standard error = " << standard_error << "\n"
                  << "normal 95% CI = [" << mean - half_width
                  << ", " << mean + half_width << "]\n"
                  << "second-half minus first-half current = "
                  << paired_difference_mean << " +/- "
                  << paired_difference_se << " (paired SE)\n"
                  << "projection count = " << projection_count
                  << " (per-update rate " << projection_rate
                  << ", events/trajectory/time "
                  << projection_events_per_trajectory_time << ")\n"
                  << "elapsed seconds = " << elapsed << "\n"
                  << "wrote prefix " << p.prefix << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}
