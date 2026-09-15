/*
 * Projection-free Cartesian BC1 stationary-profile sampler.
 *
 * Dynamics and RNG ordering are ported from the validated controlled sampler
 * experiments/spectral_gap_controlled_2026-09-07/source/
 * NLS_stationary_autocorr_controlled.cpp.  The only scientific extensions are
 * runtime chain length and per-trajectory stationary means of every action and
 * every bond sine.
 */

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int LANES = 16;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_MIDPOINT_ITERATIONS = 20;
constexpr double MIDPOINT_TOLERANCE = 2.0e-13;
constexpr const char* MODEL_VERSION =
    "nls-bc1-cartesian-controlled-profile-v1";

using A16d = Eigen::Array<double, LANES, 1>;
using AlignedVec = std::vector<A16d, Eigen::aligned_allocator<A16d>>;

AlignedVec zero_arrays(int count) {
    AlignedVec result(static_cast<std::size_t>(count));
    for (auto& value : result) value.setZero();
    return result;
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
        for (auto& v : s) v = splitmix64(value);
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

std::uint64_t stream_seed(std::uint64_t base_seed, std::uint64_t stream_id) {
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

struct State {
    explicit State(int chain_length)
        : n(chain_length),
          x(zero_arrays(n)),
          y(zero_arrays(n)),
          action(zero_arrays(n)),
          square_real(zero_arrays(n)),
          square_imag(zero_arrays(n)),
          force_real(zero_arrays(n)),
          force_imag(zero_arrays(n)) {}

    int n;
    AlignedVec x;
    AlignedVec y;
    AlignedVec action;
    AlignedVec square_real;
    AlignedVec square_imag;
    AlignedVec force_real;
    AlignedVec force_imag;
    A16d total_action = A16d::Zero();
};

struct MidpointWorkspace {
    explicit MidpointWorkspace(int n)
        : midpoint(n),
          old_x(zero_arrays(n)),
          old_y(zero_arrays(n)),
          guess_x(zero_arrays(n)),
          guess_y(zero_arrays(n)),
          candidate_x(zero_arrays(n)),
          candidate_y(zero_arrays(n)) {}

    State midpoint;
    AlignedVec old_x;
    AlignedVec old_y;
    AlignedVec guess_x;
    AlignedVec guess_y;
    AlignedVec candidate_x;
    AlignedVec candidate_y;
};

void initialize_state(State& state, double T1, double Tn) {
    const double initial_action = std::sqrt(0.5 * (T1 + Tn) / state.n);
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
        A16d nr = A16d::Zero();
        A16d ni = A16d::Zero();
        if (j > 0) {
            nr += state.square_real[j - 1];
            ni += state.square_imag[j - 1];
        }
        if (j + 1 < state.n) {
            nr += state.square_real[j + 1];
            ni += state.square_imag[j + 1];
        }
        const A16d onsite = 2.0 * state.total_action - state.action[j];
        state.force_real[j] =
            onsite * state.x[j] + 2.0 * (nr * state.x[j] + ni * state.y[j]);
        state.force_imag[j] =
            onsite * state.y[j] + 2.0 * (ni * state.x[j] - nr * state.y[j]);
    }
}

bool hamiltonian_midpoint_step(State& state,
                               MidpointWorkspace& workspace,
                               double dt,
                               int& iterations) {
    compute_force(state);
    for (int j = 0; j < state.n; ++j) {
        workspace.old_x[j] = state.x[j];
        workspace.old_y[j] = state.y[j];
        workspace.guess_x[j] = state.x[j] - dt * state.force_imag[j];
        workspace.guess_y[j] = state.y[j] + dt * state.force_real[j];
    }

    bool converged = false;
    for (iterations = 1; iterations <= MAX_MIDPOINT_ITERATIONS; ++iterations) {
        for (int j = 0; j < state.n; ++j) {
            workspace.midpoint.x[j] =
                0.5 * (workspace.old_x[j] + workspace.guess_x[j]);
            workspace.midpoint.y[j] =
                0.5 * (workspace.old_y[j] + workspace.guess_y[j]);
        }
        compute_force(workspace.midpoint);
        double maximum_change = 0.0;
        for (int j = 0; j < state.n; ++j) {
            workspace.candidate_x[j] = workspace.old_x[j] -
                                       dt * workspace.midpoint.force_imag[j];
            workspace.candidate_y[j] = workspace.old_y[j] +
                                       dt * workspace.midpoint.force_real[j];
            maximum_change = std::max(
                maximum_change,
                (workspace.candidate_x[j] - workspace.guess_x[j]).abs().maxCoeff());
            maximum_change = std::max(
                maximum_change,
                (workspace.candidate_y[j] - workspace.guess_y[j]).abs().maxCoeff());
        }
        for (int j = 0; j < state.n; ++j) {
            workspace.guess_x[j] = workspace.candidate_x[j];
            workspace.guess_y[j] = workspace.candidate_y[j];
        }
        if (maximum_change < MIDPOINT_TOLERANCE) {
            converged = true;
            break;
        }
    }
    for (int j = 0; j < state.n; ++j) {
        state.x[j] = workspace.guess_x[j];
        state.y[j] = workspace.guess_y[j];
    }
    compute_force(state);
    return converged;
}

void bath_step(State& state,
               int site,
               double temperature,
               double gamma,
               const A16d& normal_x,
               const A16d& normal_y,
               double dt) {
    const double noise_scale = std::sqrt(2.0 * gamma * temperature * dt);
    state.x[site] += -gamma * state.force_real[site] * dt + noise_scale * normal_x;
    state.y[site] += -gamma * state.force_imag[site] * dt + noise_scale * normal_y;
    compute_force(state);
}

bool finite_state(const State& state) {
    for (int j = 0; j < state.n; ++j) {
        if (!state.x[j].isFinite().all() || !state.y[j].isFinite().all()) {
            return false;
        }
    }
    return true;
}

AlignedVec observe_bond_sines(const State& state) {
    AlignedVec result = zero_arrays(state.n - 1);
    for (int j = 0; j + 1 < state.n; ++j) {
        for (int lane = 0; lane < LANES; ++lane) {
            const double phi_left = std::atan2(state.y[j](lane), state.x[j](lane));
            const double phi_right =
                std::atan2(state.y[j + 1](lane), state.x[j + 1](lane));
            result[j](lane) = std::sin(2.0 * (phi_right - phi_left));
        }
    }
    return result;
}

struct Config {
    fs::path prefix;
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    double gamma = 0.1;
    double dt = 0.0;
    double burn = 0.0;
    double measure = 0.0;
    double sample_interval = 0.01;
    std::uint64_t base_seed = 0;
    std::uint64_t stream_offset = 0;
};

std::int64_t exact_steps(double time, double dt, const char* label) {
    const auto result = static_cast<std::int64_t>(std::llround(time / dt));
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(time));
    if (result < 0 || std::abs(result * dt - time) > tolerance) {
        throw std::invalid_argument(std::string(label) + "/dt must be integral");
    }
    return result;
}

Config parse_args(int argc, char** argv) {
    if (argc != 12) {
        throw std::invalid_argument(
            std::string("Usage: ") + argv[0] +
            " OUT_PREFIX T1 Tn n gamma dt burn measure sample_interval"
            " base_seed stream_offset");
    }
    Config c;
    c.prefix = argv[1];
    c.T1 = std::stod(argv[2]);
    c.Tn = std::stod(argv[3]);
    c.n = std::stoi(argv[4]);
    c.gamma = std::stod(argv[5]);
    c.dt = std::stod(argv[6]);
    c.burn = std::stod(argv[7]);
    c.measure = std::stod(argv[8]);
    c.sample_interval = std::stod(argv[9]);
    c.base_seed = std::stoull(argv[10]);
    c.stream_offset = std::stoull(argv[11]);
    if (!(c.T1 > 0.0 && c.Tn > 0.0 && c.n >= 2 && c.gamma > 0.0 &&
          c.dt > 0.0 && c.burn >= 0.0 && c.measure > 0.0 &&
          c.sample_interval > 0.0)) {
        throw std::invalid_argument("invalid argument");
    }
    return c;
}

struct Result {
    explicit Result(int n)
        : sum_action(zero_arrays(n)),
          sum_sine(zero_arrays(n - 1)),
          final_action(zero_arrays(n)),
          final_sine(zero_arrays(n - 1)) {
        min_action.fill(std::numeric_limits<double>::infinity());
    }

    AlignedVec sum_action;
    AlignedVec sum_sine;
    AlignedVec final_action;
    AlignedVec final_sine;
    std::array<double, LANES> min_action{};
    std::array<std::uint64_t, LANES> zero_radius_events{};
    std::uint64_t samples = 0;
    std::uint64_t midpoint_failure_steps = 0;
    std::uint64_t midpoint_iteration_sum = 0;
    int midpoint_iteration_max = 0;
    int nonfinite = 0;
};

void track_actions(const State& state, Result& result) {
    for (int j = 0; j < state.n; ++j) {
        for (int lane = 0; lane < LANES; ++lane) {
            const double value = state.action[j](lane);
            result.min_action[lane] = std::min(result.min_action[lane], value);
            if (value == 0.0) ++result.zero_radius_events[lane];
        }
    }
}

void record_profile_sample(const State& state, Result& result) {
    const AlignedVec sines = observe_bond_sines(state);
    for (int j = 0; j < state.n; ++j) result.sum_action[j] += state.action[j];
    for (int j = 0; j + 1 < state.n; ++j) result.sum_sine[j] += sines[j];
    ++result.samples;
}

Result run(const Config& c) {
    const std::int64_t burn_steps = exact_steps(c.burn, c.dt, "burn");
    const std::int64_t measure_steps = exact_steps(c.measure, c.dt, "measure");
    const std::int64_t sample_steps =
        exact_steps(c.sample_interval, c.dt, "sample_interval");
    if (sample_steps <= 0 || measure_steps % sample_steps != 0) {
        throw std::invalid_argument("measurement/sample step mismatch");
    }

    State state(c.n);
    MidpointWorkspace workspace(c.n);
    Result result(c.n);
    initialize_state(state, c.T1, c.Tn);
    compute_force(state);
    track_actions(state, result);

    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        rng[lane].seed(stream_seed(c.base_seed, c.stream_offset + lane));
    }

    A16d normal_left_x, normal_left_y, normal_right_x, normal_right_y;
    const std::int64_t total_steps = burn_steps + measure_steps;
    for (std::int64_t step = 0; step < total_steps; ++step) {
        if (step == burn_steps) record_profile_sample(state, result);

        // Keep the validated controlled sampler's exact normal ordering.
        fill_gaussian_pair(rng, normal_left_x, normal_right_x);
        fill_gaussian_pair(rng, normal_left_y, normal_right_y);

        int iterations = 0;
        const bool converged =
            hamiltonian_midpoint_step(state, workspace, c.dt, iterations);
        result.midpoint_iteration_sum += static_cast<std::uint64_t>(iterations);
        result.midpoint_iteration_max =
            std::max(result.midpoint_iteration_max, iterations);
        if (!converged) ++result.midpoint_failure_steps;

        if (step % 2 == 0) {
            bath_step(state, 0, c.T1, c.gamma,
                      normal_left_x, normal_left_y, c.dt);
            bath_step(state, c.n - 1, c.Tn, c.gamma,
                      normal_right_x, normal_right_y, c.dt);
        } else {
            bath_step(state, c.n - 1, c.Tn, c.gamma,
                      normal_right_x, normal_right_y, c.dt);
            bath_step(state, 0, c.T1, c.gamma,
                      normal_left_x, normal_left_y, c.dt);
        }

        if (!finite_state(state)) {
            result.nonfinite = 1;
            break;
        }
        track_actions(state, result);
        if (step >= burn_steps &&
            (step - burn_steps + 1) % sample_steps == 0) {
            record_profile_sample(state, result);
        }
    }

    if (!result.nonfinite) {
        const auto expected =
            static_cast<std::uint64_t>(measure_steps / sample_steps + 1);
        if (result.samples != expected) {
            throw std::runtime_error("internal profile sample count mismatch");
        }
        const AlignedVec final_sines = observe_bond_sines(state);
        for (int j = 0; j < c.n; ++j) result.final_action[j] = state.action[j];
        for (int j = 0; j + 1 < c.n; ++j) result.final_sine[j] = final_sines[j];
    }
    return result;
}

struct MeanSE {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double se = std::numeric_limits<double>::quiet_NaN();
};

MeanSE mean_se(const A16d& values) {
    MeanSE result;
    result.mean = values.mean();
    const double sumsq = (values - result.mean).square().sum();
    result.se = std::sqrt(sumsq / (LANES - 1) / LANES);
    return result;
}

void write_header(std::ostream& out, const Config& c, const Result& result) {
    out << "# model=" << MODEL_VERSION << '\n'
        << "# bc=BC1 canonical Cartesian bath\n"
        << "# dynamics=dc_j=i*grad_j(H/2)*dt; boundaries add"
           " -gamma*grad_j(H/2)*dt+sqrt(2*gamma*T_j)*dW\n"
        << "# integrator=implicit-midpoint Hamiltonian then alternating-order"
           " Euler-Maruyama baths\n"
        << "# projection=none floor=none adaptive_step=no standard_trig=yes\n"
        << "# n=" << c.n << " T1=" << std::setprecision(17) << c.T1
        << " Tn=" << c.Tn << " gamma=" << c.gamma << " dt=" << c.dt
        << " burnin=" << c.burn << " measure=" << c.measure
        << " sample_interval=" << c.sample_interval
        << " profile_samples_per_trajectory=" << result.samples
        << " base_seed=" << c.base_seed
        << " stream_offset=" << c.stream_offset
        << " trajectories=" << LANES << '\n';
}

void write_outputs(const Config& c, const Result& result, double elapsed) {
    if (c.prefix.has_parent_path()) fs::create_directories(c.prefix.parent_path());
    const double inv_samples = 1.0 / static_cast<double>(result.samples);

    std::ofstream trajectory(c.prefix.string() + "_trajectory_profiles.csv");
    write_header(trajectory, c, result);
    trajectory << "stream_id,stream_seed,j,mean_I,mean_sin_theta\n";
    for (int lane = 0; lane < LANES; ++lane) {
        const auto stream = c.stream_offset + static_cast<std::uint64_t>(lane);
        const auto seed = stream_seed(c.base_seed, stream);
        for (int j = 0; j < c.n; ++j) {
            trajectory << stream << ',' << seed << ',' << j + 1 << ','
                       << std::setprecision(17)
                       << result.sum_action[j](lane) * inv_samples << ',';
            if (j + 1 < c.n) {
                trajectory << result.sum_sine[j](lane) * inv_samples;
            } else {
                trajectory << "nan";
            }
            trajectory << '\n';
        }
    }

    std::ofstream profile(c.prefix.string() + "_profile.csv");
    write_header(profile, c, result);
    profile << "j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta\n";
    for (int j = 0; j < c.n; ++j) {
        const MeanSE action = mean_se(result.sum_action[j] * inv_samples);
        profile << j + 1 << ',' << std::setprecision(17)
                << action.mean << ',' << action.se << ',';
        if (j + 1 < c.n) {
            const MeanSE sine = mean_se(result.sum_sine[j] * inv_samples);
            profile << sine.mean << ',' << sine.se;
        } else {
            profile << "nan,nan";
        }
        profile << '\n';
    }

    std::ofstream final_state(c.prefix.string() + "_final_state.csv");
    write_header(final_state, c, result);
    final_state << "stream_id,stream_seed,j,I,sin_theta\n";
    for (int lane = 0; lane < LANES; ++lane) {
        const auto stream = c.stream_offset + static_cast<std::uint64_t>(lane);
        const auto seed = stream_seed(c.base_seed, stream);
        for (int j = 0; j < c.n; ++j) {
            final_state << stream << ',' << seed << ',' << j + 1 << ','
                        << std::setprecision(17) << result.final_action[j](lane)
                        << ',';
            if (j + 1 < c.n) final_state << result.final_sine[j](lane);
            else final_state << "nan";
            final_state << '\n';
        }
    }

    std::ofstream metadata(c.prefix.string() + "_metadata.csv");
    write_header(metadata, c, result);
    metadata << "stream_id,stream_seed,projection_count,floor_count,"
                "zero_radius_events,min_sampled_action,"
                "midpoint_failure_steps,midpoint_iteration_max,nonfinite\n";
    for (int lane = 0; lane < LANES; ++lane) {
        const auto stream = c.stream_offset + static_cast<std::uint64_t>(lane);
        metadata << stream << ',' << stream_seed(c.base_seed, stream)
                 << ",0,0," << result.zero_radius_events[lane] << ','
                 << std::setprecision(17) << result.min_action[lane] << ','
                 << result.midpoint_failure_steps << ','
                 << result.midpoint_iteration_max << ',' << result.nonfinite
                 << '\n';
    }

    const double global_min =
        *std::min_element(result.min_action.begin(), result.min_action.end());
    std::uint64_t zero_total = 0;
    for (const auto value : result.zero_radius_events) zero_total += value;
    std::ofstream summary(c.prefix.string() + "_summary.csv");
    summary << "model,n,T1,Tn,gamma,dt,burnin,measure,sample_interval,"
               "profile_samples_per_trajectory,base_seed,stream_offset,"
               "trajectories,projection_count,floor_count,zero_radius_events,"
               "min_sampled_action,midpoint_failure_steps,"
               "midpoint_iteration_max,nonfinite,elapsed_seconds\n";
    summary << MODEL_VERSION << ',' << c.n << ',' << std::setprecision(17)
            << c.T1 << ',' << c.Tn << ',' << c.gamma << ',' << c.dt << ','
            << c.burn << ',' << c.measure << ',' << c.sample_interval << ','
            << result.samples << ',' << c.base_seed << ',' << c.stream_offset
            << ',' << LANES << ",0,0," << zero_total << ',' << global_min << ','
            << result.midpoint_failure_steps << ','
            << result.midpoint_iteration_max << ',' << result.nonfinite << ','
            << elapsed << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config c = parse_args(argc, argv);
        const auto start = std::chrono::steady_clock::now();
        std::cout << std::setprecision(17)
                  << "model=" << MODEL_VERSION
                  << " bc=BC1 n=" << c.n
                  << " coordinates=cartesian-double projection=none floor=none"
                  << " adaptive_step=no standard_trig=yes"
                  << " T1=" << c.T1 << " Tn=" << c.Tn
                  << " gamma=" << c.gamma << " dt=" << c.dt
                  << " burn=" << c.burn << " measure=" << c.measure
                  << " sample_interval=" << c.sample_interval
                  << " base_seed=" << c.base_seed
                  << " stream_offset=" << c.stream_offset
                  << " trajectories=" << LANES << '\n';

        const Result result = run(c);
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        write_outputs(c, result, elapsed);

        const double global_min =
            *std::min_element(result.min_action.begin(), result.min_action.end());
        std::uint64_t zero_total = 0;
        for (const auto value : result.zero_radius_events) zero_total += value;
        std::cout << "profile_samples_per_trajectory=" << result.samples
                  << " projection_count=0 floor_count=0"
                  << " zero_radius_events=" << zero_total
                  << " min_sampled_action=" << std::setprecision(17)
                  << global_min
                  << " midpoint_failure_steps="
                  << result.midpoint_failure_steps
                  << " midpoint_iteration_max="
                  << result.midpoint_iteration_max
                  << " nonfinite=" << result.nonfinite
                  << " elapsed_seconds=" << elapsed << '\n';
        return result.nonfinite || result.midpoint_failure_steps ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
