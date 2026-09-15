/*
 * Controlled Cartesian stationary-autocorrelation sampler for the n=3 NLS.
 *
 * The dynamics and fixed-step splitting are derived from
 * flux/NLS_entropy_ft.cpp.  The state is double precision and Cartesian, so
 * there is no positive-action projection, no action floor, no adaptive step,
 * and no polynomial trigonometric approximation.
 */

#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
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

constexpr int N = 3;
constexpr int LANES = 16;
constexpr int NOBS = 12;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_MIDPOINT_ITERATIONS = 20;
constexpr double MIDPOINT_TOLERANCE = 2.0e-13;
constexpr const char* MODEL_VERSION =
    "nls-n3-cartesian-controlled-autocorr-v1";

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
    AlignedVec x{N};
    AlignedVec y{N};
    AlignedVec action{N};
    AlignedVec square_real{N};
    AlignedVec square_imag{N};
    AlignedVec force_real{N};
    AlignedVec force_imag{N};
    A16d total_action = A16d::Zero();
};

struct MidpointWorkspace {
    State midpoint;
    AlignedVec old_x{N};
    AlignedVec old_y{N};
    AlignedVec guess_x{N};
    AlignedVec guess_y{N};
    AlignedVec candidate_x{N};
    AlignedVec candidate_y{N};
};

void initialize_state(State& state, double T1, double T3) {
    const double initial_action = std::sqrt(0.5 * (T1 + T3) / N);
    const double initial_amplitude = std::sqrt(initial_action);
    for (int j = 0; j < N; ++j) {
        state.x[j].setConstant(initial_amplitude);
        state.y[j].setZero();
    }
}

void compute_force(State& state) {
    state.total_action.setZero();
    for (int j = 0; j < N; ++j) {
        state.action[j] = state.x[j].square() + state.y[j].square();
        state.square_real[j] = state.x[j].square() - state.y[j].square();
        state.square_imag[j] = 2.0 * state.x[j] * state.y[j];
        state.total_action += state.action[j];
    }
    for (int j = 0; j < N; ++j) {
        A16d nr = A16d::Zero();
        A16d ni = A16d::Zero();
        if (j > 0) {
            nr += state.square_real[j - 1];
            ni += state.square_imag[j - 1];
        }
        if (j + 1 < N) {
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
    for (int j = 0; j < N; ++j) {
        workspace.old_x[j] = state.x[j];
        workspace.old_y[j] = state.y[j];
        workspace.guess_x[j] = state.x[j] - dt * state.force_imag[j];
        workspace.guess_y[j] = state.y[j] + dt * state.force_real[j];
    }

    bool converged = false;
    for (iterations = 1; iterations <= MAX_MIDPOINT_ITERATIONS; ++iterations) {
        for (int j = 0; j < N; ++j) {
            workspace.midpoint.x[j] =
                0.5 * (workspace.old_x[j] + workspace.guess_x[j]);
            workspace.midpoint.y[j] =
                0.5 * (workspace.old_y[j] + workspace.guess_y[j]);
        }
        compute_force(workspace.midpoint);
        double maximum_change = 0.0;
        for (int j = 0; j < N; ++j) {
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
        for (int j = 0; j < N; ++j) {
            workspace.guess_x[j] = workspace.candidate_x[j];
            workspace.guess_y[j] = workspace.candidate_y[j];
        }
        if (maximum_change < MIDPOINT_TOLERANCE) {
            converged = true;
            break;
        }
    }
    for (int j = 0; j < N; ++j) {
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

double wrap_angle(double value) {
    value = std::fmod(value + PI, TWO_PI);
    if (value < 0.0) value += TWO_PI;
    return value - PI;
}

std::array<float, NOBS> observe(const State& state, int lane) {
    const double phi1 = std::atan2(state.y[0](lane), state.x[0](lane));
    const double phi2 = std::atan2(state.y[1](lane), state.x[1](lane));
    const double phi3 = std::atan2(state.y[2](lane), state.x[2](lane));
    const double theta1 = wrap_angle(2.0 * (phi1 - phi2));
    const double theta3 = wrap_angle(2.0 * (phi3 - phi2));
    const double c1 = std::cos(theta1);
    const double s1 = std::sin(theta1);
    const double c3 = std::cos(theta3);
    const double s3 = std::sin(theta3);
    const double i1 = state.action[0](lane);
    const double i2 = state.action[1](lane);
    const double i3 = state.action[2](lane);
    return {
        static_cast<float>(c1),
        static_cast<float>(s1),
        static_cast<float>(c3),
        static_cast<float>(s3),
        static_cast<float>(std::cos(wrap_angle(theta1 - theta3))),
        static_cast<float>(i2),
        static_cast<float>(i1 + i3),
        static_cast<float>(i1 - i3),
        static_cast<float>(c1 + c3),
        static_cast<float>(c1 - c3),
        static_cast<float>(s1 + s3),
        static_cast<float>(s1 - s3),
    };
}

bool finite_state(const State& state) {
    for (int j = 0; j < N; ++j) {
        if (!state.x[j].isFinite().all() || !state.y[j].isFinite().all()) {
            return false;
        }
    }
    return true;
}

struct Config {
    fs::path outdir;
    double T1 = 0.0;
    double T3 = 0.0;
    double gamma = 0.1;
    double dt = 0.0;
    double burn = 500.0;
    double measure = 1000.0;
    double snapshot = 0.01;
    int batches = 4;
    std::uint64_t base_seed = 0;
    int threads = 4;
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

struct BatchResult {
    std::array<std::uint64_t, LANES> midpoint_failures{};
    std::array<std::uint64_t, LANES> zero_radius_events{};
    std::array<double, LANES> min_action{};
    std::array<int, LANES> nonfinite{};
    std::uint64_t midpoint_iteration_sum = 0;
    int midpoint_iteration_max = 0;
};

BatchResult run_batch(const Config& c,
                      int batch,
                      std::int64_t burn_steps,
                      std::int64_t measure_steps,
                      std::int64_t snapshot_steps,
                      std::int64_t nsnap) {
    BatchResult result;
    result.min_action.fill(std::numeric_limits<double>::infinity());
    State state;
    MidpointWorkspace workspace;
    initialize_state(state, c.T1, c.T3);
    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto id = static_cast<std::uint64_t>(batch * LANES + lane);
        rng[lane].seed(stream_seed(c.base_seed, id));
    }

    std::vector<float> buffer(
        static_cast<std::size_t>(nsnap) * LANES * NOBS);
    auto record_snapshot = [&](std::int64_t snapshot_index) {
        compute_force(state);
        for (int lane = 0; lane < LANES; ++lane) {
            const auto values = observe(state, lane);
            const std::size_t offset =
                (static_cast<std::size_t>(lane) * nsnap + snapshot_index) * NOBS;
            std::copy(values.begin(), values.end(), buffer.begin() + offset);
            for (int j = 0; j < N; ++j) {
                const double action = state.action[j](lane);
                result.min_action[lane] = std::min(result.min_action[lane], action);
                if (action == 0.0) ++result.zero_radius_events[lane];
            }
        }
    };

    A16d normal_left_x, normal_left_y, normal_right_x, normal_right_y;
    const std::int64_t total_steps = burn_steps + measure_steps;
    std::int64_t snapshot_index = 0;
    bool bad_batch = false;
    for (std::int64_t step = 0; step < total_steps; ++step) {
        if (step == burn_steps) record_snapshot(snapshot_index++);

        // Preserve the validated entropy sampler's exact normal-generation
        // ordering: x normals for the two baths, then y normals.
        fill_gaussian_pair(rng, normal_left_x, normal_right_x);
        fill_gaussian_pair(rng, normal_left_y, normal_right_y);
        int iterations = 0;
        const bool converged = hamiltonian_midpoint_step(state, workspace, c.dt, iterations);
        result.midpoint_iteration_sum += static_cast<std::uint64_t>(iterations);
        result.midpoint_iteration_max = std::max(result.midpoint_iteration_max, iterations);
        if (!converged) {
            for (auto& value : result.midpoint_failures) ++value;
        }
        if (step % 2 == 0) {
            bath_step(state, 0, c.T1, c.gamma, normal_left_x, normal_left_y, c.dt);
            bath_step(state, N - 1, c.T3, c.gamma, normal_right_x, normal_right_y, c.dt);
        } else {
            bath_step(state, N - 1, c.T3, c.gamma, normal_right_x, normal_right_y, c.dt);
            bath_step(state, 0, c.T1, c.gamma, normal_left_x, normal_left_y, c.dt);
        }
        if (!finite_state(state)) {
            bad_batch = true;
            for (auto& value : result.nonfinite) value = 1;
            break;
        }
        if (step >= burn_steps && (step - burn_steps + 1) % snapshot_steps == 0) {
            record_snapshot(snapshot_index++);
        }
    }

    if (!bad_batch && snapshot_index != nsnap) {
        throw std::runtime_error("internal snapshot count mismatch");
    }
    fs::create_directories(c.outdir);
    for (int lane = 0; lane < LANES; ++lane) {
        const int stream = batch * LANES + lane;
        std::ostringstream filename;
        filename << "stream_" << std::setw(3) << std::setfill('0') << stream << ".f32";
        const fs::path final_path = c.outdir / filename.str();
        const fs::path partial_path = final_path.string() + ".partial";
        std::ofstream output(partial_path, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("cannot create " + partial_path.string());
        if (!bad_batch) {
            const std::size_t offset = static_cast<std::size_t>(lane) * nsnap * NOBS;
            output.write(reinterpret_cast<const char*>(buffer.data() + offset),
                         static_cast<std::streamsize>(nsnap * NOBS * sizeof(float)));
        }
        output.close();
        if (!bad_batch) fs::rename(partial_path, final_path);
        else fs::rename(partial_path, final_path.string() + ".failed");
    }
    return result;
}

Config parse_args(int argc, char** argv) {
    if (argc != 12) {
        throw std::invalid_argument(
            std::string("Usage: ") + argv[0] +
            " OUTDIR T1 T3 gamma dt burn measure snapshot batches base_seed threads");
    }
    Config c;
    c.outdir = argv[1];
    c.T1 = std::stod(argv[2]);
    c.T3 = std::stod(argv[3]);
    c.gamma = std::stod(argv[4]);
    c.dt = std::stod(argv[5]);
    c.burn = std::stod(argv[6]);
    c.measure = std::stod(argv[7]);
    c.snapshot = std::stod(argv[8]);
    c.batches = std::stoi(argv[9]);
    c.base_seed = std::stoull(argv[10]);
    c.threads = std::stoi(argv[11]);
    if (!(c.T1 > 0.0 && c.T3 > 0.0 && c.gamma > 0.0 && c.dt > 0.0 &&
          c.burn >= 0.0 && c.measure > 0.0 && c.snapshot > 0.0 &&
          c.batches > 0 && c.threads > 0)) {
        throw std::invalid_argument("invalid nonpositive argument");
    }
    return c;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config c = parse_args(argc, argv);
        const auto burn_steps = exact_steps(c.burn, c.dt, "burn");
        const auto measure_steps = exact_steps(c.measure, c.dt, "measure");
        const auto snapshot_steps = exact_steps(c.snapshot, c.dt, "snapshot");
        if (snapshot_steps <= 0 || measure_steps % snapshot_steps != 0) {
            throw std::invalid_argument("measurement/snapshot step mismatch");
        }
        const auto nsnap = measure_steps / snapshot_steps + 1;
        const int streams = c.batches * LANES;
        fs::create_directories(c.outdir);

        std::cout << std::setprecision(17)
                  << "model=" << MODEL_VERSION << " n=3 coordinates=cartesian-double"
                  << " projection=none floor=none adaptive_step=no standard_trig=yes"
                  << " T1=" << c.T1 << " T3=" << c.T3
                  << " gamma=" << c.gamma << " dt=" << c.dt
                  << " burn=" << c.burn << " measure=" << c.measure
                  << " snapshot=" << c.snapshot << " batches=" << c.batches
                  << " streams=" << streams << " nsnap=" << nsnap
                  << " base_seed=" << c.base_seed << " threads=" << c.threads
                  << '\n';

        std::vector<BatchResult> results(c.batches);
        std::atomic<int> completed{0};
#pragma omp parallel for schedule(static) num_threads(c.threads)
        for (int batch = 0; batch < c.batches; ++batch) {
            results[batch] = run_batch(
                c, batch, burn_steps, measure_steps, snapshot_steps, nsnap);
            const int done = ++completed;
#pragma omp critical
            std::cout << "completed " << done << "/" << c.batches
                      << " batches (batch=" << batch << ")\n";
        }

        std::ofstream meta(c.outdir / "metadata.csv");
        meta << "stream_id,stream_seed,projection_count,floor_count,"
                "zero_radius_events,min_sampled_action,midpoint_failure_count,"
                "midpoint_iteration_max,nonfinite\n";
        std::uint64_t total_failures = 0;
        std::uint64_t total_zero = 0;
        int any_nonfinite = 0;
        double global_min_action = std::numeric_limits<double>::infinity();
        for (int batch = 0; batch < c.batches; ++batch) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int stream = batch * LANES + lane;
                const auto seed = stream_seed(c.base_seed, stream);
                const auto& r = results[batch];
                meta << stream << ',' << seed << ",0,0,"
                     << r.zero_radius_events[lane] << ','
                     << std::setprecision(17) << r.min_action[lane] << ','
                     << r.midpoint_failures[lane] << ','
                     << r.midpoint_iteration_max << ',' << r.nonfinite[lane] << '\n';
                total_failures += r.midpoint_failures[lane];
                total_zero += r.zero_radius_events[lane];
                any_nonfinite |= r.nonfinite[lane];
                global_min_action = std::min(global_min_action, r.min_action[lane]);
            }
        }

        std::ofstream format(c.outdir / "FORMAT.txt");
        format << "model=" << MODEL_VERSION << '\n'
               << "state_precision=float64\noutput_precision=float32\n"
               << "projection=none\nfloor=none\nadaptive_step=no\n"
               << "columns=cos_theta1,sin_theta1,cos_theta3,sin_theta3,"
                  "cos_theta1_minus_theta3,I2,I1_plus_I3,I1_minus_I3,"
                  "cos_even_sum,cos_odd_diff,sin_even_sum,sin_odd_diff\n"
               << "shape_per_stream=" << nsnap << "x" << NOBS << '\n'
               << "T1=" << c.T1 << "\nT3=" << c.T3 << "\ngamma=" << c.gamma
               << "\ndt=" << c.dt << "\nburn=" << c.burn
               << "\nmeasure=" << c.measure << "\nsnapshot=" << c.snapshot
               << "\nstreams=" << streams << "\nbase_seed=" << c.base_seed << '\n';

        std::cout << "projection_count=0 floor_count=0 zero_radius_events="
                  << total_zero << " min_sampled_action=" << global_min_action
                  << " midpoint_failure_stream_events=" << total_failures
                  << " any_nonfinite=" << any_nonfinite << '\n';
        return any_nonfinite || total_failures ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
