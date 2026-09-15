// Fixed-population Feynman--Kac / GKLT cloning sampler for the medium entropy
// production of the boundary-driven resonant NLS chain.
//
// This is a separate implementation from NLS_entropy_ft.cpp.  It uses the
// same projection-free Cartesian dynamics, physical energy, alternating bath
// splitting, and exact discrete bath-heat increments.  Existing production
// code and raw data are not modified.
//
// For each selection interval Delta, clone i accumulates Delta Sigma_i and is
// assigned the Feynman--Kac weight
//
//     w_i = exp(-k Delta Sigma_i).
//
// The population is then systematic-resampled at fixed size N_c.  The SCGF
// estimator is the accumulated log mean weight divided by elapsed time.
// Descendants receive independent fresh random-number streams after every
// selection event.  Genealogical diagnostics are written at every event.
//
// Build on Apple Silicon:
//   clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
//     -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     flux/NLS_entropy_cloning.cpp -o flux/entropy_cloning
//
// Self-test:
//   ./flux/entropy_cloning selftest
//
// Naive NLS run (path weights exp(-k Delta Sigma)):
//   ./flux/entropy_cloning clone T1 Tn n clones burnin observation_time \
//       selection_time dt k seed threads out_prefix [bond]
//
// Tilted-generator NLS run (modified drift plus Feynman--Kac potential):
//   ./flux/entropy_cloning guided T1 Tn n clones burnin observation_time \
//       selection_time dt k seed threads out_prefix [bond]
//
// Small-chain endpoint-density pilot (currently n=2 only):
//   ./flux/entropy_cloning endpoints T1 Tn n streams burnin block_time \
//       blocks_per_stream dt seed threads out_prefix

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double GAMMA = 0.1;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_MIDPOINT_ITERATIONS = 20;
constexpr double MIDPOINT_TOLERANCE = 2.0e-13;
constexpr const char* MODEL_VERSION = "gibbs-cartesian-entropy-cloning-v1";

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

    std::pair<double, double> gaussian_pair() {
        const double u1 = uniform_open();
        const double u2 = uniform_open();
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double angle = TWO_PI * u2;
        return {radius * std::cos(angle), radius * std::sin(angle)};
    }
};

std::uint64_t mixed_seed(std::uint64_t base,
                         std::uint64_t tag1,
                         std::uint64_t tag2) {
    std::uint64_t x = base ^ (0xD2B74407B1CE6E93ULL * (tag1 + 1ULL));
    x ^= 0xCA5A826395121157ULL * (tag2 + 1ULL);
    return Xoshiro256pp::splitmix64(x);
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
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> action;
    std::vector<double> square_real;
    std::vector<double> square_imag;
    std::vector<double> force_real;
    std::vector<double> force_imag;
    double total_action = 0.0;

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
    std::vector<double> old_x;
    std::vector<double> old_y;
    std::vector<double> guess_x;
    std::vector<double> guess_y;
    std::vector<double> candidate_x;
    std::vector<double> candidate_y;

    explicit MidpointWorkspace(int n)
        : midpoint(n),
          old_x(n),
          old_y(n),
          guess_x(n),
          guess_y(n),
          candidate_x(n),
          candidate_y(n) {}
};

struct Clone {
    State state;
    MidpointWorkspace workspace;
    Xoshiro256pp rng;
    std::uint64_t root_id = 0;

    explicit Clone(int n) : state(n), workspace(n) {}
};

void copy_state(const State& source, State& destination) {
    destination.x = source.x;
    destination.y = source.y;
    destination.action = source.action;
    destination.square_real = source.square_real;
    destination.square_imag = source.square_imag;
    destination.force_real = source.force_real;
    destination.force_imag = source.force_imag;
    destination.total_action = source.total_action;
}

void initialize_state(State& state, double T1, double Tn) {
    const double initial_action =
        std::sqrt(0.5 * (T1 + Tn) / static_cast<double>(state.n));
    const double initial_amplitude = std::sqrt(initial_action);
    std::fill(state.x.begin(), state.x.end(), initial_amplitude);
    std::fill(state.y.begin(), state.y.end(), 0.0);
}

void compute_force(State& state) {
    state.total_action = 0.0;
    for (int j = 0; j < state.n; ++j) {
        state.action[j] = state.x[j] * state.x[j] + state.y[j] * state.y[j];
        state.square_real[j] =
            state.x[j] * state.x[j] - state.y[j] * state.y[j];
        state.square_imag[j] = 2.0 * state.x[j] * state.y[j];
        state.total_action += state.action[j];
    }
    for (int j = 0; j < state.n; ++j) {
        double neighbor_real = 0.0;
        double neighbor_imag = 0.0;
        if (j > 0) {
            neighbor_real += state.square_real[j - 1];
            neighbor_imag += state.square_imag[j - 1];
        }
        if (j + 1 < state.n) {
            neighbor_real += state.square_real[j + 1];
            neighbor_imag += state.square_imag[j + 1];
        }
        const double onsite = 2.0 * state.total_action - state.action[j];
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

double physical_energy(const State& state) {
    double sum_action_squared = 0.0;
    double interaction = 0.0;
    for (int j = 0; j < state.n; ++j) {
        sum_action_squared += state.action[j] * state.action[j];
        if (j > 0) {
            interaction +=
                state.square_real[j] * state.square_real[j - 1] +
                state.square_imag[j] * state.square_imag[j - 1];
        }
    }
    return 0.5 * state.total_action * state.total_action -
           0.25 * sum_action_squared + interaction;
}

double action_current(const State& state, int bond) {
    return 4.0 *
           (state.square_imag[bond] * state.square_real[bond - 1] -
            state.square_real[bond] * state.square_imag[bond - 1]);
}

struct MidpointResult {
    double energy_error = 0.0;
    int iterations = 0;
    bool converged = false;
};

MidpointResult hamiltonian_midpoint_step(State& state,
                                         MidpointWorkspace& workspace,
                                         double dt) {
    compute_force(state);
    const double energy_before = physical_energy(state);
    for (int j = 0; j < state.n; ++j) {
        workspace.old_x[j] = state.x[j];
        workspace.old_y[j] = state.y[j];
        workspace.guess_x[j] = state.x[j] - dt * state.force_imag[j];
        workspace.guess_y[j] = state.y[j] + dt * state.force_real[j];
    }

    MidpointResult result;
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
                workspace.old_x[j] - dt * workspace.midpoint.force_imag[j];
            workspace.candidate_y[j] =
                workspace.old_y[j] + dt * workspace.midpoint.force_real[j];
            maximum_change = std::max(
                maximum_change,
                std::abs(workspace.candidate_x[j] - workspace.guess_x[j]));
            maximum_change = std::max(
                maximum_change,
                std::abs(workspace.candidate_y[j] - workspace.guess_y[j]));
        }
        workspace.guess_x = workspace.candidate_x;
        workspace.guess_y = workspace.candidate_y;
        result.iterations = iteration;
        if (maximum_change < MIDPOINT_TOLERANCE) {
            result.converged = true;
            break;
        }
    }
    state.x = workspace.guess_x;
    state.y = workspace.guess_y;
    compute_force(state);
    result.energy_error = physical_energy(state) - energy_before;
    return result;
}

double bath_step(State& state,
                 int site,
                 double temperature,
                 double normal_x,
                 double normal_y,
                 double dt,
                 double sqrt_dt,
                 double drift_coefficient) {
    const double energy_before = physical_energy(state);
    const double noise_scale = std::sqrt(2.0 * GAMMA * temperature);
    state.x[site] +=
        drift_coefficient * state.force_real[site] * dt +
        noise_scale * sqrt_dt * normal_x;
    state.y[site] +=
        drift_coefficient * state.force_imag[site] * dt +
        noise_scale * sqrt_dt * normal_y;
    compute_force(state);
    return physical_energy(state) - energy_before;
}

struct SegmentResult {
    double entropy = 0.0;
    double q_left = 0.0;
    double q_right = 0.0;
    double action_integral = 0.0;
    double feynman_kac_integral = 0.0;
    double hamiltonian_error = 0.0;
    std::uint64_t midpoint_failures = 0;
    std::uint64_t midpoint_iterations = 0;
    bool finite = true;
};

SegmentResult advance_segment(Clone& clone,
                              double T1,
                              double Tn,
                              double dt,
                              std::int64_t steps,
                              std::int64_t global_step_start,
                              int bond,
                              bool guided = false,
                              double tilt_k = 0.0) {
    SegmentResult result;
    const double sqrt_dt = std::sqrt(dt);
    for (std::int64_t local = 0; local < steps; ++local) {
        compute_force(clone.state);
        result.action_integral += action_current(clone.state, bond) * dt;
        if (guided) {
            // For one bath, with M=sum_j |c_j|^2,
            //   dSigma_m = [gamma |F|^2/T - 4 gamma M] dt
            //              - sqrt(2 gamma/T) F.dW.
            // Exponentially tilting by exp(-k Sigma_m) changes the bath
            // drift from -gamma F to gamma(2k-1)F and leaves the scalar
            // Feynman--Kac potential
            //   gamma(k^2-k)|F|^2/T + 4 gamma k M.
            // The final term below contains both boundary baths.
            const double force_metric =
                (clone.state.force_real.front() *
                     clone.state.force_real.front() +
                 clone.state.force_imag.front() *
                     clone.state.force_imag.front()) /
                    T1 +
                (clone.state.force_real.back() *
                     clone.state.force_real.back() +
                 clone.state.force_imag.back() *
                     clone.state.force_imag.back()) /
                    Tn;
            const double potential =
                GAMMA * (tilt_k * tilt_k - tilt_k) * force_metric +
                8.0 * GAMMA * tilt_k * clone.state.total_action;
            result.feynman_kac_integral += potential * dt;
        }
        const auto normal_a = clone.rng.gaussian_pair();
        const auto normal_b = clone.rng.gaussian_pair();
        const auto midpoint = hamiltonian_midpoint_step(
            clone.state, clone.workspace, dt);
        result.midpoint_iterations +=
            static_cast<std::uint64_t>(midpoint.iterations);
        result.hamiltonian_error += midpoint.energy_error;
        if (!midpoint.converged) {
            ++result.midpoint_failures;
        }

        double q_left = 0.0;
        double q_right = 0.0;
        const double drift_coefficient =
            guided ? GAMMA * (2.0 * tilt_k - 1.0) : -GAMMA;
        const auto global_step = global_step_start + local;
        if (global_step % 2 == 0) {
            q_left = bath_step(
                clone.state, 0, T1, normal_a.first, normal_b.first,
                dt, sqrt_dt, drift_coefficient);
            q_right = bath_step(
                clone.state, clone.state.n - 1, Tn,
                normal_a.second, normal_b.second, dt, sqrt_dt,
                drift_coefficient);
        } else {
            q_right = bath_step(
                clone.state, clone.state.n - 1, Tn,
                normal_a.second, normal_b.second, dt, sqrt_dt,
                drift_coefficient);
            q_left = bath_step(
                clone.state, 0, T1, normal_a.first, normal_b.first,
                dt, sqrt_dt, drift_coefficient);
        }
        result.q_left += q_left;
        result.q_right += q_right;
        result.entropy += -q_left / T1 - q_right / Tn;
        if (!std::isfinite(result.entropy) ||
            !std::isfinite(physical_energy(clone.state))) {
            result.finite = false;
            break;
        }
    }
    return result;
}

struct Parameters {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int clone_count = 0;
    double burnin = 0.0;
    double observation_time = 0.0;
    double selection_time = 0.0;
    double dt = 0.0;
    double k = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    int bond = 0;
    std::int64_t burn_steps = 0;
    std::int64_t selection_steps = 0;
    int selection_events = 0;
};

Parameters parse_clone(int argc, char* argv[]) {
    if (argc < 14 || argc > 15) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " clone T1 Tn n clones burnin observation_time"
                << " selection_time dt k seed threads out_prefix [bond]";
        throw std::invalid_argument(message.str());
    }
    Parameters p;
    p.T1 = std::stod(argv[2]);
    p.Tn = std::stod(argv[3]);
    p.n = std::stoi(argv[4]);
    p.clone_count = std::stoi(argv[5]);
    p.burnin = std::stod(argv[6]);
    p.observation_time = std::stod(argv[7]);
    p.selection_time = std::stod(argv[8]);
    p.dt = std::stod(argv[9]);
    p.k = std::stod(argv[10]);
    p.seed = static_cast<std::uint64_t>(std::stoull(argv[11]));
    p.threads = std::stoi(argv[12]);
    p.prefix = argv[13];
    p.bond = argc == 14 ? p.n / 2 : std::stoi(argv[14]);

    if (!(p.T1 > 0.0) || !(p.Tn > 0.0) || p.n < 2 ||
        p.clone_count < 2 || p.burnin < 0.0 ||
        !(p.observation_time > 0.0) || !(p.selection_time > 0.0) ||
        !(p.dt > 0.0) || !std::isfinite(p.k) || p.threads < 1 ||
        p.bond < 1 || p.bond >= p.n) {
        throw std::invalid_argument("invalid cloning parameters");
    }
    p.burn_steps = checked_step_count(p.burnin, p.dt, "burnin", true);
    p.selection_steps =
        checked_step_count(p.selection_time, p.dt, "selection_time");
    const double raw_events = p.observation_time / p.selection_time;
    p.selection_events = static_cast<int>(std::llround(raw_events));
    if (p.selection_events <= 0 ||
        std::abs(p.selection_events * p.selection_time - p.observation_time) >
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, p.observation_time)) {
        throw std::invalid_argument(
            "observation_time must be divisible by selection_time");
    }
    return p;
}

std::vector<int> systematic_resample(const std::vector<double>& probabilities,
                                     Xoshiro256pp& rng) {
    const int count = static_cast<int>(probabilities.size());
    std::vector<int> parents(count, 0);
    const double offset = rng.uniform_open() / static_cast<double>(count);
    int parent = 0;
    double cumulative = probabilities[0];
    for (int child = 0; child < count; ++child) {
        const double position =
            offset + static_cast<double>(child) / static_cast<double>(count);
        while (position > cumulative && parent + 1 < count) {
            ++parent;
            cumulative += probabilities[parent];
        }
        parents[child] = parent;
    }
    return parents;
}

double count_ess(const std::vector<int>& counts) {
    double sum_squares = 0.0;
    for (const int value : counts) {
        sum_squares += static_cast<double>(value) * value;
    }
    const double total =
        static_cast<double>(std::accumulate(counts.begin(), counts.end(), 0));
    return total * total / sum_squares;
}

struct RunTotals {
    std::uint64_t midpoint_failures = 0;
    std::uint64_t midpoint_iterations = 0;
    double hamiltonian_error = 0.0;
    bool finite = true;
};

void run_clone(const Parameters& p, bool guided) {
    ensure_parent_directory(p.prefix);
    omp_set_num_threads(p.threads);
    std::vector<Clone> clones;
    std::vector<State> next_states;
    clones.reserve(p.clone_count);
    next_states.reserve(p.clone_count);
    for (int i = 0; i < p.clone_count; ++i) {
        clones.emplace_back(p.n);
        next_states.emplace_back(p.n);
        initialize_state(clones.back().state, p.T1, p.Tn);
        clones.back().rng.seed(mixed_seed(p.seed, 0, i));
        clones.back().root_id = static_cast<std::uint64_t>(i);
    }

    std::cout << "Entropy cloning sampler"
              << (guided ? " (tilted-generator guided)\n" : " (naive)\n")
              << "  n=" << p.n << " clones=" << p.clone_count
              << " k=" << p.k << " burnin=" << p.burnin
              << " observation=" << p.observation_time
              << " selection=" << p.selection_time << " dt=" << p.dt
              << " seed=" << p.seed << "\n";

    RunTotals totals;
    std::vector<SegmentResult> segment_results(p.clone_count);
    if (p.burn_steps > 0) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < p.clone_count; ++i) {
            segment_results[i] = advance_segment(
                clones[i], p.T1, p.Tn, p.dt, p.burn_steps, 0, p.bond,
                false, 0.0);
        }
        for (const auto& result : segment_results) {
            totals.midpoint_failures += result.midpoint_failures;
            totals.midpoint_iterations += result.midpoint_iterations;
            totals.hamiltonian_error += result.hamiltonian_error;
            totals.finite = totals.finite && result.finite;
        }
    }
    if (!totals.finite || totals.midpoint_failures != 0) {
        throw std::runtime_error("burn-in failed numerical gates");
    }

    Xoshiro256pp resampling_rng;
    resampling_rng.seed(mixed_seed(p.seed, 0xC10EULL, 0));
    auto timeseries = open_output(p.prefix + "_timeseries.csv");
    timeseries
        << "event,time,log_mean_weight,cumulative_log_normalizer,scgf,"
        << "population_mean_entropy_rate,tilted_mean_entropy_rate,"
        << "population_mean_action_current,population_mean_potential_rate,"
        << "weight_ess,max_weight_fraction,"
        << "unique_parents,parent_count_ess,unique_roots,root_count_ess\n";

    double cumulative_log_normalizer = 0.0;
    double total_observation_entropy = 0.0;
    double total_observation_action = 0.0;
    double minimum_weight_ess = p.clone_count;
    double minimum_parent_ess = p.clone_count;
    double minimum_root_ess = p.clone_count;
    int minimum_unique_parents = p.clone_count;
    int minimum_unique_roots = p.clone_count;

    std::vector<double> log_weights(p.clone_count);
    std::vector<double> probabilities(p.clone_count);
    for (int event = 0; event < p.selection_events; ++event) {
        const std::int64_t global_start =
            p.burn_steps +
            static_cast<std::int64_t>(event) * p.selection_steps;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < p.clone_count; ++i) {
            segment_results[i] = advance_segment(
                clones[i], p.T1, p.Tn, p.dt, p.selection_steps,
                global_start, p.bond, guided, p.k);
        }

        double population_entropy = 0.0;
        double population_action = 0.0;
        double population_potential = 0.0;
        for (int i = 0; i < p.clone_count; ++i) {
            const auto& result = segment_results[i];
            totals.midpoint_failures += result.midpoint_failures;
            totals.midpoint_iterations += result.midpoint_iterations;
            totals.hamiltonian_error += result.hamiltonian_error;
            totals.finite = totals.finite && result.finite;
            population_entropy += result.entropy;
            population_action += result.action_integral;
            population_potential += result.feynman_kac_integral;
            log_weights[i] =
                guided ? result.feynman_kac_integral : -p.k * result.entropy;
        }
        if (!totals.finite || totals.midpoint_failures != 0) {
            throw std::runtime_error("measurement failed numerical gates");
        }

        const double maximum_log_weight =
            *std::max_element(log_weights.begin(), log_weights.end());
        double shifted_sum = 0.0;
        double shifted_square_sum = 0.0;
        for (int i = 0; i < p.clone_count; ++i) {
            probabilities[i] = std::exp(log_weights[i] - maximum_log_weight);
            shifted_sum += probabilities[i];
            shifted_square_sum += probabilities[i] * probabilities[i];
        }
        const double log_mean_weight =
            maximum_log_weight + std::log(shifted_sum) -
            std::log(static_cast<double>(p.clone_count));
        cumulative_log_normalizer += log_mean_weight;
        for (double& probability : probabilities) {
            probability /= shifted_sum;
        }
        const double weight_ess =
            shifted_sum * shifted_sum / shifted_square_sum;
        const double max_weight_fraction =
            *std::max_element(probabilities.begin(), probabilities.end());
        double tilted_entropy = 0.0;
        for (int i = 0; i < p.clone_count; ++i) {
            tilted_entropy += probabilities[i] * segment_results[i].entropy;
        }

        const auto parents =
            systematic_resample(probabilities, resampling_rng);
        std::vector<int> parent_counts(p.clone_count, 0);
        std::vector<int> root_counts(p.clone_count, 0);
        std::vector<std::uint64_t> next_roots(p.clone_count, 0);
        for (int child = 0; child < p.clone_count; ++child) {
            const int parent = parents[child];
            ++parent_counts[parent];
            copy_state(clones[parent].state, next_states[child]);
            next_roots[child] = clones[parent].root_id;
        }
        for (const auto root : next_roots) {
            ++root_counts.at(static_cast<std::size_t>(root));
        }
        const int unique_parents = static_cast<int>(std::count_if(
            parent_counts.begin(), parent_counts.end(),
            [](int value) { return value > 0; }));
        const int unique_roots = static_cast<int>(std::count_if(
            root_counts.begin(), root_counts.end(),
            [](int value) { return value > 0; }));
        const double parent_ess = count_ess(parent_counts);
        const double root_ess = count_ess(root_counts);

        for (int child = 0; child < p.clone_count; ++child) {
            std::swap(clones[child].state, next_states[child]);
            clones[child].root_id = next_roots[child];
            clones[child].rng.seed(mixed_seed(
                p.seed, static_cast<std::uint64_t>(event + 1),
                static_cast<std::uint64_t>(child)));
        }

        const double time = (event + 1) * p.selection_time;
        const double scgf = cumulative_log_normalizer / time;
        const double population_entropy_rate =
            population_entropy /
            (static_cast<double>(p.clone_count) * p.selection_time);
        const double tilted_entropy_rate = tilted_entropy / p.selection_time;
        const double population_action_current =
            population_action /
            (static_cast<double>(p.clone_count) * p.selection_time);
        const double population_potential_rate =
            population_potential /
            (static_cast<double>(p.clone_count) * p.selection_time);
        total_observation_entropy += population_entropy;
        total_observation_action += population_action;
        minimum_weight_ess = std::min(minimum_weight_ess, weight_ess);
        minimum_parent_ess = std::min(minimum_parent_ess, parent_ess);
        minimum_root_ess = std::min(minimum_root_ess, root_ess);
        minimum_unique_parents = std::min(minimum_unique_parents, unique_parents);
        minimum_unique_roots = std::min(minimum_unique_roots, unique_roots);

        timeseries << event + 1 << ',' << time << ',' << log_mean_weight << ','
                   << cumulative_log_normalizer << ',' << scgf << ','
                   << population_entropy_rate << ',' << tilted_entropy_rate
                   << ',' << population_action_current << ','
                   << population_potential_rate << ',' << weight_ess
                   << ',' << max_weight_fraction << ',' << unique_parents
                   << ',' << parent_ess << ',' << unique_roots << ','
                   << root_ess << '\n';
        if ((event + 1) % std::max(1, p.selection_events / 10) == 0 ||
            event + 1 == p.selection_events) {
            std::cout << "  time=" << time << " psi=" << scgf
                      << " weight ESS=" << weight_ess
                      << " roots=" << unique_roots << '\n';
        }
    }

    const double total_steps_per_clone =
        static_cast<double>(p.burn_steps) +
        static_cast<double>(p.selection_steps) * p.selection_events;
    const double mean_midpoint_iterations =
        static_cast<double>(totals.midpoint_iterations) /
        (static_cast<double>(p.clone_count) * total_steps_per_clone);
    const double hamiltonian_error_rate =
        totals.hamiltonian_error /
        (static_cast<double>(p.clone_count) *
         (p.burnin + p.observation_time));
    const double mean_population_entropy_rate =
        total_observation_entropy /
        (static_cast<double>(p.clone_count) * p.observation_time);
    const double mean_population_action_current =
        total_observation_action /
        (static_cast<double>(p.clone_count) * p.observation_time);

    auto summary = open_output(p.prefix + "_summary.csv");
    summary
        << "model_version,mode,T1,Tn,gamma,n,clone_count,burnin,observation_time,"
        << "selection_time,dt,k,seed,threads,bond,selection_events,scgf,"
        << "mean_population_entropy_rate,mean_population_action_current,"
        << "minimum_weight_ess,minimum_parent_count_ess,minimum_root_count_ess,"
        << "minimum_unique_parents,minimum_unique_roots,midpoint_failures,"
        << "mean_midpoint_iterations,mean_hamiltonian_energy_error_rate\n";
    summary << MODEL_VERSION << ',' << (guided ? "guided" : "naive")
            << ',' << p.T1 << ',' << p.Tn << ',' << GAMMA
            << ',' << p.n << ',' << p.clone_count << ',' << p.burnin << ','
            << p.observation_time << ',' << p.selection_time << ',' << p.dt
            << ',' << p.k << ',' << p.seed << ',' << p.threads << ','
            << p.bond << ',' << p.selection_events << ','
            << cumulative_log_normalizer / p.observation_time << ','
            << mean_population_entropy_rate << ','
            << mean_population_action_current << ',' << minimum_weight_ess
            << ',' << minimum_parent_ess << ',' << minimum_root_ess << ','
            << minimum_unique_parents << ',' << minimum_unique_roots << ','
            << totals.midpoint_failures << ',' << mean_midpoint_iterations
            << ',' << hamiltonian_error_rate << '\n';
}

struct EndpointParameters {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int streams = 0;
    double burnin = 0.0;
    double block_time = 0.0;
    int blocks_per_stream = 0;
    double dt = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    std::int64_t burn_steps = 0;
    std::int64_t block_steps = 0;
};

EndpointParameters parse_endpoints(int argc, char* argv[]) {
    if (argc != 13) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " endpoints T1 Tn n streams burnin block_time"
                << " blocks_per_stream dt seed threads out_prefix";
        throw std::invalid_argument(message.str());
    }
    EndpointParameters p;
    p.T1 = std::stod(argv[2]);
    p.Tn = std::stod(argv[3]);
    p.n = std::stoi(argv[4]);
    p.streams = std::stoi(argv[5]);
    p.burnin = std::stod(argv[6]);
    p.block_time = std::stod(argv[7]);
    p.blocks_per_stream = std::stoi(argv[8]);
    p.dt = std::stod(argv[9]);
    p.seed = static_cast<std::uint64_t>(std::stoull(argv[10]));
    p.threads = std::stoi(argv[11]);
    p.prefix = argv[12];
    if (!(p.T1 > 0.0) || !(p.Tn > 0.0) || p.n != 2 || p.streams < 2 ||
        p.burnin < 0.0 || !(p.block_time > 0.0) ||
        p.blocks_per_stream < 1 || !(p.dt > 0.0) || p.threads < 1) {
        throw std::invalid_argument(
            "invalid endpoint parameters (this pilot requires n=2)");
    }
    p.burn_steps = checked_step_count(p.burnin, p.dt, "burnin", true);
    p.block_steps = checked_step_count(p.block_time, p.dt, "block_time");
    return p;
}

struct ReducedState {
    double log_action_1 = 0.0;
    double log_action_2 = 0.0;
    double theta = 0.0;
    double energy = 0.0;
};

ReducedState reduced_state(State& state) {
    compute_force(state);
    if (!(state.action[0] > 0.0) || !(state.action[1] > 0.0)) {
        throw std::runtime_error("zero action in endpoint sampler");
    }
    const double phi_1 = std::atan2(state.y[0], state.x[0]);
    const double phi_2 = std::atan2(state.y[1], state.x[1]);
    return {
        std::log(state.action[0]),
        std::log(state.action[1]),
        std::remainder(2.0 * (phi_2 - phi_1), TWO_PI),
        physical_energy(state),
    };
}

struct EndpointRecord {
    ReducedState start;
    ReducedState end;
    SegmentResult segment;
};

void run_endpoints(const EndpointParameters& p) {
    ensure_parent_directory(p.prefix);
    omp_set_num_threads(p.threads);
    std::vector<Clone> streams;
    streams.reserve(p.streams);
    for (int i = 0; i < p.streams; ++i) {
        streams.emplace_back(p.n);
        initialize_state(streams.back().state, p.T1, p.Tn);
        streams.back().rng.seed(mixed_seed(p.seed, 0, i));
    }

    std::vector<SegmentResult> burn_results(p.streams);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.streams; ++i) {
        burn_results[i] = advance_segment(
            streams[i], p.T1, p.Tn, p.dt, p.burn_steps, 0, 1, false, 0.0);
    }
    std::uint64_t midpoint_failures = 0;
    std::uint64_t midpoint_iterations = 0;
    double hamiltonian_error = 0.0;
    for (const auto& result : burn_results) {
        midpoint_failures += result.midpoint_failures;
        midpoint_iterations += result.midpoint_iterations;
        hamiltonian_error += result.hamiltonian_error;
        if (!result.finite) {
            throw std::runtime_error("endpoint burn-in became nonfinite");
        }
    }
    if (midpoint_failures != 0) {
        throw std::runtime_error("endpoint burn-in failed midpoint gate");
    }

    auto output = open_output(p.prefix + "_blocks.csv");
    output
        << "stream_id,block_id,log_action_1_start,log_action_2_start,"
        << "theta_start,log_action_1_end,log_action_2_end,theta_end,"
        << "energy_start,energy_end,q_left,q_right,delta_energy,"
        << "entropy_medium,energy_balance_error,action_current\n";
    std::vector<EndpointRecord> records(p.streams);
    double entropy_sum = 0.0;
    double action_sum = 0.0;
    double balance_square_sum = 0.0;
    for (int block = 0; block < p.blocks_per_stream; ++block) {
        const std::int64_t global_start =
            p.burn_steps + static_cast<std::int64_t>(block) * p.block_steps;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < p.streams; ++i) {
            records[i].start = reduced_state(streams[i].state);
            records[i].segment = advance_segment(
                streams[i], p.T1, p.Tn, p.dt, p.block_steps,
                global_start, 1, false, 0.0);
            records[i].end = reduced_state(streams[i].state);
        }
        for (int i = 0; i < p.streams; ++i) {
            const auto& record = records[i];
            const auto& result = record.segment;
            if (!result.finite || result.midpoint_failures != 0) {
                throw std::runtime_error(
                    "endpoint measurement failed numerical gate");
            }
            midpoint_failures += result.midpoint_failures;
            midpoint_iterations += result.midpoint_iterations;
            hamiltonian_error += result.hamiltonian_error;
            const double delta_energy = record.end.energy - record.start.energy;
            const double balance_error =
                result.q_left + result.q_right - delta_energy;
            const double current = result.action_integral / p.block_time;
            entropy_sum += result.entropy;
            action_sum += current;
            balance_square_sum += balance_error * balance_error;
            output << i << ',' << block << ','
                   << record.start.log_action_1 << ','
                   << record.start.log_action_2 << ','
                   << record.start.theta << ','
                   << record.end.log_action_1 << ','
                   << record.end.log_action_2 << ','
                   << record.end.theta << ','
                   << record.start.energy << ',' << record.end.energy << ','
                   << result.q_left << ',' << result.q_right << ','
                   << delta_energy << ',' << result.entropy << ','
                   << balance_error << ',' << current << '\n';
        }
    }

    const double samples =
        static_cast<double>(p.streams) * p.blocks_per_stream;
    const double total_steps =
        static_cast<double>(p.streams) *
        (p.burn_steps +
         static_cast<std::int64_t>(p.blocks_per_stream) * p.block_steps);
    auto summary = open_output(p.prefix + "_summary.csv");
    summary
        << "model_version,T1,Tn,n,streams,burnin,block_time,"
        << "blocks_per_stream,dt,seed,threads,samples,mean_entropy_rate,"
        << "mean_action_current,rms_energy_balance_error_rate,"
        << "midpoint_failures,mean_midpoint_iterations,"
        << "mean_hamiltonian_energy_error_rate\n";
    summary << MODEL_VERSION << ',' << p.T1 << ',' << p.Tn << ',' << p.n
            << ',' << p.streams << ',' << p.burnin << ',' << p.block_time
            << ',' << p.blocks_per_stream << ',' << p.dt << ',' << p.seed
            << ',' << p.threads << ',' << static_cast<std::int64_t>(samples)
            << ',' << entropy_sum / (samples * p.block_time)
            << ',' << action_sum / samples
            << ',' << std::sqrt(balance_square_sum / samples) / p.block_time
            << ',' << midpoint_failures
            << ',' << static_cast<double>(midpoint_iterations) / total_steps
            << ',' << hamiltonian_error /
                (samples * (p.burnin / p.blocks_per_stream + p.block_time))
            << '\n';
}

double gaussian_cloning_trial(int population,
                              int events,
                              double interval,
                              double mean_rate,
                              double k,
                              std::uint64_t seed) {
    std::vector<Xoshiro256pp> rng(population);
    for (int i = 0; i < population; ++i) {
        rng[i].seed(mixed_seed(seed, 0, i));
    }
    Xoshiro256pp selection_rng;
    selection_rng.seed(mixed_seed(seed, 0xA11CEULL, 0));
    const double mean = mean_rate * interval;
    const double sigma = std::sqrt(2.0 * mean);
    double log_normalizer = 0.0;
    std::vector<double> probabilities(population);
    for (int event = 0; event < events; ++event) {
        std::vector<double> log_weights(population);
        for (int i = 0; i < population; ++i) {
            const double entropy = mean + sigma * rng[i].gaussian_pair().first;
            log_weights[i] = -k * entropy;
        }
        const double maximum =
            *std::max_element(log_weights.begin(), log_weights.end());
        double sum = 0.0;
        for (int i = 0; i < population; ++i) {
            probabilities[i] = std::exp(log_weights[i] - maximum);
            sum += probabilities[i];
        }
        log_normalizer +=
            maximum + std::log(sum) - std::log(static_cast<double>(population));
        for (double& probability : probabilities) {
            probability /= sum;
        }
        const auto parents = systematic_resample(probabilities, selection_rng);
        for (int i = 0; i < population; ++i) {
            (void)parents[i];
            rng[i].seed(mixed_seed(seed, event + 1, i));
        }
    }
    return log_normalizer / (events * interval);
}

void selftest() {
    // The Gaussian entropy model Sigma_Delta ~ N(m Delta, 2m Delta) obeys
    // psi(k)=m(k^2-k)=psi(1-k).  Multiple deterministic trials reduce the
    // stochastic self-test variance without weakening the tolerance.
    const int trials = 8;
    const int population = 4096;
    const int events = 200;
    const double interval = 0.25;
    const double mean_rate = 0.7;
    for (const double k : {0.0, 0.2, 0.5, 0.8, 1.0}) {
        double estimate = 0.0;
        for (int trial = 0; trial < trials; ++trial) {
            estimate += gaussian_cloning_trial(
                population, events, interval, mean_rate, k,
                20260827ULL + static_cast<std::uint64_t>(trial) * 1009ULL);
        }
        estimate /= trials;
        const double exact = mean_rate * (k * k - k);
        if (std::abs(estimate - exact) > 0.012) {
            std::ostringstream message;
            message << "Gaussian cloning self-test failed at k=" << k
                    << ": estimate=" << estimate << " exact=" << exact;
            throw std::runtime_error(message.str());
        }
    }
    std::cout << "Gaussian cloning self-test: PASS\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            throw std::invalid_argument(
                "mode required: selftest, clone, guided, or endpoints");
        }
        const std::string mode = argv[1];
        if (mode == "selftest") {
            selftest();
        } else if (mode == "clone") {
            run_clone(parse_clone(argc, argv), false);
        } else if (mode == "guided") {
            run_clone(parse_clone(argc, argv), true);
        } else if (mode == "endpoints") {
            run_endpoints(parse_endpoints(argc, argv));
        } else {
            throw std::invalid_argument("unknown mode: " + mode);
        }
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
