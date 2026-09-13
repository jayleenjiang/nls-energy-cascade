// Exact discrete path-probability fluctuation-relation control for the
// Cartesian boundary-driven NLS chain.
//
// This program is deliberately separate from NLS_entropy_ft.cpp and
// NLS_entropy_cloning.cpp.  It does not change the established production
// samplers.  Instead, it treats the split integrator itself as a discrete
// Markov process and accumulates the exact Gaussian transition-kernel ratio
// for every bath substep.
//
// For a forward path omega and its conjugate reverse path Theta omega,
//
//   Sigma_tot^dt = log rho_0(z_0) / rho_R(Theta z_t)
//                  + log dP_F[omega|z_0] / dP_R[Theta omega|Theta z_t].
//
// The initial forward and reverse densities are the same explicitly
// normalized, time-reversal-invariant product complex Gaussian,
//
//   rho_a(z) = product_j exp(-|c_j|^2/a)/(pi a).
//
// The Hamiltonian implicit-midpoint map is volume preserving and reversible.
// Its deterministic contribution therefore cancels.  The reverse macrostep
// applies the forward bath substeps in reverse order and then advances the
// time-reversed state with the same +dt Hamiltonian map.  Indeed,
// H_dt Theta = Theta H_-dt.  At a bath substep z -> z', only one complex site
// changes,
// and the exact local log ratio is
//
//   log K_T(z,z') / K_T(Theta z',Theta z)
//     = (|-delta + gamma F(z') dt|^2
//        - | delta + gamma F(z ) dt|^2)/(4 gamma T dt).
//
// This finite-dt quantity is compared with the continuum heat proxy -dQ/T.
// Agreement is expected only as dt -> 0; the path-ratio IFT itself is exact
// for the discrete forward/reverse pair, up to midpoint-solve tolerance and
// Monte Carlo error.
//
// Build (Apple Silicon with Homebrew OpenMP):
//   clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
//     -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     flux/NLS_discrete_path_ft.cpp -o flux/discrete_path_ft
//
// Self-test:
//   ./flux/discrete_path_ft selftest
//
// Sample forward and reverse ensembles:
//   ./flux/discrete_path_ft sample T_left T_right n trajectories horizon \
//       dt initial_mean_action seed threads out_prefix

#include <omp.h>

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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double GAMMA = 0.1;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr int MAX_MIDPOINT_ITERATIONS = 30;
constexpr double MIDPOINT_TOLERANCE = 2.0e-13;
constexpr const char* MODEL_VERSION = "nls-discrete-path-ft-v2";

struct Xoshiro256pp {
    std::array<std::uint64_t, 4> s{};

    static std::uint64_t splitmix64(std::uint64_t& x) {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE6E5B9ULL;
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
                         std::uint64_t direction,
                         std::uint64_t trajectory) {
    std::uint64_t x = base ^
        (0xD2B74407B1CE6E93ULL * (direction + 1ULL));
    x ^= 0xCA5A826395121157ULL * (trajectory + 1ULL);
    return Xoshiro256pp::splitmix64(x);
}

std::int64_t checked_step_count(double time, double dt) {
    const double raw = time / dt;
    const auto steps = static_cast<std::int64_t>(std::llround(raw));
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(time));
    if (steps <= 0 || std::abs(steps * dt - time) > tolerance) {
        throw std::invalid_argument("horizon/dt must be a positive integer");
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
        : n(n_), x(n_), y(n_), action(n_), square_real(n_),
          square_imag(n_), force_real(n_), force_imag(n_) {}
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
        : midpoint(n), old_x(n), old_y(n), guess_x(n), guess_y(n),
          candidate_x(n), candidate_y(n) {}
};

void compute_force(State& state) {
    state.total_action = 0.0;
    for (int j = 0; j < state.n; ++j) {
        state.action[j] = state.x[j] * state.x[j] +
                          state.y[j] * state.y[j];
        state.square_real[j] = state.x[j] * state.x[j] -
                               state.y[j] * state.y[j];
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

void sample_initial(State& state, double mean_action, Xoshiro256pp& rng) {
    const double scale = std::sqrt(0.5 * mean_action);
    for (int j = 0; j < state.n; ++j) {
        const auto normal = rng.gaussian_pair();
        state.x[j] = scale * normal.first;
        state.y[j] = scale * normal.second;
    }
    compute_force(state);
}

double log_initial_density(const State& state, double mean_action) {
    return -static_cast<double>(state.n) * std::log(PI * mean_action) -
           state.total_action / mean_action;
}

void time_reverse(State& state) {
    for (double& value : state.y) {
        value = -value;
    }
    compute_force(state);
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
                workspace.old_x[j] -
                dt * workspace.midpoint.force_imag[j];
            workspace.candidate_y[j] =
                workspace.old_y[j] +
                dt * workspace.midpoint.force_real[j];
            maximum_change = std::max(
                maximum_change,
                std::abs(workspace.candidate_x[j] -
                         workspace.guess_x[j]));
            maximum_change = std::max(
                maximum_change,
                std::abs(workspace.candidate_y[j] -
                         workspace.guess_y[j]));
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

struct BathResult {
    double heat = 0.0;
    double log_kernel_ratio = 0.0;
};

BathResult bath_step_with_ratio(State& state,
                                int site,
                                double temperature,
                                double normal_x,
                                double normal_y,
                                double dt) {
    compute_force(state);
    const double energy_before = physical_energy(state);
    const double x_before = state.x[site];
    const double y_before = state.y[site];
    const double fx_before = state.force_real[site];
    const double fy_before = state.force_imag[site];
    const double noise_scale = std::sqrt(2.0 * GAMMA * temperature * dt);
    state.x[site] += -GAMMA * fx_before * dt + noise_scale * normal_x;
    state.y[site] += -GAMMA * fy_before * dt + noise_scale * normal_y;
    compute_force(state);

    const double dx = state.x[site] - x_before;
    const double dy = state.y[site] - y_before;
    const double forward_x = dx + GAMMA * fx_before * dt;
    const double forward_y = dy + GAMMA * fy_before * dt;
    const double reverse_x = -dx + GAMMA * state.force_real[site] * dt;
    const double reverse_y = -dy + GAMMA * state.force_imag[site] * dt;
    const double forward_square =
        forward_x * forward_x + forward_y * forward_y;
    const double reverse_square =
        reverse_x * reverse_x + reverse_y * reverse_y;

    BathResult result;
    result.heat = physical_energy(state) - energy_before;
    result.log_kernel_ratio =
        (reverse_square - forward_square) /
        (4.0 * GAMMA * temperature * dt);
    return result;
}

struct TrajectoryRecord {
    double sigma_kernel = 0.0;
    double sigma_endpoint = 0.0;
    double sigma_total = 0.0;
    double sigma_heat = 0.0;
    double q_left = 0.0;
    double q_right = 0.0;
    double delta_energy = 0.0;
    double energy_balance_error = 0.0;
    double kernel_minus_heat = 0.0;
    double hamiltonian_energy_error = 0.0;
    std::uint64_t midpoint_failures = 0;
    std::uint64_t midpoint_iterations = 0;
    bool finite = true;
};

TrajectoryRecord simulate_trajectory(bool reverse,
                                     double T1,
                                     double Tn,
                                     int n,
                                     std::int64_t steps,
                                     double dt,
                                     double mean_action,
                                     std::uint64_t seed) {
    Xoshiro256pp rng;
    rng.seed(seed);
    State state(n);
    MidpointWorkspace workspace(n);
    sample_initial(state, mean_action, rng);
    const double log_density_start =
        log_initial_density(state, mean_action);
    const double energy_start = physical_energy(state);

    TrajectoryRecord result;
    auto apply_bath = [&](int site, double temperature, bool left) {
        const auto normal = rng.gaussian_pair();
        const auto bath = bath_step_with_ratio(
            state, site, temperature, normal.first, normal.second, dt);
        result.sigma_kernel += bath.log_kernel_ratio;
        if (left) {
            result.q_left += bath.heat;
        } else {
            result.q_right += bath.heat;
        }
    };

    for (std::int64_t local = 0; local < steps; ++local) {
        if (!reverse) {
            const auto midpoint =
                hamiltonian_midpoint_step(state, workspace, dt);
            result.hamiltonian_energy_error += midpoint.energy_error;
            result.midpoint_iterations +=
                static_cast<std::uint64_t>(midpoint.iterations);
            if (!midpoint.converged) {
                ++result.midpoint_failures;
            }
            if (local % 2 == 0) {
                apply_bath(0, T1, true);
                apply_bath(n - 1, Tn, false);
            } else {
                apply_bath(n - 1, Tn, false);
                apply_bath(0, T1, true);
            }
        } else {
            const std::int64_t forward_step = steps - 1 - local;
            if (forward_step % 2 == 0) {
                apply_bath(n - 1, Tn, false);
                apply_bath(0, T1, true);
            } else {
                apply_bath(0, T1, true);
                apply_bath(n - 1, Tn, false);
            }
            // The reverse state already represents Theta z.  Advancing it
            // with +dt gives H_dt Theta = Theta H_-dt, i.e. the conjugate of
            // the inverse forward Hamiltonian segment.  Using -dt here would
            // reverse the Hamiltonian part twice.
            const auto midpoint =
                hamiltonian_midpoint_step(state, workspace, dt);
            result.hamiltonian_energy_error += midpoint.energy_error;
            result.midpoint_iterations +=
                static_cast<std::uint64_t>(midpoint.iterations);
            if (!midpoint.converged) {
                ++result.midpoint_failures;
            }
        }
        if (!std::isfinite(result.sigma_kernel) ||
            !std::isfinite(physical_energy(state))) {
            result.finite = false;
            break;
        }
    }

    const double energy_end = physical_energy(state);
    const double log_density_end = log_initial_density(state, mean_action);
    result.sigma_endpoint = log_density_start - log_density_end;
    result.sigma_total = result.sigma_kernel + result.sigma_endpoint;
    result.sigma_heat = -result.q_left / T1 - result.q_right / Tn;
    result.delta_energy = energy_end - energy_start;
    result.energy_balance_error =
        result.q_left + result.q_right - result.delta_energy;
    result.kernel_minus_heat = result.sigma_kernel - result.sigma_heat;
    result.finite = result.finite && std::isfinite(result.sigma_total);
    return result;
}

struct Parameters {
    double T1 = 0.0;
    double Tn = 0.0;
    int n = 0;
    int trajectories = 0;
    double horizon = 0.0;
    double dt = 0.0;
    double mean_action = 0.0;
    std::uint64_t seed = 0;
    int threads = 1;
    std::string prefix;
    std::int64_t steps = 0;
};

Parameters parse_sample(int argc, char* argv[]) {
    if (argc != 12) {
        std::ostringstream message;
        message << "Usage: " << argv[0]
                << " sample T_left T_right n trajectories horizon dt"
                << " initial_mean_action seed threads out_prefix";
        throw std::invalid_argument(message.str());
    }
    Parameters p;
    p.T1 = std::stod(argv[2]);
    p.Tn = std::stod(argv[3]);
    p.n = std::stoi(argv[4]);
    p.trajectories = std::stoi(argv[5]);
    p.horizon = std::stod(argv[6]);
    p.dt = std::stod(argv[7]);
    p.mean_action = std::stod(argv[8]);
    p.seed = static_cast<std::uint64_t>(std::stoull(argv[9]));
    p.threads = std::stoi(argv[10]);
    p.prefix = argv[11];
    if (!(p.T1 > 0.0) || !(p.Tn > 0.0) || p.n < 2 ||
        p.trajectories < 2 || !(p.horizon > 0.0) || !(p.dt > 0.0) ||
        !(p.mean_action > 0.0) || p.threads < 1) {
        throw std::invalid_argument("invalid sample parameters");
    }
    p.steps = checked_step_count(p.horizon, p.dt);
    return p;
}

void write_records(const std::string& path,
                   const std::vector<TrajectoryRecord>& records) {
    auto output = open_output(path);
    output << "trajectory_id,sigma_kernel,sigma_endpoint,sigma_total,"
           << "sigma_heat,q_left,q_right,delta_energy,energy_balance_error,"
           << "kernel_minus_heat,hamiltonian_energy_error,"
           << "midpoint_failures,midpoint_iterations,finite\n";
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto& r = records[i];
        output << i << ',' << r.sigma_kernel << ',' << r.sigma_endpoint
               << ',' << r.sigma_total << ',' << r.sigma_heat << ','
               << r.q_left << ',' << r.q_right << ',' << r.delta_energy
               << ',' << r.energy_balance_error << ','
               << r.kernel_minus_heat << ','
               << r.hamiltonian_energy_error << ','
               << r.midpoint_failures << ',' << r.midpoint_iterations
               << ',' << (r.finite ? 1 : 0) << '\n';
    }
}

double mean_field(const std::vector<TrajectoryRecord>& records,
                  double TrajectoryRecord::* field) {
    double sum = 0.0;
    for (const auto& r : records) {
        sum += r.*field;
    }
    return sum / static_cast<double>(records.size());
}

void run_sample(const Parameters& p) {
    ensure_parent_directory(p.prefix);
    omp_set_num_threads(p.threads);
    std::vector<TrajectoryRecord> forward(p.trajectories);
    std::vector<TrajectoryRecord> reverse(p.trajectories);
    const auto started = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.trajectories; ++i) {
        forward[i] = simulate_trajectory(
            false, p.T1, p.Tn, p.n, p.steps, p.dt, p.mean_action,
            mixed_seed(p.seed, 0, static_cast<std::uint64_t>(i)));
        reverse[i] = simulate_trajectory(
            true, p.T1, p.Tn, p.n, p.steps, p.dt, p.mean_action,
            mixed_seed(p.seed, 1, static_cast<std::uint64_t>(i)));
    }

    const auto finished = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(finished - started).count();
    write_records(p.prefix + "_forward.csv", forward);
    write_records(p.prefix + "_reverse.csv", reverse);

    std::uint64_t failures = 0;
    std::uint64_t nonfinite = 0;
    for (const auto* ensemble : {&forward, &reverse}) {
        for (const auto& r : *ensemble) {
            failures += r.midpoint_failures;
            nonfinite += r.finite ? 0ULL : 1ULL;
        }
    }
    auto summary = open_output(p.prefix + "_summary.csv");
    summary << "model_version,T_left,T_right,n,trajectories_per_direction,"
            << "horizon,dt,steps,initial_mean_action,seed,threads,"
            << "forward_mean_sigma_total,reverse_mean_sigma_total,"
            << "forward_mean_kernel_minus_heat,"
            << "reverse_mean_kernel_minus_heat,midpoint_failures,"
            << "nonfinite_trajectories,elapsed_seconds\n";
    summary << MODEL_VERSION << ',' << p.T1 << ',' << p.Tn << ',' << p.n
            << ',' << p.trajectories << ',' << p.horizon << ',' << p.dt
            << ',' << p.steps << ',' << p.mean_action << ',' << p.seed
            << ',' << p.threads << ','
            << mean_field(forward, &TrajectoryRecord::sigma_total) << ','
            << mean_field(reverse, &TrajectoryRecord::sigma_total) << ','
            << mean_field(forward, &TrajectoryRecord::kernel_minus_heat)
            << ','
            << mean_field(reverse, &TrajectoryRecord::kernel_minus_heat)
            << ',' << failures << ',' << nonfinite << ',' << elapsed << '\n';

    std::cout << "Discrete path-ratio sampler complete\n"
              << "model=" << MODEL_VERSION << " n=" << p.n
              << " trajectories/direction=" << p.trajectories
              << " horizon=" << p.horizon << " dt=" << p.dt
              << " forward mean Sigma="
              << mean_field(forward, &TrajectoryRecord::sigma_total)
              << " reverse mean Sigma="
              << mean_field(reverse, &TrajectoryRecord::sigma_total)
              << " midpoint failures=" << failures
              << " nonfinite=" << nonfinite
              << " elapsed=" << elapsed << " s\n";
}

double max_state_difference(const State& a, const State& b) {
    double maximum = 0.0;
    for (int j = 0; j < a.n; ++j) {
        maximum = std::max(maximum, std::abs(a.x[j] - b.x[j]));
        maximum = std::max(maximum, std::abs(a.y[j] - b.y[j]));
    }
    return maximum;
}

void selftest() {
    State initial(3);
    initial.x = {0.31, -0.47, 0.22};
    initial.y = {-0.19, 0.36, 0.41};
    compute_force(initial);

    State recovered = initial;
    MidpointWorkspace forward_workspace(3);
    MidpointWorkspace backward_workspace(3);
    const auto forward = hamiltonian_midpoint_step(
        recovered, forward_workspace, 5.0e-4);
    const auto backward = hamiltonian_midpoint_step(
        recovered, backward_workspace, -5.0e-4);
    const double midpoint_reversibility =
        max_state_difference(initial, recovered);

    State conjugate_path = initial;
    time_reverse(conjugate_path);
    MidpointWorkspace conjugate_workspace(3);
    const auto conjugate_step = hamiltonian_midpoint_step(
        conjugate_path, conjugate_workspace, 5.0e-4);
    time_reverse(conjugate_path);
    State inverse_path = initial;
    MidpointWorkspace inverse_workspace(3);
    const auto inverse_step = hamiltonian_midpoint_step(
        inverse_path, inverse_workspace, -5.0e-4);
    const double time_reversal_error =
        max_state_difference(conjugate_path, inverse_path);

    State forward_end = initial;
    MidpointWorkspace path_forward_workspace(3);
    const auto path_forward = hamiltonian_midpoint_step(
        forward_end, path_forward_workspace, 5.0e-4);
    time_reverse(forward_end);
    MidpointWorkspace path_reverse_workspace(3);
    const auto path_reverse = hamiltonian_midpoint_step(
        forward_end, path_reverse_workspace, 5.0e-4);
    State conjugate_initial = initial;
    time_reverse(conjugate_initial);
    const double conjugate_path_recovery_error =
        max_state_difference(forward_end, conjugate_initial);

    State bath_state = initial;
    const int site = 0;
    const double temperature = 6.0;
    const double dt = 5.0e-4;
    compute_force(bath_state);
    const double x_before = bath_state.x[site];
    const double y_before = bath_state.y[site];
    const double fx_before = bath_state.force_real[site];
    const double fy_before = bath_state.force_imag[site];
    const auto bath = bath_step_with_ratio(
        bath_state, site, temperature, 0.37, -1.13, dt);
    const double dx = bath_state.x[site] - x_before;
    const double dy = bath_state.y[site] - y_before;
    const double variance = 2.0 * GAMMA * temperature * dt;
    const double forward_square =
        std::pow(dx + GAMMA * fx_before * dt, 2) +
        std::pow(dy + GAMMA * fy_before * dt, 2);
    const double reverse_square =
        std::pow(-dx + GAMMA * bath_state.force_real[site] * dt, 2) +
        std::pow(-dy + GAMMA * bath_state.force_imag[site] * dt, 2);
    const double log_forward =
        -std::log(2.0 * PI * variance) -
        forward_square / (2.0 * variance);
    const double log_reverse =
        -std::log(2.0 * PI * variance) -
        reverse_square / (2.0 * variance);
    const double kernel_error =
        std::abs(bath.log_kernel_ratio - (log_forward - log_reverse));

    const bool passed =
        forward.converged && backward.converged && conjugate_step.converged &&
        inverse_step.converged && path_forward.converged &&
        path_reverse.converged && midpoint_reversibility < 2.0e-11 &&
        time_reversal_error < 2.0e-11 &&
        conjugate_path_recovery_error < 2.0e-11 &&
        kernel_error < 2.0e-13;
    std::cout << std::setprecision(12)
              << "midpoint reversibility error="
              << midpoint_reversibility << '\n'
              << "Hamiltonian time-reversal error="
              << time_reversal_error << '\n'
              << "conjugate-path recovery error="
              << conjugate_path_recovery_error << '\n'
              << "bath kernel-ratio error=" << kernel_error << '\n'
              << (passed ? "SELFTEST PASS\n" : "SELFTEST FAIL\n");
    if (!passed) {
        throw std::runtime_error("discrete path FT self-test failed");
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            throw std::invalid_argument("mode required: selftest or sample");
        }
        const std::string mode = argv[1];
        if (mode == "selftest") {
            selftest();
        } else if (mode == "sample") {
            run_sample(parse_sample(argc, argv));
        } else {
            throw std::invalid_argument("unknown mode: " + mode);
        }
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
