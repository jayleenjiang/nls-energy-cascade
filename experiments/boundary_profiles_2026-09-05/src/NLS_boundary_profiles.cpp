//SIMD-Optimized Hamiltonian Energy Flux Simulation

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <Eigen/Dense>
#include <omp.h>
#include <random>
#include <array>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

// ============================================================
// Global constants
// ============================================================
static const float gamma_val = 0.1f;
static double T_burnin = 2000.0;
static double T_measure = 2000.0;
static double T_final = 4000.0;
static const double dt = 0.0005f;
static const float PI_HALF_f = 1.5707963268f;

static int N_sample = 16;           // 16*16 = 256 trajectories per seed replicate
static const int N_thread = 10;
static const int LANES = 16;      // SIMD width: 16 floats = 512 bits
static uint32_t base_seed = 20260619u;

typedef Eigen::Array<float,  LANES, 1> A16f;
typedef Eigen::Array<double,  LANES, 1> A16d;
typedef std::vector<A16f, Eigen::aligned_allocator<A16f>> AlignedVec;
typedef std::vector<A16d, Eigen::aligned_allocator<A16d>> AlignedVecD;

static const int N_CHECKPOINTS = 4;

static const char* bc_label(int bc_id) {
    switch (bc_id) {
        case 0: return "BC1_canonical";
        case 1: return "BC2_no_phase_drift";
        case 2: return "BC3b_no_global_with_phase_drift";
        case 3: return "BC3_no_global_no_phase_drift";
        default: return "INVALID";
    }
}

static const char* bc_action_equation(int bc_id) {
    switch (bc_id) {
        case 0:
        case 1:
            return "b_I=2*gamma*(2*T-(2*M*I-I^2+2*I*I_neighbor*cos(delta)))";
        case 2:
        case 3:
            return "b_I=2*gamma*(2*T-I^2)";
        default:
            return "INVALID";
    }
}

static const char* bc_phase_equation(int bc_id) {
    switch (bc_id) {
        case 0:
        case 2:
            return "b_phi=gamma*(2*I_neighbor*sin(delta))";
        case 1:
        case 3:
            return "b_phi=0";
        default:
            return "INVALID";
    }
}

// ============================================================
// Fast vectorized math (operates on 16 floats simultaneously)
// ============================================================

static inline A16f wrap_pi(const A16f& x) {
    // Wrap angle to [-pi, pi] using round-to-nearest multiple of 2pi
    const float invTwoPi = 0.159154943f;
    const float twoPi = 6.283185307f;
    return x - (x * invTwoPi).round() * twoPi;
}

static inline A16f fast_sin(const A16f& x) {
    // Polynomial approximation of sin(x) for x in [-pi, pi]
    // Max error ~0.001, sufficient for stochastic simulation
    const float B = 1.27323954f;   // 4/pi
    const float C = -0.40528473f;  // -4/pi^2
    const float P = 0.225f;
    auto y = B * x + C * x * x.abs();
    return P * (y * y.abs() - y) + y;
}

static inline A16f fast_sin_w(const A16f& x) {
    // sin(x) for arbitrary x (wraps first)
    return fast_sin(wrap_pi(x));
}

static inline A16f fast_cos_w(const A16f& x) {
    // cos(x) = sin(x + pi/2) for arbitrary x
    return fast_sin(wrap_pi(x + PI_HALF_f));
}

// ============================================================
// RNG: xoshiro128+ variant, generates 16 independent uniforms
// ============================================================

struct RNGState {
    uint32_t s[4];
};

static inline A16f next_u01_16(RNGState& state) {
    A16f res;
    for (int i = 0; i < LANES; ++i) {
        const uint32_t t0 = state.s[0] + state.s[3];
        const uint32_t result = (t0 << 7) | (t0 >> 25);
        const uint32_t t = state.s[1] << 9;
        state.s[2] ^= state.s[0];
        state.s[3] ^= state.s[1];
        state.s[1] ^= state.s[2];
        state.s[0] ^= state.s[3];
        state.s[2] ^= t;
        state.s[3] = (state.s[3] << 11) | (state.s[3] >> 21);
        res(i) = (result >> 9) * 0.00000011920929f;  // map to [0, 1)
    }
    return res;
}

static inline void box_muller_16(RNGState& state, A16f& n1, A16f& n2) {
    // Generate two arrays of 16 independent standard normals
    A16f u1 = next_u01_16(state).max(1e-10f);  // avoid log(0)
    A16f u2 = next_u01_16(state);
    const float twoPi = 6.283185307f;
    A16f radius = (-2.0f * u1.log()).sqrt();
    A16f theta = twoPi * u2;
    n1 = radius * fast_cos_w(theta);
    n2 = radius * fast_sin_w(theta);
}

static RNGState initialize_rng(int rank, int sample_idx) {
    // Deterministic seeding from (base_seed, sample_idx) via splitmix64.
    // Do not use OpenMP rank here: schedule(dynamic) may assign batches to
    // different threads across runs.
    (void)rank;
    RNGState state;
    uint64_t z = static_cast<uint64_t>(base_seed)
               ^ (0xD2B74407B1CE6E93ULL * static_cast<uint64_t>(sample_idx + 1))
               ^ 0x9E3779B97F4A7C15ULL;
    auto splitmix64 = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t result = z;
        result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9ULL;
        result = (result ^ (result >> 27)) * 0x94D049BB133111EBULL;
        return result ^ (result >> 31);
    };
    uint64_t part1 = splitmix64();
    uint64_t part2 = splitmix64();
    state.s[0] = static_cast<uint32_t>(part1);
    state.s[1] = static_cast<uint32_t>(part1 >> 32);
    state.s[2] = static_cast<uint32_t>(part2);
    state.s[3] = static_cast<uint32_t>(part2 >> 32);
    return state;
}

// ============================================================
// Simulation result (16 trajectories per batch)
// ============================================================

struct SimResult16 {
    AlignedVecD accumulated_energy;  // time integrals, one independent value per lane
    AlignedVecD accumulated_sin;     // n-1 bond-sine time integrals
    std::array<AlignedVecD, N_CHECKPOINTS> burn_snapshot_energy;
    std::array<AlignedVecD, N_CHECKPOINTS> burn_snapshot_sin;
    std::array<double, N_CHECKPOINTS> burn_snapshot_time{};
    std::array<AlignedVecD, N_CHECKPOINTS> checkpoint_energy;
    std::array<AlignedVecD, N_CHECKPOINTS> checkpoint_sin;
    std::array<double, N_CHECKPOINTS> checkpoint_time{};
    std::array<bool, LANES> valid{};
    std::array<double, LANES> max_left_action{};
    std::array<double, LANES> max_right_action{};
    A16d accumulated_flux;
    double measurement_time;
    std::uint64_t measurement_steps;
    std::uint64_t projection_count;
};

// ============================================================
// Core simulation: 16 independent trajectories via SIMD
// ============================================================

void run_simulation_16_old(int n, int case_num, float T1, float Tn, int flux_mode,
                       RNGState& rng, SimResult16& result) {
    // State arrays: I[j] and phi[j] each hold 16 trajectory values
    AlignedVec I(n), phi(n);
    for (int j = 0; j < n; j++) {
        I[j].setConstant(0.1f);
        phi[j].setZero();
    }
    I[0].setConstant(1.0f);
    

    // Preallocate temporaries (avoid allocation in hot loop)
    AlignedVec drift_I(n), drift_phi(n);
    AlignedVec sin_d(n), cos_d(n);  // sin_d[j] = sin(2*(phi[j]-phi[j-1])), j=1..n-1
    AlignedVec sin_observe(n);      // recomputed after each update for profile sampling

    
    double current_time = 0.0;
    bool measuring = false;

    // ---- Main time-stepping loop ----
    while(current_time < T_final)
    {

        // 1) Compute total action (sum of I_j across oscillators, per lane)
        A16f total_mass = A16f::Zero();
        for (int j = 0; j < n; j++) total_mass += I[j];

        // 2) Precompute sin/cos of phase differences between neighbors
        //    sin_d[j] = sin(2*(phi[j] - phi[j-1]))
        //    cos_d[j] = cos(2*(phi[j] - phi[j-1]))
        //    Only need j=1,...,n-1 (the n-1 bonds)
        for (int j = 1; j < n; j++) {
            A16f d = 2.0f * (phi[j] - phi[j-1]);
            sin_d[j] = fast_sin_w(d);
            cos_d[j] = fast_cos_w(d);
        }

        // 3) Compute drift_I and drift_phi for each oscillator
        //    Reuse precomputed sin/cos (each bond used by both neighbors):
        //      sin(2*(phi[j] - phi[j-1])) = sin_d[j]
        //      sin(2*(phi[j] - phi[j+1])) = -sin_d[j+1]  (antisymmetry of sin)
        //      cos(2*(phi[j] - phi[j+1])) = cos_d[j+1]   (symmetry of cos)
        for (int j = 0; j < n; j++) {
            A16f term_I_left, term_I_right;
            A16f term_phi_left, term_phi_right;

            if (j > 0) {
                term_I_left = I[j-1] * sin_d[j];
                term_phi_left = 2.0f * I[j-1] * cos_d[j];
            } else {
                term_I_left.setZero();
                term_phi_left.setZero();
            }

            if (j < n - 1) {
                term_I_right = I[j+1] * (-sin_d[j+1]);
                term_phi_right = 2.0f * I[j+1] * cos_d[j+1];
            } else {
                term_I_right.setZero();
                term_phi_right.setZero();
            }

            drift_I[j] = 4.0f * I[j] * (term_I_left + term_I_right);
            drift_phi[j] = (2.0f * total_mass - I[j]) + term_phi_left + term_phi_right;
        }

        // 4) Compute energy flux at specified mode (before boundary modification)
        //    flux = 4 * I[flux_mode-1] * I[flux_mode] * sin(2*(phi[flux_mode] - phi[flux_mode-1]))
        A16f flux_j;
        if (flux_mode > 0 && flux_mode < n) {
            flux_j = 4.0f * I[flux_mode - 1] * I[flux_mode] * sin_d[flux_mode];
        } else {
            flux_j.setZero();
        }
        //std::cout<<flux_j.transpose()<<std::endl;
        // 5) Boundary conditions (heat bath coupling at j=0 and j=n-1)
        A16f sin_next_0, cos_next_0, sin_prev_n, cos_prev_n;
        A16f I_next_0, I_prev_n;

        if (n > 1) {
            I_next_0 = I[1];
            sin_next_0 = -sin_d[1];   // sin(2*(phi[0]-phi[1])) = -sin(2*(phi[1]-phi[0]))
            cos_next_0 = cos_d[1];    // cos(2*(phi[0]-phi[1])) = cos(2*(phi[1]-phi[0]))
            I_prev_n = I[n - 2];
            sin_prev_n = sin_d[n - 1]; // sin(2*(phi[n-1]-phi[n-2]))
            cos_prev_n = cos_d[n - 1]; // cos(2*(phi[n-1]-phi[n-2]))
        } else {
            I_next_0.setZero();
            sin_next_0.setZero();
            cos_next_0.setOnes();
            I_prev_n.setZero();
            sin_prev_n.setZero();
            cos_prev_n.setOnes();
        }

        switch (case_num) {
            case 0: // Canonical BC (drift on both I and phi at boundaries)
                drift_I[0] += 2.0f * gamma_val * (
                    2.0f * T1 - (2.0f * total_mass * I[0] - I[0] * I[0]
                    + 2.0f * I_next_0 * I[0] * cos_next_0));
                drift_phi[0] += gamma_val * (2.0f * I_next_0 * sin_next_0);

                drift_I[n-1] += 2.0f * gamma_val * (
                    2.0f * Tn - (2.0f * total_mass * I[n-1] - I[n-1] * I[n-1]
                    + 2.0f * I_prev_n * I[n-1] * cos_prev_n));
                drift_phi[n-1] += gamma_val * (2.0f * I_prev_n * sin_prev_n);
                break;

            case 1: // Drift on I only at boundaries
                drift_I[0] += 2.0f * gamma_val * (
                    2.0f * T1 - (2.0f * total_mass * I[0] - I[0] * I[0]
                    + 2.0f * I_next_0 * I[0] * cos_next_0));

                drift_I[n-1] += 2.0f * gamma_val * (
                    2.0f * Tn - (2.0f * total_mass * I[n-1] - I[n-1] * I[n-1]
                    + 2.0f * I_prev_n * I[n-1] * cos_prev_n));
                break;

            case 2: // No global term in I drift
                drift_phi[0] += gamma_val * (2.0f * I_next_0 * sin_next_0);
                drift_phi[n-1] += gamma_val * (2.0f * I_prev_n * sin_prev_n);
                drift_I[0] += 2.0f * gamma_val * (2.0f * T1 - I[0] * I[0]);
                drift_I[n-1] += 2.0f * gamma_val * (2.0f * Tn - I[n-1] * I[n-1]);
                break;

            case 3: // BC3: no global term and no phase drift
                drift_I[0] += 2.0f * gamma_val * (2.0f * T1 - I[0] * I[0]);
                drift_I[n-1] += 2.0f * gamma_val * (2.0f * Tn - I[n-1] * I[n-1]);
                break;
        }

        // 6) Generate stochastic noise (Box-Muller: 4 independent normal arrays)
        A16f nI0, nIn, nPhi0, nPhin;
        box_muller_16(rng, nI0, nIn);
        box_muller_16(rng, nPhi0, nPhin);


        // Compute noise amplitudes BEFORE updating state (use current I values)
        A16f noise_I_0 = 2.0f * (2.0f * gamma_val * T1 * I[0].max(1e-10f)).sqrt() * nI0;
        A16f noise_I_n = 2.0f * (2.0f * gamma_val * Tn * I[n-1].max(1e-10f)).sqrt() * nIn;
        A16f noise_phi_0 = (2.0f * gamma_val * T1 / I[0].max(1e-10f)).sqrt() * nPhi0;
        A16f noise_phi_n = (2.0f * gamma_val * Tn / I[n-1].max(1e-10f)).sqrt() * nPhin;

        float max_drift = 1.0f;
        for (int j = 0; j < n; j++) {
            if(max_drift < drift_I[j].cwiseAbs().maxCoeff())
            {
                max_drift = drift_I[j].cwiseAbs().maxCoeff();
            }
            if(max_drift < drift_phi[j].cwiseAbs().maxCoeff())
            {
                max_drift = drift_phi[j].cwiseAbs().maxCoeff();
            }
        }
        double ddt = std::min(1.0/double(max_drift), dt);
        ddt = std::max(1e-5, ddt);
//        std::cout<<" time = "<<current_time<<std::endl;
//        std::cout<<"max drift = "<<max_drift<<" ddt = "<<ddt<<std::endl;
        const double sqrt_ddt_val = sqrt(ddt);
        // 7) Euler-Maruyama update: deterministic part
        for (int j = 0; j < n; j++) {
            I[j] += drift_I[j] * float(ddt);
            phi[j] += drift_phi[j] * float(ddt);
        }
        
        

        // 7b) Add stochastic part at boundaries
        I[0] += noise_I_0*float(sqrt_ddt_val);
        I[n-1] += noise_I_n*float(sqrt_ddt_val);
        phi[0] += noise_phi_0*float(sqrt_ddt_val);
        phi[n-1] += noise_phi_n*float(sqrt_ddt_val);

        // 8) Detect non-finite lanes, then clamp action variables to positive.
        // Invalid trajectories are retained only in the diagnostics count and
        // are excluded from every reported profile statistic.
        std::array<bool, LANES> newly_invalid{};
        for (int j = 0; j < n; j++) {
            if (!(I[j].isFinite().all() && phi[j].isFinite().all())) {
                for (int l = 0; l < LANES; ++l) {
                    if (!std::isfinite(I[j](l)) || !std::isfinite(phi[j](l))) {
                        newly_invalid[l] = true;
                    }
                }
            }
            result.projection_count += static_cast<std::uint64_t>((I[j] < 1e-10f).count());
            I[j] = I[j].max(1e-10f);
            phi[j] = wrap_pi(phi[j]);  // <-- THIS IS THE CRITICAL FIX
        }
        for (int l = 0; l < LANES; ++l) {
            if (newly_invalid[l]) {
                result.valid[l] = false;
                for (int j = 0; j < n; ++j) {
                    I[j](l) = 0.1f;
                    phi[j](l) = 0.0f;
                }
            }
            if (result.valid[l]) {
                result.max_left_action[l] = std::max(
                    result.max_left_action[l], static_cast<double>(I[0](l)));
                result.max_right_action[l] = std::max(
                    result.max_right_action[l], static_cast<double>(I[n - 1](l)));
            }
        }
/*
        float max_I = 0;
        for(int j = 0; j < n; j++)
        {
            if(max_I < I[j](0))
                max_I = I[j](0);
        }
        std::cout<<"max value of I in lane 0 = "<<max_I<<std::endl;
 */
        current_time += ddt;
        /*
        int flag = 0;
        for(int j = 0; j < n; j++)
        {
            if(phi[j].hasNaN() || I[j].hasNaN())
            {
                std::cout<<"time = "<<current_time<<std::endl;
                std::cout<<"mag: "<<std::endl;
                for(int k = 0; k < n; k++)
                {
                    std::cout<<I[k].transpose()<<std::endl;
                }
                std::cout<<"phase: "<<std::endl;
                for(int k = 0; k < n; k++)
                {
                    std::cout<<phi[k].transpose()<<std::endl;
                }
                flag++;
                break;
            }
        }
        if(flag != 0 )
        {
            std::cout<<" nan error!"<<std::endl;
            break;
        }
*/

        // Recompute bond sines at the post-update state. This keeps the action
        // and bond-angle profile observations at the same physical time.
        for (int j = 1; j < n; ++j) {
            const A16f d = 2.0f * (phi[j] - phi[j - 1]);
            sin_observe[j] = fast_sin_w(d);
        }

        // Save trajectory-wise snapshots during burn-in. The four fixed times
        // are predeclared fractions of the requested burn-in and are not chosen
        // after inspecting the stationary profile.
        if (!measuring && T_burnin > 0.0) {
            for (int checkpoint = 0; checkpoint < N_CHECKPOINTS; ++checkpoint) {
                const double target = T_burnin * (checkpoint + 1) / N_CHECKPOINTS;
                if (result.burn_snapshot_time[checkpoint] == 0.0 &&
                    current_time >= target) {
                    result.burn_snapshot_time[checkpoint] = current_time;
                    for (int j = 0; j < n; ++j)
                        result.burn_snapshot_energy[checkpoint][j] = I[j].cast<double>();
                    for (int j = 0; j < n - 1; ++j)
                        result.burn_snapshot_sin[checkpoint][j] =
                            sin_observe[j + 1].cast<double>();
                }
            }
        }

        // 9) Accumulate measurements after burn-in
        if (!measuring && current_time >= T_burnin) {
            measuring = true;
        }

        if (measuring) {
            for (int j = 0; j < n; j++) {
                result.accumulated_energy[j] += I[j].cast<double>() * ddt;
            }
            for (int j = 0; j < n - 1; ++j) {
                // theta_j = 2*(phi_{j+1}-phi_j).
                result.accumulated_sin[j] += sin_observe[j + 1].cast<double>() * ddt;
            }
            result.accumulated_flux += flux_j.cast<double>() * ddt;
            result.measurement_time += ddt;
            ++result.measurement_steps;

            for (int checkpoint = 0; checkpoint < N_CHECKPOINTS; ++checkpoint) {
                const double target = T_measure * (checkpoint + 1) / N_CHECKPOINTS;
                if (result.checkpoint_time[checkpoint] == 0.0 &&
                    result.measurement_time >= target) {
                    result.checkpoint_time[checkpoint] = result.measurement_time;
                    const double inv_time = 1.0 / result.measurement_time;
                    for (int j = 0; j < n; ++j) {
                        result.checkpoint_energy[checkpoint][j] =
                            result.accumulated_energy[j] * inv_time;
                    }
                    for (int j = 0; j < n - 1; ++j) {
                        result.checkpoint_sin[checkpoint][j] =
                            result.accumulated_sin[j] * inv_time;
                    }
                }
            }
        }
    }
    //std::cout<< "max drift = "<<max_drift<<std::endl;
    /*
    std::cout<<"current time = "<<current_time<<std::endl;
    for(int i = 0; i < 16; i++)
    {
        std::cout<<"local flux = "<<result.accumulated_flux(i)/(T_final - T_burnin)<<std::endl;
    }
      */

     
}




// ============================================================
// Output/statistics helpers and main
// ============================================================

struct MeanSE {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double se = std::numeric_limits<double>::quiet_NaN();
    std::size_t count = 0;
};

static MeanSE mean_se(const std::vector<double>& values,
                      const std::vector<unsigned char>& valid) {
    MeanSE result;
    double sum = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (valid[i] && std::isfinite(values[i])) {
            sum += values[i];
            ++result.count;
        }
    }
    if (result.count == 0) return result;
    result.mean = sum / static_cast<double>(result.count);
    if (result.count == 1) return result;
    double sumsq = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (valid[i] && std::isfinite(values[i])) {
            const double d = values[i] - result.mean;
            sumsq += d * d;
        }
    }
    result.se = std::sqrt(sumsq /
        (static_cast<double>(result.count - 1) * static_cast<double>(result.count)));
    return result;
}

static void write_metadata(std::ostream& out, int bc_id, int n, float T1, float Tn,
                           int total_trajectories, int valid_trajectories,
                           std::uint64_t step_samples) {
    out << "# bc_id=" << bc_id << "\n"
        << "# bc_label=" << bc_label(bc_id) << "\n"
        << "# " << bc_action_equation(bc_id) << "\n"
        << "# " << bc_phase_equation(bc_id) << "\n"
        << "# sigma_I=2*sqrt(2*gamma*T*I)\n"
        << "# sigma_phi=sqrt(2*gamma*T/I)\n"
        << "# n=" << n << " T1=" << T1 << " Tn=" << Tn
        << " gamma=" << gamma_val << " dt_max=" << dt
        << " adaptive_dt=min(dt_max,1/max_drift),floor=1e-5\n"
        << "# burnin=" << T_burnin << " measure=" << T_measure
        << " seed=" << base_seed << " requested_trajectories=" << total_trajectories
        << " valid_trajectories=" << valid_trajectories
        << " discarded_trajectories=" << total_trajectories - valid_trajectories
        << " integration_step_samples=" << step_samples << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 10) {
        std::cerr << "Usage: " << argv[0]
                  << " <bc_id:0..3> <T1> <Tn> <n> <batches> <burnin> <measure> <seed> <output_prefix>\n"
                  << "  0=BC1 canonical; 1=BC2 no phase drift;\n"
                  << "  2=BC3b no global term with phase drift;\n"
                  << "  3=BC3 no global term and no phase drift\n";
        return 1;
    }

    const int case_num = std::stoi(argv[1]);
    const float T1 = std::stof(argv[2]);
    const float Tn = std::stof(argv[3]);
    const int n = std::stoi(argv[4]);
    N_sample = std::stoi(argv[5]);
    T_burnin = std::stod(argv[6]);
    T_measure = std::stod(argv[7]);
    base_seed = static_cast<uint32_t>(std::stoul(argv[8]));
    const std::string prefix = argv[9];
    T_final = T_burnin + T_measure;
    if (case_num < 0 || case_num > 3 || n < 2 || N_sample < 1 ||
        !(T1 > 0.0f) || !(Tn > 0.0f) || !(T_burnin >= 0.0) || !(T_measure > 0.0)) {
        throw std::invalid_argument("invalid simulation parameters");
    }

    const int flux_mode = n / 2;
    const int total_trajectories = N_sample * LANES;
    std::cout << std::setprecision(17)
              << "=== Boundary-profile SIMD production ===\n"
              << "BC: " << bc_label(case_num) << " (id=" << case_num << ")\n"
              << bc_action_equation(case_num) << "\n"
              << bc_phase_equation(case_num) << "\n"
              << "sigma_I=2*sqrt(2*gamma*T*I)\n"
              << "sigma_phi=sqrt(2*gamma*T/I)\n"
              << "n=" << n << " T1=" << T1 << " Tn=" << Tn
              << " gamma=" << gamma_val << " dt_max=" << dt
              << " burnin=" << T_burnin << " measure=" << T_measure << "\n"
              << "batches=" << N_sample << " lanes=" << LANES
              << " trajectories=" << total_trajectories << " seed=" << base_seed
              << " threads=" << N_thread << " prefix=" << prefix << "\n";

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<unsigned char> valid(total_trajectories, 0);
    std::vector<double> flux_samples(total_trajectories, nan);
    std::vector<double> max_left(total_trajectories, nan);
    std::vector<double> max_right(total_trajectories, nan);
    std::vector<std::vector<double>> action_samples(
        n, std::vector<double>(total_trajectories, nan));
    std::vector<std::vector<double>> sin_samples(
        n - 1, std::vector<double>(total_trajectories, nan));
    std::array<std::vector<std::vector<double>>, N_CHECKPOINTS> burn_action;
    std::array<std::vector<std::vector<double>>, N_CHECKPOINTS> burn_sin;
    std::array<std::vector<std::vector<double>>, N_CHECKPOINTS> checkpoint_action;
    std::array<std::vector<std::vector<double>>, N_CHECKPOINTS> checkpoint_sin;
    for (int c = 0; c < N_CHECKPOINTS; ++c) {
        burn_action[c].assign(n, std::vector<double>(total_trajectories, nan));
        burn_sin[c].assign(n - 1, std::vector<double>(total_trajectories, nan));
        checkpoint_action[c].assign(n, std::vector<double>(total_trajectories, nan));
        checkpoint_sin[c].assign(n - 1, std::vector<double>(total_trajectories, nan));
    }
    std::vector<std::uint64_t> batch_steps(N_sample, 0);
    std::vector<std::uint64_t> batch_projections(N_sample, 0);

    const auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for schedule(dynamic) num_threads(N_thread)
    for (int i = 0; i < N_sample; ++i) {
        const int rank = omp_get_thread_num();
        RNGState rng = initialize_rng(rank, i);
        SimResult16 result;
        result.accumulated_energy.resize(n);
        result.accumulated_sin.resize(n - 1);
        for (int j = 0; j < n; ++j) result.accumulated_energy[j].setZero();
        for (int j = 0; j < n - 1; ++j) result.accumulated_sin[j].setZero();
        for (int c = 0; c < N_CHECKPOINTS; ++c) {
            result.burn_snapshot_energy[c].resize(n);
            result.burn_snapshot_sin[c].resize(n - 1);
            result.checkpoint_energy[c].resize(n);
            result.checkpoint_sin[c].resize(n - 1);
            for (int j = 0; j < n; ++j) result.burn_snapshot_energy[c][j].setZero();
            for (int j = 0; j < n - 1; ++j) result.burn_snapshot_sin[c][j].setZero();
            for (int j = 0; j < n; ++j) result.checkpoint_energy[c][j].setZero();
            for (int j = 0; j < n - 1; ++j) result.checkpoint_sin[c][j].setZero();
        }
        result.valid.fill(true);
        result.max_left_action.fill(0.0);
        result.max_right_action.fill(0.0);
        result.accumulated_flux.setZero();
        result.measurement_time = 0.0;
        result.measurement_steps = 0;
        result.projection_count = 0;

        run_simulation_16_old(n, case_num, T1, Tn, flux_mode, rng, result);
        batch_steps[i] = result.measurement_steps;
        batch_projections[i] = result.projection_count;
        const double inv_time = result.measurement_time > 0.0
            ? 1.0 / result.measurement_time : nan;
        for (int l = 0; l < LANES; ++l) {
            const int trajectory = i * LANES + l;
            bool lane_valid = result.valid[l] && std::isfinite(inv_time) &&
                              std::isfinite(result.accumulated_flux(l));
            for (int j = 0; lane_valid && j < n; ++j)
                lane_valid = std::isfinite(result.accumulated_energy[j](l));
            for (int j = 0; lane_valid && j < n - 1; ++j)
                lane_valid = std::isfinite(result.accumulated_sin[j](l));
            if (!lane_valid) continue;
            valid[trajectory] = 1;
            flux_samples[trajectory] = result.accumulated_flux(l) * inv_time;
            max_left[trajectory] = result.max_left_action[l];
            max_right[trajectory] = result.max_right_action[l];
            for (int j = 0; j < n; ++j)
                action_samples[j][trajectory] = result.accumulated_energy[j](l) * inv_time;
            for (int j = 0; j < n - 1; ++j)
                sin_samples[j][trajectory] = result.accumulated_sin[j](l) * inv_time;
            for (int c = 0; c < N_CHECKPOINTS; ++c) {
                for (int j = 0; j < n; ++j)
                    burn_action[c][j][trajectory] = result.burn_snapshot_energy[c][j](l);
                for (int j = 0; j < n - 1; ++j)
                    burn_sin[c][j][trajectory] = result.burn_snapshot_sin[c][j](l);
                for (int j = 0; j < n; ++j)
                    checkpoint_action[c][j][trajectory] = result.checkpoint_energy[c][j](l);
                for (int j = 0; j < n - 1; ++j)
                    checkpoint_sin[c][j][trajectory] = result.checkpoint_sin[c][j](l);
            }
        }
        #pragma omp critical
        std::cout << "completed_batch=" << i + 1 << "/" << N_sample << "\n";
    }

    const int valid_count = static_cast<int>(
        std::count(valid.begin(), valid.end(), static_cast<unsigned char>(1)));
    std::uint64_t step_samples = 0;
    std::uint64_t projection_count = 0;
    for (int i = 0; i < N_sample; ++i) {
        int valid_lanes = 0;
        for (int l = 0; l < LANES; ++l) valid_lanes += valid[i * LANES + l];
        step_samples += batch_steps[i] * static_cast<std::uint64_t>(valid_lanes);
        projection_count += batch_projections[i];
    }

    {
        std::ofstream out(prefix + "_profile.csv");
        write_metadata(out, case_num, n, T1, Tn, total_trajectories,
                       valid_count, step_samples);
        out << "j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta\n";
        out << std::setprecision(17);
        for (int j = 0; j < n; ++j) {
            const MeanSE a = mean_se(action_samples[j], valid);
            out << j + 1 << ',' << a.mean << ',' << a.se;
            if (j < n - 1) {
                const MeanSE s = mean_se(sin_samples[j], valid);
                out << ',' << s.mean << ',' << s.se;
            } else {
                out << ",,";
            }
            out << '\n';
        }
    }
    {
        std::ofstream out(prefix + "_burnin_checkpoints.csv");
        write_metadata(out, case_num, n, T1, Tn, total_trajectories,
                       valid_count, step_samples);
        out << "checkpoint_time,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta\n";
        out << std::setprecision(17);
        for (int c = 0; c < N_CHECKPOINTS; ++c) {
            const double checkpoint_time = T_burnin * (c + 1) / N_CHECKPOINTS;
            for (int j = 0; j < n; ++j) {
                const MeanSE a = mean_se(burn_action[c][j], valid);
                out << checkpoint_time << ',' << j + 1 << ',' << a.mean << ',' << a.se;
                if (j < n - 1) {
                    const MeanSE s = mean_se(burn_sin[c][j], valid);
                    out << ',' << s.mean << ',' << s.se;
                } else {
                    out << ",,";
                }
                out << '\n';
            }
        }
    }
    {
        std::ofstream out(prefix + "_checkpoints.csv");
        write_metadata(out, case_num, n, T1, Tn, total_trajectories,
                       valid_count, step_samples);
        out << "checkpoint_time,j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta\n";
        out << std::setprecision(17);
        for (int c = 0; c < N_CHECKPOINTS; ++c) {
            const double checkpoint_time = T_measure * (c + 1) / N_CHECKPOINTS;
            for (int j = 0; j < n; ++j) {
                const MeanSE a = mean_se(checkpoint_action[c][j], valid);
                out << checkpoint_time << ',' << j + 1 << ',' << a.mean << ',' << a.se;
                if (j < n - 1) {
                    const MeanSE s = mean_se(checkpoint_sin[c][j], valid);
                    out << ',' << s.mean << ',' << s.se;
                } else {
                    out << ",,";
                }
                out << '\n';
            }
        }
    }
    {
        std::ofstream out(prefix + "_trajectory_diagnostics.csv");
        write_metadata(out, case_num, n, T1, Tn, total_trajectories,
                       valid_count, step_samples);
        out << "trajectory_id,valid,time_averaged_flux,max_I_left,max_I_right\n";
        out << std::setprecision(17);
        for (int i = 0; i < total_trajectories; ++i)
            out << i << ',' << static_cast<int>(valid[i]) << ',' << flux_samples[i]
                << ',' << max_left[i] << ',' << max_right[i] << '\n';
    }

    const MeanSE flux = mean_se(flux_samples, valid);
    const auto stop = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(stop - start).count();
    {
        std::ofstream out(prefix + "_summary.csv");
        out << "bc_id,bc_label,n,T1,Tn,gamma,dt_max,burnin,measure,seed,requested_trajectories,"
               "valid_trajectories,nonfinite_trajectories,discarded_trajectories,"
               "integration_step_samples,projection_count,mean_flux,se_mean_flux,elapsed_seconds\n";
        out << std::setprecision(17)
            << case_num << ',' << bc_label(case_num) << ',' << n << ',' << T1 << ',' << Tn
            << ',' << gamma_val << ',' << dt << ',' << T_burnin << ',' << T_measure
            << ',' << base_seed << ',' << total_trajectories << ',' << valid_count
            << ',' << total_trajectories - valid_count
            << ',' << total_trajectories - valid_count << ',' << step_samples
            << ',' << projection_count << ',' << flux.mean << ',' << flux.se << ',' << elapsed << '\n';
    }

    std::cout << "valid_trajectories=" << valid_count
              << " nonfinite_trajectories=" << total_trajectories - valid_count
              << " discarded_trajectories=" << total_trajectories - valid_count
              << " projection_count=" << projection_count
              << " mean_flux=" << flux.mean << " se_flux=" << flux.se
              << " elapsed_seconds=" << elapsed << "\n"
              << "wrote_prefix=" << prefix << "\n";
    return valid_count == 0 ? 2 : 0;
}
