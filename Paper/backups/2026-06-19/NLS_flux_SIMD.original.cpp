//SIMD-Optimized Hamiltonian Energy Flux Simulation
// hamiltonian_flux_simd.cpp
// SIMD-optimized simulation of mean energy flux in a Hamiltonian system
// with action-angle coordinates coupled to two heat baths.
//
// Each thread processes 16 independent trajectories simultaneously using
// Eigen::Array<float,16,1> which maps to AVX-512 (or 2x AVX2) registers.
//
// Compile with:
//   g++ -O3 -mavx512f -mfma -fopenmp -I/path/to/eigen hamiltonian_flux_simd.cpp -o flux_sim
// or for AVX2:
//   g++ -O3 -mavx2 -mfma -fopenmp -I/path/to/eigen hamiltonian_flux_simd.cpp -o flux_sim

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>
#include <random>

// ============================================================
// Global constants
// ============================================================
static const float gamma_val = 0.1f;
static double T_final = 1200.0;     // task: burn 1000 + measure 200 (override via argv[6])
static double T_burnin = 1000.0;    // task: burn-in 1000
static const double dt = 0.0005f;
static const float PI_f = 3.14159265358979323846f;
static const float PI_HALF_f = 1.5707963268f;

static int N_sample = 625;          // 625*16 = 10^4 trajectories (override via argv[5])
static const int N_thread = 10;
static const int LANES = 16;      // SIMD width: 16 floats = 512 bits

typedef Eigen::Array<float,  LANES, 1> A16f;
typedef Eigen::Array<double,  LANES, 1> A16d;
typedef std::vector<A16f, Eigen::aligned_allocator<A16f>> AlignedVec;

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
    // Deterministic seeding from (rank, sample_idx) via splitmix64
    RNGState state;
    std::random_device rd;
    uint32_t seed = rd();
    uint64_t z = static_cast<uint64_t>(rank * 100000 + sample_idx + seed) + 0x9E3779B97F4A7C15ULL;
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
    AlignedVec accumulated_energy;  // size n, each A16f
    A16d accumulated_flux;
    double measurement_time;
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
        }

        // 6) Generate stochastic noise (Box-Muller: 4 independent normal arrays)
        A16f nI0, nIn, nPhi0, nPhin;
        box_muller_16(rng, nI0, nIn);
        box_muller_16(rng, nPhi0, nPhin);


        // Compute noise amplitudes BEFORE updating state (use current I values)
        A16f noise_I_0 = 2.0f * (gamma_val * T1 * I[0].max(1e-10f)).sqrt() * nI0;
        A16f noise_I_n = 2.0f * (gamma_val * Tn * I[n-1].max(1e-10f)).sqrt() * nIn;
        A16f noise_phi_0 = (gamma_val * T1 / I[0].max(1e-10f)).sqrt() * nPhi0;
        A16f noise_phi_n = (gamma_val * Tn / I[n-1].max(1e-10f)).sqrt() * nPhin;

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

        // 8) Clamp action variables to positive
        for (int j = 0; j < n; j++) {
            I[j] = I[j].max(1e-10f);
            phi[j] = wrap_pi(phi[j]);  // <-- THIS IS THE CRITICAL FIX
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

        // 9) Accumulate measurements after burn-in
        if (!measuring && current_time >= T_burnin) {
            measuring = true;
        }

        if (measuring) {
            for (int j = 0; j < n; j++) {
                result.accumulated_energy[j] += I[j] * ddt;
            }
            result.accumulated_flux += flux_j.cast<double>() * ddt;
            result.measurement_time += ddt;
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
// Main
// ============================================================

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <case_number> <T1> <Tn> <n>" << std::endl;
        std::cerr << "  case_number: 0 (canonical BC), 1 (I-only drift), 2 (no global term)" << std::endl;
        std::cerr << "Example: " << argv[0] << " 0 10 2 50" << std::endl;
        return 1;
    }

    int case_num = std::atoi(argv[1]);
    float T1 = static_cast<float>(std::atof(argv[2]));
    float Tn = static_cast<float>(std::atof(argv[3]));
    int n = std::atoi(argv[4]);
    if (argc >= 6) N_sample = std::atoi(argv[5]);   // batches; *16 = #trajectories
    if (argc >= 7) T_final  = std::atof(argv[6]);   // optional: override measuring end
    float deltaT = T1 - Tn;

    int flux_mode = n / 2;
    int total_trajectories = N_sample * LANES;

    std::cout << "=== SIMD-Optimized Hamiltonian Flux Simulation ===" << std::endl;
    std::cout << "Case " << case_num << " | n=" << n << " | T1=" << T1 << " | Tn=" << Tn << std::endl;
    std::cout << "deltaT = " << deltaT << std::endl;
    std::cout << "Burn-in: " << T_burnin << "s | Measurement: " << (T_final - T_burnin) << "s" << std::endl;
    std::cout << "Flux measured at mode " << flux_mode << std::endl;
    std::cout << "Threads: " << N_thread << " | SIMD lanes: " << LANES << std::endl;
    std::cout << "Batches: " << N_sample << " | Total trajectories: " << total_trajectories << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // Storage for per-trajectory flux values (for statistics)
    std::vector<double> flux_samples(total_trajectories, 0.0f);

    // Accumulator for mean energy profile (double for precision)
    std::vector<double> total_e(n, 0.0);

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        std::vector<double> local_e(n, 0.0);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N_sample; ++i) {
            // Initialize RNG with unique seed per (rank, sample)
            RNGState rng = initialize_rng(rank, i);

            // Prepare result structure
            SimResult16 result;
            result.accumulated_energy.resize(n);
            for (int j = 0; j < n; j++) result.accumulated_energy[j].setZero();
            result.accumulated_flux.setZero();
            result.measurement_time = 0.0f;

            // Run 16 trajectories in parallel via SIMD
            run_simulation_16_old(n, case_num, T1, Tn, flux_mode, rng, result);

            // Extract per-lane results
            double inv_mtime = 1.0 / result.measurement_time;
            for (int l = 0; l < LANES; l++) {
                double mf = result.accumulated_flux(l) * inv_mtime;
                flux_samples[i * LANES + l] = mf;  // unique index, no race
            }

            // Accumulate energy (local to avoid contention)
            for (int j = 0; j < n; j++) {
                for (int l = 0; l < LANES; l++) {
                    local_e[j] += static_cast<double>(result.accumulated_energy[j](l) * inv_mtime);
                }
            }

            if (rank == 0) {
                std::cout << "  Completed batch " << (i + 1) << "/" << N_sample
                          << " (" << (i + 1) * LANES << " trajectories so far)" << std::endl;
            }
        }



        // Reduce local accumulators
        #pragma omp critical
        {
            for (int j = 0; j < n; j++) total_e[j] += local_e[j];
        }
    }

    // ---- Compute statistics ----
    double inv_total = 1.0 / static_cast<double>(total_trajectories);

    // Mean energy per oscillator
    for (int j = 0; j < n; j++) total_e[j] *= inv_total;
    /*
    std::cout<< " global energy profile printed "<<std::endl;
    std::ofstream flux_file("energy_profile_summary.csv");
    for(int j = 0; j < n; j++)
    {
        flux_file<<total_e[j]<<" ";
    }
    flux_file<<std::endl;
*/
    // Mean flux
    double mean_flux = 0.0;
    for (int i = 0; i < total_trajectories; i++) {
        mean_flux += static_cast<double>(flux_samples[i]);
        //std::cout<<" local flux = "<<flux_samples[i]<<std::endl;
    }
    mean_flux *= inv_total;

    // Variance and confidence interval
    double variance = 0.0;
    for (int i = 0; i < total_trajectories; i++) {
        double diff = static_cast<double>(flux_samples[i]) - mean_flux;
        variance += diff * diff;
    }
    variance /= (total_trajectories - 1);
    double std_dev = sqrt(variance);
    double ci_95 = 1.96 * std_dev;

    // ---- Output results ----
    std::cout << "\n====== Results ======" << std::endl;
    std::cout << "Case: " << case_num << std::endl;
    std::cout << "n = " << n << ", T1 = " << T1 << ", Tn = " << Tn << std::endl;
    std::cout << "deltaT = " << deltaT << std::endl;
    std::cout << "Total trajectories: " << total_trajectories << std::endl;
    std::cout << "Flux measured at mode " << flux_mode << std::endl;
    std::cout << "Mean Flux = " << mean_flux << std::endl;
    std::cout << "Std Dev = " << std_dev << std::endl;
    std::cout << "95% CI = [" << (mean_flux - ci_95) << ", " << (mean_flux + ci_95) << "]" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // ---- Write ALL per-trajectory flux samples (for distribution / tail analysis) ----
    {
        std::string fn = "flux_n" + std::to_string(n) + ".txt";
        std::ofstream fs(fn.c_str());
        fs << "# per-trajectory time-averaged flux | n=" << n
           << " T1=" << T1 << " Tn=" << Tn
           << " burnin=" << T_burnin << " Tfinal=" << T_final
           << " window=" << (T_final - T_burnin)
           << " n_traj=" << total_trajectories << "\n";
        for (int i = 0; i < total_trajectories; i++) fs << flux_samples[i] << "\n";
        fs.close();
        std::cout << "Wrote " << fn << " (" << total_trajectories << " flux samples)\n";
    }
    // ---- Write CSV summary ----
    std::ofstream flux_file("flux_vs_length_2.csv", std::ios::app);
    flux_file.seekp(0, std::ios::end);
    if (flux_file.tellp() == 0) {
        flux_file << "case,n,T1,Tn,deltaT,flux_mode,n_traj,mean_flux,std_dev,ci_lower,ci_upper,covers_zero\n";
    }
    flux_file << case_num << ","<< n << "," << T1 << "," << Tn << "," << deltaT << ","
              << flux_mode << "," << total_trajectories << ","
              << mean_flux << "," << std_dev << ","
              << (mean_flux - ci_95) << "," << (mean_flux + ci_95) << ","
              << ((mean_flux - ci_95 <= 0 && mean_flux + ci_95 >= 0) ? 1 : 0) << "\n";
    flux_file.close();
    

    std::cout << "Results appended to flux_vs_length.csv" << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;

    return 0;
}


