// Dumps (I_j, phi_j) trajectories for a range of modes for LTE analysis.
// Fixed to Case 0 (closed chain).
//
// Usage:
//   ./lte_dump <T1> <Tn> <n> <j_lo> <j_hi> <output.csv>
// Example:
//   ./lte_dump 10 2 25 5 20 traj_N25.csv

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>

const double gamma_val = 0.1;
const double T_final   = 30000.0;
const double T_burnin  = 2000.0;
const double dt        = 0.001;
const int    dump_every_steps = 100;
const int    N_traj    = 64;
const int    N_thread  = 8;
const double PI        = 3.14159265358979323846;

using Vector = Eigen::VectorXd;

void run_one_trajectory(int n, double T1, double Tn,
                        int j_lo, int j_hi, int traj_id,
                        std::mt19937& rng,
                        std::ostream& out)
{
    std::normal_distribution<double> dist(0.0, 1.0);

    Vector I   = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);

    double current_time = 0.0;
    long step = 0;
    bool measuring = false;

    while (current_time < T_final) {
        Vector I_padded = Vector::Zero(n + 2);
        I_padded.segment(1, n) = I;
        Vector I_prev = I_padded.head(n);
        Vector I_next = I_padded.tail(n);

        Vector phi_padded = Vector::Zero(n + 2);
        phi_padded.segment(1, n) = phi;
        Vector phi_prev = phi_padded.head(n);
        Vector phi_next = phi_padded.tail(n);

        Vector d_phi_prev = 2.0 * (phi - phi_prev);
        Vector d_phi_next = 2.0 * (phi - phi_next);

        Vector drift_I = 4.0 * I.array() *
            (I_prev.array() * d_phi_prev.array().sin() +
             I_next.array() * d_phi_next.array().sin());

        double total_mass_I = I.sum();

        Vector drift_phi = (Vector::Constant(n, 2.0 * total_mass_I) - I).array()
            + 2.0 * I_prev.array() * d_phi_prev.array().cos()
            + 2.0 * I_next.array() * d_phi_next.array().cos();

        // --- Case 0 boundary terms (FIXED: phi drift sign +=) ---
        drift_I(0)   += 2.0 * gamma_val * (2.0 * T1 -
            (2.0 * total_mass_I * I(0) - pow(I(0), 2)
             + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
        drift_phi(0) += gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));   // FIX

        drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn -
            (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2)
             + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
        drift_phi(n-1) += gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));   // FIX

        // --- noise (FIXED: sqrt(2) factor) ---
        double sqrt_dt = sqrt(dt);
        Vector noise_I   = Vector::Zero(n);
        Vector noise_phi = Vector::Zero(n);
        noise_I(0)     = 2.0 * sqrt(2.0 * gamma_val * T1 * I(0))     * dist(rng);   // FIX
        noise_I(n-1)   = 2.0 * sqrt(2.0 * gamma_val * Tn * I(n-1))   * dist(rng);   // FIX
        noise_phi(0)   = sqrt(2.0 * gamma_val * T1 / I(0))           * dist(rng);   // FIX
        noise_phi(n-1) = sqrt(2.0 * gamma_val * Tn / I(n-1))         * dist(rng);   // FIX

        I   += drift_I   * dt + noise_I   * sqrt_dt;
        phi += drift_phi * dt + noise_phi * sqrt_dt;
        I = I.cwiseMax(1e-10);

        current_time += dt;
        step++;

        if (!measuring && current_time >= T_burnin) measuring = true;

        if (measuring && (step % dump_every_steps == 0)) {
            std::ostringstream oss;
            oss << traj_id << "," << current_time;
            for (int j = j_lo; j <= j_hi; ++j) oss << "," << I(j);
            for (int j = j_lo; j <= j_hi; ++j) {
                double p = phi(j);
                p = std::fmod(p + PI, 2.0 * PI);
                if (p < 0) p += 2.0 * PI;
                p -= PI;
                oss << "," << p;
            }
            oss << "\n";
            #pragma omp critical
            {
                out << oss.str();
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0]
                  << " <T1> <Tn> <n> <j_lo> <j_hi> <output.csv>\n"
                  << "Example: " << argv[0] << " 10 2 25 5 20 traj_N25.csv\n";
        return 1;
    }
    double T1   = std::atof(argv[1]);
    double Tn   = std::atof(argv[2]);
    int    n    = std::atoi(argv[3]);
    int    j_lo = std::atoi(argv[4]);
    int    j_hi = std::atoi(argv[5]);
    std::string outfile = argv[6];

    if (j_lo < 0 || j_hi >= n || j_lo > j_hi) {
        std::cerr << "Bad mode range: need 0 <= j_lo <= j_hi < n\n";
        return 1;
    }

    std::cout << "LTE trajectory dump (Case 0 closed chain) [FIXED]\n"
              << "T1=" << T1 << " Tn=" << Tn << " n=" << n
              << " modes=[" << j_lo << "," << j_hi << "]\n"
              << "T_final=" << T_final << " T_burnin=" << T_burnin
              << " dt=" << dt << " dump_every=" << dump_every_steps
              << " (=" << dump_every_steps * dt << " time units)\n"
              << "N_traj=" << N_traj << "\n"
              << "Expected samples per trajectory: "
              << (long)((T_final - T_burnin) / (dump_every_steps * dt))
              << "\n"
              << "Total expected samples: "
              << N_traj * (long)((T_final - T_burnin) / (dump_every_steps * dt))
              << "\n";

    std::ofstream out(outfile);
    if (!out) { std::cerr << "Cannot open " << outfile << "\n"; return 1; }

    out << "traj_id,time";
    for (int j = j_lo; j <= j_hi; ++j) out << ",I_" << j;
    for (int j = j_lo; j <= j_hi; ++j) out << ",phi_" << j;
    out << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        std::seed_seq seed{
            static_cast<unsigned>(rank),
            static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<unsigned>(clock())
        };
        std::mt19937 thread_rng(seed);

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N_traj; ++t) {
            run_one_trajectory(n, T1, Tn, j_lo, j_hi, t, thread_rng, out);
            if (rank == 0) {
                #pragma omp critical
                std::cout << "  trajectory " << t << " (thread "
                          << rank << ") done\n";
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
    std::cout << "Done in " << sec << "s. Output: " << outfile << "\n";
    return 0;
}
