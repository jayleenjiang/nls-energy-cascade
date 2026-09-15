#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>

// global constants
const double gamma = 0.1;
const double T1 = 10.0;
const double Tn = 0.01;
const double T_final = 2000.0;
const double dt = 0.001;
const double PI = 3.14159265358979323846;

const int N_sample = 16; 
const int N_thread = 8;  

using Vector = Eigen::VectorXd;
void run_simulation(int n, std::mt19937& rng, Vector& accumulated_energy);

int main() {
    std::vector<int> n_values = {200};  

    auto start = std::chrono::high_resolution_clock::now();

    for (int n : n_values) {
        std::cout << "Running simulation for n = " << n << "..." << std::endl;
        
        Vector total_e = Vector::Zero(n); 

        #pragma omp parallel num_threads(N_thread)
        {
            int rank = omp_get_thread_num();
            // Create a seed sequence using multiple sources of entropy
            std::seed_seq seed{
                static_cast<unsigned>(rank),
                static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                static_cast<unsigned>(clock())
            };
            std::mt19937 thread_rng(seed);

            // Thread-local accumulator
            Vector local_e = Vector::Zero(n);

            #pragma omp for
            for (int i = 0; i < N_sample; ++i) {
                run_simulation(n, thread_rng, local_e);
                
                if (rank == 0 && (i + 1) % 4 == 0) {
                    std::cout << "  Progress: " << (i + 1) * N_thread << "/" << N_sample << " samples" << std::endl;
                }
            }

            // Combine thread-local results
            #pragma omp critical
            {
                total_e += local_e;
            }
        } 

        total_e /= N_sample;

        std::ofstream outfile;
        outfile.open("energy_7n" + std::to_string(n) + ".csv");
        outfile << "mode,energy\n";
        for (int i = 0; i < n; ++i) {
            outfile << i + 1 << "," << total_e(i) << "\n";
        }
        outfile.close();
        std::cout << "Finished n = " << n << ". Output written to energy_n" << n << ".csv" << std::endl;
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "\nTotal execution time: " << duration.count() << " seconds." << std::endl;
    
    return 0;
}

// accumulate energy on the fly
void run_simulation(int n, std::mt19937& rng, Vector& accumulated_energy) {
    const int N_steps = static_cast<int>(T_final / dt);
    std::normal_distribution<double> dist(0.0, 1.0);

    Vector I = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);
    
    for (int i = 0; i < N_steps; ++i) {
        // Accumulate energy at this timestep
        accumulated_energy += I;
        
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

        Vector drift_I = 4.0 * I.array() * (I_prev.array() * d_phi_prev.array().sin() + I_next.array() * d_phi_next.array().sin());
        
        double total_mass_I = I.sum();
        
        Vector drift_phi = (Vector::Constant(n, 2 * total_mass_I) - I).array() + 2 * I_prev.array() * d_phi_prev.array().cos() + 2 * I_next.array() * d_phi_next.array().cos();

        // Original
        drift_I(0) += 2.0 * gamma * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
        drift_phi(0) += gamma * (2.0 * I_next(0) * sin(d_phi_next(0)));
        
        //drift_I(n - 1) += 2.0 * gamma * (2.0 * Tn - (2.0 * total_mass_I * I(n - 1) - pow(I(n - 1), 2) + 2.0 * I_prev(n - 1) * I(n - 1) * cos(d_phi_prev(n - 1))));
        //drift_phi(n - 1) += gamma * (2.0 * I_prev(n - 1) * sin(d_phi_prev(n - 1)));
        
        Vector dW(4);
        for(int j=0; j<4; ++j) dW(j) = sqrt(dt) * dist(rng);
        
        double diff_I1 = 2.0 * sqrt(2.0 * gamma * T1 * I(0));
        double diff_In = 2.0 * sqrt(2.0 * gamma * Tn * I(n - 1));
        double diff_phi1 = sqrt(2.0 * gamma * T1 / I(0));
        double diff_phin = sqrt(2.0 * gamma * Tn / I(n - 1));

        double milstein_I1 = 2.0 * gamma * T1 * (pow(dW(0), 2) - dt);
        double milstein_In = 2.0 * gamma * Tn * (pow(dW(1), 2) - dt);

        Vector I_new = I + drift_I * dt;
        Vector phi_new = phi + drift_phi * dt;
        
        I_new(0) += diff_I1 * dW(0) + milstein_I1;
        //I_new(n - 1) += diff_In * dW(1) + milstein_In;
        
        phi_new(0) += diff_phi1 * dW(2);
        //phi_new(n - 1) += diff_phin * dW(3);
        
        I = I_new;
        phi = phi_new;
        
        I(0) = std::max(I(0), 1e-9);
        I(n - 1) = std::max(I(n - 1), 1e-9);
        for(int j=0; j<n; ++j) {
            phi(j) = fmod(phi(j) + PI, 2.0 * PI) - PI;
        }
    }
}