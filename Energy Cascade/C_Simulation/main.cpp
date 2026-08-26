#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>
#include <cmath> 
#include <iomanip> 


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
        outfile.open("energy" + std::to_string(n) + "_no_summation.csv");
        outfile << "mode,energy\n";
        for (int i = 0; i < n; ++i) {
            outfile << i + 1 << "," << total_e(i) << "\n";
        }
        outfile.close();
        std::cout << "Finished n = " << n << ". Output written to energy_n" << n << "_no_summation.csv" << std::endl;
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "\nTotal execution time: " << duration.count() << " seconds." << std::endl;
    
    return 0;
}\


// accumulate energy on the fly
void run_simulation(int n, std::mt19937& rng, Vector& accumulated_energy) {
    const int N_steps = static_cast<int>(T_final / dt);
    std::normal_distribution<double> dist(0.0, 1.0);

    Vector I = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);
    

    // new part 
    // set precision for debug 
    std::cout << std::fixed << std::setprecision(12);

    for (int i = 0; i < N_steps; ++i) {
        bool has_blowup = false;
        for (int j = 0; j < n; ++j) {
            if (!std::isfinite(I(j)) || !std::isfinite(phi(j))) {
                has_blowup = true;
                break;
            }
        }
        
        if (has_blowup) {
            #pragma omp critical
            {
                std::cout << "Thread: " << omp_get_thread_num() << std::endl;
                std::cout << "Step: " << i << " (Time: " << i * dt << ")" << std::endl;
                std::cout << "First few I values: ";
                for (int j = 0; j < std::min(5, n); ++j) {
                    std::cout << "I(" << j << ")=" << I(j) << " ";
                }
                std::cout << std::endl;
            }
            // fill accumulated_energy with NaN to mark this sample as failed
            accumulated_energy.fill(std::numeric_limits<double>::quiet_NaN());
            return; 
        }
        // end of the new part




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

        Vector drift_I = 2.0 * I.array() * (I_prev.array() * d_phi_prev.array().sin() + I_next.array() * d_phi_next.array().sin());
        
        double total_mass_I = I.sum();
        
        Vector drift_phi = (Vector::Constant(n, total_mass_I) - 0.5 * I).array() + I_prev.array() * d_phi_prev.array().cos() + I_next.array() * d_phi_next.array().cos();
        

        // remove summation term 
        drift_I(0) += 2.0 * gamma * (2.0 * T1 - (- pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
        drift_phi(0) -= gamma * (2.0 * I_next(0) * sin(d_phi_next(0)));
        
        drift_I(n - 1) += 2.0 * gamma * (2.0 * Tn - (- pow(I(n - 1), 2) + 2.0 * I_prev(n - 1) * I(n - 1) * cos(d_phi_prev(n - 1))));
        drift_phi(n - 1) -= gamma * (2.0 * I_prev(n - 1) * sin(d_phi_prev(n - 1)));
        
        Vector dW(4);
        for(int j=0; j<4; ++j) dW(j) = sqrt(dt) * dist(rng);
        
        double diff_I1 = 2.0 * sqrt(2.0 * gamma * T1 * I(0));
        double diff_In = 2.0 * sqrt(2.0 * gamma * Tn * I(n - 1));
        double diff_phi1 = sqrt(2.0 * gamma * T1 / I(0));
        double diff_phin = sqrt(2.0 * gamma * Tn / I(n - 1));

        double milstein_I1 = 2.0 * gamma * T1 * (pow(dW(0), 2) - dt);
        double milstein_In = 2.0 * gamma * Tn * (pow(dW(1), 2) - dt);

        // new values 
        Vector I_new = I + drift_I * dt;
        Vector phi_new = phi + drift_phi * dt;
        
        // NEW
        // store stochastic terms before adding them
        double I_new_0_stochastic = diff_I1 * dW(0);
        double I_new_0_milstein = milstein_I1;
        double I_new_n_stochastic = diff_In * dW(1);
        double phi_new_0_stochastic = diff_phi1 * dW(2);
        double phi_new_n_stochastic = diff_phin * dW(3);

        // add them to the new values
        I_new(0) += I_new_0_stochastic + I_new_0_milstein;
        I_new(n - 1) += I_new_n_stochastic;
        phi_new(0) += phi_new_0_stochastic;
        phi_new(n - 1) += phi_new_n_stochastic;
        
        // check for blow up in the updated values
        if (!std::isfinite(I_new(0)) || I_new(0) < 0) {
            #pragma omp critical
            {
                std::cout << "Thread: " << omp_get_thread_num() << std::endl;
                std::cout << "Step: " << i << " (Time: " << i * dt << ")" << std::endl;
                std::cout << "  PREVIOUS I(0):        " << I(0) << std::endl;
                std::cout << "  PREVIOUS I(1):        " << I(1) << std::endl;
                std::cout << "  PREVIOUS phi(0):      " << phi(0) << std::endl;
                std::cout << "--- Components of I_new(0) ---" << std::endl;
                std::cout << "  I(0) (start):         " << I(0) << std::endl;
                std::cout << "  drift_I(0):           " << drift_I(0) << std::endl;
                std::cout << "  drift_I(0) * dt:      " << drift_I(0) * dt << std::endl;
                std::cout << "  stochastic term:      " << I_new_0_stochastic << std::endl;
                std::cout << "  milstein term:        " << I_new_0_milstein << std::endl;
                std::cout << "  Sum before stochastic:" << (I(0) + drift_I(0) * dt) << std::endl;
                std::cout << "  CALCULATED I_new(0):  " << I_new(0) << std::endl;
                std::cout << "--- Stochastic Term ---" << std::endl;
                std::cout << "  sqrt(I(0)):           " << sqrt(I(0)) << std::endl;
                std::cout << "  diff_I1:              " << diff_I1 << std::endl;
                std::cout << "  dW(0):                " << dW(0) << std::endl;
            }
            // fill with NaN and exit
            accumulated_energy.fill(std::numeric_limits<double>::quiet_NaN());
            return;
        }

        // Update state
        I = I_new;
        phi = phi_new;
        
        // safety bounds
        I(0) = std::max(I(0), 1e-9);
        I(n - 1) = std::max(I(n - 1), 1e-9);
        
        // wrap phases
        for(int j=0; j<n; ++j) {
            phi(j) = fmod(phi(j) + PI, 2.0 * PI) - PI;
        }
    }
}