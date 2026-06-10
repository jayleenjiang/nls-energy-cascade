#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>

// global constants
const double gamma_val = 0.1;
const double T_final = 10000.0;
const double T_burnin = 2000.0;
const double dt = 0.001;
const double PI = 3.14159265358979323846;

const int N_sample = 32;
const int N_thread = 8;

using Vector = Eigen::VectorXd;

struct SimulationResult {
    Vector accumulated_energy;
    double accumulated_flux;
    double measurement_time;
};

void run_simulation(int n, int case_num, double T1, double Tn, int flux_mode, std::mt19937& rng, SimulationResult& result);

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <case_number> <T1> <Tn> <n>" << std::endl;
        std::cerr << "\nExample: " << argv[0] << " 0 10 2 50" << std::endl;
        std::cerr << "\nCases 0-3: Closed chain" << std::endl;
        std::cerr << "Cases 4-7: Open chain" << std::endl;
        return 1;
    }
    
    int case_num = std::atoi(argv[1]);
    double T1 = std::atof(argv[2]);
    double Tn = std::atof(argv[3]);
    int n = std::atoi(argv[4]);
    
    double deltaT = T1 - Tn;
    
    std::cout << "Running Case " << case_num << " with T1 = " << T1 << ", Tn = " << Tn << ", n = " << n << std::endl;
    std::cout << "deltaT = " << deltaT << std::endl;
    std::cout << "Chain type: " << (case_num <= 3 ? "CLOSED" : "OPEN") << std::endl;
    std::cout << "Burn in time: " << T_burnin << ", Measurement time: " << (T_final - T_burnin) << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    int flux_mode = n / 2;
    
    std::cout << "Measuring flux at mode " << flux_mode << std::endl;
    
    Vector total_e = Vector::Zero(n);
    double total_flux = 0.0;
    std::vector<double> flux_samples(N_sample);

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        std::seed_seq seed{
            static_cast<unsigned>(rank),
            static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<unsigned>(clock())
        };
        std::mt19937 thread_rng(seed);

        Vector local_e = Vector::Zero(n);
        double local_flux = 0.0;

        #pragma omp for
        for (int i = 0; i < N_sample; ++i) {
            SimulationResult result;
            result.accumulated_energy = Vector::Zero(n);
            result.accumulated_flux = 0.0;
            result.measurement_time = 0.0;
            
            run_simulation(n, case_num, T1, Tn, flux_mode, thread_rng, result);
            
            local_e += result.accumulated_energy;
            double mean_flux = result.accumulated_flux / result.measurement_time;
            
            #pragma omp critical
            {
                flux_samples[i] = mean_flux;
            }
            local_flux += mean_flux;
            
            if (rank == 0 && (i + 1) % 4 == 0) {
                std::cout << "  Progress: " << (i + 1) * N_thread << "/" << N_sample << " samples" << std::endl;
            }
        }

        #pragma omp critical
        {
            total_e += local_e;
            total_flux += local_flux;
        }
    } 

    total_e /= N_sample;
    double mean_flux = total_flux / N_sample;
    
    double variance = 0.0;
    for (int i = 0; i < N_sample; ++i) {
        variance += (flux_samples[i] - mean_flux) * (flux_samples[i] - mean_flux);
    }
    variance /= (N_sample - 1);
    double std_dev = sqrt(variance);
    double std_error = std_dev / sqrt(N_sample);
    double ci_95 = 1.96 * std_error;

    std::cout << "Case: " << case_num << " (" << (case_num <= 3 ? "Closed" : "Open") << " chain)" << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "T1 = " << T1 << ", Tn = " << Tn << std::endl;
    std::cout << "deltaT = " << deltaT << std::endl;
    std::cout << "Flux measured at mode " << flux_mode << std::endl;
    std::cout << "Mean Flux = " << mean_flux << std::endl;
    std::cout << "Std Dev = " << std_dev << std::endl;
    std::cout << "95% CI = [" << (mean_flux - ci_95) << ", " << (mean_flux + ci_95) << "]" << std::endl;
    std::cout << "-------------------------------------\n" << std::endl;

    // summary file with T1 and Tn in filename
    std::ofstream flux_file;
    flux_file.open("flux_summary_temp.csv", std::ios::app);
    flux_file.seekp(0, std::ios::end);
    if (flux_file.tellp() == 0) {
        flux_file << "case,chain_type,n,T1,Tn,deltaT,flux_mode,mean_flux,std_dev,ci_lower,ci_upper,covers_zero\n";
    }
    flux_file << case_num << "," << (case_num <= 3 ? "closed" : "open") << "," 
              << n << "," << T1 << "," << Tn << "," << deltaT << ","
              << flux_mode << "," << mean_flux << "," << std_dev << ","
              << (mean_flux - ci_95) << "," << (mean_flux + ci_95) << ","
              << ((mean_flux - ci_95 <= 0 && mean_flux + ci_95 >= 0) ? 1 : 0) << "\n";
    flux_file.close();
    
    std::cout << "Results appended to flux_summary_temp.csv" << std::endl;
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;
    
    return 0;
}

void run_simulation(int n, int case_num, double T1, double Tn, int flux_mode, std::mt19937& rng, SimulationResult& result) {
    std::normal_distribution<double> dist(0.0, 1.0);

    Vector I = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);
    double current_time = 0;
    
    bool measuring = false;

    while(current_time < T_final)
    {
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
        
        Vector drift_phi = (Vector::Constant(n, 2.0 * total_mass_I) - I).array() + 2.0 * I_prev.array() * d_phi_prev.array().cos() + 2.0 * I_next.array() * d_phi_next.array().cos();

        double flux_j = 4.0 * I_prev(flux_mode) * I(flux_mode) * sin(d_phi_prev(flux_mode));

        switch(case_num) {
            case 0:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                break;
            case 1:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                break;
            case 2:
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                break;
            case 3:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                break;
            case 4:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2)));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                break;
            case 5:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2)));
                break;
            case 6:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2)));
                break;
            case 7:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2)));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                break;
        }

        double sqrt_dt = sqrt(dt);
        Vector noise_I = Vector::Zero(n);
        Vector noise_phi = Vector::Zero(n);

        noise_I(0) = 2.0 * sqrt(gamma_val * T1 * I(0)) * dist(rng);
        noise_I(n-1) = 2.0 * sqrt(gamma_val * Tn * I(n-1)) * dist(rng);
        noise_phi(0) = sqrt(gamma_val * T1 / I(0)) * dist(rng);
        noise_phi(n-1) = sqrt(gamma_val * Tn / I(n-1)) * dist(rng);

        I += drift_I * dt + noise_I * sqrt_dt;
        phi += drift_phi * dt + noise_phi * sqrt_dt;

        I = I.cwiseMax(1e-10);

        current_time += dt;

        if (!measuring && current_time >= T_burnin) {
            measuring = true;
        }
        
        if (measuring) {
            result.accumulated_energy += I * dt;
            result.accumulated_flux += flux_j * dt;
            result.measurement_time += dt;
        }
    }
}