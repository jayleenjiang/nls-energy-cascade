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
const double T1 = 10.0;
const double T_final = 2000.0;
const double dt = 0.001;
const double PI = 3.14159265358979323846;

const int N_sample = 16;
const int N_thread = 8;

using Vector = Eigen::VectorXd;

struct SimulationResult {
    Vector accumulated_energy;
    double accumulated_flux;
    double total_time;
};

void run_simulation(int n, int case_num, double Tn, int flux_mode, std::mt19937& rng, SimulationResult& result);

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <case_number> <Tn> <n>" << std::endl;
        std::cerr << "\nCases 0-3: Closed chain" << std::endl;
        std::cerr << "Cases 4-7: Open chain" << std::endl;
        return 1;
    }
    
    int case_num = std::atoi(argv[1]);
    double Tn = std::atof(argv[2]);
    int n = std::atoi(argv[3]);
    
    std::cout << "Running Case " << case_num << " with Tn = " << Tn << ", n = " << n << std::endl;
    std::cout << "T1 = " << T1 << ", T1-Tn = " << (T1 - Tn) << std::endl;
    std::cout << "Chain type: " << (case_num <= 3 ? "CLOSED" : "OPEN") << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // Choose middle mode for flux measurement
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
            result.total_time = 0.0;
            
            run_simulation(n, case_num, Tn, flux_mode, thread_rng, result);
            
            local_e += result.accumulated_energy;
            double mean_flux = result.accumulated_flux / result.total_time;
            
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
    
    // Calculate standard deviation and confidence interval
    double variance = 0.0;
    for (int i = 0; i < N_sample; ++i) {
        variance += (flux_samples[i] - mean_flux) * (flux_samples[i] - mean_flux);
    }
    variance /= (N_sample - 1);
    double std_dev = sqrt(variance);
    double std_error = std_dev / sqrt(N_sample);
    double ci_95 = 1.96 * std_error;

    // Output results
    std::cout << "Case: " << case_num << " (" << (case_num <= 3 ? "Closed" : "Open") << " chain)" << std::endl;
    std::cout << "Chain length n = " << n << std::endl;
    std::cout << "T1 = " << T1 << ", Tn = " << Tn << std::endl;
    std::cout << "Temperature difference (T1 - Tn) = " << (T1 - Tn) << std::endl;
    std::cout << "Flux measured at mode " << flux_mode << std::endl;
    std::cout << "Mean Flux = " << mean_flux << std::endl;
    std::cout << "Std Dev = " << std_dev << std::endl;
    std::cout << "95% CI = [" << (mean_flux - ci_95) << ", " << (mean_flux + ci_95) << "]" << std::endl;
    std::cout << "CI covers zero: " << ((mean_flux - ci_95 <= 0 && mean_flux + ci_95 >= 0) ? "YES" : "NO") << std::endl;
    std::cout << "-------------------------------------\n" << std::endl;

    // Save energy distribution
    std::ofstream outfile;
    std::string filename = "energy_C" + std::to_string(case_num) + "_n" + std::to_string(n) + "_Tn" + std::to_string(Tn) + ".csv";
    outfile.open(filename);
    outfile << "mode,energy\n";
    for (int i = 0; i < n; ++i) {
        outfile << i + 1 << "," << total_e(i) << "\n";
    }
    outfile.close();
    
    std::ofstream flux_file;
    flux_file.open("flux_summary.csv", std::ios::app);
    flux_file.seekp(0, std::ios::end);
    if (flux_file.tellp() == 0) {
        flux_file << "case,chain_type,n,T1,Tn,deltaT,flux_mode,mean_flux,std_dev,ci_lower,ci_upper,covers_zero\n";
    }
    flux_file << case_num << "," << (case_num <= 3 ? "closed" : "open") << "," 
              << n << "," << T1 << "," << Tn << "," << (T1 - Tn) << ","
              << flux_mode << "," << mean_flux << "," << std_dev << ","
              << (mean_flux - ci_95) << "," << (mean_flux + ci_95) << ","
              << ((mean_flux - ci_95 <= 0 && mean_flux + ci_95 >= 0) ? 1 : 0) << "\n";
    flux_file.close();
    
    std::cout << "Results appended to flux_summary.csv" << std::endl;
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;
    
    return 0;
}

void run_simulation(int n, int case_num, double Tn, int flux_mode, std::mt19937& rng, SimulationResult& result) {
    std::normal_distribution<double> dist(0.0, 1.0);

    Vector I = Vector::Constant(n, 0.1);
    I(0) = 1.0;
    Vector phi = Vector::Zero(n);
    double current_time = 0;

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

        // flux: J_j = 4 * I_{j-1} * I_j * sin(2(phi_j - phi_{j-1}))
        double flux_j = 4.0 * I_prev(flux_mode) * I(flux_mode) * sin(d_phi_prev(flux_mode));

        switch(case_num) {
            case 0:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * total_mass_I * I(n-1) - pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                break;
                
            case 1:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * I(0) * I_next(0) + pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - (2.0 * I(n-1) * I_prev(n-1) + pow(I(n-1), 2) + 2.0 * I_prev(n-1) * I(n-1) * cos(d_phi_prev(n-1))));
                drift_phi(n-1) -= gamma_val * (2.0 * I_prev(n-1) * sin(d_phi_prev(n-1)));
                break;
                
            case 2:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - pow(I(0), 2));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - pow(I(n-1), 2));
                break;
                
            case 3:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - pow(I(0), 2));
                drift_I(n-1) += 2.0 * gamma_val * (2.0 * Tn - pow(I(n-1), 2));
                break;
                
            case 4:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * I(0) * I_next(0) + pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                break;
                
            case 5:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - pow(I(0), 2));
                break;
                
            case 6:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - pow(I(0), 2));
                break;
                
            case 7:
                drift_I(0) += 2.0 * gamma_val * (2.0 * T1 - (2.0 * total_mass_I * I(0) - pow(I(0), 2) + 2.0 * I_next(0) * I(0) * cos(d_phi_next(0))));
                drift_phi(0) -= gamma_val * (2.0 * I_next(0) * sin(d_phi_next(0)));
                break;
        }

        // Adaptive timestep
        double dt_drift_0 = (drift_I(0) < 0) ? -I(0)/(4.0*drift_I(0)) : dt;
        double dt_drift_n = (drift_I(n-1) < 0) ? -I(n-1)/(4.0*drift_I(n-1)) : dt;
        double dt_diff_0 = I(0) / (72.0 * gamma_val * T1);
        double dt_diff_n = I(n-1) / (72.0 * gamma_val * std::max(Tn, 0.001));

        double ddt = dt; 
        ddt = std::min(ddt, dt_drift_0);
        ddt = std::min(ddt, dt_drift_n);
        ddt = std::min(ddt, dt_diff_0);
        ddt = std::min(ddt, dt_diff_n);

        ddt = std::max(ddt, 1e-12);

        Vector dW(4);
        for(int j=0; j<4; ++j) dW(j) = sqrt(ddt) * dist(rng);
    
        double diff_I1 = 2.0 * sqrt(2.0 * gamma_val * T1 * I(0));
        double diff_In = 2.0 * sqrt(2.0 * gamma_val * std::max(Tn, 1e-10) * I(n - 1));
        double diff_phi1, diff_phin;
        
        if (case_num == 3 || case_num == 6) {
            diff_phi1 = sqrt(2.0 * gamma_val * T1);
            diff_phin = sqrt(2.0 * gamma_val * std::max(Tn, 1e-10));
        } else {
            diff_phi1 = sqrt(2.0 * gamma_val * T1 / I(0));
            diff_phin = sqrt(2.0 * gamma_val * std::max(Tn, 1e-10) / (I(n-1) + 1e-14));
        }
 
        double milstein_I1 = 2.0 * gamma_val * T1 * (pow(dW(0), 2) - ddt);
        double milstein_In = 2.0 * gamma_val * Tn * (pow(dW(1), 2) - ddt);

        Vector I_new = I + drift_I * ddt;
        Vector phi_new = phi + drift_phi * ddt;
        
        if (case_num >= 0 && case_num <= 7) {
            I_new(0) += diff_I1 * dW(0) + milstein_I1;
            phi_new(0) += diff_phi1 * dW(2);
        }
        
        if (case_num >= 0 && case_num <= 3) {
            I_new(n-1) += diff_In * dW(1) + milstein_In;
            phi_new(n-1) += diff_phin * dW(3);
        }

        current_time += ddt;

        result.accumulated_energy += 0.5*(I + I_new)*ddt;
        
        result.accumulated_flux += flux_j * ddt;
        result.total_time = current_time;
        
        I = I_new;
        phi = phi_new;
        
        I(0) = std::max(I(0), 1e-14);
        I(n - 1) = std::max(I(n - 1), 1e-14);
        
        for(int j=0; j<n; ++j) {
            phi(j) = fmod(phi(j) + PI, 2.0 * PI) - PI;
        }
    }
}