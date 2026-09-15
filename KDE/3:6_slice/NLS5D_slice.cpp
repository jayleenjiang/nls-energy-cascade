#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <random>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cstring>

using namespace std;

// Physical parameters
double gamma_val = 0.1;
double T1        = 10.0;
double T3        = 2.0;

// Integrator parameters
const double dt      = 0.0001;
const double sqrt_dt = sqrt(dt);
const int    Burn_in = 2000000;

// Slice parameters: fix I1 = I2 = I3 = I_slice, with tolerance
double I_slice = 2.0;       // the slice value
double I_tol   = 0.5;       // accept samples with |I_k - I_slice| < I_tol


// 2D grid for (theta1, theta3) in [-pi, pi]^2
const double Th_lo = -M_PI, Th_hi = M_PI;

// Test 1: 50x50 simple binning
const int G_bin = 50;
const double d_bin = (Th_hi - Th_lo) / G_bin;  // ~0.1257

// Test 2: 30x30 KDE grid
const int G_kde = 30;
const double d_kde = (Th_hi - Th_lo) / G_kde;  // ~0.2094

// MC parameters
long long N_sample = 20000000;
int       N_thread = 8;

// State and helpers (same as NLS5D2.cpp)
struct State6 {
    double v[6];
    double& operator()(int i) { return v[i]; }
    const double& operator()(int i) const { return v[i]; }
};

inline double wrap(double x) {
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}

void step_EM(State6& X, mt19937& rng, normal_distribution<double>& dist)
{
    double I1 = X(0), I2 = X(1), I3 = X(2);
    double p1 = X(3), p2 = X(4), p3 = X(5);

    double I_prev[3], I_next[3];
    I_prev[0] = 0.0;   I_next[0] = I2;
    I_prev[1] = I1;    I_next[1] = I3;
    I_prev[2] = I2;    I_next[2] = 0.0;

    double d_phi_prev[3], d_phi_next[3];
    d_phi_prev[0] = 2.0*(p1 - 0.0);
    d_phi_next[0] = 2.0*(p1 - p2);
    d_phi_prev[1] = 2.0*(p2 - p1);
    d_phi_next[1] = 2.0*(p2 - p3);
    d_phi_prev[2] = 2.0*(p3 - p2);
    d_phi_next[2] = 2.0*(p3 - 0.0);

    double Ivec[3] = {I1, I2, I3};
    double M = I1 + I2 + I3;

    double dI[3], dphi[3];
    for(int j = 0; j < 3; j++) {
        dI[j] = 4.0 * Ivec[j] * (I_prev[j]*sin(d_phi_prev[j])
                                 + I_next[j]*sin(d_phi_next[j]));
        dphi[j] = 2.0*M - Ivec[j]
                  + 2.0*I_prev[j]*cos(d_phi_prev[j])
                  + 2.0*I_next[j]*cos(d_phi_next[j]);
    }

    dI[0] += 2.0*gamma_val*(2.0*T1 - (2.0*M*I1 - I1*I1 + 2.0*I_next[0]*I1*cos(d_phi_next[0])));
    dphi[0] -= gamma_val*(2.0*I_next[0]*sin(d_phi_next[0]));

    dI[2] += 2.0*gamma_val*(2.0*T3 - (2.0*M*I3 - I3*I3 + 2.0*I_prev[2]*I3*cos(d_phi_prev[2])));
    dphi[2] -= gamma_val*(2.0*I_prev[2]*sin(d_phi_prev[2]));

    double dW[4];
    for(int j = 0; j < 4; j++) dW[j] = sqrt_dt * dist(rng);

    double diff_I1   = 2.0 * sqrt(2.0 * gamma_val * T1 * max(I1, 1e-14));
    double diff_I3   = 2.0 * sqrt(2.0 * gamma_val * max(T3, 1e-10) * max(I3, 1e-14));
    double diff_phi1 = sqrt(2.0 * gamma_val * T1 / max(I1, 1e-14));
    double diff_phi3 = sqrt(2.0 * gamma_val * max(T3, 1e-10) / max(I3, 1e-14));

    X(0) = I1  + dI[0]*dt + diff_I1*dW[0];
    X(1) = I2  + dI[1]*dt;
    X(2) = I3  + dI[2]*dt + diff_I3*dW[1];
    X(3) = p1  + dphi[0]*dt + diff_phi1*dW[2];
    X(4) = p2  + dphi[1]*dt;
    X(5) = p3  + dphi[2]*dt + diff_phi3*dW[3];

    X(0) = max(X(0), 1e-14);
    X(1) = max(X(1), 1e-14);
    X(2) = max(X(2), 1e-14);
    X(3) = wrap(X(3));
    X(4) = wrap(X(4));
    X(5) = wrap(X(5));
}

// Check if sample is in the slice
inline bool in_slice(double I1, double I2, double I3) {
    return fabs(I1 - I_slice) < I_tol
        && fabs(I2 - I_slice) < I_tol
        && fabs(I3 - I_slice) < I_tol;
}

// 1D Gaussian PDF
inline double normpdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * M_PI));
}

// Periodic distance for angles
inline double periodic_dist(double a, double b) {
    double d = a - b;
    if(d >  M_PI) d -= 2.0 * M_PI;
    if(d < -M_PI) d += 2.0 * M_PI;
    return d;
}

// Periodic Gaussian PDF for angles
inline double normpdf_periodic(double x, double mu, double h) {
    double d = periodic_dist(x, mu);
    double z = d / h;
    return exp(-0.5 * z * z) / (h * sqrt(2.0 * M_PI));
}


// Main
int main(int argc, char* argv[]) {
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);

    if(argc > 1) N_thread = atoi(argv[1]);
    if(argc > 2) gamma_val = atof(argv[2]);
    if(argc > 3) T1 = atof(argv[3]);
    if(argc > 4) T3 = atof(argv[4]);
    if(argc > 5) N_sample = atoll(argv[5]);
    if(argc > 6) I_slice = atof(argv[6]);
    if(argc > 7) I_tol = atof(argv[7]);

    cout << "=== NLS 5D: 2D Slice Test ===" << endl;
    cout << "Slice: I1=I2=I3=" << I_slice << " +/- " << I_tol << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << endl;
    cout << endl;

    // ================================================================
    // Phase 1: MC sampling — collect (theta1, theta3) samples in slice
    // ================================================================
    cout << "Collecting slice samples..." << endl;

    // All threads collect into shared vector
    vector<double> all_th1, all_th3;
    long long total_in_domain = 0;
    long long total_in_slice = 0;

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        mt19937 rng(random_device{}() + rank*137);
        normal_distribution<double> dist(0.0, 1.0);

        State6 X;
        X(0) = 1.0 + 0.1*rank; X(1) = 1.0;    X(2) = 0.1;
        X(3) = 0.1*rank;       X(4) = 0.0;     X(5) = -0.1*rank;

        for(int i = 0; i < Burn_in; i++)
            step_EM(X, rng, dist);

        // Thread-local collection
        vector<double> local_th1, local_th3;
        long long local_domain = 0, local_slice = 0;

        for(long long step = 0; step < N_sample; step++) {
            step_EM(X, rng, dist);

            double I1 = X(0), I2 = X(1), I3 = X(2);
            // Check domain
            if(I1 < 0 || I1 >= 30.0 || I2 < 0 || I2 >= 30.0 || I3 < 0 || I3 >= 20.0)
                continue;
            local_domain++;

            // Check slice
            if(in_slice(I1, I2, I3)) {
                double th1 = wrap(2.0*(X(3) - X(4)));
                double th3 = wrap(2.0*(X(5) - X(4)));
                local_th1.push_back(th1);
                local_th3.push_back(th3);
                local_slice++;
            }

            if(rank == 0 && step % 5000000LL == 0 && step > 0)
                cout << "  thread 0: " << step/1000000 << "M/"
                     << N_sample/1000000 << "M" << endl;
        }

        // Merge into shared vector
        #pragma omp critical
        {
            all_th1.insert(all_th1.end(), local_th1.begin(), local_th1.end());
            all_th3.insert(all_th3.end(), local_th3.begin(), local_th3.end());
            total_in_domain += local_domain;
            total_in_slice += local_slice;
        }
    }

    long long N_slice = all_th1.size();
    cout << "Total in-domain: " << total_in_domain << endl;
    cout << "Total in slice: " << N_slice
         << " (" << 100.0*N_slice/total_in_domain << "% of in-domain)" << endl;
    cout << endl;

    if(N_slice < 100) {
        cout << "ERROR: Too few slice samples. Increase N_sample or I_tol." << endl;
        return 1;
    }

    // ================================================================
    // Test 1: Simple binning on 50x50 grid
    // ================================================================
    cout << "=== Test 1: Simple binning 50x50 ===" << endl;
    cout << "Grid spacing: " << d_bin << " rad" << endl;

    vector<double> hist_50(G_bin * G_bin, 0.0);

    for(long long s = 0; s < N_slice; s++) {
        int j1 = (int)floor((all_th1[s] - Th_lo) / d_bin);
        int j3 = (int)floor((all_th3[s] - Th_lo) / d_bin);
        j1 = max(0, min(j1, G_bin - 1));
        j3 = max(0, min(j3, G_bin - 1));
        hist_50[j1 * G_bin + j3] += 1.0;
    }

    // Normalize to density: count / (N * cell_area)
    double cell_area_50 = d_bin * d_bin;
    {
        ofstream fout("slice_test1_binning50.txt");
        fout << "# Test 1: Simple binning 50x50, N_slice=" << N_slice << endl;
        fout << "# theta1_center  theta3_center  density" << endl;
        for(int j1 = 0; j1 < G_bin; j1++) {
            double th1_c = Th_lo + (j1 + 0.5) * d_bin;
            for(int j3 = 0; j3 < G_bin; j3++) {
                double th3_c = Th_lo + (j3 + 0.5) * d_bin;
                double dens = hist_50[j1 * G_bin + j3] / ((double)N_slice * cell_area_50);
                fout << th1_c << " " << th3_c << " " << dens << endl;
            }
        }
        fout.close();
    }

    // Stats
    double max_dens1 = 0, min_dens1 = 1e99;
    int nz1 = 0;
    for(int i = 0; i < G_bin*G_bin; i++) {
        double d = hist_50[i] / ((double)N_slice * cell_area_50);
        if(d > max_dens1) max_dens1 = d;
        if(d < min_dens1) min_dens1 = d;
        if(hist_50[i] > 0) nz1++;
    }
    cout << "  Nonzero cells: " << nz1 << "/" << G_bin*G_bin << endl;
    cout << "  Density range: [" << min_dens1 << ", " << max_dens1 << "]" << endl;
    cout << "  Avg counts per cell: " << (double)N_slice / (G_bin*G_bin) << endl;
    cout << "  Output: slice_test1_binning50.txt" << endl;
    cout << endl;

    // ================================================================
    // Test 2: KDE on 30x30 grid
    // For each grid point, compute:
    //   f(th1, th3) = (1/N) sum_i normpdf(th1, th1_i, h) * normpdf(th3, th3_i, h)
    //
    // Two bandwidths: h=0.05, h=0.10
    // Two sample sets: (a) local bin only, (b) all slice samples
    // ================================================================

    // First: bin all slice samples onto 30x30 grid for local-bin lookup
    cout << "=== Test 2: KDE on 30x30 grid ===" << endl;
    cout << "Grid spacing: " << d_kde << " rad" << endl;

    // Build bin lists: for each 30x30 cell, store indices of samples in it
    vector<vector<long long>> bin_lists(G_kde * G_kde);
    for(long long s = 0; s < N_slice; s++) {
        int j1 = (int)floor((all_th1[s] - Th_lo) / d_kde);
        int j3 = (int)floor((all_th3[s] - Th_lo) / d_kde);
        j1 = max(0, min(j1, G_kde - 1));
        j3 = max(0, min(j3, G_kde - 1));
        bin_lists[j1 * G_kde + j3].push_back(s);
    }

    // Stats on bin occupancy
    long long max_bin = 0, min_bin = N_slice, empty_bins = 0;
    for(int i = 0; i < G_kde*G_kde; i++) {
        long long sz = bin_lists[i].size();
        if(sz > max_bin) max_bin = sz;
        if(sz < min_bin) min_bin = sz;
        if(sz == 0) empty_bins++;
    }
    cout << "  Avg samples per bin: " << (double)N_slice / (G_kde*G_kde) << endl;
    cout << "  Min/Max per bin: " << min_bin << " / " << max_bin << endl;
    cout << "  Empty bins: " << empty_bins << "/" << G_kde*G_kde << endl;
    cout << endl;

    double h_vals[2] = {0.05, 0.10};

    for(int ih = 0; ih < 2; ih++) {
        double h = h_vals[ih];
        cout << "  --- h = " << h << " ---" << endl;

        // ============================================================
        // (a) KDE with local bin samples only
        // ============================================================
        {
            char fname[256];
            snprintf(fname, sizeof(fname), "slice_test2_kde_local_h%.2f.txt", h);
            ofstream fout(fname);
            fout << "# Test 2a: KDE local-bin only, 30x30, h=" << h
                 << ", N_slice=" << N_slice << endl;
            fout << "# theta1_center  theta3_center  density  n_local" << endl;

            int zero_count = 0;
            for(int j1 = 0; j1 < G_kde; j1++) {
                double th1_c = Th_lo + (j1 + 0.5) * d_kde;
                for(int j3 = 0; j3 < G_kde; j3++) {
                    double th3_c = Th_lo + (j3 + 0.5) * d_kde;

                    const vector<long long>& local_samples = bin_lists[j1 * G_kde + j3];
                    int n_local = local_samples.size();

                    double sum = 0.0;
                    for(int s = 0; s < n_local; s++) {
                        long long idx = local_samples[s];
                        sum += normpdf_periodic(th1_c, all_th1[idx], h)
                             * normpdf_periodic(th3_c, all_th3[idx], h);
                    }

                    // Normalize by total N_slice (not n_local) so density integrates correctly
                    double dens = sum / (double)N_slice;
                    fout << th1_c << " " << th3_c << " " << dens << " " << n_local << endl;
                    if(n_local == 0) zero_count++;
                }
            }
            fout.close();
            cout << "    Local-bin KDE: " << fname << " (empty bins: " << zero_count << ")" << endl;
        }

        // ============================================================
        // (b) KDE with ALL slice samples (global)
        // ============================================================
        {
            char fname[256];
            snprintf(fname, sizeof(fname), "slice_test2_kde_global_h%.2f.txt", h);
            ofstream fout(fname);
            fout << "# Test 2b: KDE global (all samples), 30x30, h=" << h
                 << ", N_slice=" << N_slice << endl;
            fout << "# theta1_center  theta3_center  density" << endl;

            // This is O(G_kde^2 * N_slice) — can be slow if N_slice is large
            // Parallelize over grid points
            vector<double> kde_global(G_kde * G_kde, 0.0);

            #pragma omp parallel for num_threads(N_thread) schedule(dynamic, 4)
            for(int j1 = 0; j1 < G_kde; j1++) {
                double th1_c = Th_lo + (j1 + 0.5) * d_kde;
                for(int j3 = 0; j3 < G_kde; j3++) {
                    double th3_c = Th_lo + (j3 + 0.5) * d_kde;

                    double sum = 0.0;
                    for(long long s = 0; s < N_slice; s++) {
                        sum += normpdf_periodic(th1_c, all_th1[s], h)
                             * normpdf_periodic(th3_c, all_th3[s], h);
                    }

                    kde_global[j1 * G_kde + j3] = sum / (double)N_slice;
                }
            }

            for(int j1 = 0; j1 < G_kde; j1++) {
                double th1_c = Th_lo + (j1 + 0.5) * d_kde;
                for(int j3 = 0; j3 < G_kde; j3++) {
                    double th3_c = Th_lo + (j3 + 0.5) * d_kde;
                    fout << th1_c << " " << th3_c << " " << kde_global[j1*G_kde+j3] << endl;
                }
            }
            fout.close();
            cout << "    Global KDE:    " << fname << endl;

            // Stats
            double mx = *max_element(kde_global.begin(), kde_global.end());
            double mn = *min_element(kde_global.begin(), kde_global.end());
            double avg = 0;
            for(auto v : kde_global) avg += v;
            avg /= kde_global.size();
            cout << "    Density range: [" << mn << ", " << mx << "], mean=" << avg << endl;
        }

        cout << endl;
    }

    // Summary
    gettimeofday(&t_end, NULL);
    double sec = ((t_end.tv_sec - t_start.tv_sec)*1e6
                 + t_end.tv_usec - t_start.tv_usec) / 1e6;

    cout << "=== Output files ===" << endl;
    cout << "  slice_test1_binning50.txt           -- Test 1: simple binning 50x50" << endl;
    cout << "  slice_test2_kde_local_h0.05.txt     -- Test 2a: KDE local, h=0.05" << endl;
    cout << "  slice_test2_kde_local_h0.10.txt     -- Test 2a: KDE local, h=0.10" << endl;
    cout << "  slice_test2_kde_global_h0.05.txt    -- Test 2b: KDE global, h=0.05" << endl;
    cout << "  slice_test2_kde_global_h0.10.txt    -- Test 2b: KDE global, h=0.10" << endl;
    cout << endl;
    cout << "Total wall time = " << sec << "s" << endl;

    return 0;
}
