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


// Coarse grid 
const double I1_lo = 0.0,  I1_hi = 30.0;
const double I2_lo = 0.0,  I2_hi = 30.0;
const double I3_lo = 0.0,  I3_hi = 20.0;
const double Th_lo = -M_PI, Th_hi = M_PI;

const int G1  = 30;
const int G2  = 30;
const int G3  = 25;
const int Gt1 = 40;
const int Gt3 = 40;

const double dI1  = (I1_hi - I1_lo) / G1;
const double dI2  = (I2_hi - I2_lo) / G2;
const double dI3  = (I3_hi - I3_lo) / G3;
const double dTh1 = (Th_hi - Th_lo) / Gt1;
const double dTh3 = (Th_hi - Th_lo) / Gt3;

const long long total_cells = (long long)G1 * G2 * G3 * Gt1 * Gt3;

// KDE bandwidth
const double alpha_bw = 1.5;
const double h1  = alpha_bw * dI1;
const double h2  = alpha_bw * dI2;
const double h3  = alpha_bw * dI3;
const double ht1 = alpha_bw * dTh1;
const double ht3 = alpha_bw * dTh3;

const int R1  = (int)ceil(h1 / dI1);
const int R2  = (int)ceil(h2 / dI2);
const int R3  = (int)ceil(h3 / dI3);
const int Rt1 = (int)ceil(ht1 / dTh1);
const int Rt3 = (int)ceil(ht3 / dTh3);


// MC parameters
int       N_box    = 50000;
long long N_sample = 20000000;
int       N_thread = 8;


// State
struct State6 {
    double v[6];
    double& operator()(int i) { return v[i]; }
    const double& operator()(int i) const { return v[i]; }
};

struct Coords5D {
    double I1, I2, I3, th1, th3;
};

inline double wrap(double x) {
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}

// EM
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

Coords5D to_5d(const State6& X) {
    Coords5D c;
    c.I1  = X(0);
    c.I2  = X(1);
    c.I3  = X(2);
    c.th1 = wrap(2.0*(X(3) - X(4)));
    c.th3 = wrap(2.0*(X(5) - X(4)));
    return c;
}


// Grid index helpers
inline long long cell_index_int(int n1, int n2, int n3, int nt1, int nt3) {
    return ((long long)n1 * G2 + n2) * G3 * Gt1 * Gt3
           + (long long)n3 * Gt1 * Gt3
           + (long long)nt1 * Gt3
           + nt3;
}

inline void grid_coords(const Coords5D& c, int& n1, int& n2, int& n3, int& nt1, int& nt3) {
    n1  = max(0, min((int)floor((c.I1  - I1_lo) / dI1),  G1  - 1));
    n2  = max(0, min((int)floor((c.I2  - I2_lo) / dI2),  G2  - 1));
    n3  = max(0, min((int)floor((c.I3  - I3_lo) / dI3),  G3  - 1));
    nt1 = max(0, min((int)floor((c.th1 - Th_lo) / dTh1), Gt1 - 1));
    nt3 = max(0, min((int)floor((c.th3 - Th_lo) / dTh3), Gt3 - 1));
}

inline bool in_domain(const Coords5D& c) {
    return c.I1 >= I1_lo && c.I1 < I1_hi
        && c.I2 >= I2_lo && c.I2 < I2_hi
        && c.I3 >= I3_lo && c.I3 < I3_hi;
}


// kernel
inline double epanechnikov(double u) {
    if(u > 1.0 || u < -1.0) return 0.0;
    return 0.75 * (1.0 - u * u);
}

// ====================================================================
// Phase 1: Create collocation points + register into cell-list
//
// cell_list[cell_idx] = list of collocation point IDs in that cell
// ====================================================================
void create_collocation_points(
    vector<double>&         box_coords,
    vector<vector<int>>&    cell_list
) {
    // Track occupied cells for uniqueness
    vector<bool> cell_occupied(total_cells, false);

    State6 X;
    X(0) = 1.0; X(1) = 1.0; X(2) = 0.1;
    X(3) = 0.0; X(4) = 0.0; X(5) = 0.0;

    mt19937 rng(42);
    normal_distribution<double> dist(0.0, 1.0);

    cout << "  Burning in..." << endl;
    for(int i = 0; i < Burn_in; i++)
        step_EM(X, rng, dist);

    cout << "  Sampling collocation points..." << endl;
    int count = 0;
    while(count < N_box) {
        for(int i = 0; i < 5000; i++)
            step_EM(X, rng, dist);

        Coords5D c = to_5d(X);
        if(!in_domain(c)) continue;

        int n1, n2, n3, nt1, nt3;
        grid_coords(c, n1, n2, n3, nt1, nt3);
        long long idx = cell_index_int(n1, n2, n3, nt1, nt3);

        if(!cell_occupied[idx]) {
            cell_occupied[idx] = true;

            // Store exact coordinates
            box_coords[5*count + 0] = c.I1;
            box_coords[5*count + 1] = c.I2;
            box_coords[5*count + 2] = c.I3;
            box_coords[5*count + 3] = c.th1;
            box_coords[5*count + 4] = c.th3;

            // Register into cell-list
            cell_list[idx].push_back(count);

            count++;
            if(count % 10000 == 0)
                cout << "    " << count << "/" << N_box << " points" << endl;
        }
    }

    // Stats: how many cells have collocation points
    int occupied = 0;
    int max_per_cell = 0;
    for(long long i = 0; i < total_cells; i++) {
        if(!cell_list[i].empty()) {
            occupied++;
            max_per_cell = max(max_per_cell, (int)cell_list[i].size());
        }
    }
    cout << "  Occupied cells: " << occupied << " / " << total_cells << endl;
    cout << "  Max points per cell: " << max_per_cell << endl;
}

// ====================================================================
// Phase 2: MC with direct KDE in inner loop
// ====================================================================
void MC_direct_kde(
    const vector<double>&      box_coords,
    const vector<vector<int>>& cell_list,
    vector<double>&            density
) {
    cout << "MC direct KDE: " << N_thread << " threads, "
         << N_sample << " samples each" << endl;

    long long total_samples = (long long)N_thread * N_sample;
    double h_prod = h1 * h2 * h3 * ht1 * ht3;

    // Per-thread density accumulators to avoid atomic ops on double arrays
    vector<vector<double>> thread_density(N_thread, vector<double>(N_box, 0.0));

    long long total_in = 0;
    long long total_kernel_evals = 0;

#pragma omp parallel num_threads(N_thread) reduction(+:total_in,total_kernel_evals)
    {
        int rank = omp_get_thread_num();
        mt19937 rng(random_device{}() + rank*137);
        normal_distribution<double> dist(0.0, 1.0);

        State6 X;
        X(0) = 1.0 + 0.1*rank; X(1) = 1.0;    X(2) = 0.1;
        X(3) = 0.1*rank;       X(4) = 0.0;     X(5) = -0.1*rank;

        for(int i = 0; i < Burn_in; i++)
            step_EM(X, rng, dist);

        double* my_density = thread_density[rank].data();
        long long local_in = 0;
        long long local_evals = 0;

        for(long long step = 0; step < N_sample; step++) {
            step_EM(X, rng, dist);
            Coords5D c = to_5d(X);

            if(!in_domain(c)) continue;
            local_in++;

            int cn1, cn2, cn3, cnt1, cnt3;
            grid_coords(c, cn1, cn2, cn3, cnt1, cnt3);

            // Loop over neighbor cells
            for(int d1 = cn1 - R1; d1 <= cn1 + R1; d1++) {
                if(d1 < 0 || d1 >= G1) continue;
                double u1 = (c.I1 - (I1_lo + d1*dI1 + dI1*0.5)) / h1;

            for(int d2 = cn2 - R2; d2 <= cn2 + R2; d2++) {
                if(d2 < 0 || d2 >= G2) continue;

            for(int d3 = cn3 - R3; d3 <= cn3 + R3; d3++) {
                if(d3 < 0 || d3 >= G3) continue;

            for(int dt1 = cnt1 - Rt1; dt1 <= cnt1 + Rt1; dt1++) {
                int dt1w = ((dt1 % Gt1) + Gt1) % Gt1;

            for(int dt3 = cnt3 - Rt3; dt3 <= cnt3 + Rt3; dt3++) {
                int dt3w = ((dt3 % Gt3) + Gt3) % Gt3;

                long long cell_idx = cell_index_int(d1, d2, d3, dt1w, dt3w);
                const vector<int>& pts = cell_list[cell_idx];
                if(pts.empty()) continue;

                // For each collocation point in this cell
                for(int box_id : pts) {
                    // Exact distance from sample to collocation point
                    double xi_I1  = box_coords[5*box_id + 0];
                    double xi_I2  = box_coords[5*box_id + 1];
                    double xi_I3  = box_coords[5*box_id + 2];
                    double xi_th1 = box_coords[5*box_id + 3];
                    double xi_th3 = box_coords[5*box_id + 4];

                    double e1 = epanechnikov((c.I1 - xi_I1) / h1);
                    if(e1 == 0.0) continue;
                    double e2 = epanechnikov((c.I2 - xi_I2) / h2);
                    if(e2 == 0.0) continue;
                    double e3 = epanechnikov((c.I3 - xi_I3) / h3);
                    if(e3 == 0.0) continue;

                    double delta_t1 = c.th1 - xi_th1;
                    if(delta_t1 >  M_PI) delta_t1 -= 2.0*M_PI;
                    if(delta_t1 < -M_PI) delta_t1 += 2.0*M_PI;
                    double et1 = epanechnikov(delta_t1 / ht1);
                    if(et1 == 0.0) continue;

                    double delta_t3 = c.th3 - xi_th3;
                    if(delta_t3 >  M_PI) delta_t3 -= 2.0*M_PI;
                    if(delta_t3 < -M_PI) delta_t3 += 2.0*M_PI;
                    double et3 = epanechnikov(delta_t3 / ht3);
                    if(et3 == 0.0) continue;

                    double w = e1 * e2 * e3 * et1 * et3;
                    my_density[box_id] += w;
                    local_evals++;
                }
            }}}}}

            if(rank == 0 && step % 5000000LL == 0 && step > 0)
                cout << "  thread 0: " << step/1000000 << "M/"
                     << N_sample/1000000 << "M" << endl;
        }

        total_in += local_in;
        total_kernel_evals += local_evals;
    }

    // Merge thread results
    for(int t = 0; t < N_thread; t++)
        for(int i = 0; i < N_box; i++)
            density[i] += thread_density[t][i];

    // Normalize
    for(int i = 0; i < N_box; i++)
        density[i] /= ((double)total_samples * h_prod);

    cout << "Total in-domain samples: " << total_in << endl;
    cout << "Total kernel evaluations: " << total_kernel_evals << endl;
}

// Main
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc > 1) N_thread = atoi(argv[1]);
    if(argc > 2) gamma_val = atof(argv[2]);
    if(argc > 3) T1 = atof(argv[3]);
    if(argc > 4) T3 = atof(argv[4]);
    if(argc > 5) N_box = atoi(argv[5]);
    if(argc > 6) N_sample = atoll(argv[6]);

    cout << "=== NLS 5D: Direct KDE in inner loop (Prof's approach) ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "dt=" << dt << " Burn_in=" << Burn_in << endl;
    cout << "N_box=" << N_box << " N_sample/thread=" << N_sample << endl;
    cout << "N_thread=" << N_thread << endl;
    cout << endl;
    cout << "Coarse grid: " << G1 << "x" << G2 << "x" << G3
         << "x" << Gt1 << "x" << Gt3 << " = " << total_cells << " cells" << endl;
    cout << "KDE bandwidth: h=(" << h1 << "," << h2 << "," << h3
         << "," << ht1 << "," << ht3 << ")" << endl;
    cout << "Neighbor radius: R=(" << R1 << "," << R2 << "," << R3
         << "," << Rt1 << "," << Rt3 << ")" << endl;
    cout << "Per-thread density memory: "
         << N_box * 8 / (1024*1024) << " MB x " << N_thread << " threads" << endl;
    cout << endl;

    // ----------------------------------------------------------------
    // Phase 1: Create collocation points + cell-list
    // ----------------------------------------------------------------
    cout << "Allocating cell-list..." << endl;
    vector<vector<int>> cell_list(total_cells);
    vector<double> box_coords(5 * N_box);

    cout << "Creating collocation points..." << endl;
    create_collocation_points(box_coords, cell_list);
    cout << "Done." << endl << endl;

    // ----------------------------------------------------------------
    // Phase 2: MC with direct KDE
    // ----------------------------------------------------------------
    struct timeval mc_start, mc_end;
    gettimeofday(&mc_start, NULL);

    vector<double> density(N_box, 0.0);
    MC_direct_kde(box_coords, cell_list, density);

    gettimeofday(&mc_end, NULL);
    double mc_sec = ((mc_end.tv_sec - mc_start.tv_sec)*1e6
                    + mc_end.tv_usec - mc_start.tv_usec) / 1e6;
    cout << "MC direct KDE time: " << mc_sec << "s" << endl;

    // ----------------------------------------------------------------
    // Output
    // ----------------------------------------------------------------
    ofstream f1("NLS_FP_boxes_direct.txt");
    ofstream f2("NLS_FP_density_direct.txt");
    for(int i = 0; i < N_box; i++) {
        f1 << box_coords[5*i+0] << " " << box_coords[5*i+1] << " "
           << box_coords[5*i+2] << " " << box_coords[5*i+3] << " "
           << box_coords[5*i+4] << endl;
        f2 << density[i] << endl;
    }
    f1.close();
    f2.close();

    // Summary
    double mx = *max_element(density.begin(), density.end());
    double mn = *min_element(density.begin(), density.end());
    int nz = 0;
    for(int i = 0; i < N_box; i++)
        if(density[i] > 0) nz++;

    cout << endl;
    cout << "=== Results ===" << endl;
    cout << "Max density = " << mx << endl;
    cout << "Min density = " << mn << endl;
    cout << "Nonzero density = " << nz << "/" << N_box << endl;

    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec - t1.tv_sec)*1e6 + t2.tv_usec - t1.tv_usec) / 1e6;
    cout << "Total wall time = " << sec << "s" << endl;

    return 0;
}