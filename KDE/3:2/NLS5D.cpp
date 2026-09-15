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

// Action variables
const double I1_lo = 0.0,  I1_hi = 30.0;
const double I2_lo = 0.0,  I2_hi = 30.0;
const double I3_lo = 0.0,  I3_hi = 20.0;
// Angle variables
const double Th_lo = -M_PI, Th_hi = M_PI;

const int G1  = 30;   // cells for I1
const int G2  = 30;   // cells for I2
const int G3  = 25;   // cells for I3
const int Gt1 = 40;   // cells for theta1
const int Gt3 = 40;   // cells for theta3

const double dI1  = (I1_hi - I1_lo) / G1;    // 1.0
const double dI2  = (I2_hi - I2_lo) / G2;    // 1.0
const double dI3  = (I3_hi - I3_lo) / G3;    // 0.8
const double dTh1 = (Th_hi - Th_lo) / Gt1;   // ~0.157
const double dTh3 = (Th_hi - Th_lo) / Gt3;   // ~0.157

// Total histogram cells
const long long total_cells = (long long)G1 * G2 * G3 * Gt1 * Gt3;  // 36M

// KDE bandwidth per dimension: h_d = alpha * delta_d
const double alpha_bw = 1.5;
const double h1  = alpha_bw * dI1;
const double h2  = alpha_bw * dI2;
const double h3  = alpha_bw * dI3;
const double ht1 = alpha_bw * dTh1;
const double ht3 = alpha_bw * dTh3;

// Neighbor search radius per dimension: R_d = ceil(h_d / delta_d)
const int R1  = (int)ceil(h1 / dI1);    // 2
const int R2  = (int)ceil(h2 / dI2);    // 2
const int R3  = (int)ceil(h3 / dI3);    // 2
const int Rt1 = (int)ceil(ht1 / dTh1);  // 2
const int Rt3 = (int)ceil(ht3 / dTh3);  // 2

// MC parameters
int       N_box    = 50000;
long long N_sample = 20000000;
int       N_thread = 8;

// 6D state: (I1, I2, I3, phi1, phi2, phi3)
struct State6 {
    double v[6];
    double& operator()(int i) { return v[i]; }
    const double& operator()(int i) const { return v[i]; }
};

struct Coords5D {
    double I1, I2, I3, th1, th3;
};

// Wrap angle to [-pi, pi]
inline double wrap(double x) {
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}

// EM integrator (unchanged from original)
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

    // Heat bath on mode 1
    dI[0] += 2.0*gamma_val*(2.0*T1 - (2.0*M*I1 - I1*I1 + 2.0*I_next[0]*I1*cos(d_phi_next[0])));
    dphi[0] -= gamma_val*(2.0*I_next[0]*sin(d_phi_next[0]));

    // Heat bath on mode 3
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

    // Reflecting boundary for actions
    X(0) = max(X(0), 1e-14);
    X(1) = max(X(1), 1e-14);
    X(2) = max(X(2), 1e-14);

    // Periodic boundary for angles
    X(3) = wrap(X(3));
    X(4) = wrap(X(4));
    X(5) = wrap(X(5));
}

// 6D -> 5D projection
Coords5D to_5d(const State6& X) {
    Coords5D c;
    c.I1  = X(0);
    c.I2  = X(1);
    c.I3  = X(2);
    c.th1 = wrap(2.0*(X(3) - X(4)));
    c.th3 = wrap(2.0*(X(5) - X(4)));
    return c;
}


// Histogram cell index: 5D -> 1D flat index
// Returns -1 if out of bounds
inline long long cell_index(const Coords5D& c) {
    if(c.I1 < I1_lo || c.I1 >= I1_hi) return -1;
    if(c.I2 < I2_lo || c.I2 >= I2_hi) return -1;
    if(c.I3 < I3_lo || c.I3 >= I3_hi) return -1;

    int n1  = (int)floor((c.I1 - I1_lo) / dI1);
    int n2  = (int)floor((c.I2 - I2_lo) / dI2);
    int n3  = (int)floor((c.I3 - I3_lo) / dI3);
    int nt1 = (int)floor((c.th1 - Th_lo) / dTh1);
    int nt3 = (int)floor((c.th3 - Th_lo) / dTh3);

    // Clamp 
    n1  = max(0, min(n1, G1 - 1));
    n2  = max(0, min(n2, G2 - 1));
    n3  = max(0, min(n3, G3 - 1));
    nt1 = max(0, min(nt1, Gt1 - 1));
    nt3 = max(0, min(nt3, Gt3 - 1));

    return ((long long)n1 * G2 + n2) * G3 * Gt1 * Gt3
           + (long long)n3 * Gt1 * Gt3
           + (long long)nt1 * Gt3
           + nt3;
}


// Cell index from integer grid coordinates (no bounds check)
inline long long cell_index_int(int n1, int n2, int n3, int nt1, int nt3) {
    return ((long long)n1 * G2 + n2) * G3 * Gt1 * Gt3
           + (long long)n3 * Gt1 * Gt3
           + (long long)nt1 * Gt3
           + nt3;
}


// Phase 1: Create collocation points along trajectory
void create_collocation_points(
    vector<double>& box_coords,   // 5 * N_box
    vector<bool>&   cell_occupied  // total_cells, to ensure uniqueness
) {
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
        long long idx = cell_index(c);
        if(idx < 0) continue;

        if(!cell_occupied[idx]) {
            cell_occupied[idx] = true;

            // Store the EXACT 5D coordinates 
            // This gives exact KDE at the collocation point
            box_coords[5*count + 0] = c.I1;
            box_coords[5*count + 1] = c.I2;
            box_coords[5*count + 2] = c.I3;
            box_coords[5*count + 3] = c.th1;
            box_coords[5*count + 4] = c.th3;
            count++;

            if(count % 10000 == 0)
                cout << "    " << count << "/" << N_box << " points" << endl;
        }
    }
}


// Phase 2: Parallel MC histogram
// Inner loop: step_EM -> to_5d -> cell_index -> histogram[idx]++
void MC_histogram(double* histogram) {
    cout << "MC histogram: " << N_thread << " threads, "
         << N_sample << " samples each" << endl;

    long long total_in = 0;
    long long total_out = 0;

#pragma omp parallel num_threads(N_thread) reduction(+:total_in,total_out)
    {
        int rank = omp_get_thread_num();
        mt19937 rng(random_device{}() + rank*137);
        normal_distribution<double> dist(0.0, 1.0);

        State6 X;
        X(0) = 1.0 + 0.1*rank; X(1) = 1.0;    X(2) = 0.1;
        X(3) = 0.1*rank;       X(4) = 0.0;     X(5) = -0.1*rank;

        // Per-thread burn-in
        for(int i = 0; i < Burn_in; i++)
            step_EM(X, rng, dist);

        long long local_in = 0, local_out = 0;

        for(long long step = 0; step < N_sample; step++) {
            step_EM(X, rng, dist);
            Coords5D c = to_5d(X);
            long long idx = cell_index(c);

            if(idx >= 0) {
                // Shared histogram: rare collisions acceptable for statistics
                #pragma omp atomic
                histogram[idx] += 1.0;
                local_in++;
            } else {
                local_out++;
            }

            if(rank == 0 && step % 5000000LL == 0 && step > 0)
                cout << "  thread 0: " << step/1000000 << "M/"
                     << N_sample/1000000 << "M" << endl;
        }

        total_in += local_in;
        total_out += local_out;
    }

    long long total_samples = (long long)N_thread * N_sample;
    cout << "Total samples: " << total_samples << endl;
    cout << "In-domain: " << total_in
         << " (" << 100.0*total_in/total_samples << "%)" << endl;
    cout << "Out-of-domain: " << total_out
         << " (" << 100.0*total_out/total_samples << "%)" << endl;
}

// kernel: k(u) = (3/4)(1 - u^2) for |u| <= 1, else 0
inline double epanechnikov(double u) {
    if(u > 1.0 || u < -1.0) return 0.0;
    return 0.75 * (1.0 - u * u);
}

// Phase 3: Post-convolution (kernel smoothing on histogram)
// For each collocation point, convolve histogram with product
// kernel over neighbor cells
void post_convolution(
    const double*        histogram,
    const vector<double>& box_coords,
    vector<double>&      density,
    long long            total_samples
) {
    cout << "Post-convolution: " << N_box << " collocation points, "
         << "neighbor radius = (" << R1 << "," << R2 << "," << R3
         << "," << Rt1 << "," << Rt3 << ")" << endl;

    int neighbors_per_point = (2*R1+1)*(2*R2+1)*(2*R3+1)*(2*Rt1+1)*(2*Rt3+1);
    cout << "  Max neighbor cells per point: " << neighbors_per_point << endl;

    long long total_kernel_evals = 0;

    #pragma omp parallel for reduction(+:total_kernel_evals) schedule(dynamic, 256)
    for(int i = 0; i < N_box; i++) {
        double xi_I1  = box_coords[5*i + 0];
        double xi_I2  = box_coords[5*i + 1];
        double xi_I3  = box_coords[5*i + 2];
        double xi_th1 = box_coords[5*i + 3];
        double xi_th3 = box_coords[5*i + 4];

        // Grid coordinates of this collocation point
        int cn1  = (int)floor((xi_I1 - I1_lo) / dI1);
        int cn2  = (int)floor((xi_I2 - I2_lo) / dI2);
        int cn3  = (int)floor((xi_I3 - I3_lo) / dI3);
        int cnt1 = (int)floor((xi_th1 - Th_lo) / dTh1);
        int cnt3 = (int)floor((xi_th3 - Th_lo) / dTh3);

        cn1  = max(0, min(cn1, G1 - 1));
        cn2  = max(0, min(cn2, G2 - 1));
        cn3  = max(0, min(cn3, G3 - 1));
        cnt1 = max(0, min(cnt1, Gt1 - 1));
        cnt3 = max(0, min(cnt3, Gt3 - 1));

        double sum = 0.0;

        // Loop over neighbor cells
        for(int d1 = cn1 - R1; d1 <= cn1 + R1; d1++) {
            if(d1 < 0 || d1 >= G1) continue;  // reflecting boundary for I1
            double c1 = I1_lo + d1 * dI1 + dI1 * 0.5;  // cell center
            double u1 = (xi_I1 - c1) / h1;
            double k1 = epanechnikov(u1);
            if(k1 == 0.0) continue;

        for(int d2 = cn2 - R2; d2 <= cn2 + R2; d2++) {
            if(d2 < 0 || d2 >= G2) continue;
            double c2 = I2_lo + d2 * dI2 + dI2 * 0.5;
            double u2 = (xi_I2 - c2) / h2;
            double k2 = epanechnikov(u2);
            if(k2 == 0.0) continue;

        for(int d3 = cn3 - R3; d3 <= cn3 + R3; d3++) {
            if(d3 < 0 || d3 >= G3) continue;
            double c3 = I3_lo + d3 * dI3 + dI3 * 0.5;
            double u3 = (xi_I3 - c3) / h3;
            double k3 = epanechnikov(u3);
            if(k3 == 0.0) continue;

        for(int dt1 = cnt1 - Rt1; dt1 <= cnt1 + Rt1; dt1++) {
            // Periodic wrap for theta1
            int dt1w = ((dt1 % Gt1) + Gt1) % Gt1;
            double ct1 = Th_lo + dt1w * dTh1 + dTh1 * 0.5;
            // Periodic distance
            double delta_t1 = xi_th1 - ct1;
            if(delta_t1 >  M_PI) delta_t1 -= 2.0 * M_PI;
            if(delta_t1 < -M_PI) delta_t1 += 2.0 * M_PI;
            double u_t1 = delta_t1 / ht1;
            double kt1 = epanechnikov(u_t1);
            if(kt1 == 0.0) continue;

        for(int dt3 = cnt3 - Rt3; dt3 <= cnt3 + Rt3; dt3++) {
            // Periodic wrap for theta3
            int dt3w = ((dt3 % Gt3) + Gt3) % Gt3;
            double ct3 = Th_lo + dt3w * dTh3 + dTh3 * 0.5;
            double delta_t3 = xi_th3 - ct3;
            if(delta_t3 >  M_PI) delta_t3 -= 2.0 * M_PI;
            if(delta_t3 < -M_PI) delta_t3 += 2.0 * M_PI;
            double u_t3 = delta_t3 / ht3;
            double kt3 = epanechnikov(u_t3);
            if(kt3 == 0.0) continue;

            // Product kernel weight (unnormalized — we normalize at the end)
            double w = k1 * k2 * k3 * kt1 * kt3;

            long long cell_idx = cell_index_int(d1, d2, d3, dt1w, dt3w);
            sum += histogram[cell_idx] * w;

            total_kernel_evals++;
        }}}}}

        // Normalize: divide by (N * prod(h_d))
        double h_prod = h1 * h2 * h3 * ht1 * ht3;
        density[i] = sum / ((double)total_samples * h_prod);
    }

    cout << "  Total kernel evaluations: " << total_kernel_evals << endl;
}


// Main
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    // Parse command line arguments
    if(argc > 1) N_thread = atoi(argv[1]);
    if(argc > 2) gamma_val = atof(argv[2]);
    if(argc > 3) T1 = atof(argv[3]);
    if(argc > 4) T3 = atof(argv[4]);
    if(argc > 5) N_box = atoi(argv[5]);
    if(argc > 6) N_sample = atoll(argv[6]);

    cout << "=== NLS 5D: Pre-binning + Post-convolution KDE ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "dt=" << dt << " Burn_in=" << Burn_in << endl;
    cout << "N_box=" << N_box << " N_sample/thread=" << N_sample << endl;
    cout << "N_thread=" << N_thread << endl;
    cout << endl;
    cout << "Coarse grid: " << G1 << "x" << G2 << "x" << G3
         << "x" << Gt1 << "x" << Gt3 << " = " << total_cells << " cells" << endl;
    cout << "Grid spacing: dI1=" << dI1 << " dI2=" << dI2 << " dI3=" << dI3
         << " dTh1=" << dTh1 << " dTh3=" << dTh3 << endl;
    cout << "KDE bandwidth (alpha=" << alpha_bw << "): h1=" << h1 << " h2=" << h2
         << " h3=" << h3 << " ht1=" << ht1 << " ht3=" << ht3 << endl;
    cout << "Neighbor radius: R=(" << R1 << "," << R2 << "," << R3
         << "," << Rt1 << "," << Rt3 << ")" << endl;
    cout << "Histogram memory: " << total_cells * 8 / (1024*1024) << " MB" << endl;
    cout << endl;

    // ----------------------------------------------------------------
    // Phase 0: Allocate histogram
    // ----------------------------------------------------------------
    cout << "Allocating histogram (" << total_cells << " cells)..." << endl;
    double* histogram = new double[total_cells];
    memset(histogram, 0, total_cells * sizeof(double));

    vector<bool> cell_occupied(total_cells, false);

    // ----------------------------------------------------------------
    // Phase 1: Create collocation points
    // ----------------------------------------------------------------
    cout << "Creating collocation points..." << endl;
    vector<double> box_coords(5 * N_box);
    create_collocation_points(box_coords, cell_occupied);
    cout << "Done. " << N_box << " collocation points created." << endl;

    // Free the occupancy map (not needed anymore)
    cell_occupied.clear();
    cell_occupied.shrink_to_fit();

    // ----------------------------------------------------------------
    // Phase 2: Parallel MC histogram
    // ----------------------------------------------------------------
    struct timeval mc_start, mc_end;
    gettimeofday(&mc_start, NULL);

    MC_histogram(histogram);

    gettimeofday(&mc_end, NULL);
    double mc_sec = ((mc_end.tv_sec - mc_start.tv_sec)*1e6
                    + mc_end.tv_usec - mc_start.tv_usec) / 1e6;
    cout << "MC histogram time: " << mc_sec << "s" << endl;

    // histogram stats
    long long nonzero_cells = 0;
    double max_count = 0;
    for(long long i = 0; i < total_cells; i++) {
        if(histogram[i] > 0) nonzero_cells++;
        if(histogram[i] > max_count) max_count = histogram[i];
    }
    cout << "Nonzero histogram cells: " << nonzero_cells
         << " / " << total_cells
         << " (" << 100.0*nonzero_cells/total_cells << "%)" << endl;
    cout << "Max histogram count: " << max_count << endl;
    cout << endl;

    // ----------------------------------------------------------------
    // Phase 3: Post-convolution (kernel smoothing)
    // ----------------------------------------------------------------
    struct timeval conv_start, conv_end;
    gettimeofday(&conv_start, NULL);

    long long total_samples = (long long)N_thread * N_sample;
    vector<double> density(N_box, 0.0);
    post_convolution(histogram, box_coords, density, total_samples);

    gettimeofday(&conv_end, NULL);
    double conv_sec = ((conv_end.tv_sec - conv_start.tv_sec)*1e6
                      + conv_end.tv_usec - conv_start.tv_usec) / 1e6;
    cout << "Post-convolution time: " << conv_sec << "s" << endl;

    // ----------------------------------------------------------------
    // Phase 4: Output
    // ----------------------------------------------------------------
    ofstream f1("NLS_FP_boxes.txt");
    ofstream f2("NLS_FP_density.txt");
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
    double sum_d = 0.0;
    for(int i = 0; i < N_box; i++) {
        if(density[i] > 0) nz++;
        sum_d += density[i];
    }
    cout << endl;
    cout << "=== Results ===" << endl;
    cout << "Max density = " << mx << endl;
    cout << "Min density = " << mn << endl;
    cout << "Nonzero density = " << nz << "/" << N_box << endl;
    cout << "Mean density = " << sum_d / N_box << endl;

    delete[] histogram;

    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec - t1.tv_sec)*1e6 + t2.tv_usec - t1.tv_usec) / 1e6;
    cout << "Total wall time = " << sec << "s" << endl;

    return 0;
}