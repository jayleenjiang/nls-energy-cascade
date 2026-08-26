
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <algorithm>
#include <omp.h>

using namespace Eigen;
using namespace std;

// Physical parameters
double gamma_val = 0.1;      
double T1        = 10.0;    
double T3        = 2.0;       
const int n_modes = 3;      

// ====================================================================
// Integrator parameters
// ====================================================================
const double dt        = 0.0001;      // fixed timestep 
const double sqrt_dt   = sqrt(dt);    
const int    Burn_in   = 2000000;     

// ====================================================================
// Domain parameters for 5D box lookup
// ====================================================================
const double I_max  = 50.0;          // upper bound for action variables
const double I_lo   = 0.0;           // lower bound for action variables
const double Th_lo  = -M_PI;         // lower bound for angle variables
const double Th_hi  = M_PI;          // upper bound for angle variables

const int N_I  = 150;                // grid cells per I dimension
const int N_Th = 150;                // grid cells per theta dimension
const double hI  = (I_max - I_lo) / N_I;    // grid spacing for I
const double hTh = (Th_hi - Th_lo) / N_Th;  // grid spacing for theta


const int sizeA = N_I * N_I;
const int sizeB = N_I * N_Th;
const int sizeC = N_Th;
const int total_list = sizeA + sizeB + sizeC;

// ====================================================================
// MC parameters 
// ====================================================================
int    N_box    = 50000;             // number of reference boxes
long long N_sample = 20000000;       // MC steps per thread
int    N_thread = 8;                 // number of OpenMP threads

// ====================================================================
// 6D (I1, I2, I3, phi1, phi2, phi3)
// ====================================================================
typedef Matrix<double, 6, 1> State6;

// Wrap angle to [-pi, pi]
inline double wrap(double x) {
    double res = x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
    return res;
}

// ====================================================================
// Euler integrator 
// ====================================================================
void step_EM(State6& X, mt19937& rng, normal_distribution<double>& dist)
{
    double I1 = X(0), I2 = X(1), I3 = X(2);
    double p1 = X(3), p2 = X(4), p3 = X(5);
    
    double d_phi_prev[3], d_phi_next[3];
    double I_prev[3], I_next[3];
    
    // Neighbor arrays with boundary conditions I_0 = I_4 = 0, phi_0 = phi_4 = 0
    I_prev[0] = 0.0;   I_next[0] = I2;
    I_prev[1] = I1;    I_next[1] = I3;
    I_prev[2] = I2;    I_next[2] = 0.0;
    
    // Phase differences: d_phi_prev[j] = 2*(phi_j - phi_{j-1})
    d_phi_prev[0] = 2.0*(p1 - 0.0);     
    d_phi_next[0] = 2.0*(p1 - p2);       
    d_phi_prev[1] = 2.0*(p2 - p1);       
    d_phi_next[1] = 2.0*(p2 - p3);       
    d_phi_prev[2] = 2.0*(p3 - p2);       
    d_phi_next[2] = 2.0*(p3 - 0.0);      
    
    double Ivec[3] = {I1, I2, I3};
    double M = I1 + I2 + I3;  // total mass
    
    // Hamiltonian drift (conservative part)
    double dI[3], dphi[3];
    for(int j = 0; j < 3; j++) {
        dI[j] = 4.0 * Ivec[j] * (I_prev[j]*sin(d_phi_prev[j]) 
                                 + I_next[j]*sin(d_phi_next[j]));
        dphi[j] = 2.0*M - Ivec[j] 
                  + 2.0*I_prev[j]*cos(d_phi_prev[j]) 
                  + 2.0*I_next[j]*cos(d_phi_next[j]);
    }
    
    // Heat bath dissipation on mode 1 (coupled to T1)
    dI[0] += 2.0*gamma_val*(2.0*T1 - (2.0*M*I1 - I1*I1 + 2.0*I_next[0]*I1*cos(d_phi_next[0])));
    dphi[0] += gamma_val*(2.0*I_next[0]*sin(d_phi_next[0]));
    
    // Heat bath dissipation on mode 3 (coupled to T3)
    dI[2] += 2.0*gamma_val*(2.0*T3 - (2.0*M*I3 - I3*I3 + 2.0*I_prev[2]*I3*cos(d_phi_prev[2])));
    dphi[2] += gamma_val*(2.0*I_prev[2]*sin(d_phi_prev[2]));
    
    // Generate noise increments dW = sqrt(dt) * N(0,1)
    double dW[4];
    for(int j = 0; j < 4; j++) dW[j] = sqrt_dt * dist(rng);
    
    // Clamp at 1e-14 
    double diff_I1   = 2.0 * sqrt(2.0 * gamma_val * T1 * max(I1, 1e-14));
    double diff_I3   = 2.0 * sqrt(2.0 * gamma_val * max(T3, 1e-10) * max(I3, 1e-14));
    double diff_phi1 = sqrt(2.0 * gamma_val * T1 / max(I1, 1e-14));
    double diff_phi3 = sqrt(2.0 * gamma_val * max(T3, 1e-10) / max(I3, 1e-14));
    
    // Euler 
    X(0) = I1  + dI[0]*dt + diff_I1*dW[0];
    X(1) = I2  + dI[1]*dt;   
    X(2) = I3  + dI[2]*dt + diff_I3*dW[1];
    X(3) = p1  + dphi[0]*dt + diff_phi1*dW[2];
    X(4) = p2  + dphi[1]*dt;
    X(5) = p3  + dphi[2]*dt + diff_phi3*dW[3];
    
    // Reflecting boundary: clamp I >= 0 to prevent negative actions
    X(0) = max(X(0), 1e-14);
    X(1) = max(X(1), 1e-14);
    X(2) = max(X(2), 1e-14);
    
    // Periodic boundary: wrap angles to [-pi, pi]
    X(3) = wrap(X(3));
    X(4) = wrap(X(4));
    X(5) = wrap(X(5));
}


struct Coords5D {
    double I1, I2, I3, th1, th3;
};

Coords5D to_5d(const State6& X) {
    Coords5D c;
    c.I1  = X(0);
    c.I2  = X(1);
    c.I3  = X(2);
    c.th1 = wrap(2.0*(X(3) - X(4)));  // theta1 = 2*(phi1 - phi2)
    c.th3 = wrap(2.0*(X(5) - X(4)));  // theta3 = 2*(phi3 - phi2)
    return c;
}

// ====================================================================
// Box lookup 
// ====================================================================
int which_box(vector<vector<int>>& lob, const Coords5D& c) {
    if(c.I1 < I_lo || c.I1 >= I_max) return -1;
    if(c.I2 < I_lo || c.I2 >= I_max) return -1;
    if(c.I3 < I_lo || c.I3 >= I_max) return -1;
    
    int n1 = (int)floor((c.I1 - I_lo) / hI);
    int n2 = (int)floor((c.I2 - I_lo) / hI);
    int n3 = (int)floor((c.I3 - I_lo) / hI);
    int nt1 = (int)floor((c.th1 - Th_lo) / hTh);
    int nt3 = (int)floor((c.th3 - Th_lo) / hTh);

    if(n1 < 0 || n1 >= N_I || n2 < 0 || n2 >= N_I || n3 < 0 || n3 >= N_I) return -1;
    nt1 = max(0, min(nt1, N_Th-1));
    nt3 = max(0, min(nt3, N_Th-1));
    
    int idxA = n1*N_I + n2;               // Group A: (I1, I2)
    int idxB = sizeA + n3*N_Th + nt1;     // Group B: (I3, theta1)
    int idxC = sizeA + sizeB + nt3;       // Group C: (theta3)
    
    if(lob[idxA].empty() || lob[idxB].empty() || lob[idxC].empty()) return -1;
    
    static thread_local vector<int> ab;
    static thread_local vector<int> abc;
    
    ab.clear();
    abc.clear();

    set_intersection(lob[idxA].begin(), lob[idxA].end(),
                     lob[idxB].begin(), lob[idxB].end(), back_inserter(ab));
    if(ab.empty()) return -1;
    
    set_intersection(ab.begin(), ab.end(),
                     lob[idxC].begin(), lob[idxC].end(), back_inserter(abc));
    
    return abc.empty() ? -1 : abc[0];
}

// ====================================================================
// Register a box 
// ====================================================================
bool register_box(vector<vector<int>>& lob, const Coords5D& c, int box_id) {
    int n1 = max(0, min((int)floor((c.I1 - I_lo)/hI),  N_I-1));
    int n2 = max(0, min((int)floor((c.I2 - I_lo)/hI),  N_I-1));
    int n3 = max(0, min((int)floor((c.I3 - I_lo)/hI),  N_I-1));
    int nt1 = max(0, min((int)floor((c.th1 - Th_lo)/hTh), N_Th-1));
    int nt3 = max(0, min((int)floor((c.th3 - Th_lo)/hTh), N_Th-1));
    
    int idxA = n1*N_I + n2;
    int idxB = sizeA + n3*N_Th + nt1;
    int idxC = sizeA + sizeB + nt3;
    
    // Check if cell is already occupied
    vector<int> ab;
    set_intersection(lob[idxA].begin(), lob[idxA].end(),
                     lob[idxB].begin(), lob[idxB].end(), back_inserter(ab));
    if(!ab.empty()) {
        vector<int> abc;
        set_intersection(ab.begin(), ab.end(),
                         lob[idxC].begin(), lob[idxC].end(), back_inserter(abc));
        if(!abc.empty()) return false;  // already occupied
    }
    
    // Register in all three groups
    lob[idxA].push_back(box_id); sort(lob[idxA].begin(), lob[idxA].end());
    lob[idxB].push_back(box_id); sort(lob[idxB].begin(), lob[idxB].end());
    lob[idxC].push_back(box_id); sort(lob[idxC].begin(), lob[idxC].end());
    return true;
}

// ====================================================================
// Create reference boxes: half from trajectory, half from uniform random
// (Following Zhai-Dobson-Li paper: trajectory samples cover high-density
//  regions, uniform samples cover low-density tails)
// ====================================================================
void create_boxes(vector<vector<int>>& lob, MatrixXd& Boxes) {
    const double sample_ratio = 0.5;  // fraction from trajectory

    State6 X;
    X << 1.0, 1.0, 0.1, 0.0, 0.0, 0.0;
    
    mt19937 rng(42);
    normal_distribution<double> dist(0.0, 1.0);
    uniform_real_distribution<double> uni(0.0, 1.0);
    
    // Burn in to reach steady state
    cout << "  Burning in..." << endl;
    for(int i = 0; i < Burn_in; i++)
        step_EM(X, rng, dist);
    
    cout << "  Sampling boxes (ratio=" << sample_ratio << ")..." << endl;
    int count = 0;
    int from_traj = 0, from_unif = 0;

    while(count < N_box) {
        Coords5D c;

        if(uni(rng) < sample_ratio) {
            // --- Sample from trajectory (importance sampling) ---
            for(int i = 0; i < 5000; i++)
                step_EM(X, rng, dist);
            c = to_5d(X);
            if(c.I1 < I_lo || c.I1 >= I_max) continue;
            if(c.I2 < I_lo || c.I2 >= I_max) continue;
            if(c.I3 < I_lo || c.I3 >= I_max) continue;
        } else {
            // --- Uniform random sample in domain ---
            c.I1  = I_lo + uni(rng) * (I_max - I_lo);
            c.I2  = I_lo + uni(rng) * (I_max - I_lo);
            c.I3  = I_lo + uni(rng) * (I_max - I_lo);
            c.th1 = Th_lo + uni(rng) * (Th_hi - Th_lo);
            c.th3 = Th_lo + uni(rng) * (Th_hi - Th_lo);
        }
        
        if(register_box(lob, c, count)) {
            // Store the center of the grid cell
            int n1 = max(0, min((int)floor((c.I1-I_lo)/hI), N_I-1));
            int n2 = max(0, min((int)floor((c.I2-I_lo)/hI), N_I-1));
            int n3 = max(0, min((int)floor((c.I3-I_lo)/hI), N_I-1));
            int nt1 = max(0, min((int)floor((c.th1-Th_lo)/hTh), N_Th-1));
            int nt3 = max(0, min((int)floor((c.th3-Th_lo)/hTh), N_Th-1));
            
            Boxes(0, count) = I_lo  + n1*hI  + hI/2.0;
            Boxes(1, count) = I_lo  + n2*hI  + hI/2.0;
            Boxes(2, count) = I_lo  + n3*hI  + hI/2.0;
            Boxes(3, count) = Th_lo + nt1*hTh + hTh/2.0;
            Boxes(4, count) = Th_lo + nt3*hTh + hTh/2.0;

            if(uni(rng) < sample_ratio) from_traj++; else from_unif++;
            count++;
            
            if(count % 5000 == 0)
                cout << "    " << count << "/" << N_box 
                     << " (traj=" << from_traj << " unif=" << from_unif << ")" << endl;
        }
    }
    cout << "  Final: " << from_traj << " from trajectory, " 
         << from_unif << " from uniform" << endl;
}

// ====================================================================
// Monte Carlo sampling
// ====================================================================
void MC(VectorXd& Box_count, vector<vector<int>>& lob) {
    cout << "MC: " << N_thread << " threads, " << N_sample << " samples each" << endl;
    
    vector<long long> counts(N_thread, 0);
    
#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        mt19937 rng(random_device{}() + rank*137);
        normal_distribution<double> dist(0.0, 1.0);
        
        // Slightly different initial conditions per thread
        State6 X;
        X << 1.0+0.1*rank, 1.0, 0.1, 0.1*rank, 0.0, -0.1*rank;
        
        // Per-thread burn-in
        for(int i = 0; i < Burn_in; i++)
            step_EM(X, rng, dist);
        
        long long ct = 0;
        long long steps_done = 0;
        
        while(steps_done < N_sample) {
            step_EM(X, rng, dist);
            Coords5D c = to_5d(X);
            int idx = which_box(lob, c);
            if(idx != -1) {
                // Since dt is fixed, simple counting is correct
                Box_count(rank*N_box + idx) += 1.0;
                ct++;
            }
            steps_done++;
            
            if(rank == 0 && steps_done % 2000000LL == 0)
                cout << "  thread 0: " << steps_done/1000000 << "M/" 
                     << N_sample/1000000 << "M" << endl;
        }
        counts[rank] = ct;
    }
    
    long long total = 0;
    for(int i = 0; i < N_thread; i++) total += counts[i];
    cout << "Total hits = " << total << endl;
    
    // Merge thread results
    for(int i = 1; i < N_thread; i++)
        for(int j = 0; j < N_box; j++)
            Box_count(j) += Box_count(i*N_box + j);
}

// ====================================================================
// Main
// ====================================================================
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
    
    cout << "=== FIXED STEP EULER-MARUYAMA ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "Fixed dt=" << dt << endl;
    cout << "N_box=" << N_box << " N_sample/thread=" << N_sample << endl;
    cout << "N_thread=" << N_thread << endl;
    
    // Allocate splitting structure
    vector<vector<int>> lob(total_list);
    for(auto& v : lob) v.reserve(8);
    MatrixXd Boxes(5, N_box);
    
    // Phase 1: create reference boxes
    cout << "Creating boxes..." << endl;
    create_boxes(lob, Boxes);
    cout << "Done." << endl;
    
    // Phase 2: run parallel MC
    VectorXd Box_count(N_thread * N_box);
    Box_count.fill(0);
    MC(Box_count, lob);
    
    // Phase 3: compute density = count / (total_steps * box_volume)
    // With fixed dt and +1.0 counting, dividing by total steps is correct
    double vol = pow(hI, 3) * pow(hTh, 2);
    long long total_samples = (long long)N_thread * N_sample;
    VectorXd density = Box_count.head(N_box) / (total_samples * vol);
    
    // Write output
    ofstream f1("NLS_FP_boxes.txt");
    ofstream f2("NLS_FP_density.txt");
    for(int i = 0; i < N_box; i++) {
        f1 << Boxes(0,i) << " " << Boxes(1,i) << " " << Boxes(2,i) 
           << " " << Boxes(3,i) << " " << Boxes(4,i) << endl;
        f2 << density(i) << endl;
    }
    f1.close(); f2.close();
    
    // Summary
    double mx = density.maxCoeff();
    int nz = 0;
    for(int i = 0; i < N_box; i++) if(density(i) > 0) nz++;
    cout << "Max density = " << mx << endl;
    cout << "Nonzero = " << nz << "/" << N_box << endl;
    
    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec-t1.tv_sec)*1e6 + t2.tv_usec-t1.tv_usec) / 1e6;
    cout << "Wall time = " << sec << "s" << endl;
    
    return 0;
}