/*
 * NLS5D_MC.cpp
 * 
 * Monte Carlo sampler for the 5D NLS energy cascade system
 * in action-angle variables (I1, I2, I3, theta1, theta3)
 * where theta_j = 2(phi_j - phi_2) for j = 1, 3.
 *
 * Based on equation (4.1) from Summary_NLS.pdf.
 * Adapted from 4DVDP_omp_nomesh.cpp using the splitting method
 * from Zhai, Dobson, Li (2021) for high-dimensional MC sampling.
 *
 * ====================================================================
 * THE 5D SYSTEM
 * ====================================================================
 *
 * State: X = (I1, I2, I3, theta1, theta3)
 *   where theta1 = 2(phi1 - phi2), theta3 = 2(phi3 - phi2)
 *
 * From eq (4.1) in Summary_NLS.pdf:
 *
 * dI1 = [4*I2*I1*(sin(theta1) - gamma*(1 + cos(theta1)))
 *         + 4*gamma*T1 - 2*gamma*(I1^2 + 2*I3*I1)] dt
 *        + 2*sqrt(2*gamma*T1*I1) dW1^I
 *
 * dI2 = -4*I2*(I1*sin(theta1) + I3*sin(theta3)) dt
 *        [NO NOISE on I2]
 *
 * dI3 = [4*I2*I3*(sin(theta3) - gamma*(1 + cos(theta3)))
 *         + 4*gamma*T3 - 2*gamma*(I3^2 + 2*I3*I1)] dt
 *        + 2*sqrt(2*gamma*T3*I3) dW3^I
 *
 * For the theta equations, we use theta_j = 2(phi_j - phi_2):
 *   d(theta1) = 2*(d(phi1) - d(phi2))
 *   d(theta3) = 2*(d(phi3) - d(phi2))
 *
 * From eq (4.1):
 *   dphi1 = [2*I2*(1 + cos(theta1) - gamma*sin(theta1)) + I1 + 2*I3] dt
 *           + sqrt(2*gamma*T1) * I1^{-1/2} dW1^phi
 *   dphi2 = [I2 + 2*I1*(1 + cos(theta1)) + 2*I3*(1 + cos(theta3))] dt
 *           [NO NOISE on phi2]
 *   dphi3 = [2*I2*(1 + cos(theta3) - gamma*sin(theta3)) + I3 + 2*I1] dt
 *           + sqrt(2*gamma*T3) * I3^{-1/2} dW3^phi
 *
 * Note the sign conventions: in eq (4.1), theta1 appears as 
 *   "Delta_{1,2}phi" = phi1 - phi2, and our theta1 = 2*(phi1-phi2),
 *   so cos(2*Delta_{1,2}phi) = cos(theta1), etc.
 *   Similarly "Delta_{2,1}phi" = phi2 - phi1 = -theta1/2
 *   so cos(2*Delta_{2,1}phi) = cos(theta1) (even function).
 *
 * Therefore:
 *   d(theta1) = 2*(dphi1 - dphi2)
 *     = 2*{[2*I2*(1+cos(theta1) - gamma*sin(theta1)) + I1 + 2*I3]
 *          - [I2 + 2*I1*(1+cos(theta1)) + 2*I3*(1+cos(theta3))]} dt
 *       + 2*sqrt(2*gamma*T1)*I1^{-1/2} dW1^phi
 *
 *   Drift of theta1:
 *     = 2*{ 2*I2 + 2*I2*cos(theta1) - 2*gamma*I2*sin(theta1) + I1 + 2*I3
 *           - I2 - 2*I1 - 2*I1*cos(theta1) - 2*I3 - 2*I3*cos(theta3) }
 *     = 2*{ I2*(1 + 2*cos(theta1) - 2*gamma*sin(theta1))
 *           - I1*(1 + 2*cos(theta1)) - 2*I3*cos(theta3) }
 *     [Simplified form from page 20 of Summary_NLS.pdf]
 *     = 2*{ I2*(1 + 2*cos(theta1) + 2*gamma*sin(theta1))
 *           - I1 - 2*I1*cos(theta1) - 2*I3*cos(theta3) }
 *
 *   Wait - let me re-derive carefully. The phi equations in (4.1) use
 *   Delta_{2,1}phi notation. Let me use the explicit form from page 20:
 *
 *   dtheta1 = [I2*(1 + 2*cos(theta1) + 2*gamma*sin(theta1))
 *              - I1 - 2*I1*cos(theta1) - 2*I3*cos(theta3)] dt
 *             + sqrt(2*gamma*T1)*I1^{-1/2} dW1^phi
 *
 *   dtheta3 = [I2*(1 + 2*cos(theta3) + 2*gamma*sin(theta3))
 *              - I3 - 2*I3*cos(theta3) - 2*I1*cos(theta1)] dt
 *             + sqrt(2*gamma*T3)*I3^{-1/2} dW3^phi
 *
 *   [These match the equations on page 20 of the summary]
 *
 * Noise structure: 4 independent Wiener processes
 *   W1^I, W3^I (for I1, I3)
 *   W1^phi, W3^phi (for theta1, theta3)
 *   I2 and the deterministic part have no noise.
 *
 * ====================================================================
 * SPLITTING METHOD FOR 5D
 * ====================================================================
 *
 * We split dimensions into groups for efficient box lookup:
 *   Group A: (I1, I2)   -> index in [0, N*N)         -> n_I1 * N + n_I2
 *   Group B: (I3, theta1)-> index in [N*N, 2*N*N)    -> N*N + n_I3 * N + n_th1
 *   Group C: (theta3)    -> index in [2*N*N, 2*N*N+N) -> 2*N*N + n_th3
 *
 * A sample hits box j iff j is in the intersection of
 *   list_of_box[groupA_idx] ∩ list_of_box[groupB_idx] ∩ list_of_box[groupC_idx]
 *
 * Memory: 2*N*N + N entries (vs N^5 for full grid)
 *
 * ====================================================================
 * IMPORTANT NOTES
 * ====================================================================
 *
 * 1. I_j >= 0 always (action variables are non-negative).
 *    We need to handle the boundary I_j = 0 carefully in the integrator
 *    (reflecting boundary or clamping).
 *
 * 2. theta_j are periodic: theta_j in [-pi, pi].
 *    The grid for theta wraps around.
 *
 * 3. The noise coefficients sqrt(I_j) become singular at I_j = 0.
 *    We clamp I_j >= I_min > 0 to avoid numerical issues.
 *
 * ====================================================================
 */

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define NDEBUG
#define EIGEN_NO_DEBUG
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <Eigen/Dense>
#include <random>
#include <complex>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Eigen;
using namespace std;

// ====================================================================
// DOMAIN PARAMETERS
// ====================================================================
// I variables: [0, I_max]
// theta variables: [-pi, pi]

const double I_max    = 6.0;      // upper bound for I1, I2, I3
const double I_min    = 1e-8;     // lower clamp for I_j (avoid sqrt singularity)

const double lowI     = 0.0;      // lower bound for I grid
const double SpI      = I_max;    // span for I grid

const double lowTh    = -M_PI;    // lower bound for theta grid
const double SpTh     = 2.0*M_PI; // span for theta grid

const int N_I         = 200;      // grid points per I dimension
const int N_Th        = 200;      // grid points per theta dimension

const double hI       = SpI / N_I;    // grid spacing for I
const double hTh      = SpTh / N_Th;  // grid spacing for theta

// ====================================================================
// SIMULATION PARAMETERS
// ====================================================================
const int    N_box     = 20000;          // number of reference boxes
const long long N_sample = 1000000000LL; // MC samples per thread (1e9)
const int    Burn_in   = 1000000;        // burn-in steps
const double dt        = 0.0005;         // time step (smaller for stability)
const double ratio     = 1.0;            // fraction from trajectory vs uniform

// ====================================================================
// PHYSICAL PARAMETERS
// ====================================================================
double gamma_param = 0.1;   // coupling to heat bath
double T1          = 2.0;   // temperature of heat bath 1 (left)
double T3          = 1.0;   // temperature of heat bath 3 (right)

// ====================================================================
// SPLITTING STRUCTURE
// ====================================================================
// Group A: (I1, I2)    -> N_I * N_I entries
// Group B: (I3, theta1) -> N_I * N_Th entries
// Group C: (theta3)     -> N_Th entries
const int sizeA = N_I * N_I;
const int sizeB = N_I * N_Th;
const int sizeC = N_Th;
const int total_list_size = sizeA + sizeB + sizeC;

// ====================================================================
// 5D STATE TYPE
// ====================================================================
// X = (I1, I2, I3, theta1, theta3)
typedef Matrix<double, 5, 1> Vector5d;

// ====================================================================
// WRAP theta to [-pi, pi]
// ====================================================================
inline double wrap_angle(double theta)
{
    while(theta >  M_PI) theta -= 2.0*M_PI;
    while(theta < -M_PI) theta += 2.0*M_PI;
    return theta;
}

// ====================================================================
// DRIFT FUNCTION
// ====================================================================
Vector5d drift(const Vector5d& X)
{
    double I1 = X(0), I2 = X(1), I3 = X(2);
    double th1 = X(3), th3 = X(4);
    double g = gamma_param;
    
    double sin_th1 = sin(th1), cos_th1 = cos(th1);
    double sin_th3 = sin(th3), cos_th3 = cos(th3);
    
    Vector5d f;
    
    // dI1
    f(0) = 4.0*I2*I1*(sin_th1 - g*(1.0 + cos_th1))
           + 4.0*g*T1
           - 2.0*g*(I1*I1 + 2.0*I3*I1);
    
    // dI2 (no noise, purely deterministic)
    f(1) = -4.0*I2*(I1*sin_th1 + I3*sin_th3);
    
    // dI3
    f(2) = 4.0*I2*I3*(sin_th3 - g*(1.0 + cos_th3))
           + 4.0*g*T3
           - 2.0*g*(I3*I3 + 2.0*I3*I1);
    
    // dtheta1 (from page 20 of summary)
    f(3) = I2*(1.0 + 2.0*cos_th1 + 2.0*g*sin_th1)
           - I1 - 2.0*I1*cos_th1 - 2.0*I3*cos_th3;
    
    // dtheta3 (from page 20 of summary)
    f(4) = I2*(1.0 + 2.0*cos_th3 + 2.0*g*sin_th3)
           - I3 - 2.0*I3*cos_th3 - 2.0*I1*cos_th1;
    
    return f;
}

// ====================================================================
// DIFFUSION COEFFICIENTS (diagonal noise)
// ====================================================================
// dX = f(X)dt + sigma(X) dW
// where dW = (dW1^I, dW3^I, dW1^phi, dW3^phi)  (4 independent noises)
// 
// sigma maps R^4 -> R^5:
//   I1:     2*sqrt(2*gamma*T1*I1)  * dW1^I
//   I2:     0
//   I3:     2*sqrt(2*gamma*T3*I3)  * dW3^I
//   theta1: sqrt(2*gamma*T1)*I1^{-1/2} * dW1^phi  [factor 2 from theta=2*(phi1-phi2)]
//           Actually: d(theta1) = 2*d(phi1) - 2*d(phi2)
//           noise on phi1 = sqrt(2*gamma*T1)*I1^{-1/2} dW1^phi
//           noise on phi2 = 0
//           So noise on theta1 = 2*sqrt(2*gamma*T1)*I1^{-1/2} dW1^phi
//   theta3: 2*sqrt(2*gamma*T3)*I3^{-1/2} * dW3^phi

struct DiffCoeffs {
    double sig_I1;      // coefficient for dW1^I
    double sig_I3;      // coefficient for dW3^I  
    double sig_th1;     // coefficient for dW1^phi
    double sig_th3;     // coefficient for dW3^phi
};

DiffCoeffs diffusion(const Vector5d& X)
{
    double I1 = max(X(0), I_min);
    double I3 = max(X(2), I_min);
    double g = gamma_param;
    
    DiffCoeffs s;
    s.sig_I1  = 2.0*sqrt(2.0*g*T1*I1);
    s.sig_I3  = 2.0*sqrt(2.0*g*T3*I3);
    s.sig_th1 = 2.0*sqrt(2.0*g*T1) / sqrt(I1);
    s.sig_th3 = 2.0*sqrt(2.0*g*T3) / sqrt(I3);
    
    return s;
}

// ====================================================================
// EULER-MARUYAMA STEP
// ====================================================================
void EM_step(Vector5d& X, const double dt_step, 
             double dW1I, double dW3I, double dW1phi, double dW3phi)
{
    Vector5d f = drift(X);
    DiffCoeffs s = diffusion(X);
    
    X(0) += f(0)*dt_step + s.sig_I1  * dW1I;
    X(1) += f(1)*dt_step;  // no noise on I2
    X(2) += f(2)*dt_step + s.sig_I3  * dW3I;
    X(3) += f(3)*dt_step + s.sig_th1 * dW1phi;
    X(4) += f(4)*dt_step + s.sig_th3 * dW3phi;
    
    // Clamp I_j >= I_min (reflecting boundary at 0)
    if(X(0) < I_min) X(0) = I_min;
    if(X(1) < I_min) X(1) = I_min;
    if(X(2) < I_min) X(2) = I_min;
    
    // Wrap theta to [-pi, pi]
    X(3) = wrap_angle(X(3));
    X(4) = wrap_angle(X(4));
}

// ====================================================================
// BOX LOOKUP: which_box
// ====================================================================
// Returns the box index that the point X falls into, or -1 if none.
int which_box(vector<vector<int>>& list_of_box, const Vector5d& X)
{
    double I1 = X(0), I2 = X(1), I3 = X(2);
    double th1 = X(3), th3 = X(4);
    
    // Check bounds
    if(I1 < lowI || I1 >= lowI + SpI) return -1;
    if(I2 < lowI || I2 >= lowI + SpI) return -1;
    if(I3 < lowI || I3 >= lowI + SpI) return -1;
    // theta always in [-pi, pi] after wrapping, so always in range
    
    // Compute grid indices
    int n_I1  = (int)floor((I1 - lowI) / hI);
    int n_I2  = (int)floor((I2 - lowI) / hI);
    int n_I3  = (int)floor((I3 - lowI) / hI);
    int n_th1 = (int)floor((th1 - lowTh) / hTh);
    int n_th3 = (int)floor((th3 - lowTh) / hTh);
    
    // Clamp indices
    if(n_I1 < 0 || n_I1 >= N_I) return -1;
    if(n_I2 < 0 || n_I2 >= N_I) return -1;
    if(n_I3 < 0 || n_I3 >= N_I) return -1;
    if(n_th1 < 0) n_th1 = 0; if(n_th1 >= N_Th) n_th1 = N_Th - 1;
    if(n_th3 < 0) n_th3 = 0; if(n_th3 >= N_Th) n_th3 = N_Th - 1;
    
    // Group A: (I1, I2)
    int idxA = n_I1 * N_I + n_I2;
    // Group B: (I3, theta1)
    int idxB = sizeA + n_I3 * N_Th + n_th1;
    // Group C: (theta3)
    int idxC = sizeA + sizeB + n_th3;
    
    // Check if any box is registered in all three groups
    if(list_of_box[idxA].empty() || list_of_box[idxB].empty() || 
       list_of_box[idxC].empty())
        return -1;
    
    // 3-way intersection
    // First intersect A and B
    vector<int> tmp_AB;
    set_intersection(list_of_box[idxA].begin(), list_of_box[idxA].end(),
                     list_of_box[idxB].begin(), list_of_box[idxB].end(),
                     back_inserter(tmp_AB));
    if(tmp_AB.empty()) return -1;
    
    // Then intersect (A∩B) with C
    vector<int> tmp_ABC;
    set_intersection(tmp_AB.begin(), tmp_AB.end(),
                     list_of_box[idxC].begin(), list_of_box[idxC].end(),
                     back_inserter(tmp_ABC));
    
    if(tmp_ABC.empty()) return -1;
    return tmp_ABC[0];
}

// ====================================================================
// CREATE BOXES (reference points)
// ====================================================================
// Sample N_box reference points from the SDE trajectory, register them
// in the splitting structure.
void create_boxes(vector<vector<int>>& list_of_box, MatrixXd& Boxes)
{
    Vector5d x_old, x_new;
    // Initial condition: moderate values
    x_old << 1.0, 1.0, 1.0, 0.5, -0.5;
    
    mt19937 mt(42);  // fixed seed for reproducibility
    uniform_real_distribution<double> u(0.0, 1.0);
    normal_distribution<double> normal(0.0, 1.0);
    
    // Burn in
    for(int i = 0; i < Burn_in; i++)
    {
        double dW1I = sqrt(dt)*normal(mt);
        double dW3I = sqrt(dt)*normal(mt);
        double dW1p = sqrt(dt)*normal(mt);
        double dW3p = sqrt(dt)*normal(mt);
        EM_step(x_old, dt, dW1I, dW3I, dW1p, dW3p);
    }
    
    int count = 0;
    while(count < N_box)
    {
        // Sample from trajectory
        if(u(mt) < ratio)
        {
            bool in_domain = false;
            while(!in_domain)
            {
                for(int i = 0; i < 5000; i++)
                {
                    double dW1I = sqrt(dt)*normal(mt);
                    double dW3I = sqrt(dt)*normal(mt);
                    double dW1p = sqrt(dt)*normal(mt);
                    double dW3p = sqrt(dt)*normal(mt);
                    EM_step(x_old, dt, dW1I, dW3I, dW1p, dW3p);
                }
                x_new = x_old;
                if(x_new(0) >= lowI && x_new(0) < lowI + SpI &&
                   x_new(1) >= lowI && x_new(1) < lowI + SpI &&
                   x_new(2) >= lowI && x_new(2) < lowI + SpI)
                {
                    in_domain = true;
                }
            }
        }
        // else: random point (not used when ratio=1.0)
        
        // Compute grid indices
        int n_I1  = (int)floor((x_new(0) - lowI)  / hI);
        int n_I2  = (int)floor((x_new(1) - lowI)  / hI);
        int n_I3  = (int)floor((x_new(2) - lowI)  / hI);
        int n_th1 = (int)floor((x_new(3) - lowTh) / hTh);
        int n_th3 = (int)floor((x_new(4) - lowTh) / hTh);
        
        // Clamp
        n_I1  = max(0, min(n_I1,  N_I  - 1));
        n_I2  = max(0, min(n_I2,  N_I  - 1));
        n_I3  = max(0, min(n_I3,  N_I  - 1));
        n_th1 = max(0, min(n_th1, N_Th - 1));
        n_th3 = max(0, min(n_th3, N_Th - 1));
        
        int idxA = n_I1 * N_I + n_I2;
        int idxB = sizeA + n_I3 * N_Th + n_th1;
        int idxC = sizeA + sizeB + n_th3;
        
        // Check this box isn't already occupied (3-way intersection must be empty)
        vector<int> tmp_AB, tmp_ABC;
        set_intersection(list_of_box[idxA].begin(), list_of_box[idxA].end(),
                         list_of_box[idxB].begin(), list_of_box[idxB].end(),
                         back_inserter(tmp_AB));
        if(!tmp_AB.empty())
        {
            set_intersection(tmp_AB.begin(), tmp_AB.end(),
                             list_of_box[idxC].begin(), list_of_box[idxC].end(),
                             back_inserter(tmp_ABC));
        }
        
        if(tmp_ABC.empty())
        {
            // Register box
            list_of_box[idxA].push_back(count);
            sort(list_of_box[idxA].begin(), list_of_box[idxA].end());
            list_of_box[idxB].push_back(count);
            sort(list_of_box[idxB].begin(), list_of_box[idxB].end());
            list_of_box[idxC].push_back(count);
            sort(list_of_box[idxC].begin(), list_of_box[idxC].end());
            
            // Store box center
            Boxes(0, count) = lowI  + n_I1  * hI  + hI/2.0;
            Boxes(1, count) = lowI  + n_I2  * hI  + hI/2.0;
            Boxes(2, count) = lowI  + n_I3  * hI  + hI/2.0;
            Boxes(3, count) = lowTh + n_th1 * hTh + hTh/2.0;
            Boxes(4, count) = lowTh + n_th3 * hTh + hTh/2.0;
            count++;
            
            if(count % 5000 == 0)
                cout << "  boxes created: " << count << "/" << N_box << endl;
        }
    }
}

// ====================================================================
// MONTE CARLO SAMPLING (OpenMP parallelized)
// ====================================================================
void MC(VectorXd& Box_count, vector<vector<int>>& list_of_box, 
        const int N_thread)
{
    cout << "Begin parallel MC session (" << N_thread << " threads, "
         << N_sample << " samples each)" << endl;
    
    VectorXi count_vec(N_thread);
    
#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        
        mt19937 mt(random_device{}() + rank);
        normal_distribution<double> norm(0.0, 1.0);
        
        // Random initial condition
        Vector5d x_old;
        x_old << fabs(2.0*norm(mt)), fabs(2.0*norm(mt)), fabs(2.0*norm(mt)),
                 M_PI*norm(mt)/2.0, M_PI*norm(mt)/2.0;
        x_old(0) = max(x_old(0), I_min);
        x_old(1) = max(x_old(1), I_min);
        x_old(2) = max(x_old(2), I_min);
        x_old(3) = wrap_angle(x_old(3));
        x_old(4) = wrap_angle(x_old(4));
        
        // Burn in
        for(int i = 0; i < Burn_in; i++)
        {
            double dW1I = sqrt(dt)*norm(mt);
            double dW3I = sqrt(dt)*norm(mt);
            double dW1p = sqrt(dt)*norm(mt);
            double dW3p = sqrt(dt)*norm(mt);
            EM_step(x_old, dt, dW1I, dW3I, dW1p, dW3p);
        }
        
        int count_f = 0;
        for(long long i = 0; i < N_sample; i++)
        {
            double dW1I = sqrt(dt)*norm(mt);
            double dW3I = sqrt(dt)*norm(mt);
            double dW1p = sqrt(dt)*norm(mt);
            double dW3p = sqrt(dt)*norm(mt);
            EM_step(x_old, dt, dW1I, dW3I, dW1p, dW3p);
            
            int index = which_box(list_of_box, x_old);
            if(index != -1)
            {
                Box_count(rank * N_box + index) += 1.0;
                count_f++;
            }
            
            // Progress report
            if(rank == 0 && i % 100000000LL == 0 && i > 0)
                cout << "  thread 0: " << i/1000000 << "M / " 
                     << N_sample/1000000 << "M steps" << endl;
        }
        count_vec(rank) = count_f;
    }
    
    long long sum = 0;
    for(int i = 0; i < N_thread; i++)
        sum += count_vec(i);
    cout << "Total effective count = " << sum << endl;
    
    // Merge thread results
    for(int i = 1; i < N_thread; i++)
        for(int j = 0; j < N_box; j++)
            Box_count(j) += Box_count(i * N_box + j);
}

// ====================================================================
// CREATE BOXES ON A 2D SLICE
// ====================================================================
// Fix theta1=0, theta3=0, I2=I2_fix; vary I1, I3 on a grid.
// This is for visualizing the marginal density on a 2D slice.
void create_boxes_slice(vector<vector<int>>& list_of_box, MatrixXd& Boxes,
                        int N_slice, double I2_fix, double th1_fix, double th3_fix)
{
    double hI_slice = SpI / N_slice;
    
    for(int i = 0; i < N_slice; i++)
    {
        for(int j = 0; j < N_slice; j++)
        {
            int idx = i * N_slice + j;
            double I1_val  = lowI + i * hI_slice + hI_slice / 2.0;
            double I3_val  = lowI + j * hI_slice + hI_slice / 2.0;
            
            // Grid indices for lookup
            int n_I1  = (int)floor((I1_val  - lowI)  / hI);
            int n_I2  = (int)floor((I2_fix  - lowI)  / hI);
            int n_I3  = (int)floor((I3_val  - lowI)  / hI);
            int n_th1 = (int)floor((th1_fix - lowTh) / hTh);
            int n_th3 = (int)floor((th3_fix - lowTh) / hTh);
            
            n_I1  = max(0, min(n_I1,  N_I  - 1));
            n_I2  = max(0, min(n_I2,  N_I  - 1));
            n_I3  = max(0, min(n_I3,  N_I  - 1));
            n_th1 = max(0, min(n_th1, N_Th - 1));
            n_th3 = max(0, min(n_th3, N_Th - 1));
            
            int idxA = n_I1 * N_I + n_I2;
            int idxB = sizeA + n_I3 * N_Th + n_th1;
            int idxC = sizeA + sizeB + n_th3;
            
            // Register
            list_of_box[idxA].push_back(idx);
            sort(list_of_box[idxA].begin(), list_of_box[idxA].end());
            list_of_box[idxB].push_back(idx);
            sort(list_of_box[idxB].begin(), list_of_box[idxB].end());
            list_of_box[idxC].push_back(idx);
            sort(list_of_box[idxC].begin(), list_of_box[idxC].end());
            
            Boxes(0, idx) = I1_val;
            Boxes(1, idx) = I2_fix;
            Boxes(2, idx) = I3_val;
            Boxes(3, idx) = th1_fix;
            Boxes(4, idx) = th3_fix;
        }
    }
}

// ====================================================================
// MAIN
// ====================================================================
int main(int argc, char* argv[])
{
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    
    int N_thread = 10;
    
    // Parse command line arguments (optional)
    if(argc > 1) N_thread = atoi(argv[1]);
    if(argc > 2) gamma_param = atof(argv[2]);
    if(argc > 3) T1 = atof(argv[3]);
    if(argc > 4) T3 = atof(argv[4]);
    
    cout << "====================================" << endl;
    cout << "NLS 5D Monte Carlo Sampler" << endl;
    cout << "====================================" << endl;
    cout << "Parameters:" << endl;
    cout << "  gamma = " << gamma_param << endl;
    cout << "  T1    = " << T1 << endl;
    cout << "  T3    = " << T3 << endl;
    cout << "  dt    = " << dt << endl;
    cout << "  N_box = " << N_box << endl;
    cout << "  N_sample per thread = " << N_sample << endl;
    cout << "  N_thread = " << N_thread << endl;
    cout << "  I range  = [" << lowI << ", " << lowI+SpI << "]" << endl;
    cout << "  Th range = [" << lowTh << ", " << lowTh+SpTh << "]" << endl;
    cout << "  Grid: N_I=" << N_I << ", N_Th=" << N_Th << endl;
    cout << "  Splitting structure size = " << total_list_size << endl;
    cout << "  hI = " << hI << ", hTh = " << hTh << endl;
    cout << "====================================" << endl;
    
    // ---- Allocate splitting structure ----
    vector<vector<int>> list_of_box(total_list_size);
    for(int i = 0; i < total_list_size; i++)
        list_of_box[i].reserve(10);
    
    // ---- Create reference boxes ----
    MatrixXd Boxes(5, N_box);
    cout << "Creating reference boxes..." << endl;
    create_boxes(list_of_box, Boxes);
    cout << "Boxes generated." << endl;
    
    // ---- Run MC ----
    VectorXd Box_count(N_thread * N_box);
    Box_count.fill(0);
    MC(Box_count, list_of_box, N_thread);
    
    // Compute density: count / (total_samples * box_volume)
    // Box volume in 5D: hI^3 * hTh^2 (three I dims, two theta dims)
    double box_vol = pow(hI, 3) * pow(hTh, 2);
    long long total_samples = (long long)N_thread * N_sample;
    VectorXd density = Box_count.head(N_box) / (total_samples * box_vol);
    
    // ---- Output results ----
    cout << "Writing output files..." << endl;
    
    ofstream file_boxes("NLS5D_boxes.txt");
    ofstream file_density("NLS5D_density.txt");
    
    for(int i = 0; i < N_box; i++)
    {
        file_boxes << Boxes(0,i) << " " << Boxes(1,i) << " " 
                   << Boxes(2,i) << " " << Boxes(3,i) << " " 
                   << Boxes(4,i) << endl;
        file_density << density(i) << endl;
    }
    file_boxes.close();
    file_density.close();
    
    // ---- Summary statistics ----
    double max_density = density.maxCoeff();
    double mean_density = density.mean();
    int nonzero = 0;
    for(int i = 0; i < N_box; i++)
        if(density(i) > 0) nonzero++;
    
    cout << "====================================" << endl;
    cout << "Results summary:" << endl;
    cout << "  Max density  = " << max_density << endl;
    cout << "  Mean density = " << mean_density << endl;
    cout << "  Nonzero boxes = " << nonzero << "/" << N_box << endl;
    cout << "====================================" << endl;
    
    gettimeofday(&t2, NULL);
    double elapsed = ((t2.tv_sec - t1.tv_sec) * 1000000u + 
                      t2.tv_usec - t1.tv_usec) / 1.0e6;
    cout << "Total wall time = " << elapsed << " seconds" << endl;
    
    return 0;
}