/*
 * NLS_backward.cpp
 * 
 * For each initial point x_i (from LHS), runs N_sample MC trajectories,
 * records E_x[f(X_t)] at each time step, outputs the time series.
 * 
 * Variables: (I1, I2, I3, phi1, phi2, phi3) — 6D state
 * Observable: f(x) = max(0, c-I1) + max(0, c-I2) + max(0, c-I3)
 */

#define EIGEN_NO_DEBUG
#define NDEBUG
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <Eigen/Dense>
#include <random>
#include <omp.h>

using namespace Eigen;
using namespace std;

// ====================================================================
// Physical parameters 
// ====================================================================
float gamma_val = 0.1f;
float T1        = 5.0f;
float T3        = 5.0f;

// ====================================================================
// Integrator
// ====================================================================
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const float PI_f    = 3.14159265358979323846f;

// ====================================================================
// Backward solver parameters
// ====================================================================
int N_initials      = 4096;    // number of initial points (must be multiple of 16)
int T_steps         = 50;      // number of recorded time steps
int gap             = 100;     // record every gap*dt time units (gap*dt = 1.0)
long long N_sample  = 10000;    // MC trajectories per initial point
float obs_cutoff    = 3.0f;     // cutoff for localized observable
int N_thread        = 8;

// ====================================================================
// SIMD types and fast math 
// ====================================================================
typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

inline A16f wrap_pi_16(const A16f& x) {
    const float invTwoPi=0.159154943f, twoPi=6.283185307f;
    return x-(x*invTwoPi).round()*twoPi;
}
inline A16f fast_sin_16(const A16f& x) {
    const float B=1.27323954f, C=-0.40528473f, P=0.225f;
    auto y=B*x+C*x*x.abs(); return P*(y*y.abs()-y)+y;
}
inline A16f fast_cos_16(const A16f& x) {
    return fast_sin_16(wrap_pi_16(x+1.570796327f));
}

// ====================================================================
// RNG 
// ====================================================================
struct RNGState { uint32_t s[4]; };
RNGState init_rng(int rank) {
    RNGState st; uint64_t z=(uint64_t)rank+0x9E3779B97F4A7C15ULL;
    auto sm=[&z]()->uint64_t{z+=0x9E3779B97F4A7C15ULL;uint64_t r=z;
    r=(r^(r>>30))*0xBF58476D1CE4E5B9ULL;r=(r^(r>>27))*0x94D049BB133111EBULL;return r^(r>>31);};
    uint64_t p1=sm(),p2=sm();
    st.s[0]=(uint32_t)p1;st.s[1]=(uint32_t)(p1>>32);st.s[2]=(uint32_t)p2;st.s[3]=(uint32_t)(p2>>32);
    return st;
}
inline A16f next_u01_16(RNGState& s) {
    A16f res; for(int i=0;i<16;i++){
    uint32_t r=((s.s[0]+s.s[3])<<7)|((s.s[0]+s.s[3])>>25);uint32_t t=s.s[1]<<9;
    s.s[2]^=s.s[0];s.s[3]^=s.s[1];s.s[1]^=s.s[2];s.s[0]^=s.s[3];s.s[2]^=t;
    s.s[3]=(s.s[3]<<11)|(s.s[3]>>21);res(i)=(r>>9)*0.00000011920929f;} return res;
}
inline void gen_noise_4x16(RNGState& s, A16f& n1, A16f& n2, A16f& n3, A16f& n4){
    const float twoPi=6.283185307f, piO2=1.570796327f;
    A16f u1=next_u01_16(s),u2=next_u01_16(s);
    A16f r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt(),th=twoPi*u2;
    n1=r*fast_sin_16(wrap_pi_16(th+piO2));n2=r*fast_sin_16(wrap_pi_16(th));
    u1=next_u01_16(s);u2=next_u01_16(s);
    r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt();th=twoPi*u2;
    n3=r*fast_sin_16(wrap_pi_16(th+piO2));n4=r*fast_sin_16(wrap_pi_16(th));
}

// ====================================================================
// NLS Euler-Maruyama step 
// ====================================================================
inline void step_EM_batch16(State6x16& X, RNGState& rng) {
    A16f I1=X.row(0).array(),I2=X.row(1).array(),I3=X.row(2).array();
    A16f p1=X.row(3).array(),p2=X.row(4).array(),p3=X.row(5).array();
    A16f d12w=wrap_pi_16(2.0f*(p1-p2)),d32w=wrap_pi_16(2.0f*(p3-p2));
    A16f s12=fast_sin_16(d12w),c12=fast_cos_16(d12w);
    A16f s32=fast_sin_16(d32w),c32=fast_cos_16(d32w);
    A16f M=I1+I2+I3;
    // Hamiltonian drift
    A16f dI1=4.0f*I1*I2*s12, dI2=4.0f*I2*(-I1*s12-I3*s32), dI3=4.0f*I3*I2*s32;
    A16f dp1=2.0f*M-I1+2.0f*I2*c12, dp2=2.0f*M-I2+2.0f*I1*c12+2.0f*I3*c32, dp3=2.0f*M-I3+2.0f*I2*c32;
    // Dissipation
    dI1+=2.0f*gamma_val*(2.0f*T1-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    dI3+=2.0f*gamma_val*(2.0f*T3-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);
    // Noise
    A16f nI1,nI3,np1,np3; gen_noise_4x16(rng,nI1,nI3,np1,np3);
    A16f I1c=I1.max(1e-14f),I3c=I3.max(1e-14f);
    A16f sI1=2.0f*(2.0f*gamma_val*T1*I1c).sqrt(),sI3=2.0f*(2.0f*gamma_val*T3*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T1/I1c).sqrt(),sp3=(2.0f*gamma_val*T3/I3c).sqrt();
    // Update
    X.row(0).array()=(I1+dI1*dt+sI1*sqrt_dt*nI1).max(1e-14f);
    X.row(1).array()=(I2+dI2*dt).max(1e-14f);
    X.row(2).array()=(I3+dI3*dt+sI3*sqrt_dt*nI3).max(1e-14f);
    X.row(3).array()=p1+dp1*dt+sp1*sqrt_dt*np1;
    X.row(4).array()=p2+dp2*dt;
    X.row(5).array()=p3+dp3*dt+sp3*sqrt_dt*np3;
    X.row(3)=wrap_pi_16(X.row(3).array());
    X.row(4)=wrap_pi_16(X.row(4).array());
    X.row(5)=wrap_pi_16(X.row(5).array());
}

// ====================================================================
// Observable: localized bump in action space
// f(x) = max(0, c-I1) + max(0, c-I2) + max(0, c-I3)
// ====================================================================
inline float observable_scalar(float I1, float I2, float I3) {
    return fmax(0.0f, obs_cutoff - I1) + fmax(0.0f, obs_cutoff - I2) + fmax(0.0f, obs_cutoff - I3);
}

// ====================================================================
// Load initial points from file
// ====================================================================
MatrixXf loadCSV(const string& path, int rows, int cols) {
    MatrixXf matrix(rows, cols);
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Error: could not open " << path << endl;
        return matrix;
    }
    string line;
    int row = 0;
    while (getline(file, line) && row < rows) {
        stringstream ss(line);
        string cell;
        int col = 0;
        while (getline(ss, cell, ' ') && col < cols) {
            try { matrix(row, col) = stof(cell); }
            catch (...) { matrix(row, col) = 0.0f; }
            col++;
        }
        row++;
    }
    return matrix;
}

// ====================================================================
// Main MC loop
// ====================================================================
void MC_backward(MatrixXf& Initial, MatrixXf& Obs_traj, const int nthread) {
#pragma omp parallel num_threads(nthread)
    {
        int rank = omp_get_thread_num();
        RNGState rng = init_rng(rank);

#pragma omp for
        for (long long int n = 0; n < N_initials / 16; n++) {
            State6x16 state_batch;

            for (int m = 0; m < N_sample; m++) {
                // Load initial conditions for this batch of 16
                for (int i = 0; i < 16; i++) {
                    state_batch.col(i) = Initial.col(n * 16 + i);
                }

                // Run trajectory, record observable at each gap
                for (int j = 0; j <= T_steps * gap; j++) {
                    // Record observable
                    if (j % gap == 0) {
                        int t_idx = j / gap;
                        for (int i = 0; i < 16; i++) {
                            float I1 = state_batch(0, i);
                            float I2 = state_batch(1, i);
                            float I3 = state_batch(2, i);
                            Obs_traj(n * 16 + i, t_idx) += observable_scalar(I1, I2, I3);
                        }
                    }
                    // EM step
                    if (j < T_steps * gap) {
                        step_EM_batch16(state_batch, rng);
                    }
                }
            }
        }
    }
}

// ====================================================================
// Main
// ====================================================================
int main(int argc, char* argv[]) {
    // Parse command line
    if (argc >= 2) N_thread   = atoi(argv[1]);
    if (argc >= 3) gamma_val  = atof(argv[2]);
    if (argc >= 4) T1         = atof(argv[3]);
    if (argc >= 5) T3         = atof(argv[4]);
    if (argc >= 6) N_sample   = atoll(argv[5]);

    cout << "=== NLS Backward Solver ===" << endl;
    cout << "gamma=" << gamma_val << "  T1=" << T1 << "  T3=" << T3 << endl;
    cout << "N_initials=" << N_initials << "  N_sample=" << N_sample << endl;
    cout << "T_steps=" << T_steps << "  gap=" << gap << "  dt=" << dt << endl;
    cout << "Observable cutoff c=" << obs_cutoff << endl;
    cout << "Threads=" << N_thread << endl;

    // Load initial points (5D: I1, I2, I3, theta1, theta3)
    // Convert to 6D: (I1, I2, I3, phi1, phi2=0, phi3)
    // where theta1 = 2*(phi1 - phi2), theta3 = 2*(phi3 - phi2)
    // so phi1 = theta1/2, phi3 = theta3/2 (with phi2 = 0)
    MatrixXf X_5d = loadCSV("NLS_backward_LHS_X_train.txt", N_initials, 5);
    cout << "Loaded " << X_5d.rows() << " initial points" << endl;

    // Convert 5D → 6D
    MatrixXf X_6d(N_initials, 6);
    for (int i = 0; i < N_initials; i++) {
        X_6d(i, 0) = X_5d(i, 0);          // I1
        X_6d(i, 1) = X_5d(i, 1);          // I2
        X_6d(i, 2) = X_5d(i, 2);          // I3
        X_6d(i, 3) = X_5d(i, 3) / 2.0f;  // phi1 = theta1/2
        X_6d(i, 4) = 0.0f;                // phi2 = 0
        X_6d(i, 5) = X_5d(i, 4) / 2.0f;  // phi3 = theta3/2
    }
    MatrixXf X_6dT = X_6d.transpose();  // 6 x N_initials

    // Allocate output: N_initials x (T_steps+1)
    MatrixXf Obs_traj(N_initials, T_steps + 1);
    Obs_traj.fill(0.0f);

    // Run MC
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    MC_backward(X_6dT, Obs_traj, N_thread);
    gettimeofday(&t2, NULL);

    double wall = ((t2.tv_sec - t1.tv_sec) * 1e6 + t2.tv_usec - t1.tv_usec) / 1e6;
    cout << "MC time: " << wall << "s" << endl;

    // Normalize by N_sample
    Obs_traj /= (float)N_sample;

    // Save output
    ofstream myfile("NLS_backward_Y_train.txt");
    for (int i = 0; i < N_initials; i++) {
        for (int j = 0; j <= T_steps; j++) {
            myfile << Obs_traj(i, j) << " ";
        }
        myfile << endl;
    }
    myfile.close();

    // Also save the 5D initial points used
    ofstream xfile("NLS_backward_X_train.txt");
    for (int i = 0; i < N_initials; i++) {
        for (int j = 0; j < 5; j++) {
            xfile << X_5d(i, j) << " ";
        }
        xfile << endl;
    }
    xfile.close();

    cout << "Output written: NLS_backward_Y_train.txt (" 
         << N_initials << " x " << T_steps + 1 << ")" << endl;
    cout << "Wall time: " << wall << "s" << endl;

    return 0;
}
