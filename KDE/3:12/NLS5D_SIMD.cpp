#define EIGEN_NO_DEBUG
#define NDEBUG
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
float gamma_val = 0.1f;
float T1        = 10.0f;
float T3        = 2.0f;

// Integrator parameters
const float dt      = 0.0001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

// Domain parameters
const float I_max  = 50.0f;
const float I_lo   = 0.0f;
const float PI_f   = 3.14159265358979323846f;
const float Th_lo  = -PI_f;
const float Th_hi  = PI_f;

const int N_I  = 150;
const int N_Th = 150;
const float hI  = (I_max - I_lo) / N_I;
const float hTh = (Th_hi - Th_lo) / N_Th;

const int sizeA = N_I * N_I;
const int sizeB = N_I * N_Th;
const int sizeC = N_Th;
const int total_list = sizeA + sizeB + sizeC;

// ================
// MC parameters
// ================
int       N_box    = 50000;
long long N_sample = 20000000LL;
int       N_thread = 8;
const double sample_ratio = 0.5;

typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 5, 16, RowMajor> Block5x16;  // 5D state: (I1,I2,I3,th1,th3)

// ===========================
// wrap angle to [-pi, pi]
// ===========================
inline A16f wrap_pi_16(const A16f& x) {
    const float invTwoPi = 0.159154943f;
    const float twoPi    = 6.283185307f;
    return x - (x * invTwoPi).round() * twoPi;
}

// ========================
// Fast sin 
// =======================
inline A16f fast_sin_16(const A16f& x) {
    const float B = 1.27323954f;   // 4/pi
    const float C = -0.40528473f;  // -4/pi^2
    const float P = 0.225f;
    auto y = B * x + C * x * x.abs();
    return P * (y * y.abs() - y) + y;
}

// =========================
// Fast cos 
// ========================
inline A16f fast_cos_16(const A16f& x) {
    const float piOver2 = 1.570796327f;
    return fast_sin_16(wrap_pi_16(x + piOver2));
}

// =================
// Xoshiro128++ RNG state
// =================
struct RNGState {
    uint32_t s[4];
};

RNGState initialize_rng_for_rank(int rank) {
    RNGState state;
    uint64_t z = static_cast<uint64_t>(rank) + 0x9E3779B97F4A7C15ULL;
    auto splitmix64 = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t r = z;
        r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9ULL;
        r = (r ^ (r >> 27)) * 0x94D049BB133111EBULL;
        return r ^ (r >> 31);
    };
    uint64_t p1 = splitmix64(), p2 = splitmix64();
    state.s[0] = (uint32_t)p1;
    state.s[1] = (uint32_t)(p1 >> 32);
    state.s[2] = (uint32_t)p2;
    state.s[3] = (uint32_t)(p2 >> 32);
    return state;
}

inline A16f next_u01_16(RNGState& state) {
    A16f res;
    for(int i = 0; i < 16; i++) {
        const uint32_t result = ((state.s[0] + state.s[3]) << 7) |
                                ((state.s[0] + state.s[3]) >> 25);
        const uint32_t t = state.s[1] << 9;
        state.s[2] ^= state.s[0]; state.s[3] ^= state.s[1];
        state.s[1] ^= state.s[2]; state.s[0] ^= state.s[3];
        state.s[2] ^= t;
        state.s[3] = (state.s[3] << 11) | (state.s[3] >> 21);
        res(i) = (result >> 9) * 0.00000011920929f;
    }
    return res;
}

// ====================================================================
// Generate 4 rows of N(0,1) noise for 16 trajectories 
// Rows: dW_I1, dW_I3, dW_phi1, dW_phi3
// ====================================================================
inline void generate_noise_4x16(RNGState& state,
    A16f& nI1, A16f& nI3, A16f& np1, A16f& np3)
{
    const float twoPi = 6.283185307f;
    const float piOver2 = 1.570796327f;

    // Pair 1: dW_I1, dW_I3
    A16f u1 = next_u01_16(state);
    A16f u2 = next_u01_16(state);
    A16f radius = (-2.0f * (1.0f - u1).max(1e-30f).log()).sqrt();
    A16f theta  = twoPi * u2;
    nI1 = radius * fast_sin_16(wrap_pi_16(theta + piOver2));
    nI3 = radius * fast_sin_16(wrap_pi_16(theta));

    // Pair 2: dW_phi1, dW_phi3
    u1 = next_u01_16(state);
    u2 = next_u01_16(state);
    radius = (-2.0f * (1.0f - u1).max(1e-30f).log()).sqrt();
    theta  = twoPi * u2;
    np1 = radius * fast_sin_16(wrap_pi_16(theta + piOver2));
    np3 = radius * fast_sin_16(wrap_pi_16(theta));
}

// SIMD batch EM step for 16 trajectories
// 6D state: (I1, I2, I3, phi1, phi2, phi3)
typedef Matrix<float, 6, 16, RowMajor> State6x16;

inline void step_EM_batch16(State6x16& X, RNGState& rng_state)
{
    // Extract rows as arrays for SIMD
    A16f I1 = X.row(0).array();
    A16f I2 = X.row(1).array();
    A16f I3 = X.row(2).array();
    A16f p1 = X.row(3).array();
    A16f p2 = X.row(4).array();
    A16f p3 = X.row(5).array();

    // The two independent phase differences (from which all sin/cos derive)
    A16f dphi12 = 2.0f * (p1 - p2);  // 2*(phi1 - phi2)
    A16f dphi32 = 2.0f * (p3 - p2);  // 2*(phi3 - phi2)

    // Compute sin and cos of the two independent angles
    A16f dphi12_w = wrap_pi_16(dphi12);
    A16f dphi32_w = wrap_pi_16(dphi32);

    A16f sin12 = fast_sin_16(dphi12_w);  // sin(2(phi1-phi2))
    A16f cos12 = fast_cos_16(dphi12_w);  // cos(2(phi1-phi2))
    A16f sin32 = fast_sin_16(dphi32_w);  // sin(2(phi3-phi2))
    A16f cos32 = fast_cos_16(dphi32_w);  // cos(2(phi3-phi2))

    // Total mass
    A16f M = I1 + I2 + I3;

    // === Hamiltonian drift ===
    // dI1 = 4*I1*I2*sin(2(p1-p2))    [I_prev=0 for mode 1]
    A16f dI1 = 4.0f * I1 * I2 * sin12;
    // dI2 = 4*I2*(I1*sin(2(p2-p1)) + I3*sin(2(p2-p3)))
    //     = 4*I2*(-I1*sin12 - I3*sin32)
    A16f dI2 = 4.0f * I2 * (-I1 * sin12 - I3 * sin32);
    // dI3 = 4*I3*I2*sin(2(p3-p2))
    A16f dI3 = 4.0f * I3 * I2 * sin32;

    // dphi1 = 2M - I1 + 2*I2*cos(2(p1-p2))
    A16f dp1 = 2.0f * M - I1 + 2.0f * I2 * cos12;
    // dphi2 = 2M - I2 + 2*I1*cos(2(p2-p1)) + 2*I3*cos(2(p2-p3))
    //       = 2M - I2 + 2*I1*cos12 + 2*I3*cos32
    A16f dp2 = 2.0f * M - I2 + 2.0f * I1 * cos12 + 2.0f * I3 * cos32;
    // dphi3 = 2M - I3 + 2*I2*cos(2(p3-p2))
    A16f dp3 = 2.0f * M - I3 + 2.0f * I2 * cos32;

    // === Heat bath dissipation (mode 1, coupled to T1) ===
    // dI1 += 2*gamma*(2*T1 - (2*M*I1 - I1^2 + 2*I2*I1*cos12))
    dI1 += 2.0f * gamma_val * (2.0f * T1 - (2.0f * M * I1 - I1 * I1 + 2.0f * I2 * I1 * cos12));
    // dphi1 += gamma*(2*I2*sin12)  [CORRECTED SIGN: += not -=]
    dp1 += gamma_val * (2.0f * I2 * sin12);

    // === Heat bath dissipation (mode 3, coupled to T3) ===
    dI3 += 2.0f * gamma_val * (2.0f * T3 - (2.0f * M * I3 - I3 * I3 + 2.0f * I2 * I3 * cos32));
    dp3 += gamma_val * (2.0f * I2 * sin32);

    // === Noise ===
    A16f nI1, nI3, np1, np3;
    generate_noise_4x16(rng_state, nI1, nI3, np1, np3);

    // Diffusion coefficients (vectorized)
    A16f I1c = I1.max(1e-14f);
    A16f I3c = I3.max(1e-14f);
    A16f sigI1 = 2.0f * (2.0f * gamma_val * T1 * I1c).sqrt();
    A16f sigI3 = 2.0f * (2.0f * gamma_val * T3 * I3c).sqrt();
    A16f sigp1 = (2.0f * gamma_val * T1 / I1c).sqrt();
    A16f sigp3 = (2.0f * gamma_val * T3 / I3c).sqrt();

    // === Euler-Maruyama update ===
    X.row(0).array() = (I1 + dI1 * dt + sigI1 * sqrt_dt * nI1).max(1e-14f);
    X.row(1).array() = (I2 + dI2 * dt).max(1e-14f);
    X.row(2).array() = (I3 + dI3 * dt + sigI3 * sqrt_dt * nI3).max(1e-14f);
    X.row(3).array() = p1 + dp1 * dt + sigp1 * sqrt_dt * np1;
    X.row(4).array() = p2 + dp2 * dt;
    X.row(5).array() = p3 + dp3 * dt + sigp3 * sqrt_dt * np3;

    // Wrap angles
    X.row(3) = wrap_pi_16(X.row(3).array());
    X.row(4) = wrap_pi_16(X.row(4).array());
    X.row(5) = wrap_pi_16(X.row(5).array());
}

// =====================================================
// Convert 6D state to 5D coordinates for box lookup
// ====================================================
struct Coords5D { float I1, I2, I3, th1, th3; };

inline Coords5D to_5d_col(const State6x16& X, int col) {
    Coords5D c;
    c.I1  = X(0, col);
    c.I2  = X(1, col);
    c.I3  = X(2, col);
    float p1 = X(3, col), p2 = X(4, col), p3 = X(5, col);
    // theta1 = 2*(phi1-phi2), theta3 = 2*(phi3-phi2)
    float th1 = 2.0f * (p1 - p2);
    float th3 = 2.0f * (p3 - p2);
    // wrap to [-pi, pi]
    th1 = th1 - roundf(th1 * 0.159154943f) * 6.283185307f;
    th3 = th3 - roundf(th3 * 0.159154943f) * 6.283185307f;
    c.th1 = th1;
    c.th3 = th3;
    return c;
}

// =============
// Box lookup 
// ==============
int which_box(vector<vector<int>>& lob, const Coords5D& c) {
    if(c.I1 < I_lo || c.I1 >= I_max) return -1;
    if(c.I2 < I_lo || c.I2 >= I_max) return -1;
    if(c.I3 < I_lo || c.I3 >= I_max) return -1;

    int n1  = (int)floorf((c.I1 - I_lo) / hI);
    int n2  = (int)floorf((c.I2 - I_lo) / hI);
    int n3  = (int)floorf((c.I3 - I_lo) / hI);
    int nt1 = (int)floorf((c.th1 - Th_lo) / hTh);
    int nt3 = (int)floorf((c.th3 - Th_lo) / hTh);

    if(n1 < 0 || n1 >= N_I || n2 < 0 || n2 >= N_I || n3 < 0 || n3 >= N_I) return -1;
    nt1 = max(0, min(nt1, N_Th-1));
    nt3 = max(0, min(nt3, N_Th-1));

    int idxA = n1*N_I + n2;
    int idxB = sizeA + n3*N_Th + nt1;
    int idxC = sizeA + sizeB + nt3;

    if(lob[idxA].empty() || lob[idxB].empty() || lob[idxC].empty()) return -1;

    static thread_local vector<int> ab, abc;
    ab.clear(); abc.clear();
    set_intersection(lob[idxA].begin(), lob[idxA].end(),
                     lob[idxB].begin(), lob[idxB].end(), back_inserter(ab));
    if(ab.empty()) return -1;
    set_intersection(ab.begin(), ab.end(),
                     lob[idxC].begin(), lob[idxC].end(), back_inserter(abc));
    return abc.empty() ? -1 : abc[0];
}

// ===============
// Register box 
// ==============
bool register_box(vector<vector<int>>& lob, const Coords5D& c, int box_id) {
    int n1  = max(0, min((int)floorf((c.I1 - I_lo)/hI),  N_I-1));
    int n2  = max(0, min((int)floorf((c.I2 - I_lo)/hI),  N_I-1));
    int n3  = max(0, min((int)floorf((c.I3 - I_lo)/hI),  N_I-1));
    int nt1 = max(0, min((int)floorf((c.th1 - Th_lo)/hTh), N_Th-1));
    int nt3 = max(0, min((int)floorf((c.th3 - Th_lo)/hTh), N_Th-1));

    int idxA = n1*N_I + n2;
    int idxB = sizeA + n3*N_Th + nt1;
    int idxC = sizeA + sizeB + nt3;

    vector<int> ab;
    set_intersection(lob[idxA].begin(), lob[idxA].end(),
                     lob[idxB].begin(), lob[idxB].end(), back_inserter(ab));
    if(!ab.empty()) {
        vector<int> abc;
        set_intersection(ab.begin(), ab.end(),
                         lob[idxC].begin(), lob[idxC].end(), back_inserter(abc));
        if(!abc.empty()) return false;
    }
    lob[idxA].push_back(box_id); sort(lob[idxA].begin(), lob[idxA].end());
    lob[idxB].push_back(box_id); sort(lob[idxB].begin(), lob[idxB].end());
    lob[idxC].push_back(box_id); sort(lob[idxC].begin(), lob[idxC].end());
    return true;
}

// =================
// Scalar EM step (for create_boxes burn-in, uses double precision)
// =================
typedef Matrix<double, 6, 1> State6d;

inline double wrap_d(double x) {
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}

void step_EM_scalar(State6d& X, mt19937& rng, normal_distribution<double>& dist) {
    double I1=X(0), I2=X(1), I3=X(2), p1=X(3), p2=X(4), p3=X(5);
    double dpn0=2*(p1-p2), dpn1=2*(p2-p3), dpp1=2*(p2-p1), dpp2=2*(p3-p2);
    double M=I1+I2+I3;
    double s12=sin(dpn0), c12=cos(dpn0), s32=sin(dpp2), c32=cos(dpp2);

    double dI1=4*I1*I2*s12, dI2=4*I2*(-I1*s12-I3*s32), dI3=4*I3*I2*s32;
    double dp1=2*M-I1+2*I2*c12, dp2_=2*M-I2+2*I1*c12+2*I3*c32, dp3=2*M-I3+2*I2*c32;

    dI1 += 2*gamma_val*(2*T1-(2*M*I1-I1*I1+2*I2*I1*c12));
    dp1 += gamma_val*(2*I2*s12);
    dI3 += 2*gamma_val*(2*T3-(2*M*I3-I3*I3+2*I2*I3*c32));
    dp3 += gamma_val*(2*I2*s32);

    double dW[4]; for(int j=0;j<4;j++) dW[j]=sqrt(dt)*dist(rng);
    double sI1=2*sqrt(2.0*gamma_val*T1*max(I1,1e-14));
    double sI3=2*sqrt(2.0*gamma_val*T3*max(I3,1e-14));
    double sp1=sqrt(2.0*gamma_val*T1/max(I1,1e-14));
    double sp3=sqrt(2.0*gamma_val*T3/max(I3,1e-14));

    X(0)=max(I1+dI1*dt+sI1*dW[0],1e-14);
    X(1)=max(I2+dI2*dt,1e-14);
    X(2)=max(I3+dI3*dt+sI3*dW[1],1e-14);
    X(3)=wrap_d(p1+dp1*dt+sp1*dW[2]);
    X(4)=wrap_d(p2+dp2_*dt);
    X(5)=wrap_d(p3+dp3*dt+sp3*dW[3]);
}


// Create reference boxes (ratio=0.5: half trajectory, half uniform)
// Uses scalar integrator (double precision) for accuracy
void create_boxes(vector<vector<int>>& lob, MatrixXf& Boxes) {
    State6d X; X << 1.0, 1.0, 0.1, 0.0, 0.0, 0.0;
    mt19937 rng(42);
    normal_distribution<double> dist(0.0, 1.0);
    uniform_real_distribution<float> uni(0.0f, 1.0f);

    cout << "  Burning in..." << endl;
    for(int i = 0; i < Burn_in; i++) step_EM_scalar(X, rng, dist);

    cout << "  Sampling boxes (ratio=" << sample_ratio << ")..." << endl;
    int count = 0, from_traj = 0, from_unif = 0;

    while(count < N_box) {
        Coords5D c;
        if(uni(rng) < sample_ratio) {
            for(int i = 0; i < 5000; i++) step_EM_scalar(X, rng, dist);
            c.I1=(float)X(0); c.I2=(float)X(1); c.I3=(float)X(2);
            c.th1=(float)wrap_d(2*(X(3)-X(4)));
            c.th3=(float)wrap_d(2*(X(5)-X(4)));
            if(c.I1<I_lo||c.I1>=I_max||c.I2<I_lo||c.I2>=I_max||c.I3<I_lo||c.I3>=I_max) continue;
        } else {
            c.I1 = I_lo + uni(rng)*(I_max-I_lo);
            c.I2 = I_lo + uni(rng)*(I_max-I_lo);
            c.I3 = I_lo + uni(rng)*(I_max-I_lo);
            c.th1 = Th_lo + uni(rng)*(Th_hi-Th_lo);
            c.th3 = Th_lo + uni(rng)*(Th_hi-Th_lo);
        }

        if(register_box(lob, c, count)) {
            int n1=max(0,min((int)floorf((c.I1-I_lo)/hI),N_I-1));
            int n2=max(0,min((int)floorf((c.I2-I_lo)/hI),N_I-1));
            int n3=max(0,min((int)floorf((c.I3-I_lo)/hI),N_I-1));
            int nt1=max(0,min((int)floorf((c.th1-Th_lo)/hTh),N_Th-1));
            int nt3=max(0,min((int)floorf((c.th3-Th_lo)/hTh),N_Th-1));
            Boxes(0,count)=I_lo+n1*hI+hI/2; Boxes(1,count)=I_lo+n2*hI+hI/2;
            Boxes(2,count)=I_lo+n3*hI+hI/2;
            Boxes(3,count)=Th_lo+nt1*hTh+hTh/2; Boxes(4,count)=Th_lo+nt3*hTh+hTh/2;
            if(uni(rng)<sample_ratio) from_traj++; else from_unif++;
            count++;
            if(count%10000==0)
                cout<<"    "<<count<<"/"<<N_box<<" (traj="<<from_traj<<" unif="<<from_unif<<")"<<endl;
        }
    }
    cout<<"  Final: "<<from_traj<<" traj, "<<from_unif<<" unif"<<endl;
}

// ====================================================================
// Monte Carlo with SIMD batch processing
// Each thread runs 16 long trajectories in parallel (SIMD lanes).
// N_sample = total steps per trajectory
// Total sample points = N_thread * 16 * N_sample.
// ====================================================================
void MC(VectorXd& Box_count, vector<vector<int>>& lob) {
    long long total_points = (long long)N_thread * 16LL * N_sample;
    cout << "MC SIMD: " << N_thread << " threads x 16 lanes, "
         << N_sample << " steps/traj" << endl;
    cout << "Total sample points = " << total_points << endl;

    vector<long long> counts(N_thread, 0);

#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        RNGState rng_state = initialize_rng_for_rank(rank * 1000 + 42);
        mt19937 mt_init(rank * 137 + 99);
        normal_distribution<float> ndist(0.0f, 1.0f);

        // Initialize 16 trajectories
        State6x16 X;
        for(int col = 0; col < 16; col++) {
            X(0, col) = 1.0f + 0.1f * ndist(mt_init);
            X(1, col) = 1.0f + 0.1f * ndist(mt_init);
            X(2, col) = 0.1f + 0.05f * fabsf(ndist(mt_init));
            X(3, col) = 0.5f * ndist(mt_init);
            X(4, col) = 0.0f;
            X(5, col) = 0.5f * ndist(mt_init);
        }

        // Burn-in
        for(int s = 0; s < Burn_in; s++)
            step_EM_batch16(X, rng_state);

        // Sampling: N_sample steps, each step checks all 16 trajectories
        long long local_hits = 0;
        for(long long step = 0; step < N_sample; step++) {
            step_EM_batch16(X, rng_state);

            for(int col = 0; col < 16; col++) {
                Coords5D c = to_5d_col(X, col);
                int idx = which_box(lob, c);
                if(idx >= 0 && idx < N_box) {
                    Box_count(rank * N_box + idx) += 1.0;
                    local_hits++;
                }
            }

            if(rank == 0 && step % 2000000LL == 0 && step > 0)
                cout << "  thread 0: " << step/1000000 << "M/" 
                     << N_sample/1000000 << "M hits=" << local_hits << endl;
        }
        counts[rank] = local_hits;
    }

    long long total = 0;
    for(int i = 0; i < N_thread; i++) total += counts[i];
    cout << "Total hits = " << total << endl;

    // Merge thread results
    for(int i = 1; i < N_thread; i++)
        for(int j = 0; j < N_box; j++)
            Box_count(j) += Box_count(i * N_box + j);
}

// ==========
// Main
// ==========
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc > 1) N_thread  = atoi(argv[1]);
    if(argc > 2) gamma_val = atof(argv[2]);
    if(argc > 3) T1        = atof(argv[3]);
    if(argc > 4) T3        = atof(argv[4]);
    if(argc > 5) N_box     = atoi(argv[5]);
    if(argc > 6) N_sample  = atoll(argv[6]);

    cout << "=== NLS5D SIMD (16-wide batch) ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "dt=" << dt << " N_box=" << N_box << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << endl;

    vector<vector<int>> lob(total_list);
    for(auto& v : lob) v.reserve(8);
    MatrixXf Boxes(5, N_box);

    cout << "Creating boxes..." << endl;
    create_boxes(lob, Boxes);

    VectorXd Box_count(N_thread * N_box);
    Box_count.fill(0);
    MC(Box_count, lob);

    // Density: total samples = N_thread * 16 * N_sample
    double vol = pow(hI, 3) * pow(hTh, 2);
    long long total_samples = (long long)N_thread * 16LL * N_sample;
    VectorXd density = Box_count.head(N_box) / (total_samples * vol);

    ofstream f1("NLS_FP_boxes.txt"), f2("NLS_FP_density.txt");
    for(int i = 0; i < N_box; i++) {
        f1 << Boxes(0,i) << " " << Boxes(1,i) << " " << Boxes(2,i)
           << " " << Boxes(3,i) << " " << Boxes(4,i) << endl;
        f2 << density(i) << endl;
    }
    f1.close(); f2.close();

    double mx = density.maxCoeff();
    int nz = 0; for(int i=0;i<N_box;i++) if(density(i)>0) nz++;
    cout << "Max density=" << mx << " Nonzero=" << nz << "/" << N_box << endl;

    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec-t1.tv_sec)*1e6 + t2.tv_usec-t1.tv_usec)/1e6;
    cout << "Wall time = " << sec << "s" << endl;
    return 0;
}
