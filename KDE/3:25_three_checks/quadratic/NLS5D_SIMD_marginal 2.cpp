/*
 * NLS5D_SIMD_marginal.cpp
 *
 * Compute marginal distribution of (theta1, theta3) from Monte Carlo,
 * do quadratic fitting: log(PDF) ~ C + A*theta1^2 + B*theta3^2,
 * and compare A*2T with <I1*I2>.
 *
 * No slice restriction — ALL samples contribute to the marginal.
 *
 * Compile (Mac):
 *   clang++ -O3 -std=c++17 -I/opt/homebrew/include/eigen3 \
 *       -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/opt/libomp/lib -lomp NLS5D_SIMD_marginal.cpp -o nls_marginal
 *
 * Usage:
 *   ./nls_marginal [N_thread] [gamma] [T] [N_sample]
 *   e.g.: ./nls_marginal 8 0.1 5.0 200000000
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
#include <vector>
#include <omp.h>

using namespace Eigen;
using namespace std;

// Physical parameters
float gamma_val = 0.1f;
float T_eq      = 5.0f;

// Integrator
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

// Grid: 30x30 on (theta1, theta3)
const float PI_f = 3.14159265358979323846f;
const int   G = 30;
const float H_box = 2.0f * PI_f / G;
const float Th_lo = -PI_f;

// MC parameters
long long N_sample = 200000000LL;
int       N_thread = 8;

// SIMD types
typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

// Fast math
inline A16f wrap_pi_16(const A16f& x) {
    const float invTwoPi = 0.159154943f, twoPi = 6.283185307f;
    return x - (x * invTwoPi).round() * twoPi;
}
inline A16f fast_sin_16(const A16f& x) {
    const float B=1.27323954f, C=-0.40528473f, P=0.225f;
    auto y = B*x + C*x*x.abs();
    return P*(y*y.abs()-y)+y;
}
inline A16f fast_cos_16(const A16f& x) {
    return fast_sin_16(wrap_pi_16(x + 1.570796327f));
}

// RNG
struct RNGState { uint32_t s[4]; };
RNGState init_rng(int rank) {
    RNGState st;
    uint64_t z = (uint64_t)rank + 0x9E3779B97F4A7C15ULL;
    auto sm = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL; uint64_t r = z;
        r = (r^(r>>30))*0xBF58476D1CE4E5B9ULL;
        r = (r^(r>>27))*0x94D049BB133111EBULL;
        return r^(r>>31);
    };
    uint64_t p1=sm(), p2=sm();
    st.s[0]=(uint32_t)p1; st.s[1]=(uint32_t)(p1>>32);
    st.s[2]=(uint32_t)p2; st.s[3]=(uint32_t)(p2>>32);
    return st;
}
inline A16f next_u01_16(RNGState& s) {
    A16f res;
    for(int i=0;i<16;i++){
        uint32_t r=((s.s[0]+s.s[3])<<7)|((s.s[0]+s.s[3])>>25);
        uint32_t t=s.s[1]<<9;
        s.s[2]^=s.s[0]; s.s[3]^=s.s[1];
        s.s[1]^=s.s[2]; s.s[0]^=s.s[3];
        s.s[2]^=t; s.s[3]=(s.s[3]<<11)|(s.s[3]>>21);
        res(i)=(r>>9)*0.00000011920929f;
    }
    return res;
}
inline void gen_noise_4x16(RNGState& s, A16f& n1, A16f& n2, A16f& n3, A16f& n4) {
    const float twoPi=6.283185307f, piO2=1.570796327f;
    A16f u1=next_u01_16(s), u2=next_u01_16(s);
    A16f r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt(), th=twoPi*u2;
    n1=r*fast_sin_16(wrap_pi_16(th+piO2)); n2=r*fast_sin_16(wrap_pi_16(th));
    u1=next_u01_16(s); u2=next_u01_16(s);
    r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt(); th=twoPi*u2;
    n3=r*fast_sin_16(wrap_pi_16(th+piO2)); n4=r*fast_sin_16(wrap_pi_16(th));
}

// EM step
inline void step_EM_batch16(State6x16& X, RNGState& rng) {
    A16f I1=X.row(0).array(), I2=X.row(1).array(), I3=X.row(2).array();
    A16f p1=X.row(3).array(), p2=X.row(4).array(), p3=X.row(5).array();
    A16f d12w=wrap_pi_16(2.0f*(p1-p2)), d32w=wrap_pi_16(2.0f*(p3-p2));
    A16f s12=fast_sin_16(d12w), c12=fast_cos_16(d12w);
    A16f s32=fast_sin_16(d32w), c32=fast_cos_16(d32w);
    A16f M=I1+I2+I3;

    A16f dI1=4.0f*I1*I2*s12;
    A16f dI2=4.0f*I2*(-I1*s12-I3*s32);
    A16f dI3=4.0f*I3*I2*s32;
    A16f dp1=2.0f*M-I1+2.0f*I2*c12;
    A16f dp2=2.0f*M-I2+2.0f*I1*c12+2.0f*I3*c32;
    A16f dp3=2.0f*M-I3+2.0f*I2*c32;

    dI1+=2.0f*gamma_val*(2.0f*T_eq-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    dI3+=2.0f*gamma_val*(2.0f*T_eq-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);

    A16f nI1,nI3,np1,np3;
    gen_noise_4x16(rng, nI1,nI3,np1,np3);
    A16f I1c=I1.max(1e-14f), I3c=I3.max(1e-14f);
    A16f sI1=2.0f*(2.0f*gamma_val*T_eq*I1c).sqrt();
    A16f sI3=2.0f*(2.0f*gamma_val*T_eq*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T_eq/I1c).sqrt();
    A16f sp3=(2.0f*gamma_val*T_eq/I3c).sqrt();

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

inline float wrap_f(float x) {
    return x - roundf(x*0.159154943f)*6.283185307f;
}

int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T_eq      = atof(argv[3]);
    if(argc>4) N_sample  = atoll(argv[4]);

    cout << "=== Marginal Distribution + Quadratic Fit ===" << endl;
    cout << "gamma=" << gamma_val << " T=" << T_eq << endl;
    cout << "Grid: " << G << "x" << G << " H=" << H_box << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << " (x16 SIMD)" << endl;
    cout << endl;

    // Accumulate: 30x30 histogram + <I1*I2> and <I2*I3>
    vector<vector<double>> thread_hist(N_thread, vector<double>(G*G, 0.0));
    vector<double> thread_I1I2(N_thread, 0.0);
    vector<double> thread_I2I3(N_thread, 0.0);
    vector<long long> thread_counts(N_thread, 0);

#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        RNGState rng = init_rng(rank*1000+42);
        mt19937 mt_init(rank*137+99);
        normal_distribution<float> nd(0.0f, 1.0f);

        State6x16 X;
        for(int c=0; c<16; c++){
            X(0,c)=1.0f+0.1f*nd(mt_init);
            X(1,c)=1.0f+0.1f*nd(mt_init);
            X(2,c)=0.1f+0.05f*fabsf(nd(mt_init));
            X(3,c)=0.5f*nd(mt_init); X(4,c)=0; X(5,c)=0.5f*nd(mt_init);
        }
        for(int s=0; s<Burn_in; s++) step_EM_batch16(X, rng);

        long long loc_count = 0;
        double loc_I1I2 = 0, loc_I2I3 = 0;

        for(long long step=0; step<N_sample; step++) {
            step_EM_batch16(X, rng);

            for(int c=0; c<16; c++) {
                float I1=X(0,c), I2=X(1,c), I3=X(2,c);
                float p1=X(3,c), p2=X(4,c), p3=X(5,c);
                float th1 = wrap_f(2.0f*(p1-p2));
                float th3 = wrap_f(2.0f*(p3-p2));

                // Bin into 30x30 grid
                int i1 = (int)floorf((th1 - Th_lo) / H_box);
                int i3 = (int)floorf((th3 - Th_lo) / H_box);
                if(i1 < 0) i1 = 0; if(i1 >= G) i1 = G-1;
                if(i3 < 0) i3 = 0; if(i3 >= G) i3 = G-1;

                thread_hist[rank][i1*G + i3] += 1.0;

                // Accumulate <I1*I2> and <I2*I3>
                loc_I1I2 += (double)(I1 * I2);
                loc_I2I3 += (double)(I2 * I3);
                loc_count++;
            }

            if(rank==0 && step%50000000LL==0 && step>0)
                cout << "  thread 0: " << step/1000000 << "M/" << N_sample/1000000 << "M" << endl;
        }

        thread_I1I2[rank] = loc_I1I2;
        thread_I2I3[rank] = loc_I2I3;
        thread_counts[rank] = loc_count;
    }

    // Merge
    vector<double> hist(G*G, 0.0);
    double total_I1I2 = 0, total_I2I3 = 0;
    long long total_count = 0;
    for(int r=0; r<N_thread; r++) {
        for(int k=0; k<G*G; k++) hist[k] += thread_hist[r][k];
        total_I1I2 += thread_I1I2[r];
        total_I2I3 += thread_I2I3[r];
        total_count += thread_counts[r];
    }

    double mean_I1I2 = total_I1I2 / total_count;
    double mean_I2I3 = total_I2I3 / total_count;

    cout << "Total samples: " << total_count << endl;
    cout << "<I1*I2> = " << mean_I1I2 << endl;
    cout << "<I2*I3> = " << mean_I2I3 << endl;
    cout << endl;

    // Normalize histogram to PDF
    vector<double> pdf(G*G);
    for(int k=0; k<G*G; k++) pdf[k] = hist[k] / (total_count * H_box * H_box);

    // === Quadratic fit: log(PDF) ~ C + A*theta1^2 + B*theta3^2 ===
    // Use bins near center (|theta| < pi/2) for fitting
    // Since density has minimum at center, A and B should be positive
    // Theory predicts A = <I1*I2>/(2T)

    vector<double> y_fit, x1sq_fit, x3sq_fit;
    for(int i=0; i<G; i++) {
        double th1 = Th_lo + (i+0.5)*H_box;
        for(int j=0; j<G; j++) {
            double th3 = Th_lo + (j+0.5)*H_box;
            double p = pdf[i*G+j];
            if(p > 0 && fabs(th1) < PI_f/2 && fabs(th3) < PI_f/2) {
                y_fit.push_back(log(p));
                x1sq_fit.push_back(th1*th1);
                x3sq_fit.push_back(th3*th3);
            }
        }
    }

    int N_fit = y_fit.size();
    cout << "Fitting " << N_fit << " bins with |theta| < pi/2" << endl;

    // Solve: log(pdf) = C + A*th1^2 + B*th3^2
    // Since density has MINIMUM at center, log(pdf) curves upward → A > 0
    // Using normal equations
    MatrixXd X_mat(N_fit, 3);
    VectorXd y_vec(N_fit);
    for(int k=0; k<N_fit; k++) {
        X_mat(k, 0) = 1.0;
        X_mat(k, 1) = x1sq_fit[k];
        X_mat(k, 2) = x3sq_fit[k];
        y_vec(k) = y_fit[k];
    }

    Vector3d beta = (X_mat.transpose() * X_mat).ldlt().solve(X_mat.transpose() * y_vec);
    double C_fit = beta(0);
    double A_fit = beta(1);
    double B_fit = beta(2);

    cout << endl;
    cout << "=== Quadratic Fit Results ===" << endl;
    cout << "log(PDF) ~ C + A*theta1^2 + B*theta3^2" << endl;
    cout << "C = " << C_fit << endl;
    cout << "A = " << A_fit << "  (expect > 0, density minimum at center)" << endl;
    cout << "B = " << B_fit << endl;
    cout << endl;

    // Theory: log rho ~ const + <I1I2>/(2T) * theta1^2 + <I2I3>/(2T) * theta3^2
    // So A should equal <I1I2>/(2T), i.e. A*2T = <I1I2>
    cout << "=== Comparison ===" << endl;
    cout << "Theory predicts: A = <I1*I2>/(2T), so A*2T = <I1*I2>" << endl;
    cout << endl;
    cout << "A * 2T  = " << A_fit * 2 * T_eq << endl;
    cout << "<I1*I2> = " << mean_I1I2 << endl;
    cout << "Ratio A*2T / <I1*I2> = " << (A_fit * 2 * T_eq) / mean_I1I2 << endl;
    cout << endl;
    cout << "B * 2T  = " << B_fit * 2 * T_eq << endl;
    cout << "<I2*I3> = " << mean_I2I3 << endl;
    cout << "Ratio B*2T / <I2*I3> = " << (B_fit * 2 * T_eq) / mean_I2I3 << endl;

    // Write marginal PDF to CSV
    {
        ofstream f("marginal_pdf.csv");
        f << "theta1,theta3,pdf,log_pdf\n";
        for(int i=0; i<G; i++) {
            double th1 = Th_lo + (i+0.5)*H_box;
            for(int j=0; j<G; j++) {
                double th3 = Th_lo + (j+0.5)*H_box;
                double p = pdf[i*G+j];
                f << th1 << "," << th3 << "," << p << "," << (p>0?log(p):0) << "\n";
            }
        }
    }

    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << endl << "Wall time = " << sec << "s" << endl;
    cout << "Output: marginal_pdf.csv" << endl;

    return 0;
}
