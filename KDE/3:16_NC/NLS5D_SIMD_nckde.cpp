
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

// ====================================================================
// Physical parameters
// ====================================================================
float gamma_val = 0.1f;
float T_eq      = 5.0f;   // T1 = T3 = T (equilibrium)

// ====================================================================
// Integrator
// ====================================================================
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

// ====================================================================
// Slice parameters
// ====================================================================
float I_slice = 2.0f;
float I_tol   = 0.5f;

// ====================================================================
// Grid: 30x30 on (theta1, theta3), box size H = 2pi/30
// ====================================================================
const float PI_f = 3.14159265358979323846f;
const int   G = 30;
const float H_box = 2.0f * PI_f / G;        // box edge length
const float sigma_kde = H_box / 3.0f;        // KDE bandwidth
const float Th_lo = -PI_f;

// ====================================================================
// MC parameters
// ====================================================================
long long N_sample = 200000000LL;
int       N_thread = 8;

// ====================================================================
// SIMD types and fast math (same as NLS5D_SIMD.cpp)
// ====================================================================
typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

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

// ====================================================================
// RNG
// ====================================================================
struct RNGState { uint32_t s[4]; };

RNGState init_rng(int rank) {
    RNGState st;
    uint64_t z = (uint64_t)rank + 0x9E3779B97F4A7C15ULL;
    auto sm = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t r = z;
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

// ====================================================================
// SIMD EM step
// ====================================================================
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

// ====================================================================
// Helpers
// ====================================================================
inline bool in_slice(const State6x16& X, int col) {
    return fabsf(X(0,col)-I_slice)<I_tol
        && fabsf(X(1,col)-I_slice)<I_tol
        && fabsf(X(2,col)-I_slice)<I_tol;
}

inline float wrap_f(float x) {
    return x - roundf(x*0.159154943f)*6.283185307f;
}

// 1D Gaussian PDF
inline double normpdf(double x, double mu, double sig) {
    double d = x - mu;
    return exp(-0.5*d*d/(sig*sig)) / (sig * sqrt(2.0*M_PI));
}

// 2D Gaussian PDF (product of two 1D, since theta1 and theta3 independent)
inline double normpdf2d(double x1, double x2, double mu1, double mu2, double sig) {
    return normpdf(x1, mu1, sig) * normpdf(x2, mu2, sig);
}

double bessel_I0(double x) {
    double s=1,t=1;
    for(int k=1;k<200;k++){t*=(x/(2.0*k))*(x/(2.0*k));s+=t;if(t<1e-16*s)break;}
    return s;
}

// ====================================================================
// Main
// ====================================================================
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T_eq      = atof(argv[3]);
    if(argc>4) N_sample  = atoll(argv[4]);
    if(argc>5) I_slice   = atof(argv[5]);
    if(argc>6) I_tol     = atof(argv[6]);

    double a = (double)I_slice * I_slice / T_eq;

    cout << "=== Noise-Compensating KDE Slice Test ===" << endl;
    cout << "gamma=" << gamma_val << " T=" << T_eq << " I_slice=" << I_slice << endl;
    cout << "Grid: " << G << "x" << G << " H=" << H_box << " sigma=" << sigma_kde << endl;
    cout << "a = I^2/T = " << a << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << " (x16 SIMD)" << endl;
    cout << endl;

    // === Theoretical density on 30x30 grid ===
    vector<double> theory(G*G);
    double I0a = bessel_I0(a);
    double Z_theory = (2*M_PI*I0a) * (2*M_PI*I0a);
    for(int i=0; i<G; i++) {
        double tc1 = Th_lo + (i+0.5)*H_box;
        for(int j=0; j<G; j++) {
            double tc3 = Th_lo + (j+0.5)*H_box;
            theory[i*G+j] = exp(-a*(cos(tc1)+cos(tc3))) / Z_theory;
        }
    }

    // === MC with noise-compensating KDE ===
    // Each thread accumulates into its own 30x30 grid
    vector<vector<double>> thread_grids(N_thread, vector<double>(G*G, 0.0));
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

        for(long long step=0; step<N_sample; step++) {
            step_EM_batch16(X, rng);

            for(int c=0; c<16; c++) {
                if(!in_slice(X, c)) continue;

                float p1=X(3,c), p2=X(4,c), p3=X(5,c);
                float th1 = wrap_f(2.0f*(p1-p2));
                float th3 = wrap_f(2.0f*(p3-p2));

                // Which box does this sample fall in?
                int i1 = (int)floorf((th1 - Th_lo) / H_box);
                int i3 = (int)floorf((th3 - Th_lo) / H_box);
                if(i1 < 0) i1 = 0; if(i1 >= G) i1 = G-1;
                if(i3 < 0) i3 = 0; if(i3 >= G) i3 = G-1;

                // Box center
                double tc1 = Th_lo + (i1 + 0.5) * H_box;
                double tc3 = Th_lo + (i3 + 0.5) * H_box;

                // Noise-compensating KDE weight
                double w = normpdf2d((double)th1, (double)th3, tc1, tc3, sigma_kde);

                thread_grids[rank][i1*G + i3] += w;
                loc_count++;
            }

            if(rank==0 && step%50000000LL==0 && step>0)
                cout << "  thread 0: " << step/1000000 << "M/" << N_sample/1000000 << "M" << endl;
        }
        thread_counts[rank] = loc_count;
    }

    // Merge threads
    vector<double> mc_density(G*G, 0.0);
    long long total_sl = 0;
    for(int r=0; r<N_thread; r++) {
        for(int k=0; k<G*G; k++) mc_density[k] += thread_grids[r][k];
        total_sl += thread_counts[r];
    }

    cout << "Total slice samples: " << total_sl << endl;

    // Normalize MC: sum * H^2 = 1
    double mc_sum = 0;
    for(int k=0; k<G*G; k++) mc_sum += mc_density[k];
    double mc_norm = mc_sum * H_box * H_box;
    for(int k=0; k<G*G; k++) mc_density[k] /= mc_norm;

    // Normalize theory: sum * H^2 = 1
    double th_sum = 0;
    for(int k=0; k<G*G; k++) th_sum += theory[k];
    double th_norm = th_sum * H_box * H_box;
    for(int k=0; k<G*G; k++) theory[k] /= th_norm;

    // === Compare ===
    double rmse = 0, max_err = 0;
    for(int k=0; k<G*G; k++) {
        double d = fabs(mc_density[k] - theory[k]);
        rmse += d*d;
        if(d > max_err) max_err = d;
    }
    rmse = sqrt(rmse / (G*G));

    cout << endl;
    cout << "RMSE  = " << rmse << endl;
    cout << "MaxErr = " << max_err << endl;

    // Check corners vs center
    // Corners: (0,0), (0,G-1), (G-1,0), (G-1,G-1) → theta ≈ ±pi → high density
    // Center: (G/2, G/2) → theta ≈ 0 → low density
    cout << endl << "Corner/center comparison:" << endl;
    printf("  %-20s %10s %10s %10s\n", "Location", "Theory", "MC", "Ratio");
    
    int corners[][2] = {{0,0},{0,G-1},{G-1,0},{G-1,G-1}};
    for(auto& c : corners) {
        int k = c[0]*G + c[1];
        double th1 = Th_lo + (c[0]+0.5)*H_box;
        double th3 = Th_lo + (c[1]+0.5)*H_box;
        printf("  (%.2f, %.2f)       %10.6f %10.6f %10.4f\n",
               th1, th3, theory[k], mc_density[k], mc_density[k]/theory[k]);
    }
    {
        int k = (G/2)*G + G/2;
        double th1 = Th_lo + (G/2+0.5)*H_box;
        double th3 = Th_lo + (G/2+0.5)*H_box;
        printf("  (%.2f, %.2f) center %10.6f %10.6f %10.4f\n",
               th1, th3, theory[k], mc_density[k], mc_density[k]/theory[k]);
    }

    // Cross-section at theta3 ~ 0
    cout << endl << "Cross-section theta3~0:" << endl;
    printf("  %8s %10s %10s %10s\n", "theta1", "Theory", "MC_NCKDE", "Ratio");
    int j3m = G/2;
    for(int i=0; i<G; i++) {
        double t1 = Th_lo + (i+0.5)*H_box;
        int k = i*G + j3m;
        printf("  %8.3f %10.6f %10.6f %10.4f\n",
               t1, theory[k], mc_density[k], 
               theory[k]>0 ? mc_density[k]/theory[k] : 0.0);
    }

    // === Write CSV for MATLAB ===
    {
        ofstream f("nckde_mc.csv");
        f << "theta1,theta3,density\n";
        for(int i=0; i<G; i++) {
            double t1 = Th_lo + (i+0.5)*H_box;
            for(int j=0; j<G; j++) {
                double t3 = Th_lo + (j+0.5)*H_box;
                f << t1 << "," << t3 << "," << mc_density[i*G+j] << "\n";
            }
        }
    }
    {
        ofstream f("nckde_theory.csv");
        f << "theta1,theta3,density\n";
        for(int i=0; i<G; i++) {
            double t1 = Th_lo + (i+0.5)*H_box;
            for(int j=0; j<G; j++) {
                double t3 = Th_lo + (j+0.5)*H_box;
                f << t1 << "," << t3 << "," << theory[i*G+j] << "\n";
            }
        }
    }
    {
        ofstream f("nckde_diff.csv");
        f << "theta1,theta3,diff\n";
        for(int i=0; i<G; i++) {
            double t1 = Th_lo + (i+0.5)*H_box;
            for(int j=0; j<G; j++) {
                double t3 = Th_lo + (j+0.5)*H_box;
                f << t1 << "," << t3 << "," << mc_density[i*G+j]-theory[i*G+j] << "\n";
            }
        }
    }

    gettimeofday(&t2, NULL);
    double sec = ((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << endl << "Wall time = " << sec << "s" << endl;

    cout << endl << "Output files: nckde_mc.csv, nckde_theory.csv, nckde_diff.csv" << endl;
    cout << "MATLAB: mesh(reshape(readtable('nckde_diff.csv').diff, 30, 30))" << endl;

    return 0;
}
