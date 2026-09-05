/*
 * NLS5D_SIMD_nckde_Tscan.cpp
 *
 * Based on NLS5D_SIMD_nckde.cpp (advisor's NC-KDE code).
 * Modification: sweep T_run, run MC+NCKDE at each T_run,
 *   compare with fixed theory(T_target), output RMSE vs T_run.
 *
 * If T* correction applies to NC-KDE, RMSE minimum is at T* < T_target.
 * If NC-KDE has no smoothing bias, RMSE minimum is at T_target.
 *
 * Compile (Mac):
 *   clang++ -O3 -std=c++17 -I/opt/homebrew/include/eigen3 \
 *       -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/opt/libomp/lib -lomp NLS5D_SIMD_nckde_Tscan.cpp -o nckde_tscan
 *
 * Usage:
 *   ./nckde_tscan [N_thread] [gamma] [T_target] [N_sample] [I_slice] [I_tol]
 *   e.g.: ./nckde_tscan 8 0.1 5.0 30000000 2.0 0.5
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

// ====================================================================
// Physical parameters
// ====================================================================
float gamma_val = 0.1f;
float T_run_global = 5.0f;   // will be set before each MC run

// ====================================================================
// Integrator
// ====================================================================
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 1000000;

// ====================================================================
// Slice parameters
// ====================================================================
float I_slice = 2.0f;
float I_tol   = 0.5f;

// ====================================================================
// Grid: 30x30 on (theta1, theta3)
// ====================================================================
const float PI_f = 3.14159265358979323846f;
const int   G = 30;
const float H_box = 2.0f * PI_f / G;
const float sigma_kde = H_box / 3.0f;   // NC-KDE bandwidth (same as advisor's code)
const float Th_lo = -PI_f;

// ====================================================================
// MC parameters
// ====================================================================
long long N_sample = 30000000LL;
int       N_thread = 8;

// ====================================================================
// SIMD types and fast math (unchanged from advisor's code)
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
// RNG (unchanged)
// ====================================================================
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

// ====================================================================
// SIMD EM step — uses T_run_global (changed from T_eq)
// ====================================================================
inline void step_EM_batch16(State6x16& X, RNGState& rng) {
    float T = T_run_global;
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
    dI1+=2.0f*gamma_val*(2.0f*T-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    dI3+=2.0f*gamma_val*(2.0f*T-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);
    A16f nI1,nI3,np1,np3;
    gen_noise_4x16(rng, nI1,nI3,np1,np3);
    A16f I1c=I1.max(1e-14f), I3c=I3.max(1e-14f);
    A16f sI1=2.0f*(2.0f*gamma_val*T*I1c).sqrt();
    A16f sI3=2.0f*(2.0f*gamma_val*T*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T/I1c).sqrt();
    A16f sp3=(2.0f*gamma_val*T/I3c).sqrt();
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
// Helpers (unchanged from advisor's code)
// ====================================================================
inline bool in_slice(const State6x16& X, int col) {
    return fabsf(X(0,col)-I_slice)<I_tol
        && fabsf(X(1,col)-I_slice)<I_tol
        && fabsf(X(2,col)-I_slice)<I_tol;
}
inline float wrap_f(float x) {
    return x - roundf(x*0.159154943f)*6.283185307f;
}
inline double normpdf(double x, double mu, double sig) {
    double d = x - mu;
    return exp(-0.5*d*d/(sig*sig)) / (sig * sqrt(2.0*M_PI));
}
inline double normpdf2d(double x1, double x2, double mu1, double mu2, double sig) {
    return normpdf(x1, mu1, sig) * normpdf(x2, mu2, sig);
}
double bessel_I0(double x) {
    double s=1,t=1;
    for(int k=1;k<200;k++){t*=(x/(2.0*k))*(x/(2.0*k));s+=t;if(t<1e-16*s)break;}
    return s;
}

// ====================================================================
// Run MC + NC-KDE at given T_run, return normalized density on 30x30 grid
// ====================================================================
void run_nckde(float T_run, vector<double>& mc_density, long long& total_sl) {
    T_run_global = T_run;

    vector<vector<double>> thread_grids(N_thread, vector<double>(G*G, 0.0));
    vector<long long> thread_counts(N_thread, 0);

#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        RNGState rng = init_rng(rank*1000+42+(int)(T_run*1000));
        mt19937 mt_init(rank*137+99+(int)(T_run*1000));
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

                // === NC-KDE: same as advisor's code ===
                int i1 = (int)floorf((th1 - Th_lo) / H_box);
                int i3 = (int)floorf((th3 - Th_lo) / H_box);
                if(i1 < 0) i1 = 0; if(i1 >= G) i1 = G-1;
                if(i3 < 0) i3 = 0; if(i3 >= G) i3 = G-1;
                double tc1 = Th_lo + (i1 + 0.5) * H_box;
                double tc3 = Th_lo + (i3 + 0.5) * H_box;
                double w = normpdf2d((double)th1, (double)th3, tc1, tc3, sigma_kde);
                thread_grids[rank][i1*G + i3] += w;
                loc_count++;
            }
        }
        thread_counts[rank] = loc_count;
    }

    mc_density.assign(G*G, 0.0);
    total_sl = 0;
    for(int r=0; r<N_thread; r++) {
        for(int k=0; k<G*G; k++) mc_density[k] += thread_grids[r][k];
        total_sl += thread_counts[r];
    }
    // Normalize: sum * H^2 = 1
    double s=0;
    for(int k=0;k<G*G;k++) s+=mc_density[k];
    for(int k=0;k<G*G;k++) mc_density[k]/=(s*H_box*H_box);
}

// ====================================================================
// Main
// ====================================================================
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    float T_target = 5.0f;
    double I1I2_measured = 0;
    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T_target  = atof(argv[3]);
    if(argc>4) N_sample  = atoll(argv[4]);
    if(argc>5) I_slice   = atof(argv[5]);
    if(argc>6) I_tol     = atof(argv[6]);

    // Step 0: Run MC at T_target to measure actual <I1I2>
    cout << "=== NC-KDE T* Verification ===" << endl;
    cout << "Step 0: Measuring <I1I2> from MC at T_target..." << endl;
    {
        T_run_global = T_target;
        // Quick run to measure <I1I2>
        vector<vector<double>> thr_I1I2(N_thread, vector<double>(1, 0.0));
        vector<long long> thr_cnt(N_thread, 0);
#pragma omp parallel num_threads(N_thread)
        {
            int rank = omp_get_thread_num();
            RNGState rng = init_rng(rank*1000+7777);
            mt19937 mt_init(rank*137+7777);
            normal_distribution<float> nd(0.0f, 1.0f);
            State6x16 X;
            for(int c=0;c<16;c++){
                X(0,c)=1.0f+0.1f*nd(mt_init); X(1,c)=1.0f+0.1f*nd(mt_init);
                X(2,c)=0.1f+0.05f*fabsf(nd(mt_init));
                X(3,c)=0.5f*nd(mt_init); X(4,c)=0; X(5,c)=0.5f*nd(mt_init);
            }
            for(int s=0;s<Burn_in;s++) step_EM_batch16(X,rng);
            for(long long step=0;step<N_sample;step++){
                step_EM_batch16(X,rng);
                for(int c=0;c<16;c++){
                    if(!in_slice(X,c)) continue;
                    thr_I1I2[rank][0] += (double)X(0,c)*(double)X(1,c);
                    thr_cnt[rank]++;
                }
            }
        }
        double sum_I1I2=0; long long sum_cnt=0;
        for(int r=0;r<N_thread;r++){sum_I1I2+=thr_I1I2[r][0]; sum_cnt+=thr_cnt[r];}
        I1I2_measured = sum_I1I2/sum_cnt;
        cout << "  <I1*I2> = " << I1I2_measured << " (from " << sum_cnt << " samples)" << endl;
        cout << "  (vs I_slice^2 = " << I_slice*I_slice << ")" << endl;
    }

    double a_target = I1I2_measured / T_target;
    cout << "Sweep T_run, run advisor's NC-KDE, compare with theory(T_target)." << endl;
    cout << "gamma=" << gamma_val << " T_target=" << T_target << endl;
    cout << "I_slice=" << I_slice << " +/- " << I_tol << endl;
    cout << "sigma_kde=" << sigma_kde << " (H_box/3)" << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << " (x16 SIMD)" << endl;
    cout << endl;

    // Theory density at T_target (fixed reference)
    vector<double> theory(G*G);
    double I0a = bessel_I0(a_target);
    double Z = (2*M_PI*I0a) * (2*M_PI*I0a);
    for(int i=0; i<G; i++) {
        double tc1 = Th_lo + (i+0.5)*H_box;
        for(int j=0; j<G; j++) {
            double tc3 = Th_lo + (j+0.5)*H_box;
            theory[i*G+j] = exp(-a_target*(cos(tc1)+cos(tc3))) / Z;
        }
    }
    double th_sum=0;
    for(int k=0;k<G*G;k++) th_sum+=theory[k];
    for(int k=0;k<G*G;k++) theory[k]/=(th_sum*H_box*H_box);

    // T* prediction
    double T_star = T_target - sigma_kde*sigma_kde*I1I2_measured;
    cout << "T* = T_target - sigma^2*<I1I2> = " << T_star << endl;
    cout << endl;

    // T_run scan
    int N_T = 15;
    double T_lo = T_star - 0.6;
    double T_hi = T_target + 0.6;
    if(T_lo < 1.0) T_lo = 1.0;

    ofstream fout("nckde_tstar_curves.csv");
    fout << "T_run,rmse" << endl;

    printf("%-10s  %-10s  %-8s  %s\n", "T_run", "RMSE", "samples", "");
    printf("----------  ----------  --------  ----\n");

    double best_rmse = 1e30, best_T = T_target;

    for(int it=0; it<N_T; it++) {
        double T_run = T_lo + (T_hi - T_lo) * it / (N_T - 1);

        vector<double> mc_dens;
        long long nsl;
        run_nckde((float)T_run, mc_dens, nsl);

        double rmse = 0;
        for(int k=0;k<G*G;k++){double d=mc_dens[k]-theory[k]; rmse+=d*d;}
        rmse = sqrt(rmse/(G*G));

        fout << T_run << "," << rmse << endl;

        string note = "";
        if(fabs(T_run - T_target) < (T_hi-T_lo)/(N_T-1)*0.6) note = "<-- T_target";
        if(fabs(T_run - T_star) < (T_hi-T_lo)/(N_T-1)*0.6)   note = "<-- T*";

        printf("%-10.4f  %-10.6f  %-8lld  %s\n", T_run, rmse, nsl, note.c_str());

        if(rmse < best_rmse) { best_rmse = rmse; best_T = T_run; }
    }

    printf("\n>> BEST: T_run=%.4f  RMSE=%.6f\n", best_T, best_rmse);
    printf(">> T*=%.4f  T_target=%.4f\n", T_star, T_target);

    fout.close();
    cout << endl << "Output: nckde_tstar_curves.csv" << endl;

    gettimeofday(&t2, NULL);
    double sec=((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << "Wall time = " << sec << "s" << endl;
    return 0;
}
