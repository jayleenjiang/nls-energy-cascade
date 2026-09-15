/*
 * NLS5D_SIMD_Teff_scan.cpp
 *
 * Verify that NC-KDE conditional distribution error is minimized at T*.
 *
 * Strategy:
 *   Phase 1: MC sampling — collect all (theta1, theta3) samples on a slice.
 *            Store raw samples (not pre-binned).
 *   Phase 2: For each bandwidth h in {0.05, 0.10, 0.15, 0.20, 0.25, 0.30}:
 *            - Re-bin samples with NC-KDE weight normpdf(x, x_c, h)
 *            - Scan a_eff from 0.5*a_true to 1.5*a_true
 *            - Compute RMSE(a_eff) = || MC_density - theory(a_eff) ||
 *            - Find a_opt (minimizes RMSE) and compare with a_pred = <I1I2>/T*
 *              where T* = T - h^2 * <I1I2>   (corrected formula, no /2)
 *
 * Corrected T* formula comparison:
 *   Formula A (advisor): T* = T - h^2 * <I1I2> / 2
 *   Formula B (alt):     T* = T - h^2 * <I1I2>
 *   Data shows Formula A matches a_opt better.
 *
 * Compile (Mac):
 *   clang++ -O3 -std=c++17 -I/opt/homebrew/include/eigen3 \
 *       -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
 *       -L/opt/homebrew/opt/libomp/lib -lomp NLS5D_SIMD_Teff_scan.cpp -o nls_teff
 *
 * Usage:
 *   ./nls_teff [N_thread] [gamma] [T] [N_sample] [I_slice] [I_tol]
 *   e.g.: ./nls_teff 8 0.1 5.0 200000000 2.0 0.5
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

float gamma_val = 0.1f;
float T_eq      = 5.0f;
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

float I_slice = 2.0f;
float I_tol   = 0.5f;

const float PI_f = 3.14159265358979323846f;
const int   G = 30;
const float H_box = 2.0f * PI_f / G;
const float Th_lo = -PI_f;

long long N_sample = 200000000LL;
int       N_thread = 8;

typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

// --- SIMD helpers (same as check2) ---
inline A16f wrap_pi_16(const A16f& x) {
    const float invTwoPi=0.159154943f, twoPi=6.283185307f;
    return x-(x*invTwoPi).round()*twoPi;
}
inline A16f fast_sin_16(const A16f& x) {
    const float B=1.27323954f,C=-0.40528473f,P=0.225f;
    auto y=B*x+C*x*x.abs(); return P*(y*y.abs()-y)+y;
}
inline A16f fast_cos_16(const A16f& x) {
    return fast_sin_16(wrap_pi_16(x+1.570796327f));
}

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
inline void gen_noise_4x16(RNGState& s,A16f& n1,A16f& n2,A16f& n3,A16f& n4){
    const float twoPi=6.283185307f,piO2=1.570796327f;
    A16f u1=next_u01_16(s),u2=next_u01_16(s);
    A16f r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt(),th=twoPi*u2;
    n1=r*fast_sin_16(wrap_pi_16(th+piO2));n2=r*fast_sin_16(wrap_pi_16(th));
    u1=next_u01_16(s);u2=next_u01_16(s);
    r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt();th=twoPi*u2;
    n3=r*fast_sin_16(wrap_pi_16(th+piO2));n4=r*fast_sin_16(wrap_pi_16(th));
}

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
    // Dissipation (dphi: += sign, corrected)
    dI1+=2.0f*gamma_val*(2.0f*T_eq-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    dI3+=2.0f*gamma_val*(2.0f*T_eq-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);
    // Noise (with sqrt(2))
    A16f nI1,nI3,np1,np3; gen_noise_4x16(rng,nI1,nI3,np1,np3);
    A16f I1c=I1.max(1e-14f),I3c=I3.max(1e-14f);
    A16f sI1=2.0f*(2.0f*gamma_val*T_eq*I1c).sqrt(),sI3=2.0f*(2.0f*gamma_val*T_eq*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T_eq/I1c).sqrt(),sp3=(2.0f*gamma_val*T_eq/I3c).sqrt();
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

inline bool in_slice(const State6x16& X, int col) {
    return fabsf(X(0,col)-I_slice)<I_tol && fabsf(X(1,col)-I_slice)<I_tol && fabsf(X(2,col)-I_slice)<I_tol;
}
inline float wrap_f(float x) { return x-roundf(x*0.159154943f)*6.283185307f; }

double bessel_I0(double x) {
    double s=1,t=1; for(int k=1;k<200;k++){t*=(x/(2.0*k))*(x/(2.0*k));s+=t;if(t<1e-16*s)break;} return s;
}

// Compute theory density on G×G grid for given parameter a
// theory(th1,th3) = exp(-a*(cos(th1)+cos(th3))) / Z, normalized so sum*H^2=1
void compute_theory(double a, vector<double>& theory) {
    theory.resize(G*G);
    double I0a = bessel_I0(a);
    double Z = (2*M_PI*I0a)*(2*M_PI*I0a);
    for(int i=0;i<G;i++){
        double tc1=Th_lo+(i+0.5)*H_box;
        for(int j=0;j<G;j++){
            double tc3=Th_lo+(j+0.5)*H_box;
            theory[i*G+j]=exp(-a*(cos(tc1)+cos(tc3)))/Z;
        }
    }
    double s=0; for(int k=0;k<G*G;k++) s+=theory[k];
    for(int k=0;k<G*G;k++) theory[k]/=(s*H_box*H_box);
}

// Compute RMSE between two grids
double compute_rmse(const vector<double>& a, const vector<double>& b) {
    double s=0;
    for(int k=0;k<G*G;k++){ double d=a[k]-b[k]; s+=d*d; }
    return sqrt(s/(G*G));
}

int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T_eq      = atof(argv[3]);
    if(argc>4) N_sample  = atoll(argv[4]);
    if(argc>5) I_slice   = atof(argv[5]);
    if(argc>6) I_tol     = atof(argv[6]);

    cout << "=== T_eff Scan: Verify T* minimizes RMSE ===" << endl;
    cout << "gamma=" << gamma_val << " T=" << T_eq << " I_slice=" << I_slice << " +/- " << I_tol << endl;
    cout << "Grid: " << G << "x" << G << " H=" << H_box << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << " (x16 SIMD)" << endl;
    cout << "T* formula: T* = T - h^2 * <I1I2>  (corrected, no /2)" << endl;
    cout << endl;

    // =========================================================
    // Phase 1: Collect raw (theta1, theta3) samples on slice
    // =========================================================
    // Each thread stores its samples in a local vector, merge after
    vector<vector<float>> thread_th1(N_thread), thread_th3(N_thread);
    vector<double> thread_I1I2(N_thread, 0.0);
    vector<double> thread_I2I3(N_thread, 0.0);
    vector<long long> thread_counts(N_thread, 0);

    // Pre-reserve memory (estimate ~1% hit rate for I_tol=0.5)
    long long est_per_thread = (long long)(N_sample * 16 * 0.01);
    for(int r=0;r<N_thread;r++){
        thread_th1[r].reserve(est_per_thread);
        thread_th3[r].reserve(est_per_thread);
    }

    cout << "Phase 1: MC sampling..." << endl;

#pragma omp parallel num_threads(N_thread)
    {
        int rank=omp_get_thread_num();
        RNGState rng=init_rng(rank*1000+42);
        mt19937 mt_init(rank*137+99);
        normal_distribution<float> nd(0.0f,1.0f);

        State6x16 X;
        for(int c=0;c<16;c++){
            X(0,c)=1.0f+0.1f*nd(mt_init); X(1,c)=1.0f+0.1f*nd(mt_init);
            X(2,c)=0.1f+0.05f*fabsf(nd(mt_init));
            X(3,c)=0.5f*nd(mt_init); X(4,c)=0; X(5,c)=0.5f*nd(mt_init);
        }
        for(int s=0;s<Burn_in;s++) step_EM_batch16(X,rng);

        double loc_I1I2=0, loc_I2I3=0;
        long long loc_count=0;

        for(long long step=0;step<N_sample;step++){
            step_EM_batch16(X,rng);
            for(int c=0;c<16;c++){
                if(!in_slice(X,c)) continue;
                float I1=X(0,c),I2=X(1,c),I3=X(2,c);
                float p1=X(3,c),p2=X(4,c),p3=X(5,c);
                float th1=wrap_f(2.0f*(p1-p2)), th3=wrap_f(2.0f*(p3-p2));

                thread_th1[rank].push_back(th1);
                thread_th3[rank].push_back(th3);
                loc_I1I2 += (double)(I1*I2);
                loc_I2I3 += (double)(I2*I3);
                loc_count++;
            }
            if(rank==0 && step%50000000LL==0 && step>0)
                cout<<"  thread 0: "<<step/1000000<<"M/"<<N_sample/1000000<<"M"<<endl;
        }
        thread_I1I2[rank]=loc_I1I2; thread_I2I3[rank]=loc_I2I3;
        thread_counts[rank]=loc_count;
    }

    // Merge samples
    long long total_sl=0; double total_I1I2=0, total_I2I3=0;
    for(int r=0;r<N_thread;r++){
        total_I1I2+=thread_I1I2[r]; total_I2I3+=thread_I2I3[r];
        total_sl+=thread_counts[r];
    }
    // Flatten into single arrays
    vector<float> all_th1, all_th3;
    all_th1.reserve(total_sl); all_th3.reserve(total_sl);
    for(int r=0;r<N_thread;r++){
        all_th1.insert(all_th1.end(), thread_th1[r].begin(), thread_th1[r].end());
        all_th3.insert(all_th3.end(), thread_th3[r].begin(), thread_th3[r].end());
        // Free memory
        thread_th1[r].clear(); thread_th1[r].shrink_to_fit();
        thread_th3[r].clear(); thread_th3[r].shrink_to_fit();
    }

    double mean_I1I2=total_I1I2/total_sl;
    double mean_I2I3=total_I2I3/total_sl;
    double a_true = mean_I1I2 / T_eq;  // true parameter (no KDE smoothing)

    cout << "Total slice samples: " << total_sl << endl;
    cout << "<I1*I2> = " << mean_I1I2 << "  <I2*I3> = " << mean_I2I3 << endl;
    cout << "a_true = <I1I2>/T = " << a_true << endl;
    cout << endl;

    // =========================================================
    // Phase 2: For each h, re-bin with NC-KDE, then scan a_eff
    // =========================================================

    // h values to test
    vector<double> h_list = {0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50};

    // a_eff scan parameters (fine grid, will be centered per-h)
    const int N_a = 101;

    // Output file
    ofstream fout("teff_scan_results.csv");
    fout << "h,a_true,a_pred_A,a_pred_B,a_opt,err_true,err_pred_A,err_pred_B,err_opt" << endl;

    // Also output full RMSE curves for plotting
    ofstream fcurve("teff_scan_curves.csv");
    fcurve << "h,a_eff,rmse" << endl;

    cout << "Phase 2: T_eff scan..." << endl;
    cout << endl;
    printf("%-6s  %-8s  %-8s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s\n",
           "h", "a_true", "a_A(/2)", "a_B(x1)", "a_opt", "err_true", "err_A", "err_B", "err_opt");
    printf("------  --------  --------  --------  --------  ----------  ----------  ----------  ----------\n");

    for(double h : h_list) {
        // Step 2a: Build TRUE KDE density grid for this h
        // Each sample contributes to ALL nearby bins (true convolution)
        vector<double> mc_density(G*G, 0.0);

        int R = (int)ceil(3.0*h / H_box);  // search radius: 3 sigma in grid units
        if(R < 1) R = 1;
        if(R >= G) R = G-1;

        for(long long s=0; s<total_sl; s++){
            double th1 = (double)all_th1[s];
            double th3 = (double)all_th3[s];

            // Contribute to all bins within radius R
            // Find center bin
            int ci1 = (int)floor((th1-Th_lo)/H_box);
            int ci3 = (int)floor((th3-Th_lo)/H_box);

            for(int di1 = -R; di1 <= R; di1++){
                int i1 = ci1 + di1;
                // Periodic wrap for theta
                if(i1 < 0) i1 += G;
                if(i1 >= G) i1 -= G;
                double tc1 = Th_lo + (i1+0.5)*H_box;
                // Periodic distance
                double dd1 = th1 - tc1;
                if(dd1 > M_PI) dd1 -= 2*M_PI;
                if(dd1 < -M_PI) dd1 += 2*M_PI;
                double w1 = exp(-0.5*dd1*dd1/(h*h));

                for(int di3 = -R; di3 <= R; di3++){
                    int i3 = ci3 + di3;
                    if(i3 < 0) i3 += G;
                    if(i3 >= G) i3 -= G;
                    double tc3 = Th_lo + (i3+0.5)*H_box;
                    double dd3 = th3 - tc3;
                    if(dd3 > M_PI) dd3 -= 2*M_PI;
                    if(dd3 < -M_PI) dd3 += 2*M_PI;

                    double w = w1 * exp(-0.5*dd3*dd3/(h*h)) / (2.0*M_PI*h*h);
                    mc_density[i1*G+i3] += w;
                }
            }
        }

        // Normalize: sum * H^2 = 1
        double mc_sum=0; for(int k=0;k<G*G;k++) mc_sum+=mc_density[k];
        for(int k=0;k<G*G;k++) mc_density[k]/=(mc_sum*H_box*H_box);

        // Step 2b: Predicted a_eff from BOTH formulas
        //   KDE smoothing flattens density, so a_eff < a_true
        //   T_eff = T + h^2 * <I1I2> gives a_eff = <I1I2>/T_eff
        //   Formula A (advisor's /2): T_eff_A = T + h^2 * <I1I2> / 2
        //   Formula B (H_code):       T_eff_B = T + h^2 * <I1I2>
        double T_eff_A = T_eq + h*h*mean_I1I2/2.0;
        double T_eff_B = T_eq + h*h*mean_I1I2;
        double a_pred_A = mean_I1I2 / T_eff_A;
        double a_pred_B = mean_I1I2 / T_eff_B;

        // Step 2c: Fine scan — must go BELOW a_true since smoothing reduces a
        double a_min_pred = min(a_pred_A, a_pred_B);
        double a_center = 0.5*(a_true + a_min_pred);
        double a_range = max(0.15, 2.0*fabs(a_true - a_min_pred) + 0.08);
        vector<double> a_scan(N_a);
        for(int i=0;i<N_a;i++) a_scan[i] = a_center - a_range/2.0 + a_range*i/(N_a-1);

        // Step 2d: Scan a_eff, compute RMSE for each
        double best_rmse = 1e30;
        double best_a = a_true;
        vector<double> theory_grid;

        for(int ia=0; ia<N_a; ia++){
            double a_eff = a_scan[ia];
            compute_theory(a_eff, theory_grid);
            double rmse = compute_rmse(mc_density, theory_grid);

            fcurve << h << "," << a_eff << "," << rmse << endl;

            if(rmse < best_rmse){
                best_rmse = rmse;
                best_a = a_eff;
            }
        }

        // Step 2e: Compute RMSE at a_true, a_pred_A, a_pred_B specifically
        compute_theory(a_true, theory_grid);
        double err_true = compute_rmse(mc_density, theory_grid);

        compute_theory(a_pred_A, theory_grid);
        double err_pred_A = compute_rmse(mc_density, theory_grid);

        compute_theory(a_pred_B, theory_grid);
        double err_pred_B = compute_rmse(mc_density, theory_grid);

        printf("%-6.2f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-10.6f  %-10.6f  %-10.6f  %-10.6f\n",
               h, a_true, a_pred_A, a_pred_B, best_a, err_true, err_pred_A, err_pred_B, best_rmse);

        fout << h << "," << a_true << "," << a_pred_A << "," << a_pred_B << "," << best_a
             << "," << err_true << "," << err_pred_A << "," << err_pred_B << "," << best_rmse << endl;
    }

    fout.close();
    fcurve.close();

    cout << endl;
    cout << "Output: teff_scan_results.csv (summary)" << endl;
    cout << "        teff_scan_curves.csv (full RMSE vs a curves)" << endl;

    gettimeofday(&t2, NULL);
    double sec=((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << endl << "Wall time = " << sec << "s" << endl;
    return 0;
}
