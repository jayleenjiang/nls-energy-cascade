
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
float T1_val    = 10.0f;
float T3_val    = 2.0f;
const float dt      = 0.001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

const float PI_f = 3.14159265358979323846f;
const int   G = 30;
const float H_box = 2.0f * PI_f / G;
const float sigma_kde = H_box / 3.0f;
const float Th_lo = -PI_f;

long long N_sample = 500000000LL;
int       N_thread = 8;

typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

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

// NON-EQUILIBRIUM EM step: T1 and T3 are different
inline void step_EM_batch16(State6x16& X, RNGState& rng) {
    A16f I1=X.row(0).array(),I2=X.row(1).array(),I3=X.row(2).array();
    A16f p1=X.row(3).array(),p2=X.row(4).array(),p3=X.row(5).array();
    A16f d12w=wrap_pi_16(2.0f*(p1-p2)),d32w=wrap_pi_16(2.0f*(p3-p2));
    A16f s12=fast_sin_16(d12w),c12=fast_cos_16(d12w);
    A16f s32=fast_sin_16(d32w),c32=fast_cos_16(d32w);
    A16f M=I1+I2+I3;
    A16f dI1=4.0f*I1*I2*s12, dI2=4.0f*I2*(-I1*s12-I3*s32), dI3=4.0f*I3*I2*s32;
    A16f dp1=2.0f*M-I1+2.0f*I2*c12, dp2=2.0f*M-I2+2.0f*I1*c12+2.0f*I3*c32, dp3=2.0f*M-I3+2.0f*I2*c32;

    // Mode 1: coupled to T1
    dI1+=2.0f*gamma_val*(2.0f*T1_val-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    // Mode 3: coupled to T3
    dI3+=2.0f*gamma_val*(2.0f*T3_val-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);

    A16f nI1,nI3,np1,np3; gen_noise_4x16(rng,nI1,nI3,np1,np3);
    A16f I1c=I1.max(1e-14f),I3c=I3.max(1e-14f);
    // Noise uses T1 for mode 1, T3 for mode 3
    A16f sI1=2.0f*(2.0f*gamma_val*T1_val*I1c).sqrt();
    A16f sI3=2.0f*(2.0f*gamma_val*T3_val*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T1_val/I1c).sqrt();
    A16f sp3=(2.0f*gamma_val*T3_val/I3c).sqrt();

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

inline float wrap_f(float x) { return x-roundf(x*0.159154943f)*6.283185307f; }

inline double normpdf(double x, double mu, double sig) {
    double d=x-mu; return exp(-0.5*d*d/(sig*sig))/(sig*sqrt(2.0*M_PI));
}
inline double normpdf2d(double x1, double x2, double mu1, double mu2, double sig) {
    return normpdf(x1,mu1,sig)*normpdf(x2,mu2,sig);
}

double bessel_I0(double x) {
    double s=1,t=1; for(int k=1;k<200;k++){t*=(x/(2.0*k))*(x/(2.0*k));s+=t;if(t<1e-16*s)break;} return s;
}

// Slice definition
struct SliceDef {
    float I_center;
    float I_tol;
};

void run_slice(const SliceDef& sl, const vector<vector<double>>& all_grids,
               const vector<double>& all_I1I2, const vector<double>& all_I2I3,
               const vector<long long>& all_counts, int N_thread) {
    // This function is not used — slices are processed inline below
}

int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T1_val    = atof(argv[3]);
    if(argc>4) T3_val    = atof(argv[4]);
    if(argc>5) N_sample  = atoll(argv[5]);

    cout << "=== Check 3: Non-equilibrium NC-KDE with T1*, T3* ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1_val << " T3=" << T3_val << endl;
    cout << "Grid: " << G << "x" << G << " H=" << H_box << " sigma_kde=" << sigma_kde << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << " (x16 SIMD)" << endl;
    cout << endl;

    // Define multiple slices
    const int N_slices = 3;
    float slice_centers[N_slices] = {1.0f, 2.0f, 3.0f};
    float slice_tol = 0.1f;

    // Per-slice accumulators
    vector<vector<double>> slice_grids(N_slices, vector<double>(G*G, 0.0));
    vector<double> slice_I1I2(N_slices, 0.0);
    vector<double> slice_I2I3(N_slices, 0.0);
    vector<long long> slice_counts(N_slices, 0);

    // Thread-local copies
    vector<vector<vector<double>>> t_grids(N_thread, vector<vector<double>>(N_slices, vector<double>(G*G, 0.0)));
    vector<vector<double>> t_I1I2(N_thread, vector<double>(N_slices, 0.0));
    vector<vector<double>> t_I2I3(N_thread, vector<double>(N_slices, 0.0));
    vector<vector<long long>> t_counts(N_thread, vector<long long>(N_slices, 0));

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

        for(long long step=0;step<N_sample;step++){
            step_EM_batch16(X,rng);
            for(int c=0;c<16;c++){
                float I1=X(0,c),I2=X(1,c),I3=X(2,c);

                // Check each slice
                for(int sl=0;sl<N_slices;sl++){
                    float Ic = slice_centers[sl];
                    if(fabsf(I1-Ic)<slice_tol && fabsf(I2-Ic)<slice_tol && fabsf(I3-Ic)<slice_tol){
                        float p1=X(3,c),p2=X(4,c),p3=X(5,c);
                        float th1=wrap_f(2.0f*(p1-p2)), th3=wrap_f(2.0f*(p3-p2));

                        int i1=(int)floorf((th1-Th_lo)/H_box);
                        int i3=(int)floorf((th3-Th_lo)/H_box);
                        if(i1<0)i1=0; if(i1>=G)i1=G-1;
                        if(i3<0)i3=0; if(i3>=G)i3=G-1;

                        double tc1=Th_lo+(i1+0.5)*H_box, tc3=Th_lo+(i3+0.5)*H_box;
                        double w=normpdf2d((double)th1,(double)th3,tc1,tc3,sigma_kde);
                        t_grids[rank][sl][i1*G+i3]+=w;
                        t_I1I2[rank][sl]+=(double)(I1*I2);
                        t_I2I3[rank][sl]+=(double)(I2*I3);
                        t_counts[rank][sl]++;
                    }
                }
            }
            if(rank==0 && step%100000000LL==0 && step>0)
                cout<<"  thread 0: "<<step/1000000<<"M/"<<N_sample/1000000<<"M"<<endl;
        }
    }

    // Merge threads
    for(int r=0;r<N_thread;r++){
        for(int sl=0;sl<N_slices;sl++){
            for(int k=0;k<G*G;k++) slice_grids[sl][k]+=t_grids[r][sl][k];
            slice_I1I2[sl]+=t_I1I2[r][sl];
            slice_I2I3[sl]+=t_I2I3[r][sl];
            slice_counts[sl]+=t_counts[r][sl];
        }
    }

    // Process each slice
    double h = (double)sigma_kde;
    for(int sl=0;sl<N_slices;sl++){
        float Ic = slice_centers[sl];
        if(slice_counts[sl] < 1000){
            cout << "=== Slice I=" << Ic << " : only " << slice_counts[sl] << " samples, skipping ===" << endl << endl;
            continue;
        }

        double mI1I2 = slice_I1I2[sl]/slice_counts[sl];
        double mI2I3 = slice_I2I3[sl]/slice_counts[sl];

        // T1*, T3*
        double T1s = T1_val - h*h*mI1I2/2.0;
        double T3s = T3_val - h*h*mI2I3/2.0;

        // a values: conditional density ∝ exp(-I1I2 cosθ1/T1 - I2I3 cosθ3/T3)
        // Using <I1I2> as the fixed I1I2 on the narrow slice
        double a1_orig = mI1I2 / T1_val;
        double a3_orig = mI2I3 / T3_val;
        double a1_star = mI1I2 / T1s;
        double a3_star = mI2I3 / T3s;

        // Normalize MC
        vector<double> mc(G*G);
        double mc_sum=0;
        for(int k=0;k<G*G;k++){ mc[k]=slice_grids[sl][k]; mc_sum+=mc[k]; }
        for(int k=0;k<G*G;k++) mc[k]/=(mc_sum*H_box*H_box);

        // Theory (T1, T3) — non-symmetric!
        // ρ ∝ exp(-a1 cosθ1 - a3 cosθ3)
        vector<double> th_orig(G*G), th_star(G*G);
        double Z_o=0, Z_s=0;
        for(int i=0;i<G;i++){
            double tc1=Th_lo+(i+0.5)*H_box;
            for(int j=0;j<G;j++){
                double tc3=Th_lo+(j+0.5)*H_box;
                th_orig[i*G+j]=exp(-a1_orig*cos(tc1)-a3_orig*cos(tc3));
                th_star[i*G+j]=exp(-a1_star*cos(tc1)-a3_star*cos(tc3));
                Z_o+=th_orig[i*G+j]; Z_s+=th_star[i*G+j];
            }
        }
        for(int k=0;k<G*G;k++){ th_orig[k]/=(Z_o*H_box*H_box); th_star[k]/=(Z_s*H_box*H_box); }

        // RMSE
        double rmse_o=0, rmse_s=0;
        for(int k=0;k<G*G;k++){
            double d1=mc[k]-th_orig[k], d2=mc[k]-th_star[k];
            rmse_o+=d1*d1; rmse_s+=d2*d2;
        }
        rmse_o=sqrt(rmse_o/(G*G)); rmse_s=sqrt(rmse_s/(G*G));

        cout << "=== Slice I=" << Ic << " +/- " << slice_tol << " ===" << endl;
        cout << "Samples: " << slice_counts[sl] << endl;
        cout << "<I1I2> = " << mI1I2 << "  <I2I3> = " << mI2I3 << endl;
        cout << "T1=" << T1_val << "  T1*=" << T1s << endl;
        cout << "T3=" << T3_val << "  T3*=" << T3s << endl;
        cout << "a1=" << a1_orig << "  a1*=" << a1_star << endl;
        cout << "a3=" << a3_orig << "  a3*=" << a3_star << endl;
        cout << "RMSE(T)=" << rmse_o << "  RMSE(T*)=" << rmse_s;
        if(rmse_o>0) cout << "  Improvement: " << (rmse_o-rmse_s)/rmse_o*100 << "%";
        cout << endl;

        // Corner/center
        int k_corner=0, k_center=(G/2)*G+G/2;
        printf("  Corner: MC=%.6f  T=%.6f  T*=%.6f  MC/T=%.4f  MC/T*=%.4f\n",
            mc[k_corner], th_orig[k_corner], th_star[k_corner],
            mc[k_corner]/th_orig[k_corner], mc[k_corner]/th_star[k_corner]);
        printf("  Center: MC=%.6f  T=%.6f  T*=%.6f  MC/T=%.4f  MC/T*=%.4f\n",
            mc[k_center], th_orig[k_center], th_star[k_center],
            mc[k_center]/th_orig[k_center], mc[k_center]/th_star[k_center]);
        cout << endl;
    }

    gettimeofday(&t2, NULL);
    double sec=((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << "Wall time = " << sec << "s" << endl;
    return 0;
}