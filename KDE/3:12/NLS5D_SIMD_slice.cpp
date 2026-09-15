
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
float T1        = 5.0f;
float T3        = 5.0f;

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
// Histogram / KDE parameters
// ====================================================================
const float PI_f = 3.14159265358979323846f;
const float Th_lo = -PI_f, Th_hi = PI_f;
const int G_bin = 50;
const float d_bin = (Th_hi - Th_lo) / G_bin;
const int G_kde = 30;
const float d_kde = (Th_hi - Th_lo) / G_kde;

// ====================================================================
// MC parameters
// ====================================================================
long long N_sample = 200000000LL;
int       N_thread = 8;

// ====================================================================
// SIMD types
// ====================================================================
typedef Array<float, 16, 1> A16f;
typedef Matrix<float, 6, 16, RowMajor> State6x16;

// ====================================================================
// Fast math
// ====================================================================
inline A16f wrap_pi_16(const A16f& x) {
    const float invTwoPi = 0.159154943f;
    const float twoPi    = 6.283185307f;
    return x - (x * invTwoPi).round() * twoPi;
}

inline A16f fast_sin_16(const A16f& x) {
    const float B = 1.27323954f, C = -0.40528473f, P = 0.225f;
    auto y = B * x + C * x * x.abs();
    return P * (y * y.abs() - y) + y;
}

inline A16f fast_cos_16(const A16f& x) {
    return fast_sin_16(wrap_pi_16(x + 1.570796327f));
}

// ====================================================================
// RNG (Xoshiro128++)
// ====================================================================
struct RNGState { uint32_t s[4]; };

RNGState init_rng(int rank) {
    RNGState st;
    uint64_t z = (uint64_t)rank + 0x9E3779B97F4A7C15ULL;
    auto sm = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t r = z;
        r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9ULL;
        r = (r ^ (r >> 27)) * 0x94D049BB133111EBULL;
        return r ^ (r >> 31);
    };
    uint64_t p1 = sm(), p2 = sm();
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
// SIMD EM step (identical physics to NLS5D_SIMD.cpp)
// ====================================================================
inline void step_EM_batch16(State6x16& X, RNGState& rng) {
    A16f I1=X.row(0).array(), I2=X.row(1).array(), I3=X.row(2).array();
    A16f p1=X.row(3).array(), p2=X.row(4).array(), p3=X.row(5).array();

    A16f d12w=wrap_pi_16(2.0f*(p1-p2)), d32w=wrap_pi_16(2.0f*(p3-p2));
    A16f s12=fast_sin_16(d12w), c12=fast_cos_16(d12w);
    A16f s32=fast_sin_16(d32w), c32=fast_cos_16(d32w);
    A16f M=I1+I2+I3;

    // Hamiltonian
    A16f dI1=4.0f*I1*I2*s12;
    A16f dI2=4.0f*I2*(-I1*s12-I3*s32);
    A16f dI3=4.0f*I3*I2*s32;
    A16f dp1=2.0f*M-I1+2.0f*I2*c12;
    A16f dp2=2.0f*M-I2+2.0f*I1*c12+2.0f*I3*c32;
    A16f dp3=2.0f*M-I3+2.0f*I2*c32;

    // Dissipation (mode 1)
    dI1+=2.0f*gamma_val*(2.0f*T1-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    // Dissipation (mode 3)
    dI3+=2.0f*gamma_val*(2.0f*T3-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);

    // Noise
    A16f nI1,nI3,np1,np3;
    gen_noise_4x16(rng, nI1,nI3,np1,np3);

    A16f I1c=I1.max(1e-14f), I3c=I3.max(1e-14f);
    A16f sI1=2.0f*(2.0f*gamma_val*T1*I1c).sqrt();
    A16f sI3=2.0f*(2.0f*gamma_val*T3*I3c).sqrt();
    A16f sp1=(2.0f*gamma_val*T1/I1c).sqrt();
    A16f sp3=(2.0f*gamma_val*T3/I3c).sqrt();

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
// Check if column is in slice
// ====================================================================
inline bool in_slice(const State6x16& X, int col) {
    return fabsf(X(0,col)-I_slice)<I_tol
        && fabsf(X(1,col)-I_slice)<I_tol
        && fabsf(X(2,col)-I_slice)<I_tol;
}

inline float wrap_f(float x) {
    return x - roundf(x*0.159154943f)*6.283185307f;
}

// ====================================================================
// Bessel I_0
// ====================================================================
double bessel_I0(double x) {
    double s=1,t=1;
    for(int k=1;k<200;k++){t*=(x/(2.0*k))*(x/(2.0*k));s+=t;if(t<1e-16*s)break;}
    return s;
}

// ====================================================================
// Periodic Gaussian for KDE
// ====================================================================
inline double normpdf_p(double x, double mu, double h) {
    double d=x-mu; if(d>M_PI)d-=2*M_PI; if(d<-M_PI)d+=2*M_PI;
    double z=d/h; return exp(-0.5*z*z)/(h*sqrt(2*M_PI));
}

// ====================================================================
// Main
// ====================================================================
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    if(argc>1) N_thread  = atoi(argv[1]);
    if(argc>2) gamma_val = atof(argv[2]);
    if(argc>3) T1        = atof(argv[3]);
    if(argc>4) T3        = atof(argv[4]);
    if(argc>5) N_sample  = atoll(argv[5]);
    if(argc>6) I_slice   = atof(argv[6]);
    if(argc>7) I_tol     = atof(argv[7]);

    bool eq = (fabsf(T1-T3)<1e-6f);
    float a = I_slice*I_slice/T1;  // theoretical parameter

    cout<<"=== NLS SIMD Slice Test ==="<<endl;
    cout<<"gamma="<<gamma_val<<" T1="<<T1<<" T3="<<T3;
    if(eq) cout<<" [EQUILIBRIUM]";
    cout<<endl;
    cout<<"I_slice="<<I_slice<<" +/- "<<I_tol<<" dt="<<dt<<endl;
    cout<<"N_sample/thread="<<N_sample<<" N_thread="<<N_thread<<" (x16 SIMD)"<<endl;
    if(eq) cout<<"Theory: a = I^2/T = "<<a<<endl;
    cout<<endl;

    // === Collect samples ===
    vector<float> all_th1, all_th3;
    long long total_steps=0, total_sl=0;

#pragma omp parallel num_threads(N_thread)
    {
        int rank=omp_get_thread_num();
        RNGState rng=init_rng(rank*1000+42);
        mt19937 mt_init(rank*137+99);
        normal_distribution<float> nd(0.0f,1.0f);

        State6x16 X;
        for(int c=0;c<16;c++){
            X(0,c)=1.0f+0.1f*nd(mt_init);
            X(1,c)=1.0f+0.1f*nd(mt_init);
            X(2,c)=0.1f+0.05f*fabsf(nd(mt_init));
            X(3,c)=0.5f*nd(mt_init); X(4,c)=0; X(5,c)=0.5f*nd(mt_init);
        }

        for(int s=0;s<Burn_in;s++) step_EM_batch16(X,rng);

        vector<float> loc_th1, loc_th3;
        long long loc_steps=0, loc_sl=0;

        for(long long step=0; step<N_sample; step++){
            step_EM_batch16(X,rng);
            loc_steps+=16;

            for(int c=0;c<16;c++){
                if(in_slice(X,c)){
                    float p1=X(3,c), p2=X(4,c), p3=X(5,c);
                    loc_th1.push_back(wrap_f(2.0f*(p1-p2)));
                    loc_th3.push_back(wrap_f(2.0f*(p3-p2)));
                    loc_sl++;
                }
            }

            if(rank==0 && step%50000000LL==0 && step>0)
                cout<<"  thread 0: "<<step/1000000<<"M/"<<N_sample/1000000<<"M"<<endl;
        }

#pragma omp critical
        {
            all_th1.insert(all_th1.end(), loc_th1.begin(), loc_th1.end());
            all_th3.insert(all_th3.end(), loc_th3.begin(), loc_th3.end());
            total_steps+=loc_steps; total_sl+=loc_sl;
        }
    }

    long long N_sl=all_th1.size();
    cout<<"Total steps: "<<total_steps<<endl;
    cout<<"In slice: "<<N_sl<<" ("<<100.0*N_sl/max(total_steps,1LL)<<"%)"<<endl<<endl;

    if(N_sl<100){cout<<"Too few samples.\n";return 1;}

    // === <cos(theta)> ===
    double mc1=0,mc3=0;
    for(long long s=0;s<N_sl;s++){mc1+=cos(all_th1[s]);mc3+=cos(all_th3[s]);}
    mc1/=N_sl; mc3/=N_sl;
    cout<<"Measured <cos(th1)> = "<<mc1<<endl;
    cout<<"Measured <cos(th3)> = "<<mc3<<endl;
    if(eq){
        double I0a=bessel_I0(a), I1a=0;
        for(int k=0;k<200;k++){
            double t=1; for(int j=1;j<=k;j++) t*=(a/2.0)/j;
            for(int j=1;j<=k+1;j++) t*=(a/2.0)/j; I1a+=t;
        }
        cout<<"Theory  <cos(th)>  = "<<(-I1a/I0a)<<" (a="<<a<<")"<<endl;
    }
    cout<<endl;

    // === Binning ===
    vector<double> hist(G_bin*G_bin, 0.0);
    for(long long s=0;s<N_sl;s++){
        int j1=max(0,min((int)floor((all_th1[s]-Th_lo)/d_bin),G_bin-1));
        int j3=max(0,min((int)floor((all_th3[s]-Th_lo)/d_bin),G_bin-1));
        hist[j1*G_bin+j3]+=1.0;
    }
    vector<double> dens_bin(G_bin*G_bin);
    for(int i=0;i<G_bin*G_bin;i++) dens_bin[i]=hist[i]/(N_sl*d_bin*d_bin);

    // === KDE ===
    double h_vals[2]={0.05,0.10};
    vector<vector<double>> kde_res(2);
    for(int ih=0;ih<2;ih++){
        double h=h_vals[ih];
        cout<<"KDE h="<<h<<"..."<<endl;
        kde_res[ih].resize(G_kde*G_kde,0);
#pragma omp parallel for num_threads(N_thread) schedule(dynamic,4)
        for(int j1=0;j1<G_kde;j1++){
            double tc=Th_lo+(j1+0.5)*d_kde;
            for(int j3=0;j3<G_kde;j3++){
                double t3c=Th_lo+(j3+0.5)*d_kde;
                double sum=0;
                for(long long s=0;s<N_sl;s++)
                    sum+=normpdf_p(tc,all_th1[s],h)*normpdf_p(t3c,all_th3[s],h);
                kde_res[ih][j1*G_kde+j3]=sum/(double)N_sl;
            }
        }
    }

    // === Theory comparison ===
    if(eq){
        cout<<endl<<"=== Theory comparison ==="<<endl;
        double I0a=bessel_I0(a);
        double Z=(2*M_PI*I0a)*(2*M_PI*I0a);
        cout<<"a="<<a<<" I0(a)="<<I0a<<" Z="<<Z<<endl;

        // Theory on KDE grid
        vector<double> theory(G_kde*G_kde);
        for(int j1=0;j1<G_kde;j1++){
            double t1=Th_lo+(j1+0.5)*d_kde;
            for(int j3=0;j3<G_kde;j3++){
                double t3=Th_lo+(j3+0.5)*d_kde;
                theory[j1*G_kde+j3]=exp(-a*(cos(t1)+cos(t3)))/Z;
            }
        }

        // Error metrics
        for(int ih=0;ih<2;ih++){
            double l2=0,linf=0;
            for(int i=0;i<G_kde*G_kde;i++){
                double d=fabs(kde_res[ih][i]-theory[i]);
                l2+=d*d; if(d>linf)linf=d;
            }
            l2=sqrt(l2/(G_kde*G_kde));
            printf("  KDE h=%.2f: RMSE=%.6f  MaxErr=%.6f\n",h_vals[ih],l2,linf);
        }

        // Theory on bin grid
        vector<double> theory_bin(G_bin*G_bin);
        for(int j1=0;j1<G_bin;j1++){
            double t1=Th_lo+(j1+0.5)*d_bin;
            for(int j3=0;j3<G_bin;j3++){
                double t3=Th_lo+(j3+0.5)*d_bin;
                theory_bin[j1*G_bin+j3]=exp(-a*(cos(t1)+cos(t3)))/Z;
            }
        }
        {
            double l2=0,linf=0;
            for(int i=0;i<G_bin*G_bin;i++){
                double d=fabs(dens_bin[i]-theory_bin[i]);
                l2+=d*d; if(d>linf)linf=d;
            }
            l2=sqrt(l2/(G_bin*G_bin));
            printf("  Binning %dx%d: RMSE=%.6f  MaxErr=%.6f\n",G_bin,G_bin,l2,linf);
        }

        // Cross-section
        cout<<endl<<"  Cross-section theta3~0:"<<endl;
        int j3m=G_kde/2;
        printf("  %8s %10s %10s %10s\n","theta1","Theory","KDE_0.05","KDE_0.10");
        for(int j1=0;j1<G_kde;j1++){
            double t1=Th_lo+(j1+0.5)*d_kde;
            printf("  %8.3f %10.6f %10.6f %10.6f\n",
                t1,theory[j1*G_kde+j3m],kde_res[0][j1*G_kde+j3m],kde_res[1][j1*G_kde+j3m]);
        }
    }

    // === Write CSV ===
    {
        ofstream f("density_bin.csv");
        f<<"theta1,theta3,density\n";
        for(int j1=0;j1<G_bin;j1++){
            float t1=Th_lo+(j1+0.5f)*d_bin;
            for(int j3=0;j3<G_bin;j3++){
                float t3=Th_lo+(j3+0.5f)*d_bin;
                f<<t1<<","<<t3<<","<<dens_bin[j1*G_bin+j3]<<"\n";
            }
        }
    }
    for(int ih=0;ih<2;ih++){
        char fn[64]; snprintf(fn,sizeof(fn),"density_kde_h%.2f.csv",h_vals[ih]);
        ofstream f(fn);
        f<<"theta1,theta3,density\n";
        for(int j1=0;j1<G_kde;j1++){
            float t1=Th_lo+(j1+0.5f)*d_kde;
            for(int j3=0;j3<G_kde;j3++){
                float t3=Th_lo+(j3+0.5f)*d_kde;
                f<<t1<<","<<t3<<","<<kde_res[ih][j1*G_kde+j3]<<"\n";
            }
        }
    }
    if(eq){
        ofstream f("density_theory.csv");
        f<<"theta1,theta3,density\n";
        double I0a=bessel_I0(a); double Z=(2*M_PI*I0a)*(2*M_PI*I0a);
        for(int j1=0;j1<G_kde;j1++){
            double t1=Th_lo+(j1+0.5)*d_kde;
            for(int j3=0;j3<G_kde;j3++){
                double t3=Th_lo+(j3+0.5)*d_kde;
                f<<t1<<","<<t3<<","<<exp(-a*(cos(t1)+cos(t3)))/Z<<"\n";
            }
        }
    }

    gettimeofday(&t2,NULL);
    double sec=((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout<<endl<<"Wall time = "<<sec<<"s"<<endl;
    return 0;
}
