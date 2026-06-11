
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

// Integrator
const float dt      = 0.0001f;
const float sqrt_dt = sqrtf(dt);
const int   Burn_in = 2000000;

// Domain
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

// NC-KDE bandwidths (one per dimension type)
// sigma_I for action variables, sigma_Th for angle variables
const float sigma_I  = hI / 3.0f;    // ~0.111
const float sigma_Th = hTh / 3.0f;   // ~0.014

// MC parameters
int       N_box    = 50000;
long long N_sample = 20000000LL;
int       N_thread = 8;
const double sample_ratio = 0.5;

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
inline void gen_noise_4x16(RNGState& s,A16f& n1,A16f& n2,A16f& n3,A16f& n4){
    const float twoPi=6.283185307f,piO2=1.570796327f;
    A16f u1=next_u01_16(s),u2=next_u01_16(s);
    A16f r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt(),th=twoPi*u2;
    n1=r*fast_sin_16(wrap_pi_16(th+piO2));n2=r*fast_sin_16(wrap_pi_16(th));
    u1=next_u01_16(s);u2=next_u01_16(s);
    r=(-2.0f*(1.0f-u1).max(1e-30f).log()).sqrt();th=twoPi*u2;
    n3=r*fast_sin_16(wrap_pi_16(th+piO2));n4=r*fast_sin_16(wrap_pi_16(th));
}

// ====================================================================
// SIMD EM step (corrected: dphi +=, noise with sqrt(2))
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
    // Dissipation (corrected sign: +=)
    dI1+=2.0f*gamma_val*(2.0f*T1-(2.0f*M*I1-I1*I1+2.0f*I2*I1*c12));
    dp1+=gamma_val*(2.0f*I2*s12);
    dI3+=2.0f*gamma_val*(2.0f*T3-(2.0f*M*I3-I3*I3+2.0f*I2*I3*c32));
    dp3+=gamma_val*(2.0f*I2*s32);
    // Noise (with sqrt(2))
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
// 5D coords and box lookup (unchanged from NLS5D_SIMD.cpp)
// ====================================================================
struct Coords5D { float I1, I2, I3, th1, th3; };

inline Coords5D to_5d_col(const State6x16& X, int col) {
    Coords5D c;
    c.I1=X(0,col); c.I2=X(1,col); c.I3=X(2,col);
    float p1=X(3,col),p2=X(4,col),p3=X(5,col);
    float th1=2.0f*(p1-p2), th3=2.0f*(p3-p2);
    th1=th1-roundf(th1*0.159154943f)*6.283185307f;
    th3=th3-roundf(th3*0.159154943f)*6.283185307f;
    c.th1=th1; c.th3=th3;
    return c;
}

int which_box(vector<vector<int>>& lob, const Coords5D& c) {
    if(c.I1<I_lo||c.I1>=I_max||c.I2<I_lo||c.I2>=I_max||c.I3<I_lo||c.I3>=I_max) return -1;
    int n1=(int)floorf((c.I1-I_lo)/hI), n2=(int)floorf((c.I2-I_lo)/hI);
    int n3=(int)floorf((c.I3-I_lo)/hI);
    int nt1=(int)floorf((c.th1-Th_lo)/hTh), nt3=(int)floorf((c.th3-Th_lo)/hTh);
    if(n1<0||n1>=N_I||n2<0||n2>=N_I||n3<0||n3>=N_I) return -1;
    nt1=max(0,min(nt1,N_Th-1)); nt3=max(0,min(nt3,N_Th-1));
    int idxA=n1*N_I+n2, idxB=sizeA+n3*N_Th+nt1, idxC=sizeA+sizeB+nt3;
    if(lob[idxA].empty()||lob[idxB].empty()||lob[idxC].empty()) return -1;
    static thread_local vector<int> ab,abc;
    ab.clear(); abc.clear();
    set_intersection(lob[idxA].begin(),lob[idxA].end(),lob[idxB].begin(),lob[idxB].end(),back_inserter(ab));
    if(ab.empty()) return -1;
    set_intersection(ab.begin(),ab.end(),lob[idxC].begin(),lob[idxC].end(),back_inserter(abc));
    return abc.empty()?-1:abc[0];
}

bool register_box(vector<vector<int>>& lob, const Coords5D& c, int box_id) {
    int n1=max(0,min((int)floorf((c.I1-I_lo)/hI),N_I-1));
    int n2=max(0,min((int)floorf((c.I2-I_lo)/hI),N_I-1));
    int n3=max(0,min((int)floorf((c.I3-I_lo)/hI),N_I-1));
    int nt1=max(0,min((int)floorf((c.th1-Th_lo)/hTh),N_Th-1));
    int nt3=max(0,min((int)floorf((c.th3-Th_lo)/hTh),N_Th-1));
    int idxA=n1*N_I+n2, idxB=sizeA+n3*N_Th+nt1, idxC=sizeA+sizeB+nt3;
    vector<int> ab;
    set_intersection(lob[idxA].begin(),lob[idxA].end(),lob[idxB].begin(),lob[idxB].end(),back_inserter(ab));
    if(!ab.empty()){vector<int> abc;
    set_intersection(ab.begin(),ab.end(),lob[idxC].begin(),lob[idxC].end(),back_inserter(abc));
    if(!abc.empty()) return false;}
    lob[idxA].push_back(box_id); sort(lob[idxA].begin(),lob[idxA].end());
    lob[idxB].push_back(box_id); sort(lob[idxB].begin(),lob[idxB].end());
    lob[idxC].push_back(box_id); sort(lob[idxC].begin(),lob[idxC].end());
    return true;
}

// ====================================================================
// Scalar EM (for create_boxes, double precision)
// ====================================================================
typedef Matrix<double, 6, 1> State6d;
inline double wrap_d(double x) { return x-2.0*M_PI*floor((x+M_PI)/(2.0*M_PI)); }

void step_EM_scalar(State6d& X, mt19937& rng, normal_distribution<double>& dist) {
    double I1=X(0),I2=X(1),I3=X(2),p1=X(3),p2=X(4),p3=X(5);
    double dpn0=2*(p1-p2),dpp2=2*(p3-p2);
    double M=I1+I2+I3;
    double s12=sin(dpn0),c12=cos(dpn0),s32=sin(dpp2),c32=cos(dpp2);
    double dI1=4*I1*I2*s12, dI2=4*I2*(-I1*s12-I3*s32), dI3=4*I3*I2*s32;
    double dp1=2*M-I1+2*I2*c12, dp2_=2*M-I2+2*I1*c12+2*I3*c32, dp3=2*M-I3+2*I2*c32;
    dI1+=2*gamma_val*(2*T1-(2*M*I1-I1*I1+2*I2*I1*c12));
    dp1+=gamma_val*(2*I2*s12);
    dI3+=2*gamma_val*(2*T3-(2*M*I3-I3*I3+2*I2*I3*c32));
    dp3+=gamma_val*(2*I2*s32);
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

// ====================================================================
// Create boxes (same as NLS5D_SIMD.cpp: ratio=0.5)
// ====================================================================
void create_boxes(vector<vector<int>>& lob, MatrixXf& Boxes) {
    State6d X; X<<1.0,1.0,0.1,0.0,0.0,0.0;
    mt19937 rng(42);
    normal_distribution<double> dist(0.0,1.0);
    uniform_real_distribution<float> uni(0.0f,1.0f);
    cout<<"  Burning in..."<<endl;
    for(int i=0;i<Burn_in;i++) step_EM_scalar(X,rng,dist);
    cout<<"  Sampling boxes (ratio="<<sample_ratio<<")..."<<endl;
    int count=0, from_traj=0, from_unif=0;
    while(count<N_box){
        Coords5D c;
        if(uni(rng)<sample_ratio){
            for(int i=0;i<5000;i++) step_EM_scalar(X,rng,dist);
            c.I1=(float)X(0);c.I2=(float)X(1);c.I3=(float)X(2);
            c.th1=(float)wrap_d(2*(X(3)-X(4))); c.th3=(float)wrap_d(2*(X(5)-X(4)));
            if(c.I1<I_lo||c.I1>=I_max||c.I2<I_lo||c.I2>=I_max||c.I3<I_lo||c.I3>=I_max) continue;
        } else {
            c.I1=I_lo+uni(rng)*(I_max-I_lo); c.I2=I_lo+uni(rng)*(I_max-I_lo);
            c.I3=I_lo+uni(rng)*(I_max-I_lo);
            c.th1=Th_lo+uni(rng)*(Th_hi-Th_lo); c.th3=Th_lo+uni(rng)*(Th_hi-Th_lo);
        }
        if(register_box(lob,c,count)){
            // Store EXACT box center coordinates (grid cell center)
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
            if(count%10000==0) cout<<"    "<<count<<"/"<<N_box
                <<" (traj="<<from_traj<<" unif="<<from_unif<<")"<<endl;
        }
    }
    cout<<"  Final: "<<from_traj<<" traj, "<<from_unif<<" unif"<<endl;
}

// ====================================================================
// 5D NC-KDE weight: product of 5 univariate Gaussian PDFs
// Uses separate bandwidths for I (sigma_I) and theta (sigma_Th)
// ====================================================================
inline double nckde_weight_5d(const Coords5D& sample, const MatrixXf& Boxes, int box_id) {
    double d1 = (double)sample.I1  - (double)Boxes(0, box_id);
    double d2 = (double)sample.I2  - (double)Boxes(1, box_id);
    double d3 = (double)sample.I3  - (double)Boxes(2, box_id);
    double d4 = (double)sample.th1 - (double)Boxes(3, box_id);
    double d5 = (double)sample.th3 - (double)Boxes(4, box_id);

    // Product kernel: normpdf in each dimension
    double exponent = -0.5 * (d1*d1/(sigma_I*sigma_I) + d2*d2/(sigma_I*sigma_I)
                             + d3*d3/(sigma_I*sigma_I)
                             + d4*d4/(sigma_Th*sigma_Th) + d5*d5/(sigma_Th*sigma_Th));
    // Normalization constant for 5D product kernel
    double norm = pow(sigma_I, 3) * pow(sigma_Th, 2) * pow(2.0*M_PI, 2.5);
    return exp(exponent) / norm;
}

// ====================================================================
// MC with NC-KDE
// ====================================================================
void MC_NCKDE(VectorXd& Box_count, vector<vector<int>>& lob, const MatrixXf& Boxes) {
    long long total_points = (long long)N_thread * 16LL * N_sample;
    cout << "MC NC-KDE: " << N_thread << " threads x 16 lanes, "
         << N_sample << " steps/traj" << endl;
    cout << "Total sample points = " << total_points << endl;
    cout << "NC-KDE bandwidths: sigma_I=" << sigma_I << " sigma_Th=" << sigma_Th << endl;

    vector<long long> counts(N_thread, 0);

#pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        RNGState rng = init_rng(rank*1000+42);
        mt19937 mt_init(rank*137+99);
        normal_distribution<float> nd(0.0f,1.0f);

        State6x16 X;
        for(int col=0;col<16;col++){
            X(0,col)=1.0f+0.1f*nd(mt_init); X(1,col)=1.0f+0.1f*nd(mt_init);
            X(2,col)=0.1f+0.05f*fabsf(nd(mt_init));
            X(3,col)=0.5f*nd(mt_init); X(4,col)=0; X(5,col)=0.5f*nd(mt_init);
        }
        for(int s=0;s<Burn_in;s++) step_EM_batch16(X,rng);

        long long local_hits=0;
        for(long long step=0; step<N_sample; step++){
            step_EM_batch16(X,rng);
            for(int col=0;col<16;col++){
                Coords5D c = to_5d_col(X,col);
                int idx = which_box(lob, c);
                if(idx >= 0 && idx < N_box) {
                    // NC-KDE: weight by normpdf instead of +1
                    double w = nckde_weight_5d(c, Boxes, idx);
                    Box_count(rank*N_box + idx) += w;
                    local_hits++;
                }
            }
            if(rank==0 && step%2000000LL==0 && step>0)
                cout<<"  thread 0: "<<step/1000000<<"M/"<<N_sample/1000000<<"M hits="<<local_hits<<endl;
        }
        counts[rank] = local_hits;
    }

    long long total=0;
    for(int i=0;i<N_thread;i++) total+=counts[i];
    cout << "Total hits = " << total << endl;

    // Merge
    for(int i=1;i<N_thread;i++)
        for(int j=0;j<N_box;j++)
            Box_count(j) += Box_count(i*N_box+j);
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
    if(argc>5) N_box     = atoi(argv[5]);
    if(argc>6) N_sample  = atoll(argv[6]);

    cout << "=== NLS5D SIMD NC-KDE (5D density estimation) ===" << endl;
    cout << "gamma=" << gamma_val << " T1=" << T1 << " T3=" << T3 << endl;
    cout << "dt=" << dt << " N_box=" << N_box << endl;
    cout << "N_sample/thread=" << N_sample << " N_thread=" << N_thread << endl;
    cout << "NC-KDE: sigma_I=" << sigma_I << " sigma_Th=" << sigma_Th << endl;
    cout << endl;

    // Phase 1: Create boxes
    vector<vector<int>> lob(total_list);
    for(auto& v:lob) v.reserve(8);
    MatrixXf Boxes(5, N_box);
    cout << "Creating boxes..." << endl;
    create_boxes(lob, Boxes);

    // Phase 2: MC with NC-KDE
    VectorXd Box_count(N_thread * N_box);
    Box_count.fill(0);
    MC_NCKDE(Box_count, lob, Boxes);

    // Phase 3: Density normalization
    // density = total_weight / (N_total_samples)
    // (NC-KDE weights already include the kernel normalization)
    long long total_samples = (long long)N_thread * 16LL * N_sample;
    VectorXd density = Box_count.head(N_box) / (double)total_samples;

    // Write output
    ofstream f1("NLS_FP_boxes.txt"), f2("NLS_FP_density.txt");
    for(int i=0;i<N_box;i++){
        f1<<Boxes(0,i)<<" "<<Boxes(1,i)<<" "<<Boxes(2,i)
          <<" "<<Boxes(3,i)<<" "<<Boxes(4,i)<<endl;
        f2<<density(i)<<endl;
    }
    f1.close(); f2.close();

    // Summary
    double mx=density.maxCoeff();
    int nz=0; for(int i=0;i<N_box;i++) if(density(i)>0) nz++;
    cout << endl;
    cout << "Max density = " << mx << endl;
    cout << "Nonzero = " << nz << "/" << N_box << endl;

    // T* info for reference
    // Note: NC-KDE (single-bin weighting) does not produce full convolution,
    // so T* correction is small. Record for completeness.
    cout << endl << "T* correction info (for reference):" << endl;
    cout << "  T* = T - h^2 * <I1I2>" << endl;
    cout << "  sigma_I=" << sigma_I << " sigma_Th=" << sigma_Th << endl;
    cout << "  (T* correction is small for single-bin NC-KDE)" << endl;

    gettimeofday(&t2, NULL);
    double sec=((t2.tv_sec-t1.tv_sec)*1e6+t2.tv_usec-t1.tv_usec)/1e6;
    cout << "Wall time = " << sec << "s" << endl;

    cout << endl << "Output: NLS_FP_boxes.txt, NLS_FP_density.txt" << endl;
    return 0;
}
