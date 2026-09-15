#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define NDEBUG
#define EIGEN_NO_DEBUG
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <complex>
#include <omp.h>
using namespace Eigen;
using namespace std;

// domain
const float M_PI_f = 3.14159265358979323846f;
float Sp = 2*M_PI_f;
int N = 40;
float h = Sp/N;
int N_box = 200000;
// time
int T = 100;
int gap = 20;
int Burn_in = 0;
float dt = 0.01;
float sig = 2.0;
float drift = 1.0f;
//int gap = 10;

// num samples
long long int N_sample = 1e6;
double ratio = 0.5;

// parameters

inline Eigen::Array<float, 16, 1> wrap_pi_16(const Eigen::Array<float, 16, 1>& x) {
    const float invTwoPi = 0.159154943f; // 1.0 / (2.0 * PI)
    const float twoPi = 6.283185307f;    // 2.0 * PI

    // 1. x / 2pi
    // 2. round to nearest integer (This is one instruction in AVX-512)
    // 3. x - (rounded * 2pi)
    return x - (x * invTwoPi).round() * twoPi;
}

inline Eigen::Array<float, 16, 1> fast_sin_16(const Eigen::Array<float, 16, 1>& x) {
    const float B = 1.27323954f; // 4/pi
    const float C = -0.40528473f; // -4/pi^2
    const float P = 0.225f;

    // x.abs() and x*x are mapped to vabsps and vmulps instructions
    auto y = B * x + C * x * x.abs();
    return P * (y * y.abs() - y) + y;
}

void compute_batch_16(const Eigen::Ref<const Eigen::Matrix<float, 6, 16, Eigen::RowMajor>>& x_block, Eigen::Ref<Eigen::Matrix<float, 6, 16, Eigen::RowMajor>> results)
{
    // Initialize all results to 1.0 (the constant in your formula)
    results.setConstant(drift);

    // We only need to compute 15 unique pairs (i, j) where i < j
    // because sin(x_i - x_j) = -sin(x_j - x_i)
    
    // Pairwise interaction indices for 6 variables
    static const int pairs_i[] = {0,0,0,0,0, 1,1,1,1, 2,2,2, 3,3, 4};
    static const int pairs_j[] = {1,2,3,4,5, 2,3,4,5, 3,4,5, 4,5, 5};

    for (int p = 0; p < 15; ++p) {
        int i = pairs_i[p];
        int j = pairs_j[p];
        
        //cout<<"row "<<i<<" = "<<x_block.row(i)<<endl;
        //cout<<"row "<<j<<" = "<<x_block.row(j)<<endl;
        
        // 1. Get x_i and x_j for all 16 trajectories (16 floats = 1 AVX-512 register)
        // 2. Compute the difference
        Eigen::Array<float, 16, 1> diff = x_block.row(i).array() - x_block.row(j).array();

        Eigen::Array<float, 16, 1> wrapped_diff = wrap_pi_16(diff);
        
        // 3. Compute sin(diff) for all 16 trajectories in one go
        Eigen::Array<float, 16, 1> s_val = fast_sin_16(wrapped_diff);

        // 4. Update y(i) and y(j) using symmetry
        results.row(i).array() -= s_val;
        results.row(j).array() += s_val;
    }
}


struct RNGState {
    uint32_t s[4]; // State for Xoshiro128++
};

inline Eigen::Array<float, 16, 1> next_u01_16(RNGState& state) {
    Eigen::Array<float, 16, 1> res;
    for (int i = 0; i < 16; ++i) {
        // Xoshiro128++ 1.0 logic
        const uint32_t result = ((state.s[0] + state.s[3]) << 7) | ((state.s[0] + state.s[3]) >> 25);
        const uint32_t t = state.s[1] << 9;
        state.s[2] ^= state.s[0];
        state.s[3] ^= state.s[1];
        state.s[1] ^= state.s[2];
        state.s[0] ^= state.s[3];
        state.s[2] ^= t;
        state.s[3] = (state.s[3] << 11) | (state.s[3] >> 21);
        
        // Convert to float in [0, 1)
        res(i) = (result >> 9) * 0.00000011920929f; // result / 2^23
    }
    return res;
}

typedef Eigen::Array<float, 16, 1> Array16f;

void generate_normal_batch_16(RNGState& state, Eigen::Ref<Eigen::Matrix<float, 6, 16, Eigen::RowMajor>> noise_out) {
    const float twoPi = 6.283185307f;
    const float piOverTwo = 1.570796327f;

    for (int r = 0; r < 6; r += 2) {
        // 1. Generate 16 uniforms (This is already your fast XOR logic)
        Array16f u1 = next_u01_16(state);
        Array16f u2 = next_u01_16(state);

        // 2. VECTORIZED LOG AND SQRT
        // On Intel, ensure you compile with -fveclib=SVML or use Intel Compiler.
        // Without SVML, Eigen uses a slower internal approximation or scalar fallback.
        Array16f radius = (-2.0f * (1.0f - u1).log()).sqrt();
        Array16f theta = twoPi * u2;

        // 3. YOUR FAST SIN/COS (Already vectorized)
        Array16f cos_theta = fast_sin_16(wrap_pi_16(theta + piOverTwo));
        Array16f sin_theta = fast_sin_16(wrap_pi_16(theta));

        // 4. Store (16-wide store is very fast)
        noise_out.row(r).array() = radius * cos_theta;
        noise_out.row(r+1).array() = radius * sin_theta;
    }
}

RNGState initialize_rng_for_rank(int rank) {
    RNGState state;
    
    // 1. Initial seed based on rank.
    // Adding a large constant (like a golden ratio) ensures rank 0 and 1
    // are far apart in the bit-space.
    uint64_t z = static_cast<uint64_t>(rank) + 0x9E3779B97F4A7C15ULL;
    
    auto splitmix64 = [&z]() -> uint64_t {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t result = z;
        result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9ULL;
        result = (result ^ (result >> 27)) * 0x94D049BB133111EBULL;
        return result ^ (result >> 31);
    };

    // 2. Fill the 128-bit Xoshiro state using two 64-bit SplitMix pulls
    uint64_t part1 = splitmix64();
    uint64_t part2 = splitmix64();

    state.s[0] = static_cast<uint32_t>(part1);
    state.s[1] = static_cast<uint32_t>(part1 >> 32);
    state.s[2] = static_cast<uint32_t>(part2);
    state.s[3] = static_cast<uint32_t>(part2 >> 32);

    // 3. Safety Check: Xoshiro states must not be all zeros.
    // SplitMix64 is guaranteed not to return 0 for these constants.
    
    return state;
}

float wrap(float angle)
{
    float wrapped_angle = fmod(angle + M_PI_f, 2*M_PI_f);
    if (wrapped_angle < 0)
    {
        wrapped_angle += 2*M_PI_f;
    }
    return wrapped_angle - M_PI_f;
}

VectorXf wrap(VectorXf x)
{
    VectorXf result = x;
    for(int i = 0; i < x.size(); i++)
    {
        result(i) = wrap(x(i));
    }
    return result;
}

VectorXf Kuramoto(const VectorXf xx)
{
    VectorXf tmp(6);
    for(int i = 0; i < 6; i++)
    {
        tmp(i) = drift;//drift
        for(int j = 0; j < 6; j++)
        {
            if(j != i)
            {
                tmp(i) += sin(xx(j) - xx(i));
            }
        }
    }
    return tmp;
}

int which_box(vector<vector<int>>& list_of_box, const Eigen::Ref<const Eigen::VectorXf>& xx)
{
    int ret_val = -1;
    //cout<<"xx = "<<xx.transpose()<<endl;
    int th1 = int(floor((xx(0) + M_PI_f)/h));
    int th2 = int(floor((xx(1) + M_PI_f)/h));
    int th3 = int(floor((xx(2) + M_PI_f)/h));
    int th4 = int(floor((xx(3) + M_PI_f)/h));
    int th5 = int(floor((xx(4) + M_PI_f)/h));
    int th6 = int(floor((xx(5) + M_PI_f)/h));
    

        //order: th1 to th6
    int th123 = th1*N*N + th2*N + th3;
    int th456 = th4*N*N + th5*N + th6;
    int Num_123 = N*N*N;
        //cout<<"nxy = "<<nxy<<" nzw = "<<nzw<<endl;
    if(th123 >= 0 && th123 < Num_123-1 && th456 >= 0 && th456 < Num_123 - 1)
    {
        if(list_of_box[th123].size()!=0 && list_of_box[Num_123 + th456].size()!=0)
        {
            vector<int> tmp;
            set_intersection(list_of_box[th123].begin(), list_of_box[th123].end(),list_of_box[Num_123 + th456].begin(), list_of_box[Num_123 + th456].end(), back_inserter(tmp));
            if(tmp.size() != 0)
            {
                ret_val = tmp[0];
            }
        }
    }
    
    return ret_val;
}


void create_boxes(vector<vector<int>>& list_of_box, MatrixXf& Boxes)
{
    VectorXf x_old(6);
    VectorXf x_new(6);
    x_old.fill(0);
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<float> u(0, 1.0);
    normal_distribution<float> normal(0.0, 1.0);
    VectorXd rnd(6);
    for(int i = 0; i < Burn_in; i++)
    {
        VectorXf rnd(6);
        for(int i = 0; i < 6; i++)
        {
            rnd(i) = sig*normal(mt);
        }
        x_new = x_old + dt*Kuramoto(x_old) + sqrt(dt)*rnd;
        x_new = wrap(x_new);
        x_old = x_new;
    }
    int th1,th2,th3,th4,th5,th6,th123,th456;
    int count = 0;
    int flag = 0;
    while(count < N_box)
    {
        flag = 0;
        if(u(mt) < ratio)
        {
            for(int i = 0; i < 500; i++)
            {
                VectorXf rnd(6);
                for(int i = 0; i < 6; i++)
                {
                    rnd(i) = sig*normal(mt);
                }
                x_new = x_old + dt*Kuramoto(x_old) + sqrt(dt)*rnd;
                x_new = wrap(x_new);
                x_old = x_new;
            }
            th1 = int(floor((x_old(0) + M_PI_f)/h));
            th2 = int(floor((x_old(1) + M_PI_f)/h));
            th3 = int(floor((x_old(2) + M_PI_f)/h));
            th4 = int(floor((x_old(3) + M_PI_f)/h));
            th5 = int(floor((x_old(4) + M_PI_f)/h));
            th6 = int(floor((x_old(5) + M_PI_f)/h));
            //order: th1 to th6
        }
        else
        {
            th1 = int(N*u(mt));
            th2 = int(N*u(mt));
            th3 = int(N*u(mt));
            th4 = int(N*u(mt));
            th5 = int(N*u(mt));
            th6 = int(N*u(mt));
        }
        int th123 = th1*N*N + th2*N + th3;
        int th456 = th4*N*N + th5*N + th6;
        int Num_123 = N*N*N;
        vector<int> tmp;
        if(th123 >= 0 && th123 < Num_123-1 && th456 >= 0 && th456 < Num_123 - 1)
        {
            set_intersection(list_of_box[th123].begin(), list_of_box[th123].end(),list_of_box[Num_123 + th456].begin(), list_of_box[Num_123 + th456].end(), back_inserter(tmp));
            if(tmp.size() == 0)
            {
                list_of_box[th123].push_back(count);
                sort(list_of_box[th123].begin(), list_of_box[th123].end());
                list_of_box[Num_123 + th456].push_back(count);
                sort(list_of_box[Num_123 + th456].begin(), list_of_box[Num_123 + th456].end());
                Boxes(0, count) = -M_PI_f + th1*h + h/2.0;
                Boxes(1, count) = -M_PI_f + th2*h + h/2.0;
                Boxes(2, count) = -M_PI_f + th3*h + h/2.0;
                Boxes(3, count) = -M_PI_f + th4*h + h/2.0;
                Boxes(4, count) = -M_PI_f + th5*h + h/2.0;
                Boxes(5, count) = -M_PI_f + th6*h + h/2.0;
                count++;
            }
        }
    }
    //cout<<Boxes;

    
}




void MC(vector<vector<int>> &list_of_box, MatrixXf &Box_count, const int N_thread)
{
    VectorXi count(N_thread);
    count.fill(0);
#pragma omp parallel num_threads(N_thread)
    {
        // setup
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        random_device rd;
        uint32_t seed = rank;
        mt19937 mt(rd() + rank);
        normal_distribution<double> norm(0.0, 1.0);
        RNGState rng_state = initialize_rng_for_rank(seed);
        
#pragma omp for
        for (long long int n = 0; n < N_sample/16; n++)
        {
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> state_batch;
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> rnd;
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> drift;
            // initial conditions

            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 16; j++)
                {
                    state_batch(i,j) = wrap(2.0 + norm(mt));
                }
            }
            
            VectorXf x_old = state_batch.col(0);
            VectorXf x_new(6);

            // run the MC sims in a batch of 16
            for (int j = 0; j < T*gap + 1 ; j++)
            {
                // noise and vector field
                generate_normal_batch_16(rng_state, rnd);
                compute_batch_16(state_batch, drift);
                state_batch = state_batch + dt*drift + sqrt(dt)* rnd;
                for(int i = 0; i < 6; i++)
                {
                    state_batch.row(i) =  wrap_pi_16(state_batch.row(i));
                }
                if(j% gap == 0)
                {
                    for(int i = 0; i < 16; i++)
                    {
                        int index = which_box(list_of_box, state_batch.col(i));
                        if(index>=0 && index < N_box)
                        {
                            //cout<<"T = "<<j/gap<<" index = "<<index<<endl;
                            Box_count(index, j/gap  + rank*(T + 1)) += 1.0;
                            count(rank) ++;
                        }
                    }
                }

            }

            // keep track
//            cout << i << endl;
        }
    }
    cout<<" total number of sample = "<<count.sum()<<endl;
    //ofstream myfile;
    //myfile.open("Kuramoto.txt");
    //myfile<<trajs;
    


}



int main()
{
    // setup
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    random_device rd;
    uint32_t seed = rd();
    
    
    int N_thread = 10;
    int Num = N*N*N;
    vector<vector<int>> list_of_box(2*Num);//indices of boxes in each projection
    for(int i = 0; i < 2*Num; i++)
    {
        list_of_box[i].reserve(10);
    }
    MatrixXf Boxes(6, N_box);//coordinate of boxes
    MatrixXf Box_count(N_box, (T+1)*N_thread);
    Box_count.fill(0);
    //cout<<data<<endl;
    create_boxes(list_of_box, Boxes);
    cout<<"box generated"<<endl;
    gettimeofday(&t2, NULL);
    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;

    cout<<"Box generation time = " <<delta<<endl;
    /*
    for(int i = 0; i < 2*Num; i++)
    {
        if(list_of_box[i].size() > 0)
        {
            cout<<" at index "<<i<<" ";
            for(int j = 0; j < list_of_box[i].size(); j++)
            {
                cout<<list_of_box[i][j]<<" ";
            }
            cout<<endl;
        }
    }
    */
    // generate data
    
    gettimeofday(&t1,NULL);
    MatrixXd data(N_thread, N * N);
    MC(list_of_box, Box_count, N_thread);
    cout << "data generated" << endl;
    for (int i = 0; i < T + 1; i++)
    {
        for (int j = 1; j < N_thread; j++)
        {
                Box_count.col(i) += Box_count.col(j*(T + 1) + i);
        }
    }
    gettimeofday(&t2, NULL);
    ofstream myfile, myfile2;
    myfile.open("Kuramoto_Y_train.txt");
    for (int i = 0; i < N_box; i++)
    {
        for (int j = 0; j < T; j++)
        {
                myfile << Box_count(i, j) << " ";
        }
        myfile << endl;
    }
    myfile2.open("Kuramoto_X_train.txt");
    myfile2<<Boxes.transpose()<<endl;
    delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
    cout<<"total CPU time = " <<delta<<endl;

    
    return 0;

}
