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
int N_initials = 64000;
// time
int T = 200;
int gap = 20;
int Burn_in = 0;
float dt = 0.01;
float sig = 1.5;
float drift = 1.0f;
//int gap = 10;

// num samples
long long int N_sample = 10000;
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


float observable(VectorXf X)
{
    return (2.0f - X.array().abs()).cwiseMax(0.0f).sum();
}


void MC(MatrixXf &Initial, MatrixXf &Obs_traj, const int N_thread)
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
        for (long long int n = 0; n < N_initials/16; n++)
        {
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> state_batch;
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> rnd;
            Eigen::Matrix<float, 6, 16, Eigen::RowMajor> drift;
            //inner loop

            for(int m = 0; m < N_sample; m++)
            {

                // initial conditions
                for(int i = 0; i < 16; i++)
                {
                    state_batch.col(i) = Initial.col(n*16 + i);
                }
                //cout<<"state_batch = "<<state_batch<<endl;

            // run the MC sims in a batch of 16
                for (int j = 0; j < T*gap + 1 ; j++)
                {
                // noise and vector field
                    generate_normal_batch_16(rng_state, rnd);
                    compute_batch_16(state_batch, drift);
                    state_batch = state_batch + dt*drift + sig*sqrt(dt)* rnd;
                    for(int i = 0; i < 6; i++)
                    {
                        state_batch.row(i) =  wrap_pi_16(state_batch.row(i));
                    }
                    if(j% gap == 0)
                    {
                        for(int i = 0; i < 16; i++)
                        {
                            Obs_traj(n*16+i, j/gap) += observable(state_batch.col(i));
                        }
                    }

                }
            }
            // keep track
//            cout << i << endl;
        }
    }

    


}

MatrixXf loadCSV(const std::string& path, int rows, int cols) {
    MatrixXf matrix(rows, cols);
    std::ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << path << std::endl;
        return matrix; // Returns empty/unitialized matrix
    }

    std::string line;
    int row = 0;
    while (std::getline(file, line) && row < rows) {
        std::stringstream lineStream(line);
        std::string cell;
        int col = 0;

        while (std::getline(lineStream, cell, ' ') && col < cols) {
            try {
                matrix(row, col) = std::stod(cell);
            } catch (const std::invalid_argument& e) {
                matrix(row, col) = 0.0; // Handle non-numeric junk
            }
            col++;
        }
        row++;
    }

    return matrix;
}

int main()
{
    // setup
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);
    random_device rd;
    uint32_t seed = rd();
    
    
    int N_thread = 10;
    MatrixXf myData = loadCSV("Kuramoto_backward_LHS_X_train.txt", N_initials, 6);
    
    MatrixXf myDataT = myData.transpose();
    //cout<<myDataT<<endl;
    cout<<"data loaded"<<endl;
    MatrixXf Obs_traj(N_initials, (T+1));
    Obs_traj.fill(0);
    ofstream myfile;
    
    
    gettimeofday(&t1,NULL);
    MC(myDataT, Obs_traj, N_thread);
    cout << "data generated" << endl;

    gettimeofday(&t2, NULL);
    myfile.open("Kuramoto_backward_Y_train2.txt");
    for (int i = 0; i < N_initials; i++)
    {
        for (int j = 0; j < T + 1; j++)
        {
                myfile << Obs_traj(i, j) << " ";
        }
        myfile << endl;
    }


    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
    cout<<"total CPU time = " <<delta<<endl;

    
    return 0;

}
