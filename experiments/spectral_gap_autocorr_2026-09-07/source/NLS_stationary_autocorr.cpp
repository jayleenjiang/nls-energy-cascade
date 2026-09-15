/*
 * Stationary-trajectory sampler for the n=3 NLS autocorrelation experiment.
 *
 * The Euler--Maruyama map is copied from cpp/backward/NLS_backward.cpp.
 * This file changes stream management and output only: every trajectory has
 * an explicit independent seed and writes its seven observables separately.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>

namespace fs = std::filesystem;

constexpr float PI = 3.14159265358979323846f;
constexpr float TWO_PI = 6.28318530717958647692f;
constexpr float ACTION_FLOOR = 1.0e-14f;
constexpr int NOBS = 7;

struct Config {
    std::string outdir;
    float T1 = 2.0f;
    float T3 = 8.0f;
    float gamma = 0.1f;
    float dt = 0.001f;
    double burn = 500.0;
    double measure = 1000.0;
    double snapshot = 0.01;
    int streams = 64;
    uint64_t base_seed = 2026090701ULL;
    int threads = 8;
};

struct State {
    float I1, I2, I3;
    float p1, p2, p3;
};

struct RNGState { uint32_t s[4]; };

static uint64_t splitmix64(uint64_t &z) {
    z += 0x9E3779B97F4A7C15ULL;
    uint64_t r = z;
    r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9ULL;
    r = (r ^ (r >> 27)) * 0x94D049BB133111EBULL;
    return r ^ (r >> 31);
}

static RNGState init_rng(uint64_t seed) {
    uint64_t z = seed;
    uint64_t p1 = splitmix64(z), p2 = splitmix64(z);
    RNGState st{{static_cast<uint32_t>(p1), static_cast<uint32_t>(p1 >> 32),
                 static_cast<uint32_t>(p2), static_cast<uint32_t>(p2 >> 32)}};
    if ((st.s[0] | st.s[1] | st.s[2] | st.s[3]) == 0U) st.s[0] = 1U;
    return st;
}

static uint32_t next_u32(RNGState &s) {
    uint32_t sum = s.s[0] + s.s[3];
    uint32_t r = (sum << 7) | (sum >> 25);
    uint32_t t = s.s[1] << 9;
    s.s[2] ^= s.s[0]; s.s[3] ^= s.s[1]; s.s[1] ^= s.s[2]; s.s[0] ^= s.s[3];
    s.s[2] ^= t;
    s.s[3] = (s.s[3] << 11) | (s.s[3] >> 21);
    return r;
}

static float uniform01(RNGState &s) {
    return static_cast<float>(next_u32(s) >> 9) * 0.00000011920928955078125f;
}

static float wrap_pi(float x) {
    return x - std::round(x / TWO_PI) * TWO_PI;
}

static float fast_sin(float x) {
    constexpr float B = 1.27323954f, C = -0.40528473f, P = 0.225f;
    float y = B * x + C * x * std::fabs(x);
    return P * (y * std::fabs(y) - y) + y;
}

static float fast_cos(float x) {
    return fast_sin(wrap_pi(x + 1.570796327f));
}

static std::array<float,4> normals4(RNGState &rng) {
    float u1 = uniform01(rng), u2 = uniform01(rng);
    float r1 = std::sqrt(-2.0f * std::log(std::max(1.0f - u1, 1.0e-30f)));
    float th1 = TWO_PI * u2;
    float u3 = uniform01(rng), u4 = uniform01(rng);
    float r2 = std::sqrt(-2.0f * std::log(std::max(1.0f - u3, 1.0e-30f)));
    float th2 = TWO_PI * u4;
    return {r1 * fast_cos(wrap_pi(th1)), r1 * fast_sin(wrap_pi(th1)),
            r2 * fast_cos(wrap_pi(th2)), r2 * fast_sin(wrap_pi(th2))};
}

static void step_em(State &x, RNGState &rng, const Config &c,
                    uint64_t &projection_count, float &min_proposed) {
    const float I1 = x.I1, I2 = x.I2, I3 = x.I3;
    const float th1 = wrap_pi(2.0f * (x.p1 - x.p2));
    const float th3 = wrap_pi(2.0f * (x.p3 - x.p2));
    const float s1 = fast_sin(th1), c1 = fast_cos(th1);
    const float s3 = fast_sin(th3), c3 = fast_cos(th3);
    const float M = I1 + I2 + I3;

    float dI1 = 4.0f * I1 * I2 * s1;
    float dI2 = 4.0f * I2 * (-I1 * s1 - I3 * s3);
    float dI3 = 4.0f * I3 * I2 * s3;
    float dp1 = 2.0f * M - I1 + 2.0f * I2 * c1;
    float dp2 = 2.0f * M - I2 + 2.0f * I1 * c1 + 2.0f * I3 * c3;
    float dp3 = 2.0f * M - I3 + 2.0f * I2 * c3;

    dI1 += 2.0f * c.gamma *
           (2.0f * c.T1 - (2.0f * M * I1 - I1 * I1 + 2.0f * I2 * I1 * c1));
    dp1 += c.gamma * (2.0f * I2 * s1);
    dI3 += 2.0f * c.gamma *
           (2.0f * c.T3 - (2.0f * M * I3 - I3 * I3 + 2.0f * I2 * I3 * c3));
    dp3 += c.gamma * (2.0f * I2 * s3);

    const auto z = normals4(rng);
    const float I1c = std::max(I1, ACTION_FLOOR);
    const float I3c = std::max(I3, ACTION_FLOOR);
    const float sI1 = 2.0f * std::sqrt(2.0f * c.gamma * c.T1 * I1c);
    const float sI3 = 2.0f * std::sqrt(2.0f * c.gamma * c.T3 * I3c);
    const float sp1 = std::sqrt(2.0f * c.gamma * c.T1 / I1c);
    const float sp3 = std::sqrt(2.0f * c.gamma * c.T3 / I3c);
    const float sqrt_dt = std::sqrt(c.dt);

    const float pI1 = I1 + dI1 * c.dt + sI1 * sqrt_dt * z[0];
    const float pI2 = I2 + dI2 * c.dt;
    const float pI3 = I3 + dI3 * c.dt + sI3 * sqrt_dt * z[1];
    min_proposed = std::min(min_proposed, std::min(pI1, std::min(pI2, pI3)));
    projection_count += static_cast<uint64_t>(pI1 < ACTION_FLOOR);
    projection_count += static_cast<uint64_t>(pI2 < ACTION_FLOOR);
    projection_count += static_cast<uint64_t>(pI3 < ACTION_FLOOR);

    x.I1 = std::max(pI1, ACTION_FLOOR);
    x.I2 = std::max(pI2, ACTION_FLOOR);
    x.I3 = std::max(pI3, ACTION_FLOOR);
    x.p1 = wrap_pi(x.p1 + dp1 * c.dt + sp1 * sqrt_dt * z[2]);
    x.p2 = wrap_pi(x.p2 + dp2 * c.dt);
    x.p3 = wrap_pi(x.p3 + dp3 * c.dt + sp3 * sqrt_dt * z[3]);
}

static std::array<float,NOBS> observe(const State &x) {
    const float th1 = wrap_pi(2.0f * (x.p1 - x.p2));
    const float th3 = wrap_pi(2.0f * (x.p3 - x.p2));
    return {std::cos(th1), std::sin(th1), std::cos(th3),
            std::cos(wrap_pi(th1 - th3)), x.I2, x.I1 + x.I3, x.I1 - x.I3};
}

static bool finite_state(const State &x) {
    return std::isfinite(x.I1) && std::isfinite(x.I2) && std::isfinite(x.I3) &&
           std::isfinite(x.p1) && std::isfinite(x.p2) && std::isfinite(x.p3);
}

int main(int argc, char **argv) {
    // Eleven user arguments plus argv[0].
    if (argc != 12) {
        std::cerr << "Usage: " << argv[0]
                  << " OUTDIR T1 T3 gamma dt burn measure snapshot streams base_seed threads\n";
        return 2;
    }
    Config c;
    c.outdir = argv[1]; c.T1 = std::stof(argv[2]); c.T3 = std::stof(argv[3]);
    c.gamma = std::stof(argv[4]); c.dt = std::stof(argv[5]);
    c.burn = std::stod(argv[6]); c.measure = std::stod(argv[7]);
    c.snapshot = std::stod(argv[8]); c.streams = std::stoi(argv[9]);
    c.base_seed = std::stoull(argv[10]); c.threads = std::stoi(argv[11]);

    const long long burn_steps = std::llround(c.burn / c.dt);
    const long long measure_steps = std::llround(c.measure / c.dt);
    const long long snapshot_steps = std::llround(c.snapshot / c.dt);
    if (snapshot_steps <= 0 || std::fabs(snapshot_steps * c.dt - c.snapshot) > 1e-7 ||
        measure_steps % snapshot_steps != 0) {
        std::cerr << "dt, measurement duration, and snapshot interval are incompatible\n";
        return 2;
    }
    const long long nsnap = measure_steps / snapshot_steps + 1;
    fs::create_directories(c.outdir);

    std::cout << std::setprecision(10)
              << "n=3 T1=" << c.T1 << " T3=" << c.T3 << " gamma=" << c.gamma
              << " dt=" << c.dt << " burn=" << c.burn << " measure=" << c.measure
              << " snapshot=" << c.snapshot << " streams=" << c.streams
              << " nsnap=" << nsnap << " base_seed=" << c.base_seed
              << " threads=" << c.threads << "\n";

    std::vector<uint64_t> projections(c.streams, 0);
    std::vector<float> min_proposed(c.streams, std::numeric_limits<float>::infinity());
    std::vector<int> nonfinite(c.streams, 0);
    std::atomic<int> completed{0};

#pragma omp parallel for schedule(dynamic) num_threads(c.threads)
    for (int stream = 0; stream < c.streams; ++stream) {
        RNGState rng = init_rng(c.base_seed + static_cast<uint64_t>(stream));
        State x{1.0f, 1.0f, 1.0f,
                (2.0f * uniform01(rng) - 1.0f) * PI,
                (2.0f * uniform01(rng) - 1.0f) * PI,
                (2.0f * uniform01(rng) - 1.0f) * PI};
        uint64_t proj = 0;
        float minp = std::numeric_limits<float>::infinity();
        bool bad = false;
        for (long long i = 0; i < burn_steps; ++i) {
            step_em(x, rng, c, proj, minp);
            if (!finite_state(x)) { bad = true; break; }
        }

        std::ostringstream name;
        name << c.outdir << "/stream_" << std::setw(3) << std::setfill('0') << stream << ".f32";
        const std::string final_path = name.str();
        const std::string tmp_path = final_path + ".partial";
        std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
        if (!out) bad = true;

        if (!bad) {
            auto values = observe(x);
            out.write(reinterpret_cast<const char*>(values.data()), sizeof(float) * NOBS);
            for (long long step = 1; step <= measure_steps; ++step) {
                step_em(x, rng, c, proj, minp);
                if (!finite_state(x)) { bad = true; break; }
                if (step % snapshot_steps == 0) {
                    values = observe(x);
                    out.write(reinterpret_cast<const char*>(values.data()), sizeof(float) * NOBS);
                }
            }
        }
        out.close();
        if (!bad) fs::rename(tmp_path, final_path);
        else {
            nonfinite[stream] = 1;
            std::error_code ec; fs::rename(tmp_path, final_path + ".failed", ec);
        }
        projections[stream] = proj;
        min_proposed[stream] = minp;
        int done = ++completed;
#pragma omp critical
        std::cout << "completed " << done << "/" << c.streams
                  << " stream=" << stream << " projections=" << proj
                  << " min_proposed=" << minp << " bad=" << bad << "\n";
    }

    std::ofstream meta(c.outdir + "/metadata.csv");
    meta << "stream_id,seed,projection_count,min_proposed_action,nonfinite\n";
    for (int i = 0; i < c.streams; ++i) {
        meta << i << ',' << (c.base_seed + static_cast<uint64_t>(i)) << ','
             << projections[i] << ',' << std::setprecision(17) << min_proposed[i]
             << ',' << nonfinite[i] << '\n';
    }
    return std::any_of(nonfinite.begin(), nonfinite.end(), [](int x){ return x != 0; }) ? 1 : 0;
}
