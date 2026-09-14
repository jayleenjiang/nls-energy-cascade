// Lossless reduced-state writer around the frozen controlled n=3 integrator.
// The included source supplies the validated Cartesian dynamics verbatim.
#define main nls_original_f32_main
#include "../../spectral_gap_controlled_2026-09-07/source/NLS_stationary_autocorr_controlled.cpp"
#undef main

namespace {

constexpr int NSTATE5 = 5;

std::array<double, NSTATE5> observe_state5(const State& state, int lane) {
    const double phi1 = std::atan2(state.y[0](lane), state.x[0](lane));
    const double phi2 = std::atan2(state.y[1](lane), state.x[1](lane));
    const double phi3 = std::atan2(state.y[2](lane), state.x[2](lane));
    return {
        state.action[0](lane),
        state.action[1](lane),
        state.action[2](lane),
        wrap_angle(2.0 * (phi1 - phi2)),
        wrap_angle(2.0 * (phi3 - phi2)),
    };
}

BatchResult run_batch_state5(const Config& c,
                             int batch,
                             std::int64_t burn_steps,
                             std::int64_t measure_steps,
                             std::int64_t snapshot_steps,
                             std::int64_t nsnap) {
    BatchResult result;
    result.min_action.fill(std::numeric_limits<double>::infinity());
    State state;
    MidpointWorkspace workspace;
    initialize_state(state, c.T1, c.T3);
    std::array<Xoshiro256pp, LANES> rng;
    for (int lane = 0; lane < LANES; ++lane) {
        const auto id = static_cast<std::uint64_t>(batch * LANES + lane);
        rng[lane].seed(stream_seed(c.base_seed, id));
    }

    std::vector<double> buffer(
        static_cast<std::size_t>(nsnap) * LANES * NSTATE5);
    auto record_snapshot = [&](std::int64_t snapshot_index) {
        compute_force(state);
        for (int lane = 0; lane < LANES; ++lane) {
            const auto values = observe_state5(state, lane);
            const std::size_t offset =
                (static_cast<std::size_t>(lane) * nsnap + snapshot_index) * NSTATE5;
            std::copy(values.begin(), values.end(), buffer.begin() + offset);
            for (int j = 0; j < N; ++j) {
                const double action = state.action[j](lane);
                result.min_action[lane] = std::min(result.min_action[lane], action);
                if (action == 0.0) ++result.zero_radius_events[lane];
            }
        }
    };

    A16d normal_left_x, normal_left_y, normal_right_x, normal_right_y;
    const std::int64_t total_steps = burn_steps + measure_steps;
    std::int64_t snapshot_index = 0;
    bool bad_batch = false;
    for (std::int64_t step = 0; step < total_steps; ++step) {
        if (step == burn_steps) record_snapshot(snapshot_index++);
        fill_gaussian_pair(rng, normal_left_x, normal_right_x);
        fill_gaussian_pair(rng, normal_left_y, normal_right_y);
        int iterations = 0;
        const bool converged =
            hamiltonian_midpoint_step(state, workspace, c.dt, iterations);
        result.midpoint_iteration_sum += static_cast<std::uint64_t>(iterations);
        result.midpoint_iteration_max = std::max(result.midpoint_iteration_max, iterations);
        if (!converged) {
            for (auto& value : result.midpoint_failures) ++value;
        }
        if (step % 2 == 0) {
            bath_step(state, 0, c.T1, c.gamma, normal_left_x, normal_left_y, c.dt);
            bath_step(state, N - 1, c.T3, c.gamma, normal_right_x, normal_right_y, c.dt);
        } else {
            bath_step(state, N - 1, c.T3, c.gamma, normal_right_x, normal_right_y, c.dt);
            bath_step(state, 0, c.T1, c.gamma, normal_left_x, normal_left_y, c.dt);
        }
        if (!finite_state(state)) {
            bad_batch = true;
            for (auto& value : result.nonfinite) value = 1;
            break;
        }
        if (step >= burn_steps &&
            (step - burn_steps + 1) % snapshot_steps == 0) {
            record_snapshot(snapshot_index++);
        }
    }

    if (!bad_batch && snapshot_index != nsnap) {
        throw std::runtime_error("internal snapshot count mismatch");
    }
    fs::create_directories(c.outdir);
    for (int lane = 0; lane < LANES; ++lane) {
        const int stream = batch * LANES + lane;
        std::ostringstream filename;
        filename << "stream_" << std::setw(3) << std::setfill('0') << stream << ".f64";
        const fs::path final_path = c.outdir / filename.str();
        const fs::path partial_path = final_path.string() + ".partial";
        std::ofstream output(partial_path, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("cannot create " + partial_path.string());
        if (!bad_batch) {
            const std::size_t offset =
                static_cast<std::size_t>(lane) * nsnap * NSTATE5;
            output.write(reinterpret_cast<const char*>(buffer.data() + offset),
                         static_cast<std::streamsize>(
                             nsnap * NSTATE5 * sizeof(double)));
        }
        output.close();
        if (!bad_batch) fs::rename(partial_path, final_path);
        else fs::rename(partial_path, final_path.string() + ".failed");
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Config c = parse_args(argc, argv);
        const auto burn_steps = exact_steps(c.burn, c.dt, "burn");
        const auto measure_steps = exact_steps(c.measure, c.dt, "measure");
        const auto snapshot_steps = exact_steps(c.snapshot, c.dt, "snapshot");
        if (snapshot_steps <= 0 || measure_steps % snapshot_steps != 0) {
            throw std::invalid_argument("measurement/snapshot step mismatch");
        }
        const auto nsnap = measure_steps / snapshot_steps + 1;
        const int streams = c.batches * LANES;
        fs::create_directories(c.outdir);

        std::cout << std::setprecision(17)
                  << "model=nls-n3-cartesian-controlled-state5-f64-v1"
                  << " parent_model=" << MODEL_VERSION
                  << " n=3 coordinates=cartesian-double"
                  << " projection=none floor=none adaptive_step=no standard_trig=yes"
                  << " T1=" << c.T1 << " T3=" << c.T3
                  << " gamma=" << c.gamma << " dt=" << c.dt
                  << " burn=" << c.burn << " measure=" << c.measure
                  << " snapshot=" << c.snapshot << " batches=" << c.batches
                  << " streams=" << streams << " nsnap=" << nsnap
                  << " base_seed=" << c.base_seed << " threads=" << c.threads
                  << '\n';

        std::vector<BatchResult> results(c.batches);
        std::atomic<int> completed{0};
#pragma omp parallel for schedule(static) num_threads(c.threads)
        for (int batch = 0; batch < c.batches; ++batch) {
            results[batch] = run_batch_state5(
                c, batch, burn_steps, measure_steps, snapshot_steps, nsnap);
            const int done = ++completed;
#pragma omp critical
            std::cout << "completed " << done << "/" << c.batches
                      << " batches (batch=" << batch << ")\n";
        }

        std::ofstream meta(c.outdir / "metadata.csv");
        meta << "stream_id,stream_seed,projection_count,floor_count,"
                "zero_radius_events,min_sampled_action,midpoint_failure_count,"
                "midpoint_iteration_max,nonfinite\n";
        std::uint64_t total_failures = 0;
        std::uint64_t total_zero = 0;
        int any_nonfinite = 0;
        double global_min_action = std::numeric_limits<double>::infinity();
        for (int batch = 0; batch < c.batches; ++batch) {
            for (int lane = 0; lane < LANES; ++lane) {
                const int stream = batch * LANES + lane;
                const auto seed = stream_seed(c.base_seed, stream);
                const auto& r = results[batch];
                meta << stream << ',' << seed << ",0,0,"
                     << r.zero_radius_events[lane] << ','
                     << std::setprecision(17) << r.min_action[lane] << ','
                     << r.midpoint_failures[lane] << ','
                     << r.midpoint_iteration_max << ',' << r.nonfinite[lane] << '\n';
                total_failures += r.midpoint_failures[lane];
                total_zero += r.zero_radius_events[lane];
                any_nonfinite |= r.nonfinite[lane];
                global_min_action = std::min(global_min_action, r.min_action[lane]);
            }
        }

        std::ofstream format(c.outdir / "FORMAT.txt");
        format << "model=nls-n3-cartesian-controlled-state5-f64-v1\n"
               << "parent_model=" << MODEL_VERSION << '\n'
               << "state_precision=float64\noutput_precision=float64\n"
               << "projection=none\nfloor=none\nadaptive_step=no\n"
               << "columns=I1,I2,I3,theta1,theta3\n"
               << "shape_per_stream=" << nsnap << "x" << NSTATE5 << '\n'
               << "T1=" << c.T1 << "\nT3=" << c.T3 << "\ngamma=" << c.gamma
               << "\ndt=" << c.dt << "\nburn=" << c.burn
               << "\nmeasure=" << c.measure << "\nsnapshot=" << c.snapshot
               << "\nstreams=" << streams << "\nbase_seed=" << c.base_seed << '\n';

        std::cout << "projection_count=0 floor_count=0 zero_radius_events="
                  << total_zero << " min_sampled_action=" << global_min_action
                  << " midpoint_failure_stream_events=" << total_failures
                  << " any_nonfinite=" << any_nonfinite << '\n';
        return any_nonfinite || total_failures || total_zero ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}

