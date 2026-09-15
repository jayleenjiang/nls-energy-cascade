#!/usr/bin/env bash
set -euo pipefail

MODULES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${MODULES_DIR}/../.." && pwd)"

link_one() {
  local target="$1"
  local link_path="$2"
  mkdir -p "$(dirname "${link_path}")"
  if [[ -L "${link_path}" ]]; then
    if [[ "$(readlink "${link_path}")" == "${target}" ]]; then
      return
    fi
    echo "Refusing to replace existing symlink: ${link_path}" >&2
    exit 2
  fi
  if [[ -e "${link_path}" ]]; then
    echo "Refusing to replace existing path: ${link_path}" >&2
    exit 2
  fi
  ln -s "${target}" "${link_path}"
}

# 00 manuscript versions
link_one "${REPO_ROOT}/Paper/revision/draft.tex" "${MODULES_DIR}/00_manuscript_versions/local/sources/current_local_enhanced.tex"
link_one "${REPO_ROOT}/Paper/backups/2026-08-25_repo_organization/draft_advisor_synced_2026-08-25.tex" "${MODULES_DIR}/00_manuscript_versions/local/sources/advisor_synced_2026-08-25.tex"
link_one "${REPO_ROOT}/Paper/backups/2026-08-25_repo_organization/main_enhanced_reference_2026-08-25.tex" "${MODULES_DIR}/00_manuscript_versions/local/sources/enhanced_reference_2026-08-25.tex"
link_one "${REPO_ROOT}/Paper/revision/references.bib" "${MODULES_DIR}/00_manuscript_versions/local/sources/references.bib"
link_one "${REPO_ROOT}/Paper/packages/NLS_best_current_2026-07-21" "${MODULES_DIR}/00_manuscript_versions/local/build/best_current_package"

# 01 model and theory
link_one "${REPO_ROOT}/Energy Cascade/Rect Plots/main.cpp" "${MODULES_DIR}/01_model_and_theory/local/code/cascade_geometry.cpp"
link_one "${REPO_ROOT}/Energy Cascade/Rect Plots/plot_rectangles.py" "${MODULES_DIR}/01_model_and_theory/local/code/plot_rectangles.py"
link_one "${REPO_ROOT}/Energy Cascade/Rect Plots/verify_rectangles.py" "${MODULES_DIR}/01_model_and_theory/local/code/verify_rectangles.py"
link_one "${REPO_ROOT}/Energy Cascade/Rect Plots" "${MODULES_DIR}/01_model_and_theory/local/data/cascade_geometry_workdir"
link_one "${REPO_ROOT}/Paper/revision/cascade_embedding.pdf" "${MODULES_DIR}/01_model_and_theory/local/figures/cascade_embedding.pdf"
link_one "${REPO_ROOT}/Paper/revision/cascade_embedding.png" "${MODULES_DIR}/01_model_and_theory/local/figures/cascade_embedding.png"
link_one "${REPO_ROOT}/Paper/CKSTT.pdf" "${MODULES_DIR}/01_model_and_theory/local/references/CKSTT.pdf"
link_one "${REPO_ROOT}/Paper/Summary_NLS.pdf" "${MODULES_DIR}/01_model_and_theory/local/references/Summary_NLS.pdf"
link_one "${REPO_ROOT}/Paper/Non-equilibrium steady state for a three-mode energy cascade model.pdf" "${MODULES_DIR}/01_model_and_theory/local/references/three_mode_NESS.pdf"

# 02 long-chain NESS and LTE
link_one "${REPO_ROOT}/cpp/simulation/lte_histogram_simd.cpp" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/lte_histogram_simd.cpp"
link_one "${REPO_ROOT}/cpp/simulation/lte_histogram.cpp" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/lte_histogram.cpp"
link_one "${REPO_ROOT}/cpp/simulation/lte_simd_conditioned.cpp" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/lte_simd_conditioned.cpp"
link_one "${REPO_ROOT}/cpp/simulation/lte_simd_dt25.cpp" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/lte_simd_dt25.cpp"
link_one "${REPO_ROOT}/python/analysis/hist_analysis.py" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/hist_analysis.py"
link_one "${REPO_ROOT}/python/analysis/lte_full_analysis.py" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/lte_full_analysis.py"
link_one "${REPO_ROOT}/Paper/revision/scripts/generate_manuscript_figures.py" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/generate_manuscript_figures.py"
link_one "${REPO_ROOT}/Paper/revision/scripts/export_compare_residual_mesh_metrics.py" "${MODULES_DIR}/02_long_chain_ness_lte/local/code/export_compare_residual_mesh_metrics.py"
link_one "${REPO_ROOT}/lte/n15 data" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/n15"
link_one "${REPO_ROOT}/lte/n25 data" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/n25"
link_one "${REPO_ROOT}/lte/n50 data" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/n50"
link_one "${REPO_ROOT}/lte/n100 data" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/n100"
link_one "${REPO_ROOT}/lte/T712:8:9 data" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/temperature_controls"
link_one "${REPO_ROOT}/lte/cond_hist" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/conditioned_histograms"
link_one "${REPO_ROOT}/experiments/lte" "${MODULES_DIR}/02_long_chain_ness_lte/local/data/curated_profiles_and_analysis"
link_one "${REPO_ROOT}/Paper/revision/action_profiles.pdf" "${MODULES_DIR}/02_long_chain_ness_lte/local/figures/action_profiles.pdf"
link_one "${REPO_ROOT}/Paper/revision/lte_residual_midchain.pdf" "${MODULES_DIR}/02_long_chain_ness_lte/local/figures/lte_residual_midchain.pdf"
link_one "${REPO_ROOT}/Paper/revision/report_assets/compare_residual_mesh.pdf" "${MODULES_DIR}/02_long_chain_ness_lte/local/figures/compare_residual_mesh.pdf"
link_one "${REPO_ROOT}/Paper/revision/report_assets/compare_residual_matlab.pdf" "${MODULES_DIR}/02_long_chain_ness_lte/local/figures/compare_residual_matlab.pdf"
link_one "${REPO_ROOT}/Paper/revision/source_trace_metrics.json" "${MODULES_DIR}/02_long_chain_ness_lte/local/reports/source_trace_metrics.json"
link_one "${REPO_ROOT}/Paper/revision/manuscript_figure_metrics.json" "${MODULES_DIR}/02_long_chain_ness_lte/local/reports/manuscript_figure_metrics.json"
link_one "${REPO_ROOT}/Paper/revision/report_assets/compare_residual_mesh_metrics.json" "${MODULES_DIR}/02_long_chain_ness_lte/local/reports/compare_residual_mesh_metrics.json"

# 03 action-current scaling
link_one "${REPO_ROOT}/flux/NLS_flux_canonical.cpp" "${MODULES_DIR}/03_action_current_scaling/local/code/canonical/NLS_flux_canonical.cpp"
link_one "${REPO_ROOT}/flux/analyze_canonical_flux.py" "${MODULES_DIR}/03_action_current_scaling/local/code/canonical/analyze_canonical_flux.py"
link_one "${REPO_ROOT}/flux/analyze_current_windows.py" "${MODULES_DIR}/03_action_current_scaling/local/code/canonical/analyze_current_windows.py"
link_one "${REPO_ROOT}/flux/flux_data/NLS_flux_SIMD.cpp" "${MODULES_DIR}/03_action_current_scaling/local/code/simd_comparison/NLS_flux_SIMD.cpp"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/production_dt5e-4" "${MODULES_DIR}/03_action_current_scaling/local/data/canonical_production"
link_one "${REPO_ROOT}/flux_data/non_simd_flux_data" "${MODULES_DIR}/03_action_current_scaling/local/data/canonical_distribution_export"
link_one "${REPO_ROOT}/flux/flux_data/fixed_simd_2026-06-19" "${MODULES_DIR}/03_action_current_scaling/local/data/fixed_simd_comparison"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.pdf" "${MODULES_DIR}/03_action_current_scaling/local/figures/flux_primary_scaling.pdf"
link_one "${REPO_ROOT}/Paper/revision/report_assets/simd_canonical_scaling_overlay.pdf" "${MODULES_DIR}/03_action_current_scaling/local/figures/simd_canonical_scaling_overlay.pdf"
link_one "${REPO_ROOT}/Paper/revision/report_assets/canonical_flux_distribution.pdf" "${MODULES_DIR}/03_action_current_scaling/local/figures/canonical_flux_distribution.pdf"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/production_manifest.md" "${MODULES_DIR}/03_action_current_scaling/local/reports/production_manifest.md"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/validation_report.md" "${MODULES_DIR}/03_action_current_scaling/local/reports/validation_report.md"

# 04 burn-in and finite-time distributions
link_one "${REPO_ROOT}/flux/NLS_flux_relaxation_tau.cpp" "${MODULES_DIR}/04_burnin_flux_distribution/local/code/NLS_flux_relaxation_tau.cpp"
link_one "${REPO_ROOT}/flux/analyze_burnin_ld.py" "${MODULES_DIR}/04_burnin_flux_distribution/local/code/analyze_burnin_ld.py"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/burnin_ld_2026-07-10" "${MODULES_DIR}/04_burnin_flux_distribution/local/data/pilot_1024"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/burnin_ld_full100k_2026-08-25" "${MODULES_DIR}/04_burnin_flux_distribution/local/data/full100k_active"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/burnin_ld_2026-07-10/analysis" "${MODULES_DIR}/04_burnin_flux_distribution/local/figures/pilot_analysis"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/burnin_ld_2026-07-10/README.md" "${MODULES_DIR}/04_burnin_flux_distribution/local/reports/pilot_readme.md"

# 05 short-chain Fokker--Planck density
link_one "${REPO_ROOT}/cpp/fp5d" "${MODULES_DIR}/05_short_chain_fp_density/local/code/fp5d_cpp"
link_one "${REPO_ROOT}/NN notebooks/FKE_5d_NLS.ipynb" "${MODULES_DIR}/05_short_chain_fp_density/local/code/FKE_5d_NLS.ipynb"
link_one "${REPO_ROOT}/Paper/revision/scripts/recompute_short_chain_nn_metrics.py" "${MODULES_DIR}/05_short_chain_fp_density/local/code/recompute_short_chain_nn_metrics.py"
link_one "${REPO_ROOT}/KDE/4:15_NN" "${MODULES_DIR}/05_short_chain_fp_density/local/data/fp5d_archived_run"
link_one "${REPO_ROOT}/Paper/revision/eq_validation.png" "${MODULES_DIR}/05_short_chain_fp_density/local/figures/eq_validation.png"
link_one "${REPO_ROOT}/Paper/revision/neq_density.png" "${MODULES_DIR}/05_short_chain_fp_density/local/figures/neq_density.png"
link_one "${REPO_ROOT}/Paper/revision/short_chain_nn_rerun_metrics.json" "${MODULES_DIR}/05_short_chain_fp_density/local/reports/short_chain_nn_rerun_metrics.json"
link_one "${REPO_ROOT}/Paper/revision/source_trace_metrics.json" "${MODULES_DIR}/05_short_chain_fp_density/local/reports/source_trace_metrics.json"

# 06 short-chain stabilization
link_one "${REPO_ROOT}/NN notebooks/FKE_4d.ipynb" "${MODULES_DIR}/06_short_chain_stabilization/local/code/FKE_4d.ipynb"
link_one "${REPO_ROOT}/NN notebooks/FKE_5d_NLS.ipynb" "${MODULES_DIR}/06_short_chain_stabilization/local/code/FKE_5d_NLS.ipynb"
link_one "${REPO_ROOT}/python/analysis/find_peak.py" "${MODULES_DIR}/06_short_chain_stabilization/local/code/find_peak.py"
link_one "${REPO_ROOT}/KDE/4:15_NN/NLS_FP_density.txt" "${MODULES_DIR}/06_short_chain_stabilization/local/data/NLS_FP_density.txt"
link_one "${REPO_ROOT}/KDE/4:15_NN/h5_files/final.keras" "${MODULES_DIR}/06_short_chain_stabilization/local/data/final_density_model.keras"
link_one "${REPO_ROOT}/Paper/revision/symmetry_breaking.png" "${MODULES_DIR}/06_short_chain_stabilization/local/figures/symmetry_breaking.png"
link_one "${REPO_ROOT}/Paper/revision/short_chain_nn_rerun_metrics.json" "${MODULES_DIR}/06_short_chain_stabilization/local/reports/short_chain_nn_rerun_metrics.json"

# 07 generator spectrum
link_one "${REPO_ROOT}/cpp/backward" "${MODULES_DIR}/07_generator_spectrum/local/code/backward_cpp"
link_one "${REPO_ROOT}/NN notebooks/FKE_eigen.ipynb" "${MODULES_DIR}/07_generator_spectrum/local/code/FKE_eigen.ipynb"
link_one "${REPO_ROOT}/python/data_gen/generate_LHS_backward.py" "${MODULES_DIR}/07_generator_spectrum/local/code/generate_LHS_backward.py"
link_one "${REPO_ROOT}/python/data_gen/get_train.py" "${MODULES_DIR}/07_generator_spectrum/local/code/get_train.py"
link_one "${REPO_ROOT}/Paper/revision/scripts/analyze_eigen_fit_windows.py" "${MODULES_DIR}/07_generator_spectrum/local/code/analyze_eigen_fit_windows.py"
link_one "${REPO_ROOT}/KDE/backward_NLS_X.txt" "${MODULES_DIR}/07_generator_spectrum/local/data/backward_NLS_X.txt"
link_one "${REPO_ROOT}/KDE/backward_NLS_Q1.txt" "${MODULES_DIR}/07_generator_spectrum/local/data/backward_NLS_Q1.txt"
link_one "${REPO_ROOT}/KDE/4:23_eigen/NLS_backward_Y_train.txt" "${MODULES_DIR}/07_generator_spectrum/local/data/NLS_backward_Y_train.txt"
link_one "${REPO_ROOT}/KDE/h5_files_eigen/final.keras" "${MODULES_DIR}/07_generator_spectrum/local/data/final_eigen_model.keras"
link_one "${REPO_ROOT}/Paper/revision/eigenvalue_scatter.png" "${MODULES_DIR}/07_generator_spectrum/local/figures/eigenvalue_scatter.png"
link_one "${REPO_ROOT}/Paper/revision/Q1_slices.png" "${MODULES_DIR}/07_generator_spectrum/local/figures/Q1_slices.png"
link_one "${REPO_ROOT}/KDE/eigen_3_histograms.png" "${MODULES_DIR}/07_generator_spectrum/local/figures/eigen_3_histograms.png"
link_one "${REPO_ROOT}/KDE/eigen_4_example_fits.png" "${MODULES_DIR}/07_generator_spectrum/local/figures/eigen_4_example_fits.png"
link_one "${REPO_ROOT}/Paper/revision/eigen_fit_sensitivity.json" "${MODULES_DIR}/07_generator_spectrum/local/reports/eigen_fit_sensitivity.json"
link_one "${REPO_ROOT}/Paper/revision/short_chain_nn_rerun_metrics.json" "${MODULES_DIR}/07_generator_spectrum/local/reports/short_chain_nn_rerun_metrics.json"

# 08 validation and robustness
link_one "${REPO_ROOT}/Paper/revision/scripts/analyze_flux_scaling_sensitivity.py" "${MODULES_DIR}/08_validation_robustness/local/code/analyze_flux_scaling_sensitivity.py"
link_one "${REPO_ROOT}/Paper/revision/scripts/run_gamma_robustness_production.py" "${MODULES_DIR}/08_validation_robustness/local/code/run_gamma_robustness_production.py"
link_one "${REPO_ROOT}/flux/compare_gibbs_sde.py" "${MODULES_DIR}/08_validation_robustness/local/code/compare_gibbs_sde.py"
link_one "${REPO_ROOT}/flux/gibbs_mcmc_reference.py" "${MODULES_DIR}/08_validation_robustness/local/code/gibbs_mcmc_reference.py"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/pilot" "${MODULES_DIR}/08_validation_robustness/local/data/timestep_and_stationarity_pilot"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/larger_n_pilot_2026-06-20" "${MODULES_DIR}/08_validation_robustness/local/data/larger_n50"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/larger_n60_pilot_2026-06-20" "${MODULES_DIR}/08_validation_robustness/local/data/larger_n60"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/parameter_robustness_2026-06-20" "${MODULES_DIR}/08_validation_robustness/local/data/bath_temperature_robustness"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/gamma_robustness_2026-06-21" "${MODULES_DIR}/08_validation_robustness/local/data/gamma_robustness"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/gibbs_mcmc" "${MODULES_DIR}/08_validation_robustness/local/data/gibbs_mcmc"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/gibbs_sde" "${MODULES_DIR}/08_validation_robustness/local/data/gibbs_sde"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/validation_report.md" "${MODULES_DIR}/08_validation_robustness/local/reports/validation_report.md"
link_one "${REPO_ROOT}/Paper/revision/experiments/flux_validation/production_manifest.md" "${MODULES_DIR}/08_validation_robustness/local/reports/production_manifest.md"
link_one "${REPO_ROOT}/Paper/revision/pre_submission_reviewer_audit_2026-06-21.md" "${MODULES_DIR}/08_validation_robustness/local/reports/pre_submission_reviewer_audit.md"

# 09 reproducibility and submission
link_one "${REPO_ROOT}/Paper/revision/scripts" "${MODULES_DIR}/09_submission_reproducibility/local/code/audit_and_packaging_scripts"
link_one "${REPO_ROOT}/Paper/revision/source_trace_metrics.json" "${MODULES_DIR}/09_submission_reproducibility/local/data/source_trace_metrics.json"
link_one "${REPO_ROOT}/Paper/revision/raw_data_archive_manifest.md" "${MODULES_DIR}/09_submission_reproducibility/local/data/raw_data_archive_manifest.md"
link_one "${REPO_ROOT}/Paper/revision/manuscript_claim_audit.md" "${MODULES_DIR}/09_submission_reproducibility/local/reports/manuscript_claim_audit.md"
link_one "${REPO_ROOT}/Paper/revision/reference_integrity_audit.md" "${MODULES_DIR}/09_submission_reproducibility/local/reports/reference_integrity_audit.md"
link_one "${REPO_ROOT}/Paper/revision/submission_checks_summary.md" "${MODULES_DIR}/09_submission_reproducibility/local/reports/submission_checks_summary.md"
link_one "${REPO_ROOT}/Paper/packages/NLS_best_current_2026-07-21" "${MODULES_DIR}/09_submission_reproducibility/local/build/best_current_package"

echo "Paper module navigation links are ready under ${MODULES_DIR}."
