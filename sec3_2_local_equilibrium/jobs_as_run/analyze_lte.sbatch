#!/usr/bin/env bash
# Frozen Sec. 3.2 analysis (same estimators, bootstrap seed and gates as fixed_reference_2026-09-18)
# applied to raw/cartesian_v2.  Submit with --dependency=afterok:<run_lte job>.
#SBATCH --job-name=lte_v2_an
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_lte_analysis_%j.log
set -euo pipefail
source /home/jiang37j_mtholyoke_edu/NLS/nls7d_stream_2026-10-03/venv/bin/activate
python3 analyze_lte_v2.py
ls -l analysis
