#!/usr/bin/env bash
# Untrained-moment check against the long trajectories, for the v2 baseline models (../unified_v2/density)
# and the v5 models (models/), plus the n=4 conditional sin-theta / bounded weak-residual diagnostic.
# Submit with --dependency=afterok:<train_v5 job> after the v2 density training has finished.
#SBATCH --job-name=v5eval
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm_eval_%j.log
set -euo pipefail
export TF_CPP_MIN_LOG_LEVEL=3
source /home/jiang37j_mtholyoke_edu/NLS/nls7d_stream_2026-10-03/venv/bin/activate
L=$(cd ../unified_v2/density && pwd)
mkdir -p eval
# v5 models laid out like a results directory (models + code + Gibbs references + old short validation files)
mkdir -p eval/v5_res
for f in n3_common.py n3_model.py n4_common.py n4_model.py gibbs_n3_T5.bin gibbs_n4_T5.bin \
         n3_eq_val1.bin n3_eq_val2.bin n3_neq_val1.bin n3_neq_val2.bin n4_eq_val1.bin n4_eq_val2.bin n4_neq_val1.bin n4_neq_val2.bin; do
  ln -sfn "$L/$f" eval/v5_res/$f
done
for m in n3_eq n3_neq_s1 n3_neq_s2 n4_eq n4_neq_s1 n4_neq_s2; do ln -sfn "$(pwd)/models/$m" eval/v5_res/$m; done
cd eval
python3 ../eval_check_long.py "$L" ../val_long check_long_v2_baseline.json > check_long_v2_baseline.log 2>&1
python3 ../eval_check_long.py v5_res ../val_long check_long_v5.json > check_long_v5.log 2>&1
RES=$L VAL=../val_long python3 ../diag_sin_n4.py "$L/n4_neq_s1" "$L/n4_neq_s2" "$(pwd)/v5_res/n4_neq_s1" "$(pwd)/v5_res/n4_neq_s2" > diag_sin_n4.log 2>&1
grep -h "z_rms" check_long_v2_baseline.log check_long_v5.log
grep -h "weak z" diag_sin_n4.log
