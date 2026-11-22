#!/usr/bin/env bash
# Density training v5 on the v2 (co-rotating Cartesian) labels: the same configuration as the v2 density
# training plus --weak-set bounded (bounded angular stationarity test functions) and the long validation
# trajectories (--val-layout streams).  Submit with --dependency=afterok:<val_long job>.
#SBATCH --job-name=v5train
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_train_%A_%a.log
set -euo pipefail
export TF_CPP_MIN_LOG_LEVEL=3
source /home/jiang37j_mtholyoke_edu/NLS/nls7d_stream_2026-10-03/venv/bin/activate
NCPU=${SLURM_CPUS_PER_TASK:-16}
L=../unified_v2/density
for f in n3_eq_kde_kde.bin n3_neq_kde_kde.bin n4_eq_kde_kde.bin n4_neq_kde_kde.bin gibbs_n3_T5.bin gibbs_n4_T5.bin \
         n3_queries_eq.bin n3_queries_neq.bin n4_queries_eq.bin n4_queries_neq.bin n3_common.py n3_model.py n4_common.py n4_model.py train_kde_nn.py; do
  [ -e "$L/$f" ] || { echo "missing $L/$f (v2 density prep not finished?)"; exit 1; }
done
for f in n3_common.py n3_model.py n4_common.py n4_model.py train_kde_nn.py; do [ -e "$f" ] || cp "$L/$f" .; done
V5="--weak-set bounded --val-layout streams --ref-eval 200000 --pde-weight 0.3 --weak-weight 30 --moment-weight 30 --steps 50"
N3="--labels $L/n3_neq_kde_kde.bin --queries $L/n3_queries_neq.bin --reference $L/gibbs_n3_T5.bin --validation val_long/n3_neq_valA.bin val_long/n3_neq_valB.bin $V5"
N4="--labels $L/n4_neq_kde_kde.bin --queries $L/n4_queries_neq.bin --reference $L/gibbs_n4_T5.bin --validation val_long/n4_neq_valA.bin val_long/n4_neq_valB.bin $V5"
JOBS=(
  "n3_eq|train_kde_nn_n3_v5.py|--case eq --labels $L/n3_eq_kde_kde.bin --queries $L/n3_queries_eq.bin --reference $L/gibbs_n3_T5.bin --validation val_long/n3_eq_valA.bin val_long/n3_eq_valB.bin $V5 --epochs 100"
  "n3_neq_s1|train_kde_nn_n3_v5.py|--case neq $N3 --epochs 600 --seed 1"
  "n3_neq_s2|train_kde_nn_n3_v5.py|--case neq $N3 --epochs 600 --seed 2"
  "n4_eq|train_kde_nn_v5.py|--case eq --labels $L/n4_eq_kde_kde.bin --queries $L/n4_queries_eq.bin --reference $L/gibbs_n4_T5.bin --validation val_long/n4_eq_valA.bin val_long/n4_eq_valB.bin $V5 --epochs 100"
  "n4_neq_s1|train_kde_nn_v5.py|--case neq $N4 --epochs 600 --seed 1"
  "n4_neq_s2|train_kde_nn_v5.py|--case neq $N4 --epochs 600 --seed 2"
)
J=${JOBS[$SLURM_ARRAY_TASK_ID]}
OUT=models/${J%%|*}; R=${J#*|}; SCRIPT=${R%%|*}; ARGS=${R#*|}
mkdir -p models
echo "task $SLURM_ARRAY_TASK_ID -> $OUT ($SCRIPT) on $(hostname), cpus=$NCPU, $(date)"
# shellcheck disable=SC2086
python3 $SCRIPT $ARGS --output "$OUT" --threads "$NCPU" --jit 1 > "$OUT.log" 2>&1
echo "done $OUT $(date)"
python3 -c "import json;r=json.load(open('$OUT/result.json'));f=r['final'];print('runtime_h=%.2f'%(r['runtime_seconds']/3600),'heldout=%.4f'%f['heldout_label_rmse'],'fp=%.3f'%f['fp_abs_median'],'weak=%.4f'%f['weak_rms'],'mom=%.4f'%f['moment_rms'],'check_z=%.2f'%f['check_z_rms'],'I=',[round(v,3) for v in f['model_mean_I']],'J=',[round(v,3) for v in f['model_currents']])"
