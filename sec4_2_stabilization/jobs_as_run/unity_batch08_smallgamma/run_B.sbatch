#!/usr/bin/env bash
#SBATCH --job-name=b08_B
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-10
#SBATCH --output=slurm_runB_%A_%a.log
# One line of tasks_B.tsv per array task: <folder> <name> <arguments>. Output: out/<folder>/<name>_{summary.txt,dump.f32,log}.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks_B.tsv)
IFS=$'\t' read -r dir name args <<<"$line"
mkdir -p "out/$dir"
cd "out/$dir"
if ls "${name}"_* "${name}".log >/dev/null 2>&1; then echo "output exists for $dir/$name, not overwriting" >&2; exit 1; fi
cmd="../../bin/rk4a_chain_v2 $args --threads 16 --out $name"
echo "$cmd" > "${name}.command.txt"
t0=$(date +%s)
$cmd > "${name}.log" 2>&1
echo "wall_seconds $(( $(date +%s) - t0 ))" >> "${name}.command.txt"
echo "task $SLURM_ARRAY_TASK_ID ($dir/$name) done $(date): $(tail -1 ${name}.log)"
