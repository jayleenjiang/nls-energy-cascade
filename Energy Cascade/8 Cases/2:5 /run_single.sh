#!/bin/bash
# job_id: 1-100
# job_id = (T1-1)*10 + Tn
#   for i in $(seq 1 100); do sbatch --wrap="./run_single_grid_100.sh $i"; done

JOB_ID=$1

if [ -z "$JOB_ID" ]; then
    echo "Usage: $0 <job_id (1-100)>"
    exit 1
fi

T1=$(( (JOB_ID - 1) / 10 + 1 ))
Tn=$(( (JOB_ID - 1) % 10 + 1 ))
n=50
case_num=0

echo "Job $JOB_ID: T1=$T1, Tn=$Tn, n=$n"
./EnergyCascade $case_num $T1 $Tn $n
echo "Job $JOB_ID done"
