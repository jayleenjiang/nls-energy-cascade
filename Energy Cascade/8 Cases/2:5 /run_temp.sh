#!/bin/bash
# run_temp.sh
# T1 = 1-10, Tn = 1-10, n = 50

echo "=== Starting 100 parallel jobs ==="

n=50
case_num=0

for T1 in 1 2 3 4 5 6 7 8 9 10; do
    for Tn in 1 2 3 4 5 6 7 8 9 10; do
        ./EnergyCascade $case_num $T1 $Tn $n > "log_T1_${T1}_Tn_${Tn}.txt" 2>&1 &
    done
done

echo "All jobs submitted"
echo "Waiting for completion..."
wait

echo "Completed"
