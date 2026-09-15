#!/bin/bash
# T1 = 1, 2, 3, ..., 10
# Tn = 1, 2, 3, ..., 10
# n = 50, Case 0

echo "===T1=[1-10], Tn=[1-10], n=50 ==="

n=50
case_num=0

job=0
for T1 in 1 2 3 4 5 6 7 8 9 10; do
    for Tn in 1 2 3 4 5 6 7 8 9 10; do
        job=$((job + 1))
        deltaT=$((T1 - Tn))
        echo "Run $job/100: T1=$T1, Tn=$Tn, deltaT=$deltaT"
        ./EnergyCascade $case_num $T1 $Tn $n
    done
done

echo ""
echo "=== completed ==="
echo "Results: flux_summary_temp.csv"