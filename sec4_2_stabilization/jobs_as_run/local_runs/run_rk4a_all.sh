#!/bin/bash
# §4.3 action-angle runs redone with rk4a_chain (schemes4 RK4_A + Cartesian bath). Same cases, parameters and seeds as
# run_aa_core.sh + run_aa_angles_sweep.sh (aa_chain); the floor test is dropped (no floor in rk4a_chain).
set -e
B=$(cd "$(dirname "$0")" && pwd)/rk4a_chain
R="$1"; mkdir -p "$R"/{val,case1,drain,lock,ness,angles,sweep}
run() { local d=$1 name=$2; shift 2; (cd "$R/$d" && $B --threads 2 --out $name "$@" > $name.log); echo "$d/$name $(tail -1 $R/$d/$name.log)"; }
run val n4_eqT5_full --n 4 --T1 5 --Tn 5 --dt 1e-3 --burn 50 --measure 1000 --streams 128 --seed 11
run val n3_hold2_T5_5 --n 3 --T1 5 --Tn 5 --hold 2 --dt 1e-3 --burn 50 --measure 2000 --streams 256 --seed 12
run val n3_hold2_T2_8_dt5e-4 --n 3 --T1 2 --Tn 8 --hold 2 --dt 5e-4 --burn 50 --measure 1000 --streams 256 --seed 13
run val n3_hold2_T2_8_dt1e-3 --n 3 --T1 2 --Tn 8 --hold 2 --dt 1e-3 --burn 50 --measure 1000 --streams 256 --seed 13
run val n4_hold2_T2_8_dt5e-4 --n 4 --T1 2 --Tn 8 --hold 2 --dt 5e-4 --burn 50 --measure 1000 --streams 256 --seed 14
for TT in "2 8" "5 5" "8 2"; do set -- $TT
  run case1 n3_hold2_T$1_$2 --n 3 --T1 $1 --Tn $2 --hold 2 --dt 1e-3 --burn 50 --measure 2000 --streams 256 --seed 21
  run case1 n4_hold2_T$1_$2 --n 4 --T1 $1 --Tn $2 --hold 2 --dt 1e-3 --burn 50 --measure 2000 --streams 256 --seed 22
  run case1 n4_hold23_T$1_$2 --n 4 --T1 $1 --Tn $2 --hold 2,3 --dt 1e-3 --burn 50 --measure 2000 --streams 256 --seed 23
done
run ness n3_T2_8 --n 3 --T1 2 --Tn 8 --dt 1e-3 --burn 50 --measure 1000 --streams 64 --seed 31
for p in "10 1e-3" "30 2e-4" "100 1e-4" "200 5e-5"; do set -- $p
  run drain n3_E$1 --n 3 --T1 2 --Tn 8 --init 0.01,$1,0.01 --dt $2 --burn 0 --measure 10 --interval 0.002 --streams 256 --seed 32
done
for g in 0.01 0.1 0.2 0.5 1.0; do
  run lock n3_g$g --n 3 --T1 0.1 --Tn 0.1 --gamma $g --init 1e-6,20,1e-6 --dt 2e-5 --burn 0 --measure 0.6 --interval 2e-4 --streams 256 --seed 33
done
for TT in "2 8" "5 5" "8 2"; do set -- $TT
  run angles n3_hold2_T$1_$2 --n 3 --T1 $1 --Tn $2 --hold 2 --dt 1e-3 --burn 50 --measure 200 --interval 0.02 --streams 64 --seed 41
  run angles n4_hold2_T$1_$2 --n 4 --T1 $1 --Tn $2 --hold 2 --dt 1e-3 --burn 50 --measure 200 --interval 0.02 --streams 64 --seed 42
  run angles n4_hold23_T$1_$2 --n 4 --T1 $1 --Tn $2 --hold 2,3 --dt 1e-3 --burn 50 --measure 200 --interval 0.02 --streams 64 --seed 43
done
for a in 0.5 1 2 5 10; do for b in 0.5 1 2 5 10; do
  run sweep n3_hold2_T${a}_${b} --n 3 --T1 $a --Tn $b --hold 2 --dt 1e-3 --burn 50 --measure 500 --streams 128 --seed 51
done; done
echo ALLDONE
