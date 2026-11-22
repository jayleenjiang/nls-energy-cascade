G=${G:-$PWD}   # as run: a scratch folder of the cloud workspace; loc/ holds the local_test_v2 build (now sec4_3_return)
L=$G/loc
cd $G/g005/runs
$L/local_test_v2 ex ex_c0_neq_g0.05.json 2 8 0.05 1e-3 0 0 0.01 501 4096 161 2 1 1 1 0 0
$L/local_test_v2 ex ex_c0_eq_g0.05.json  5 5 0.05 1e-3 0 0 0.01 501 4096 162 2 1 1 1 0 0
$L/local_test_v2 ex ex_c2_neq_g0.05.json 2 8 0.05 1e-3 0 0 0.01 501 4096 165 2 2 1 1 0 0
$L/local_test_v2 ex ex_c2_eq_g0.05.json  5 5 0.05 1e-3 0 0 0.01 501 4096 166 2 2 1 1 0 0
$G/sampler_n3 neq_g005.bin 2 8 0.05 1e-3 200 5000 0.2 4 7501 2
$G/sampler_n3 eq_g005.bin  5 5 0.05 1e-3 200 5000 0.2 4 7502 2
cd $G/g001/runs
$L/local_test stat stat_eq_g0.01.json  5 5 0.01 1e-3 1000 5000 0.01 500 64 152 2
cd $G/g005/runs
$L/local_test stat stat_neq_g0.05.json 2 8 0.05 1e-3 200 5000 0.01 500 64 171 2
$L/local_test stat stat_eq_g0.05.json  5 5 0.05 1e-3 200 5000 0.01 500 64 172 2
echo ALLDONE
