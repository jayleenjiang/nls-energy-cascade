L=../../loc
$L/local_test_v2 ex ex_c0_neq_g0.02.json 2 8 0.02 1e-3 0 0 0.01 501 4096 121 2 1 1 1 0 0
$L/local_test_v2 ex ex_c0_eq_g0.02.json  5 5 0.02 1e-3 0 0 0.01 501 4096 122 2 1 1 1 0 0
$L/local_test_v2 ex ex_c2_neq_g0.02.json 2 8 0.02 1e-3 0 0 0.01 501 4096 125 2 2 1 1 0 0
$L/local_test_v2 ex ex_c2_eq_g0.02.json  5 5 0.02 1e-3 0 0 0.01 501 4096 126 2 2 1 1 0 0
../../sampler_n3 neq_g002.bin 2 8 0.02 1e-3 500 5000 0.5 4 7201 2
../../sampler_n3 eq_g002.bin  5 5 0.02 1e-3 500 5000 0.5 4 7202 2
$L/local_test stat stat_neq_g0.02.json 2 8 0.02 1e-3 500 5000 0.01 500 64 131 2
$L/local_test stat stat_eq_g0.02.json  5 5 0.02 1e-3 500 5000 0.01 500 64 132 2
echo ALLDONE
