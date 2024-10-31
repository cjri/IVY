../../run_ivy --input Sequences.fa --times Times.dat --model V --error 0.207 --set 62 --uncertainty 1 --verb 1 --systematic 0 > Run_ModelVE_U.out
grep Index Run_ModelVE_U.out > Index.data

