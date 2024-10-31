../../run_ivy --input Sequences.fa --times Times.dat --model Y --error 0.207 --set 2421 --uncertainty 1 --verb 1 --systematic 0 > Run_ModelY_U.out
grep Index Run_ModelY_U.out > Index.data
