../../run_ivy --input Sequences.fa --verb 1 --times Times.dat --model I --error 0.207 > Run_ModelIE.out
grep Index Run_ModelIE.out > Index.data
