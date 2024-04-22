for s in Sim*; do
  echo $s
  cd $s
  /Users/christopher.illingworth/Documents/GitHub/IVY/run_ivy --input PopSeqs.dat --times Times.in --model I --error 0.207 > Run_ModelIE.out
  cd ../
done
