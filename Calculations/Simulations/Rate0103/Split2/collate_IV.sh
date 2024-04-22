for s in Sim*; do
  echo $s
  cd $s
  grep 'Best likelihood' Run_ModelIE.out > Likelihoods.out
  grep 'Best set' Run_ModelVE.out >> Likelihoods.out
  grep Rate Run_ModelIE.out | tail -n2 > RatesIE.out
  cd ../
done
