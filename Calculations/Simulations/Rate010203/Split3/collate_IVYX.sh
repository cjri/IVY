for s in Sim*; do
  echo $s
  cd $s
  grep 'Best likelihood' Run_ModelIE.out > Likelihoods.out
  grep 'Best set' Run_ModelVE.out >> Likelihoods.out
  grep 'Best set' Run_ModelYE.out >> Likelihoods.out
  grep 'Best set' Run_ModelXE.out >> Likelihoods.out
  grep Rate Run_ModelIE.out | tail -n2 > RatesIE.out
  grep -A1 'Best parameters' Run_ModelVE.out | tail -n1 > RatesVE.out
  grep -A1 'Best parameters' Run_ModelYE.out | tail -n1 > RatesYE.out
  grep -A1 'Best parameters' Run_ModelXE.out | tail -n1 > RatesXE.out
  cd ../
done
