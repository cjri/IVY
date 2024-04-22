for p in Patient*; do
  echo $p
  cd $p
  ./run_calcs.sh
  cd ../
done
