while read line
  do
  echo $line
  cd $line
  ./run_grid.sh &
  cd ../
done < Grid_calcs.in
