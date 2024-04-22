for p in Patient*; do
  echo $p
  cd $p
  ./final_run.sh
  cd ../
done
