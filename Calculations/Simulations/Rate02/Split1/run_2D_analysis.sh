for s in Sim*; do
  echo $s
  cd $s
  /Users/chris/Documents/Projects/MPhil_WithinHost/Writeup/MappingCode/run_align --input PopSeqs.dat --dim 2
  mv Output_points.dat Output_points2.dat
  mv Distances.dat Distances2.dat
  /Users/chris/Documents/Projects/MPhil_WithinHost/Writeup/MappingCode/run_align --input PopSeqs.dat --dim 3
  mv Output_points.dat Output_points3.dat
  mv Distances.dat Distances3.dat
  cd ../
done
