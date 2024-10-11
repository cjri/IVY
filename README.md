# IVY

The IVY code fits models of evolution to consensus sequence data describing viral infection in a single host.  The results of IVY can be used to identify the minimum number of viral subpopulations within an individual that are needed to explain the sequences observed from that individual.

Documentation for the use of this code is found in the file Instructions.txt

This code accompanies a manuscript, accessible here: https://www.biorxiv.org/content/10.1101/2024.05.31.596818v1

**Compilation and use**

IVY is written in C++ and makes use of the GSL software library.  The CC_FLAGS and LD_FLAGS lines of the accompanying Makefile will likely need to be changed on your machine to facilitate use, specifically the path to where the GSL library has been installed on your machine.  Running the command

make

in the directory containing the code will carry out compilation, once everything is in place.

**Reproducibility and Example data**

The Calculations directory contains data and scripts facilitating the reproduction of the calculations run in the manuscript.  Within this directory the file Notes.txt provides details about running these scripts.  We note that a full reproduction of the results may take a few days of computational time on a laptop.

The directory Calculations/Simulations contains data and scripts facilitating the reproduction of the generation and analysis of simulated data.  The Mathematica notebooks in this directory generate sequence data.  While use of these notebooks requires a Mathematica license, the simulated data generated by these notebooks is provided in subdirectories.  For example the directory Rate0102/Split2 contains simulated data from cases of infection with two subpopulations, with respective rates of evolution 0.1 and 0.2 substitutions per day.  The script run_modelV.sh in this directory runs the inference model V on the data.

Scripts as provided need some basic editing in order to work:

1.  You will need to modify the path to the run_ivy executable to fit with its location on your machine.

2.  The scripts make use of the GNU Parallel package.  You can either install this, or modify the script to run the code in serial as follows:

for s in Sim*; do
    cd $s
    <PATH>/run_ivy --input PopSeqs.dat --times Times.in --model V --error 0.207 > Run_ModelVE.out
    cd ../
done


