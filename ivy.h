
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


struct run_params {
    int verb; //Verbose output
    int seed; //For random number generator
    int fix; //Include ambiguous nucleotides
    int sets; //Number of sets into which we divide the sequences
    string model;
    int systematic; //Flag for systematic search
    int uncertainty; //Flag to calculate uncertainty in parameters
    int chosen_set; //Specify a set to look at e.g. in ModelV
    double fix_error;
    string seqs_file;
    string times_file;
};

struct sparseseq {
    vector<int> locus;
    vector<char> allele;
};

struct sample {
    int dt;
    double nfix;
    double nfluc;
    vector<double> xfluc;
};

struct modelstore {
    vector<int> start_seq;  //Initial sequence for this model
    vector<double> rates; //Rates of evolution
    double error;
    double lL;
    double lL_rel; //Relative likelihood.  Not used at the moment
    int index; //Index within the set of data uncertainties
    int splittime; //Time of split of populations (where used)
};




