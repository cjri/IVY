#include "ivy.h"
#include "modelI.h"
#include "model.h"
#include "io.h"
#include "model_tools.h"
#include <iostream>
#include <string>
#include <cstring>

void ModelSinglePopulation (run_params& p, vector<int> times, vector< vector<double> > varbin, gsl_rng *rgen) {
    
    //Set up a starting point - at time zero
    AddSingleOrigin (times,varbin);
    
    //Thoughts here about the initial point.  It could be equal to the first time point. If we assume an infinite sites
    //model then any fixations will be in sites not observed here as polymorphic: The consequence of this will be more
    //fixations in the time before the first observed sample.
    
    //Identify vectors against which to compare which variants are fixations and which are fluctuations
    vector< vector<double> > constant;
    vector< vector<double> > fixes;
    MakeConstantFixes(p,varbin,constant,fixes);

    //Test variant data against constant and fix vectors
    vector<int> fixpos;
    vector<int> flucpos;
    FindFixationSites (p,varbin,constant,fixes,fixpos,flucpos);

    //Calculate times of fixations - first observation of a 1
    vector<int> fixtimes;
    FindFixationTimes (p,fixpos,varbin,fixtimes);

    //Use the fixation times to find a new category of sites that fix in the last time-point
    //Not clear whether these are fixations or fluctuations: Need to adjust optimisation to account for this
    //What we do is consider the possible numbers of which are which.  Precise allocation doesn't matter.
    vector<int> qfixpos;
    FindQFixPositions (p,fixpos,fixtimes,varbin,qfixpos);
    
    //Find numbers of fluctuations and fixations in each sample point
    //Note: The numbers of these depends upon the ambiguity in the qfixpos.
    //Note 2: The precise location of the qfixpos doesn't matter: There is just uncertainty in the number of fluctuations or fixations
    
    //How do we incorporate uncertainty into the seq_data?
    //The number of fixes or fluctuations has to be an integer
    //1.  Make vectors of non-binary probabilities of fixes and flucs.
    //2.  Use the vector to calculate probabilities of 0, 1, 2, 3 fixes/flucs
    //Make a composite likelihood of P(#fluc) x L(rate,error|#fluc)

    //Note: Where an ambigous nucleotide occurs in an otherwise non-polymorphic site we ignore it.  Need an ACGT to call a variant.
    
    int nq=qfixpos.size();
    vector<sample> seq_data;
    SetupSequenceData (p,times,fixtimes,fixpos,flucpos,varbin,seq_data);
    //Note: We could get nq from the basic definition of seq_data but it might be useful to know which are which later
    cout << "Size " << fixpos.size() << " " << qfixpos.size() << "\n";
    cout << "Nq " << nq << "\n";

    //Updated to here: The routine SetupSequence data now sets up the data to account for uncertain positions:
    //The probability of k additional fluctuations is given by seq_data[i].xfluc[k]
    
    //We want to optimise over the number of fixations prior to the first observed time point.
    //Do this in a simple way for the single-line model: Assume infinite sites and just add #fixations
    
    vector<double> model_parameters;
    SetupModelParameters(model_parameters);
    vector<double> model_parameters_best=model_parameters;
    vector<sample> seq_data_best;
    vector< vector<sample> > seq_data_record;
    vector< vector<double> > model_parameters_record;
    vector<double> likelihood_record;
    FindBestModelFirstLast (p,nq,seq_data,seq_data_best,model_parameters,model_parameters_best,seq_data_record,fixpos,qfixpos,model_parameters_record,rgen);  //Accounts for uncertainty in original sequence and final sample

    
    cout << "Best sequence input\n";
    //N.B. xfluc data not actually accounted for in the likelihood calculation: Commented out in likelihood.cpp
    for (int j=0;j<seq_data_best.size();j++) {
        cout << j << " " << seq_data_best[j].dt << " " << seq_data_best[j].nfix << " " << seq_data_best[j].nfluc << " " << seq_data_best[j].xfluc.size() << " ";
        for (int j=0;j<seq_data_best[j].xfluc.size();j++) {
            cout << seq_data_best[j].xfluc[j] << " ";
        }
        cout << "\n";
    }

    //Uncertainty calculation in the single model.  Use the optimum dataset for this?  Could try this initially but should really consider the whole dataset
    //Perhaps if we make a vector of datasets that are within 2 log likelihood units of the best one at their optimum.
    //Can do this by recording all of the locally optimal likelihoods.
    //Do this recording in FindBestModelLast
    
    
    //Filter sequence inputs by closeness to best likelihood.  Keep within two units of best likelihood
    //We also want details of the optimal parameters for each dataset
    double best_likelihood=-1e9;
    GetOptimisationData (best_likelihood,seq_data_record,model_parameters_record);
    
    //Calculate parameter uncertainties for each of the sequence inputs
    cout << "Best likelihood " << best_likelihood << "\n";
    cout << "Rate " << model_parameters_best[0] << "\n";
    cout << "Error " << model_parameters_best[1] << "\n";
    vector< vector<double> > limits;
    InitialiseLimits (p,model_parameters_best,limits);
    vector<double> extreme_model_parameters;
    SingleRateModelExtremes (p,best_likelihood,seq_data_record,model_parameters_record,extreme_model_parameters,limits,rgen);
    WriteLimits(limits);
    
}

void AddSingleOrigin (vector<int>& times, vector< vector<double> >& varbin) {
    vector<double> varbin_init=varbin[0];  //Start point is first sample; allow for a number of substitutions to be unobserved
    //Accounting for these is done later...
    //Edit varbin information
    vector< vector<double> > varbin_temp;
    varbin_temp.push_back(varbin_init);
    for (int i=0;i<varbin.size();i++) {
        varbin_temp.push_back(varbin[i]);
    }
    varbin=varbin_temp;
    //Edit times vector
    vector<int> times_temp;
    times_temp.push_back(0);
    for (int i=0;i<times.size();i++) {
        times_temp.push_back(times[i]);
    }
    times=times_temp;
}


void FindBestModelFirstLast (run_params& p, int nq, const vector<sample>& seq_data, vector<sample>& seq_data_best, vector<double>& model_parameters, vector<double>& model_parameters_best, vector< vector<sample> >& seq_data_record, vector<int>& fixpos, vector<int>& qfixpos, vector< vector<double> >& model_parameters_record, gsl_rng *rgen) {
    int pre_subs=-1;
    int better=1;
    int first=1;
    double lL_best=-1e6;
    do { //Assume a smooth pattern in the likelihood function - keep going until decrease in likelihood
        if (first==0) {
            better=0;
            if (model_parameters_best[2]>lL_best) {
                better=1;
                lL_best=model_parameters_best[2];
            }
        }
        first=0;
        pre_subs++;

        if (better==1) {
            //Construct an array of seq_data objects, spanning all of the possibilities for the last time point...
            vector< vector<sample> > seq_data_array;
            vector<int> removed;
            ConstructSequenceDataArray(nq,seq_data,removed,seq_data_array);
            for (int i=0;i<seq_data_array.size();i++) {
                seq_data_array[i][1].nfix=seq_data_array[i][1].nfix+pre_subs;
            }
            
            //Find best model parameters
            int index_best=-1;
            FindBestModelParameters (p,pre_subs,seq_data_best,seq_data_array,index_best,model_parameters,model_parameters_best,seq_data_record,fixpos,qfixpos,removed,model_parameters_record,rgen); //Rename as FindBestSingleRateModelParameters?
        }
    }
    while (better==1);
}

void OptimiseSingleRateModel (run_params& p, const vector<sample>& seq_data, vector<double>& model_parameters, gsl_rng *rgen) {
    double rate=0.1;
    double error=0.1;
    double rate_best=rate;
    double error_best=error;
    if (p.fix_error>-1) {
        error=p.fix_error;
        error_best=error;
    }
    double dx=0.001;
    double lL=-1e9;
    double lL_best=lL;
    int first=1;
    
    cout << "Data\n";
    for (int i=0;i<seq_data.size();i++) {
        cout << seq_data[i].dt << " " << seq_data[i].nfix << " " << seq_data[i].nfluc << "\n";
    }
    
    for (int it=0;it<100000;it++) {
        if (first==0) {
            if (lL>lL_best) {
                rate_best=rate;
                error_best=error;
                lL_best=lL;
                //cout << "Better " << rate << " " << error << " " << lL << "\n";
            } else {
                rate=rate_best;
                error=error_best;
            }
        }
        first=0;
        //Make a change to the parameters
        rate=rate+(gsl_rng_uniform(rgen)*dx)-(dx/2);
        if (rate<0) {
            rate=1e-9;
        }
        if (p.fix_error==-1) {
            error=error+(gsl_rng_uniform(rgen)*dx)-(dx/2);
        }
        
        //Evaluate the likelihood
        
        
        //Original code
        //Evaluate the likelihood
        lL=GetSingleLikelihood (rate,error,seq_data);
    }
    model_parameters[0]=rate_best;
    model_parameters[1]=error_best;
    model_parameters[2]=lL_best;
    cout << "Optimised rate " << rate_best << " error " << error_best << " lL " << lL_best << "\n";

}

double GetSingleLikelihood (double& rate, double& error, const vector<sample>& seq_data) {
    double lL=0;
    for (int i=0;i<seq_data.size();i++) {
        if (seq_data[i].nfix>=0) {
            int done=0;
            double r=rate*seq_data[i].dt;
            //cout << "Here " << rate << " " << seq_data[i].dt << " " << seq_data[i].nfix << " ";
            
            if (r==0&&seq_data[i].nfix>0) {//dt is zero: Will happen if we put an initial time-point at the beginning and the first time is zero
                lL=lL-1e9;
                //cout << "Now -1e9" << "\n";
                done=1;
            }
            if (done==0) {
                if (r>0&&seq_data[i].nfix>=0) {
                    lL=lL+log(gsl_ran_poisson_pdf(seq_data[i].nfix,r));
                    // cout << gsl_ran_poisson_pdf(seq_data[i].nfix,r) << " " << log(gsl_ran_poisson_pdf(seq_data[i].nfix,r)) << "\n";
                    
                }
                lL=lL+log(gsl_ran_poisson_pdf(seq_data[i].nfluc,error));
                //cout << "F " << seq_data[i].nfluc << " " << error << " " << log(gsl_ran_poisson_pdf(seq_data[i].nfluc,error)) << " " << lL << "\n";
                
                /*
                 if (seq_data[i].xfluc.size()==0) {
                 lL=lL+log(gsl_ran_poisson_pdf(seq_data[i].nfluc,error));
                 } else {
                 //Account for uncertainty in the number of flucutations arising from ambiguous nucleotides
                 for (int j=0;j<seq_data[i].xfluc.size();j++) {
                 lL=lL+seq_data[i].xfluc[j]*log(gsl_ran_poisson_pdf(seq_data[i].nfluc+j,error));
                 }
                 }*/
            }
        }
    }
    return lL;
}

void GetOptimisationData (double& best_likelihood, vector< vector<sample> >& seq_data_record, vector< vector<double> >& model_parameters_record) {
    //Get best likelihood
    for (int i=0;i<model_parameters_record.size();i++) {
        if (model_parameters_record[i][2]>best_likelihood) {
            best_likelihood=model_parameters_record[i][2];
        }
    }
    //Get sequence inputs which give likelihoods close to this
    vector< vector<sample> > seq_data_record_temp;
    for (int i=0;i<model_parameters_record.size();i++) {
        if (model_parameters_record[i][2]>best_likelihood-2) {
            seq_data_record_temp.push_back(seq_data_record[i]);
        }
    }
    seq_data_record=seq_data_record_temp;
}


void SingleRateModelExtremes (run_params& p, double best_likelihood, vector< vector<sample> >& seq_data_record, vector< vector<double> >& model_parameters_record, vector<double>& extreme_model_parameters, vector< vector<double> >& limits,gsl_rng *rgen) {
    for (int i=0;i<seq_data_record.size();i++) {
        UncertaintySingleRateModel (p,0,0,best_likelihood,seq_data_record[i],model_parameters_record[i], model_parameters_record[i],extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]<limits[0][0]) {
            limits[0][0]=extreme_model_parameters[0];
        }
        UncertaintySingleRateModel (p,0,1,best_likelihood,seq_data_record[i],model_parameters_record[i], model_parameters_record[i],extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]>limits[0][1]) {
            limits[0][1]=extreme_model_parameters[0];
        }
        UncertaintySingleRateModel (p,1,0,best_likelihood,seq_data_record[i],model_parameters_record[i], model_parameters_record[i],extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]<limits[1][0]) {
            limits[1][0]=extreme_model_parameters[1];
        }
        UncertaintySingleRateModel (p,1,1,best_likelihood,seq_data_record[i],model_parameters_record[i], model_parameters_record[i],extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]>limits[1][1]) {
            limits[1][1]=extreme_model_parameters[1];
        }

    }

}

void UncertaintySingleRateModel (run_params& p, int parameter, int direction, double best_likelihood, const vector<sample>& seq_data, vector<double> initial_model_parameters, vector<double>& model_parameters, vector<double>& extreme_model_parameters, gsl_rng *rgen) {
    //cout << " Uncertainty model " << best_likelihood << "\n";
    double rate=initial_model_parameters[0];
    double error=initial_model_parameters[1];
    extreme_model_parameters=model_parameters;
    double rate_previous=rate;
    double error_previous=error;
    if (p.fix_error>-1) {
        error=p.fix_error;
        error_previous=error;
    }
    double dx=0.001;
    double lL=-1e9;
    double lL_best=best_likelihood;
    int first=1;
    for (int it=0;it<100000;it++) {
        if (first==0) {
            if (lL>lL_best-2) {
                //cout << lL << " " << rate << " " << error << "\n";
                model_parameters[0]=rate;
                model_parameters[1]=error;
                model_parameters[2]=lL;
                if (direction==1&&model_parameters[parameter]>extreme_model_parameters[parameter]) {
                    extreme_model_parameters[0]=rate;
                    extreme_model_parameters[1]=error;
                    extreme_model_parameters[2]=lL;
                    rate_previous=rate;
                    error_previous=error;
                }
                if (direction==0&&model_parameters[parameter]<extreme_model_parameters[parameter]) {
                    extreme_model_parameters[0]=rate;
                    extreme_model_parameters[1]=error;
                    extreme_model_parameters[2]=lL;
                }
            } else {
                rate=rate_previous;
                error=error_previous;
            }
        }
        first=0;
        //Make a change to the parameters
        if (p.fix_error==-1) {
            if (parameter==0&&direction==1) {
                rate=rate+(gsl_rng_uniform(rgen)*dx/2);
                error=error+(gsl_rng_uniform(rgen)*dx)-(dx/2);
            }
            if (parameter==0&&direction==0) {
                rate=rate-(gsl_rng_uniform(rgen)*dx/2);
                error=error+(gsl_rng_uniform(rgen)*dx)-(dx/2);
            }

            if (parameter==1&&direction==1) {
                rate=rate+(gsl_rng_uniform(rgen)*dx)-(dx/2);
                error=error+(gsl_rng_uniform(rgen)*dx/2);
            }
            if (parameter==1&&direction==0) {
                rate=rate+(gsl_rng_uniform(rgen)*dx)-(dx/2);
                error=error-(gsl_rng_uniform(rgen)*dx/2);
            }
        } else {
            if (parameter==0) {
                if (direction==1) {
                    rate=rate+(gsl_rng_uniform(rgen)*dx/2);
                } else {
                    rate=rate-(gsl_rng_uniform(rgen)*dx/2);
                    if (rate<0) {
                        rate=-rate;
                    }
                }
            }
        }
        
        //Evaluate the likelihood
        lL=GetSingleLikelihood (rate,error,seq_data);
        //cout << "Likelihood " << lL << "\n";
    }
  //  cout << "Extreme parameters " << extreme_model_parameters[0] << " " << extreme_model_parameters[1] << " " << extreme_model_parameters[2] << "\n";
    
}


void SetupModelParameters (vector<double>& model_parameters) {
    for (int i=0;i<3;i++) {
        model_parameters.push_back(-1e6);
    }
}

void FindBestModelParameters (run_params& p, int pre_subs, vector<sample>& seq_data_best, const vector< vector<sample> >& seq_data_array, int& index_best, vector<double>& model_parameters, vector<double>& model_parameters_best, vector< vector<sample> >& seq_data_record, vector<int>& fixpos, vector<int>& qfixpos, vector<int>& removed, vector< vector<double> >& model_parameters_record, gsl_rng *rgen) {
    //Finds the best across an array of points
    for (int i=0;i<seq_data_array.size();i++) {
        /*if (p.verb==1) {
            cout << "Sequence input\n";
            for (int j=0;j<seq_data_array[i].size();j++) {
                cout << j << " " << seq_data_array[i][j].dt << " " << seq_data_array[i][j].nfix << " " << seq_data_array[i][j].nfluc << " " << seq_data_array[i][j].xfluc.size() << " ";
                for (int j=0;j<seq_data_array[i][j].xfluc.size();j++) {
                    cout << seq_data_array[i][j].xfluc[j] << " ";
                }
                cout << "\n";
            }
        }*/
        OptimiseSingleRateModel (p,seq_data_array[i],model_parameters,rgen);
        seq_data_record.push_back(seq_data_array[i]);
        model_parameters_record.push_back(model_parameters);
        if (model_parameters[2]>-1e+9) {
            cout << "Index " << i << " Rate " << model_parameters[0] << " Prior Substitutions " << pre_subs << " Fixes " << fixpos.size()-qfixpos.size() << " ";
            for (int k=0;k<fixpos.size();k++) {
                int match=0;
                for (int l=0;l<qfixpos.size();l++) {
                    if (qfixpos[l]==fixpos[k]) {
                        match=1;
                    }
                }
                if (match==0) {
                    cout << fixpos[k] << " ";
                }
            }
            cout << "Provisional " << qfixpos.size() << " ";
            for (int l=0;l<qfixpos.size();l++) {
                cout << qfixpos[l] << " ";
            }

            //Want to output how many are removed...
            cout << " Removed " << removed[i];
            cout << " Likelihood " << model_parameters[2] << "\n";
        }
        if (model_parameters[2]>model_parameters_best[2]) {
            model_parameters_best=model_parameters;
            index_best=i;
            seq_data_best=seq_data_array[index_best];
        }
    }
    if (p.verb==1) {
        if (index_best>-1) {
            cout << "Best index " << index_best << " ";
            for (int i=0;i<model_parameters_best.size();i++) {
                cout << model_parameters_best[i] << " ";
            }
            cout << "\n";
            
            /*cout << "Corresponding sequence input\n";
            for (int j=0;j<seq_data_array[index_best].size();j++) {
                cout << j << " " << seq_data_array[index_best][j].dt << " " << seq_data_array[index_best][j].nfix << " " << seq_data_array[index_best][j].nfluc << " " << seq_data_array[index_best][j].xfluc.size() << " ";
                for (int j=0;j<seq_data_array[index_best][j].xfluc.size();j++) {
                    cout << seq_data_array[index_best][j].xfluc[j] << " ";
                }
                cout << "\n";
            }*/

        }
    }
}



