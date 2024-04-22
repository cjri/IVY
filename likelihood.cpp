#include "ivy.h"
#include "io.h"
#include "likelihood.h"
#include "modelI.h"
#include "statespace.h"
#include <iostream>
#include <string>
#include <cstring>

void OptimiseMultiRateModel (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters, gsl_rng *rgen) {
    cout << "Optimise " << gseq_data.size() << "\n";
    vector<double> rates;
    for (int i=0;i<n_rates;i++) {
        rates.push_back(0.1);
    }
    vector<double> rates_best=rates;
    double error=0.1;
    double error_best=error;
    if (p.fix_error>-1) {
        error=p.fix_error;
        error_best=error;
    }
    double dx=0.001;
    double lL=-1e9;
    double lL_best=lL;
    int first=1;
    int reduced=0;
    int fail=0;
    if (gseq_data.size()>0) {
        
        vector<double> track;
        for (int it=0;it<100000;it++) {
            
            if (first==0) {
                if (lL>lL_best) {
                    rates_best=rates;
                    error_best=error;
                    double diff=lL-lL_best;
                    DoTrack(diff,track);
                    CheckTrack (reduced,dx,track);
                    lL_best=lL;
                    fail=0;
                    /*cout << "Better " << it << " ";
                    for (int i=0;i<rates.size();i++) {
                        cout << rates[i] << " ";
                    }
                    cout << lL << "\n";*/
                    
                } else {
                    fail++;
                    rates=rates_best;
                    error=error_best;
                }
            }
            if (fail==10000) {
                //cout << "End at " << it << "\n";
                break;
            }
            first=0;
            //Make a change to the parameters
            for (int i=0;i<rates.size();i++) {
                rates[i]=rates[i]+(gsl_rng_uniform(rgen)*dx)-(dx/2);
                if (rates[i]<0) {
                    rates[i]=-rates[i];
                }
            }
            
            if (p.fix_error==-1) {
                error=error+(gsl_rng_uniform(rgen)*dx)-(dx/2);
            }
            //Evaluate the likelihood
            lL=FindLikelihood (n_rates,rates,error,gseq_data);
            
        }
        //cout << "Done loop\n";
    }
    for (int i=0;i<n_rates;i++) {
        model_parameters.push_back(rates_best[i]);
    }
    model_parameters.push_back(error_best);
    model_parameters.push_back(lL_best);
}

void DoTrack (double diff, vector<double>& track) {
    if (track.size()<10) {
        track.push_back(diff);
    } else {
        vector<double> ntrack;
        for (int i=1;i<track.size();i++) {
            ntrack.push_back(track[i]);
        }
        track=ntrack;
        track.push_back(diff);
    }
}

void CheckTrack (int& reduced, double& dx, vector<double>& track) {
    if (track.size()>9&&reduced==0) {
        double tot=0;
        for (int i=0;i<track.size();i++) {
            tot=tot+track[i];
        }
        if (tot<0.01) {
            dx=dx/10;
            reduced++;
        }
    }
}

double FindLikelihood (int n_rates, vector<double>& rates, double error, const vector< vector<sample> >& gseq_data) {
    double lL=0;
    for (int i=0;i<gseq_data.size();i++) {
        for (int j=0;j<gseq_data[i].size();j++) {
            if (gseq_data[i][j].nfix>=0) {
                //cout << i << " " << j << "\n";
                double r=rates[i]*gseq_data[i][j].dt;
                //cout << "Here " << rates[i] << " " << gseq_data[i][j].dt << " " << gseq_data[i][j].nfix << " ";
                int done=0;
                if (r==0&&gseq_data[i][j].nfix>0) {
                    lL=lL-1e9;
                    //cout << "Now -1e9" << "\n";
                    done=1;
                }
                if (done==0) {
                    if (r>0&&gseq_data[i][j].nfix>=0) {
                        lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfix,r));
                        //cout << gsl_ran_poisson_pdf(gseq_data[i][j].nfix,r) << " " << log(gsl_ran_poisson_pdf(gseq_data[i][j].nfix,r)) << "\n";
                    }
                    lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,error));
                    //cout << "F " << gseq_data[i][j].nfluc << " " << error << " " << log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,error)) << " " << lL << "\n";
                    
                    /*if (gseq_data[i][j].nfluc>=0||gseq_data[i][j].xfluc.size()>0) {
                     if (gseq_data[i][j].xfluc.size()==0) {
                     lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,error));
                     } else {
                     //Account for uncertainty in the number of flucutations arising from ambiguous nucleotides
                     for (int j=0;j<gseq_data[i][j].xfluc.size();j++) {
                     lL=lL+gseq_data[i][j].xfluc[j]*log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc+j,error));
                     }
                     }
                     }*/
                }
            }
        }
    }
    return lL;
}


void UncertaintyMultiRateModel (run_params& p, int n_rates, int parameter, int direction, double best_likelihood, const vector< vector<sample> >& gseq_data, const vector<double>& initial_model_parameters, vector<double>& model_parameters, vector<double>& extreme_model_parameters, gsl_rng *rgen) {
    
    //Set up model
    for (int i=0;i<initial_model_parameters.size();i++) {
        model_parameters[i]=initial_model_parameters[i];
    }
    extreme_model_parameters=model_parameters;
    vector<double> model_parameters_previous=model_parameters;
    if (p.fix_error>-1) {
        model_parameters[n_rates]=p.fix_error;
        extreme_model_parameters=model_parameters;
    }

    double dx=0.001;
    double lL=-1e9;
    int first=1;
    for (int it=0;it<100000;it++) {
        if (first==0) {
            if (lL>best_likelihood-2) {
                if (direction==1&&model_parameters[parameter]>extreme_model_parameters[parameter]) {
                    extreme_model_parameters=model_parameters;
                }
                
                if (direction==0&&model_parameters[parameter]<extreme_model_parameters[parameter]) {
                    extreme_model_parameters=model_parameters;
                }
            } else {
                model_parameters=model_parameters_previous;
            }
        }
        first=0;
        //Make a change to the parameters
        for (int i=0;i<model_parameters.size()-1;i++) {
            if (i!=n_rates||p.fix_error==-1) {
                if (parameter==i) {
                    if (direction==1) {
                        model_parameters[i]=model_parameters[i]+(gsl_rng_uniform(rgen)*dx/2);
                    } else {
                        model_parameters[i]=model_parameters[i]-(gsl_rng_uniform(rgen)*dx/2);
                        if (model_parameters[i]<0) {
                            model_parameters[i]=-model_parameters[i];
                        }
                    }
                } else {
                    model_parameters[i]=model_parameters[i]+(gsl_rng_uniform(rgen)*dx)-(dx/2);
                }
            }
        }
        
        //Evaluate the likelihood
        lL=FindLikelihood (n_rates,model_parameters,model_parameters[n_rates],gseq_data);
        
    }
}

