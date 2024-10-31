#include "ivy.h"
#include "model.h"
#include "clustering.h"
#include "io.h"
#include "likelihood.h"
#include "model_tools.h"
#include "process_sequences.h"
#include "statespace.h"
#include <iostream>
#include <string>
#include <cstring>

void RunModel (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen) {

    WriteVarbin(varbin);
    
    //Find samples that share k variants - put into clusters.  In the below k=2: Allows for a single error position
    vector< vector<int> > clusters;
    FindClusters(p,2,varbin,clusters);

    cout << "Clusters\n";
    
    //Note here: We find all of the sets made up of clusters, keeping clustered data points together.
    vector< vector< vector<int> > > sets;
    CompileSets (p,clusters,sets);
    cout << "Number of sets " << sets.size() << "\n";
    
    if (p.systematic==1) {
        //Systematic search of all possible sets
        int best_set=-1;
        double best_log=-1e9;
        for (int i=0;i<sets.size();i++) { //Systematic search through all possible sets.  Want option to explore one set out of all of them, with uncertainty
            cout << "Run systematic search of sets\n";
            vector<modelstore> outputs;
            ExploreSet (p,i,best_set,best_log,times,outputs,varbin,sets,rgen);
        }
        cout << "Best set " << best_set << " likelihood " << best_log << "\n";
    } else {
        //Look at specific set   N.B. Come back to this as there are two uncertainty calculations.  Should make these consistent.
        double best_log=-1e9;
        cout << "Explore specified set\n";
        vector<modelstore> outputs;
        ExploreSet (p,p.chosen_set,p.chosen_set,best_log,times,outputs,varbin,sets,rgen);
        
        //Find uncertainty in parameters
        if (p.uncertainty==1) {
            cout << "Do uncertainty here\n";
            double maxL=best_log;
            RefineModels(maxL,outputs);
            
            //Set up best model parameters
            vector<double> model_parameters_best;
            GetBestModelParameters (maxL,outputs,model_parameters_best);

            vector< vector<int> > gtimes;
            vector< vector< vector<double> > > gvarbin;
            AssignDataToSetsG (p,p.chosen_set,sets,times,varbin,gtimes,gvarbin);
            vector< vector< vector<double> > > gvarbin_orig=gvarbin;
            vector< vector<int> > gtimes_orig=gtimes;

            //Calculate uncertainty in each of the models
            vector< vector<double> > limits;
            InitialiseLimits (p,model_parameters_best,limits);
            //Find limits
            ModelRateExtremes (p,maxL,gvarbin_orig,gtimes_orig,outputs,limits,rgen);
            
            WriteLimits(limits);
        }
        
        if (p.uncertainty==2) {
            double maxL=FindMaxL(outputs);
            vector< vector<double> > acceptable;
            
            ModelRateExploration (p,p.chosen_set,maxL,times,varbin,sets,outputs,acceptable);
            ofstream acc_file;
            acc_file.open("Accepted_points.dat");
            for (int j=0;j<acceptable.size();j++) {
                for (int k=0;k<acceptable[j].size();k++) {
                    acc_file << acceptable[j][k] << " ";
                }
                acc_file << "\n";
            }
        }
    }

}

void ExploreSet (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, vector<modelstore>& outputs, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen) {
    //Go through all sets and find the one producing the maximum likelihood reconstruction.
    //We don't need to do the uncertainty calculation in this case.
    
    cout << "Set " << i << "\n";
    
    //Calculate penalty for variants found in more than one population.
    double Lpen=FindPenalty(i,times,varbin,sets);
    cout << "Sharing Penalty " << Lpen << "\n";
    
    double maxL=-1e6;

    //Store starting data in a vector of sets
    vector< vector<int> > gtimes;
    vector< vector< vector<double> > > gvarbin;
    AssignDataToSetsG (p,i,sets,times,varbin,gtimes,gvarbin);
    vector< vector< vector<double> > > gvarbin_orig=gvarbin;
    vector< vector<int> > gtimes_orig=gtimes;
    WriteGVarbin(gtimes,gvarbin);
        
    //Check times compliance - max set size >= 3
    int max=CheckTimes(times,gtimes);
    if (max>2) {
        cout << "Running inferences...\n";
        vector< vector<int> > start_seqs;
        GetStartSeqs (p,gvarbin,start_seqs);
        
        //Set up object to store parameters
        //We want a likelihood, set of parameters, start vector
        vector<double> varbin_init;
        for (int st=0;st<start_seqs.size();st++){  //Loop over potential start points
            if (p.verb==1) {
                cout << "Start sequence " << st << "\n";
            }
            vector<double> varbin_init;
            GetOriginalSeq (p,st,start_seqs,varbin_init);

            //Make temporary samples with this point added
            gvarbin=gvarbin_orig;
            gtimes=gtimes_orig;
            vector< vector<int> > nremoved;
            vector< vector<int> > gfixpos;
            vector< vector<int> > gqfixpos;
            vector< vector< vector<sample> > > gseq_data_array;
            //Series of routines to compile an array containing all of the data neded for the optimisation
            GetGSeqDataArray (p,varbin_init,gvarbin,gtimes,gfixpos,gqfixpos,nremoved,gseq_data_array);

            //Calculate likelihoods for each combination of end-point uncertainty.
            CalculateBestModels (p,st,Lpen,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,rgen);
        }
            
        int max_index=-1;
        maxL=FindBestModel(max_index,outputs);
        OutputBestModelParameters (max_index,outputs);
        if (maxL>best_log) {
            cout << "Better " << maxL << " " << best_log << " at " << i << "\n";
            best_log=maxL;
            best_set=i;
        }
    } else {
        cout << "Set excluded\n";
    }
    gvarbin=gvarbin_orig;
    gtimes=gtimes_orig;
}

double FindPenalty (const int i, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets) {
    double Lpen=0;
    vector<int> overlap; //Number of repeats of a variant across populations + 1
    for (int l=0;l<varbin[sets[i][0][0]].size();l++) {
        overlap.push_back(0);
    }
    vector<int> inset; //Samples that have been considered
    for (int l=0;l<times.size();l++) {
        inset.push_back(1);
    }
    double Lfrac=log(varbin[sets[i][0][0]].size()/29782.);  //Rough probability of a repeated variant
    
    //Find the number of times variants have been found
    for (int j=0;j<sets[i].size();j++) { //Subpopulations
        vector<int> done = overlap;
        for (int k=0;k<done.size();k++) {
            done[k]=0;
        }
        for (int k=0;k<sets[i][j].size();k++) {//Individuals in subpopulation
            inset[sets[i][j][k]]=0;
            for (int l=0;l<varbin[sets[i][j][k]].size();l++) { //Variant positions
                if (varbin[sets[i][j][k]][l]==1&&done[l]==0) {
                    overlap[l]++;
                    done[l]=1;
                }
            }
        }
    }
    
    //Repeat for the individuals not in the main subpopulation.
    vector<int> done;
    for (int j=0;j<varbin[0].size();j++) {
        done.push_back(0);
    }
    for (int j=0;j<inset.size();j++) {
        if (inset[j]==1) {
            for (int l=0;l<varbin[j].size();l++) { //Variant positions
                if (varbin[j][l]==1&&done[l]==0) {
                    overlap[l]++;
                    done[l]=1;
                }
            }
        }
    }
        
    for (int j=0;j<overlap.size();j++) {
        if (overlap[j]>1) {
            Lpen=Lpen+Lfrac;
        }
    }
    return Lpen;
}


void CalculateBestModels (run_params& p, int st, const double Lpen, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen) {
    vector<double> model_parameters;
    for (int i=0;i<gseq_data_array.size();i++) {
        modelstore m;
        m.start_seq=start_seqs[st];
        if (p.verb==1) {
            cout << "i= " << i << " " << gseq_data_array.size() << "\n";
        }
        model_parameters.clear();
        //Check model
        int check=1;
        for (int j=0;j<gseq_data_array[i].size();j++) {
            for (int k=0;k<gseq_data_array[i][j].size();k++){
                if (gseq_data_array[i][j][k].dt==0&&gseq_data_array[i][j][k].nfix>0) {
                    check=0;
                    for (int l=0;l<p.sets;l++) {
                        m.rates.push_back(-1);
                    }
                    m.error=-1;
                    m.lL=-1e9;
                    m.index=i;
                    outputs.push_back(m);
                    break;
                }
            }
        }
        //cout << "Check " << check << "\n";
        if (check!=0) { //Don't have a gain of variants in zero time.
            //cout << "Index " << i << " Problem with seq_data\n";
        //} else {
            OptimiseMultiRateModel (p,p.sets,Lpen,gseq_data_array[i],model_parameters,rgen);
            for (int j=0;j<p.sets;j++) {
                m.rates.push_back(model_parameters[j]);
            }
            m.error=model_parameters[p.sets];
            m.lL=model_parameters[p.sets+1];
            m.index=i;
            outputs.push_back(m);
            //Print out details of the individual model - likelihood and which variants fix or are provisional
            //cout << "Go to indexing\n";
            Indexing (p,i,m,model_parameters,gfixpos,gqfixpos,nremoved);
        }
    }
}

double FindBestModel (int& max_index, const vector<modelstore>& outputs) {
    double maxL=-1e10;
    for (int j=0;j<outputs.size();j++) {
        if (outputs[j].lL>maxL) {
            maxL=outputs[j].lL;
            max_index=j;
        }
    }
    return maxL;
}

void RefineModels (double maxL, vector<modelstore>& outputs) {
    vector<modelstore> outputs_best;
    for (int j=0;j<outputs.size();j++) {
        if (outputs[j].lL>maxL-2) {
            outputs_best.push_back(outputs[j]);
        }
    }
    outputs=outputs_best;
}

void GetBestModelParameters (double maxL, const vector<modelstore>& outputs, vector<double>& model_parameters_best) {
    for (int j=0;j<outputs.size();j++) {
        if (outputs[j].lL==maxL) {
            for (int k=0;k<outputs[j].rates.size();k++) {
                model_parameters_best.push_back(outputs[j].rates[k]);
            }
            model_parameters_best.push_back(outputs[j].error);
            model_parameters_best.push_back(maxL);
            break;
        }
    }
}

void InitialiseLimits (run_params& p, const vector<double>& model_parameters_best, vector< vector<double> >& limits) {
    vector<double> l;
    for (int i=0;i<=p.sets;i++) {
        l.push_back(model_parameters_best[i]);
        l.push_back(model_parameters_best[i]);
        limits.push_back(l);
        l.clear();
    }
}

void ModelRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen) {
    cout << "Model Rate Extremes\n";
    vector<double> extreme_model_parameters;
    for (int j=0;j<outputs.size();j++) {
        //cout << "Check output " << j << "\n";
        p.verb=0;
        vector<double> varbin_init;
        GetOriginalSeq2 (p,outputs[j].start_seq,varbin_init);
        //Set up data to generate uncertainty
        
        
        vector< vector< vector<double> > > gvarbin=gvarbin_orig;
        vector< vector<int> > gtimes=gtimes_orig;
        vector< vector<int> > nremoved;
        vector< vector<int> > gfixpos;
        vector< vector<int> > gqfixpos;
        vector< vector< vector<sample> > > gseq_data_array;
        //Series of routines to compile an array containing all of the data neded for the optimisation
        GetGSeqDataArray (p,varbin_init,gvarbin,gtimes,gfixpos,gqfixpos,nremoved,gseq_data_array);


        vector<double> initial_model_parameters;
        for (int k=0;k<outputs[j].rates.size();k++) {
            initial_model_parameters.push_back(outputs[j].rates[k]);
        }
        initial_model_parameters.push_back(outputs[j].error);
        initial_model_parameters.push_back(maxL);
        vector<double> model_parameters=initial_model_parameters;
        
        for (int m=0;m<=p.sets;m++) {
            for (int n=0;n<=1;n++) {
                //cout << "M " << m << " N " << n << "\n";
                UncertaintyMultiRateModel (p,p.sets,m,n,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
                /*cout << "Limits\n";
                for (int ii=0;ii<limits.size();ii++) {
                    for (int jj=0;jj<limits[ii].size();jj++) {
                        cout << limits[ii][jj] << " ";
                    }
                    cout << "\n";
                }*/
                if (n==1&&extreme_model_parameters[m]>limits[m][1-n]) {
                    limits[m][1-n]=extreme_model_parameters[m];
                }
                if (n==0&&extreme_model_parameters[m]<limits[m][1-n]) {
                    limits[m][1-n]=extreme_model_parameters[m];
                }
            }
        }
        
        if (p.verb==1) {
            cout << "New limits\n";
            for (int k=0;k<limits.size();k++) {
                for (int l=0;l<limits[k].size();l++) {
                    cout << limits[k][l] << " ";
                }
                cout << "\n";
            }
        }
    }
}

void GetOriginalSeq2 (run_params& p, const vector<int>& start_seqs, vector<double>& varbin_init) {
    for (int j=0;j<start_seqs.size();j++) {
        varbin_init.push_back(start_seqs[j]);
    }
    if (p.verb==1) {
        cout << "Original\n";
        for (int j=0;j<varbin_init.size();j++) {
            cout << varbin_init[j] << " ";
        }
        cout << "\n";
    }
}


double FindMaxL (vector<modelstore>& outputs) {
    double maxL=-1e9;
    for (int j=0;j<outputs.size();j++) {
        if (outputs[j].lL>maxL) {
            maxL=outputs[j].lL;
        }
    }
    return maxL;
}


void ModelVRateExploration (run_params& p, int i, const double maxL, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, vector<modelstore>& outputs, vector< vector<double> >& acceptable) {
    cout << "Rate exploration\n";
    
    cout << "Set " << i << "\n";

    //Store starting data in a vector of sets
    vector< vector<int> > gtimes;
    vector< vector< vector<double> > > gvarbin;
    AssignDataToSetsG (p,i,sets,times,varbin,gtimes,gvarbin);
    vector< vector< vector<double> > > gvarbin_orig=gvarbin;
    vector< vector<int> > gtimes_orig=gtimes;
    WriteGVarbin(gtimes,gvarbin);
        
    //Check times compliance - max set size >= 3
    int max=CheckTimes(times,gtimes);
    if (max>2) {
        cout << "Running inferences...\n";
        vector< vector<int> > start_seqs;
        GetStartSeqs (p,gvarbin,start_seqs);
        vector<double> varbin_init;
        for (int st=0;st<start_seqs.size();st++){  //Loop over potential start points
            vector<double> varbin_init;
            GetOriginalSeq (p,st,start_seqs,varbin_init);
            //Make temporary samples with this point added
            gvarbin=gvarbin_orig;
            gtimes=gtimes_orig;
            vector< vector<int> > nremoved;
            vector< vector<int> > gfixpos;
            vector< vector<int> > gqfixpos;
            vector< vector< vector<sample> > > gseq_data_array;
            //Series of routines to compile an array containing all of the data neded for the optimisation
            GetGSeqDataArray (p,varbin_init,gvarbin,gtimes,gfixpos,gqfixpos,nremoved,gseq_data_array);

            CalculateStateSpace(p,2,st,maxL,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,acceptable);
        }
            
    } else {
        cout << "Set excluded\n";
    }
}







void AssignDataToSetsG (run_params& p, int i, const vector< vector< vector<int> > >& sets, const vector<int>& times, const vector< vector<double> >& varbin, vector< vector<int> >& gtimes, vector< vector< vector<double> > >& gvarbin) {
    vector<int> pushed;
    //Generalise the following to multiple sets
    for (int j=0;j<p.sets-1;j++) {  //Set j of p.sets sets
        vector<int> t;
        vector< vector<double> > v;
        for (int k=0;k<sets[i][j].size();k++) {
            t.push_back(times[sets[i][j][k]]);
            v.push_back(varbin[sets[i][j][k]]);
            pushed.push_back(sets[i][j][k]);
        }
        gtimes.push_back(t);
        gvarbin.push_back(v);
    }
    sort(pushed.begin(),pushed.end());
    int index=0;  //Final set p.sets of p.sets
    vector<int> t;
    vector< vector<double> > v;
    for (int j=0;j<varbin.size();j++) {
        if (index<pushed.size()&&pushed[index]==j) {
            index++;
        } else {
            t.push_back(times[j]);
            v.push_back(varbin[j]);
        }
    }
    gtimes.push_back(t);
    gvarbin.push_back(v);
    if (p.verb==1) {
        cout << "Times\n";
        for (int j=0;j<gtimes.size();j++) {
            for (int k=0;k<gtimes[j].size();k++) {
                cout << gtimes[j][k] << " ";
            }
            cout << "\n";
        }
    }
}

















