#include "seqmodel.h"
#include "io.h"
#include "modelI.h"
#include "modelV.h"
#include "modelX.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void RunModelX (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen) {
    cout << "Model X\n";
    p.sets=4;
    cout << "Varbin\n";
    for (int i=0;i<varbin.size();i++) {
        for (int j=0;j<varbin[i].size();j++) {
            cout << varbin[i][j] << " ";
        }
        cout << "\n";
    }
    //Find samples that share k variants - put into clusters.  In the below k=2: Allows for an error position
    vector< vector<int> > clusters;
    FindClusters(p,2,varbin,clusters);

    //Note here: We find all of the sets made up of clusters, keeping clustered data points together.
    vector< vector< vector<int> > > sets;
    CompileSetsX (p,clusters,sets);
    
    cout << "Sets " << sets.size() << "\n";
    
    //How are we going to find all of the sets?  Could do this in a sytematic way?
    
    //General thought: Adding in substitutions before the first observed sample will bring the substitution rates closer to one another.  Can try this anyway...
    //In effect we are optimising the sequence at time zero
    
    //Split the samples given a particular division
    
    //Can't have just the first one as a distinct set
    
    //Could do systematic search of sets or focus on a specific set
    if (p.systematic==1) {
        int best_set=-1;
        double best_log=-1e9;
        for (int i=0;i<sets.size();i++) { //Systematic search through all possible sets.  Want option to explore one set out of all of them, with uncertainty
            cout << "Run systematic search of sets\n";
            ExploreSetX (p,i,best_set,best_log,times,varbin,sets,rgen);
        }
        cout << "Best set " << best_set << " likelihood " << best_log << "\n";
    } else {
        double best_log=-1e9;
        cout << "Explore specified set\n";
        ExploreSetX (p,p.chosen_set,p.chosen_set,best_log,times,varbin,sets,rgen);
    }

}

void ExploreSetX (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen) {
    //Go through all sets and find the one producing the maximum likelihood reconstruction.
    //We don't need to do the uncertainty calculation in this case.
    
    cout << "Set " << i << "\n";
    double maxL=-1e6;
    vector<modelstore> outputs;

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
            
        //Problem here: Need to propose a constant start point to calculate fixations.
        //Definitely needed if we have only one sample in one of the clusters
        //Default is to use the first time point
        //How do we change this?
            
        //Strategy for changing the nominal first time point: Given the default we will see fixations in the
        //population that doesn't contain the first sample.  We could propose that the reverse fixation to this
        //occurred in the first population before the time of the first sample.
            
        //For the moment we will propose a first sample identical to the earliest observed sample.
        //Note that we need to fix this in the single-population measure: Currently not accounting for the
        //first sample in the likelihood at all.  Uncertainty in the single-population sample - add fixes to the first
        //sample i.e. more than 0.
            
        //Find a set of original vectors.  Can be done in the two-set case by finding the first sequence of each of the two sets
            
        //Differences between first sequences - any combination of these could be the initial sequence
        vector< vector<int> > start_seqs;
        GetStartSeqsX (p,gvarbin,start_seqs);
        
        
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
            if (p.verb==1) {
                cout << "Here GVarbin\n";
                WriteGVarbin(gtimes,gvarbin);
            }
            AddStartToGVarbin (varbin_init,gvarbin);
                
            //Edit GVarbin so that the first time point is zero - revert according to the initial point
            RevertGVarbin(gvarbin);
            
            //Make temporary times vectors with new point at time zero
            AddStartToGTimes (gtimes);
                
            if (p.verb==1) {
                cout << "Now GVarbin\n";

                WriteGVarbin(gtimes,gvarbin);
            }

            //Set up fixation and constant vectors.  Use general form
            vector< vector< vector<double> > > gconstant;
            vector< vector< vector<double> > > gfixes;
            MakeGConstantFix (p,gvarbin,gconstant,gfixes);
                
            //Test variant data against constant and fix vectors
            vector< vector<int> > gfixpos;
            vector< vector<int> > gflucpos;
            MakeGFixFluc (p,gvarbin,gconstant,gfixes,gfixpos,gflucpos);
                                    
            //Find fixation times
            vector< vector<int> > gfixtimes;
            MakeGFixTimes (p,gvarbin,gfixpos,gfixtimes);
                
            //Find positions of fixation events
            vector< vector<int> > gqfixpos;
            vector<int> gnq;
            MakeGQFixPos (p,gvarbin,gfixpos,gfixtimes,gqfixpos,gnq);
                
            //Generate sequence data
            vector< vector<sample> > gseq_data;
            MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
           
            //Make array over error in the final time point.
            vector< vector< vector<sample> > > gseq_data_array;
            MakeGSDataArrayX (p,gnq,gseq_data,gseq_data_array);
               

            //Calculate likelihoods for each combination of end-point uncertainty.
            //Need to embed all of the above in a loop over the zero time point.
            //Need to work out what is the scope for the first unobserved time point - substitutions between first-observed samples?  Does this work for three-way split?
            CalculateBestModelsX (p,st,start_seqs,gseq_data_array,outputs,rgen);
            
        }
            
        int max_index=-1;
        maxL=FindBestModel(max_index,outputs);
        OutputBestModelParametersV (max_index,outputs);
        if (maxL>best_log) {
            cout << "Better " << maxL << " " << best_log << " at " << i << "\n";
            best_log=maxL;
            best_set=i;
        }
    } else {
        cout << "Set excluded\n";
    }
    
    //Find uncertainty in parameters
    if (p.uncertainty==1) {
        RefineModels(maxL,outputs);
        /*cout << "Refined models " << outputs.size() << "\n";
        for (int j=0;j<outputs.size();j++) {
            cout << outputs[j].rates[0] << " " << outputs[j].rates[1] << "\n";
        }*/
        
        //Set up best model parameters
        vector<double> model_parameters_best;
        GetBestModelParameters (maxL,outputs,model_parameters_best);

        //Calculate uncertainty in each of the models
        vector< vector<double> > limits;
        InitialiseLimitsX (model_parameters_best,limits);
        //Find limits
        ModelXRateExtremes (p,maxL,gvarbin_orig,gtimes_orig,outputs,limits,rgen);
        
        WriteLimits(limits);
    }

    
    //In theory could be any combination of the variants between the two.
    //Need a routine to find all of these.
    
    
    

    //In the optimisation here, we can have more than one rate of evolution, but only one rate of error
    //Key is that a faster rate of evolution allows for a lower rate of error in the other population
    
    

    //NB We have an array of possibilities that will potential grow as n^k.
    //Might need to think about how we do the optimisation on this grid.
    
    
    //Thought: Could just not allow ambiguous data to contribute to the likelihood at all.

}



void CompileSetsX (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    CalculateSetsSystematic4 (clusters,sets);
    ConvertSetsClusters (clusters,sets);
    if (p.verb==1) {
        cout << "Sets size " << sets.size() << "\n";
        for (int i=0;i<sets.size();i++) {
            cout << i << " ";
            for (int j=0;j<sets[i].size();j++) {
                cout << "(";
                for (int k=0;k<sets[i][j].size();k++) {
                    cout << sets[i][j][k] << " ";
                }
                cout << ") ";
            }
            cout << "\n";
        }
    }
}



void GetStartSeqsX (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs) {
    cout << "GetStartSeqsX\n";
    start_seqs.clear();
    vector<int> start;
    for (int k=0;k<gvarbin[0][0].size();k++) {
        start.push_back(0);
    }
    start_seqs.push_back(start);
    
    cout << "Number of start seqs " << start_seqs.size() << "\n";
}

void MakeGSDataArrayX (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector< vector<sample> > >& gseq_data_array) {
    vector< vector< vector<sample> > > gseq_data_temp;
    for (int i=0;i<gseq_data.size();i++) {
        vector< vector<sample> > seq_data_array;
        vector<int> temp;
        ConstructSequenceDataArray(gnq[i],gseq_data[i],temp,seq_data_array);
        gseq_data_temp.push_back(seq_data_array);
        if (p.verb==1) {
            cout << "Size of array " << seq_data_array.size() << "\n";
        }
    }

    //All combinations of uncertainty in last time point
    if (p.verb==1) {
        cout << "Interpretation of i value\n";
        cout << "Index Set1_Fix- Set2_Fix- Set3_Fix- Set4_Fix-\n";
    }
    int index=0;
    for (int i0=0;i0<gseq_data_temp[0].size();i0++) {
        for (int i1=0;i1<gseq_data_temp[1].size();i1++) {
            for (int i2=0;i2<gseq_data_temp[2].size();i2++) {
                for (int i3=0;i3<gseq_data_temp[3].size();i3++) {
                    if (p.verb==1) {
                        cout << index << " " << i0 << " " << i1 << " " << i2 << " " << i3 << "\n";
                    }
                    vector< vector<sample> > seq_data_array;
                    seq_data_array.push_back(gseq_data_temp[0][i0]);
                    seq_data_array.push_back(gseq_data_temp[1][i1]);
                    seq_data_array.push_back(gseq_data_temp[2][i2]);
                    seq_data_array.push_back(gseq_data_temp[3][i3]);
                    gseq_data_array.push_back(seq_data_array);
                    index++;
                }
            }
        }
    }
}


void CalculateBestModelsX (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen) {
    vector<double> model_parameters;
    int index_best=0;
    for (int i=0;i<gseq_data_array.size();i++) {
        modelstore m;
        m.start_seq=start_seqs[st];
        if (p.verb==1) {
            cout << "i= " << i << "\n";
        }
        model_parameters.clear();
        //Check model
        int check=1;
        for (int j=0;j<gseq_data_array[i].size();j++) {
            for (int k=0;k<gseq_data_array[i][j].size();k++){
                if (gseq_data_array[i][j][k].dt==0&&gseq_data_array[i][j][k].nfix>0) {
                    check=0;
                    break;
                }
            }
        }
        if (check==0) {
            cout << "Problem with seq_data\n";
        } else {
            OptimiseMultiRateModel (p,4,gseq_data_array[i],model_parameters,rgen);
            m.rates.push_back(model_parameters[0]);
            m.rates.push_back(model_parameters[1]);
            m.rates.push_back(model_parameters[2]);
            m.rates.push_back(model_parameters[3]);
            m.error=model_parameters[4];
            m.lL=model_parameters[5];
            m.index=i;
            outputs.push_back(m);
        }
    }
}

void InitialiseLimitsX (const vector<double>& model_parameters_best, vector< vector<double> >& limits) {
    vector<double> l;
    l.push_back(model_parameters_best[0]);
    l.push_back(model_parameters_best[0]);
    limits.push_back(l);
    l.clear();
    l.push_back(model_parameters_best[1]);
    l.push_back(model_parameters_best[1]);
    limits.push_back(l);
    l.clear();
    l.push_back(model_parameters_best[2]);
    l.push_back(model_parameters_best[2]);
    limits.push_back(l);
    l.clear();
    l.push_back(model_parameters_best[3]);
    l.push_back(model_parameters_best[3]);
    limits.push_back(l);
    l.clear();
    l.push_back(model_parameters_best[4]);
    l.push_back(model_parameters_best[4]);
    limits.push_back(l);
}

void ModelXRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen) {
    vector<double> extreme_model_parameters;
    for (int j=0;j<outputs.size();j++) {
        //cout << "Check output " << j << "\n";
        p.verb=0;
        vector<double> varbin_init;
        GetOriginalSeq2 (p,outputs[j].start_seq,varbin_init);
        //Set up data to generate uncertainty
        vector< vector< vector<double> > > gvarbin=gvarbin_orig;
        AddStartToGVarbin (varbin_init,gvarbin);
        RevertGVarbin(gvarbin);
        vector< vector<int> > gtimes=gtimes_orig;
        AddStartToGTimes (gtimes);
        vector< vector< vector<double> > > gconstant;
        vector< vector< vector<double> > > gfixes;
        MakeGConstantFix (p,gvarbin,gconstant,gfixes);
        vector< vector<int> > gfixpos;
        vector< vector<int> > gflucpos;
        MakeGFixFluc (p,gvarbin,gconstant,gfixes,gfixpos,gflucpos);
        vector< vector<int> > gfixtimes;
        MakeGFixTimes (p,gvarbin,gfixpos,gfixtimes);
        vector< vector<int> > gqfixpos;
        vector<int> gnq;
        MakeGQFixPos (p,gvarbin,gfixpos,gfixtimes,gqfixpos,gnq);
        vector< vector<sample> > gseq_data;
        MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
        vector< vector< vector<sample> > > gseq_data_array;
        MakeGSDataArrayX (p,gnq,gseq_data,gseq_data_array);

        vector<double> initial_model_parameters;
        for (int k=0;k<outputs[j].rates.size();k++) {
            initial_model_parameters.push_back(outputs[j].rates[k]);
        }
        initial_model_parameters.push_back(outputs[j].error);
        initial_model_parameters.push_back(maxL);
        vector<double> model_parameters=initial_model_parameters;
        
        for (int i=0;i<4;i++) {
            UncertaintyMultiRateModel (p,4,i,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
            if (extreme_model_parameters[i]>limits[i][0]) {
                limits[i][0]=extreme_model_parameters[i];
            }
            UncertaintyMultiRateModel (p,4,i,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
            if (extreme_model_parameters[i]<limits[i][1]) {
                limits[i][1]=extreme_model_parameters[i];
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

