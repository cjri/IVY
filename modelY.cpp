#include "seqmodel.h"
#include "io.h"
#include "modelI.h"
#include "modelV.h"
#include "modelY.h"
#include "statespace.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void RunModelY (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen) {
    cout << "Model Y\n";
    p.sets=3;
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
    CompileSetsY (p,clusters,sets);
    
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
            vector<modelstore> outputs;
            ExploreSetY (p,i,best_set,best_log,times,outputs,varbin,sets,rgen);
        }
        cout << "Best set " << best_set << " likelihood " << best_log << "\n";
    } else {
        double best_log=-1e9;
        cout << "Explore specified set\n";
        vector<modelstore> outputs;
        ExploreSetY (p,p.chosen_set,p.chosen_set,best_log,times,outputs,varbin,sets,rgen);
        if (p.uncertainty==2) {
            double maxL=-1e9;
            for (int j=0;j<outputs.size();j++) {
                if (outputs[j].lL>maxL) {
                    maxL=outputs[j].lL;
                }
            }
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

void ExploreSetY (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, vector<modelstore>& outputs, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen) {
    //Go through all sets and find the one producing the maximum likelihood reconstruction.
    //We don't need to do the uncertainty calculation in this case.
    
    cout << "Set " << i << "\n";
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
        GetStartSeqsY (p,gvarbin,start_seqs);
        
        
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
            cout << "Gfixtimes\n";
            for (int k=0;k<gfixtimes.size();k++) {
                for (int l=0;l<gfixtimes[k].size();l++){
                    cout << gfixtimes[k][l] << " ";
                }
                cout << "\n";
            }

            //Find positions of fixation events
            vector< vector<int> > gqfixpos;
            vector<int> gnq;
            MakeGQFixPos (p,gvarbin,gfixpos,gfixtimes,gqfixpos,gnq);
                
            //Generate sequence data
            vector< vector<sample> > gseq_data;
            MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
            cout << "Here GS Size " << gseq_data.size() << "\n";

            //Make array over error in the final time point.
            vector< vector<int> > nremoved;  //Count of how many fixations are removed from the final timepoint
            vector< vector< vector<sample> > > gseq_data_array;
            MakeGSDataArrayY (p,gnq,gseq_data,nremoved,gseq_data_array);
               
            cout << "Size of array " << gseq_data_array.size() << "\n";

            //Calculate likelihoods for each combination of end-point uncertainty.
            //Need to embed all of the above in a loop over the zero time point.
            //Need to work out what is the scope for the first unobserved time point - substitutions between first-observed samples?  Does this work for three-way split?
            CalculateBestModelsY(p,st,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,rgen);
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
        InitialiseLimitsY (model_parameters_best,limits);
        //Find limits
        ModelYRateExtremes (p,maxL,gvarbin_orig,gtimes_orig,outputs,limits,rgen);
        
        WriteLimits(limits);
    }
    
  

}


void UnwindObject (vector<objx> new_input_unique, vector< vector<double> >& new_input) {
    for (int l=0;l<new_input_unique.size();l++) {
        for (int m=0;m<new_input_unique[l].y.size();m++) {
            for (int n=0;n<new_input_unique[l].y[m].z.size();n++) {
                vector<double> temp;
                temp.push_back(new_input_unique[l].x);
                temp.push_back(new_input_unique[l].y[m].y);
                temp.push_back(new_input_unique[l].y[m].z[n]);
                new_input.push_back(temp);
            }
        }
    }
}


void CompileSetsY (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    CalculateSetsSystematic3 (clusters,sets);
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

void GetStartSeqsY (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs) {
    cout << "GetStartSeqsY\n";
    start_seqs.clear();
    vector<int> start;
    for (int k=0;k<gvarbin[0][0].size();k++) {
        start.push_back(0);
    }
    start_seqs.push_back(start);
    
    cout << "Number of start seqs " << start_seqs.size() << "\n";
}

void MakeGSDataArrayY (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array) {
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
       // cout << "Interpretation of i value\n";
       // cout << "Index Set1_Fix- Set2_Fix- Set3_Fix-\n";
    }
    //int index=0;
    for (int i0=0;i0<gseq_data_temp[0].size();i0++) {
        for (int i1=0;i1<gseq_data_temp[1].size();i1++) {
            for (int i2=0;i2<gseq_data_temp[2].size();i2++) {
                vector<int> nr;
                nr.push_back(i0);
                nr.push_back(i1);
                nr.push_back(i2);
                nremoved.push_back(nr);

                /*if (p.verb==1) {
                    cout << index << " " << i0 << " " << i1 << " " << i2 << "\n";
                }*/
                vector< vector<sample> > seq_data_array;
                seq_data_array.push_back(gseq_data_temp[0][i0]);
                seq_data_array.push_back(gseq_data_temp[1][i1]);
                seq_data_array.push_back(gseq_data_temp[2][i2]);
                gseq_data_array.push_back(seq_data_array);
                //index++;
            }
        }
    }
}


void CalculateBestModelsY (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen) {
    vector<double> model_parameters;
    //int index_best=0;
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
                    m.rates.push_back(-1);
                    m.rates.push_back(-1);
                    m.rates.push_back(-1);
                    m.error=-1;
                    m.lL=-1e9;
                    m.index=i;
                    outputs.push_back(m);
                    break;
                }
            }
        }
        if (check!=0) { //Don't have a gain of variants in zero time.
            //cout << "Index " << i << " Problem with seq_data\n";
        //} else {
            OptimiseMultiRateModel (p,3,gseq_data_array[i],model_parameters,rgen);
            m.rates.push_back(model_parameters[0]);
            m.rates.push_back(model_parameters[1]);
            m.rates.push_back(model_parameters[2]);
            m.error=model_parameters[3];
            m.lL=model_parameters[4];
            m.index=i;
            outputs.push_back(m);
            cout << "Index " << i << " ";
            for (int p=0;p<=2;p++) {
                cout << "Population " << p << " Rate " << model_parameters[p] << " Fixes " << gfixpos[p].size()-gqfixpos[p].size() << " ";
                for (int k=0;k<gfixpos[p].size();k++) {
                    int match=0;
                    for (int l=0;l<gqfixpos[p].size();l++) {
                        if (gqfixpos[p][l]==gfixpos[p][k]) {
                            match=1;
                        }
                    }
                    if (match==0) {
                        cout << gfixpos[p][k] << " ";
                    }
                }
                cout << "Provisional " << gqfixpos[p].size() << " ";
                for (int l=0;l<gqfixpos[p].size();l++) {
                    cout << gqfixpos[p][l] << " ";
                }
                cout << "Removed " << nremoved[i][p] << " ";
            }
            cout << "Likelihood " << m.lL << "\n";
        }
    }
}

void InitialiseLimitsY (const vector<double>& model_parameters_best, vector< vector<double> >& limits) {
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
}
    
    

void ModelYRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen) {
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
        vector< vector<int> > nremoved;
        vector< vector< vector<sample> > > gseq_data_array;
        MakeGSDataArrayY (p,gnq,gseq_data,nremoved,gseq_data_array);

        vector<double> initial_model_parameters;
        for (int k=0;k<outputs[j].rates.size();k++) {
            initial_model_parameters.push_back(outputs[j].rates[k]);
        }
        initial_model_parameters.push_back(outputs[j].error);
        initial_model_parameters.push_back(maxL);
        vector<double> model_parameters=initial_model_parameters;
        
        UncertaintyMultiRateModel (p,3,0,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]>limits[0][0]) {
            limits[0][0]=extreme_model_parameters[0];
        }
        UncertaintyMultiRateModel (p,3,0,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]<limits[0][1]) {
            limits[0][1]=extreme_model_parameters[0];
        }
        UncertaintyMultiRateModel (p,3,1,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]>limits[1][0]) {
            limits[1][0]=extreme_model_parameters[1];
        }
        UncertaintyMultiRateModel (p,3,1,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]<limits[1][1]) {
            limits[1][1]=extreme_model_parameters[1];
        }
        UncertaintyMultiRateModel (p,3,2,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[2]>limits[2][0]) {
            limits[2][0]=extreme_model_parameters[2];
        }
        UncertaintyMultiRateModel (p,3,2,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[2]<limits[2][1]) {
            limits[2][1]=extreme_model_parameters[2];
        }
        UncertaintyMultiRateModel (p,3,3,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[3]>limits[3][0]) {
            limits[3][0]=extreme_model_parameters[3];
        }
        UncertaintyMultiRateModel (p,3,3,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[3]<limits[3][1]) {
            limits[3][1]=extreme_model_parameters[3];
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









