#include "seqmodel.h"
#include "io.h"
#include "modelI.h"
#include "modelV.h"
#include "statespace.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void RunModelV (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen) {
    p.sets=2;
    cout << "Varbin\n";
    for (int i=0;i<varbin.size();i++) {
        for (int j=0;j<varbin[i].size();j++) {
            cout << varbin[i][j] << " ";
        }
        cout << "\n";
    }
    //Find samples that share k variants - put into clusters.  In the below k=2: Allows for a single error position
    vector< vector<int> > clusters;
    FindClusters(p,2,varbin,clusters);

    //Note here: We find all of the sets made up of clusters, keeping clustered data points together.
    vector< vector< vector<int> > > sets;
    CompileSets (p,clusters,sets);
    
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
            ExploreSet (p,i,best_set,best_log,times,outputs,varbin,sets,rgen);
        }
        cout << "Best set " << best_set << " likelihood " << best_log << "\n";
    } else {
        double best_log=-1e9;
        cout << "Explore specified set\n";
        vector<modelstore> outputs;
        ExploreSet (p,p.chosen_set,p.chosen_set,best_log,times,outputs,varbin,sets,rgen);
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

void ExploreSet (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, vector<modelstore>& outputs, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen) {
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
        GetStartSeqsV (p,gvarbin,start_seqs);
        
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

            //Output fixation information here
            cout << "Fixation details\n";
            for (int i=0;i<gfixpos.size();i++) {
                cout << "Population " << i << "\n";
                for (int j=0;j<gfixpos[i].size();j++) {
                    cout << "(Pos " << gfixpos[i][j] << ", Time " << gfixtimes[i][j];
                    int match=0;
                    for (int k=0;k<gqfixpos[i].size();k++) {
                        if (gqfixpos[i][k]==gfixpos[i][j]) {
                            match=1;
                        }
                    }
                    if (match==1) {
                        cout << "*)\n";
                    } else {
                        cout << ")\n";
                    }
                }
            }
            
            
            //Generate sequence data
            vector< vector<sample> > gseq_data;
            MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
                
            //Make array over error in the final time point.
            vector< vector<int> > nremoved;  //Count of how many fixations are removed from the final timepoint
            vector< vector< vector<sample> > > gseq_data_array;
            MakeGSDataArrayV (p,gnq,gseq_data,nremoved,gseq_data_array);
               

            //Calculate likelihoods for each combination of end-point uncertainty.
            //Need to embed all of the above in a loop over the zero time point.
            //Need to work out what is the scope for the first unobserved time point - substitutions between first-observed samples?  Does this work for three-way split?
            CalculateBestModelsV (p,st,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,rgen);
            
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
        InitialiseLimitsV (model_parameters_best,limits);
        //Find limits
        ModelVRateExtremes (p,maxL,gvarbin_orig,gtimes_orig,outputs,limits,rgen);
        
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
        GetStartSeqsV (p,gvarbin,start_seqs);
        vector<double> varbin_init;
        for (int st=0;st<start_seqs.size();st++){  //Loop over potential start points
            vector<double> varbin_init;
            GetOriginalSeq (p,st,start_seqs,varbin_init);
            //Make temporary samples with this point added
            gvarbin=gvarbin_orig;
            gtimes=gtimes_orig;
            AddStartToGVarbin (varbin_init,gvarbin);
            RevertGVarbin(gvarbin);
            AddStartToGTimes (gtimes);
            vector< vector< vector<double> > > gconstant;
            vector< vector< vector<double> > > gfixes;
            MakeGConstantFix (p,gvarbin,gconstant,gfixes);
            //Test variant data against constant and fix vectors
            vector< vector<int> > gfixpos;
            vector< vector<int> > gflucpos;
            MakeGFixFluc (p,gvarbin,gconstant,gfixes,gfixpos,gflucpos);
            vector< vector<int> > gfixtimes;
            MakeGFixTimes (p,gvarbin,gfixpos,gfixtimes);
            //Find positions of fixation events
            vector< vector<int> > gqfixpos;
            vector<int> gnq;
            MakeGQFixPos (p,gvarbin,gfixpos,gfixtimes,gqfixpos,gnq);
            vector< vector<sample> > gseq_data;
            MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
            vector< vector<int> > nremoved;  //Count of how many fixations are removed from the final timepoint
            vector< vector< vector<sample> > > gseq_data_array;
            MakeGSDataArrayV (p,gnq,gseq_data,nremoved,gseq_data_array);

            //Change the code here: Use data from previous outputs and do the calculation in discrete space.  NB Need to import the outputs.


            CalculateStateSpace(p,2,st,maxL,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,acceptable);
        }
            
    } else {
        cout << "Set excluded\n";
    }
}



void FindClusters(run_params& p, int threshold, const vector< vector<double> >& varbin, vector< vector<int> >& clusters) {
    vector<int> accounted;
    for (int i=0;i<varbin.size();i++) {
        accounted.push_back(0);
    }
    for (int i=0;i<varbin.size();i++) {
        vector<int> c;
        if (accounted[i]==0) {
            c.push_back(i);
            accounted[i]=1;
            vector<int> added;
            added.clear();
            for (int j=i+1;j<varbin.size();j++) {  //Problem here: We add j but not things attached to j
                if (accounted[j]==0) {
                    int same=0;
                    for (int k=0;k<varbin[i].size();k++) {
                        if (varbin[i][k]==1&&varbin[j][k]==1) {
                            same++;
                        }
                    }
                    if (same>=threshold) {
                        accounted[j]=1;
                        c.push_back(j);
                        added.push_back(j);

                    }
                }
            }
            while (added.size()>0) {
                vector<int> added_temp;
                //Go through the things that have been added.  Add links to c and links to added.
                for (int a=0;a<added.size();a++) {
                    for (int j=i+1;j<varbin.size();j++) {
                        if (accounted[j]==0) {
                            int same=0;
                            for (int k=0;k<varbin[added[a]].size();k++) {
                                if (varbin[added[a]][k]==1&&varbin[j][k]==1) {
                                    same++;
                                }
                            }
                            if (same>=threshold) {
                                accounted[j]=1;
                                c.push_back(j);
                                added_temp.push_back(j);
                            }
                        }
                    }
                }
                added.clear();
                added=added_temp;
            }
            clusters.push_back(c);
        }
    }
    //Sort clusters - find unique elements
    for (int i=0;i<clusters.size();i++) {
        sort(clusters[i].begin(),clusters[i].end());
        clusters[i].erase(unique(clusters[i].begin(),clusters[i].end()),clusters[i].end());
    }
    
    if (p.verb==1) {
        cout << "Clusters\n";
        for (int i=0;i<clusters.size();i++) {
            for (int j=0;j<clusters[i].size();j++) {
                cout << clusters[i][j] << " ";
            }
            cout << "\n";
            for (int j=0;j<clusters[i].size();j++) {
                for (int k=0;k<varbin[clusters[i][j]].size();k++) {
                    cout << varbin[clusters[i][j]][k] << " ";
                }
                cout << "\n";
            }
        }
    }
}

void CompileSets (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    CalculateSetsSystematic2 (clusters,sets);
    ConvertSetsClusters (clusters,sets);
    if (p.verb==1) {
        cout << "Sets size " << sets.size() << "\n";
        for (int i=0;i<sets.size();i++) {
            cout << i << " ";
            for (int j=0;j<sets[i][0].size();j++) {
                cout << sets[i][0][j] << " ";
            }
            cout << "\n";
        }
    }
}

void ConvertSetsClusters (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    vector< vector< vector<int> > > sets_temp;
    for (int i=0;i<sets.size();i++) {
        vector< vector<int> > ss;
        for (int j=0;j<sets[i].size();j++) {
            vector<int> s;
            for (int k=0;k<sets[i][j].size();k++) {
                for (int l=0;l<clusters[sets[i][j][k]].size();l++) {
                    s.push_back(clusters[sets[i][j][k]][l]);
                }
            }
            sort(s.begin(),s.end());
            ss.push_back(s);
        }
        sets_temp.push_back(ss);
    }
    sets=sets_temp;
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

void GetStartSeqsV (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs) {
    vector<int> start;
    for (int k=0;k<gvarbin[0][0].size();k++) {
        start.push_back(0);
    }
    vector<int> start_orig=start;

    //Find differences between initial sequences - first one in each varbin
    vector<int> diffs;
    for (int k=0;k<gvarbin[0][0].size();k++) {
        if (gvarbin[0][0][k]!=gvarbin[1][0][k]) {
            diffs.push_back(k);
        }
    }

    cout << "Differences " << diffs.size() << "\n";
    for (int i=0;i<diffs.size();i++) {
        cout << diffs[i] << "\n";
    }
    
    //Find identical columns i.e. in the two sets of gvarbin
    vector< vector<int> > ident;
    vector<int> found;
    for (int j=0;j<gvarbin[0][0].size();j++) {
        found.push_back(0);
    }
    
    //Find identical columns in two sets - just need one instance of n of these rather than considering every combination
    int index=0;
    while (index<gvarbin[0][0].size()) {
        if (found[index]==0) {
            vector<int> id;
            //Check if part of difference
            int part=0;
            for (int n=0;n<diffs.size();n++) {
                if (diffs[n]==index) {
                    part=1;
                    break;
                }
            }
            if (part==1) {
                id.push_back(index);
            }
            found[index]=1;
            for (int k=index+1;k<gvarbin[0][0].size();k++) {
                if (found[k]==0) {
                    int same=1;
                    for (int j=0;j<gvarbin[0].size();j++) {
                        if (gvarbin[0][j][index]!=gvarbin[0][j][k]) {
                            //cout << "Difference0 " << index << " " << k << " at " << j << "\n";
                            same=0;
                            break;
                        }
                    }
                    if (same==1) {
                        for (int j=0;j<gvarbin[1].size();j++) {
                            if (gvarbin[1][j][index]!=gvarbin[1][j][k]) {
                                //cout << "Difference1 " << index << " " << k << " at " << j << "\n";
                                same=0;
                                break;

                            }
                        }
                    }
                    if (same==1) {
                        int part=0;
                        for (int n=0;n<diffs.size();n++) {
                            if (diffs[n]==k) {
                                part=1;
                                break;
                            }
                        }
                        if (part==1) {
                            id.push_back(k);
                        }
                        found[k]=1;
                    }
                }
            }
            if (id.size()>0) {
                ident.push_back(id);
            }
        }
        index++;
    }
    
    if (p.verb==1) {
        cout << "Equivalences\n";
        for (int i=0;i<ident.size();i++) {
            for (int j=0;j<ident[i].size();j++) {
                cout << ident[i][j] << " ";
            }
            cout << "\n";
        }
    }
    
    //How the next bit works: Do the binary bit for the number of sets.
    //Within each binary part, have from 1 to n identical variants

    //Binary counting routine - binary strings of length ident.size()
    vector< vector<int> > allbin;
    vector<int> bin;
    for (int i=0;i<ident.size();i++) {
        bin.push_back(0);
    }
    index=0;
    int tot=pow(2,ident.size());
    allbin.push_back(bin);
    for (int j=0;j<tot-1;j++) {
        index=0;
        if (bin[index]==0) {
            bin[index]=1;
        } else {
            while (bin[index]==1) {
                bin[index]=0;
                index++;
            }
            bin[index]=1;
        }
        allbin.push_back(bin);
    }
    /*cout << "Made allbin\n";
    for (int i=0;i<allbin.size();i++) {
        for (int j=0;j<allbin[i].size();j++) {
            cout << allbin[i][j] << " ";
        }
        cout << "\n";
    }*/
    //Generate a set of the number of variants in each set included if there is a 1 at allbin[i][j]
    vector< vector<int> > numbers;
    vector<int> nc;
    tot=1;
    for (int s=0;s<ident.size();s++) {
        nc.push_back(1);
        tot=tot*ident[s].size();
    }
    numbers.push_back(nc);
    index=0;
    for (int j=1;j<tot;j++) {
        index=0;
        if (nc[index]<ident[index].size()) {
            nc[index]++;
        } else {
            while (nc[index]==ident[index].size()) {
                nc[index]=1;
                index++;
            }
            nc[index]++;
        }
        numbers.push_back(nc);
    }
    
    /*cout << "Made numbers\n";
    for (int i=0;i<numbers.size();i++) {
        for (int j=0;j<numbers[i].size();j++) {
            cout << numbers[i][j] << "\n";
        }
    }*/
    
    for (int i=0;i<numbers.size();i++) { //Loop over the number involved if allbin[i][j] is a 1
        for (int j=0;j<allbin.size();j++) {
            //cout << "j= " << j << "\n";
            start=start_orig;
            for (int k=0;k<allbin[j].size();k++) {
                //cout << "k= " << k << "\n";
                if (allbin[j][k]==1) {
                    for (int l=0;l<numbers[i][k];l++) {
                        start[ident[k][l]]=1;
                    }
                }
            }
            /*for (int j=0;j<start.size();j++) {
                cout << start[j] << " ";
            }
            cout << "\n";*/
            start_seqs.push_back(start);
        }
    }
    
    cout << "Number of start seqs " << start_seqs.size() << "\n";
}

void GetOriginalSeq (run_params& p, int st, const vector< vector<int> >& start_seqs, vector<double>& varbin_init) {
    for (int j=0;j<start_seqs[st].size();j++) {
        varbin_init.push_back(start_seqs[st][j]);
    }
    if (p.verb==1) {
        cout << "Original\n";
        for (int j=0;j<varbin_init.size();j++) {
            cout << varbin_init[j] << " ";
        }
        cout << "\n";
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


void AddStartToGVarbin (const vector<double>& varbin_init, vector< vector< vector<double> > >& gvarbin) {
    vector< vector< vector<double> > > gvarbin_temp;
    for (int j=0;j<gvarbin.size();j++) {
        vector< vector<double> > vv;
        vv.push_back(varbin_init);
        for (int k=0;k<gvarbin[j].size();k++) {
            vv.push_back(gvarbin[j][k]);
        }
        gvarbin_temp.push_back(vv);
    }
    gvarbin=gvarbin_temp;
}

void RevertGVarbin(vector< vector< vector<double> > >& gvarbin) {
    for (int j=0;j<gvarbin.size();j++) {
        for (int k=0;k<gvarbin[j][0].size();k++) {
            if (gvarbin[j][0][k]==1) {
                for (int l=0;l<gvarbin[j].size();l++) {
                    gvarbin[j][l][k]=1-gvarbin[j][l][k];
                }
            }
        }
    }
}

void AddStartToGTimes (vector< vector<int> >& gtimes) {
    vector< vector<int> > gtimes_temp;
    for (int j=0;j<gtimes.size();j++) {
        vector<int> t;
        t.push_back(0);
        for (int k=0;k<gtimes[j].size();k++) {
            t.push_back(gtimes[j][k]);
        }
        gtimes_temp.push_back(t);
    }
    gtimes=gtimes_temp;
}

void MakeGConstantFix (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes) {
    for (int j=0;j<gvarbin.size();j++) {
        vector< vector<double> > constant;
        vector< vector<double> > fixes;
        MakeConstantFixes(p,gvarbin[j],constant,fixes);
        gconstant.push_back(constant);
        gfixes.push_back(fixes);
    }
}

void MakeGFixFluc (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<int> fixpos;
        vector<int> flucpos;
        //NB This routine now accounts for positions with no variation: Are neither fixations or fluctuations
        FindFixationSites (p,gvarbin[j],gconstant[j],gfixes[j],fixpos,flucpos);
        gfixpos.push_back(fixpos);
        gflucpos.push_back(flucpos);
    }
}

void MakeGFixTimes (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, vector< vector<int> >& gfixtimes) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<int> fixtimes;
        FindFixationTimes (p,gfixpos[j],gvarbin[j],fixtimes);
        gfixtimes.push_back(fixtimes);
    }
}

void MakeGQFixPos (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gqfixpos, vector<int>& gnq) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<int> qfixpos;
        FindQFixPositions (p,gfixpos[j],gfixtimes[j],gvarbin[j],qfixpos);
        gqfixpos.push_back(qfixpos);
    }
    for (int j=0;j<gqfixpos.size();j++) {
        gnq.push_back(gqfixpos[j].size());
    }
}

void MakeGSeqData (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gtimes, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos, vector< vector<sample> >& gseq_data) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<sample> seq_data;
        SetupSequenceData (p,gtimes[j],gfixtimes[j],gfixpos[j],gflucpos[j],gvarbin[j],seq_data);
        gseq_data.push_back(seq_data);
    }
}

void MakeGSDataArrayV (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array) {
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
    
    //Note: In this model we don't specify which fixations occur in the final time-point, only the number of fixations
    //Identifying which fixations occur is a matter of combinatorics
    if (p.verb==1) {
        //cout << "Interpretation of i value\n";
        //cout << "Index Set1_Fix- Set2_Fix-\n";
    }
    for (int i0=0;i0<gseq_data_temp[0].size();i0++) {
        for (int i1=0;i1<gseq_data_temp[1].size();i1++) {
            vector<int> nr;
            nr.push_back(i0);
            nr.push_back(i1);
            nremoved.push_back(nr);
            /*if (p.verb==1) {
                cout << index << " " << i0 << " " << i1 << "\n";
            }*/
            vector< vector<sample> > seq_data_array;
            seq_data_array.push_back(gseq_data_temp[0][i0]);
            seq_data_array.push_back(gseq_data_temp[1][i1]);
            gseq_data_array.push_back(seq_data_array);
        }
    }
}

void CalculateBestModelsV (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen) {
    vector<double> model_parameters;
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
        cout << "Check " << check << "\n";
        if (check!=0) { //Don't have a gain of variants in zero time.
            //cout << "Index " << i << " Problem with seq_data\n";
        //} else {
            OptimiseMultiRateModel (p,2,gseq_data_array[i],model_parameters,rgen);
            m.rates.push_back(model_parameters[0]);
            m.rates.push_back(model_parameters[1]);
            m.error=model_parameters[2];
            m.lL=model_parameters[3];
            m.index=i;
            outputs.push_back(m);
            //1.  What happened to the start seq?
            //2.  Which fixations occur in each population?
            //We want to print out which variants fix in each population in this line, along with which are provisional
            cout << "Index " << i << " ";
            for (int p=0;p<=1;p++) {
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


void OptimiseMultiRateModel (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters, gsl_rng *rgen) {
    //cout << "Optimise " << gseq_data.size() << "\n";
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
    for (int it=0;it<100000;it++) {
        if (first==0) {
            if (lL>lL_best) {
                rates_best=rates;
                error_best=error;
                lL_best=lL;
               // cout << "Better " << rate << " " << error << " " << lL << "\n";
            } else {
                rates=rates_best;
                error=error_best;
            }
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
        lL=0;
        for (int i=0;i<gseq_data.size();i++) {
            for (int j=0;j<gseq_data[i].size();j++) {
                double r=rates[i]*gseq_data[i][j].dt;
                int done=0;
                if (r==0&&gseq_data[i][j].dt>0) {
                    lL=lL-1e9;
                    done=1;
                }
                if (done==0) {
                    if (r>0||gseq_data[i][j].nfix>0) {
                        lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfix,r));
                    }
                    if (gseq_data[i][j].nfluc>=0||gseq_data[i][j].xfluc.size()>0) {
                        if (gseq_data[i][j].xfluc.size()==0) {
                            lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,error));
                        } else {
                            //Account for uncertainty in the number of flucutations arising from ambiguous nucleotides
                            for (int j=0;j<gseq_data[i][j].xfluc.size();j++) {
                                lL=lL+gseq_data[i][j].xfluc[j]*log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc+j,error));
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i=0;i<n_rates;i++) {
        model_parameters.push_back(rates_best[i]);
    }
    model_parameters.push_back(error_best);
    model_parameters.push_back(lL_best);
    /*if (p.verb==1) {
        cout << "Parameters here\n";
        for (int i=0;i<model_parameters.size();i++) {
            cout << model_parameters[i] << " ";
        }
        cout << "\n";
        cout << "Check zeros\n";
        for (int i=0;i<gseq_data.size();i++) {
            for (int j=0;j<gseq_data[i].size();j++) {
                double r=rates[i]*gseq_data[i][j].dt;
                if (r==0) {
                    cout << i << " " << j << " R=0 nf= " << gseq_data[i][j].nfix << "\n";
                }
            }
        }

    }*/
}

void InitialiseLimitsV (const vector<double>& model_parameters_best, vector< vector<double> >& limits) {
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
}

void ModelVRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen) {
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
        MakeGSDataArrayV (p,gnq,gseq_data,nremoved,gseq_data_array);

        vector<double> initial_model_parameters;
        for (int k=0;k<outputs[j].rates.size();k++) {
            initial_model_parameters.push_back(outputs[j].rates[k]);
        }
        initial_model_parameters.push_back(outputs[j].error);
        initial_model_parameters.push_back(maxL);
        vector<double> model_parameters=initial_model_parameters;
        
        UncertaintyMultiRateModel (p,2,0,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]>limits[0][0]) {
            limits[0][0]=extreme_model_parameters[0];
        }
        UncertaintyMultiRateModel (p,2,0,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[0]<limits[0][1]) {
            limits[0][1]=extreme_model_parameters[0];
        }
        UncertaintyMultiRateModel (p,2,1,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]>limits[1][0]) {
            limits[1][0]=extreme_model_parameters[1];
        }
        UncertaintyMultiRateModel (p,2,1,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[1]<limits[1][1]) {
            limits[1][1]=extreme_model_parameters[1];
        }
        UncertaintyMultiRateModel (p,2,2,1,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[2]>limits[2][0]) {
            limits[2][0]=extreme_model_parameters[2];
        }
        UncertaintyMultiRateModel (p,2,2,0,maxL,gseq_data_array[outputs[j].index],initial_model_parameters,model_parameters,extreme_model_parameters,rgen);
        if (extreme_model_parameters[2]<limits[2][1]) {
            limits[2][1]=extreme_model_parameters[2];
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
                    }
                } else {
                    model_parameters[i]=model_parameters[i]+(gsl_rng_uniform(rgen)*dx)-(dx/2);
                }
            }
        }
        
        //Evaluate the likelihood
        lL=0;
        for (int i=0;i<gseq_data.size();i++) {
            for (int j=0;j<gseq_data[i].size();j++) {
                double r=model_parameters[i]*gseq_data[i][j].dt;
                int done=0;
                if (r==0&&gseq_data[i][j].dt>0) {
                        lL=lL-1e9;
                        done=1;
                }
                if (done==0) {
                    if (r>0||gseq_data[i][j].nfix>0) {
                        lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfix,r));
                    }
                    lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,model_parameters[n_rates]));

                    
                    /*if (gseq_data[i][j].xfluc.size()==0) {
                        lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,model_parameters[n_rates]));
                    } else {
                        //Account for uncertainty in the number of flucutations arising from ambiguous nucleotides
                        for (int j=0;j<gseq_data[i][j].xfluc.size();j++) {
                            lL=lL+gseq_data[i][j].xfluc[j]*log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc+j,model_parameters[n_rates]));
                        }
                    }*/
                }
            }
        }
    }
}
