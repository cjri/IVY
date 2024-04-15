#include "seqmodel.h"
#include "io.h"
#include "modelI.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void FindConsensus (string& consensus, vector<string>& seqs) {
    //Here we are finding the initial nucleotide
    consensus=seqs[0];
    if (seqs.size()<2) {
        cout << "Error: Need at least two sequences\n";
    }
    for (int seq=1;seq<seqs.size();seq++) {
        for (int pos=0;pos<seqs[0].size();pos++) {
            if (consensus[pos]!='A'&&consensus[pos]!='C'&&consensus[pos]!='G'&&consensus[pos]!='T') {
                consensus[pos]=seqs[seq][pos];
            }
        }
    }
    for (int pos=0;pos<seqs[0].size();pos++) {
        if (consensus[pos]!='A'&&consensus[pos]!='C'&&consensus[pos]!='G'&&consensus[pos]!='T') {
            consensus[pos]='-';
        }
    }
    ofstream cons_file;
    cons_file.open("Consensus.fa");
    cons_file << ">Consensus_sequence\n";
    cons_file << consensus << "\n";

}

void FindVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //cout << "Find variants\n";
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<seqs[i].size();pos++) {
            if (seqs[i].compare(pos,1,consensus,pos,1)!=0) {
                if (seqs[i].compare(pos,1,"A")==0||seqs[i].compare(pos,1,"C")==0||seqs[i].compare(pos,1,"G")==0||seqs[i].compare(pos,1,"T")==0) {
                    //cout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << seqs[i][pos] << "\n";
                    s.locus.push_back(pos);
                    s.allele.push_back(seqs[i][pos]);
                }
            }
        }
        variants.push_back(s);
    }
}

void FindAmbiguousVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //Variants which are not A, C, G, or T
    vector<int> vpos;
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            vpos.push_back(variants[i].locus[j]);
        }
    }
    sort(vpos.begin(),vpos.end());
    vpos.erase(unique(vpos.begin(),vpos.end()),vpos.end());
    FindVariants2 (vpos,variants,consensus,seqs);
}

void FindVariants2 (vector<int>& vpos, vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //cout << "Find variants\n";
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<vpos.size();pos++) {
            if (seqs[i].compare(vpos[pos],1,consensus,vpos[pos],1)!=0) {
                int found=0;
                for (int j=0;j<variants[i].locus.size();j++) {
                    if (variants[i].locus[j]==vpos[pos]) {
                        found=1;
                    }
                }
                if (found==0) {
                    variants[i].locus.push_back(vpos[pos]);
                    variants[i].allele.push_back(seqs[i][vpos[pos]]);
                }
            }
        }
    }
}


void FixVariants (vector<string>& seqs, vector<sparseseq>& variants) {
    //Fix for R
    vector<int> rpos;
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='R') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'G','A','R',seqs,variants);
    //Fix for K
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='K') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'T','G','K',seqs,variants);

    //Fix for Y
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='Y') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'C','T','Y',seqs,variants);
    
    //Fix for W
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='W') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'A','T','W',seqs,variants);
}

void RevertToN (vector<int>& xpos, char a1, char a2, char ax, vector<string>& seqs, vector<sparseseq>& variants) {
    //cout << "Revert to N\n";
    for (int i=0;i<xpos.size();i++) {
        int n1=0;
        int n2=0;
        for (int j=0;j<seqs.size();j++) {
            //cout << "Here " << j << " " << xpos[i] << " " << seqs[j][xpos[i]] << "\n";
            if (seqs[j][xpos[i]]==a1) {
                n1++;
            } else if (seqs[j][xpos[i]]==a2) {
                n2++;
            }
        }
        //cout << "Count " << n1 << " " << n2 << "\n";
        if (n1>0&&n2>0) {
            //Go through variants and find N values
            for (int j=0;j<variants.size();j++) {
                for (int k=0;k<variants[j].locus.size();k++) {
                    if (variants[j].allele[k]==ax) {
                        variants[j].allele[k]='N';
                    }
                }
            }
        }
    }
}

void FindVariantPositions (run_params& p, vector<sparseseq>& variants, vector<int>& varpos) {
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            varpos.push_back(variants[i].locus[j]);
        }
    }
    sort(varpos.begin(),varpos.end());
    varpos.erase(unique(varpos.begin(),varpos.end()),varpos.end());
}

void MakeVarmatVarbin (run_params& p, const vector<string>& seqs, const vector<int>& varpos, vector< vector<char> >& varmat, vector< vector<double> >& varbin) {
    //Sequence at each variant position
    for (int i=0;i<seqs.size();i++) {
        vector<char> v;
        for (int j=0;j<varpos.size();j++) {
            v.push_back(seqs[i][varpos[j]]);
        }
        varmat.push_back(v);
    }
    //Convert to binary code, where 0 is initial consensus
    for (int i=0;i<varmat.size();i++) {
        vector<double> v;
        for (int j=0;j<varmat[i].size();j++) {
            if (varmat[i][j]==varmat[0][j]) {
                v.push_back(0);
            } else {
                //Find what variant this is
                if (varmat[i][j]=='A'||varmat[i][j]=='C'||varmat[i][j]=='G'||varmat[i][j]=='T') {
                    v.push_back(1);
                } else {
                    //Do basic calculation.  Number of cases that are ACGT but not 0 as a fraction of ACGT
                    double same=0;
                    double diff=0;
                    for (int k=0;k<varmat.size();k++) {
                        if (varmat[k][j]=='A'||varmat[k][j]=='C'||varmat[k][j]=='G'||varmat[k][j]=='T') {
                            if (varmat[k][j]==varmat[0][j]) {
                                same++;
                            } else {
                                diff++;
                            }
                        }
                    }
                    double r=diff/(same+diff);
                    //If it is one half, introduce a bias of \epsilon towards the consensus sequence
                    if (r==0.5) {
                        r=r-1e-10;
                    }
                    v.push_back(r);
                }
            }
        }
        varbin.push_back(v);
    }
    if (p.verb==1) {
        cout << "Variant matrix\n";
        for (int i=0;i<varmat.size();i++) {
            for (int j=0;j<varmat[i].size();j++) {
                cout << varmat[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "Variant binary\n";
        for (int i=0;i<varbin.size();i++) {
            for (int j=0;j<varbin[i].size();j++) {
                cout << varbin[i][j] << " ";
            }
            cout << "\n";
        }

    }
}

void MakeVarmatVarbin2 (run_params& p, const vector<string>& seqs, const vector<int>& varpos, vector<char>& consensus, vector< vector<char> >& varmat, vector< vector<double> >& varbin) {
    //New version of this algorithm.  Now we are not using uncertainties other than just as N we don't care what the fractions are any more
    //Added in correction for when the first variant could be an N
    //Sequence at each variant position
    for (int i=0;i<seqs.size();i++) {
        vector<char> v;
        for (int j=0;j<varpos.size();j++) {
            v.push_back(seqs[i][varpos[j]]);
        }
        varmat.push_back(v);
    }
    for (int j=0;j<varmat[0].size();j++) {
        consensus.push_back(varmat[0][j]);
    }
    for (int j=0;j<consensus.size();j++) {
        if (consensus[j]=='N') {
            for (int i=0;i<varmat.size();i++) {
                if (varmat[i][j]!='N') {
                    consensus[j]=varmat[i][j];
                    break;
                }
            }
        }
    }
    cout << "Consensus\n";
    for (int j=0;j<consensus.size();j++) {
        cout << consensus[j] << " ";
    }
    cout << "\n";
    
    //Convert to binary code, where 0 is initial consensus
    for (int i=0;i<varmat.size();i++) {
        vector<double> v;
        for (int j=0;j<varmat[i].size();j++) {
            if (varmat[i][j]==consensus[j]) {
                v.push_back(0);
            } else {
                //Find what variant this is
                if (varmat[i][j]=='A'||varmat[i][j]=='C'||varmat[i][j]=='G'||varmat[i][j]=='T') {
                    v.push_back(1);
                } else {
                    v.push_back(0.5);
                }
            }
        }
        varbin.push_back(v);
    }
    if (p.verb==1) {
        cout << "Variant matrix\n";
        for (int i=0;i<varmat.size();i++) {
            for (int j=0;j<varmat[i].size();j++) {
                cout << varmat[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "Variant binary\n";
        for (int i=0;i<varbin.size();i++) {
            for (int j=0;j<varbin[i].size();j++) {
                cout << varbin[i][j] << " ";
            }
            cout << "\n";
        }

    }
}



void MakeConstantFixes (run_params& p, const vector< vector<double> >& varbin, vector< vector<double> >& constant, vector< vector<double> >& fixes) {
    //Vectors to test for fit for fixation of flucutation at a site.
    for (int i=0;i<varbin.size();i++) {
        vector<double> c;
        for (int j=0;j<2;j++) {
            c.push_back(j);
        }
        constant.push_back(c);
    }

    for (int i=0;i<varbin.size();i++) {
        vector<double> c;
        for (int j=0;j<varbin.size()-1;j++) {
            if (j>=i) {
                c.push_back(0);
            } else {
                c.push_back(1);
            }
        }
        fixes.push_back(c);
    }
    /*if (p.verb==1) {
        cout << "Constant\n";
        for (int i=0;i<constant.size();i++) {
            for (int j=0;j<constant[i].size();j++) {
                cout << constant[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "Fixed\n";
        for (int i=0;i<fixes.size();i++) {
            for (int j=0;j<fixes[i].size();j++) {
                cout << fixes[i][j] << " ";
            }
            cout << "\n";
        }
    }*/
}

void FindFixationSites (run_params& p, const vector< vector<double> >& varbin, const vector< vector<double> >& constant, const vector< vector<double> >& fixes, vector<int>& fixpos, vector<int>& flucpos) {
    //Compare seuqence trajectories to constant type
    vector<int> min_comp;
    vector<int> zeros;
    for (int i=0;i<varbin[0].size();i++) {
        int z=0;
        double min_c=1000;
        for (int j=0;j<constant[0].size();j++) {
            double c=0;
            for (int k=0;k<varbin.size();k++) {
                if (varbin[k][i]!=constant[k][j]) {
                    c=c+abs(varbin[k][i]-constant[k][j]);
                }
            }
            if (c<min_c) {
                min_c=c;
            }
            if (c==0&&j==0) {
                z=1;
            }

        }
        min_comp.push_back(min_c);
        zeros.push_back(z);
    }
    //Compare sequence trajectories to fixation types
    vector<int> min_fix;
    for (int i=0;i<varbin[0].size();i++) {
        double min_c=1000;
        for (int j=0;j<fixes[0].size();j++) {
            double c=0;
            for (int k=0;k<varbin.size();k++) {
                if (varbin[k][i]!=fixes[k][j]) {
                    c=c+abs(varbin[k][i]-fixes[k][j]);
                }
            }
            if (c<min_c) {
                min_c=c;
                //cout << "New minimum " << min_c << "\n";
            }
        }
        min_fix.push_back(min_c);
    }

    //Identify fixation sites
    for (int i=0;i<min_comp.size();i++) {
        if (min_fix[i]<min_comp[i]) {
            fixpos.push_back(i);
        } else {
            if (zeros[i]==0) {
                flucpos.push_back(i);
            }
        }
    }
    
    if (p.verb==1) {
        cout << "Fixations\n";
        for (int i=0;i<fixpos.size();i++) {
            cout << fixpos[i] << " ";
        }
        cout << "\n";
        cout << "Fluctuations\n";
        for (int i=0;i<flucpos.size();i++) {
            cout << flucpos[i] << " ";
        }
        cout << "\n";
    }

}

void FindFixationTimes (run_params& p, const vector<int>& fixpos, const vector< vector<double> >& varbin, vector<int>& fixtimes) {
    for (int i=0;i<fixpos.size();i++) {
        for (int j=0;j<varbin.size();j++) {
            if (varbin[j][fixpos[i]]==1) {
                fixtimes.push_back(j);
                break;
            }
        }
    }
    if (p.verb==1) {
        cout << "Fixation times\n";
        for (int i=0;i<fixtimes.size();i++) {
            cout << fixtimes[i] << " ";
        }
        cout << "\n";
    }
}

void FindQFixPositions (run_params& p, const vector<int>& fixpos, const vector<int>& fixtimes, const vector< vector<double> >& varbin, vector<int>& qfixpos) {
    //Sites which fix in the final time-point
    for (int i=0;i<fixtimes.size();i++) {
        if (fixtimes[i]==varbin.size()-1) {
            qfixpos.push_back(fixpos[i]);
        }
    }
    if (p.verb==1) {
        cout << "Qfix positions\n";
        for (int i=0;i<qfixpos.size();i++) {
            cout << qfixpos[i] << " ";
        }
        cout << "\n";
    }
}

void SetupSequenceData (run_params& p, const vector<int>& times, const vector<int>& fixtimes, const vector<int>& fixpos, const vector<int>& flucpos, const vector< vector<double> >& varbin, vector<sample>& seq_data) {
    //Set up an initial sample.  We don't know anything about it.  If we have multiple clusters we might want to alter the basic parameters:  Added fixations before the first sample could detract distance from a weird sample, implying fixations in the main population prior to the first sample.
  //  cout << "SetupSequenceData\n";
  //  cout << varbin.size() << "\n";
    sample s;
    s.dt=times[0];
    s.nfix=-1;
    s.nfluc=-1;
    seq_data.push_back(s);
    if (p.verb==1) {
        cout << "Sequence data\n";
        cout << "0 " << s.dt << " " << s.nfix << " " << s.nfluc << "\n";
    }
    for (int i=1;i<times.size();i++) {
        sample s;
        //Difference in times
        s.dt=times[i]-times[i-1];
        //Number of fixation events.  N.B. we define a fixation time by the first observation of a non-ambiguous nucleotide.  Implies that all ambiguous nucleotides in fixation sites are fluctuations
        s.nfix=0;
        for (int j=0;j<fixtimes.size();j++) {
            if (fixtimes[j]==i) {
                s.nfix++;
            }
        }
        //Number of fluctuations
        s.nfluc=0;
        vector<double> xf;
        for (int k=0;k<flucpos.size();k++) {  //Is there a fluctuation at this site?
            if (varbin[i][flucpos[k]]==1) {
                s.nfluc++;
            }
            if (varbin[i][flucpos[k]]>0&&varbin[i][flucpos[k]]<1) {
                xf.push_back(varbin[i][flucpos[k]]);
            }
        }
        for (int k=0;k<fixpos.size();k++) {  //Is there a reverse fluctuation at a fixation site?
            if (fixtimes[k]<i) {
                if (varbin[i][fixpos[k]]==0) {
                    s.nfluc++;
                }
                if (varbin[i][fixpos[k]]>0&&varbin[i][fixpos[k]]<1) { //Reverse fluctuations
                    xf.push_back(1.-varbin[i][fixpos[k]]);
                }
            }
            if (fixtimes[k]>=i) {
                if (varbin[i][fixpos[k]]>0&&varbin[i][fixpos[k]]<1) {
                    xf.push_back(varbin[i][fixpos[k]]);  //Forward fluctuations previous to sequencing
                }
            }
        }
        /*cout << "Time point " << i << "\n";
        for (int j=0;j<xf.size();j++) {
            cout << xf[j] << " ";
        }
        cout << "\n";*/
        
        //Make a vector of binary vectors of size 2^n where n is the size of xf
        if (xf.size()>0) {
            CalculateSxFluc(p,xf,s);
        }
        //Calculate the probabilites of each outcome - number of stochastic fluctuations
        
        if (p.verb==1) {
            cout << i << " " << s.dt << " " << s.nfix << " " << s.nfluc << " " << s.xfluc.size() << " ";
            for (int j=0;j<s.xfluc.size();j++) {
                cout << s.xfluc[j] << " ";
            }
            cout << "\n";
        }
        seq_data.push_back(s);
    }
}

void CalculateSxFluc (run_params& p, const vector<double>& xf, sample& s) {
    vector< vector<int> > binset;
    MakeBinaries(xf.size(),binset);
    vector<double> binprobs;
    vector<int> bincounts;
    for (int i=0;i<binset.size();i++) {
        double p=1;
        int c=0;
        for (int j=0;j<binset[i].size();j++) {
            if (binset[i][j]==1) {
                p=p*xf[j];
                c++;
            } else {
                p=p*(1-xf[j]);
            }
        }
        binprobs.push_back(p);
        bincounts.push_back(c);
    }
    
    for (int j=0;j<xf.size()+1;j++) {
        s.xfluc.push_back(0);
    }
    
    for (int j=0;j<binprobs.size();j++) {
        s.xfluc[bincounts[j]]=s.xfluc[bincounts[j]]+binprobs[j];
    }
}

void MakeBinaries (int n_size, vector< vector<int> >& binset) {
    //Intervals of 1, 2, 4, 8, etc
    int vsize=pow(2,n_size);
    vector<int> v;
    for (int i=0;i<n_size;i++) {
        v.push_back(0);
    }
    binset.push_back(v);
    for (int i=1;i<vsize;i++) {
        v[0]=1-v[0];
        if (i%2==0) {
            v[1]=1-v[1];
        }
        if (i%4==0) {
            v[2]=1-v[2];
        }
        if (i%8==0) {
            v[3]=1-v[3];
        }
        if (i%16==0) {
            v[4]=1-v[4];
        }
        if (i%32==0) {
            v[5]=1-v[5];
        }
        if (i%64==0) {
            v[6]=1-v[6];
        }
        if (i%128==0) {
            v[7]=1-v[7];
        }
        binset.push_back(v);
    }
}

void ConstructSequenceDataArray(int nq, const vector<sample>& seq_data, vector<int>& removed, vector< vector<sample> >& seq_data_array) {
    //Array includes all uncertainty about fluctuations and fixations in the final time point
    for (int i=0;i<=nq;i++) {
        vector<sample> sd=seq_data;
        sd[sd.size()-1].nfix=sd[sd.size()-1].nfix-i;
        sd[sd.size()-1].nfluc=sd[sd.size()-1].nfluc+i;
        removed.push_back(i);
        seq_data_array.push_back(sd);
    }
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
            cout << "Provisional " << qfixpos.size();
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


void InitialiseLimits (const vector<double>& model_parameters_best, vector< vector<double> >& limits) {
    vector<double> l;
    l.push_back(model_parameters_best[0]);
    l.push_back(model_parameters_best[0]);
    limits.push_back(l);
    l.clear();
    l.push_back(model_parameters_best[1]);
    l.push_back(model_parameters_best[1]);
    limits.push_back(l);
}

/*void AssignDataToSets2 (const vector< vector<int> >& sets, const vector<int>& times, const vector< vector<double> >& varbin, vector<int>& times1, vector<int>& times2, vector< vector<double> >& varbin1, vector< vector<double> >& varbin2) {
    vector<int> pushed;
    //Generalise the following to multiple sets
    for (int j=0;j<sets[i].size();j++) {
        times1.push_back(times[sets[i][j]]);
        varbin1.push_back(varbin[sets[i][j]]);
        pushed.push_back(sets[i][j]);
    }
    int index=0;
    for (int j=0;j<varbin.size();j++) {
        if (index<pushed.size()&&pushed[index]==j) {
            index++;
        } else {
            times2.push_back(times[j]);
            varbin2.push_back(varbin[j]);
        }
    }
}*/




int CheckTimes (const vector<int>& times, vector< vector<int> >& gtimes) {
    int max=0;
    for (int i=0;i<gtimes.size();i++) {
        if (gtimes[i].size()>max) {
            max=gtimes[i].size();
        }
    }
    for (int i=0;i<gtimes.size();i++) {
        if (gtimes[i].size()==1&&gtimes[i][0]==times[0]) {
            max=0;
        }
    }
    return max;
}




void CalculateSetsSystematic2 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/2;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector< vector<int> > sets1;
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
                sets.push_back(sets1);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector< vector<int> > sets1;
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                    sets.push_back(sets1);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector< vector<int> > sets1;
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                        sets.push_back(sets1);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector< vector<int> > sets1;
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                            sets.push_back(sets1);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector< vector<int> > sets1;
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                                sets.push_back(sets1);
                            }
                        }
                    }
                }
            }
        }
        if (n==6) {
            for (int i1=0;i1<clusters.size()-5;i1++) {
                for (int i2=i1+1;i2<clusters.size()-4;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-3;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-2;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-1;i5++) {
                                for (int i6=i5+1;i6<clusters.size();i6++) {
                                    vector< vector<int> > sets1;
                                    vector<int> s;
                                    s.push_back(i1);
                                    s.push_back(i2);
                                    s.push_back(i3);
                                    s.push_back(i4);
                                    s.push_back(i5);
                                    s.push_back(i6);
                                    sets1.push_back(s);
                                    sets.push_back(sets1);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==7) {
            for (int i1=0;i1<clusters.size()-6;i1++) {
                for (int i2=i1+1;i2<clusters.size()-5;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-4;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-3;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-2;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-1;i6++) {
                                    for (int i7=i6+1;i7<clusters.size();i7++) {
                                        vector< vector<int> > sets1;
                                        vector<int> s;
                                        s.push_back(i1);
                                        s.push_back(i2);
                                        s.push_back(i3);
                                        s.push_back(i4);
                                        s.push_back(i5);
                                        s.push_back(i6);
                                        s.push_back(i7);
                                        sets1.push_back(s);
                                        sets.push_back(sets1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==8) {
            for (int i1=0;i1<clusters.size()-7;i1++) {
                for (int i2=i1+1;i2<clusters.size()-6;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-5;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-4;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-3;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-2;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-1;i7++) {
                                        for (int i8=i7+1;i8<clusters.size();i8++) {
                                            vector< vector<int> > sets1;
                                            vector<int> s;
                                            s.push_back(i1);
                                            s.push_back(i2);
                                            s.push_back(i3);
                                            s.push_back(i4);
                                            s.push_back(i5);
                                            s.push_back(i6);
                                            s.push_back(i7);
                                            s.push_back(i8);
                                            sets1.push_back(s);
                                            sets.push_back(sets1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==9) {
            for (int i1=0;i1<clusters.size()-8;i1++) {
                for (int i2=i1+1;i2<clusters.size()-7;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-6;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-5;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-4;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-3;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-2;i7++) {
                                        for (int i8=i7+1;i8<clusters.size()-1;i8++) {
                                            for (int i9=i8+1;i9<clusters.size();i9++) {
                                                vector< vector<int> > sets1;
                                                vector<int> s;
                                                s.push_back(i1);
                                                s.push_back(i2);
                                                s.push_back(i3);
                                                s.push_back(i4);
                                                s.push_back(i5);
                                                s.push_back(i6);
                                                s.push_back(i7);
                                                s.push_back(i8);
                                                s.push_back(i9);
                                                sets1.push_back(s);
                                                sets.push_back(sets1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==10) {
            for (int i1=0;i1<clusters.size()-9;i1++) {
                for (int i2=i1+1;i2<clusters.size()-8;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-7;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-6;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-5;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-4;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-3;i7++) {
                                        for (int i8=i7+1;i8<clusters.size()-2;i8++) {
                                            for (int i9=i8+1;i9<clusters.size()-1;i9++) {
                                                for (int ia=i9+1;ia<clusters.size();ia++) {
                                                    vector< vector<int> > sets1;
                                                    vector<int> s;
                                                    s.push_back(i1);
                                                    s.push_back(i2);
                                                    s.push_back(i3);
                                                    s.push_back(i4);
                                                    s.push_back(i5);
                                                    s.push_back(i6);
                                                    s.push_back(i7);
                                                    s.push_back(i8);
                                                    s.push_back(i9);
                                                    s.push_back(ia);
                                                    sets1.push_back(s);
                                                    sets.push_back(sets1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void CalculateSetsSystematic3 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/3;
    vector< vector<int> > sets1;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                            }
                        }
                    }
                }
            }
        }
        if (n==6) {
            for (int i1=0;i1<clusters.size()-5;i1++) {
                for (int i2=i1+1;i2<clusters.size()-4;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-3;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-2;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-1;i5++) {
                                for (int i6=i5+1;i6<clusters.size();i6++) {
                                    vector<int> s;
                                    s.push_back(i1);
                                    s.push_back(i2);
                                    s.push_back(i3);
                                    s.push_back(i4);
                                    s.push_back(i5);
                                    s.push_back(i6);
                                    sets1.push_back(s);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "First set " << sets1.size() << " possibilities\n";
    for (int i=0;i<sets1.size();i++) {
        int min=sets1[i].size();
        int r=clusters.size()-min;
        int max=r/2;
        //Construct remainder
        vector<int> remainder;
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==sets1[i][index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2;
                    s2.push_back(sets1[i]);
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2;
                        s2.push_back(sets1[i]);
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-2;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-1;i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2;
                            s2.push_back(sets1[i]);
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-3;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-2;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-1;i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2;
                                s2.push_back(sets1[i]);
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-4;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-3;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-2;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-1;i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2;
                                    s2.push_back(sets1[i]);
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-5;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-4;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-3;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-2;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-1;i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2;
                                        s2.push_back(sets1[i]);
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-6;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-5;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-4;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-3;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-2;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-1;i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2;
                                            s2.push_back(sets1[i]);
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==8) {
                for (int i1=0;i1<remainder.size()-7;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-6;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-5;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-4;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-3;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-2;i6++) {
                                        for (int i7=i6+1;i7<remainder.size()-1;i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                vector< vector<int> > s2;
                                                s2.push_back(sets1[i]);
                                                vector<int> s;
                                                s.push_back(remainder[i1]);
                                                s.push_back(remainder[i2]);
                                                s.push_back(remainder[i3]);
                                                s.push_back(remainder[i4]);
                                                s.push_back(remainder[i5]);
                                                s.push_back(remainder[i6]);
                                                s.push_back(remainder[i7]);
                                                s.push_back(remainder[i8]);
                                                s2.push_back(s);
                                                sets.push_back(s2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==9) {
                for (int i1=0;i1<remainder.size()-8;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-7;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-6;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-5;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-4;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-3;i6++) {
                                        for (int i7=i6+1;i7<remainder.size()-2;i7++) {
                                            for (int i8=i7+1;i8<remainder.size()-1;i8++) {
                                                for (int i9=i8+1;i9<remainder.size();i9++) {
                                                    vector< vector<int> > s2;
                                                    s2.push_back(sets1[i]);
                                                    vector<int> s;
                                                    s.push_back(remainder[i1]);
                                                    s.push_back(remainder[i2]);
                                                    s.push_back(remainder[i3]);
                                                    s.push_back(remainder[i4]);
                                                    s.push_back(remainder[i5]);
                                                    s.push_back(remainder[i6]);
                                                    s.push_back(remainder[i7]);
                                                    s.push_back(remainder[i8]);
                                                    s.push_back(remainder[i9]);
                                                    s2.push_back(s);
                                                    sets.push_back(s2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void CalculateSetsSystematic4 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/4; //Smallest
    vector< vector<int> > sets1;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "First set " << sets1.size() << " possibilities\n";
    for (int i=0;i<sets1.size();i++) {
        int min=sets1[i].size(); //Set 1 is the smallest
        int r=clusters.size()-min;
        int max=r/3;
        //Construct remainder - clusters not in a set
        vector<int> remainder;
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==sets1[i][index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2;
                    s2.push_back(sets1[i]);
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2;
                        s2.push_back(sets1[i]);
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-2;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-1;i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2;
                            s2.push_back(sets1[i]);
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-3;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-2;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-1;i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2;
                                s2.push_back(sets1[i]);
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-4;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-3;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-2;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-1;i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2;
                                    s2.push_back(sets1[i]);
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-5;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-4;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-3;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-2;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-1;i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2;
                                        s2.push_back(sets1[i]);
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-6;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-5;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-4;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-3;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-2;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-1;i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2;
                                            s2.push_back(sets1[i]);
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    vector< vector< vector<int> > > sets_new;
    //For each of the sets consider the next division
    for (int i=0;i<sets.size();i++) {
        int min=sets[i][1].size(); //Set 2 is the next smallest
        int r=clusters.size()-(sets[i][0].size()+sets[i][1].size());
        int max=r/2;
        vector<int> remainder;
        vector<int> used=sets[i][0];
        for (int j=0;j<sets[i][1].size();j++) {
            used.push_back(sets[i][1][j]);
        }
        sort(used.begin(),used.end());
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==used[index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2=sets[i];
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets_new.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2=sets[i];
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets_new.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2=sets[i];
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets_new.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2=sets[i];
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets_new.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2=sets[i];
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets_new.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2=sets[i];
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets_new.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2=sets[i];
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets_new.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==8) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                vector< vector<int> > s2=sets[i];
                                                vector<int> s;
                                                s.push_back(remainder[i1]);
                                                s.push_back(remainder[i2]);
                                                s.push_back(remainder[i3]);
                                                s.push_back(remainder[i4]);
                                                s.push_back(remainder[i5]);
                                                s.push_back(remainder[i6]);
                                                s.push_back(remainder[i7]);
                                                s.push_back(remainder[i8]);
                                                s2.push_back(s);
                                                sets_new.push_back(s2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==9) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                for (int i9=i8+1;i9<remainder.size();i9++) {
                                                    vector< vector<int> > s2=sets[i];
                                                    vector<int> s;
                                                    s.push_back(remainder[i1]);
                                                    s.push_back(remainder[i2]);
                                                    s.push_back(remainder[i3]);
                                                    s.push_back(remainder[i4]);
                                                    s.push_back(remainder[i5]);
                                                    s.push_back(remainder[i6]);
                                                    s.push_back(remainder[i7]);
                                                    s.push_back(remainder[i8]);
                                                    s.push_back(remainder[i9]);
                                                    s2.push_back(s);
                                                    sets_new.push_back(s2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    sets=sets_new;
}
