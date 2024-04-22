#include "ivy.h"
#include "io.h"
#include "model_tools.h"
#include <iostream>
#include <string>
#include <cstring>

void GetGSeqDataArray (run_params& p, vector<double> varbin_init, vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& gtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gqfixpos, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array) {
    //Overarching code which gets the data all in order
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
        vector< vector<int> > gflucpos;
    MakeGFixFluc (p,gvarbin,gconstant,gfixes,gfixpos,gflucpos);
                            
    //Find fixation times
    vector< vector<int> > gfixtimes;
    MakeGFixTimes (p,gvarbin,gfixpos,gfixtimes);
    
    //Find positions of fixation events
    vector<int> gnq;
    MakeGQFixPos (p,gvarbin,gfixpos,gfixtimes,gqfixpos,gnq);

    //Output fixation information here
    if (p.verb==1) {
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
    }
    
    //Generate sequence data
    vector< vector<sample> > gseq_data;
    MakeGSeqData (p,gvarbin,gtimes,gfixtimes,gfixpos,gflucpos,gseq_data);
        
    //Make array over error in the final time point.
    
    //Now construct object
    MakeGSDataArray (p,gnq,gseq_data,nremoved,gseq_data_array);
    
}



void GetStartSeqs (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs) {
    start_seqs.clear();
    if (p.sets==2) {  //Start sequence is quite complicated to find
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
        
        for (int i=0;i<numbers.size();i++) { //Loop over the number involved if allbin[i][j] is a 1
            for (int j=0;j<allbin.size();j++) {
                start=start_orig;
                for (int k=0;k<allbin[j].size();k++) {
                    if (allbin[j][k]==1) {
                        for (int l=0;l<numbers[i][k];l++) {
                            start[ident[k][l]]=1;
                        }
                    }
                }
                start_seqs.push_back(start);
            }
        }
    }
    if (p.sets>=3) {  //Easier to find start sequences
        vector<int> start;
        for (int k=0;k<gvarbin[0][0].size();k++) {
            start.push_back(0);
        }
        start_seqs.push_back(start);
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


void MakeGFixTimes (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, vector< vector<int> >& gfixtimes) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<int> fixtimes;
        FindFixationTimes (p,gfixpos[j],gvarbin[j],fixtimes);
        gfixtimes.push_back(fixtimes);
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


void MakeGSeqData (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gtimes, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos, vector< vector<sample> >& gseq_data) {
    for (int j=0;j<gvarbin.size();j++) {
        vector<sample> seq_data;
        SetupSequenceData (p,gtimes[j],gfixtimes[j],gfixpos[j],gflucpos[j],gvarbin[j],seq_data);
        gseq_data.push_back(seq_data);
    }
}

void SetupSequenceData (run_params& p, const vector<int>& times, const vector<int>& fixtimes, const vector<int>& fixpos, const vector<int>& flucpos, const vector< vector<double> >& varbin, vector<sample>& seq_data) {
    //Set up an initial sample.  We don't know anything about it.  If we have multiple clusters we might want to alter the basic parameters:  Added fixations before the first sample could detract distance from a weird sample, implying fixations in the main population prior to the first sample.
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



void MakeGSDataArray (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array) {
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
    if (p.sets==2) {
        for (int i0=0;i0<gseq_data_temp[0].size();i0++) {
            for (int i1=0;i1<gseq_data_temp[1].size();i1++) {
                vector<int> nr;
                nr.push_back(i0);
                nr.push_back(i1);
                nremoved.push_back(nr);
                /*if (p.verb==1) {
                 cout << index << " " << i0 << " " << i1 << "\n";  NB Can initialise index variable and increase by one each time
                 }*/
                vector< vector<sample> > seq_data_array;
                seq_data_array.push_back(gseq_data_temp[0][i0]);
                seq_data_array.push_back(gseq_data_temp[1][i1]);
                gseq_data_array.push_back(seq_data_array);
            }
        }
    }
    if (p.sets==3) {
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
                }
            }
        }
    }
    if (p.sets==4) {
        for (int i0=0;i0<gseq_data_temp[0].size();i0++) {
            for (int i1=0;i1<gseq_data_temp[1].size();i1++) {
                for (int i2=0;i2<gseq_data_temp[2].size();i2++) {
                    for (int i3=0;i3<gseq_data_temp[3].size();i3++) {
                        vector<int> nr;
                        nr.push_back(i0);
                        nr.push_back(i1);
                        nr.push_back(i2);
                        nr.push_back(i3);
                        nremoved.push_back(nr);
                        /*if (p.verb==1) {
                            cout << index << " " << i0 << " " << i1 << " " << i2 << " " << i3 << "\n";
                        }*/
                        vector< vector<sample> > seq_data_array;
                        seq_data_array.push_back(gseq_data_temp[0][i0]);
                        seq_data_array.push_back(gseq_data_temp[1][i1]);
                        seq_data_array.push_back(gseq_data_temp[2][i2]);
                        seq_data_array.push_back(gseq_data_temp[3][i3]);
                        gseq_data_array.push_back(seq_data_array);
                    }
                }
            }
        }
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

