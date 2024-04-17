#include "ivy.h"
#include "io.h"
#include "model.h"
#include "modelI.h"
#include <iostream>
#include <string>
#include <cstring>

void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
    p.verb=0;
    p.sets=1;
    p.systematic=1;
    p.uncertainty=0;
    p.model="I";
    p.chosen_set=-1;
    p.fix_error=-1; //Fixed error rate.  -1 implies unfixed
    p.times_file="Times.dat";
    p.seqs_file="Sets.in";
    p.fix=0;
    while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--input")==0) {
                x++;
            p.seqs_file=argv[x];
        } else if (p_switch.compare("--model")==0) {
                x++;
                p.model=argv[x];
        } else if (p_switch.compare("--times")==0) {
                x++;
                p.times_file=argv[x];
        } else if (p_switch.compare("--verb")==0) {
                x++;
                p.verb=atoi(argv[x]);
        } else if (p_switch.compare("--sets")==0) {
                x++;
                p.sets=atoi(argv[x]);
        } else if (p_switch.compare("--uncertainty")==0) {
                x++;
                p.uncertainty=atoi(argv[x]);
        } else if (p_switch.compare("--error")==0) {
                x++;
                p.fix_error=atof(argv[x]);
        } else if (p_switch.compare("--set")==0) {
                x++;
                p.chosen_set=atoi(argv[x]);
        } else if (p_switch.compare("--systematic")==0) {
                x++;
                p.systematic=atoi(argv[x]);
       } else if (p_switch.compare("--fix")==0) {
                x++;
                p.fix=atoi(argv[x]);
        } else {
			cout << "Incorrect usage\n ";
            cout << p_switch << "\n";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ReadFastaAli (run_params p, vector<string>& seqs) {
    ifstream ali_file;
    ali_file.open(p.seqs_file.c_str());
    vector<string> names;
    string seq;
    string str;
    for (int i=0;i<1000000;i++) {
        if (!(ali_file >> str)) break;
        if (str.at(0)=='>') {
            names.push_back(str);
            if (p.verb==1) {
                cout << "Read seqname " << str << "\n";
            }
            if (seq.size()>0) {
                seqs.push_back(seq);
                seq.clear();
            }
        } else {
            seq=seq+str;
        }
    }
    if (seq.size()>0) {
        seqs.push_back(seq);
    }
}

void ReadTimes (run_params p, vector<int>& times) {
    ifstream times_file;
    times_file.open(p.times_file.c_str());
    int n;
    for (int i=0;i<1000000;i++) {
        if (!(times_file >> n)) break;
        times.push_back(n);
    }
}


void WriteVariants (run_params& p, const vector<sparseseq>& variants) {
    if (p.verb==1) {
        for (int i=0;i<variants.size();i++) {
            cout << "i= " << i << "\n";
            for (int j=0;j<variants[i].locus.size();j++) {
                cout << variants[i].locus[j] << " " << variants[i].allele[j] << "\n";
            }
        }
    }
}

void WriteVariantsToFile (run_params& p, const vector<char>& consensus, const vector<sparseseq>& variants) {
    sparseseq allvars;
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            int match=0;
            for (int k=0;k<allvars.locus.size();k++) {
                if (variants[i].locus[j]==allvars.locus[k]) {
                    if (variants[i].allele[j]==allvars.allele[k]) {
                        match=1;
                    }
                }
            }
            if (match==0) {
                allvars.locus.push_back(variants[i].locus[j]);
                allvars.allele.push_back(variants[i].allele[j]);
            }
        }
    }
    SortVariantData(allvars);
    ofstream var_file;
    var_file.open("Variant_list.out");
    for (int i=0;i<allvars.locus.size();i++) {
        var_file << allvars.locus[i] << " " << consensus[i] << " " << allvars.allele[i] << "\n";
    }
    var_file.close();
}

void WriteVarbin(const vector< vector<double> >& varbin) {
    cout << "Varbin\n";
    for (int i=0;i<varbin.size();i++) {
        for (int j=0;j<varbin[i].size();j++) {
            cout << varbin[i][j] << " ";
        }
        cout << "\n";
    }
}

void Indexing (run_params& p, int i, modelstore m, vector<double> model_parameters, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector<int> >& nremoved) {
    cout << "Index " << i << " ";
    for (int pp=0;pp<p.sets;pp++) {
        cout << "Population " << pp << " Rate " << model_parameters[pp] << " Fixes " << gfixpos[pp].size()-gqfixpos[pp].size() << " ";
        for (int k=0;k<gfixpos[pp].size();k++) {
            int match=0;
            for (int l=0;l<gqfixpos[pp].size();l++) {
                if (gqfixpos[pp][l]==gfixpos[pp][k]) {
                    match=1;
                }
            }
            if (match==0) {
                cout << gfixpos[pp][k] << " ";
            }
        }
        cout << "Provisional " << gqfixpos[pp].size() << " ";
        for (int l=0;l<gqfixpos[pp].size();l++) {
            cout << gqfixpos[pp][l] << " ";
        }
        cout << "Removed " << nremoved[i][pp] << " ";
    }
    cout << "Likelihood " << m.lL << "\n";
}

void OutputBestModelParameters (int max_index, vector<modelstore>& outputs) {
    cout << "Best model " << max_index << "\n";
    cout << "Best parameters\n";
    for (int j=0;j<outputs[max_index].rates.size();j++) {
        cout << outputs[max_index].rates[j] << " ";
    }
    cout << outputs[max_index].error << " " << outputs[max_index].lL << "\n";
}


void WriteLimits(const vector< vector<double> >& limits) {
    cout << "Parameter uncertainties\n";
    for (int i=0;i<limits.size();i++) {
        if (i<limits.size()-1) {
            cout << "Rate  ";
        } else {
            cout << "Error ";
        }
        for (int j=0;j<limits[i].size();j++) {
            cout << limits[i][j] << " ";
        }
        cout << "\n";
    }
}






void SortVariantData (sparseseq& allvars) {
    sparseseq allvars_sort;
    vector<int> index;
    vector<int> done;
    for (int i=0;i<allvars.locus.size();i++) {
        done.push_back(0);
    }
    for (int i=0;i<allvars.locus.size();i++) {
        int min=1e7;
        int pos=-1;
        //Find first locus
        for (int j=0;j<allvars.locus.size();j++) {
            if (allvars.locus[j]<min&&done[j]==0) {
                min=allvars.locus[j];
                pos=j;
            }
        }
        allvars_sort.locus.push_back(allvars.locus[pos]);
        allvars_sort.allele.push_back(allvars.allele[pos]);
        done[pos]=1;
    }
    allvars=allvars_sort;
}


void WriteAcceptable (vector< vector<double> >& acceptable) {
    for (int i=0;i<acceptable.size();i++) {
        for (int j=0;j<acceptable[i].size();j++) {
            cout << acceptable[i][j] << " ";
        }
        cout << "\n";
    }
}

void WriteGVarbin (const vector< vector<int> >& gtimes, const vector< vector< vector<double> > >& gvarbin) {
    for (int j=0;j<gvarbin.size();j++) {
        cout << "Varbin " << j << "\n";
        for (int k=0;k<gvarbin[j].size();k++) {
            cout << gtimes[j][k] << " ";
            if (gtimes[j][k]<10) {
                cout << " ";
            }
            if (gtimes[j][k]<100) {
                cout << " ";
            }
            for (int l=0;l<gvarbin[j][k].size();l++) {
                if (gvarbin[j][k][l]<1&&gvarbin[j][k][l]>0) {
                    cout << "N ";
                } else {
                    cout << gvarbin[j][k][l] << " ";
                }
            }
            cout << "\n";
        }
    }
}


