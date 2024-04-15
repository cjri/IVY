#include "seqmodel.h"
#include "io.h"
#include "modelV.h"
#include "modelY.h"
#include "statespace.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>



void ModelRateExploration (run_params& p, int i, const double maxL, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, vector<modelstore>& outputs, vector< vector<double> >& acceptable) {
    
    cout << "Rate exploration\n";
    
    //cout << "Set " << i << "\n";
    cout << "Max likelihood " << maxL << "\n";

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
        //cout << "Running inferences...\n";
        vector< vector<int> > start_seqs;
        if (p.model.compare("Y")==0) {
            //cout << "Diagnosed model Y\n";
            GetStartSeqsY (p,gvarbin,start_seqs);
        }
        if (p.model.compare("V")==0) {
            //cout << "Diagnosed model V\n";
            GetStartSeqsV (p,gvarbin,start_seqs);
        }
        
        vector<double> varbin_init;
        for (int st=0;st<start_seqs.size();st++){  //Loop over potential start points
            cout << "Start point " << st << " " << start_seqs.size() << "\n";
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
            if (p.model.compare("Y")==0) {
                MakeGSDataArrayY (p,gnq,gseq_data,nremoved,gseq_data_array);
            } else if (p.model.compare("V")==0) {
                MakeGSDataArrayV (p,gnq,gseq_data,nremoved,gseq_data_array);
            }

            //Change the code here: Use data from previous outputs and do the calculation in discrete space.  NB Need to import the outputs.

            if (p.model.compare("Y")==0) {
                CalculateStateSpace(p,3,st,maxL,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,acceptable);
            } else if (p.model.compare("V")==0) {
                CalculateStateSpace(p,2,st,maxL,start_seqs,nremoved,gfixpos,gqfixpos,gseq_data_array,outputs,acceptable);
            }
            varbin_init.clear();
            gconstant.clear();
            gfixes.clear();
            gfixpos.clear();
            gflucpos.clear();
            gfixtimes.clear();
            gqfixpos.clear();
            gnq.clear();
            gseq_data.clear();
            nremoved.clear();
            gseq_data_array.clear();
        }
            
    } else {
        cout << "Set excluded\n";
    }
}


void CalculateStateSpace (run_params& p, int dim, int st, const double maxL, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, vector< vector<double> >& acceptable) {
    vector<double> model_parameters;
    vector<double> a;
    vector<objx> acceptable_temp;
    for (int i=0;i<gseq_data_array.size();i++) {
        cout << "i " << i << " " << gseq_data_array.size() << "\n";
        modelstore m;
        m.start_seq=start_seqs[st];
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
        if (check!=0) { //Don't have a gain of variants in zero time.
            //Parameters as rounded versions of the previous parameters
            for (int j=0;j<=dim-1;j++) {
                model_parameters.push_back(outputs[i].rates[j]);
                model_parameters[j]=floor((model_parameters[j]*500)+0.5);
                model_parameters[j]=(model_parameters[j]+0.)/500.;
                model_parameters[j]=max(0.002,model_parameters[j]);
            }
            model_parameters.push_back(outputs[i].error);
            vector<double> initial_model_parameters(3);
            for (int j=0;j<model_parameters.size();j++) {
                initial_model_parameters[j]=model_parameters[j];
            }
            //First step.  Set up input object
            vector< vector<double> > input; //Parameters to test
            vector<double> in(dim);
            for (int j=0;j<=dim-1;j++) {
                in[j]=model_parameters[j];
            }
            input.push_back(in);
            vector< vector<double> > new_input; //Parameters to test
            int added=1;
            while (added==1) {
                added=0;
                for (int k=0;k<input.size();k++) {
                    for (int l=0;l<input[k].size();l++) {
                        model_parameters[l]=input[k][l];
                    }
                    double lL=EvaluateLikelihood (p,dim,gseq_data_array[i],model_parameters);
                    for (int m=0;m<model_parameters.size();m++) {
                        cout << model_parameters[m] << " ";
                    }
                    cout << lL << "\n";
                    if (lL>maxL-2) {
                        added=1;
                        a.clear();
                        for (int l=0;l<input[k].size();l++) {
                            a.push_back(input[k][l]);
                        }
                        new_input.push_back(a);
                        //Push a into acceptable_temp
                        PushVec(dim,a,acceptable_temp);
                        
                    }
                }
                input.clear();
                vector<objx> new_input_unique;
                for (int k=0;k<new_input.size();k++) {
                    vector<double> temp=new_input[k];
                    int unq=IsUniqueVec(dim,temp,new_input_unique);
                    if (unq==1) {
                        PushVec(dim,temp,new_input_unique);
                    }
                }
                new_input.clear();
                UnwindObject (dim,new_input_unique,new_input);
                new_input_unique.clear();
                for (int k=0;k<new_input.size();k++) {
                    AddPoint (dim,k,initial_model_parameters,new_input,acceptable_temp,input);
                }
                new_input.clear();
            }
        }
    }
    //Transfer acceptable_temp to acceptable
    UnwindObject (dim,acceptable_temp,acceptable);
    acceptable_temp.clear();
}

int IsUniqueVec (int dim, vector<double>& temp, const vector<objx>& new_input_unique) {
    int uniq=1;
    if (dim==3) {
        for (int l=0;l<new_input_unique.size();l++) {
            if (new_input_unique[l].x==temp[0]) {
                for (int m=0;m<new_input_unique[l].y.size();m++) {
                    if (new_input_unique[l].y[m].y==temp[1]) {
                        for (int n=0;n<new_input_unique[l].y[m].z.size();n++) {
                            if (new_input_unique[l].y[m].z[n]==temp[2]) {
                                uniq=0;
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else if (dim==2) {
        for (int l=0;l<new_input_unique.size();l++) {
            if (new_input_unique[l].x==temp[0]) {
                for (int m=0;m<new_input_unique[l].y.size();m++) {
                    if (new_input_unique[l].y[m].y==temp[1]) {
                        uniq=0;
                        break;
                    }
                }
            }
        }
    }
    return uniq;
}

void PushVec(int dim, vector<double>& a, vector<objx>& acceptable) {
    if (dim==2) {
        //cout << "Push here\n";
        int foundx=0;
        for (int l=0;l<acceptable.size();l++) {
            if (acceptable[l].x==a[0]) {
                foundx=1;
                objy y_new;
                y_new.y=a[1];
                y_new.z.push_back(0);
                acceptable[l].y.push_back(y_new);
            }
        }
        //cout << "Foundx " << foundx << "\n";
        if (foundx==0) {
            objx newxo;
            newxo.x=a[0];
            objy y_new;
            y_new.y=a[1];
            y_new.z.push_back(0);
            newxo.y.push_back(y_new);
            acceptable.push_back(newxo);
        }
    } else if (dim==3) {
        int foundx=0;
        for (int l=0;l<acceptable.size();l++) {
            if (acceptable[l].x==a[0]) {
                foundx=1;
                int foundy=0;
                for (int m=0;m<acceptable[l].y.size();m++) {
                    if (acceptable[l].y[m].y==a[1]) {
                        foundy=1;
                        acceptable[l].y[m].z.push_back(a[2]);
                    }
                }
                if (foundy==0) {
                    objy newo;
                    newo.y=a[1];
                    newo.z.push_back(a[2]);
                    acceptable[l].y.push_back(newo);
                }
            }
        }
        if (foundx==0) {
            objx newxo;
            newxo.x=a[0];
            objy newo;
            newo.y=a[1];
            newo.z.push_back(a[2]);
            newxo.y.push_back(newo);
            acceptable.push_back(newxo);
        }
    }
}

void UnwindObject (int dim, vector<objx> new_input_unique, vector< vector<double> >& new_input) {
    if (dim==3) {
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
    } else if (dim==2) {
        for (int l=0;l<new_input_unique.size();l++) {
            for (int m=0;m<new_input_unique[l].y.size();m++) {
                vector<double> temp;
                temp.push_back(new_input_unique[l].x);
                temp.push_back(new_input_unique[l].y[m].y);
                new_input.push_back(temp);
            }
        }
    }
}


void AddPoint (int dim, int k, const vector<double>& initial_model_parameters, const vector< vector<double> >& new_input, const vector<objx>& acceptable, vector< vector<double> >& input) {
    if (dim==3) {
        for (int i=0;i<6;i++) {
            vector<double> temp=new_input[k];
            if (i==0) {
                temp[0]=temp[0]+0.002;
            } else if (i==1) {
                temp[0]=max(0.002,temp[0]-0.002);
            } else if (i==2) {
                temp[1]=temp[1]+0.002;
            } else if (i==3) {
                temp[1]=max(0.002,temp[1]-0.002);
            } else if (i==4) {
                temp[2]=temp[2]+0.002;
            } else if (i==5) {
                temp[2]=max(0.002,temp[2]-0.002);
            }
            int isf=IsFurther(temp,new_input[k],initial_model_parameters);
            if (isf==1) {
                int unq=IsUniqueVec (dim,temp,acceptable);
                if (unq==1) {
                    input.push_back(temp);
                }
            }
        }
    } else if (dim==2) {
        for (int i=0;i<4;i++) {
            //cout << "i= " << i << "\n";
            //cout << "Size " << new_input[k].size() << "\n";
            vector<double> temp;
            for (int j=0;j<new_input[k].size();j++) {
                temp.push_back(new_input[k][j]);
            }
            if (i==0) {
                temp[0]=temp[0]+0.002;
            } else if (i==1) {
                temp[0]=max(0.002,temp[0]-0.002);
            } else if (i==2) {
                temp[1]=temp[1]+0.002;
            } else if (i==3) {
                temp[1]=max(0.002,temp[1]-0.002);
            }
            //cout << "Call IsFurther\n";
            int isf=IsFurther(temp,new_input[k],initial_model_parameters);
            //cout << "Is F " << isf << "\n";
            if (isf==1) {
                int unq=IsUniqueVec (dim,temp,acceptable);
                if (unq==1) {
                    input.push_back(temp);
                }
            }
        }
    }
}


int IsFurther(vector<double>& newpoint, const vector<double>& original, vector<double> opt) {
    double metric1=0;
    double metric2=0;
    //cout << "ISF size " << newpoint.size() << " " << original.size() << " " << opt.size() << "\n";
    for (int i=0;i<newpoint.size();i++) {
        metric1=metric1+abs(original[i]-opt[i]);
    }
    for (int i=0;i<newpoint.size();i++) {
        metric2=metric2+abs(newpoint[i]-opt[i]);
    }
    if (metric2>metric1) {
        return 1;
    } else {
        return 0;
    }
}

double EvaluateLikelihood (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters) {
    //cout << "FE " << p.fix_error << "\n";
    if (p.fix_error>-1) {
        model_parameters[n_rates]=p.fix_error;
    }
    //Evaluate the likelihood
    double lL=0;
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

                /*if (gseq_data[i][j].nfluc>=0||gseq_data[i][j].xfluc.size()>0) {
                    if (gseq_data[i][j].xfluc.size()==0) {
                        lL=lL+log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc,model_parameters[n_rates]));
                    } else {
                        //Account for uncertainty in the number of flucutations arising from ambiguous nucleotides
                        for (int j=0;j<gseq_data[i][j].xfluc.size();j++) {
                            lL=lL+gseq_data[i][j].xfluc[j]*log(gsl_ran_poisson_pdf(gseq_data[i][j].nfluc+j,model_parameters[n_rates]));
                        }
                    }
                }*/
            }
        }
    }
    
    
    

    return(lL);
}
