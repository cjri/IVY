
void GetGSeqDataArray (run_params& p, vector<double> varbin_init, vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& gtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gqfixpos, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array);

void GetStartSeqs (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs);
void GetOriginalSeq (run_params& p, int st, const vector< vector<int> >& start_seqs, vector<double>& varbin_init);
void AddStartToGVarbin (const vector<double>& varbin_init, vector< vector< vector<double> > >& gvarbin);
void RevertGVarbin(vector< vector< vector<double> > >& gvarbin);
void AddStartToGTimes (vector< vector<int> >& gtimes);

void MakeGConstantFix (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes);
void MakeConstantFixes (run_params& p, const vector< vector<double> >& varbin, vector< vector<double> >& constant, vector< vector<double> >& fixes);

void MakeGFixFluc (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos);
void FindFixationSites (run_params& p, const vector< vector<double> >& varbin, const vector< vector<double> >& constant, const vector< vector<double> >& fixes, vector<int>& fixpos, vector<int>& flucpos);

void MakeGFixTimes (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, vector< vector<int> >& gfixtimes);
void FindFixationTimes (run_params& p, const vector<int>& fixpos, const vector< vector<double> >& varbin, vector<int>& fixtimes);

void MakeGQFixPos (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gqfixpos, vector<int>& gnq);
void FindQFixPositions (run_params& p, const vector<int>& fixpos, const vector<int>& fixtimes, const vector< vector<double> >& varbin, vector<int>& qfixpos);

void MakeGSeqData (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gtimes, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos, vector< vector<sample> >& gseq_data);
void SetupSequenceData (run_params& p, const vector<int>& times, const vector<int>& fixtimes, const vector<int>& fixpos, const vector<int>& flucpos, const vector< vector<double> >& varbin, vector<sample>& seq_data);
void CalculateSxFluc (run_params& p, const vector<double>& xf, sample& s);
void MakeBinaries (int n_size, vector< vector<int> >& binset);

void MakeGSDataArray (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array);
void ConstructSequenceDataArray(int nq, const vector<sample>& seq_data, vector<int>& removed, vector< vector<sample> >& seq_data_array);







