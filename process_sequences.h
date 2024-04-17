void FindConsensus (string& consensus, vector<string>& seqs);
void FindVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindAmbiguousVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindVariants2 (vector<int>& vpos, vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FixVariants (vector<string>& seqs, vector<sparseseq>& variants);
void RevertToN (vector<int>& xpos, char a1, char a2, char ax, vector<string>& seqs, vector<sparseseq>& variants);

void FindVariantPositions (run_params& p, vector<sparseseq>& variants, vector<int>& varpos);
void MakeVarmatVarbin (run_params& p, const vector<string>& seqs, const vector<int>& varpos, vector< vector<char> >& varmat, vector< vector<double> >& varbin);
void MakeVarmatVarbin2 (run_params& p, const vector<string>& seqs, const vector<int>& varpos, vector<char>& consensus, vector< vector<char> >& varmat, vector< vector<double> >& varbin);

void MakeConstantFixes (run_params& p, const vector< vector<double> >& varbin, vector< vector<double> >& constant, vector< vector<double> >& fixes);
void FindFixationSites (run_params& p, const vector< vector<double> >& varbin, const vector< vector<double> >& constant, const vector< vector<double> >& fixes, vector<int>& fixpos, vector<int>& flucpos);
void FindFixationTimes (run_params& p, const vector<int>& fixpos, const vector< vector<double> >& varbin, vector<int>& fixtimes);
void FindQFixPositions (run_params& p, const vector<int>& fixpos, const vector<int>& fixtimes, const vector< vector<double> >& varbin, vector<int>& qfixpos);
void SetupSequenceData (run_params& p, const vector<int>& times, const vector<int>& fixtimes, const vector<int>& fixpos, const vector<int>& flucpos, const vector< vector<double> >& varbin, vector<sample>& seq_data);
void CalculateSxFluc (run_params& p, const vector<double>& xf, sample& s);
void MakeBinaries (int n_size, vector< vector<int> >& binset);
void ConstructSequenceDataArray(int nq, const vector<sample>& seq_data, vector<int>& removed, vector< vector<sample> >& seq_data_array);
void SetupModelParameters (vector<double>& model_parameters);
void FindBestModelParameters (run_params& p, int pre_subs, vector<sample>& seq_data_best, const vector< vector<sample> >& seq_data_array, int& index_best, vector<double>& model_parameters, vector<double>& model_parameters_best, vector< vector<sample> >& seq_data_record, vector<int>& fixpos, vector<int>& qfixpos, vector<int>& removed, vector< vector<double> >& model_parameters_record, gsl_rng *rgen);
void InitialiseLimits (const vector<double>& model_parameters_best, vector< vector<double> >& limits);




int CheckTimes (const vector<int>& times, vector< vector<int> >& gtimes);


void CalculateSetsSystematic2 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void CalculateSetsSystematic3 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void CalculateSetsSystematic4 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
