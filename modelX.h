
void RunModelX (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen);
void ExploreSetX (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen);

void CompileSetsX (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);

    
void GetStartSeqsX (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs);
void MakeGSDataArrayX (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector< vector<sample> > >& gseq_data_array);
void CalculateBestModelsX (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen);
void InitialiseLimitsX (const vector<double>& model_parameters_best, vector< vector<double> >& limits);
void ModelXRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen);
