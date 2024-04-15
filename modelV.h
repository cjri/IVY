
void RunModelV (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen);
void ExploreSet (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, vector<modelstore>& outputs, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen);

void FindClusters(run_params& p, int threshold, const vector< vector<double> >& varbin, vector< vector<int> >& clusters);
void CompileSets (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void ConvertSetsClusters (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void AssignDataToSetsG (run_params& p, int i, const vector< vector< vector<int> > >& sets, const vector<int>& times, const vector< vector<double> >& varbin, vector< vector<int> >& gtimes, vector< vector< vector<double> > >& gvarbin);
void GetStartSeqsV (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs);
void AddStartToGVarbin (const vector<double>& varbin_init, vector< vector< vector<double> > >& gvarbin);
void RevertGVarbin(vector< vector< vector<double> > >& gvarbin);
void GetOriginalSeq (run_params& p, int st, const vector< vector<int> >& start_seqs, vector<double>& varbin_init);
void GetOriginalSeq2 (run_params& p, const vector<int>& start_seqs, vector<double>& varbin_init);
void AddStartToGTimes (vector< vector<int> >& gtimes);
void MakeGConstantFix (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes);
void MakeGFixFluc (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector< vector<double> > >& gconstant, vector< vector< vector<double> > >& gfixes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos);
void MakeGFixTimes (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, vector< vector<int> >& gfixtimes);
void MakeGQFixPos (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gqfixpos, vector<int>& gnq);
void MakeGSeqData (run_params& p, const vector< vector< vector<double> > >& gvarbin, const vector< vector<int> >& gtimes, const vector< vector<int> >& gfixtimes, vector< vector<int> >& gfixpos, vector< vector<int> >& gflucpos, vector< vector<sample> >& gseq_data);
void MakeGSDataArrayV (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array);
void CalculateBestModelsV (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen);


double FindBestModel (int& max_index, const vector<modelstore>& outputs);
void RefineModels (double maxL, vector<modelstore>& outputs);
void GetBestModelParameters (double maxL, const vector<modelstore>& outputs, vector<double>& model_parameters_best);

void OptimiseMultiRateModel (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters, gsl_rng *rgen);
void InitialiseLimitsV (const vector<double>& model_parameters_best, vector< vector<double> >& limits);
void ModelVRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen);

void UncertaintyMultiRateModel (run_params& p, int n_rates, int parameter, int direction, double best_likelihood, const vector< vector<sample> >& gseq_data, const vector<double>& initial_model_parameters, vector<double>& model_parameters, vector<double>& extreme_model_parameters, gsl_rng *rgen);


