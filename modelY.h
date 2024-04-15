
void RunModelY (run_params& p, const vector<int>& times, const vector< vector<double> >& varbin, gsl_rng *rgen);
void ExploreSetY (run_params& p, int i, int& best_set, double& best_log, const vector<int>& times, vector<modelstore>& outputs, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, gsl_rng *rgen);

void CompileSetsY (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void GetStartSeqsY (run_params& p, const vector< vector< vector<double> > >& gvarbin, vector< vector<int> >& start_seqs);
void MakeGSDataArrayY (run_params& p, const vector<int>& gnq, const vector< vector<sample> >& gseq_data, vector< vector<int> >& nremoved, vector< vector< vector<sample> > >& gseq_data_array);
void CalculateBestModelsY (run_params& p, int st, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, gsl_rng *rgen);
void InitialiseLimitsY (const vector<double>& model_parameters_best, vector< vector<double> >& limits);
void ModelYRateExtremes (run_params& p, const double maxL, const vector< vector< vector<double> > >& gvarbin_orig, const vector< vector<int> >& gtimes_orig, const vector<modelstore>& outputs, vector< vector<double> >& limits, gsl_rng *rgen);
void ModelYRateExploration (run_params& p, int i, const double maxL, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, vector<modelstore>& outputs, vector< vector<double> >& acceptable);
int IsFurther(vector<double>& newpoint, const vector<double>& original, vector<double> opt);




int GetUnique (vector<double>& temp, const vector< vector<double> >& acceptable);
