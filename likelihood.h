
void OptimiseMultiRateModel (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters, gsl_rng *rgen);
void DoTrack (double diff, vector<double>& track);
void CheckTrack (int& reduced, double& dx, vector<double>& track);
double FindLikelihood (int n_rates, vector<double>& rates, double error, const vector< vector<sample> >& gseq_data);
void UncertaintyMultiRateModel (run_params& p, int n_rates, int parameter, int direction, double best_likelihood, const vector< vector<sample> >& gseq_data, const vector<double>& initial_model_parameters, vector<double>& model_parameters, vector<double>& extreme_model_parameters, gsl_rng *rgen);

