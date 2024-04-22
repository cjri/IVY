
void AddSingleOrigin (vector<int>& times, vector< vector<double> >& varbin);

void FindBestModelFirstLast (run_params& p, int nq, const vector<sample>& seq_data, vector<sample>& seq_data_best, vector<double>& model_parameters, vector<double>& model_parameters_best, vector< vector<sample> >& seq_data_record, vector<int>& fixpos, vector<int>& qfixpos, vector< vector<double> >& model_parameters_record, gsl_rng *rgen);
void OptimiseSingleRateModel (run_params& p, const vector<sample>& seq_data, vector<double>& model_parameters, gsl_rng *rgen);
double GetSingleLikelihood (double& rate, double& error, const vector<sample>& seq_data);

void ModelSinglePopulation (run_params& p, vector<int> times, vector< vector<double> > varbin, gsl_rng *rgen);
void GetOptimisationData (double& best_likelihood, vector< vector<sample> >& seq_data_record, vector< vector<double> >& model_parameters_record);
void SingleRateModelExtremes (run_params& p, double best_likelihood, vector< vector<sample> >& seq_data_record, vector< vector<double> >& model_parameters_record, vector<double>& extreme_model_parameters, vector< vector<double> >& limits,gsl_rng *rgen);
void UncertaintySingleRateModel (run_params& p, int parameter, int direction, double best_likelihood, const vector<sample>& seq_data, vector<double> initial_model_parameters, vector<double>& model_parameters, vector<double>& extreme_model_parameters, gsl_rng *rgen);



void SetupModelParameters (vector<double>& model_parameters);
void FindBestModelParameters (run_params& p, int pre_subs, vector<sample>& seq_data_best, const vector< vector<sample> >& seq_data_array, int& index_best, vector<double>& model_parameters, vector<double>& model_parameters_best, vector< vector<sample> >& seq_data_record, vector<int>& fixpos, vector<int>& qfixpos, vector<int>& removed, vector< vector<double> >& model_parameters_record, gsl_rng *rgen);
