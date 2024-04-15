

struct objy {
    double y;
    vector<double> z;
};

struct objx {
    double x;
    vector<objy> y;
};

void ModelRateExploration (run_params& p, int i, const double maxL, const vector<int>& times, const vector< vector<double> >& varbin, const vector< vector< vector<int> > >& sets, vector<modelstore>& outputs, vector< vector<double> >& acceptable);
void CalculateStateSpace (run_params& p, int dim, int st, const double maxL, const vector< vector<int> >& start_seqs, const vector< vector<int> >& nremoved, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector< vector<sample> > >& gseq_data_array, vector<modelstore>& outputs, vector< vector<double> >& acceptable);
void PushVec(int dim, vector<double>& a, vector<objx>& acceptable);
void UnwindObject (int dim, vector<objx> new_input_unique, vector< vector<double> >& new_input);
int IsUniqueVec (int dim, vector<double>& temp, const vector<objx>& new_input_unique);
void AddPoint (int dim, int k, const vector<double>& initial_model_parameters, const vector< vector<double> >& new_input, const vector<objx>& acceptable, vector< vector<double> >& input);
int IsFurther(vector<double>& newpoint, const vector<double>& original, vector<double> opt);
double EvaluateLikelihood (run_params& p, int n_rates, const vector< vector<sample> >& gseq_data, vector<double>& model_parameters);




