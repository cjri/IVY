
void GetParameters (run_params& p, int argc, const char **argv);
void ReadFastaAli (run_params p, vector<string>& seqs);
void ReadTimes (run_params p, vector<int>& times);
void WriteVariants (run_params& p, const vector<sparseseq>& variants);
void WriteVariantsToFile (run_params& p, const vector<char>& consensus, const vector<sparseseq>& variants);
void WriteVarbin(const vector< vector<double> >& varbin);
void Indexing (run_params& p, int i, modelstore m, vector<double> model_parameters, const vector< vector<int> >& gfixpos, const vector< vector<int> >& gqfixpos, const vector< vector<int> >& nremoved);
void OutputBestModelParameters (int max_index, vector<modelstore>& outputs);
void WriteLimits(const vector< vector<double> >& limits);





void SortVariantData (sparseseq& allvars);
void WriteAcceptable (vector< vector<double> >& acceptable);
void WriteGVarbin (const vector< vector<int> >& gtimes, const vector< vector< vector<double> > >& gvarbin);




