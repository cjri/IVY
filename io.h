
void GetParameters (run_params& p, int argc, const char **argv);
void ReadFastaAli (run_params p, vector<string>& seqs);
void ReadTimes (run_params p, vector<int>& times);
void WriteVariants (run_params& p, const vector<sparseseq>& variants);
void WriteVariantsToFile (run_params& p, const vector<char>& consensus, const vector<sparseseq>& variants);
void SortVariantData (sparseseq& allvars);
void WriteLimits(const vector< vector<double> >& limits);
void WriteAcceptable (vector< vector<double> >& acceptable);
void WriteGVarbin (const vector< vector<int> >& gtimes, const vector< vector< vector<double> > >& gvarbin);
void OutputBestModelParametersV (int max_index, vector<modelstore>& outputs);

