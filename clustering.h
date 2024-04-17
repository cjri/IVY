
void FindClusters(run_params& p, int threshold, const vector< vector<double> >& varbin, vector< vector<int> >& clusters);
void CompileSets (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void CalculateSetsSystematic2 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void CalculateSetsSystematic3 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void CalculateSetsSystematic4 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
void ConvertSetsClusters (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets);
