#include "ivy.h"
#include <iostream>
#include <string>
#include <cstring>

int CheckInputs (vector<string>& seqs, vector<int>& times) {
    if (seqs.size()==0) {
        cout << "Error: No sequence data\n";
        return 1;
    }
    if (times.size()!=seqs.size()) {
        cout << "Error: Mismatch between time and sequence data\n";
        return 1;
    }
    return 0;
}
