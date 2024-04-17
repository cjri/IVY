#include "ivy.h"
#include "clustering.h"
#include "io.h"
#include "modelI.h"
#include "statespace.h"
#include <iostream>
#include <string>
#include <cstring>

void FindClusters(run_params& p, int threshold, const vector< vector<double> >& varbin, vector< vector<int> >& clusters) {
    vector<int> accounted;
    for (int i=0;i<varbin.size();i++) {
        accounted.push_back(0);
    }
    for (int i=0;i<varbin.size();i++) {
        vector<int> c;
        if (accounted[i]==0) {
            c.push_back(i);
            accounted[i]=1;
            vector<int> added;
            added.clear();
            for (int j=i+1;j<varbin.size();j++) {  //Problem here: We add j but not things attached to j
                if (accounted[j]==0) {
                    int same=0;
                    for (int k=0;k<varbin[i].size();k++) {
                        if (varbin[i][k]==1&&varbin[j][k]==1) {
                            same++;
                        }
                    }
                    if (same>=threshold) {
                        accounted[j]=1;
                        c.push_back(j);
                        added.push_back(j);

                    }
                }
            }
            while (added.size()>0) {
                vector<int> added_temp;
                //Go through the things that have been added.  Add links to c and links to added.
                for (int a=0;a<added.size();a++) {
                    for (int j=i+1;j<varbin.size();j++) {
                        if (accounted[j]==0) {
                            int same=0;
                            for (int k=0;k<varbin[added[a]].size();k++) {
                                if (varbin[added[a]][k]==1&&varbin[j][k]==1) {
                                    same++;
                                }
                            }
                            if (same>=threshold) {
                                accounted[j]=1;
                                c.push_back(j);
                                added_temp.push_back(j);
                            }
                        }
                    }
                }
                added.clear();
                added=added_temp;
            }
            clusters.push_back(c);
        }
    }
    //Sort clusters - find unique elements
    for (int i=0;i<clusters.size();i++) {
        sort(clusters[i].begin(),clusters[i].end());
        clusters[i].erase(unique(clusters[i].begin(),clusters[i].end()),clusters[i].end());
    }
    
    if (p.verb==1) {
        cout << "Clusters\n";
        for (int i=0;i<clusters.size();i++) {
            for (int j=0;j<clusters[i].size();j++) {
                cout << clusters[i][j] << " ";
            }
            cout << "\n";
            for (int j=0;j<clusters[i].size();j++) {
                for (int k=0;k<varbin[clusters[i][j]].size();k++) {
                    cout << varbin[clusters[i][j]][k] << " ";
                }
                cout << "\n";
            }
        }
    }
}

void CompileSets (run_params& p, const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    if (p.sets==2) {
        CalculateSetsSystematic2 (clusters,sets);
    }
    if (p.sets==3) {
        CalculateSetsSystematic3 (clusters,sets);
    }
    if (p.sets==4) {
        CalculateSetsSystematic4 (clusters,sets);
    }
    ConvertSetsClusters (clusters,sets);
    if (p.verb==1) {
        if (p.sets==2) {
            cout << "Sets size " << sets.size() << "\n";
            for (int i=0;i<sets.size();i++) {
                cout << i << " ";
                for (int j=0;j<sets[i][0].size();j++) {
                    cout << sets[i][0][j] << " ";
                }
                cout << "\n";
            }
        }
        if (p.sets>=3) {
            cout << "Sets size " << sets.size() << "\n";
            for (int i=0;i<sets.size();i++) {
                cout << i << " ";
                for (int j=0;j<sets[i].size();j++) {
                    cout << "(";
                    for (int k=0;k<sets[i][j].size();k++) {
                        cout << sets[i][j][k] << " ";
                    }
                    cout << ") ";
                }
                cout << "\n";
            }
        }
    }
}

void ConvertSetsClusters (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    vector< vector< vector<int> > > sets_temp;
    for (int i=0;i<sets.size();i++) {
        vector< vector<int> > ss;
        for (int j=0;j<sets[i].size();j++) {
            vector<int> s;
            for (int k=0;k<sets[i][j].size();k++) {
                for (int l=0;l<clusters[sets[i][j][k]].size();l++) {
                    s.push_back(clusters[sets[i][j][k]][l]);
                }
            }
            sort(s.begin(),s.end());
            ss.push_back(s);
        }
        sets_temp.push_back(ss);
    }
    sets=sets_temp;
}






















void CalculateSetsSystematic2 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/2;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector< vector<int> > sets1;
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
                sets.push_back(sets1);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector< vector<int> > sets1;
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                    sets.push_back(sets1);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector< vector<int> > sets1;
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                        sets.push_back(sets1);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector< vector<int> > sets1;
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                            sets.push_back(sets1);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector< vector<int> > sets1;
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                                sets.push_back(sets1);
                            }
                        }
                    }
                }
            }
        }
        if (n==6) {
            for (int i1=0;i1<clusters.size()-5;i1++) {
                for (int i2=i1+1;i2<clusters.size()-4;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-3;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-2;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-1;i5++) {
                                for (int i6=i5+1;i6<clusters.size();i6++) {
                                    vector< vector<int> > sets1;
                                    vector<int> s;
                                    s.push_back(i1);
                                    s.push_back(i2);
                                    s.push_back(i3);
                                    s.push_back(i4);
                                    s.push_back(i5);
                                    s.push_back(i6);
                                    sets1.push_back(s);
                                    sets.push_back(sets1);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==7) {
            for (int i1=0;i1<clusters.size()-6;i1++) {
                for (int i2=i1+1;i2<clusters.size()-5;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-4;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-3;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-2;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-1;i6++) {
                                    for (int i7=i6+1;i7<clusters.size();i7++) {
                                        vector< vector<int> > sets1;
                                        vector<int> s;
                                        s.push_back(i1);
                                        s.push_back(i2);
                                        s.push_back(i3);
                                        s.push_back(i4);
                                        s.push_back(i5);
                                        s.push_back(i6);
                                        s.push_back(i7);
                                        sets1.push_back(s);
                                        sets.push_back(sets1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==8) {
            for (int i1=0;i1<clusters.size()-7;i1++) {
                for (int i2=i1+1;i2<clusters.size()-6;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-5;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-4;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-3;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-2;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-1;i7++) {
                                        for (int i8=i7+1;i8<clusters.size();i8++) {
                                            vector< vector<int> > sets1;
                                            vector<int> s;
                                            s.push_back(i1);
                                            s.push_back(i2);
                                            s.push_back(i3);
                                            s.push_back(i4);
                                            s.push_back(i5);
                                            s.push_back(i6);
                                            s.push_back(i7);
                                            s.push_back(i8);
                                            sets1.push_back(s);
                                            sets.push_back(sets1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==9) {
            for (int i1=0;i1<clusters.size()-8;i1++) {
                for (int i2=i1+1;i2<clusters.size()-7;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-6;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-5;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-4;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-3;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-2;i7++) {
                                        for (int i8=i7+1;i8<clusters.size()-1;i8++) {
                                            for (int i9=i8+1;i9<clusters.size();i9++) {
                                                vector< vector<int> > sets1;
                                                vector<int> s;
                                                s.push_back(i1);
                                                s.push_back(i2);
                                                s.push_back(i3);
                                                s.push_back(i4);
                                                s.push_back(i5);
                                                s.push_back(i6);
                                                s.push_back(i7);
                                                s.push_back(i8);
                                                s.push_back(i9);
                                                sets1.push_back(s);
                                                sets.push_back(sets1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (n==10) {
            for (int i1=0;i1<clusters.size()-9;i1++) {
                for (int i2=i1+1;i2<clusters.size()-8;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-7;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-6;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-5;i5++) {
                                for (int i6=i5+1;i6<clusters.size()-4;i6++) {
                                    for (int i7=i6+1;i7<clusters.size()-3;i7++) {
                                        for (int i8=i7+1;i8<clusters.size()-2;i8++) {
                                            for (int i9=i8+1;i9<clusters.size()-1;i9++) {
                                                for (int ia=i9+1;ia<clusters.size();ia++) {
                                                    vector< vector<int> > sets1;
                                                    vector<int> s;
                                                    s.push_back(i1);
                                                    s.push_back(i2);
                                                    s.push_back(i3);
                                                    s.push_back(i4);
                                                    s.push_back(i5);
                                                    s.push_back(i6);
                                                    s.push_back(i7);
                                                    s.push_back(i8);
                                                    s.push_back(i9);
                                                    s.push_back(ia);
                                                    sets1.push_back(s);
                                                    sets.push_back(sets1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void CalculateSetsSystematic3 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/3;
    vector< vector<int> > sets1;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                            }
                        }
                    }
                }
            }
        }
        if (n==6) {
            for (int i1=0;i1<clusters.size()-5;i1++) {
                for (int i2=i1+1;i2<clusters.size()-4;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-3;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-2;i4++) {
                            for (int i5=i4+1;i5<clusters.size()-1;i5++) {
                                for (int i6=i5+1;i6<clusters.size();i6++) {
                                    vector<int> s;
                                    s.push_back(i1);
                                    s.push_back(i2);
                                    s.push_back(i3);
                                    s.push_back(i4);
                                    s.push_back(i5);
                                    s.push_back(i6);
                                    sets1.push_back(s);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "First set " << sets1.size() << " possibilities\n";
    for (int i=0;i<sets1.size();i++) {
        int min=sets1[i].size();
        int r=clusters.size()-min;
        int max=r/2;
        //Construct remainder
        vector<int> remainder;
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==sets1[i][index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2;
                    s2.push_back(sets1[i]);
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2;
                        s2.push_back(sets1[i]);
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-2;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-1;i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2;
                            s2.push_back(sets1[i]);
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-3;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-2;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-1;i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2;
                                s2.push_back(sets1[i]);
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-4;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-3;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-2;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-1;i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2;
                                    s2.push_back(sets1[i]);
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-5;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-4;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-3;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-2;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-1;i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2;
                                        s2.push_back(sets1[i]);
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-6;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-5;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-4;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-3;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-2;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-1;i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2;
                                            s2.push_back(sets1[i]);
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==8) {
                for (int i1=0;i1<remainder.size()-7;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-6;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-5;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-4;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-3;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-2;i6++) {
                                        for (int i7=i6+1;i7<remainder.size()-1;i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                vector< vector<int> > s2;
                                                s2.push_back(sets1[i]);
                                                vector<int> s;
                                                s.push_back(remainder[i1]);
                                                s.push_back(remainder[i2]);
                                                s.push_back(remainder[i3]);
                                                s.push_back(remainder[i4]);
                                                s.push_back(remainder[i5]);
                                                s.push_back(remainder[i6]);
                                                s.push_back(remainder[i7]);
                                                s.push_back(remainder[i8]);
                                                s2.push_back(s);
                                                sets.push_back(s2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==9) {
                for (int i1=0;i1<remainder.size()-8;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-7;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-6;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-5;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-4;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-3;i6++) {
                                        for (int i7=i6+1;i7<remainder.size()-2;i7++) {
                                            for (int i8=i7+1;i8<remainder.size()-1;i8++) {
                                                for (int i9=i8+1;i9<remainder.size();i9++) {
                                                    vector< vector<int> > s2;
                                                    s2.push_back(sets1[i]);
                                                    vector<int> s;
                                                    s.push_back(remainder[i1]);
                                                    s.push_back(remainder[i2]);
                                                    s.push_back(remainder[i3]);
                                                    s.push_back(remainder[i4]);
                                                    s.push_back(remainder[i5]);
                                                    s.push_back(remainder[i6]);
                                                    s.push_back(remainder[i7]);
                                                    s.push_back(remainder[i8]);
                                                    s.push_back(remainder[i9]);
                                                    s2.push_back(s);
                                                    sets.push_back(s2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void CalculateSetsSystematic4 (const vector< vector<int> >& clusters, vector< vector< vector<int> > >& sets) {
    int max=clusters.size()/4; //Smallest
    vector< vector<int> > sets1;
    for (int n=1;n<=max;n++) {
        if (n==1) {
            for (int i1=0;i1<clusters.size();i1++) {
                vector<int> s;
                s.push_back(i1);
                sets1.push_back(s);
            }
        }
        if (n==2) {
            for (int i1=0;i1<clusters.size()-1;i1++) {
                for (int i2=i1+1;i2<clusters.size();i2++) {
                    vector<int> s;
                    s.push_back(i1);
                    s.push_back(i2);
                    sets1.push_back(s);
                }
            }
        }
        if (n==3) {
            for (int i1=0;i1<clusters.size()-2;i1++) {
                for (int i2=i1+1;i2<clusters.size()-1;i2++) {
                    for (int i3=i2+1;i3<clusters.size();i3++) {
                        vector<int> s;
                        s.push_back(i1);
                        s.push_back(i2);
                        s.push_back(i3);
                        sets1.push_back(s);
                    }
                }
            }
        }
        if (n==4) {
            for (int i1=0;i1<clusters.size()-3;i1++) {
                for (int i2=i1+1;i2<clusters.size()-2;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-1;i3++) {
                        for (int i4=i3+1;i4<clusters.size();i4++) {
                            vector<int> s;
                            s.push_back(i1);
                            s.push_back(i2);
                            s.push_back(i3);
                            s.push_back(i4);
                            sets1.push_back(s);
                        }
                    }
                }
            }
        }
        if (n==5) {
            for (int i1=0;i1<clusters.size()-4;i1++) {
                for (int i2=i1+1;i2<clusters.size()-3;i2++) {
                    for (int i3=i2+1;i3<clusters.size()-2;i3++) {
                        for (int i4=i3+1;i4<clusters.size()-1;i4++) {
                            for (int i5=i4+1;i5<clusters.size();i5++) {
                                vector<int> s;
                                s.push_back(i1);
                                s.push_back(i2);
                                s.push_back(i3);
                                s.push_back(i4);
                                s.push_back(i5);
                                sets1.push_back(s);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "First set " << sets1.size() << " possibilities\n";
    for (int i=0;i<sets1.size();i++) {
        int min=sets1[i].size(); //Set 1 is the smallest
        int r=clusters.size()-min;
        int max=r/3;
        //Construct remainder - clusters not in a set
        vector<int> remainder;
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==sets1[i][index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2;
                    s2.push_back(sets1[i]);
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2;
                        s2.push_back(sets1[i]);
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-2;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-1;i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2;
                            s2.push_back(sets1[i]);
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-3;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-2;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-1;i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2;
                                s2.push_back(sets1[i]);
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-4;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-3;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-2;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-1;i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2;
                                    s2.push_back(sets1[i]);
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-5;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-4;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-3;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-2;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-1;i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2;
                                        s2.push_back(sets1[i]);
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-6;i1++) {
                    for (int i2=i1+1;i2<remainder.size()-5;i2++) {
                        for (int i3=i2+1;i3<remainder.size()-4;i3++) {
                            for (int i4=i3+1;i4<remainder.size()-3;i4++) {
                                for (int i5=i4+1;i5<remainder.size()-2;i5++) {
                                    for (int i6=i5+1;i6<remainder.size()-1;i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2;
                                            s2.push_back(sets1[i]);
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    vector< vector< vector<int> > > sets_new;
    //For each of the sets consider the next division
    for (int i=0;i<sets.size();i++) {
        int min=sets[i][1].size(); //Set 2 is the next smallest
        int r=clusters.size()-(sets[i][0].size()+sets[i][1].size());
        int max=r/2;
        vector<int> remainder;
        vector<int> used=sets[i][0];
        for (int j=0;j<sets[i][1].size();j++) {
            used.push_back(sets[i][1][j]);
        }
        sort(used.begin(),used.end());
        int index=0;
        for (int j=0;j<clusters.size();j++) {
            if (j==used[index]) {
                index++;
            } else {
                remainder.push_back(j);
            }
        }
        for (int n=min;n<=max;n++) {
            if (n==1) {
                for (int i1=0;i1<remainder.size();i1++) {
                    vector< vector<int> > s2=sets[i];
                    vector<int> s;
                    s.push_back(remainder[i1]);
                    s2.push_back(s);
                    sets_new.push_back(s2);
                }
            }
            if (n==2) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        vector< vector<int> > s2=sets[i];
                        vector<int> s;
                        s.push_back(remainder[i1]);
                        s.push_back(remainder[i2]);
                        s2.push_back(s);
                        sets_new.push_back(s2);
                    }
                }
            }
            if (n==3) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            vector< vector<int> > s2=sets[i];
                            vector<int> s;
                            s.push_back(remainder[i1]);
                            s.push_back(remainder[i2]);
                            s.push_back(remainder[i3]);
                            s2.push_back(s);
                            sets_new.push_back(s2);
                        }
                    }
                }
            }
            if (n==4) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                vector< vector<int> > s2=sets[i];
                                vector<int> s;
                                s.push_back(remainder[i1]);
                                s.push_back(remainder[i2]);
                                s.push_back(remainder[i3]);
                                s.push_back(remainder[i4]);
                                s2.push_back(s);
                                sets_new.push_back(s2);
                            }
                        }
                    }
                }
            }
            if (n==5) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    vector< vector<int> > s2=sets[i];
                                    vector<int> s;
                                    s.push_back(remainder[i1]);
                                    s.push_back(remainder[i2]);
                                    s.push_back(remainder[i3]);
                                    s.push_back(remainder[i4]);
                                    s.push_back(remainder[i5]);
                                    s2.push_back(s);
                                    sets_new.push_back(s2);
                                }
                            }
                        }
                    }
                }
            }
            if (n==6) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        vector< vector<int> > s2=sets[i];
                                        vector<int> s;
                                        s.push_back(remainder[i1]);
                                        s.push_back(remainder[i2]);
                                        s.push_back(remainder[i3]);
                                        s.push_back(remainder[i4]);
                                        s.push_back(remainder[i5]);
                                        s.push_back(remainder[i6]);
                                        s2.push_back(s);
                                        sets_new.push_back(s2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==7) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            vector< vector<int> > s2=sets[i];
                                            vector<int> s;
                                            s.push_back(remainder[i1]);
                                            s.push_back(remainder[i2]);
                                            s.push_back(remainder[i3]);
                                            s.push_back(remainder[i4]);
                                            s.push_back(remainder[i5]);
                                            s.push_back(remainder[i6]);
                                            s.push_back(remainder[i7]);
                                            s2.push_back(s);
                                            sets_new.push_back(s2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==8) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                vector< vector<int> > s2=sets[i];
                                                vector<int> s;
                                                s.push_back(remainder[i1]);
                                                s.push_back(remainder[i2]);
                                                s.push_back(remainder[i3]);
                                                s.push_back(remainder[i4]);
                                                s.push_back(remainder[i5]);
                                                s.push_back(remainder[i6]);
                                                s.push_back(remainder[i7]);
                                                s.push_back(remainder[i8]);
                                                s2.push_back(s);
                                                sets_new.push_back(s2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (n==9) {
                for (int i1=0;i1<remainder.size()-1;i1++) {
                    for (int i2=i1+1;i2<remainder.size();i2++) {
                        for (int i3=i2+1;i3<remainder.size();i3++) {
                            for (int i4=i3+1;i4<remainder.size();i4++) {
                                for (int i5=i4+1;i5<remainder.size();i5++) {
                                    for (int i6=i5+1;i6<remainder.size();i6++) {
                                        for (int i7=i6+1;i7<remainder.size();i7++) {
                                            for (int i8=i7+1;i8<remainder.size();i8++) {
                                                for (int i9=i8+1;i9<remainder.size();i9++) {
                                                    vector< vector<int> > s2=sets[i];
                                                    vector<int> s;
                                                    s.push_back(remainder[i1]);
                                                    s.push_back(remainder[i2]);
                                                    s.push_back(remainder[i3]);
                                                    s.push_back(remainder[i4]);
                                                    s.push_back(remainder[i5]);
                                                    s.push_back(remainder[i6]);
                                                    s.push_back(remainder[i7]);
                                                    s.push_back(remainder[i8]);
                                                    s.push_back(remainder[i9]);
                                                    s2.push_back(s);
                                                    sets_new.push_back(s2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    sets=sets_new;
}
