#ifndef LAPACK_WRAPPER_H
#define LAPACK_WRAPPER_H

#include <vector>
using namespace std;

pair<vector<double>, vector<vector<double> > > eigen_symm(const vector<vector<double> > to_diag);

#endif