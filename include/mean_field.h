#ifndef MEAN_FIELD_H
#define MEAN_FIELD_H
#include <vector>
using namespace std;

vector<double> main_mf_loop(const double& n1_guess, const double& n2_guess, const double&m1_guess, const double&m2_guess, const int&max_kx, const int&max_ky, const double&u, const double&mu, const double&temperature);
vector<vector<double> > hamiltonian(const double& n1, const double& n2, const double&m1, const double&m2, const double&kx, const double&ky, const int&sigma, const double&u, const double&mu);
#endif