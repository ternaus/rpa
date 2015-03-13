#ifndef COMMON_H
#define COMMON_H
#include <vector>
using namespace std;

long double fermi(const double& band, const double &temperature);
double gamma(const int&x, const int&y, const double& tkx, const double&tky);
vector<double> krange(const int&max_k, const double&a);
double random_double(const double&min, const double&max);
double band_square(const double& tkx, const double& tky, const double& mu);
double band_Lieb(const int& band_number, const double& tkx, const double& tky, const double& mu);

#endif