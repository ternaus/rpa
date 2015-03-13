#ifndef SUSEPTIBILITY_H
#define SUSEPTIBILITY_H
using namespace std;

double susceptibility2(const int& alpha, 
                      const int& beta, 
                      const double &temperature, 
                      const double &a, 
                      const double& mu, 
                      const int& max_kx, 
                      const int&max_ky, 
                      const double&qx, 
                      const double&qy);

double susceptibility(const double &temperature, 
                    const double &a, 
                    const double &mu, 
                    const int&max_kx, 
                    const int& max_ky, 
                    const double&qx, 
                    const double&qy);

#endif
