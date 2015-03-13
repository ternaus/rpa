#ifndef SUSCEPTIBILITY_SC_H
#define SUSCEPTIBILITY_SC_H
using namespace std;


double susceptibility2_sc(const int& alpha, 
                      const int& beta, 
                      const double &temperature, 
                      const double &a, 
                      const double& mu, 
                      const int& max_kx, 
                      const int&max_ky, 
                      const double&qx, 
                      const double&qy);

double susceptibility_sc(const double &temperature, 
                    const double &a, 
                    const double &mu, 
                    const int&max_kx, 
                    const int& max_ky, 
                    const double&qx, 
                    const double&qy);

#endif
