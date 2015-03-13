//Magnetic susceptiblity.

using namespace std;

#include "susceptibility_square.h"
#include "common.h"
#include <cmath>
#include <iostream>

double susceptibility_square(const double& temperature, 
                      const double& a, 
                      const double& mu, 
                      const int& max_kx, 
                      const int& max_ky, 
                      const double&qx, 
                      const double&qy) {
  double result = 0;

  double kx_step = 2 * M_PI / (a * max_kx);
  double ky_step = 2 * M_PI / (a * max_ky);
  double k_min = -M_PI / a;
  double k_max = M_PI / a;

  int num_sites = max_kx * max_ky;
 
  double factor = 0;

  for (double tkx = k_min; tkx < k_max; tkx += kx_step) {
    for (double tky = k_min; tky < k_max; tky += ky_step) {
      double band1 = band_square(tkx, tky, mu);
      double band2 = band_square(tkx + qx, tky + qy, mu);                
      double fm1 = fermi(band1, temperature);
      double fm2 = fermi(band2, temperature);          
      if (abs(band1 - band2) > 1e-10) {

        result += -(fm1 - fm2) / (band1 - band2);        
      }                    
      else {            
        result += fm1 * (1 - fm1) / temperature;
        
      }                
    }
  }
  return result / num_sites;
}
