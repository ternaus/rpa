using namespace std;

#include "susceptibility_sc_square.h"
#include "common.h"
#include <cmath>
#include <iostream>


static const int num_bands = 1;

double susceptibility_sc_square(const double& temperature, 
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

  int num_sites = num_bands * max_kx * max_ky;
 
  double factor = 0;

  for (double tkx = k_min; tkx < k_max; tkx += kx_step) {
    for (double tky = k_min; tky < k_max; tky += ky_step) {
      double band1 = band_square(tkx, tky, mu);
      double band2 = band_square(-tkx -qx, -tky - qy, mu);                
      double fm1 = fermi(band1, temperature);
      double fm2 = fermi(band2, temperature);          
      if (abs(band1 + band2) > 1e-10) {
        double C1 = (exp( (band1 + band2) / temperature) - 1) / (band1 + band2);
        // double C2 = -(exp(-(band1 + band2) / temperature) - 1) / (band1 + band2);
        // factor = (C1 * fm1 * fm2 + C2 * (1 - fm1) * (1 - fm2));            
        factor = 2 * C1 * fm1 * fm2;
      }                    
      else {            
        factor = (fm1 * fm2 + (1 - fm1) * (1 - fm2)) / temperature;  
      }          
      result += factor;
    }
  }
  return result / (2 * num_sites);
}


