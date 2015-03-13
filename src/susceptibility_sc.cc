using namespace std;

#include "susceptibility_sc.h"
#include "common.h"
#include <cmath>
#include <iostream>


static const int num_bands = 3;

double susceptibility_sc(const double&temperature, 
                    const double&a, 
                    const double &mu, 
                    const int& max_kx, 
                    const int& max_ky, 
                    const double&qx, 
                    const double&qy) {
  double result = 0;
  
  if ((abs(qx) < 1e-10) || (abs(qy) < 1e-10)) {
    result = 2 * susceptibility2_sc(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy) +
    + 2 * susceptibility2_sc(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + susceptibility2_sc(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + 2 * susceptibility2_sc(1, 0, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + 2 * susceptibility2_sc(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);
  }
  else {
    for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
      for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
        result += susceptibility2_sc(band_number1, band_number2, temperature, a, mu, max_kx, max_ky, qx, qy);
      }  
    }
  }
return result;  
}

double susceptibility2_sc(const int& alpha, 
                      const int& beta, 
                      const double& temperature, 
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
      for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
        double band1 = band_Lieb(band_number1, tkx, tky, mu);
        double fm1 = fermi(band1, temperature);


        for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
          double band2 = band_Lieb(band_number2, -tkx -qx, -tky - qy, mu);          
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
          result += factor * gamma(alpha, band_number1, tkx, tky) * gamma(alpha, band_number2, -tkx - qx, - tky - qy) * gamma(beta, band_number2, -tkx - qx, -tky - qy) * gamma(beta, band_number1, tkx, tky);
        }
      }
    }
  }

  return result / (2 * num_sites);
}
