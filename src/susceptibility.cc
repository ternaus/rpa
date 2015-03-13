using namespace std;

#include "susceptibility.h"
#include "common.h"
#include <cmath>
#include <iostream>


static const int num_bands = 3;

double susceptibility(const double&temperature, 
                    const double&a, 
                    const double &mu, 
                    const int& max_kx, 
                    const int& max_ky, 
                    const double&qx, 
                    const double&qy) {
  double result = 0;
  
  if ((abs(qx) < 1e-10) || (abs(qy) < 1e-10)) {
    result = 2 * susceptibility2(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy) +
    + 2 * susceptibility2(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + susceptibility2(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + 2 * susceptibility2(1, 0, temperature, a, mu, max_kx, max_ky, qx, qy) + 
    + 2 * susceptibility2(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);
  }
  else {
    for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
      for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
        result += susceptibility2(band_number1, band_number2, temperature, a, mu, max_kx, max_ky, qx, qy);
      }  
    }
  }
return result;  
}

double susceptibility2(const int& alpha, 
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

  double factor;
  if (max_kx != max_ky) {    
    //Without symmetries
    for (double tkx = k_min; tkx < k_max; tkx += kx_step) {
      for (double tky = k_min; tky < k_max; tky += ky_step) {
        for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
          double band1 = band_Lieb(band_number1, tkx, tky, mu);
          for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
            double band2 = band_Lieb(band_number2, tkx + qx, tky + qy, mu);

            if (abs(band1 - band2) < 1e-10) {
              double fm = fermi(band1, temperature);
              factor = -fm * (1 - fm) / temperature;
            }
            else {
              factor = (fermi(band1, temperature) - fermi(band2, temperature)) / (band1 - band2);
            } 
            result += factor * gamma(alpha, band_number1, tkx, tky) * gamma(alpha, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number1, tkx, tky);
          }
        }
      }
    }
  }
  else {
  //With symmetries
    
  //corner      
  for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
      double band1 = band_Lieb(band_number1, k_min, k_min, mu);
      for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
        double band2 = band_Lieb(band_number2, k_min + qx, k_min + qy, mu);

        if (abs(band1 - band2) < 1e-10) {
          double fm = fermi(band1, temperature);
          factor = -fm * (1 - fm) / temperature;
        }
        else {
          factor = (fermi(band1, temperature) - fermi(band2, temperature)) / (band1 - band2);
        } 
        result += factor * gamma(alpha, band_number1, k_min, k_min) * gamma(alpha, band_number2, k_min + qx, k_min + qy) * gamma(beta, band_number2, k_min + qx, k_min + qy) * gamma(beta, band_number1, k_min, k_min);
      }
    }

  //along the edge  
  for (double tkx = k_min + kx_step; tkx < k_max; tkx += kx_step) {    
    for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
      double band1 = band_Lieb(band_number1, tkx, k_min, mu);
      for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
        double band2 = band_Lieb(band_number2, tkx + qx, k_min + qy, mu);

        if (abs(band1 - band2) < 1e-10) {
          double fm = fermi(band1, temperature);
          factor = -fm * (1 - fm) / temperature;
        }
        else {
          factor = (fermi(band1, temperature) - fermi(band2, temperature)) / (band1 - band2);
        } 
        result += 2 * factor * gamma(alpha, band_number1, tkx, k_min) * gamma(alpha, band_number2, tkx + qx, k_min + qy) * gamma(beta, band_number2, tkx + qx, k_min + qy) * gamma(beta, band_number1, tkx, k_min);
      }
    }
  }
  //In the bulk above diagonal
  for (double tkx = k_min + kx_step; tkx < k_max; tkx += kx_step) {
      for (double tky = tkx + kx_step; tky < k_max; tky += ky_step) {
        for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
          double band1 = band_Lieb(band_number1, tkx, tky, mu);
          for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
            double band2 = band_Lieb(band_number2, tkx + qx, tky + qy, mu);

            if (abs(band1 - band2) < 1e-10) {
              double fm = fermi(band1, temperature);
              factor = -fm * (1 - fm) / temperature;
            }
            else {
              factor = (fermi(band1, temperature) - fermi(band2, temperature)) / (band1 - band2);
            } 
            result += 2 * factor * gamma(alpha, band_number1, tkx, tky) * gamma(alpha, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number1, tkx, tky);
          }
        }
      }
    }    
  //In the bulk along diagonal
    for (double tkx = k_min + kx_step; tkx < k_max; tkx += kx_step) {
      double tky = tkx;
  
        for (int band_number1 = 0; band_number1 < num_bands; band_number1++) {
          double band1 = band_Lieb(band_number1, tkx, tky, mu);
          for (int band_number2 = 0; band_number2 < num_bands; band_number2++) {
            double band2 = band_Lieb(band_number2, tkx + qx, tky + qy, mu);

            if (abs(band1 - band2) < 1e-10) {
              double fm = fermi(band1, temperature);
              factor = -fm * (1 - fm) / temperature;
            }
            else {
              factor = (fermi(band1, temperature) - fermi(band2, temperature)) / (band1 - band2);
            } 
            result += factor * gamma(alpha, band_number1, tkx, tky) * gamma(alpha, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number2, tkx + qx, tky + qy) * gamma(beta, band_number1, tkx, tky);
          }
        }      
    } 
  }
  //In the bulk along diagonal

  return -result / num_sites;
}
