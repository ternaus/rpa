#include "common.h"
#include "rho.h"
#include <cmath>
#include <iostream>

using namespace std;

static const int num_bands = 3;

double rho(const int& max_kx, const int& max_ky, const double&a, const double&mu, const double&temperature ) { 
  double result = 0;

  double kx_step = 2 * M_PI / (a * max_kx);
  double ky_step = 2 * M_PI / (a * max_ky);
  double k_min = -M_PI / a;
  double k_max = M_PI / a;

  int num_sites = num_bands * max_kx * max_ky;
  
  for (double tkx = k_min; tkx < k_max; tkx += kx_step) {
    for (double tky = k_min; tky < k_max; tky += ky_step) {
      for (int band_number = 0; band_number < num_bands; band_number++) {
        double band = band_Lieb(band_number, tkx, tky, mu);        
        result += fermi(band, temperature);
      }
    }
  }  
  return 2 * result / num_sites;
}