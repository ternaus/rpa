#include "common.h"
#include "lapack_wrapper.h"
#include "mean_field.h"
#include <vector>
#include <cmath>


static const int num_bands = 3;
static const double a = 2;
static const double t = 1;

vector<double> main_mf_loop(const double& n1_guess, 
  const double& n2_guess, 
  const double&m1_guess, 
  const double&m2_guess, 
  const int&max_kx, 
  const int&max_ky, 
  const double&u, 
  const double&mu, 
  const double&temperature) {
  
  double n1_up = 0;
  double n2_up = 0;
  double n3_up = 0;

  double n1_down = 0;
  double n2_down = 0;
  double n3_down = 0;  

  vector<double> kx_range = krange(max_kx, a);
  vector<double> ky_range = krange(max_ky, a);

  int num_momentum = max_kx * max_ky;

  for (double tkx: kx_range) {
    for (double tky: ky_range) {
      vector<vector<double> > hamiltonain_up = hamiltonian(n1_guess, n2_guess, m1_guess, m2_guess, tkx, tky, 1, u, mu);
      vector<vector<double> > hamiltonain_down = hamiltonian(n1_guess, n2_guess, m1_guess, m2_guess, tkx, tky, -1, u, mu);

      pair<vector<double>, vector<vector<double> > > pair_up = eigen_symm(hamiltonain_up);
      pair<vector<double>, vector<vector<double> > > pair_down = eigen_symm(hamiltonain_down);

      //here we get bands_up, bands_down, transform_up, transform_down
      for (int band_number = 0; band_number < num_bands; band_number++) {

        double band_up = pair_up.first[band_number];
        double band_down = pair_down.first[band_number];

        long double fm_up = fermi(band_up, temperature);
        long double fm_down = fermi(band_down, temperature);

        n1_up += pow(pair_up.second[band_number][0], 2) * fm_up;
        n2_up += pow(pair_up.second[band_number][1], 2) * fm_up;
        n3_up += pow(pair_up.second[band_number][2], 2) * fm_up;

        n1_down += pow(pair_down.second[band_number][0], 2) * fm_down;
        n2_down += pow(pair_down.second[band_number][1], 2) * fm_down;
        n3_down += pow(pair_down.second[band_number][2], 2) * fm_down;
      }
    }
  } 

  vector<double> result(6);

  result[0] = n1_up / num_momentum;
  result[1] = n2_up / num_momentum;
  result[2] = n3_up / num_momentum;
  result[3] = n1_down / num_momentum;
  result[4] = n2_down / num_momentum;
  result[5] = n3_down / num_momentum;

  return result;
}

vector<vector<double> > hamiltonian(const double& n1, const double& n2, const double&m1, const double&m2, const double&kx, const double&ky, const int&sigma, const double&u, const double&mu) {  
  vector<vector<double> > result(num_bands, vector<double>(num_bands));  
  result[0][0] = u * (n1 - sigma * m1) - mu - u / 2.0;
  result[1][1] = result[2][2] = u * (n2 - sigma * m2) - mu - u / 2.0;
  result[0][1] = result[1][0]  = -2 * t * cos(a * kx / 2);
  result[0][2] = result[2][0]  = -2 * t * cos(a * ky / 2);
  result[1][2] = result[2][1] = 0;
  return result;
}

