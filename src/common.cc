#include <vector>
#include "common.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

long double fermi(const double& band, const double &temperature) {
  if (temperature > 0) {
    return 1 / (exp(band / temperature) + 1);
  }
  else if (temperature == 0) {
    if (band < 0) {
      return 1;
    } 
    else if (band > 0) {
      return 0;
    }
    else {
      return 0.5;
    }  
  }
  return -1;
}

double gamma(const int&x, const int&y, const double& kx, const double&ky) {
  double result = 0;
  if (x == 0) {
    if (y == 0) {
      result = 0;
    }
    else if (y == 1) {
      result = -cos(ky)/(sqrt(1 + pow(cos(ky), 2)/pow(cos(kx), 2))*cos(kx));
    }
    else if (y == 2) {
      result = pow(1 + pow(cos(ky), 2)/pow(cos(kx), 2), -1.0/2.0);
    }
  }
  else if (x == 1) {
    if (y == 0) {
      result = (-pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) - cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)))/sqrt(pow(-pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) - cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1);
    }
    else if (y == 1) {
      result = cos(kx)/(sqrt(pow(-pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) - cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1)*cos(ky));
    }
    else if (y == 2) {
      result = pow(pow(-pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) - cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1, -1.0/2.0);
    }
  }
  else if (x == 2) {
    if (y == 0) {
      result = (pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) + cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)))/sqrt(pow(pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) + cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1);      
    }
    else if (y == 1) {
      result = cos(kx)/(sqrt(pow(pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) + cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1)*cos(ky));
    }
    else if (y == 2) {
      result = pow(pow(pow(cos(kx), 2)/(sqrt(pow(cos(kx), 2) + pow(cos(ky), 2))*cos(ky)) + cos(ky)/sqrt(pow(cos(kx), 2) + pow(cos(ky), 2)), 2) + pow(cos(kx), 2)/pow(cos(ky), 2) + 1, -1.0/2.0);
    }
  }
  return result;
}

vector<double> krange(const int&max_k, const double&a) {
  vector<double> result;
  double k_step = 2 * M_PI / (a * max_k);
  double k_min = -M_PI / a;
  double k_max = M_PI / a;
  for (double k = k_min; k < k_max; k += k_step) {
    result.push_back(k);
  }
  return result;
}

double random_double(const double&min, const double&max) {
  srand(time(NULL));
  double r = (double)rand() / (double)RAND_MAX;
  return min + r * (max - min);
}

double band_square(const double& tkx, 
  const double& tky, 
  const double& mu) {
  return 2 * (cos(tkx) + cos(tky)) - mu;
}

double band_Lieb(const int& band_number, 
  const double& tkx, 
  const double& tky, 
  const double& mu) {
  if (band_number == 0) {
    return -mu;
  }
  else if (band_number == 1) {
    return -2 * sqrt(pow(cos(tkx), 2) + pow(cos(tky), 2)) - mu;
  }
  else if (band_number == 2) {
    return 2 * sqrt(pow(cos(tkx), 2) + pow(cos(tky), 2)) - mu;
  }
  return 0;
}


