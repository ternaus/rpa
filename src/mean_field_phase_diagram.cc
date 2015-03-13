#include "common.h"
#include "lapack_wrapper.h"
#include "mean_field_phase_diagram.h"
#include "mean_field.h"
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

int main() {  
  double mu_min = -3;
  double mu_max = 3;
  double mu_step = 0.1;
  double u_min = 0;
  double u_max = 8;
  double u_step = 1;
  int iter_max = 10000;
  int kx_max = 1000;
  int ky_max = 1000;
  double temperature = 0.1;
  string modelname = "Lieb";
  ofstream myfile;

  string temp = modelname +  
    string("_k_")  
    + to_string(kx_max)     
    + "_T_" 
    + to_string(temperature) 
    + string(".txt");

  const char *fName = temp.c_str();
  myfile.open(fName); 

  for (double mu = mu_min; mu < mu_max; mu += mu_step) {
    for (double u = u_min; u < u_max; u++) {
      
      double n1_guess = random_double(0, 2);
      double n2_guess = random_double(0, 2);
      double m1_guess = random_double(-1, 1);
      double m2_guess = random_double(-1, 1);
      cout << "mu = " << mu << "  u = " << u << "  n1_guess = " << n1_guess
      << "  n2_guess = " << n2_guess << "  m1_guess = " << m1_guess << "  m2_guess = " << m2_guess << endl;
      int ind = 0;

      for (ind = 0; ind < iter_max; ind++) {
        vector<double> mf = main_mf_loop(n1_guess, 
          n2_guess, 
          m1_guess, 
          m2_guess, 
          kx_max, 
          ky_max,
          u,
          mu,
          temperature);
        double n1_up = mf[0];
        double n2_up = mf[1];
        double n3_up = mf[2];
        double n1_down = mf[3];
        double n2_down = mf[4];
        double n3_down = mf[5];

        assert (abs(n2_up - n3_up) < 1e-10);
        assert (abs(n2_down - n3_down) < 1e-10);

        double n1_guess_new = n1_up + n1_down;
        double n2_guess_new = n2_up + n2_down;
        double m1_guess_new = 0.5 * (n1_up - n1_down);
        double m2_guess_new = 0.5 * (n2_up - n2_down);

        if ((abs(n1_guess - n1_guess_new) < 1e-3) and (abs(n2_guess - n2_guess_new) < 1e-3) and (
        abs(m1_guess - m1_guess_new) < 1e-3) and (abs(m2_guess - m2_guess_new) < 1e-3)) {
          break;
        }

        n1_guess = normalized_n(0.5 * (n1_guess_new + n1_guess));
        n2_guess = normalized_n(0.5 * (n2_guess_new + n2_guess));
        m1_guess = normalized_m(0.5 * (m1_guess_new + m1_guess));
        m2_guess = normalized_m(0.5 * (m2_guess_new + m2_guess));
      }
      myfile << ind << "\t" << mu << "\t" << u << "\t" << n1_guess << "\t" << n2_guess << "\t" << m1_guess << "\t" << m2_guess << endl;
    }
  }
  myfile.close();
  return 0;
}

double normalized_n(const double&n) {
  double result = min(2.0, n);
  result = max(result, 0.0);
  return result;
}

double normalized_m(const double&m) {
  double result = min(1.0, m);
  result = max(result, -1.0);
  return result;
}
