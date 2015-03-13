//This file will generate suseptability vs temperature for the Lieb lattice

using namespace std;
#include "susceptibility_sc.h"
#include <string>
#include <cmath>
#include <fstream>
#include <vector>

int main(int argc, char const *argv[]) {
  
  const double a = 2;

  // const double temperature = 1.0 / 36.0;
  
  string modelname = "Lieb";
  double qx = 0;
  double qy = 0;
  // double qx = M_PI / a;
  // double qy = M_PI / a;
  int max_kx = 1000;
  int max_ky = 1000;

  ofstream myfile;  

  // string temp = modelname +
  //   string("_sc_qx_") 
  //   + to_string(qx) 
  //   + string("_qy_") 
  //   + to_string(qy) 
  //   + string("_k_") 
  //   + to_string(max_kx)     
  //   + "_T_" 
  //   + to_string(temperature) 
  //   + string(".txt");

  string temp = modelname +
    string("_sc_qx_") 
    + to_string(qx) 
    + string("_qy_") 
    + to_string(qy) 
    + string("_k_") 
    + to_string(max_kx)     
    + "_mu_" 
    + to_string(0) 
    + string(".txt");


  const char *fName = temp.c_str();

  
  myfile.open(fName);  

  // double mu_initial = -3;
  // double mu_final = 3;
  // double mu_step = 0.01; 
  // double mu_initial = 0;
  // double mu_final = 0;
  // double mu_step = 0; 
  double mu = 0;
  double beta_initial = 0.01;
  double beta_final = 1;
  double beta_step = 0.01;



  // for (double mu = mu_initial; mu <= mu_final; mu += mu_step) {        
  for (double beta = beta_initial; beta <= beta_final; beta += beta_step) {
    double temperature = 1.0 / beta;
    double x00 = susceptibility2_sc(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x01 = susceptibility2_sc(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x11 = susceptibility2_sc(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x12 = susceptibility2_sc(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);
    
    // myfile << mu << "\t" << x00 << "\t" << x01 << "\t" << x11 << "\t" << x12 << endl;
    myfile << beta << "\t" << x00 << "\t" << x01 << "\t" << x11 << "\t" << x12 << endl;
  }
  // }
  myfile.close();
  return 0;
}

