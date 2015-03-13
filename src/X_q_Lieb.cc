//This file will generate suseptability vs temperature for the Lieb lattice

using namespace std;
#include "susceptibility.h"
#include <string>
#include <cmath>
#include <fstream>
#include <vector>

int main(int argc, char const *argv[]) {
  
  const double a = 2;

  const double temperature = 1.0 / 36.0;
  
  string modelname = "Lieb";
  double qx = 0;
  double qy = 0;
  // double qx = M_PI / a;
  // double qy = M_PI / a;
  int max_kx = 1000;
  int max_ky = 1000;

  ofstream myfile;  

  string temp = modelname +
    string("_qx_") 
    + to_string(qx) 
    + string("_qy_") 
    + to_string(qy) 
    + string("_k_") 
    + to_string(max_kx)     
    + "_T_" 
    + to_string(temperature) 
    + string(".txt");

  const char *fName = temp.c_str();

  
  myfile.open(fName);  

  double mu_initial = -3;
  double mu_final = 3;
  double mu_step = 0.01; 


  for (double mu = mu_initial; mu <= mu_final; mu += mu_step) {        
    double x00 = susceptibility2(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x01 = susceptibility2(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x11 = susceptibility2(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
    double x12 = susceptibility2(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);
    
    myfile << mu << "\t" << x00 << "\t" << x01 << "\t" << x11 << "\t" << x12 << endl;
  }
  myfile.close();
  return 0;
}

