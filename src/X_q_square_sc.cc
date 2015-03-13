//This file will generate suseptability vs temperature for the Lieb lattice

using namespace std;
#include "susceptibility_sc_square.h"
#include <string>
#include <cmath>
#include <fstream>
#include <vector>

int main(int argc, char const *argv[]) {
  
  const double a = 1;

  const double temperature = 1.0 / 36.0;
  
  string modelname = "square";
  // double qx = M_PI;
  // double qy = M_PI;
  double qx = 0;
  double qy = 0;
  int max_kx = 1000;
  int max_ky = 1000;

  ofstream myfile;  

  string temp = modelname +
    string("_sc_qx_") 
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
    double x = susceptibility_sc_square(temperature, a, mu, max_kx, max_ky, qx, qy);
    myfile << mu << "\t" << x << endl;
  }
  myfile.close();
  return 0;
}

