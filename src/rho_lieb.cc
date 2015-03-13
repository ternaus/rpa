//This file will generate suseptability vs temperature for the Lieb lattice

using namespace std;
#include "rho.h"
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
  
  const double a = 2;

  //const double temperature = 1.0 / 36.0;
  const double temperature = 0;
  
  string modelname = "Lieb";
  
  int max_kx = 1000;
  int max_ky = 1000;

  ofstream myfile;  

  string temp = modelname +
    string("_rho")     
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
    double rh =  rho(max_kx, max_ky, a, mu, temperature);    
    myfile << mu << "\t" << rh << endl;
  }
  myfile.close();  

  return 0;
}

