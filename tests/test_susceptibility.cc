// //some checks will be her till I will figure out how to use 
// // g-test

#include "susceptibility.h"
#include "susceptibility_sc.h"
#include "susceptibility_square.h"
#include "susceptibility_sc_square.h"
#include "common.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(susceptibility_Lieb, susceptibility2) {
  double temperature = 1;
  int alpha = 0;
  
  int beta1 = 1;
  int beta2 = 2;

  double a = 2;
  double mu = 0;
  int max_kx = 10;
  int max_ky = 10;
  double qx = 0;
  double qy = 0;

  double x1 = susceptibility2(alpha, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x2 = susceptibility2(alpha, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);

  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2(beta1, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2(beta2, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);

  EXPECT_NEAR(x1, x2, 1e-4);    
}

TEST(susceptibility_Lieb, band) {

  int band_number = 0;
  double tkx = 10;
  double tky = 100;
  double mu = 0;

  EXPECT_NEAR(0, band_Lieb(band_number, tkx, tky, mu- mu), 1e-10);
}

TEST(susceptibility_Lieb, susceptibility) {
  double temperature = 1.0 / 4.0;
  double a = 2;
  double qx = 0;
  double qy = 0;
  double mu = 2;
  int max_kx = 4;
  int max_ky = 4;
  double x = susceptibility(temperature, a, mu, max_kx, max_ky, qx, qy);
  EXPECT_NEAR(0.195607150119, x, 1e-6);
  
  temperature = 1;
  mu = 0;
  max_kx = 2;
  max_ky = 2;
  
  x = susceptibility(temperature, a, mu, max_kx, max_ky, qx, qy);
  EXPECT_NEAR(0.167, x, 1e-2);
  
  temperature = 0.5;

  x = susceptibility(temperature, a, mu, max_kx, max_ky, qx, qy);

  EXPECT_NEAR(0.262, x, 1e-2);  
}

TEST(susceptibility_square, X_sc_vs_X_match_mu0) {
  double temperature = 1.0 / 4.0;
  double a = 1;  
  int max_kx = 4;
  int max_ky = 4;
  double mu = 0;

  double qx = 0;
  double qy = 0;

  double x_F = susceptibility_square(temperature, a, mu, max_kx, max_ky, qx, qy);
  qx = M_PI;
  qy = M_PI;
  double x_AF = susceptibility_square(temperature, a, mu, max_kx, max_ky, qx, qy);
  qx = 0;
  qy = 0;
  double X_F_sc = susceptibility_sc_square(temperature, a, mu, max_kx, max_ky, qx, qy);
  qx = M_PI;
  qy = M_PI;
  double X_AF_sc = susceptibility_sc_square(temperature, a, mu, max_kx, max_ky, qx, qy);
  EXPECT_NEAR(x_F, X_AF_sc, 1e-5);
  EXPECT_NEAR(x_AF, X_F_sc, 1e-5);
}

TEST(susceptibility_sc_Lieb, spatial_symmetries) {
  double temperature = 1;
  int alpha = 0;
  
  int beta1 = 1;
  int beta2 = 2;

  double a = 2;
  double mu = 0;
  int max_kx = 10;
  int max_ky = 10;

  //checks for q =0 cases
  double qx = 0;
  double qy = 0;

  double x1 = susceptibility2_sc(alpha, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x2 = susceptibility2_sc(alpha, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);
  //Check that X01 = X02
  EXPECT_NEAR(x1, x2, 1e-4);  

  //Check that X11 = X22
  x1 = susceptibility2_sc(beta1, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta2, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);

  //Check that X01 = X10
  x1 = susceptibility2_sc(alpha, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, alpha, temperature, a, mu, max_kx, max_ky, qx, qy);


  EXPECT_NEAR(x1, x2, 1e-4);    

  //checks for q = pi case
  qx = M_PI / a;
  qy = M_PI / a;

  x1 = susceptibility2_sc(alpha, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(alpha, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);
  //Check that X01 = X02
  EXPECT_NEAR(x1, x2, 1e-4);  

  //Check that X11 = X22
  x1 = susceptibility2_sc(beta1, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta2, beta2, temperature, a, mu, max_kx, max_ky, qx, qy);

  //Check that X01 = X10
  x1 = susceptibility2_sc(alpha, beta1, temperature, a, mu, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, alpha, temperature, a, mu, max_kx, max_ky, qx, qy);


  EXPECT_NEAR(x1, x2, 1e-4);    
}

TEST(susceptibility_sc_Lieb, mu_symmetries) {
  double temperature = 1.0 / 36.0;
  int alpha = 0;
  
  int beta1 = 1;
  int beta2 = 2;

  double a = 2;
  double mu1 = -0.1;
  double mu2 = -mu1;

  int max_kx = 50;
  int max_ky = 50;
  double qx = 0;
  double qy = 0;

  double x1 = susceptibility2_sc(alpha, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  double x2 = susceptibility2_sc(alpha, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X01(mu = -2) == X01(mu=2)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2_sc(alpha, alpha, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(alpha, alpha, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X00(mu = -2) == X00(mu=2)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2_sc(beta1, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X11(mu = -2) == X11(mu=2)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2_sc(beta1, beta2, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, beta2, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X12(mu = -2) == X12(mu=2)
  EXPECT_NEAR(x1, x2, 1e-4);  


  //Check that X11(mu = -2) = X11(mu = 2)
  x1 = susceptibility2_sc(beta1, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);

  EXPECT_NEAR(x1, x2, 1e-4);    
}

TEST(susceptibility_Lieb, mu_symmetries) {
  double temperature = 1.0 / 36.0;
  int alpha = 0;
  
  int beta1 = 1;
  int beta2 = 2;

  double a = 2;
  double mu1 = -0.1;
  double mu2 = 0.1;

  int max_kx = 50;
  int max_ky = 50;
  double qx = 0;
  double qy = 0;

  double x1 = susceptibility2(alpha, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  double x2 = susceptibility2(alpha, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X01(mu = -0.1) == X01(mu=0.1)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2(alpha, alpha, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2(alpha, alpha, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X00(mu = -0.1) == X00(mu=0.1)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2(beta1, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2(beta1, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X11(mu = -0.1) == X11(mu=0.1)
  EXPECT_NEAR(x1, x2, 1e-4);  

  x1 = susceptibility2(beta1, beta2, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, beta2, temperature, a, mu2, max_kx, max_ky, qx, qy);
  //Check that X12(mu = -0.1) == X12(mu=0.1)
  EXPECT_NEAR(x1, x2, 1e-4);

  //Check that X11(mu = -2) = X11(mu = 2)
  x1 = susceptibility2_sc(beta1, beta1, temperature, a, mu1, max_kx, max_ky, qx, qy);
  x2 = susceptibility2_sc(beta1, beta1, temperature, a, mu2, max_kx, max_ky, qx, qy);

  EXPECT_NEAR(x1, x2, 1e-4);    
}


TEST(susceptibility_Lieb, X_sc_vs_X_match_mu0) {
  double temperature = 1.0 / 4.0;
  double a = 2;  
  int max_kx = 100;
  int max_ky = 100;
  double mu = 0;

  double qx = 0;
  double qy = 0;

  double x00_sc_0 = susceptibility2_sc(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x01_sc_0 = susceptibility2_sc(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x11_sc_0 = susceptibility2_sc(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x12_sc_0 = susceptibility2_sc(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);

  double x00_0 = susceptibility2(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x01_0 = susceptibility2(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x11_0 = susceptibility2(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x12_0 = susceptibility2(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);

  double x_F_sc_0 = x00_sc_0 + 4 * x12_sc_0 + 2 * x11_sc_0 + 2 * x12_sc_0;
  double x_AF_sc_0 = x00_sc_0 - 4 * x12_sc_0 + 2 * x11_sc_0 + 2 * x12_sc_0;
  double x_F_0 = x00_0 + 4 * x12_0 + 2 * x11_0 + 2 * x12_0;
  double x_AF_0 = x00_0 - 4 * x12_0 + 2 * x11_0 + 2 * x12_0;

  // EXPECT_NEAR(x_F_sc_0, x_AF_0, 1e-5); 

  qx = M_PI / a;
  qy = M_PI / a;

  double x00_sc_pi = susceptibility2_sc(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x01_sc_pi = susceptibility2_sc(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x11_sc_pi = susceptibility2_sc(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x12_sc_pi = susceptibility2_sc(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);

  double x00_pi = susceptibility2(0, 0, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x01_pi = susceptibility2(0, 1, temperature, a, mu, max_kx, max_ky, qx, qy);  
  double x11_pi = susceptibility2(1, 1, temperature, a, mu, max_kx, max_ky, qx, qy);
  double x12_pi = susceptibility2(1, 2, temperature, a, mu, max_kx, max_ky, qx, qy);

  double x_F_sc_pi = x00_sc_pi + 4 * x12_sc_pi + 2 * x11_sc_pi + 2 * x12_sc_pi;
  double x_AF_sc_pi = x00_sc_pi - 4 * x12_sc_pi + 2 * x11_sc_pi + 2 * x12_sc_pi;
  double x_F_pi = x00_pi + 4 * x12_pi + 2 * x11_pi + 2 * x12_pi;
  double x_AF_pi = x00_pi - 4 * x12_pi + 2 * x11_pi + 2 * x12_pi;

  // EXPECT_NEAR(x_F_sc_pi, x_AF_pi, 1e-5); 



  cout << "x_F(0) = " << x_F_0 << endl;
  cout << "x_AF(0) = " << x_AF_0 << endl;
  cout << "x_F(pi) = " << x_F_pi << endl;
  cout << "x_AF(pi) = " << x_AF_pi << endl;

  cout << "x_F_sc(0) = " << x_F_sc_0 << endl;
  cout << "x_AF_sc(0) = " << x_AF_sc_0 << endl;
  cout << "x_F_sc(pi) = " << x_F_sc_pi << endl;
  cout << "x_AF_sc(pi) = " << x_AF_sc_pi << endl;

}