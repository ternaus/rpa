#include "common.h"
#include "gtest/gtest.h"
#include <cmath>


TEST(common, fermi) {
  double temperature = 0;
  EXPECT_EQ(0, fermi(1e-5, temperature));
  EXPECT_EQ(1, fermi(-1e-5, temperature));
}

TEST(common, krange) {
  int max_k = 4;
  double a = 2;
  vector<double> krange_x = krange(max_k, a);
  EXPECT_NEAR(-M_PI / a, krange_x[0], 1e-10);
  EXPECT_NEAR(-M_PI / (2 * a), krange_x[1], 1e-10);
  EXPECT_NEAR(0, krange_x[2], 1e-10);
  EXPECT_NEAR(M_PI /(2 * a), krange_x[3], 1e-10);
  EXPECT_EQ(max_k, krange_x.size());
}

//Every column in the diagonalization matrix should be normalized.
TEST(Lieb_gamma, normalization) {
  int num_bands = 3;

  double kx = 4;
  double ky = 6;
  
  for (int column = 0; column < num_bands; column++) {
    double temp = 0;
    for (int row = 0; row < num_bands; row++) {
      temp += pow(gamma(0, row, kx, ky), 2);
    } 
    EXPECT_NEAR(1, temp, 1e-5);
  }
}
//Check that diagonalization matrix is unitary
TEST(Lieb_gamma, unitary) {
  int num_bands = 3;

  double kx = 4;
  double ky = 6;

  for (int x = 0; x < num_bands; x++) {
    for (int y = 0; y < num_bands; y++) {
      double temp = 0;
      for (int k = 0; k < num_bands; k++) {
        temp += gamma(k, y, kx, ky) * gamma(k, x, kx, ky);
      }
      if (x == y) {
        EXPECT_NEAR(1, temp, 1e-5);
      }
      else {
        EXPECT_NEAR(0, temp, 1e-5);
      }
    }
  }
}
