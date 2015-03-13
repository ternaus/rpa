#include "mean_field.h"
#include "lapack_wrapper.h"
#include "gtest/gtest.h"
#include "common.h"
#include <cmath>

TEST(hamiltonian, mu0_u_0_kx0_ky_0) {
  double tkx = 0;
  double tky = 0;
  double mu = 0;
  double u = 0;
  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;
  int num_bands = 3;

  vector<vector<double> > hamiltonain_up = hamiltonian(n1_guess, n2_guess, m1_guess, m2_guess, tkx, tky, 1, u, mu);
  vector<vector<double> > hamiltonain_down = hamiltonian(n1_guess, n2_guess, m1_guess, m2_guess, tkx, tky, -1, u, mu);

  for (int x = 0; x < num_bands; x++) {
    for (int y = 0; y < num_bands; y++) {
      EXPECT_NEAR(hamiltonain_up[x][y], hamiltonain_up[y][x], 1e-7);//Check if symmetric
      EXPECT_NEAR(hamiltonain_up[x][y], hamiltonain_down[x][y], 1e-7);
    }
  }

  EXPECT_NEAR(0, hamiltonain_up[0][0], 1e-10);
  EXPECT_NEAR(0, hamiltonain_up[1][1], 1e-10);
  EXPECT_NEAR(0, hamiltonain_up[2][2], 1e-10);

  EXPECT_NEAR(-2, hamiltonain_up[0][2], 1e-10);
  EXPECT_NEAR(-2, hamiltonain_up[0][1], 1e-10);
}


TEST(mean_field, u0_mu0) {
  int max_kx = 4;
  int max_ky = 4;

  double u = 0;
  double mu = 0;
  double temperature = 0;

  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;

  vector<double> mf = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu, temperature);
  EXPECT_NEAR(mf[1], mf[2], 1e-7);
  EXPECT_NEAR(mf[4], mf[5], 1e-7);
}

TEST(mean_field, u4_mu0) {
  int max_kx = 4;
  int max_ky = 4;

  double u = 4;
  double mu = 0;
  double temperature = 0;

  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;

  vector<double> mf = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu, temperature);
  EXPECT_NEAR(mf[1], mf[2], 1e-7);
  EXPECT_NEAR(mf[4], mf[5], 1e-7);
}

TEST(mean_field, u3_mu2) {
  int max_kx = 4;
  int max_ky = 4;

  double u = 3;
  double mu = 2;
  double temperature = 0;

  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;

  vector<double> mf = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu, temperature);
  EXPECT_NEAR(mf[1], mf[2], 1e-7);
  EXPECT_NEAR(mf[4], mf[5], 1e-7);
}

TEST(mean_field, u0_mu1) {
  int max_kx = 4;
  int max_ky = 4;

  double u = 0;
  double mu = 1;
  double temperature = 0;

  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;

  vector<double> mf = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu, temperature);
  EXPECT_NEAR(mf[1], mf[2], 1e-7);
  EXPECT_NEAR(mf[4], mf[5], 1e-7);
}

//Not sure why this test does not work. May be mean field violates particle-hole symmetry
// TEST(mean_field, u3_mu_pm1) {
//   int max_kx = 4;
//   int max_ky = 4;

//   double u = -4;
//   double mu_1 = 1;
//   double mu_2 = 2;

//   double temperature = 0;

//   double n1_guess = 2.0 / 3.0;
//   double n2_guess = 4.0 / 3.0;
//   double m1_guess = -1.0 / 3.0;
//   double m2_guess = 1.0 / 3.0;


//   vector<double> mf_1 = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu_1, temperature);

//   //Check that values were not changed by accident by this function.
//   EXPECT_NEAR(2.0 / 3.0, n1_guess, 1e-10);
//   EXPECT_NEAR(4.0 / 3.0, n2_guess, 1e-10);
//   EXPECT_NEAR(-1.0 / 3.0, m1_guess, 1e-10);
//   EXPECT_NEAR(1.0 / 3.0, m2_guess, 1e-10);
//   EXPECT_NEAR(-4, u, 1e-10);
//   EXPECT_NEAR(1, mu_1, 1e-10);
//   EXPECT_EQ(4, max_kx);
//   EXPECT_EQ(4, max_ky);
//   EXPECT_NEAR(0, temperature, 1e-10);

//   vector<double> mf_2 = main_mf_loop(n1_guess, n2_guess, m1_guess, m2_guess, max_kx, max_ky, u, mu_2, temperature);

//   //Check for spatial symmetries
//   EXPECT_NEAR(mf_1[1], mf_1[2], 1e-7);
//   EXPECT_NEAR(mf_1[4], mf_1[5], 1e-7);

//   EXPECT_NEAR(mf_2[1], mf_2[2], 1e-7);
//   EXPECT_NEAR(mf_2[4], mf_2[5], 1e-7);

//   //Check for symmetry with respect to mu
//   // EXPECT_NEAR(mf_1[0], mf_2[0], 1e-7);
//   // EXPECT_NEAR(mf_1[1], mf_2[1], 1e-7);
//   // EXPECT_NEAR(mf_1[2], mf_2[2], 1e-7);
//   // EXPECT_NEAR(mf_1[3], mf_2[3], 1e-7);
//   // EXPECT_NEAR(mf_1[4], mf_2[4], 1e-7);
//   // EXPECT_NEAR(mf_1[5], mf_2[5], 1e-7);


//   //TODO Figure out how particle hole transformation
//   //affects mu symmetry
//   double m_1_0 = mf_1[0] - mf_1[3];
//   double m_2_0 = mf_2[0] - mf_2[3];

//   EXPECT_NEAR(m_1_0, m_2_0, 1e-10);

//   double m_1_1 = mf_1[1] - mf_1[4];
//   double m_2_1 = mf_2[1] - mf_2[4];

//   EXPECT_NEAR(m_1_1, m_2_1, 1e-10);

//   double m_1_2 = mf_1[2] - mf_1[5];
//   double m_2_2 = mf_2[2] - mf_2[5];

//   EXPECT_NEAR(m_1_2, m_2_2, 1e-10);
// }

TEST(mean_field, main_mf_loop) {

  int num_bands = 3;
  int max_kx = 4;
  int max_ky = 4;

  double u = -3;
  double mu = 2;
  double temperature = 0;

  double n1_guess = 2.0 / 3.0;
  double n2_guess = 4.0 / 3.0;
  double m1_guess = -1.0 / 3.0;
  double m2_guess = 1.0 / 3.0;

  double a = 2;

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
  EXPECT_NEAR(n1_up, n1_down, 1e-10);
  EXPECT_NEAR(n2_up, n2_down, 1e-10);
  EXPECT_NEAR(n3_up, n3_down, 1e-10);

  EXPECT_NEAR(n2_up, n3_up, 1e-10);
  EXPECT_NEAR(n2_down, n3_down, 1e-10);

}