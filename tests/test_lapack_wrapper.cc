#include "lapack_wrapper.h"
#include "gtest/gtest.h"
#include <vector>
#include <cmath>
using namespace std;

TEST(common, diagonalization) {

  vector<vector<double> > victim(5, vector<double>(5));
  victim[0][0] = 1.96;
  victim[0][1] = victim[1][0] = -6.49;
  victim[0][2] = victim[2][0] = -0.47;
  victim[0][3] = victim[3][0] = -7.20;
  victim[0][4] = victim[4][0] = -0.65;
  victim[1][1] = 3.80;
  victim[1][2] = victim[2][1] = -6.39;
  victim[1][3] = victim[3][1] = 1.50;
  victim[1][4] = victim[4][1] = -6.34;
  victim[2][2] = 4.17;
  victim[2][3] = victim[3][2] = -1.51;
  victim[2][4] = victim[4][2] = 2.67;
  victim[3][3] = 5.70;
  victim[3][4] = victim[4][3] = 1.80;
  victim[4][4] = -7.10;

  pair<vector<double>, vector<vector<double> > > result = eigen_symm(victim);  

  //eigenvalues check
  EXPECT_NEAR(-11.0656, result.first[0], 1e-2);  
  EXPECT_NEAR(-6.23, result.first[1], 1e-2);
  EXPECT_NEAR(0.86, result.first[2], 1e-2);
  EXPECT_NEAR(8.87, result.first[3], 1e-2);
  EXPECT_NEAR(16.09, result.first[4], 1e-2);  

  //eigenvectors check
  EXPECT_NEAR(-0.30, result.second[0][0], 1e-2);
  EXPECT_NEAR(-0.51, result.second[0][1], 1e-2);
  EXPECT_NEAR(-0.08, result.second[0][2], 1e-2);
  EXPECT_NEAR(0    , result.second[0][3], 1e-2);
  EXPECT_NEAR(-0.80, result.second[0][4], 1e-2);

  EXPECT_NEAR(-0.61, result.second[1][0], 1e-2);
  EXPECT_NEAR(-0.29, result.second[1][1], 1e-2);
  EXPECT_NEAR(-0.38, result.second[1][2], 1e-2);
  EXPECT_NEAR(-0.45, result.second[1][3], 1e-2);
  EXPECT_NEAR(0.45, result.second[1][4], 1e-2);

  EXPECT_NEAR(0.40, result.second[2][0], 1e-2);
  EXPECT_NEAR(-0.41, result.second[2][1], 1e-2);
  EXPECT_NEAR(-0.66, result.second[2][2], 1e-2);
  EXPECT_NEAR(0.46    , result.second[2][3], 1e-2);
  EXPECT_NEAR(0.17, result.second[2][4], 1e-2);

  EXPECT_NEAR(-0.37, result.second[3][0], 1e-2);
  EXPECT_NEAR(-0.36, result.second[3][1], 1e-2);
  EXPECT_NEAR(0.50, result.second[3][2], 1e-2);
  EXPECT_NEAR(0.62    , result.second[3][3], 1e-2);
  EXPECT_NEAR(0.31, result.second[3][4], 1e-2);

  EXPECT_NEAR(0.49, result.second[4][0], 1e-2);
  EXPECT_NEAR(-0.61, result.second[4][1], 1e-2);
  EXPECT_NEAR(0.40, result.second[4][2], 1e-2);
  EXPECT_NEAR(-0.46, result.second[4][3], 1e-2);
  EXPECT_NEAR(0.16, result.second[4][4], 1e-2);

  //Check normalization

  for (auto x: result.second){
    double temp = 0;
    for (auto y: x) {
      temp += y * y;
    }
    EXPECT_NEAR(1, temp, 1e-4);
  }
} 

TEST(common, diagonalization3x3) {

  vector<vector<double> > victim(3, vector<double>(3));
  victim[0][0] = 0;
  victim[0][1] = victim[1][0] = -2;
  victim[0][2] = victim[2][0] = -2;
  victim[1][1] = 0;
  victim[1][2] = victim[2][1] = 0;
  victim[2][2] = 0;
  
  pair<vector<double>, vector<vector<double> > > result = eigen_symm(victim);  

  //eigenvalues check
  EXPECT_NEAR(-2 * sqrt(2), result.first[0], 1e-2);  
  EXPECT_NEAR(0 * sqrt(2), result.first[1], 1e-2);
  EXPECT_NEAR(2 * sqrt(2), result.first[2], 1e-2);
  //Check normalization

  for (auto x: result.second){
    double temp = 0;
    for (auto y: x) {
      temp += y * y;
    }
    EXPECT_NEAR(1, temp, 1e-4);
  }

  //eigenvectors check
  EXPECT_NEAR(-sqrt(2) / 2, result.second[0][0], 1e-2);
  EXPECT_NEAR(-0.5, result.second[0][1], 1e-2);
  EXPECT_NEAR(-0.5, result.second[0][2], 1e-2);
  
  EXPECT_NEAR(0, result.second[1][0], 1e-2);
  EXPECT_NEAR(-1.0 / sqrt(2), result.second[1][1], 1e-2);
  EXPECT_NEAR(1.0 / sqrt(2), result.second[1][2], 1e-2);
  
  EXPECT_NEAR(sqrt(2) / 2, result.second[2][0], 1e-2);
  EXPECT_NEAR(-0.5, result.second[2][1], 1e-2);
  EXPECT_NEAR(-0.5, result.second[2][2], 1e-2);

}
