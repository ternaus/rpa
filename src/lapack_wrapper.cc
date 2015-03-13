#include "lapack_wrapper.h"
#include "mkl.h"
#include <iostream>

//computes all eigenvalues and eigenvectors of a
//real symmetric matrix
//to_diag - vector<vector<double>> to diagonolize

pair<vector<double>, vector<vector<double> > > eigen_symm(const vector<vector<double> > to_diag) {  
  
  int num_columns = to_diag.size();
  int num_rows = num_columns;
  double to_diag_array[num_rows * num_columns];

  //TODO: this should be separate function.
  for (int x = 0; x < num_rows; x++) {
    for (int y = 0; y < num_columns; y++) {
      to_diag_array[x + y * num_rows] = to_diag[x][y];
    }
  }

  //no idea what do these variables mean.
  double wkopt;
  double* work;  //WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) 
          //On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

  double w[num_rows];//array with eigenvalues


  //LWORK is INTEGER
          // The length of the array WORK.  LWORK >= max(1,3*N-1).
          // For optimal efficiency, LWORK >= (NB+2)*N,
          // where NB is the blocksize for DSYTRD returned by ILAENV.

          // If LWORK = -1, then a workspace query is assumed; the routine
          // only calculates the optimal size of the WORK array, returns
          // this value as the first entry of the WORK array, and no error
          // message related to LWORK is issued by XERBLA.
  int lwork;
  int info; //If info == 0 after evaluation, 

  /* Query and allocate the optimal workspace */
  lwork = -1;

  dsyev("V", "Upper", &num_rows, to_diag_array, &num_columns, w, &wkopt, &lwork, &info );

  //TODO This should be changed to the static cast if this will work.
  lwork = (int)wkopt;
  //TODO this should be changed to using new
  work = new double[lwork];

  dsyev("V", "Upper", &num_rows, to_diag_array, &num_columns, w, work, &lwork, &info );

  delete[] work;
  /* Check for convergence */
  if( info > 0 ) {
    cout <<( "The algorithm failed to compute eigenvalues." ) << endl;
    exit(1);
  }
  
 vector<double> eigenvalues(num_columns);

  for (int x = 0; x < num_rows; x++) {
    eigenvalues[x] = w[x];
  }

  vector<vector<double> > eigenvectors(num_columns, vector<double>(num_columns));

  for (int x = 0; x < num_columns; x++) {
    for (int y = 0; y < num_rows; y++) {      
      eigenvectors[y][x] = to_diag_array[x + num_rows * y];
    }
  }

  pair<vector<double>, vector<vector<double> > > result;
  result = make_pair(eigenvalues, eigenvectors);
  return result;
}  