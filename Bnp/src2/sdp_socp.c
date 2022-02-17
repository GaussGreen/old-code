/*
Feb 2001. LTQLM Inc.
Wrapper Author: AdA (Alex d'Aspremont  , aspremon@stanford.edu  , +1 415 864
2962) Base code: Miguel Sousa Lobo  , Lieven Vandenberghe  , and Stephen Boyd.
See code website: www.stanford.edu/~boyd/SOCP.html

WARNING: Matrixes are stored in flat BLAS-LAPACK (FORTRAN) format  , column
major order  , arrays starting at 0  , the function ijtok(i  ,j ,number_of_rows)
makes the conversion.

This solves a quadratic optimisation program under LP constraints using an SOCP
solver. See header for source description.
*/

#include "float.h"
#include "math.h"

#include "sdp_sdplib.h"
#include "sdp_socp.h"
#include "utallhdr.h"

/************************************************************************/
/*				WRAPPER CODE (Phase One and Two) */
/************************************************************************/

// Set up a demo vol surface smoothing problem
Err socp_smoothvol(int **index_set, double *upperbounds, double *lowerbounds,
                   int numpoints, int num_dates, double centering_weight,
                   double **result, int printlevel) {
  int dim = (int)((num_dates * (num_dates + 1)) / 2);
  int dimobj = dim - num_dates;
  double **constraint_matrix = dmatrix(1, 2 * numpoints + 2 * dim, 1, dim);
  double **objective_matrix = dmatrix(1, 2 * dimobj, 1, dim);
  double *constraint_vector = dvector(1, 2 * numpoints + 2 * dim);
  double *resbuf = dvector(1, dim + 2);
  int i, j, k, pos;
  Err err;
  // Upper bound on vol.
  double SOCP_MAXVOL =
      1.0; /// *** Warning: this must be set to an appropriate default value ***

  // Fill the objective matrix
  k = 0;
  for (j = 2; j <= num_dates; j++)
    for (i = 1; i <= j - 1; i++) {
      k++;
      objective_matrix[k][toflat(i, j)] = 1.0;
      objective_matrix[k][toflat(i + 1, j)] = -1.0;
      k++;
      objective_matrix[k][toflat(i, j)] = 1.0;
      objective_matrix[k][toflat(i, j - 1)] = -1.0;
    }
  // Fill the Bid-Ask constraints
  for (i = 1; i <= numpoints; i++) {
    pos = toflat(index_set[1][i], index_set[2][i]);
    constraint_matrix[2 * (i - 1) + 1][pos] = -1.0;
    constraint_vector[2 * (i - 1) + 1] = upperbounds[i];
    constraint_matrix[2 * i][pos] = 1.0;
    constraint_vector[2 * i] = -lowerbounds[i];
  }
  // Fill the upper bounds and positivity constraints
  for (i = 2 * numpoints + 1; i <= 2 * numpoints + dim; i++)
    constraint_matrix[i][i - 2 * numpoints] = 1.0;
  for (i = 2 * numpoints + dim + 1; i <= 2 * numpoints + 2 * dim; i++) {
    constraint_matrix[i][i - (2 * numpoints + dim)] = -1.0;
    constraint_vector[i] = SOCP_MAXVOL;
  }
  // Launch the optimization
  err = socp_smooth_lp(objective_matrix, constraint_matrix, constraint_vector,
                       2 * dimobj, dim, 2 * numpoints + 2 * dim,
                       centering_weight, resbuf, printlevel);
  // return result
  for (j = 1; j <= num_dates; j++)
    for (i = 1; i <= j; i++)
      result[i][j] = resbuf[toflat(i, j)];
  for (j = 1; j <= num_dates; j++)
    for (i = j; i <= num_dates; i++)
      result[i][j] = resbuf[toflat(j, i)];
  // free and return
  free_dmatrix(constraint_matrix, 1, 2 * numpoints, 1, numpoints);
  free_dmatrix(objective_matrix, 1, 2 * dimobj, 1, dim);
  free_dvector(constraint_vector, 1, 2 * numpoints);
  free_dvector(resbuf, 1, dim + 2);
  return err;
}

// Set up a demo Bid_Ask problem
Err socp_bidask(double *upperbounds, double *lowerbounds, int numpoints,
                double centering_weight, double *result, int printlevel) {
  double **constraint_matrix = dmatrix(1, 2 * numpoints, 1, numpoints);
  double **objective_matrix = dmatrix(1, numpoints - 1, 1, numpoints);
  double *constraint_vector = dvector(1, 2 * numpoints);
  int i;
  Err err;

  // Set the Bid-Ask constraints.
  for (i = 1; i <= numpoints; i++) {
    constraint_matrix[2 * (i - 1) + 1][i] = -1.0;
    constraint_vector[2 * (i - 1) + 1] = upperbounds[i];
    constraint_matrix[2 * i][i] = 1.0;
    constraint_vector[2 * i] = -lowerbounds[i];
  }
  // The objective matrix computes the gradient  ,
  fill_gradient_matrix(objective_matrix, numpoints, numpoints - 1);
  // Call the main function
  err = socp_smooth_lp(objective_matrix, constraint_matrix, constraint_vector,
                       numpoints - 1, numpoints, 2 * numpoints,
                       centering_weight, result, printlevel);
  // free and return
  free_dmatrix(objective_matrix, 1, numpoints - 1, 1, numpoints);
  free_dmatrix(constraint_matrix, 1, 2 * numpoints, 1, numpoints);
  free_dvector(constraint_vector, 1, 2 * numpoints);
  return err;
}

//	Main function for vol surface smoothing  , solves the following program:
//
//		Minimize ||Ax||^2
//		subject to M.x+b>=0 where M is a (dim  ,numconst) matrix
//
//	"center" is the l-inf center of the inequalities and beta>0 avoids
//	possible degenerate solutions (parallel lines).
//	"center" is computed by the phase one problem.
Err socp_smooth_lp(double **objective_matrix, double **constraints_matrix,
                   double *constraints_value, int nobjrows, int dim,
                   int numconst, double alpha, double *result, int printlevel) {
  Err err = NULL;
  double *centersol = dvector(0, dim + 2);
  double *dualsol = dvector(0, nobjrows + numconst + 3);
  double defaut_alpha = .5;
  int i;

  // Check alpha  , if incorrect input  , use default (this is mostly harmless).
  if ((alpha < 0.0) || (alpha > 1.0)) {
    alpha = defaut_alpha;
    if (printlevel >= 2)
      smessage("SdpSOCP warning: Centering should be in [0  ,1]  , using "
               "default value... ");
  }
  // Make some tests
  if (numconst < 2)
    err = serror("SdpSOCP error: Not enough constraints ....");
  if (dim < 2)
    err = serror("SdpSOCP error: Problem is scalar.... ");
  // Start phase one and compute the center solution.
  if (err == NULL)
    err = phase_one(objective_matrix, constraints_matrix, constraints_value,
                    nobjrows, dim, numconst, centersol, dualsol, printlevel);
  for (i = 1; i <= dim; i++)
    result[i] = alpha * centersol[i - 1];
  // Test for infeasibility
  if (centersol[dim + 1] >= 0)
    err = serror("SdpSOCP error: Problem is infeasible.... ");
  // Start phase two if phase one OK and feasible  , compute the smooth
  // solution.
  if (err == NULL)
    err = phase_two(objective_matrix, constraints_matrix, constraints_value,
                    nobjrows, dim, numconst, centersol, dualsol, printlevel);
  for (i = 1; i <= dim; i++)
    result[i] += (1 - alpha) * centersol[i - 1];
  // Free and return
  free_dvector(centersol, 0, dim + 2);
  free_dvector(dualsol, 0, nobjrows + numconst + 3);
  return err;
}

// converts between (i  ,j) matrix storage and flat LAPACK matrix storage  ,
// column major order.
int ijtok(int i, int j, int numrow) { return (i - 1) + (j - 1) * numrow; }

// Fills the matrix that turns the dim vector x into the dim-1 vector dx.
Err fill_gradient_matrix(double **gradmat, int dimcol, int rowstofill) {
  int i, j;
  for (i = 1; i <= rowstofill; i++)
    for (j = 1; j <= dimcol; j++) {
      if (i == j)
        gradmat[i][j] = -1.0;
      if (j == i + 1)
        gradmat[i][j] = 1.0;
    }
  return NULL;
}

// Converts between upper triangular storage index and flat storage.
int toflat(int i, int j) { return (int)((j * (j - 1)) / 2) + i; }

// Looks for the center of the inequality by solving phase one.
// Wrapper for the core phase one function "phase_one_go"
Err phase_one(double **objective_matrix, double **constraints_matrix,
              double *constraints_value, int nobjrows, int dim, int numconst,
              double *centersol, double *dualsol, int printlevel) {
  // Algo params
  int iter = 80;
  int out_mode = 1;
  int *vecdimconstin = ivector(0, numconst + 2);
  int ndbl, nint, mhist, nhist;
  // Working vars
  int i;
  Err err;

  // Defines the problem size
  vecdimconstin[0] = nobjrows + 1;
  for (i = 2; i <= numconst + 2; i++)
    vecdimconstin[i - 1] = 1;
  // Set up the working arrays
  socp_getwork(numconst + 2, vecdimconstin, dim + 2, iter, out_mode, &mhist,
               &nhist, &ndbl, &nint);
  // Calls the core procedure
  err = phase_one_go(objective_matrix, constraints_matrix, constraints_value,
                     nobjrows, dim, numconst, centersol, dualsol, vecdimconstin,
                     ndbl, nint, mhist, nhist, iter, out_mode, printlevel);
  // Free and return
  free_ivector(vecdimconstin, 0, numconst + 2);
  return err;
}

// Looks for the smoothest solution inside the inequalitites.
// Wrapper for the core phase one function "phase_two_go"
Err phase_two(double **objective_matrix, double **constraints_matrix,
              double *constraints_value, int nobjrows, int dim, int numconst,
              double *centersol, double *dualsol, int printlevel) {
  // Algo params
  int iter = 50;
  int out_mode = 1;
  int *vecdimconstin = ivector(0, numconst + 2);
  int ndbl, nint, mhist, nhist;
  // Working vars
  int i;
  Err err;

  // Defines the problem size
  vecdimconstin[0] = nobjrows + 1;
  for (i = 2; i <= numconst + 2; i++)
    vecdimconstin[i - 1] = 1;
  // Set up the working arrays
  socp_getwork(numconst + 2, vecdimconstin, dim + 1, iter, out_mode, &mhist,
               &nhist, &ndbl, &nint);
  // Calls the core procedure
  err = phase_two_go(objective_matrix, constraints_matrix, constraints_value,
                     nobjrows, dim, numconst, centersol, dualsol, vecdimconstin,
                     ndbl, nint, mhist, nhist, iter, out_mode, printlevel);
  // Free and return
  free_ivector(vecdimconstin, 0, numconst + 2);
  return err;
}

// Looks for the smoothest solution inside the inequalitites.
Err phase_one_go(double **objective_matrix, double **constraints_matrix,
                 double *constraints_value, int nobjrows, int dim, int numconst,
                 double *centersol, double *dualsol, int *vecdimconstin,
                 int ndbl, int nint, int mhist, int nhist, int iter,
                 int out_mode, int printlevel) {
  // Algorithm params
  double abs_tol = -1.0;
  double rel_tol = 1e-6;
  double target = 0.0;
  double nu = 6.0;
  double big_M = 1000;
  double feasibility_margin = 1.5;
  // Problem variables
  double *fvecin = dvector(0, dim + 2);
  double *amatin = dvector(0, (dim + 2) * (numconst + nobjrows + 2));
  double *bvecin = dvector(0, numconst + nobjrows + 2);
  double *xvar = dvector(0, dim + 2);
  double *zvar = dvector(0, numconst + nobjrows + 2);
  // Working variables
  double *dblwork = dvector(0, ndbl);
  int *intwork = ivector(0, nint);
  double *hist = dvector(0, iter + 1);
  int info;
  int i, j;
  double bufa = 0.0;
  int output = 0;
  Err err = NULL;

  // Start filling the problem data.... the big A matrix first...
  // The objective matrix
  for (i = 1; i <= nobjrows; i++)
    for (j = 1; j <= dim; j++)
      amatin[ijtok(i, j, numconst + nobjrows + 2)] = objective_matrix[i][j];
  // the minimizer for the norm and the feasibility variable
  amatin[ijtok(nobjrows + 1, dim + 1, numconst + nobjrows + 2)] = 1.0;
  amatin[ijtok(nobjrows + 1, dim + 2, numconst + nobjrows + 2)] = 1.0;
  // then the LP constraints
  for (i = nobjrows + 2; i <= nobjrows + 1 + numconst; i++)
    for (j = 1; j <= dim; j++) {
      amatin[ijtok(i, j, numconst + nobjrows + 2)] =
          constraints_matrix[i - (nobjrows + 1)][j];
      // the big M constraint.
      amatin[ijtok(numconst + nobjrows + 2, j, numconst + nobjrows + 2)] -=
          constraints_matrix[i - (nobjrows + 1)][j];
    }
  amatin[ijtok(numconst + nobjrows + 2, dim + 1, numconst + nobjrows + 2)] =
      -1.0;
  // The feasibility variable column.
  for (i = nobjrows + 2; i <= nobjrows + 1 + numconst; i++)
    amatin[ijtok(i, dim + 2, numconst + nobjrows + 2)] = 1.0;
  // Big_M again.
  amatin[ijtok(numconst + nobjrows + 2, dim + 2, numconst + nobjrows + 2)] =
      -(double)(numconst + 1);
  // The constraints vector
  for (i = nobjrows + 2; i <= nobjrows + 1 + numconst; i++) {
    bvecin[i - 1] = constraints_value[i - (nobjrows + 1)];
    // the big M constraint
    bvecin[numconst + nobjrows + 1] -= constraints_value[i - (nobjrows + 1)];
  }
  bvecin[numconst + nobjrows + 1] += big_M;
  // The objective vector
  fvecin[nobjrows + 2] = 1.0;
  // The primal variables
  for (i = 1; i <= numconst; i++)
    if (bufa < fabs(constraints_value[i]))
      bufa = fabs(constraints_value[i]);
  xvar[dim + 1] = feasibility_margin * bufa + 1;
  // The dual variables ...
  err = feasible_dual_vecs_one(fvecin, amatin, nobjrows, dim, numconst, zvar);
  // Check feasibility
  info = 0;
  info = check_feasibility(amatin, fvecin, nobjrows, dim, numconst, zvar,
                           vecdimconstin, bvecin, xvar);
  if (info > 0)
    err = serror(
        "SdpSOCP error: No strictly feasible solution has been found... ");
  // Calls the SOCP procedure
  if (err == NULL)
    output = socp(numconst + 2, vecdimconstin, dim + 2, fvecin, amatin, bvecin,
                  xvar, zvar, abs_tol, rel_tol, target, &iter, nu, &info,
                  out_mode, hist, dblwork, intwork);
  // Process the output
  for (i = 1; i <= nobjrows + numconst + 3; i++)
    dualsol[i - 1] = zvar[i - 1];
  for (i = 1; i <= dim + 2; i++)
    centersol[i - 1] = xvar[i - 1];
  if (output != 0)
    err = serror("SdpSOCP error: Error in socp core... ");
  // free and return
  free_dvector(fvecin, 0, dim + 2);
  free_dvector(amatin, 0, (dim + 2) * (numconst + nobjrows + 2));
  free_dvector(bvecin, 0, numconst + nobjrows + 2);
  free_dvector(xvar, 0, dim + 2);
  free_dvector(zvar, 0, numconst + nobjrows + 2);
  free_dvector(dblwork, 0, ndbl);
  free_ivector(intwork, 0, nint);
  free_dvector(hist, 0, iter + 1);
  return err;
}

// Looks for the smoothest solution inside the inequalitites.
Err phase_two_go(double **objective_matrix, double **constraints_matrix,
                 double *constraints_value, int nobjrows, int dim, int numconst,
                 double *centersol, double *dualsol, int *vecdimconstin,
                 int ndbl, int nint, int mhist, int nhist, int iter,
                 int out_mode, int printlevel) {
  // Algorithm params
  double abs_tol = -1.0;
  double rel_tol = 1e-6;
  double target = 0.0;
  double nu = 6.0;
  double big_M = 10000;
  double feasibility_margin = 1.5;
  // Problem variables
  double *fvecin = dvector(0, dim + 1);
  double *amatin = dvector(0, (dim + 1) * (numconst + nobjrows + 2));
  double *bvecin = dvector(0, numconst + nobjrows + 2);
  double *xvar = dvector(0, dim + 1);
  double *zvar = dvector(0, numconst + nobjrows + 2);
  // Working variables
  double *dblwork = dvector(0, ndbl);
  int *intwork = ivector(0, nint);
  double *hist = dvector(0, iter + 1);
  int info;
  int i, j;
  double bufa = 0.0;
  int output = 0;
  Err err = NULL;

  // Start filling the problem data.... the big A matrix first...
  // The objective matrix
  for (i = 1; i <= nobjrows; i++)
    for (j = 1; j <= dim; j++)
      amatin[ijtok(i, j, numconst + nobjrows + 2)] = objective_matrix[i][j];
  // the minimizer for the gradient norm and the feasibility variable
  amatin[ijtok(nobjrows + 1, dim + 1, numconst + nobjrows + 2)] = 1.0;
  // then the LP constraints
  for (i = nobjrows + 2; i <= nobjrows + 1 + numconst; i++)
    for (j = 1; j <= dim; j++) {
      amatin[ijtok(i, j, numconst + nobjrows + 2)] =
          constraints_matrix[i - (nobjrows + 1)][j];
      // the big M constraint.
      amatin[ijtok(numconst + nobjrows + 2, j, numconst + nobjrows + 2)] -=
          constraints_matrix[i - (nobjrows + 1)][j];
    }
  amatin[ijtok(numconst + dim + 1, dim + 1, numconst + dim + 1)] = -1.0;
  // The constraints vector
  for (i = nobjrows + 2; i <= nobjrows + 1 + numconst; i++) {
    bvecin[i - 1] = constraints_value[i - (nobjrows + 1)];
    // the big M constraint
    bvecin[numconst + nobjrows + 1] -= constraints_value[i - (nobjrows + 1)];
  }
  bvecin[numconst + nobjrows + 1] += big_M;
  // The objective vector
  fvecin[dim + 1] = 1.0;
  // The dual variables ...
  err =
      feasible_dual_vecs_two(fvecin, amatin, nobjrows, dim, numconst, dualsol);
  // Calls the SOCP procedure
  if (err == NULL)
    output = socp(numconst + 2, vecdimconstin, dim + 1, fvecin, amatin, bvecin,
                  centersol, dualsol, abs_tol, rel_tol, target, &iter, nu,
                  &info, out_mode, hist, dblwork, intwork);
  if (output != 0)
    err = serror("SdpSOCP error: Error in socp core... ");
  // free and return
  free_dvector(fvecin, 0, dim + 1);
  free_dvector(amatin, 0, (dim + 1) * (numconst + nobjrows + 2));
  free_dvector(bvecin, 0, numconst + nobjrows + 2);
  free_dvector(xvar, 0, dim + 1);
  free_dvector(zvar, 0, numconst + nobjrows + 2);
  free_dvector(dblwork, 0, ndbl);
  free_ivector(intwork, 0, nint);
  free_dvector(hist, 0, iter + 1);
  return err;
}

// Euclidian norm of a vector
double normvec(double *vec, int ilow, int ihigh) {
  int i;
  double buf = 0;
  for (i = ilow; i <= ihigh; i++)
    buf += vec[i] * vec[i];
  return sqrt(buf);
}

// square
double ddsqr(double a) { return a * a; }

// square
int isqr(int a) { return a * a; }

// Computes a dual feasible starting point for phase one. z=inv(A)(f-c2)  , ....
Err feasible_dual_vecs_one(double *fvecin, double *amatin, int nobjrows,
                           int dim, int numconst, double *zvar) {
  int i, j;
  double *abuf = dvector(0, (dim + 2) * (numconst + nobjrows + 2));
  double *bvecbuf = dvector(0, numconst + nobjrows + 1);
  int mp = dim + 2, np = numconst + nobjrows + 1, nrhsp = 1, ldap = dim + 2,
      ldbp = numconst + nobjrows + 1, infop;
  double *sp = dvector(0, dim + 2);
  double rcond = -1.0;
  int dimlwork = 2 * (dim + 2 + max(2 * (dim + 2), numconst + nobjrows + 1));
  double *dwork = dvector(0, dimlwork);
  double bufa = 0.0;
  double feasibility_margin = 1.5; // Must be greater than 1...

  // Copy A-> At in abuf and f into bvec
  for (i = 1; i <= numconst + nobjrows + 2; i++)
    for (j = 1; j <= dim + 2; j++)
      abuf[ijtok(j, i, dim + 2)] = amatin[ijtok(i, j, numconst + nobjrows + 2)];
  for (i = 0; i <= dim + 1; i++)
    bvecbuf[i] = fvecin[i];
  // Calls the LAPACK function.
  dgels_("N", &mp, &np, &nrhsp, abuf, &ldap, bvecbuf, &ldbp, dwork, &dimlwork,
         &infop);
  // Process the result
  bufa = normvec(bvecbuf, 0, numconst + nobjrows);
  for (i = nobjrows + 1; i <= numconst + nobjrows + 1; i++)
    if (bvecbuf[i - 1] < -bufa)
      bufa = -bvecbuf[i - 1];
  for (i = nobjrows + 1; i <= numconst + nobjrows + 1; i++)
    bvecbuf[i - 1] += feasibility_margin * bufa + 1;
  // return the result
  for (i = 1; i <= numconst + nobjrows + 1; i++)
    zvar[i - 1] = bvecbuf[i - 1];
  zvar[numconst + nobjrows + 1] = feasibility_margin * bufa + 1;
  // free and return
  free_dvector(abuf, 0, (dim + 2) * (numconst + nobjrows + 2));
  free_dvector(bvecbuf, 0, numconst + nobjrows + 1);
  free_dvector(sp, 0, dim + 2);
  free_dvector(dwork, 0, dimlwork);
  return NULL;
}

// Computes a dual feasible starting point for phase one. z=inv(A)(f-c2)  , ....
Err feasible_dual_vecs_two(double *fvecin, double *amatin, int nobjrows,
                           int dim, int numconst, double *zvar) {
  int i, j;
  double *abuf = dvector(0, (dim + 1) * (numconst + nobjrows + 2));
  double *bvecbuf = dvector(0, numconst + nobjrows + 1);
  int mp = dim + 1, np = numconst + nobjrows + 1, nrhsp = 1, ldap = dim + 1,
      ldbp = numconst + nobjrows + 1, infop;
  double *sp = dvector(0, dim + 1);
  double rcond = -1.0;
  int dimlwork =
      2 * (dim + 2 +
           max(2 * (dim + 2),
               numconst + nobjrows + 1)); // This is twice the minimum size in
                                          // dgelss but improves speed.
  double *dwork = dvector(0, dimlwork);
  double bufa = 0.0;
  double feasibility_margin = 1.5; // Must be greater than 1...

  // Copy A-> At in abuf and f into bvec
  for (i = 1; i <= numconst + dim + 1; i++)
    for (j = 1; j <= dim + 1; j++)
      abuf[ijtok(j, i, dim + 1)] = amatin[ijtok(i, j, numconst + nobjrows + 2)];
  for (i = 0; i <= dim; i++)
    bvecbuf[i] = fvecin[i];
  // Calls the LAPACK function.
  dgels_("N", &mp, &np, &nrhsp, abuf, &ldap, bvecbuf, &ldbp, dwork, &dimlwork,
         &infop);
  // Process the result
  bufa = normvec(bvecbuf, 0, numconst + nobjrows);
  for (i = nobjrows + 1; i <= numconst + nobjrows + 1; i++)
    if (bvecbuf[i - 1] < -bufa)
      bufa = -bvecbuf[i - 1];
  for (i = nobjrows + 1; i <= numconst + nobjrows + 1; i++)
    bvecbuf[i - 1] += feasibility_margin * bufa + 1;
  // return the result
  for (i = 1; i <= numconst + nobjrows + 1; i++)
    zvar[i - 1] = bvecbuf[i - 1];
  zvar[numconst + nobjrows + 1] = feasibility_margin * bufa + 1;
  // free and return
  free_dvector(abuf, 0, (dim + 1) * (numconst + nobjrows + 2));
  free_dvector(bvecbuf, 0, numconst + nobjrows + 1);
  free_dvector(sp, 0, dim + 1);
  free_dvector(dwork, 0, dimlwork);
  return NULL;
}

// Checks for strict feasibility... (Possible trouble if initial point is to
// close to infeasible).
int check_feasibility(double *amatin, double *fvecin, int nobjrows, int dim,
                      int numconst, double *zvar, int *vecdimconstin,
                      double *bvecin, double *xvar) {
  int i, j;
  double *abuf = dvector(0, (dim + 2) * (numconst + nobjrows + 2));
  double *vecbuf1 = dvector(0, dim + 2);
  double *vecbuf2 = dvector(0, numconst + nobjrows + 2);
  double alpha = 1.0, beta = -1.0, bufa = 0.0, bufb = 0.0;
  int int1 = 1, mp = dim + 2, np = numconst + nobjrows + 2, ldap = dim + 2;
  int info = 0; // for DEBUG
  double minmargin = 1e-8;

  // Dual feasibility...
  for (i = 1; i <= numconst + nobjrows + 2; i++)
    for (j = 1; j <= dim + 2; j++)
      abuf[ijtok(j, i, dim + 2)] = amatin[ijtok(i, j, numconst + nobjrows + 2)];
  for (i = 1; i <= dim + 2; i++)
    vecbuf1[i - 1] = fvecin[i - 1];
  bufa = normvec(zvar, 0, nobjrows - 1);
  if (zvar[nobjrows] - bufa < minmargin)
    info = 1;
  if (zvar[nobjrows + numconst + 1] < minmargin)
    info = 2;
  for (i = nobjrows + 1; i <= numconst + nobjrows + 2; i++)
    if (zvar[i - 1] < minmargin)
      info = 3;
  dgemv_("N", &mp, &np, &alpha, abuf, &ldap, zvar, &int1, &beta, vecbuf1,
         &int1);
  bufa = normvec(vecbuf1, 0, dim + 1);
  if ((bufa / (dim + 1)) > minmargin)
    info = 4;
  // Primal feasibility ...
  beta = 0.0;
  alpha = 1.0;
  mp = nobjrows + numconst + 2;
  np = dim + 2;
  ldap = nobjrows + numconst + 2;
  dgemv_("N", &mp, &np, &alpha, amatin, &ldap, xvar, &int1, &beta, vecbuf2,
         &int1);
  bufa = normvec(vecbuf2, 0, nobjrows - 1);
  if (vecbuf2[nobjrows] - bufa < minmargin)
    info = 5;
  for (i = nobjrows + 2; i <= nobjrows + numconst + 1; i++)
    if ((vecbuf2[i - 1] + bvecin[i - 1]) < minmargin)
      info = 6;
  if ((vecbuf2[nobjrows + numconst + 1] + bvecin[nobjrows + numconst + 1]) <
      minmargin)
    info = 7;
  // free and return
  free_dvector(abuf, 0, (dim + 2) * (numconst + nobjrows + 2));
  free_dvector(vecbuf1, 0, dim + 2);
  free_dvector(vecbuf2, 0, numconst + nobjrows + 2);
  return info;
}

/************************************************************************/
/*				ORIGINAL CODE
*
/************************************************************************/
/*
 *   socp.c
 *
 *   Second-Order Cone Programming
 *   C implementation
 *
 *   mlobo@isl.stanford.edu -- 96/97
 */
void dzero(int n, double *p)
/*
 *  set vector of n doubles to zero
 */
{
  /* faster version:
    n*=sizeof(double)/sizeof(char);
    memset(p  ,0  ,n);
  */
  /* alternative: (using BLAS to avoid need for memory.h) */
  int int1 = 1;
  double double0 = 0.0;
  dscal_(&n, &double0, p, &int1);
}

double dsumdiv(int n, double *x, double *y)
/*
 * returns sum of elementwise division of x by y
 *   in matlab notation: sum(x./y)
 */
{
  int i;
  double sumdiv = 0.0;
  for (i = 0; i < n; ++i)
    sumdiv += x[i] / y[i];
  return sumdiv;
}

void dupge(int n, double *A)
/*
 *  Convert symmetric matrix from upper storage to full storage
 */
{
  int i, j, k;
  int int1 = 1;

  for (i = n - 1, j = 1, k = n; i > 0; --i, j += n + 1, k += n + 1)
    dcopy_(&i, A + k, &n, A + j, &int1);
}

void dgapdev(int m, int L, int *N, double *u, double *z, double *pgap,
             double *pdev)
/*
 *  note: this fct currently not used
 *
 *  compute duality gap and deviation from centrality
 *
 *  output:
 *   *pgap = duality gap
 *   *pdev = deviation from centrality
 */
{
  int i, j, k;
  double fu, fz;

  int int1 = 1;

  /* gap = u'*z; */
  *pgap = ddot_(&m, u, &int1, z, &int1);

  /* dev = -sum([log(SF'*(u.^2));log(SF'*(z.^2))]) + 2*L*(log(gap)-log(L)); */
  *pdev = 2 * L * (log(*pgap) - log(L));
  for (i = 0, k = 0; i < L; ++i) {
    for (j = 0, fu = 0.0, fz = 0.0; j < N[i] - 1; ++j, ++k) {
      fu -= ddsqr(u[k]);
      fz -= ddsqr(z[k]);
    }
    fu += ddsqr(u[k]);
    fz += ddsqr(z[k]);
    ++k;
    *pdev -= log(fu) + log(fz);
  }
}

void dgrad(
    /* input args. */
    double w, int L,
    double c1, /* pre-computed constants  , dependent on: u  , z  , du  , dz */
    double c2, double c3, double *d1, double *d2, double *d3, double *e1,
    double *e2, double *e3, double p, /* du scaling */
    double q,                         /* dz scaling */
    /* output args. */
    double *pgp, /* return gradient */
    double *pgq,
    double *t1, /* intermediate values that will be re-used in dhess() */
    double *t2, double *t3, double *t4)
/*
 *  compute gradient of potential fct wrt. p and q
 *   at u+p*du  , z+q*dz
 *
 *  output:
 *   *pgp = gradient wrt. p
 *   *pgq = gradient wrt. q
 */
{
  double ptwo = 2 * p;
  double psqr = ddsqr(p);
  double qtwo = 2 * q;
  double qsqr = ddsqr(q);

  int int1 = 1;

  /* t1 = d2 + p*d3 */
  dcopy_(&L, d2, &int1, t1, &int1);
  daxpy_(&L, &p, d3, &int1, t1, &int1);

  /* t2 = d1 + 2*p*d2 + p*p*d3 */
  dcopy_(&L, d1, &int1, t2, &int1);
  daxpy_(&L, &ptwo, d2, &int1, t2, &int1);
  daxpy_(&L, &psqr, d3, &int1, t2, &int1);

  /* t3 = e2 + q*e3 */
  dcopy_(&L, e2, &int1, t3, &int1);
  daxpy_(&L, &q, e3, &int1, t3, &int1);

  /* t4 = e1 + 2*q*e2 + q*q*e3 */
  dcopy_(&L, e1, &int1, t4, &int1);
  daxpy_(&L, &qtwo, e2, &int1, t4, &int1);
  daxpy_(&L, &qsqr, e3, &int1, t4, &int1);

  *pgp = w * c2 / (c1 + p * c2 + q * c3) - 2 * dsumdiv(L, t1, t2);
  *pgq = w * c3 / (c1 + p * c2 + q * c3) - 2 * dsumdiv(L, t3, t4);
}

void dhess(
    /* input args. */
    int L, double *d3,      /* pre-computed constants */
    double *e3, double *t1, /* intermediate values from dgrad() */
    double *t2, double *t3, double *t4,
    /* output args. */
    double *php, double *phq)
/*
 *  compute hessian of primal barrier wrt. p and q
 *   at u+p*du  , z+q*dz
 *
 *  output:
 *   *php = hessian wrt. p
 *   *phq = hessian wrt. q
 */
{
  int i;

  for (i = 0, *php = 0.0, *phq = 0.0; i < L; ++i) {
    *php += (d3[i] * t2[i] - 2 * ddsqr(t1[i])) / ddsqr(t2[i]);
    *phq += (e3[i] * t4[i] - 2 * ddsqr(t3[i])) / ddsqr(t4[i]);
  }
  *php *= -2;
  *phq *= -2;
}

void socp_getwork(
    /* input args.: */
    int L, int *N, int n, int max_iter, int out_mode,
    /* output args.: */
    int *mhist, int *nhist, int *ndbl, int *nint)
/*
 *  returns workspace size in number of doubles  , pointers and ints
 *  required by socp()
 *  (use a macro in socp.h instead?)
 *
 *  input arguments:
 *   use the same values as for socp()
 *
 *  output arguments:
 *   *mhist-- number of lines in *hist (elements are doubles)
 *   *nhist-- number of columns in *hist
 *   *ndbl -- number of doubles in *dblwork
 *   *nptr -- number of pointers in *ptrwork
 *   *nint -- number of integers in *intwork
 */
{
  int i;
  int m;

  /* compute m */
  for (i = 0, m = 0; i < L; m += N[i++])
    ;

  *mhist = out_mode;     /* 0: none;  1: gap only;  2: gap and deviation */
  *nhist = max_iter + 1; /* for initial point and after each iteration */

  *ndbl = 7 * m + 2 * n + 10 * n + 2 * n * n + 11 * L;
  *nint = n;
}

int socp(int L, int *N, int n,

         double *f, double *A, double *b,

         double *x, double *z,

         double abs_tol, double rel_tol, double target, int *iter,

         double Nu,

         int *info, int out_mode, double *hist,

         double *dblwork, int *intwork)
/*
 *  solve second order cone problem
 *
 */
{
  int m;
  double w;
  double gap;
  double dev;

  int iter_out;
  double lbound;
  double *u, *fc, *gup, *gu; /* u[m]  , fc[L]  , gup[m]  , gu[m] */
  double *du, *dz;           /* du[m]  , dz[m] */
  double *gx, *dx;           /* gx[n]  , dx[n] */
  double *Hx;                /* Hx[n*(n+1)] */
  double s;

  /* to record reason for exiting main loop */
  int xabs = 0;
  int xrel = 0;
  int xtarget = 0;
  int xiter = 0;

  /* for plane search */
  double c1, c2, c3;
  double *d1, *d2, *d3, *e1, *e2, *e3; /* [L] */
  double p, q, gp, gq, hp, hq, dp, dq;
  double *t1, *t2, *t3, *t4; /* [L] */
  double lambda2;
  int iter_plane;

  /* for line search */
  double p0, q0;
  double alpha, galpha;
  double *ua, *za; /* ua[m]  , za[m] */
  double fu, fz;
  int out_barrier;

  /* general */
  int i, j, k;

  /* for linear system solver */
  /* for both dposvx_() and dgelss_() */
  int use_ls = 0; /* will be changed to 1 if dgelss_ is to be used */
  int la_info;    /* exit info from dposvx_ and dgelss_ */
  double *Hwork;  /* Hwork[SQR(n)] (dgelss_ only needs n) */
  double rcond;
  /* for dposvx_() */
  int equed;
  double scale, ferr, berr;
  /* for dgelss_() */
  int rank;
  int lwork = 10 * n; /* size of workspace  , minimum is 5*n */

  /* constant variables for BLAS */
  int int1 = 1;
  double double1 = 1.0, double0 = 0.0, double_1 = -1.0;

  /* compute m: size of dual variable */
  for (i = 0, m = 0; i < L; m += N[i++])
    ;

  /*
   * organize workspace
   *
   * total space needed =
   *   doubles:   7*m + 2*n + max(3*n  ,lwork) + 2*SQR(n) + 11*L
   *   ints:      n
   */

  /* dblwork is used for dposvx_() and dgelss_()  , */
  /*   need 3*n and lwork doubles respectively. */
  u = dblwork + lwork; /* lwork=10*n */
  fc = u + m;
  gup = fc + L;
  gu = gup + m;
  du = gu + m;
  dz = du + m;
  gx = dz + m, dx = gx + n;
  Hx = dx + n;
  d1 = Hx + ((n) * (n));
  d2 = d1 + L;
  d3 = d2 + L;
  e1 = d3 + L;
  e2 = e1 + L;
  e3 = e2 + L;
  t1 = e3 + L;
  t2 = t1 + L;
  t3 = t2 + L;
  t4 = t3 + L;
  ua = t4 + L;
  za = ua + m;
  Hwork = za + m;
  /* Hwork needs SQR(n) doubles */

  /* gap reduction vs. centering */
  w = 2 * L + Nu * sqrt(2 * L);

  /* u=A*x+b */
  dcopy_(&m, b, &int1, u, &int1); /* u=b */
  dgemv_("N", &m, &n, &double1, A, &m, x, &int1, &double1, u,
         &int1); /* u=A*x+u */

  /* compute gap (and deviation from centrality and store in hist) */
  if (out_mode == 2) {
    dgapdev(m, L, N, u, z, &gap, &dev);
    hist[0] = gap;
    hist[1] = dev;
  } else {
    gap = ddot_(&m, u, &int1, z, &int1); /* gap = u'*z; */
    if (out_mode == 1)
      hist[0] = gap;
  }

  /* outer loop */
  iter_out = 0;
  while (!((rel_tol < 0.0 &&
            (xtarget = (ddot_(&n, f, &int1, x, &int1) < target ||
                        -ddot_(&m, b, &int1, z, &int1) >= target))) ||
           (abs_tol > 0.0 && (xabs = (gap <= abs_tol))) ||
           (rel_tol > 0.0 &&
            (xrel = (((lbound = -ddot_(&m, b, &int1, z, &int1)) > 0.0 &&
                      gap / lbound <= rel_tol) ||
                     ((lbound = -ddot_(&n, f, &int1, x, &int1)) > 0.0 &&
                      gap / lbound <= rel_tol)))) ||
           (xiter = (iter_out >= *iter)))) {
    ++iter_out;

    /* compute gup (gradient of primal barrier wrt. u) */
    /* also  , compute fc(i)=2/(t(i)^2-u(i)'*u(i)) */
    for (i = 0, k = 0; i < L; ++i) {
      for (j = 0, fc[i] = 0.0; j < N[i] - 1; ++j)
        fc[i] -= ddsqr(u[k + j]);
      fc[i] += ddsqr(u[k + j]);
      fc[i] = 2.0 / fc[i];
      for (j = 0; j < N[i] - 1; ++j, ++k)
        gup[k] = fc[i] * u[k];
      gup[k] = -fc[i] * u[k];
      ++k;
    }

    /* compute gu (gradient of potential wrt. u) */
    s = w / gap;
    dcopy_(&m, gup, &int1, gu, &int1);   /* gu=gup */
    daxpy_(&m, &s, z, &int1, gu, &int1); /* gu=gu+(w/gap)*z */

    /* compute Hx = A'*Hu*A */
    /* where   Hu(i) = fc(i)*diag([1 1 ... 1 -1]) + gup(i)*gup(i)' */
    /* Hx=0 */
    dzero((n) * (n), Hx);

    for (i = 0, k = 0; i < L; k += N[i++]) { /* for each constraint */
      /* n. of rows of A(i) */
      j = N[i] - 1;

      /* Hx = Hx + fc(i)*A(i)'*A(i) */
      if (j > 0)
        dsyrk_("U", "T", &n, &j, &(fc[i]), &(A[k]), &m, &double1, Hx, &n);

      /* Hx = Hx - fc(i)*c(i)*c(i)'  (rank one update) */
      s = -fc[i];
      dsyr_("U", &n, &s, &(A[k + j]), &m, Hx, &n);

      /* gx = [A(i);c(i)']'*gup(i)  (this is not gx  , just used as aux.) */
      dgemv_("T", &(N[i]), &n, &double1, &(A[k]), &m, &(gup[k]), &int1,
             &double0, gx, &int1);

      /* Hx = Hx + gx*gx'  (rank one update) */
      dsyr_("U", &n, &double1, gx, &int1, Hx, &n);
    }

    /* solve linear system: dx = -Hx\(A'*gu) */
    /* gx = -A'*gu */
    dgemv_("T", &m, &n, &double_1, A, &m, gu, &int1, &double0, gx, &int1);
    /* dx = Hx\gx */
    if (!use_ls) { /* solve linear system by QR fact. */
      dposvx_("N", "U", &n, &int1, Hx, &n, Hwork, &n, &equed, &scale, gx, &n,
              dx, &n, &rcond, &ferr, &berr, dblwork, intwork, &la_info);
      if (la_info > 0) /* from 1 to n  , Hessian not positive def.; */
        use_ls = 1;    /* n+1  , Hessian badly conditioned; */
                       /* do SVD now and switch to SVD for all */
                       /* future iterations */
    }
    if (use_ls) {   /* solve linear system in least squares sense using SVD */
      dupge(n, Hx); /* convert to general storage */
      rcond = -1;   /* keep singular values down to machine precision */
      dgelss_(&n, &n, &int1, Hx, &n, gx, &n, Hwork, &rcond,
              &rank, /* (only first n of Hwork used  , for S) */
              dblwork, &lwork, &la_info); /* dblwork: lwork doubles (>= 5*n) */
      dcopy_(&n, gx, &int1, dx, &int1);   /* dx=gx (dgelss_ overwrites gx) */
    }
    if (la_info)
      return la_info; /* abort: return error in lapack solver */

    /* du = A*dx */
    dgemv_("N", &m, &n, &double1, A, &m, dx, &int1, &double0, du, &int1);

    /* dz = -(gu+Hu*du) */
    /* computed one constraint at a time: */
    /* dz(i)= -(gu(i)+fc(i)*[du(i);-dt(i)]+gup(i)*(gup(i)'*[du(i);dt(i)])) */
    /* dz = gu */
    dcopy_(&m, gu, &int1, dz, &int1);
    for (i = 0, k = 0; i < L; k += N[i++]) {
      j = N[i] - 1;
      /* dz(i) = dz(i) + fc[i]*[du(i);-dt(i)] */
      daxpy_(&j, &(fc[i]), &(du[k]), &int1, &(dz[k]), &int1);
      dz[k + j] -= fc[i] * du[k + j];
      /* s = gup(i)'*du(i) */
      s = ddot_(&(N[i]), &(gup[k]), &int1, &(du[k]), &int1);
      /* dz(i) = dz(i) + s*gup(i) */
      daxpy_(&(N[i]), &s, &(gup[k]), &int1, &(dz[k]), &int1);
    }
    /* dz = - dz */
    dscal_(&m, &double_1, dz, &int1);
    /* optional: scale dz by 1/rho=gap/w
      s=gap/w;
      dscal_(&m  ,&s  ,dz  ,&int1);
    */

    /*
     *  constants for plane search
     */
    c1 = gap;
    c2 = ddot_(&m, du, &int1, z, &int1);
    c3 = ddot_(&m, u, &int1, dz, &int1);
    for (i = 0, k = 0; i < L; ++i) {
      d1[i] = 0.0;
      d2[i] = 0.0;
      d3[i] = 0.0;
      e1[i] = 0.0;
      e2[i] = 0.0;
      e3[i] = 0.0;
      for (j = 0; j < N[i] - 1; ++j, ++k) {
        d1[i] -= ddsqr(u[k]);
        d2[i] -= u[k] * du[k];
        d3[i] -= ddsqr(du[k]);
        e1[i] -= ddsqr(z[k]);
        e2[i] -= z[k] * dz[k];
        e3[i] -= ddsqr(dz[k]);
      }
      d1[i] += ddsqr(u[k]);
      d2[i] += u[k] * du[k];
      d3[i] += ddsqr(du[k]);
      e1[i] += ddsqr(z[k]);
      e2[i] += z[k] * dz[k];
      e3[i] += ddsqr(dz[k]);
      ++k;
    }

    /* plane search loop */
    p = 0.0;
    q = 0.0;
    for (lambda2 = 2 * MAX_LAMBDA2, iter_plane = 0;
         lambda2 > MAX_LAMBDA2 && iter_plane < MAX_ITER_PLANE; ++iter_plane) {

      /* compute gradient and Hessian wrt. p and q */
      /*   at u+p*du  , z+q*dz */
      dgrad(w, L, c1, c2, c3, d1, d2, d3, e1, e2, e3, p, q, &gp, &gq, t1, t2,
            t3, t4);
      dhess(L, d3, e3, t1, t2, t3, t4, &hp, &hq);

      /* Newton step */
      dp = -gp / hp;
      dq = -gq / hq;

      /* line search loop: scale down step */
      /*   until inside feasible region */
      /*   (and until between previous point and line minimum) */
      alpha = 1; /* scaling factor for dp  , dq */
      p0 = p;
      q0 = q;
      while (1) {
        p = p0 + alpha * dp;
        q = q0 + alpha * dq;

        /* ua=u+p*du */
        dcopy_(&m, u, &int1, ua, &int1);
        daxpy_(&m, &p, du, &int1, ua, &int1);
        /* za=z+q*dz */
        dcopy_(&m, z, &int1, za, &int1);
        daxpy_(&m, &q, dz, &int1, za, &int1);

        /* check constraints: */
        /*   fu = t(i)^2 - u(i)'*u(i) > 0  and  t(i) > 0 */
        /*   fz = s(i)^2 - z(i)'*z(i) > 0  and  s(i) > 0 */
        for (i = 0, k = 0, out_barrier = 0; i < L && !out_barrier; ++i) {
          for (j = 0, fu = 0.0, fz = 0.0; j < N[i] - 1; ++j, ++k) {
            fu -= ddsqr(ua[k]);
            fz -= ddsqr(za[k]);
          }
          fu += ddsqr(ua[k]);
          fz += ddsqr(za[k]);
          if (fu <= 0 || fz <= 0 || ua[k] <= 0 || za[k] <= 0) {
            out_barrier = 1; /* (replace <=0 with <TINY ?) */
            break;
          }
          ++k;
        }

        if (!out_barrier) {
          /* compute gradient along search line wrt. alpha */
          dgrad(w, L, c1, c2, c3, d1, d2, d3, e1, e2, e3, p, q, &gp, &gq, t1,
                t2, t3, t4);
          galpha = dp * gp + dq * gq;
          /* exit if: all barriers ok (i.e. out_barrier==0) */
          /*  and gradient negative (=> between initial and minimum) */
          if (galpha <= 0)
            break; /* EXIT line search */
        }
        alpha /= DIV_ALPHA;      /* scale down the p  , q step */
        if (alpha < MIN_ALPHA) { /* line search failed */
          alpha = 0.0;
          p = p0;
          q = q0;
          break; /* EXIT line search */
        }
      } /* end of line search loop */

      /* plane search convergence criterium */
      lambda2 = hp * ddsqr(dp) + hq * ddsqr(dq);

    } /* end of plane search loop */

    /* x=x+p*dx */
    daxpy_(&n, &p, dx, &int1, x, &int1);
    /* z=z+q*dz */
    daxpy_(&m, &q, dz, &int1, z, &int1);
    /* u=A*x+b*/
    dcopy_(&m, b, &int1, u, &int1);
    dgemv_("N", &m, &n, &double1, A, &m, x, &int1, &double1, u, &int1);

    /* update gap (and deviation from centrality and store in hist) */
    if (out_mode == 2) {
      dgapdev(m, L, N, u, z, &gap, &dev);
      hist[2 * iter_out] = gap;
      hist[2 * iter_out + 1] = dev;
    } else {
      gap = ddot_(&m, u, &int1, z, &int1); /* gap = u'*z; */
      if (out_mode == 1)
        hist[iter_out] = gap;
    }
  } /* end of outer loop */
  *iter = iter_out;

  /* report reason for exit */
  *info = xabs + 2 * xrel + 3 * xtarget + 4 * xiter;
  /* normal exit */
  return 0;
} /* end of socp() */
