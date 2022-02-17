/*  Author: AdA */
/*  Date:   22/06/1999 */

/*  Associated header: sdplib.h */
/*  Uses: some nr routines and some BLAS routines.  */

/*  BLAS from the Intel Math Kernel library can be linked directly        ,  */
/*  or the linear algebra fucntions provided for use in griffin can be used
 * instead. */
/*  See sdplib.h for BLAS library choice. */

/*  ------------------------------------------------------------------------------------------------------
 */
/*  This is the main toolbox used by sdp */
/*  ------------------------------------------------------------------------------------------------------
 */

#include "math.h"

#include "sdp_sdplib.h"

/*  **************************************************************** */
/*  Part One: Linear algebra basic toolbox for use in SDP. --------- */
/*  ---------------------------------------------------------------- */
/*  **************************************************************** */

/*  Copy A into B */
Err copy_mat(double **mata, int dim, double **matres) {
  int i, j;
  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      matres[i][j] = mata[i][j];
  return NULL;
}

/*  zeros a square matrix */
Err zero_mat(double **mata, int dim) {
  int i, j;

  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      mata[i][j] = 0;
  return NULL;
}

/*  zeros a vector */
Err zero_vec(double *veca, int dim) {
  int i;

  for (i = 1; i <= dim; i++)
    veca[i] = 0;
  return NULL;
}

/*  Stack a matrix of dimension dim into a [1        ,...        ,(dim*dim)]
 * vector. */
Err stack_mat(double **mata, int dim, double *vecres) {
  int i, j;
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      vecres[(j - 1) * dim + (i - 1) + 1] = mata[i][j];
  return NULL;
}

/*  Adds aplha*A+beta*B put it into C (if transpose="T"        , B is
 * transposed) */
Err add_mat(double alpha, double **mata, double beta, double **matb, int dim,
            double **matc, char *transpose) {
  int i, j;
  if (strncmp(transpose, "T", 1) == 0) {
    for (i = 1; i <= (dim); i++)
      for (j = 1; j <= (dim); j++)
        matc[i][j] = alpha * mata[i][j] + beta * matb[j][i];
  } else {
    for (i = 1; i <= (dim); i++)
      for (j = 1; j <= (dim); j++)
        matc[i][j] = alpha * mata[i][j] + beta * matb[i][j];
  }
  return NULL;
}

/*  Max of two doubles */
double d_max(double a, double b) {
  if (a >= b)
    return a;
  else
    return b;
}

#ifndef CBLAS
/*  Simple product of standard [1..n][1..n] matrixes when speed is not an issue.
 */
/*  All matrixes not symmetric. Stand-Alone. */
/*  result can be same matrix as mata or matb. */
/*  It also computes B*A instead of A*B because it proves to */
/*  improve the precision of the algo in the Lyapounov equation. */
Err mat_mult_nobuf(double **mata, double **matb, int dim, double **result) {
  double buf = 0.0;
  int i, j, k;
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++) {
      buf = 0;
      for (k = 1; k <= dim; k++)
        buf += mata[k][i] * matb[j][k];
      result[j][i] = buf;
    }
  return NULL;
}

/*  Simple product of standard [1..n][1..n] matrixes when speed is not an issue.
 */
/*  All matrixes not symmetric. Stand-Alone. */
/*  result can be same matrix as mata or matb. */
Err mat_mult_buf(double **mata, double **matb, int dim, double **result) {
  double buf = 0.0;
  int i, j, k;
  double **matbuf1 = dmatrix(1, dim, 1, dim);

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++) {
      buf = 0;
      for (k = 1; k <= dim; k++)
        buf += mata[i][k] * matb[k][j];
      matbuf1[i][j] = buf;
    }
  copy_mat(matbuf1, dim, result);
  /*  Free */
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  return NULL;
}

#else

/*  Simple product of standard [1..n][1..n] matrixes. */
/*  All matrixes not symmetric. Using CBLAS from MKL */
Err mat_mult_buf(double **mata, double **matb, int dim, double **result) {
  double buf = 0.0;
  int i, j;
  double *tempmat1 = dvector(0, dim * dim);
  double *tempmat2 = dvector(0, dim * dim);
  double *tempmat3 = dvector(0, dim * dim);

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      tempmat1[(j - 1) * dim + (i - 1)] = mata[i][j];
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      tempmat2[(j - 1) * dim + (i - 1)] = matb[i][j];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, 1.0,
              tempmat1, dim, tempmat2, dim, 0.0, tempmat3, dim);
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      result[i][j] = tempmat3[(j - 1) * dim + (i - 1)];
  /*  free */
  free_dvector(tempmat1, 0, dim * dim);
  free_dvector(tempmat2, 0, dim * dim);
  free_dvector(tempmat3, 0, dim * dim);
  return NULL;
}

#define mat_mult_nobuf mat_mult_buf

#endif

/*  Transpose a square matrix */
Err transp(double **mata, int dim, double **result) {
  int i, j;
  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      result[i][j] = mata[j][i];
  return NULL;
}

/*  Fast image of a standard [1..n][1..n] matrix by a big operator
 * [1..n^2][1..n^2] on matrixes. */
/*  Works only with CBLAS. */
#ifdef CBLAS
Err fast_imat(double **opera, double **mata, int dim, double **result) {
  double *opbuf = dvector(0, (dim * dim * dim * dim - 1));
  double *matbuf = dvector(0, dim * dim - 1);
  double *resbuf = dvector(0, dim * dim - 1);
  int i, j;

  /*  Put the arrays in CBLAS form. */
  for (i = 1; i <= (dim * dim); i++)
    for (j = 1; j <= (dim * dim); j++)
      opbuf[(j - 1) * dim * dim + (i - 1)] = opera[i][j];
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      matbuf[(j - 1) * dim + (i - 1)] = mata[i][j];
  /*  Computes the image. */
  cblas_dgemv(CblasRowMajor, CblasNoTrans, dim * dim, dim * dim, 1.0, opbuf,
              dim * dim, matbuf, 1, 0.0, resbuf, 1);
  /*  Put the resbuf in the standard form result matrix. */
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      result[i][j] = resbuf[(j - 1) * dim + (i - 1)];

  /*  free */
  free_dvector(opbuf, 0, dim * dim * dim * dim - 1);
  free_dvector(matbuf, 0, dim * dim - 1);
  free_dvector(resbuf, 0, dim * dim - 1);
  return NULL;
}
#endif

/*  Defines the A(*)B tensor product of matrixes of size dim (A        , B ,
 * must be symetric). */
/*  Returns a matrix of dimension dim^2 */
/*  This is the most intensive computation in the algorithm. */
/*  (ToDo) Look for optimized formulation. */
Err tensor_prod(double **mata, double **matb, int dim, double **result) {
  double **basmat = dmatrix(1, dim, 1, dim);
  double **tempmat1 = dmatrix(1, dim, 1, dim);
  double **tempmat2 = dmatrix(1, dim, 1, dim);
  double **tempmat3 = dmatrix(1, dim, 1, dim);
  int i, j;

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++) {
      /*  Computes the image of the e(i        ,j) baisi matrix. */
      basmat[i][j] = 1;
      mat_mult_buf(mata, basmat, dim, tempmat1);
      mat_mult_buf(tempmat1, matb, dim, tempmat2);
      mat_mult_buf(matb, basmat, dim, tempmat1);
      mat_mult_buf(tempmat1, mata, dim, tempmat3);
      /*  Update the big matrix. */
      add_mat(0.5, tempmat2, 0.5, tempmat3, dim, tempmat1, "No");
      stack_mat(tempmat1, dim, result[(j - 1) * dim + (i - 1) + 1]);
      basmat[i][j] = 0;
    }
  /*  free */
  free_dmatrix(basmat, 1, dim, 1, dim);
  free_dmatrix(tempmat1, 1, dim, 1, dim);
  free_dmatrix(tempmat2, 1, dim, 1, dim);
  free_dmatrix(tempmat3, 1, dim, 1, dim);
  return NULL;
}

/*  Computes (A(*)B)(O)        , with A        , B        , O symmetric. */
/*  Returns a matrix of dimension dim */
Err tensor_oper(double **mata, double **matb, double **mato, int dim,
                double **result) {
  double **tempmat1 = dmatrix(1, dim, 1, dim);
  double **tempmat2 = dmatrix(1, dim, 1, dim);
  double **tempmat3 = dmatrix(1, dim, 1, dim);

  /*  Computes the image of the e(i        ,j) baisi matrix. */
  mat_mult_buf(mata, mato, dim, tempmat1);
  mat_mult_buf(tempmat1, matb, dim, tempmat2);
  mat_mult_buf(matb, mato, dim, tempmat1);
  mat_mult_buf(tempmat1, mata, dim, tempmat3);
  /*  Update the big matrix. */
  add_mat(0.5, tempmat2, 0.5, tempmat3, dim, result, "No");
  /*  Symmetrize */
  add_mat(0.5, result, 0.5, result, dim, result, "T");
  /*  free */
  free_dmatrix(tempmat3, 1, dim, 1, dim);
  free_dmatrix(tempmat2, 1, dim, 1, dim);
  free_dmatrix(tempmat1, 1, dim, 1, dim);
  return NULL;
}

/*  Scalar product of two non-symmetric matrixes */
double scal_mat_nosym(double **mata, double **matb, int dim) {
  double buf = 0;
  int i, j;
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      buf += mata[i][j] * matb[i][j];
  return buf;
}

/*  Scalar product of two symmetric matrixes */
double scal_mat(double **mata, double **matb, int dim) {
  double buf = 0;
  int i, j;
  for (i = 1; i <= (dim); i++)
    for (j = i + 1; j <= (dim); j++)
      buf += mata[i][j] * matb[i][j];
  buf *= 2;
  for (i = 1; i <= (dim); i++)
    buf += mata[i][i] * matb[i][i];
  return buf;
}

/*  Scalar product of two vectors. */
double scal_vec(double *veca, double *vecb, int dim) {
  double buf = 0;
  int i;

  for (i = 1; i <= dim; i++)
    buf += veca[i] * vecb[i];
  return buf;
}

/*  Adds two scaled vectors alpha*VECA + beta*VECB = result */
Err add_vec(double alpha, double *veca, double beta, double *vecb, int dim,
            double *vecres) {
  int i;
  for (i = 1; i <= dim; i++)
    vecres[i] = alpha * veca[i] + beta * vecb[i];
  return NULL;
}

/*  Image of a small vector by a matrix */
Err slow_im_vec(double **mata, double *veca, int dim, double *vecres) {
  int i, j;
  double buf;

  for (i = 1; i <= dim; i++) {
    buf = 0;
    for (j = 1; j <= dim; j++)
      buf += mata[i][j] * (veca[j]);
    vecres[i] = buf;
  }
  return NULL;
}

/*  Test if a matrix is Positive Definite using Cholesky decomposition. */
/*  Returns 0 if PD        , or the size of the first  */
/*  non PD submatrix. */
int pd_test(double **mata, int dim) {
  double *pbuf = dvector(1, dim);
  double **bufm = dmatrix(1, dim, 1, dim);
  int i, j, k;
  double sum = 0;
  int testpd = 1, ret = 0;

  copy_mat(mata, dim, bufm);
  for (i = 1; i <= dim; i++) {
    for (j = i; j <= dim; j++) {
      sum = bufm[i][j];
      for (k = i - 1; k >= 1; k--)
        sum -= bufm[i][k] * bufm[j][k];
      if (i == j) {
        if (sum < 0.0 && ret == 0) {
          testpd = 0;
          ret = i;
        } else
          pbuf[i] = sqrt(sum);
      } else
        bufm[j][i] = sum / pbuf[i];
    }
  }
  /*  free */
  free_dvector(pbuf, 1, dim);
  free_dmatrix(bufm, 1, dim, 1, dim);
  return ret;
}

/*  Solves a linear system A.x=b        , returns inv(A) and x. if algo_type is
 * 0        , x otherwise. */
/*  Algo_type gives the choice between 0:Gauss elimination        , 1:Gauss
 * Siedel. */
Err lsolve(double **mata, double *vecb, int dim, double *vecres,
           double **matres, int algo_type, double prec_obj, int num_improve) {
  Err err = NULL;
  int i, j;
  double bufra = 0, bufrb = 0;
  double **matbuf = dmatrix(1, dim, 1, 1);
  double **matbuf1 = dmatrix(1, dim, 1, dim);
  double *vecbuf = dvector(1, dim);
  double *vecbuf1 = dvector(0, dim);

  if (algo_type == 0) {
    copy_mat(mata, dim, matres);
    for (i = 1; i <= (dim); i++)
      matbuf[i][1] = vecb[i];
    err = gaussj(matres, dim, matbuf, 1);
    for (i = 1; i <= (dim); i++)
      vecres[i] = matbuf[i][1];
    /*  Start improvement */
    bufrb = 2;
    for (i = 1; i <= 4 * num_improve; i++) {
      slow_im_vec(mata, vecres, dim, vecbuf);
      for (j = 1; j <= dim; j++)
        vecbuf1[j] = (vecb[j] - vecbuf[j]);
      if (i == 1) {
        bufra = sqrt(scal_vec(vecbuf1, vecbuf1, dim));
        if (bufra == 0)
          break;
      }
      slow_im_vec(matres, vecbuf1, dim, vecbuf);
      for (j = 1; j <= dim; j++)
        vecres[j] += vecbuf[j];
    }
    if (sqrt(scal_vec(vecbuf1, vecbuf1, dim)) > prec_obj)
      err = serror("Not enough precision in lsolve...");
  } else {
    err = serror("He He He...");
  }
  /*  free */
  free_dmatrix(matbuf, 1, dim, 1, 1);
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  free_dvector(vecbuf, 1, dim);
  free_dvector(vecbuf1, 0, dim);
  return err;
}

/*  Returns the inverse of a matrix */
Err inv_mat(double **mata, int dim, double **result) {
  double *vecbuf1 = dvector(1, dim);
  double *vecbuf2 = dvector(1, dim);
  Err err = NULL;
  int i;

  for (i = 1; i <= dim; i++)
    vecbuf1[i] = i;
  err = lsolve(mata, vecbuf1, dim, vecbuf2, result, 0, 0, 0);
  /* free */
  free_dvector(vecbuf1, 1, dim);
  free_dvector(vecbuf2, 1, dim);
  return err;
}

/*  Computes the minimum eigenvalue of a symmetric matrix. */
Err min_eigenval(double **mata, int dim, double *result, double *condition) {
  double bufmin = 0, bufmax;
  int i, j;
  Err err = NULL;
  double *evals = dvector(0, dim - 1);
  double **evects = dmatrix(0, dim - 1, 0, dim - 1);
  double **matbuf = dmatrix(0, dim - 1, 0, dim - 1);

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      matbuf[i - 1][j - 1] = mata[i][j];
  err = diagonalise_symmetric_matrix(matbuf, dim, evals, evects);
  bufmin = evals[0];
  bufmax = evals[0];
  for (i = 1; i < (dim - 1); i++)
    if (bufmin > evals[i])
      bufmin = evals[i];
  for (i = 1; i <= dim; i++)
    if (bufmax < evals[i])
      bufmax = evals[i];
  /*  free */
  free_dvector(evals, 0, dim - 1);
  free_dmatrix(evects, 0, dim - 1, 0, dim - 1);
  free_dmatrix(matbuf, 0, dim - 1, 0, dim - 1);
  *result = bufmin;
  if (bufmin != 0)
    *condition = bufmax / bufmin;
  else
    *condition = 1e+308;
  return err;
}

/*  Computes the minimum eigenvalue of a non-symmetric matrix. */
Err min_eigenval_nosym(double **mata, int dim, double *result, int max1_min0) {
  double buf = 0;
  int i, j;
  Err err = NULL;
  double *evals = dvector(0, dim);
  double *ivals = dvector(0, dim);
  double **evects = dmatrix(0, dim, 0, dim);
  double **matbuf = dmatrix(0, dim, 0, dim);

  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      matbuf[i][j] = mata[i][j];
  balanc(matbuf, dim);
  elmhes(matbuf, dim);
  err = hqr(matbuf, dim, evals, ivals);
  buf = evals[1];
  if (max1_min0 == 0) {
    for (i = 1; i <= dim; i++)
      if (buf > evals[i])
        buf = evals[i];
  } else {
    for (i = 1; i <= dim; i++)
      if (buf < evals[i])
        buf = evals[i];
  }
  /*  free */
  free_dvector(evals, 0, dim);
  free_dvector(ivals, 0, dim);
  free_dmatrix(evects, 0, dim, 0, dim);
  free_dmatrix(matbuf, 0, dim, 0, dim);
  *result = buf;
  return err;
}

/*  diagonalisation of a symmetric matrix. */
/*  Everything in [1..dim] format. */
Err diago(double **mata, int dim, double **eigenvec, double *eigenval) {
  double buf = 0;
  int i, j;
  Err err = NULL;
  double *evals = dvector(0, dim - 1);
  double **evects = dmatrix(0, dim - 1, 0, dim - 1);
  double **matbuf = dmatrix(0, dim - 1, 0, dim - 1);

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      matbuf[i - 1][j - 1] = mata[i][j];
  err = diagonalise_symmetric_matrix(matbuf, dim, evals, evects);
  for (i = 1; i <= (dim); i++)
    eigenval[i] = evals[i - 1];
  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      eigenvec[i][j] = evects[j - 1][i - 1];
  /*  free */
  free_dvector(evals, 0, dim - 1);
  free_dmatrix(evects, 0, dim - 1, 0, dim - 1);
  free_dmatrix(matbuf, 0, dim - 1, 0, dim - 1);
  return err;
}

double condition_mat(double **mata, int dim) {
  double res = 0, buf;
  Err err = NULL;
  double **matbuf = dmatrix(1, dim, 1, dim);

  transp(mata, dim, matbuf);
  mat_mult_buf(matbuf, mata, dim, matbuf);
  err = min_eigenval(matbuf, dim, &buf, &res);
  res = sqrt(fabs(res));
  /*  Free. */
  free_dmatrix(matbuf, 1, dim, 1, dim);
  if ((err == NULL) && (res * res != 1e+308))
    return res;
  else
    return 1.0e+308;
}

Err fudge_eigenvals(double *evals, int dim) {
  int i;
  double bufmax;
  double machine_epsilon = 2e-16;

  bufmax = evals[1];
  for (i = 2; i <= dim; i++)
    if (evals[i] > bufmax)
      bufmax = evals[i];
  for (i = 1; i <= dim; i++)
    evals[i] += 5.0 * bufmax * machine_epsilon;

  return NULL;
}

/*  *********************************************************************************************************
 */
/*  Part Two: SDP specific routines.
 * ------------------------------------------------------------------------ */
/*  ---------------------------------------------------------------------------------------------------------
 */
/*  *********************************************************************************************************
 */

/*  Computes the primal constraints operator (Tr(A1.X)        , ... ,Tr(Ai.X) ,
 * ....) operating on a matrix mata. */
Err op_a(double ***con_mat, int num_const, double **mata, int dim,
         double *vecres) {
  int k;
  for (k = 1; k <= num_const; k++)
    vecres[k] = scal_mat(con_mat[k], mata, dim);
  return NULL;
}

/*  Computes the augmented constraints operator (A        ,-C). */
/*  Returns a num_const+1 vector. */
Err op_a_augm(double ***con_mat, int num_const, double **cmat, double **mata,
              int dim, double *vecres) {
  int k;
  for (k = 1; k <= num_const; k++)
    vecres[k] = scal_mat(con_mat[k], mata, dim);
  vecres[num_const + 1] = -scal_mat(cmat, mata, dim);
  return NULL;
}

/*  Computes the transposition of the constraints operator operating on a vector
 * veca. */
Err op_ta(double ***con_mat, int num_const, double *veca, int dim,
          double **matres) {
  int i, j, k;

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      matres[i][j] = 0;
  for (k = 1; k <= num_const; k++) {
    add_mat(veca[k], con_mat[k], 1.0, matres, dim, matres, "No");
  }
  return NULL;
}

/*  Computes the transposition of the augmented constraints operator operating
 * on a vector veca. */
Err op_ta_augm(double ***con_mat, int num_const, double **cmat, double *veca,
               int dim, double **matres) {
  int i, j, k;

  for (i = 1; i <= (dim); i++)
    for (j = 1; j <= (dim); j++)
      matres[i][j] = 0;
  for (k = 1; k <= num_const; k++) {
    add_mat(veca[k], con_mat[k], 1.0, matres, dim, matres, "No");
  }
  add_mat(-veca[num_const + 1], cmat, 1.0, matres, dim, matres, "No");
  return NULL;
}

/*  Computes the basic step length. */
/*  Returns 0 if Z or X matrix singular */
double step_length(double **mata, double **dmata, int dim, double tau,
                   double dtau, double kappa, double dkappa, double gamma) {
  double res, buf = 0;
  double **matbuf1 = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  Err err = NULL;

  err = inv_mat(mata, dim, matbuf1);
  res = 1.0;
  if (err != NULL)
    res = 1.0;
  if (err == NULL) {
    mat_mult_buf(matbuf1, dmata, dim, matbuf2);
    err = min_eigenval_nosym(matbuf2, dim, &buf, 0);
    if ((err == NULL) && (buf < 0)) {
      if ((-gamma / buf) < res)
        res = (-gamma / buf);
    }
  }

  if (tau * dtau < 0) {
    if (((-gamma * tau) / (dtau)) < res)
      res = (-gamma * tau) / (dtau);
  }

  if (kappa * dkappa < 0) {
    if (((-gamma * kappa) / (dkappa)) < res)
      res = (-gamma * kappa) / (dkappa);
  }

  res = -d_max(-gamma, -res);
  /*  free */
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  return res;
}

/*  This computes phi        , the infeasibility measure. */
double compute_infeas(double ***const_mat, double *const_val, int num_const,
                      double **xmat, double **zmat, double **cmat, int dim,
                      double *yvec, double tau, double *prim, double *dual) {
  double bufa = 0, bufb = 0, resa, resb;
  double *vecbuf1 = dvector(1, num_const);
  double **matbuf1 = dmatrix(1, dim, 1, dim);

  /*  Relative Primal Infeasibility. */
  op_a(const_mat, num_const, xmat, dim, vecbuf1);
  add_vec(tau, const_val, -1.0, vecbuf1, num_const, vecbuf1);
  bufa = sqrt(scal_vec(vecbuf1, vecbuf1, num_const));
  resa = bufa / (tau * (1 + sqrt(scal_vec(const_val, const_val, num_const))));
  /*  Relative Dual Infeasibility. */
  op_ta(const_mat, num_const, yvec, dim, matbuf1);
  add_mat(-1.0, zmat, -1.0, matbuf1, dim, matbuf1, "No");
  add_mat(tau, cmat, 1.0, matbuf1, dim, matbuf1, "No");
  bufb = sqrt(scal_mat(matbuf1, matbuf1, dim));
  resb = bufb / (tau * (1 + sqrt(scal_mat(cmat, cmat, dim))));
  /*  free */
  free_dvector(vecbuf1, 1, num_const);
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  /*  Returns the max of the two. */
  *prim = resa;
  *dual = resb;
  return d_max(resa, resb);
}

/*  This computes sigma */
double compute_sigma(double **xmat, double **dxmat, double **zmat,
                     double **dzmat, int dim, double alpha, double beta,
                     double tau, double dtau, double kappa, double dkappa,
                     double expon) {
  double **matbuf1 = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  double buf = 0, bufb;

  add_mat(1.0, xmat, alpha, dxmat, dim, matbuf1, "No");
  add_mat(1.0, zmat, beta, dzmat, dim, matbuf2, "No");
  buf = scal_mat(matbuf1, matbuf2, dim);
  buf += (tau + alpha * dtau) * (kappa + beta * dkappa);
  buf /= (scal_mat(xmat, zmat, dim) + tau * kappa);
  buf = fabs(buf);
  buf = pow(buf, expon);
  buf = -d_max(-1.0, -buf);
  bufb = -d_max(-alpha, -beta);
  if (bufb >= 0.9) {
    buf = d_max(0.05, buf);
  } else {
    if (bufb >= 0.7) {
      buf = d_max(0.1, buf);
    } else {
      buf = d_max(0.2, buf);
    }
  }
  /*  free */
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  return buf;
}

/*  Computes mu. */
double compute_mu(double **xmat, double **zmat, int dim, double tau,
                  double kappa) {
  return (scal_mat(xmat, zmat, dim) + tau * kappa) / (dim + 1);
}

/*  This computes the inv(E) operator.  */
Err compute_inv_E(double **pmat, double **zmat, int dim, double **big_inv_E) {
  double **bigbuf1 = dmatrix(1, dim * dim, 1, dim * dim);
  double **bufmat = dmatrix(1, dim, 1, dim);
  Err err = NULL;

  mat_mult_buf(pmat, zmat, dim, bufmat);
  tensor_prod(pmat, bufmat, dim, bigbuf1);
  err = inv_mat(bigbuf1, dim * dim, big_inv_E);

  free_dmatrix(bigbuf1, 1, dim * dim, 1, dim * dim);
  free_dmatrix(bufmat, 1, dim, 1, dim);
  return err;
}

/*  This computes the inv(E) operator with the fast Kronecker product algo.  */
/*  P and Z are supposed to be commuting symmetric. (It is the case for AHO and
 * HKM). */
Err fast_compute_inv_E(double **pmat, double **zmat, int dim,
                       double **big_inv_E) {

  double **bufmat1 = dmatrix(1, dim, 1, dim);
  double **bufmat2 = dmatrix(1, dim, 1, dim);
  double **bufmat3 = dmatrix(1, dim, 1, dim);
  double **bufmat4 = dmatrix(1, dim, 1, dim);
  double **basmat = dmatrix(1, dim, 1, dim);
  double *bufvec1 = dvector(1, dim);
  double *bufvec2 = dvector(1, dim);
  int i, j, k, l;
  Err err = NULL, errbuf = NULL;

  /*  First diagonalise the two matrixes.  */
  /*  Store the common basis of eigenvectors in bufmat2. */
  /*  (for AHO        , in general we would have to find a simultanous diag
   * basis). */
  errbuf = diago(pmat, dim, bufmat1, bufvec1);
  if (err == NULL)
    err = errbuf;
  errbuf = diago(zmat, dim, bufmat2, bufvec2);
  if (err == NULL)
    err = errbuf;
  transp(bufmat2, dim, bufmat3);
  /*  Compute the image of a basis matrix. */
  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++) {
      basmat[i][j] = 1.0;
      mat_mult_nobuf(bufmat3, basmat, dim, bufmat1);
      mat_mult_nobuf(bufmat1, bufmat2, dim, bufmat4);
      for (k = 1; k <= dim; k++)
        for (l = 1; l <= dim; l++) {
          if ((0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k])) != 0)
            bufmat4[k][l] /=
                (0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k]));
          else
            err = serror("Singular system matrix building");
        }
      mat_mult_nobuf(bufmat2, bufmat4, dim, bufmat1);
      mat_mult_nobuf(bufmat1, bufmat3, dim, bufmat4);
      stack_mat(bufmat4, dim, big_inv_E[(j - 1) * dim + (i - 1) + 1]);
      basmat[i][j] = 0;
    }
  /*  Free. */
  free_dmatrix(bufmat1, 1, dim, 1, dim);
  free_dmatrix(bufmat2, 1, dim, 1, dim);
  free_dmatrix(bufmat3, 1, dim, 1, dim);
  free_dmatrix(bufmat4, 1, dim, 1, dim);
  free_dmatrix(basmat, 1, dim, 1, dim);
  free_dvector(bufvec1, 1, dim);
  free_dvector(bufvec2, 1, dim);
  return err;
}

/*  This computes the inv(E) operator image of a matrix with the fast Kronecker
 * product algo.  */
/*  P and Z are supposed to be commuting symmetric. (It is the case for AHO and
 * HKM). */
Err fast_oper_inv_E(double **pmat, double **zmat, double **mata, int dim,
                    double **result, double *error_ratio, int num_improve,
                    double **bufmat2, double *bufvec2, int diag_switch) {

  double **bufmat1 = dmatrix(1, dim, 1, dim);
  double **bufmat3 = dmatrix(1, dim, 1, dim);
  double **bufmat4 = dmatrix(1, dim, 1, dim);
  double **bufmat5 = dmatrix(1, dim, 1, dim);
  double *bufvec1 = dvector(0, dim);
  double error_improvement, buf;
  int k, l, m;
  Err err = NULL, errbuf = NULL;
  /*  This sets the precsion goal in solution improvement */
  double max_prec_improve = 10e-10;

  /*  First diagonalise the two matrixes.  */
  /*  Stores the common basis of eigenvectors in bufmat2. */
  /*  (for AHO here. In general we would have to find a simultanous diag basis).
   */
  for (k = 1; k <= dim; k++)
    bufvec1[k] = 1;
  if (diag_switch == 0) {
    errbuf = diago(zmat, dim, bufmat2, bufvec2);
    /*  Adjusting eigenvalues... */
    fudge_eigenvals(bufvec2, dim);
    if (err == NULL)
      err = errbuf;
  }
  /*  ------------------------------------ */
  transp(bufmat2, dim, bufmat3);
  /*  Compute the image of a basis matrix. */
  mat_mult_nobuf(bufmat3, mata, dim, bufmat1);
  mat_mult_nobuf(bufmat1, bufmat2, dim, bufmat4);
  for (k = 1; k <= dim; k++)
    for (l = 1; l <= dim; l++) {
      if ((0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k])) != 0)
        bufmat4[k][l] /=
            (0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k]));
      else
        err = serror("Singular system matrix operation");
    }
  mat_mult_nobuf(bufmat2, bufmat4, dim, bufmat1);
  mat_mult_nobuf(bufmat1, bufmat3, dim, result);
  /*  Result should be symmetric but  */
  /*  symmetrize it anyway to improve precision. (CPU cost is zero). */
  add_mat(0.5, result, 0.5, result, dim, result, "T");
  /*  Improve the solution by iterative method. */
  for (m = 1; m <= num_improve; m++) {
    tensor_oper(pmat, zmat, result, dim, bufmat4);
    add_mat(-1.0, bufmat4, 1.0, mata, dim, bufmat5, "No");
    /*  for error improvement tracking. */
    if (m == 1) {
      buf = sqrt(scal_mat(bufmat5, bufmat5, dim));
      if (buf != 0)
        error_improvement = 1 / buf;
      else
        break;
    }
    /*  Breaks if improved enough... */
    buf = sqrt(scal_mat(bufmat5, bufmat5, dim));
    if (error_improvement * buf < max_prec_improve)
      break;
    /*  ------------------------------- */
    mat_mult_nobuf(bufmat3, bufmat5, dim, bufmat1);
    mat_mult_nobuf(bufmat1, bufmat2, dim, bufmat4);
    for (k = 1; k <= dim; k++)
      for (l = 1; l <= dim; l++) {
        if ((0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k])) != 0)
          bufmat4[k][l] /=
              (0.5 * (bufvec1[k] * bufvec2[l] + bufvec1[l] * bufvec2[k]));
        else
          err = serror("Singular system matrix operation");
      }
    mat_mult_nobuf(bufmat2, bufmat4, dim, bufmat1);
    mat_mult_nobuf(bufmat1, bufmat3, dim, bufmat4);
    add_mat(1.0, bufmat4, 1.0, result, dim, result, "No");
    add_mat(0.5, result, 0.5, result, dim, result, "T");
  }
  /*  Result should be symmetric but  */
  /*  symmetrize it anyway to improve precision. (CPU cost is zero) */
  add_mat(0.5, result, 0.5, result, dim, result, "T");
  /*  Compute error improvement. */
  if (num_improve > 0) {
    tensor_oper(pmat, zmat, result, dim, bufmat4);
    add_mat(-1.0, bufmat4, 1.0, mata, dim, bufmat5, "No");
    buf = sqrt(scal_mat(bufmat5, bufmat5, dim));
    error_improvement *= buf;
    *error_ratio = error_improvement;
  } else
    *error_ratio = error_improvement = 1;
  /*  Free. */
  free_dmatrix(bufmat1, 1, dim, 1, dim);
  free_dmatrix(bufmat3, 1, dim, 1, dim);
  free_dmatrix(bufmat4, 1, dim, 1, dim);
  free_dmatrix(bufmat5, 1, dim, 1, dim);
  free_dvector(bufvec1, 0, dim);
  return err;
}

/*  This computes the sytem A.inv(E).F.At matrix. */
/*  (ToDo) Choice between AHO and HKM search directions. */
/*  It returns also inv(E) for further usage. */
/*  Puts the result in a [1..(num_const+1)][1..(num_const+1)] matrix. */
Err system_matrix(double ***const_mat, int num_const, double **xmat,
                  double **zmat, double **cmat, int dim, char *srch,
                  double **result, double **big_inv_E, double *errtrack,
                  int num_improve) {
  double **pmat = dmatrix(1, dim, 1, dim);
  double **bufmat = dmatrix(1, dim, 1, dim);
  double **bufmat1 = dmatrix(1, dim, 1, dim);
  double **bufmat2 = dmatrix(1, dim, 1, dim);
  double *bufvec2 = dvector(0, dim);
  double *vecbas = dvector(1, num_const + 1);
  int i;
  double errval, errmax = 0.0;
  Err err = NULL, errbuf = NULL;

  /*  Initialize pmat for AHO */
  for (i = 1; i <= (dim); i++)
    pmat[i][i] = 1.0;
/*  Computes the big operator inv(E). */
/*  This is the most CPU consuming part of the algo. */
/*  Slow direct method: */
/*  err=compute_inv_E(pmat        ,zmat        ,dim        ,big_inv_E); */
#ifndef fast_method
#ifndef very_slow
  err = fast_compute_inv_E(pmat, zmat, dim, big_inv_E);
#else
  err = compute_inv_E(pmat, zmat, dim, big_inv_E);
#endif
#endif
  /*  Computes the image of every basis vector. */
  if (err == NULL) {
    for (i = 1; i <= (num_const + 1); i++) {
      vecbas[i] = 1.0;
      op_ta_augm(const_mat, num_const, cmat, vecbas, dim, bufmat1);
      tensor_oper(xmat, pmat, bufmat1, dim, bufmat);
/*  For direct method: */
#ifndef fast_method
      fast_imat(big_inv_E, bufmat, dim, bufmat1);
#else
      /*  For "on demand" method. */
      errbuf = fast_oper_inv_E(pmat, zmat, bufmat, dim, bufmat1, &errval,
                               num_improve, bufmat2, bufvec2, (i - 1));
      if (errval > errmax)
        errmax = errval;
      if (err == NULL)
        err = errbuf;
#endif
      op_a_augm(const_mat, num_const, cmat, bufmat1, dim, result[i]);
      vecbas[i] = 0.0;
    }
  }
  /*  Return numerical error. */
  *errtrack = errmax;
  /*  free */
  free_dmatrix(pmat, 1, dim, 1, dim);
  free_dmatrix(bufmat, 1, dim, 1, dim);
  free_dmatrix(bufmat1, 1, dim, 1, dim);
  free_dmatrix(bufmat2, 1, dim, 1, dim);
  free_dvector(vecbas, 1, num_const + 1);
  free_dvector(bufvec2, 0, dim);
  return err;
}

/*  This computes the sytem A.At matrix. */
Err compute_gram(double ***const_mat, int num_const, int dim, double **result) {
  double **bufmat1 = dmatrix(1, dim, 1, dim);
  double *vecbas = dvector(1, num_const);
  int i;
  Err err = NULL;

  /*  Computes the image of every basis vector. */
  for (i = 1; i <= (num_const); i++) {
    vecbas[i] = 1.0;
    op_ta(const_mat, num_const, vecbas, dim, bufmat1);
    op_a(const_mat, num_const, bufmat1, dim, result[i]);
    vecbas[i] = 0.0;
  }
  /*  free */
  free_dmatrix(bufmat1, 1, dim, 1, dim);
  free_dvector(vecbas, 1, num_const);
  return err;
}

/*  This computes the left hand final sytem matrix. */
Err left_matrix(double **matm, int num_const, double *bvec, double kappa,
                double tau, double **result) {
  int i, j;
  for (i = 1; i <= (num_const); i++)
    for (j = 1; j <= (num_const); j++)
      result[i][j] = matm[i][j];
  for (i = 1; i <= (num_const); i++)
    result[i][num_const + 1] = matm[i][num_const + 1] - bvec[i];
  for (i = 1; i <= (num_const); i++)
    result[num_const + 1][i] = matm[num_const + 1][i] + bvec[i];
  result[num_const + 1][num_const + 1] =
      matm[num_const + 1][num_const + 1] + kappa / tau;
  return NULL;
}

/*  Computes the Rd matrix. */
Err Rd_mat(double ***const_mat, int num_const, double **zmat, double **cmat,
           int dim, double *yvec, double tau, double **result) {
  double *vecbuf = dvector(1, num_const + 1);
  double **bufmat = dmatrix(1, dim, 1, dim);
  int i;

  for (i = 1; i <= (num_const); i++)
    vecbuf[i] = yvec[i];
  vecbuf[num_const + 1] = tau;
  op_ta_augm(const_mat, num_const, cmat, vecbuf, dim, bufmat);
  add_mat(-1.0, bufmat, -1.0, zmat, dim, result, "No");
  /*  free */
  free_dvector(vecbuf, 1, num_const + 1);
  free_dmatrix(bufmat, 1, dim, 1, dim);
  return NULL;
}

/*  Computes the rp [1..num_const+1] vector. */
Err rp_vec(double ***const_mat, int num_const, double **xmat, double **cmat,
           int dim, double *yvec, double *bvec, double tau, double kappa,
           double *vecres) {
  double buf = 0;
  int i;

  op_a_augm(const_mat, num_const, cmat, xmat, dim, vecres);
  for (i = 1; i <= (num_const + 1); i++)
    vecres[i] = -vecres[i];
  for (i = 1; i <= (num_const); i++)
    vecres[i] += bvec[i] * tau;
  for (i = 1; i <= (num_const); i++)
    buf += bvec[i] * yvec[i];
  vecres[num_const + 1] += (kappa - buf);
  return NULL;
}

/*  Computes the rc scalar */
double rc_val(double sigma, double mu, double tau, double kappa) {
  return sigma * mu - tau * kappa;
}

/*  Computes the Rc matrix. */
Err Rc_mat(double **xmat, double **zmat, int dim, double sigma, double mu,
           char *srch_type, double **result) {
  double **idbuf = dmatrix(1, dim, 1, dim);
  double **matbuf = dmatrix(1, dim, 1, dim);
  int i, j;

  for (i = 1; i <= dim; i++)
    idbuf[i][i] = 1.0;
  mat_mult_buf(xmat, zmat, dim, matbuf);
  /*  (ToDo) implement the HKM search direction */
  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      result[i][j] = 0.5 * (matbuf[i][j] + matbuf[j][i]);
  add_mat(sigma * mu, idbuf, -1.0, result, dim, result, "No");
  /*  free */
  free_dmatrix(idbuf, 1, dim, 1, dim);
  free_dmatrix(matbuf, 1, dim, 1, dim);
  return NULL;
}

/*  Computes the Rq matrix. */
Err Rq_mat(double **xmat, double **dxmat, double **zmat, double **dzmat,
           int dim, double sigma, double mu, char *srch, double **result) {
  double **idbuf = dmatrix(1, dim, 1, dim);
  double **matbuf = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  int i, j;

  for (i = 1; i <= dim; i++)
    idbuf[i][i] = 1.0;
  mat_mult_buf(xmat, zmat, dim, matbuf);
  mat_mult_buf(dxmat, dzmat, dim, matbuf2);
  add_mat(1.0, matbuf2, 1.0, matbuf, dim, matbuf, "No");
  /*  (ToDo) implement the HKM search direction */
  for (i = 1; i <= dim; i++)
    for (j = 1; j <= dim; j++)
      result[i][j] = 0.5 * (matbuf[i][j] + matbuf[j][i]);
  add_mat(sigma * mu, idbuf, -1.0, result, dim, result, "No");
  /*  free */
  free_dmatrix(idbuf, 1, dim, 1, dim);
  free_dmatrix(matbuf, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  return NULL;
}

/*  This computes the right hand side of the (dy        ,dtau) system equation.
 */
/*  It uses the big_inv_E computed by system_matrix. */
Err right_vector(double ***const_mat, int num_const, double **big_inv_E,
                 double **xmat, double **zmat, double **cmat, double **dpxmat,
                 double **dpzmat, char *step_type, int dim, double *yvec,
                 double *bvec, double eta, double tau, double sigma, double mu,
                 double kappa, double dptau, double dpkappa, double *vecres,
                 int num_improve) {
  double **matbuf = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  double **matbuf3 = dmatrix(1, dim, 1, dim);
  double **matbuf4 = dmatrix(1, dim, 1, dim);
  double **matbuf8 = dmatrix(1, dim, 1, dim);
  double **idbuf = dmatrix(1, dim, 1, dim);
  double *vecbuf = dvector(1, num_const + 1);
  double *vecbuf2 = dvector(0, dim);
  double errval;
  int i;

  for (i = 1; i <= dim; i++)
    idbuf[i][i] = 1.0;
  /*  ----- */
  rp_vec(const_mat, num_const, xmat, cmat, dim, yvec, bvec, tau, kappa, vecres);
  for (i = 1; i <= (num_const + 1); i++)
    vecres[i] *= eta;
  vecres[num_const + 1] += rc_val(sigma, mu, tau, kappa) / tau;
  /*  For Corrector step. */
  if (strncmp(step_type, "Cor", 3) == 0) {
    vecres[num_const + 1] -= (dptau * dpkappa) / tau;
  }
  /*  ----- */
  Rd_mat(const_mat, num_const, zmat, cmat, dim, yvec, tau, matbuf2);
  tensor_oper(xmat, idbuf, matbuf2, dim, matbuf);
  /*  ----- */
  Rc_mat(xmat, zmat, dim, sigma, mu, "AHO", matbuf4);
  /*  For Corrector step. */
  if (strncmp(step_type, "Cor", 3) == 0) {
    mat_mult_buf(dpxmat, dpzmat, dim, matbuf3);
    mat_mult_buf(dpzmat, dpxmat, dim, matbuf2);
    add_mat(-0.5, matbuf2, 1.0, matbuf4, dim, matbuf4, "No");
    add_mat(-0.5, matbuf3, 1.0, matbuf4, dim, matbuf4, "No");
  }
  /*  ----- */
  add_mat(-1.0, matbuf4, eta, matbuf, dim, matbuf, "No");
#ifndef fast_method
  /*  For direct method */
  fast_imat(big_inv_E, matbuf, dim, matbuf4);
#else
  /*  For "on demand" method. for AHO method. */
  fast_oper_inv_E(idbuf, zmat, matbuf, dim, matbuf4, &errval, num_improve,
                  matbuf8, vecbuf2, 0);
#endif
  op_a_augm(const_mat, num_const, cmat, matbuf4, dim, vecbuf);
  for (i = 1; i <= (num_const + 1); i++)
    vecres[i] += vecbuf[i];
  /*  free */
  free_dmatrix(matbuf, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  free_dmatrix(matbuf3, 1, dim, 1, dim);
  free_dmatrix(matbuf4, 1, dim, 1, dim);
  free_dmatrix(matbuf8, 1, dim, 1, dim);
  free_dmatrix(idbuf, 1, dim, 1, dim);
  free_dvector(vecbuf, 1, num_const + 1);
  free_dvector(vecbuf2, 0, dim);
  return NULL;
}

/*  Computes dZ. */
Err delta_z(double ***const_mat, int num_const, double **zmat, double **cmat,
            int dim, double *yvec, double *dyvec, double tau, double eta,
            double dtau, double **result) {
  double *vecbuf = dvector(1, num_const + 1);
  double **matbuf9 = dmatrix(1, dim, 1, dim);
  int i;

  for (i = 1; i <= (num_const); i++)
    vecbuf[i] = yvec[i];
  vecbuf[num_const + 1] = tau;
  Rd_mat(const_mat, num_const, zmat, cmat, dim, vecbuf, tau, matbuf9);
  for (i = 1; i <= (num_const); i++)
    vecbuf[i] = dyvec[i];
  vecbuf[num_const + 1] = dtau;
  op_ta_augm(const_mat, num_const, cmat, vecbuf, dim, result);
  add_mat(eta, matbuf9, -1.0, result, dim, result, "No");
  /*  free */
  free_dmatrix(matbuf9, 1, dim, 1, dim);
  free_dvector(vecbuf, 1, num_const + 1);
  return NULL;
}

/*  Computes dX */
Err delta_x(double **big_inv_E, double **xmat, double **dpxmat, double **zmat,
            double **dpzmat, double **dczmat, int dim, double sigma, double mu,
            char *step_type, double **result, int num_improve) {
  double **matbuf = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  double **matbuf8 = dmatrix(1, dim, 1, dim);
  double *vecbuf2 = dvector(0, dim);
  double **idbuf = dmatrix(1, dim, 1, dim);
  int i;
  double errval = 0;

  for (i = 1; i <= dim; i++)
    idbuf[i][i] = 1.0;
  Rc_mat(xmat, zmat, dim, sigma, mu, "AHO", matbuf);
  /*  Predictor or Corrector... */
  if (strncmp(step_type, "Cor", 3) == 0) {
    tensor_oper(xmat, idbuf, dczmat, dim, matbuf2);
  } else {
    tensor_oper(xmat, idbuf, dpzmat, dim, matbuf2);
  }
  add_mat(1.0, matbuf, -1.0, matbuf2, dim, matbuf2, "No");
  /*  For Corrector step. */
  if (strncmp(step_type, "Cor", 3) == 0) {
    mat_mult_buf(dpxmat, dpzmat, dim, matbuf);
    mat_mult_buf(dpzmat, dpxmat, dim, matbuf8);
    add_mat(-0.5, matbuf8, 1.0, matbuf2, dim, matbuf2, "No");
    add_mat(-0.5, matbuf, 1.0, matbuf2, dim, matbuf2, "No");
  }
#ifndef fast_method
  /*  For direct method */
  fast_imat(big_inv_E, matbuf2, dim, result);
#else
  /*  For "on demand" method. for AHO method. */
  fast_oper_inv_E(idbuf, zmat, matbuf2, dim, result, &errval, num_improve,
                  matbuf8, vecbuf2, 0);
#endif
  /*  Result should be symmetric but  */
  /*  symmetrize it anyway to improve precision. (CPU cost is zero) */
  add_mat(0.5, result, 0.5, result, dim, result, "T");
  /*  free */
  free_dmatrix(matbuf, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  free_dmatrix(matbuf8, 1, dim, 1, dim);
  free_dvector(vecbuf2, 0, dim);
  free_dmatrix(idbuf, 1, dim, 1, dim);
  return NULL;
}

/*  Creates an initial scaled iterate for X and Z. */
Err init_sol(double ***const_mat, int num_const, double *const_val,
             double **cmat, int dim, double **startx, double **startz) {
  double **idbuf = dmatrix(1, dim, 1, dim);
  double xhi, eta, bufa = 0, bufb = 0;
  int k, i;

  for (i = 1; i <= dim; i++)
    idbuf[i][i] = 1.0;
  for (k = 1; k <= num_const; k++) {
    if (sqrt(scal_mat(const_mat[k], const_mat[k], dim)) > bufa)
      bufa = sqrt(scal_mat(const_mat[k], const_mat[k], dim));
  }
  if (sqrt(scal_mat(cmat, cmat, dim)) > bufa)
    bufa = sqrt(scal_mat(cmat, cmat, dim));
  eta = bufa / sqrt(dim);
  for (k = 1; k <= num_const; k++) {
    if ((1 + fabs(const_val[k])) /
            (1 + sqrt(scal_mat(const_mat[k], const_mat[k], dim))) >
        bufb)
      bufb = (1 + fabs(const_val[k])) /
             (1 + sqrt(scal_mat(const_mat[k], const_mat[k], dim)));
  }
  xhi = dim * bufb;
  add_mat(xhi, idbuf, 0.0, startx, dim, startx, "No");
  add_mat(eta, idbuf, 0.0, startz, dim, startz, "No");
  /*  Start with Id to match Mathematica package. *** DEBUG ******* */
  /* add_mat(1.0        ,idbuf        ,0.0        ,startx        ,dim ,startx
   * ,"No"); */
  /* add_mat(1.0        ,idbuf        ,0.0        ,startz        ,dim ,startz
   * ,"No"); */
  /*  free */
  free_dmatrix(idbuf, 1, dim, 1, dim);
  return NULL;
}

/*  Computes the centrality measure.  (1-min_eigenval(X.Z)/mu) */
double centrality(double **xmat, double **zmat, int dim, double tau,
                  double kappa) {
  double buf = 0, bufb = 0;
  double **matbuf = dmatrix(1, dim, 1, dim);
  Err err = NULL;

  buf = 1 / scal_mat(xmat, zmat, dim);
  mat_mult_buf(xmat, zmat, dim, matbuf);
  err = min_eigenval_nosym(matbuf, dim, &bufb, 0);
  if (err == NULL)
    buf = 1 - bufb * buf;
  else
    buf = 0;
  /*  Free */
  free_dmatrix(matbuf, 1, dim, 1, dim);
  return buf;
}

/*  Tests if steps actually go towards global progress. */
/*  Set them both to the minimum if not. */
Err set_steps(double ***const_mat, double *const_val, int num_const,
              double *yvec, double *dyvec, double **xmat, double **dxmat,
              double **zmat, double **dzmat, double **cmat, int dim,
              double kappa, double dkappa, double tau, double dtau, double eta,
              double *alpha, double *beta, double frac_search,
              char *step_type) {
  double **matbuf1 = dmatrix(1, dim, 1, dim);
  double **matbuf2 = dmatrix(1, dim, 1, dim);
  double *vecbuf1 = dvector(1, num_const);
  double taubuf1, taubuf2, newtau, newkappa, buf;
  double infa, infb, bufa, bufb;
  int test_step = 1, test_pos = 1, boucl = 0;
  int max_reduc = 60;

  buf = compute_mu(xmat, zmat, dim, tau, kappa) / (tau * tau);
  infa = compute_infeas(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                        yvec, tau, &bufa, &bufb);
  while (test_pos != 0) {
    boucl++;
    /*  First        , set new tau and kappa. */
    taubuf1 = tau + (*alpha) * dtau;
    taubuf2 = tau + (*beta) * dtau;
    if (taubuf1 > taubuf2) {
      newtau = taubuf1;
      newkappa = kappa + (*alpha) * dkappa;
      add_mat(1.0, xmat, *alpha, dxmat, dim, matbuf1, "No");
      add_mat((taubuf1 / taubuf2), zmat, (taubuf1 / taubuf2) * (*beta), dzmat,
              dim, matbuf2, "No");
      add_vec((taubuf1 / taubuf2), yvec, (taubuf1 / taubuf2) * (*beta), dyvec,
              num_const, vecbuf1);
    } else {
      newtau = taubuf2;
      newkappa = kappa + (*beta) * dkappa;
      add_mat((taubuf2 / taubuf1), xmat, (taubuf2 / taubuf1) * (*alpha), dxmat,
              dim, matbuf1, "No");
      add_mat(1.0, zmat, (*beta), dzmat, dim, matbuf2, "No");
      add_vec(1.0, yvec, (*beta), dyvec, num_const, vecbuf1);
    }
    /*  Test for positiveness. If not        , retry with smaller        , equal
     * steps (Line search). */
    test_pos = pd_test(matbuf1, dim) + pd_test(matbuf2, dim);
    /*  Uncomment to test only for corrector step */
    /*  if (strncmp(step_type        ,"Pre"        ,3)==0) test_pos=0; */
    if (test_pos != 0) {
      *alpha *= frac_search;
      *beta = *alpha;
    }
    if (boucl > max_reduc)
      test_pos = 0;
  }
  /*  Test for decreasing total complementarity. */
  if ((compute_mu(matbuf1, matbuf2, dim, newtau, newkappa) /
       (newtau * newtau)) > buf)
    test_step = 0;
  /*  Test for decreasing infesibility. */
  infb = compute_infeas(const_mat, const_val, num_const, matbuf1, matbuf2, cmat,
                        dim, vecbuf1, newtau, &bufa, &bufb);
  if (infa < infb)
    test_step = 0;
  /*  Finally        , update alpha and beta to be equal if conditions failed
   * , */
  /*  different if they hold. */
  /*  Uncomment to set to always equal steps (Less adventurous). */
  /*  test_step=0; */
  /*  ------------------------------------------------ */
  if (test_step == 0) {
    *alpha = -d_max(-*alpha, -*beta);
    *beta = *alpha;
  }
  if (boucl > max_reduc)
    *alpha = 0.0;
  /* Free. */
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  free_dmatrix(matbuf2, 1, dim, 1, dim);
  free_dvector(vecbuf1, 1, num_const);
  return NULL;
}

/*  This computes the absolute global error for the problem. */
double global_error(double ***const_mat, double *const_val, int num_const,
                    double **xmat, double **zmat, double **cmat, int dim,
                    double *yvec, double tau, double kappa) {
  double bufa, bufb, phib, mub, res;

  /*  Compute status. */
  phib = compute_infeas(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                        yvec, tau, &bufa, &bufb);
  mub = (compute_mu(xmat, zmat, dim, tau, kappa) * (dim + 1) - kappa * tau) /
        (tau * tau);
  res = mub;
  bufa *= (1 + sqrt(scal_vec(const_val, const_val, num_const)));
  bufb *= (1 + sqrt(scal_mat(cmat, cmat, dim)));
  if (bufa > mub)
    res = bufa;
  if (bufb > res)
    res = bufb;
  return res;
}

/*  Defines a step that only decrease primal infeasibility. */
Err take_feasibility_step(double ***const_mat, double *const_val, int num_const,
                          double **xmat, double **dcxmat, double **zmat,
                          double **dczmat, int dim, double *yvec,
                          double *dcyvec, double *dctau, double tau,
                          double *dckappa, double **inverse_gram) {
  double *vecbuf1 = dvector(1, num_const);
  double *vecbuf2 = dvector(1, num_const);
  double **matbuf1 = dmatrix(1, dim, 1, dim);
  double buf = 0;

  /*  Compute primal residual vector. */
  op_a(const_mat, num_const, xmat, dim, vecbuf1);
  add_vec(tau, const_val, -1.0, vecbuf1, num_const, vecbuf1);
  slow_im_vec(inverse_gram, vecbuf1, num_const, vecbuf2);
  op_ta(const_mat, num_const, vecbuf2, dim, dcxmat);
  zero_mat(dczmat, dim);
  zero_vec(dcyvec, num_const);
  *dctau = 0.0;
  *dckappa = 0.0;
  /*  Checking the new point. */
  add_mat(1.0, xmat, 1.0, dcxmat, dim, matbuf1, "No");
  op_a(const_mat, num_const, matbuf1, dim, vecbuf2);
  add_vec(tau, const_val, -1.0, vecbuf2, num_const, vecbuf1);
  buf = sqrt(scal_vec(vecbuf1, vecbuf1, num_const));
  /*  Free */
  free_dmatrix(matbuf1, 1, dim, 1, dim);
  free_dvector(vecbuf1, 1, num_const);
  free_dvector(vecbuf2, 1, num_const);
  return NULL;
}

/*  **************************************************************************************
 */
/*  Some routines for status tracking and display in the trace and Monitor
 * --------------- */
/*  --------------------------------------------------------------------------------------
 */
/*  **************************************************************************************
 */

/*  This stores an array in row major order form */
Err convert_array(double **in_array, int dim1, int dim2, double *out_array) {
  int i, j;
  for (i = 1; i <= dim1; i++)
    for (j = 1; j <= dim2; j++)
      out_array[(j - 1) * dim1 + (i - 1)] = in_array[i][j];
  return NULL;
}

/*  Stores and outputs basic convergence indicators */
/*  Mu is tored in log scale. */
/*  Status tracker is of dimension [1..max_iter][1..4]. */
Err status(int iter, double ***const_mat, double *const_val, int num_const,
           double **systmat, double **xmat, double **zmat, double **cmat,
           int dim, double *yvec, double tau, double kappa, double sigma,
           double eta, double alpha, double beta, double **status_tracker,
           int printlevel, int max_it, double num_error) {
  char buffer[512];
  double buferr = 0, bufphi = 0, bufcent = 0, pinf, dinf;
  static double scaleerr, scalephi, scaletau, scalepinf, scaledinf;
  int k, l;

  buferr = (compute_mu(xmat, zmat, dim, tau, kappa) * (dim + 1) - kappa * tau) /
           (tau * tau);
  bufphi = compute_infeas(const_mat, const_val, num_const, xmat, zmat, cmat,
                          dim, yvec, tau, &pinf, &dinf);
  bufcent = centrality(xmat, zmat, dim, tau, kappa);
  /*  Stores some scaling factors to plot everything in a homogenous way. */
  if (iter == 1) {
    for (k = 1; k <= max_it; k++)
      for (l = 1; l <= 11; l++)
        status_tracker[k][l] = 1e+308;
    if (buferr != 0)
      scaleerr = buferr;
    else
      scaleerr = 1.0;
    if (bufphi != 0)
      scalephi = bufphi;
    else
      scalephi = 1.0;
    if (tau != 0)
      scaletau = tau;
    else
      scaletau = 1.0;
    if (pinf != 0)
      scalepinf = pinf;
    else
      scalepinf = 1.0;
    if (dinf != 0)
      scaledinf = dinf;
    else
      scaledinf = 1.0;
  }
  /*  Scale and store verything */
  status_tracker[iter][1] =
      d_max(-12, (1 / log(10)) * log(fabs(buferr / scaleerr)));
  status_tracker[iter][2] =
      d_max(-12, (1 / log(10)) * log(fabs(bufphi / scalephi)));
  status_tracker[iter][3] =
      d_max(-12, (1 / log(10)) * log(fabs(tau / scaletau)));
  status_tracker[iter][4] = d_max(bufcent, 4);
  status_tracker[iter][5] = atan(eta);
  status_tracker[iter][6] = alpha;
  status_tracker[iter][7] = beta;
  status_tracker[iter][8] = d_max(-12, (1 / log(10)) * log(fabs(num_error)));
  status_tracker[iter][9] =
      d_max(-12, (1 / log(10)) * log(fabs(pinf / scalepinf)));
  status_tracker[iter][10] =
      d_max(-12, (1 / log(10)) * log(fabs(dinf / scaledinf)));
  status_tracker[iter][11] = d_max(
      -12, -d_max(0, (1 / log(10)) *
                         log(fabs(condition_mat(systmat, num_const + 1)))));
  /*  Output */
  if (printlevel == 3) {
    sprintf(buffer, "%d  Mu: %.2e  Infeas: %.2e \n", iter, buferr, bufphi);
    smessage(buffer);
  } else if (printlevel >= 4) {
    sprintf(buffer, "%d  Mu: %.2e  Infeas: %.2e  Hmgen: %.2e  Centr: %.2e \n",
            iter, buferr, bufphi, tau, bufcent);
    smessage(buffer);
  }

  return NULL;
}

/*  **************************************************************************************
 */
/*  Export to file and monitor for DEBUG
 * ------------------------------------------------------------- */
/*  --------------------------------------------------------------------------------------
 */
/*  **************************************************************************************
 */

void exportmat(char *fich, double **mat, int dimrow, int dimcol) {
  FILE *stream;
  int i, j;

  stream = fopen(fich, "w");
  for (i = 1; i <= dimrow; i++) {
    for (j = 1; j <= dimcol; j++) {
      fprintf(stream, " %e", mat[i][j]);
    }
    fprintf(stream, " \n");
  }
  fclose(stream);
}

void exportflatmat(char *fich, double *mat, int dim) {
  FILE *stream;
  int i, j;

  stream = fopen(fich, "w");
  for (i = 1; i <= dim; i++) {
    for (j = 1; j <= dim; j++) {
      fprintf(stream, " %e", mat[(j - 1) * dim + (i - 1)]);
    }
    fprintf(stream, " \n");
  }
  fclose(stream);
}

void visu_mat(double **mata, int nrows, int ncols, int monitorId) {
#ifdef _DEBUG
#ifdef _MONITOR
  double *flat_mat = dvector(0, ncols * nrows);
  long int monit;
  monit = (long int)(monitorId);
  convert_array(mata, nrows, ncols, flat_mat);
  SetData(&monit, flat_mat, nrows, ncols);
  free_dvector(flat_mat, 0, ncols * nrows);
#endif
#endif
}

void visu_vec(double *veca, int nrows, int monitorId) {
#ifdef _DEBUG
#ifdef _MONITOR
  double *flat_vec = dvector(0, nrows);
  int i;
  long int monit;
  monit = (long int)(monitorId);
  for (i = 1; i <= nrows; i++)
    flat_vec[i - 1] = veca[i];
  SetData(&monit, flat_vec, nrows, 1);
  free_dvector(flat_vec, 0, nrows);
#endif
#endif
}