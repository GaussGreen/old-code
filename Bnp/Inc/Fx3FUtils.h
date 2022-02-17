
#ifndef Fx3FUtils_h
#define Fx3FUtils_h

#include "srt_h_all.h"

#define NSTD 4
#define MAXNODE3D 400

#define MESH 1.732050808
#define HALF_MESH 0.866025404
#define DOUBLE_MESH 3.464101615
#define MESH_SQUARE 3
#define DOUBLE_MESH_SQUARE 6

long Get_Index(double T, double *Maturity, long nbrMat);
long Get_Index_Rec(double T, double *Maturity, long nbrMat, long StartIndex);

/*	Product of two matrices */
#define prod_matrix(left, right, res, i, j, k, sum)                            \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      for (j = 0; j < 3; j++) {                                                \
        sum = 0.0;                                                             \
        for (k = 0; k < 3; k++) {                                              \
          sum += left[i][k] * right[k][j];                                     \
        }                                                                      \
        res[i][j] = sum;                                                       \
      }                                                                        \
    }                                                                          \
  }

/*	Product matrix x vector */
#define prod_matrix_vector(matrix, vector, res, i, k, sum)                     \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      sum = 0.0;                                                               \
      for (k = 0; k < 3; k++) {                                                \
        sum += matrix[i][k] * vector[k];                                       \
      }                                                                        \
      res[i] = sum;                                                            \
    }                                                                          \
  }

/*	Transpose */
#define transpo_matrix(matrix, res, i, j)                                      \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      res[i][i] = matrix[i][i];                                                \
      for (j = i + 1; j < 3; j++) {                                            \
        res[i][j] = matrix[j][i];                                              \
        res[j][i] = matrix[i][j];                                              \
      }                                                                        \
    }                                                                          \
  }

/*	Product matrix x constant */
#define prod_matrix_constant(left, right, cons, res, i, j, k, sum)             \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      for (j = 0; j < 3; j++) {                                                \
        sum = 0.0;                                                             \
        for (k = 0; k < 3; k++) {                                              \
          sum += left[i][k] * right[k][j];                                     \
        }                                                                      \
        res[i][j] = sum * cons;                                                \
      }                                                                        \
    }                                                                          \
  }

/*	Product matrix x vector x constant */
#define prod_matrix_vector_constant(matrix, vector, cons, res, i, k, sum)      \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      sum = 0.0;                                                               \
      for (k = 0; k < 3; k++) {                                                \
        sum += matrix[i][k] * vector[k];                                       \
      }                                                                        \
      res[i] = sum * cons;                                                     \
    }                                                                          \
  }

/*  Product matrix x diagonal matrix */
#define prod_matrix_diago(matrix, diago, res, i, j, sum)                       \
  {                                                                            \
    for (j = 0; j < 3; j++) {                                                  \
      sum = diago[j];                                                          \
      for (i = 0; i < 3; i++) {                                                \
        res[i][j] = sum * matrix[i][j];                                        \
      }                                                                        \
    }                                                                          \
  }

/*	Product diagonal matrix x matrix */
#define prod_diago_matrix(matrix, diago, res, i, j, sum)                       \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      sum = diago[i];                                                          \
      for (j = 0; j < 3; j++) {                                                \
        res[i][j] = sum * matrix[i][j];                                        \
      }                                                                        \
    }                                                                          \
  }

/*	Product of two matrices 3x3 * 3x8 */
#define prod_matrix_special(left, right, res, i, j, k, sum)                    \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      for (j = 0; j < 8; j++) {                                                \
        sum = 0.0;                                                             \
        for (k = 0; k < 3; k++) {                                              \
          sum += left[i][k] * right[k][j];                                     \
        }                                                                      \
        res[i][j] = sum;                                                       \
      }                                                                        \
    }                                                                          \
  }

/*	Product of two matrices 3x3 * 3x8 times constant */
#define prod_matrix_constant_special(left, right, cons, res, i, j, k, sum)     \
  {                                                                            \
    for (i = 0; i < 3; i++) {                                                  \
      for (j = 0; j < 8; j++) {                                                \
        sum = 0.0;                                                             \
        for (k = 0; k < 3; k++) {                                              \
          sum += left[i][k] * right[k][j];                                     \
        }                                                                      \
        res[i][j] = sum * cons;                                                \
      }                                                                        \
    }                                                                          \
  }

Err get_transf_matrix(double *sigma, double *correl, double **transf_matrix,
                      double **transf_matrix_inv, double **covar_matrix,
                      double **eigen_vector, double **eigen_vector_transp,
                      double *eigen_val);

Err get_lambda_from_ir_ts(TermStruct *ts, double *lambda);

Err compute_vol_times(char *und3dfx, int *num_vol_times, double **vol_times,
                      double last_time);

Err fill_time_vector(double **time, int *nstp, int num_bar_times,
                     double *bar_times, int num_vol_times, double *vol_times,
                     int target_nstp);

Err find_phi(double *date, double *phi, int nstp, double theDate,
             double *phi_at_t);

/*	Adjust empirical covariance */
void adjust_emp_covar(double **g, int npth, int n);

/*	Generate BalSam numbers */
Err balsam_generation(long nbPaths, /* Must be odd */
                      long nbSteps, double **matrix);

Err correlate_variable_corr(double *dom_std, double *for_std, double *fx_std,
                            double *corr_dom_for, double *corr_dom_fx,
                            double *corr_for_fx, long nbPaths, long nbSteps,
                            double **matrix);

Err correlate_variable_covar(double *dom_std, double *for_std, double *fx_std,
                             double *covar_dom_for, double *covar_dom_fx,
                             double *covar_for_fx, long nbPaths, long nbSteps,
                             double **matrix);
#endif