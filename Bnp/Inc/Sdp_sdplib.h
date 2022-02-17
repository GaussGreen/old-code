/*  Author: AdA */
/*  Date:   22/06/1999 */

/*  Associated header: sdplib.h */
/*  Uses: some nr routines and some BLAS routines.  */

/*  BLAS from the Intel Math Kernel library can be linked directly  ,  */
/*  or the MiniBLAS provided for use in griffin can be used instead. */
/*  See sdp.h for BLAS library choice. */

/*  Comment CBLAS to use the autonomous version. */
/* #define CBLAS */
/*  ------------ */

/*  -----------------------------------------------------------------------------------------------
 */
/*  This is the main toolbox used by sdp */
/*  -----------------------------------------------------------------------------------------------
 */

#ifndef _SDPLIB_H
#define _SDPLIB_H

/*  Header for Intel MKL CBLAS */
#ifdef CBLAS
#include "cblas.h"
#endif
/*  --------------------------
#ifdef _DEBUG
#include "MonitorInterface.h"
#endif  */
//#include "math.h"
//#include "time.h"
//#include "utallhdr.h"
#include "num_h_diagonalise.h"
#include "num_h_gaussj.h"
#include "sdp_more_nr.h"

/*  These are algo options. ---------------------------------------------- */
/*  Uncomment for fast "on demand" method for inv(E) kronecker operator */
#define fast_method
/*  ---------------------------------------------------------------------- */
#ifdef CBLAS
/*  This is slower and not really more precise. Works only when using CBLAS. */
/* #define very_slow */
/*  ---------------------------------------------------------------------- */
#endif

/*  Copy A into B */
Err copy_mat(double **mata, int dim, double **matres);

/*  zeros a square matrix */
Err zero_mat(double **mata, int dim);

/*  zeros a vector */
Err zero_vec(double *veca, int dim);

/*  Stack a matrix of dimension dim into a [0  ,...  ,(dim*dim-1);] vector. */
Err stack_mat(double **mata, int dim, double *vecres);

/*  Adds aplha*A+beta*B put it into C (if transpose="T"  , B is transposed); */
Err add_mat(double alpha, double **mata, double beta, double **matb, int dim,
            double **matc, char *transpose);

/*  Max of two doubles */
double d_max(double a, double b);

/*  Simple product of standard [1..n][1..n] matrixes. */
/*  All matrixes not symmetric. */
/*  result CANNOT be same matrix as mata or matb. (Speeds-up) */
/*  It also computes B*A instead of A*B because it proves to */
/*  improve the precision of the algo in the Lyapounov equation. */
Err mat_mult_nobuf(double **mata, double **matb, int dim, double **result);

/*  Simple product of standard [1..n][1..n] matrixes. */
/*  All matrixes not symmetric. */
Err mat_mult_buf(double **mata, double **matb, int dim, double **result);

/*  Transpose a sqaure matrix */
Err transp(double **mata, int dim, double **result);

/*  Fast image of a standard [1..n][1..n] matrix by a big operator
 * [1..n^2][1..n^2] on matrixes. */
Err fast_imat(double **opera, double **mata, int dim, double **result);

/*  Defines the A(*);B tensor product of matrixes of size dim (A  , B  , must be
 * symetric);. */
/*  Returns a matrix of dimension dim^2 */
/*  This is the most intensive computation in the algorithm. */
/*  (ToDo); Look for optimized formulation. */
Err tensor_prod(double **mata, double **matb, int dim, double **result);

/*  Computes (A(*);B);(O);  , with A  , B  , O symmetric. */
/*  Returns a matrix of dimension dim */
Err tensor_oper(double **mata, double **matb, double **mato, int dim,
                double **result);

/*  Scalar product of two matrixes */
double scal_mat(double **mata, double **matb, int dim);

/*  Scalar product of two vectors. */
double scal_vec(double *veca, double *vecb, int dim);

/*  Adds two scaled vectors alpha*VECA + beta*VECB = result */
Err add_vec(double alpha, double *veca, double beta, double *vecb, int dim,
            double *vecres);

/*  Image of a small vector by a matrix */
Err slow_im_vec(double **mata, double *veca, int dim, double *vecres);

/*  Test if a matrix is Positive Definite using Cholesky decomposition. */
/*  Returns 0 if PD  , or the size of the first  */
/*  non PD submatrix. */
int pd_test(double **mata, int dim);

/*  Solves a linear system A.x=b  , returns inv(A) and x. if algo_type is 0  , x
 * otherwise. */
/*  Algo_type gives the choice between 0:Gauss elimination  , 1:Gauss Siedel. */
Err lsolve(double **mata, double *vecb, int dim, double *vecres,
           double **matres, int algo_type, double prec_obj, int num_improve);

/*  Returns the inverse of a matrix */
Err inv_mat(double **mata, int dim, double **result);

/*  Computes the minimum eigenvalue of a symmetric matrix. */
Err min_eigenval(double **mata, int dim, double *result, double *condition);

/*  Computes the minimum eigenvalue of a non-symmetric matrix. */
Err min_eigenval_nosym(double **mata, int dim, double *result, int max1_min0);

/*  diagonalisation of a symmetric matrix. */
/*  Everythinng in [1..dim] format. */
Err diago(double **mata, int dim, double **eigenvec, double *eigenval);

/*  Computes the condition number of a matrix. */
double condition_mat(double **mata, int dim);

/*  Adds 5.0*(max eigenval)*machine_epsilon to the eigenvalues. */
Err fudge_eigenvals(double *evals, int dim);

/*  *************************************** */
/*  For the msclib Gauss-Siedel solver. --- */
/*  *************************************** */

/*  WARNING: vectors are handled as [1..dim] but should be allocated [0..dim]
 * !!!!!!!!!!!!!!! */
/*  (ToDo) Correct this shit I introduced. */
void solve_linear(double **A, int row_start, int col_start, int rows, int cols,
                  double *B, double *X, double prec);

/*  *******************************************************************************
 */
/*  Part Two: SDP specific routines.
 * ---------------------------------------------- */
/*  -------------------------------------------------------------------------------
 */
/*  *******************************************************************************
 */

/*  Computes the primal constraints operator (Tr(A1.X);  , ...  ,Tr(Ai.X);  ,
 * ....); operating on a matrix mata. */
Err op_a(double ***con_mat, int num_const, double **mata, int dim,
         double *vecres);

/*  Computes the augmented constraints operator (A  ,-C);. */
/*  Returns a dim+1 vector. */
Err op_a_augm(double ***con_mat, int num_const, double **cmat, double **mata,
              int dim, double *vecres);

/*  Computes the transposition of the constraints operator operating on a vector
 * veca. */
Err op_ta(double ***con_mat, int num_const, double *veca, int dim,
          double **matres);

/*  Computes the transposition of the augmented constraints operator operating
 * on a vector veca. */
Err op_ta_augm(double ***con_mat, int num_const, double **cmat, double *veca,
               int dim, double **matres);

/*  Computes the basic step length. */
double step_length(double **mata, double **dmata, int dim, double tau,
                   double dtau, double kappa, double dkappa, double gamma);

/*  This computes phi  , the infeasibility measure. */
double compute_infeas(double ***const_mat, double *const_val, int num_const,
                      double **xmat, double **zmat, double **cmat, int dim,
                      double *yvec, double tau, double *prim, double *dual);

/*  This computes sigma */
double compute_sigma(double **xmat, double **dxmat, double **zmat,
                     double **dzmat, int dim, double alpha, double beta,
                     double tau, double dtau, double kappa, double dkappa,
                     double expon);

/*  Computes mu. */
double compute_mu(double **xmat, double **zmat, int dim, double tau,
                  double kappa);

/*  This computes the sytem A.inv(E);.F.At matrix. */
/*  (ToDo); Choice between AHO and HKM search directions. */
/*  It returns also inv(E); for further usage. */
/*  Puts the result in a [1..(dim+1);][1..(dim+1);] matrix. */
Err system_matrix(double ***const_mat, int num_const, double **xmat,
                  double **zmat, double **cmat, int dim, char *srch,
                  double **result, double **big_inv_E, double *errtrack,
                  int num_improve);

/*  This computes the sytem A.At matrix. */
Err compute_gram(double ***const_mat, int num_const, int dim, double **result);

/*  This computes the left hand final sytem matrix. */
Err left_matrix(double **matm, int num_const, double *bvec, double kappa,
                double tau, double **result);

/*  Computes the Rd matrix. */
Err Rd_mat(double ***const_mat, int num_const, double **zmat, double **cmat,
           int dim, double *yvec, double tau, double **result);

/*  Computes the rp [1..dim+1] vector. */
Err rp_vec(double ***const_mat, int num_const, double **xmat, double **cmat,
           int dim, double *yvec, double *bvec, double tau, double kappa,
           double *vecres);

/*  Computes the rc scalar */
double rc_val(double sigma, double mu, double tau, double kappa);

/*  Computes the Rc matrix. */
Err Rc_mat(double **xmat, double **zmat, int dim, double sigma, double mu,
           char *srch, double **result);

/*  Computes the Rq matrix. */
Err Rq_mat(double **xmat, double **dxmat, double **zmat, double **dzmat,
           int dim, double sigma, double mu, char *srch, double **result);

/*  This computes the right hand side of the (dy  ,dtau); system equation. */
/*  It uses the big_inv_E computed by system_matrix. */
Err right_vector(double ***const_mat, int num_const, double **big_inv_E,
                 double **xmat, double **zmat, double **cmat, double **dpxmat,
                 double **dpzmat, char *step_type, int dim, double *yvec,
                 double *bvec, double eta, double tau, double sigma, double mu,
                 double kappa, double dptau, double dpkappa, double *vecres,
                 int num_improve);

/*  Computes dZ. */
Err delta_z(double ***const_mat, int num_const, double **zmat, double **cmat,
            int dim, double *yvec, double *dyvec, double tau, double eta,
            double dtau, double **result);

/*  Computes dX */
Err delta_x(double **big_inv_E, double **xmat, double **dpxmat, double **zmat,
            double **dpzmat, double **dczmat, int dim, double sigma, double mu,
            char *pred_corr, double **result, int num_improve);

/*  Creates an initial scaled iterate for X and Z. */
Err init_sol(double ***const_mat, int num_const, double *const_val,
             double **cmat, int dim, double **startx, double **startz);

/*  Computes the centrality measure.  (1-min_eigenval(X.Z)/mu) */
double centrality(double **xmat, double **zmat, int dim, double tau,
                  double kappa);

/*  Tests if steps actually go towards global progress. */
/*  Set them both to the minimum if not. */
Err set_steps(double ***const_mat, double *const_val, int num_const,
              double *yvec, double *dyvec, double **xmat, double **dxmat,
              double **zmat, double **dzmat, double **cmat, int dim,
              double kappa, double dkappa, double tau, double dtau, double eta,
              double *alpha, double *beta, double frac_search, char *step_type);

/*  This computes the inv(E) operator image of a matrix with the fast Kronecker
 * product algo.  */
/*  P and Z are supposed to be commuting symmetric. (It is the case for AHO and
 * HKM). */
Err fast_oper_inv_E(double **pmat, double **zmat, double **mata, int dim,
                    double **result, double *error_ratio, int num_improve,
                    double **bufmat2, double *bufvec2, int diag_switch);

/*  This Computes the global absolute error. */
double global_error(double ***const_mat, double *const_val, int num_const,
                    double **xmat, double **zmat, double **cmat, int dim,
                    double *yvec, double tau, double kappa);

/*  Defines a step that only decrease primal infeasibility. */
Err take_feasibility_step(double ***const_mat, double *const_val, int num_const,
                          double **xmat, double **dcxmat, double **zmat,
                          double **dczmat, int dim, double *yvec,
                          double *dcyvec, double *dctau, double tau,
                          double *dckappa, double **inverse_gram);

/*  **************************************************************************************
 */
/*  Some routines for status display in the trace
 * ---------------------------------------- */
/*  --------------------------------------------------------------------------------------
 */
/*  **************************************************************************************
 */

/*  This stores an array in row major order form */
Err convert_array(double **in_array, int dim1, int dim2, double *out_array);

/*  Sends the status to the trace window and keep it in the status_tracker. */
/*  Status tracker is of dimension [1..max_iter][1..4]. */
Err status(int iter, double ***const_mat, double *const_val, int num_const,
           double **systmat, double **xmat, double **zmat, double **cmat,
           int dim, double *yvec, double tau, double kappa, double sigma,
           double eta, double alpha, double beta, double **status_tracker,
           int printlevel, int max_it, double num_error);

/*  **************************************************************************************
 */
/*  Export to file for DEBUG
 * ------------------------------------------------------------- */
/*  --------------------------------------------------------------------------------------
 */
/*  **************************************************************************************
 */

void exportmat(char *fich, double **mat, int dimrow, int dimcol);

void visu_mat(double **mata, int nrows, int ncols, int mointorId);

void visu_vec(double *veca, int nrows, int monitorId);

#endif