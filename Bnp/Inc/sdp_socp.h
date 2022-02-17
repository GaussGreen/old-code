/*
Feb 2001. LTQLM Inc.
Wrapper Author: AdA (Alex d'Aspremont      , aspremon@stanford.edu      , +1 415
864 2962) Base code: Miguel Sousa Lobo      , Lieven Vandenberghe      , and
Stephen Boyd. See code website: www.stanford.edu/~boyd/SOCP.html

WARNING: Matrixes are stored in flat BLAS-LAPACK (FORTRAN) format      , column
major order      , arrays starting at 0      , the function ijtok(i      ,j
,number_of_rows) makes the conversion.

This solves a quadratic optimisation program under LP constraints using an SOCP
solver. See header for source description.
*/

/******************************************************************/
/*	This is the header for the SDP_SOCP (quadratic programming functions
        It can optimize a quadratic function under LP constraints as
        a particular case (Smoothness under Bid-Ask for example)
        The original code was made by M Lobo      , see comments below.
        See code website: www.stanford.edu/~boyd/SOCP.html */
/*
******** Don't forget to link with the lapack and F77 libs... *****************
*/

/************************************************************************/
/*				WRAPPER CODE (Phase One and Two) */
/************************************************************************/

// Set up a demo vol surface smoothing problem
Err socp_smoothvol(int **index_set, double *upperbounds, double *lowerbounds,
                   int numpoints, int num_dates, double centering_weight,
                   double **result, int printlevel);

// Set up a demo Bid_Ask problem
Err socp_bidask(double *upperbounds, double *lowerbounds, int numpoints,
                double centering_weight, double *result, int printlevel);

//	Main function for smoothing under LP constraints      ,
//	solves the following program:
//
//		Minimize ||Ax||^2
//		subject to M.x+b>=0 where M is a (dim      ,numconst) matrix
//
//	"center" is the l-inf center of the inequalities and alpha>0 avoids
//	possible degenerate solutions (parallel lines).
//	"center" is computed by the phase one problem.
Err socp_smooth_lp(double **objective_matrix, double **constraints_matrix,
                   double *constraints_value, int nobjrows, int dim,
                   int numconst, double alpha, double *result, int printlevel);

// Looks for the center of the inequality by solving phase one.
// Wrapper for the core phase one function "phase_one_go"
Err phase_one(double **objective_matrix, double **constraints_matrix,
              double *constraints_value, int nobjrows, int dim, int numconst,
              double *centersol, double *dualsol, int printlevel);

// Looks for the center of the inequality by solving phase one.
Err phase_one_go(double **objective_matrix, double **constraints_matrix,
                 double *constraints_value, int nobjrows, int dim, int numconst,
                 double *centersol, double *dualsol, int *vecdimconstin,
                 int ndbl, int nint, int mhist, int nhist, int iter,
                 int out_mode, int printlevel);

// Computes a dual feasible starting point for phase one (big_M method):
// z=inv(A)(f-c2)      , ....
Err feasible_dual_vecs_one(double *fvecin, double *amatin, int nobjrows,
                           int dim, int numconst, double *zvar);

// Computes a dual feasible starting point for phase two (big_M method).
// z=inv(A)(f-c2)      , ....
Err feasible_dual_vecs_two(double *fvecin, double *amatin, int nobjrows,
                           int dim, int numconst, double *zvar);

// Euclidian norm of a vector
double normvec(double *vec, int ilow, int ihigh);

// Looks for the smoothest solution inside the inequalitites.
// Wrapper for the core phase one function "phase_two_go"
Err phase_two(double **objective_matrix, double **constraints_matrix,
              double *constraints_value, int nobjrows, int dim, int numconst,
              double *centersol, double *dualsol, int printlevel);

// Looks for the smoothest solution inside the inequalitites.
Err phase_two_go(double **objective_matrix, double **constraints_matrix,
                 double *constraints_value, int nobjrows, int dim, int numconst,
                 double *centersol, double *dualsol, int *vecdimconstin,
                 int ndbl, int nint, int mhist, int nhist, int iter,
                 int out_mode, int printlevel);

// converts between (i      ,j) matrix storage and flat LAPACK matrix storage ,
// FORTRAN style      , column major order.
int ijtok(int i, int j, int numcol);

// Converts between upper triangular storage index and flat storage.
int toflat(int i, int j);

// Fills the matrix that turns the dim vector x into the dim-1 vector dx.
Err fill_gradient_matrix(double **gradmat, int dimcol, int rowstofill);

// Checks for strict feasibility... (Possible trouble if initial point is to
// close to infeasible).
int check_feasibility(double *amatin, double *fvecin, int nobjrows, int dim,
                      int numconst, double *zvar, int *vecdimconstin,
                      double *bvecin, double *xvar);

/************************************************************************/
/*				ORIGINAL CODE
*
/************************************************************************/

/*
 *   sdp_socp.h
 *
 *   Second-Order Cone Programming
 *   Header file for socp.c and socp_mex.c
 *
 *   mlobo@isl.stanford.edu -- 96/97
 */

/*
Copyright c 1997 by Miguel Sousa Lobo      , Lieven Vandenberghe      , and
Stephen Boyd. Permis- sion to use      , copy      , modify      , and
distribute this software for any purpose without fee is hereby granted      ,
provided that this entire notice is included in all copies of any software which
is or includes a copy or modifcation of this software and in all copies of the
supporting doc- umentation for such software. This software is being provided
"as is"      , without any express or implied warranty. In particular      , the
authors do not make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any particular purpose.
*/

/* basic macros */

//#define SOCPSQR(x)    ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

/* constants for socp algorithm */

#define MAX_ITER_PLANE 20
#define MAX_LAMBDA2 1e-2 /* plane search stopping crit. */
#define DIV_ALPHA 2
#define MIN_ALPHA 1e-6 /* max. of about 20 line search iterations */

#ifdef nounderscores
#define ddot_ ddot
#define dcopy_ dcopy
#define daxpy_ daxpy
#define dscal_ dscal
#define dgemv_ dgemv
#define dsyr_ dsyr
#define dsyrk_ dsyrk
#define dposvx_ dposvx
#define dgelss_ dgelss
#endif

/* BLAS 1 */
double ddot_();
void dcopy_();
void daxpy_();
void dscal_();

/* BLAS 2 */
void dgemv_();
void dsyr_();

/* BLAS 3 */
void dsyrk_();

/* LAPACK */
void dposvx_();
void dgelss_();
void dgels_();

/* socp.c */

void socp_getwork(
    /* input args.:  problem dimensions and max. num. of iterations */
    int L, int *N, int n, int max_iter, int out_mode,
    /* output args.:  dimensions of history output matrix      , */
    /*   number of doubles      , pointers and ints required for workspace */
    int *mhist, int *nhist, int *ndbl, int *nint);

int socp(
    /* problem dimensions */
    int L, int *N, int n,
    /* problem data */
    double *f, double *A, double *b,
    /* in:  initial primal and dual strictly feasible points */
    /* out: final points */
    double *x, double *z,
    /* stopping criteria      , */
    /*   *iter on entry is max. number of iterations      , */
    /*   on exit actual number performed */
    double abs_tol, double rel_tol, double target, int *iter,
    /* algorithm parameter */
    double Nu,
    /* reason for exit      , output matrix with extra info */
    /* out_mode specifies what will be stored in *hist */
    int *info, int out_mode, double *hist,
    /* workspace      , use socp_getwork() to determine required sizes */
    double *dblwork, int *intwork);
