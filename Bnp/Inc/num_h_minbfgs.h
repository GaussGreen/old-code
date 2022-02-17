/* -------------------------------------------------------------

   AUTHOR: E. FOURNIE

   DATE : MAY 99

   FILENAME: num_h_minbfgs.h

   PURPOSE:	Quasi-Newton and Conjugate gradient as in numerical recipes

   MODIFICATION: JAN 00

 ----------------------------------------------------------------------- */

/*
        Quasi-newton DFP or BFGS :
        - BFGS is usually better
*/

typedef enum {
  ZEROGRADIENTCRITERIUM = 0,
  FUNCCRITERIUM = 1,
  DEFAULTCRITERIUM = 0
} DFPMinCriteriumType;

Err dfpmin(DFPMinCriteriumType returnCriterium, double p[], int n,
           double func_gtol, double dfunc_tol, int *iter, double *fret,
           Err (*func)(double[], double *), Err (*dfunc)(double[], double[]));
/*
        Conjugate gradient Fletcher-Reeves or Polak Ribiere :
    it can be switched between the two methods
    - Polak Ribiere is proven to always converge
    - Fletcher-reeves is usually the fastest (but can fail to converge)
*/

Err frprmn(double p[], int n, double ftol, int *iter, double *fret,
           Err (*func)(double[], double *), Err (*dfunc)(double[], double[]));

/* ==================================================================== */
