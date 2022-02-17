/* ------------------------------------------------------------------------
        FILENAME: 	num_f_levenberg.c

        FUNCTION:   levenberg_marquardt

        PURPOSE:	provide a full function that minimises an error
                                criteria using the Levenberg-Marquardt
   algorithm. This function is intended as a generic wrapper to the mrqmin
   function as described in Numerical Recipes in C. INPUTS: data:     a vector
   of data (points of evaluation) target:   a vector of target for the function
                                weight:   the relative weight of the associated
   data ndata:    number of data param:    a vector of parameters (they will be
   optimised) nparam:   number of parameters used in the optimisation niter:
   number of successfull iterations funcs():  the f(x  ,a) function as described
   in mrqmin

                                chisq:    the error criteria after optimisation

   ------------------------------------------------------------------------ */

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_levenberg.h"

#define free_levenberg_memory                                                  \
  {                                                                            \
    if (covar)                                                                 \
      free_dmatrix(covar, 1, nparam, 1, nparam);                               \
    covar = NULL;                                                              \
    if (alpha)                                                                 \
      free_dmatrix(alpha, 1, nparam, 1, nparam);                               \
    alpha = NULL;                                                              \
  }

Err levenberg_marquardt(double *data,              /* From [1] to [ndata] */
                        double *target,            /* From [1] to [ndata] */
                        double *weight,            /* From [1] to [ndata] */
                        long ndata, double *param, /* From [1] to [nparam] */
                        long nparam, long niter,
                        Err (*funcs)(double, double[], double *, double[], int),
                        double *chisq) {
  Err err = NULL;
  int i;
  long *use_param = NULL;

  /* Sets use_param to 1: optimisation is done on all */
  use_param = lngvector(1, nparam);
  for (i = 1; i <= nparam; i++)
    use_param[i] = 1;

  err = levenberg_marquardt_select(data,         /* From [1] to [ndata] */
                                   target,       /* From [1] to [ndata] */
                                   weight,       /* From [1] to [ndata] */
                                   ndata, param, /* From [1] to [nparam] */
                                   use_param,    /* From [1] to [nparam] */
                                   nparam, niter, funcs, chisq);

  /* Free what has been allocated */
  free_lngvector(use_param, 1, nparam);

  /* Return a success message */
  return err;
}

/* --------------------------------------------------------------------- */

Err levenberg_marquardt_select(
    double *data,              /* From [1] to [ndata] */
    double *target,            /* From [1] to [ndata] */
    double *weight,            /* From [1] to [ndata] */
    long ndata, double *param, /* From [1] to [nparam] */
    long *use_param,           /* From [1] to [nparam] */
    long nparam, long niter,
    Err (*funcs)(double, double[], double *, double[], int), double *chisq) {
  Err err = NULL;
  double **covar;
  double **alpha;
  double alambda, alambda_old;
  double chisq_old;
  long success_iter;
  long true_iter;
  long failed_iter = 0;
  long finish_chisq = 0;

/* Warn the user: Levenberg-Marquardt algorithm */
#ifdef _DEBUG
  smessage(" initialising...");
  smessage(" please wait...");
  smessage("");
#endif _DEBUG

  /* Defines a covar and an alpha matrixes */
  covar = dmatrix(1, nparam, 1, nparam);
  alpha = dmatrix(1, nparam, 1, nparam);

  /* Initialises parameters for search */
  alambda = -1.0;
  err = mrqmin(data, target, weight, (int)ndata, param, (int *)use_param,
               (int)nparam, covar, alpha, chisq, funcs, &alambda);
  if (err) {
    free_levenberg_memory;
    return serror("Error in LM minimisation: %s", err);
  }
#ifdef _DEBUG
  smessage(" starting minimisation  , criteria: %.2f", 100 * sqrt(*chisq));
#endif _DEBUG

  /* Starts iterations */
  success_iter = 1;
  true_iter = 1;
  /*	while (success_iter <= niter)
   */
  while ((finish_chisq <= 4) && (success_iter <= niter) && (failed_iter <= 8)) {
    alambda_old = alambda;
    chisq_old = *chisq;
#ifdef _DEBUG
    smessage(" start iteration: %d", success_iter);
#endif _DEBUG
    mrqmin(data, target, weight, (int)ndata, param, (int *)use_param,
           (int)nparam, covar, alpha, chisq, funcs, &alambda);
    if (err) {
      free_levenberg_memory;
      return serror("Error in LM minimisation: %s", err);
    }

#ifdef _DEBUG
    smessage(" end iteration: %d", success_iter);
    smessage(" criteria: %.2f", 100 * sqrt(*chisq));
#endif _DEBUG

    /* Check if the search has been succesful (lambda has reduced) */
    if (alambda <= alambda_old) {
      success_iter++;
      failed_iter = 0;
#ifdef _DEBUG
      smessage(" iteration succeeded");
      smessage("");
#endif _DEBUG
      /*			if (true_iter >= max_true_iter)
                                      break;
      */
      if (*chisq < 1e-16)
        finish_chisq++;
/*			else
			if (chisq_old - (*chisq) < 0.0001)
				finish_chisq++;
			if (finish_chisq >= 4)
				break;
*/		}
		else 
		{
#ifdef _DEBUG
  smessage(" iteration %d failed: try again", success_iter);
  smessage("");
#endif _DEBUG
  failed_iter++;
  if (*chisq < 1e-16)
    finish_chisq++;
}
true_iter++;
  }

  /* Make the last iteration with quadratic approximation: lambda = 0.0 */
  alambda = 0.0;
  mrqmin(data, target, weight, (int)ndata, param, (int *)use_param, (int)nparam,
         covar, alpha, chisq, funcs, &alambda);
  if (err) {
    free_levenberg_memory;
    return serror("Error in LM minimisation: %s", err);
  }

  /* Free what has been allocated */
  free_levenberg_memory;

  /* Return a success message */
  return NULL;
}
