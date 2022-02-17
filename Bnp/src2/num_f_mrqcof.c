/* ==========================================================
   FILENAME :        num_f_mrqcof.c

   PURPOSE:          for Levenberg-Marquardt
   ========================================================== */

#include "utallhdr.h"
#ifdef PVMPI
#include "parallel.h"
void __stdcall Sleep(unsigned short int dwMilliseconds);
#endif
#define NRANSI
#include "math.h"
#include "time.h"

Err mrqcof(double x[], double y[], double sig[], int ndata, double a[],
           int ia[], int ma, double **alpha, double beta[], double *chisq,
           Err (*funcs)(double, double[], double *, double[], int)) {
  Err err = NULL;
  int i, j, k, l, m, mfit = 0;
  double ymod, wt, sig2i, dy, *dyda;
  static int iteration = 0; /* for debug purpose */
#ifdef PVMPI
  static double max_time = 0;
  char *still_data_to_get;
  double *y_vect = NULL;
  double **dyda_mtrx = NULL; /*To receive the array of partial derivatives.*/
  int nb_local_data;
  int *local_indexes = NULL;
  int *rcvd = NULL; /*which rank was sent back*/
#endif
#ifdef DEBUG1
  FILE *fp;
  FILE *fp2;
  fp = fopen("d:\\temp\\partial_results.txt", "a");
  fp2 = fopen("d:\\temp\\matrix.txt", "a");
  fprintf(fp, "\n\n iteration : %d \n", iteration);
  fprintf(fp2, "\n\n iteration : %d \n", iteration++);
#endif

  for (j = 1; j <= ma; j++)
    if (ia[j])
      mfit++;
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= j; k++)
      alpha[j][k] = 0.0;
    beta[j] = 0.0;
  }
  *chisq = 0.0;

#ifdef PVMPI
  if (PRL_CONTEXT.bAggregate) {
#endif
    dyda = dvector(1, ma);
    for (i = 1; i <= ndata; i++) {
      /*compute y(xi  ,a) and dy/da_1...dy/da_ma*/
      err = (*funcs)(x[i], a, &ymod, dyda, ma);
      if (err) {
        free_dvector(dyda, 1, ma);
        return err;
      }
      sig2i = 1.0 / (sig[i] * sig[i]);
      dy = y[i] - ymod;
      /*Compute alpha  , beta and chi^2*/
      for (j = 0, l = 1; l <= ma; l++) {
        if (ia[l]) {
          wt = dyda[l] * sig2i;
          for (j++, k = 0, m = 1; m <= l; m++)
            if (ia[m])
              alpha[j][++k] += wt * dyda[m];
          beta[j] += dy * wt;
        }
      }
      *chisq += dy * dy * sig2i;
    }
    free_dvector(dyda, 1, ma);
    if (!_finite(*chisq)) {
      *chisq = HUGE_VAL;
      return "Overflow in mrqcof";
    }
#ifdef PVMPI
  } else {
    SendCalibData(x, ndata, a, ma);
    CheckRcvCalibData(0, NULL, &rcvd, (double **)NULL, (double ***)NULL,
                      (int **)NULL, &max_time, 1);
    still_data_to_get = (Err)calloc(256, sizeof(char));
    strcpy(still_data_to_get, "true");
    while (!strcmp(still_data_to_get, "true")) {
      /* Till all data chuncks have been received*/
      strcpy(still_data_to_get,
             CheckRcvCalibData(ma, &nb_local_data, &rcvd, &y_vect, &dyda_mtrx,
                               &local_indexes, &max_time, 0));
      if (strcmp(still_data_to_get, "true") &&
          strcmp(still_data_to_get, "false")) {
        /*this is an error*/
        /*This level deallocation done in CheckRcvCalibData*/
        smessage("Error on remote machine:");
        smessage(still_data_to_get);
        smessage("aborting calibration...");
        return still_data_to_get;
      }
      if (nb_local_data == 0)
        continue;
      /* Compute alpha  , beta  , chi square */
      for (i = 0; i < nb_local_data; i++) {
        ymod = y_vect[i];
        dyda = dyda_mtrx[i] -
               1; /* NR algorithms are for indexes starting at 1... */
                  /* This is definitely not optimal
                  2*  , 1/ & 2 [] where only 1 [] is necessary (ndata*niter times)
sig should be defined in mrqmin by sig[i] = 1.0/(sig[i]*sig[i])*/
        sig2i = 1.0 / (sig[local_indexes[i] + 1] * sig[local_indexes[i] + 1]);
        dy = y[local_indexes[i] + 1] - ymod;
        /* Compute alpha  , beta and chi^2 */
        for (j = 0, l = 1; l <= ma; l++) {
          if (ia[l]) {
            wt = dyda[l] * sig2i;
            for (j++, k = 0, m = 1; m <= l; m++)
              if (ia[m])
                alpha[j][++k] += wt * dyda[m];
            beta[j] += dy * wt;
          }
        }
        *chisq += dy * dy * sig2i;
      }
      Sleep(10); /*allow other tasks to do some job*/
    }

    /* free memory */
    TagTimes(max_time);
    free(still_data_to_get);
    free(y_vect);
    for (i = 0; i < nb_local_data; i++)
      free(dyda_mtrx[i]);
    free(dyda_mtrx);
  } /* end of if (Version...) */
#endif
  for (j = 2; j <= mfit; j++)
    for (k = 1; k < j; k++)
      alpha[k][j] = alpha[j][k];
  return NULL;
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0>)"?. */
