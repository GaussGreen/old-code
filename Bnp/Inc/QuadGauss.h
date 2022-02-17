#ifndef QUADGAUSS
#define QUADGAUSS

#include "math.h"
#include "srt_h_all.h"

typedef struct {
  /* Today */
  long today;

  /* TS Dates and Times */
  int num_dates;
  long *dates;
  double *times;

  /* Number of factors */
  int num_factor;

  /* X1 volatility */
  double *sigma1;
  /* X1 weight */
  double *e1_1;
  /* X1^2 weight */
  double *e1_2;

  /* X2 volatility */
  double *sigma2;
  /* X2 weight */
  double *e2_1;
  /* X2^2 weight */
  double *e2_2;

  /* X1 / X2 correlation */
  double *rho;

} qg_str, *QG_STR;

Err qg_free_struct(qg_str *qgstr);

Err qg_free_und_struct(SrtUndPtr pUndDesc);

Err qg_get_struct_from_und(char *und, qg_str **qgstr);

Err qg_init_struct(long today,

                   /*	TS Dates	*/
                   int num_dates, long *dates,

                   int num_factor,

                   /*	X1 TS	*/
                   double *sigma1, double *e1_1, double *e1_2,

                   /*	X2 TS	*/
                   double *sigma2, double *e2_1, double *e2_2,

                   /*	Correl TS	*/
                   double *rho,

                   /*	Structure	*/
                   qg_str *qgstr);

Err SrtInitQGUnd(/* und name */
                 char *undName,

                 /* mkt name */
                 char *ycname,

                 /*	TS Dates	*/
                 int num_dates, long *dates,

                 int num_factor,

                 /*	X1 TS	*/
                 double *sigma1, double *e1_1, double *e1_2,

                 /*	X2 TS	*/
                 double *sigma2, double *e2_1, double *e2_2,

                 /*	Correl TS	*/
                 double *rho);

Err qg_get_ts_from_und(/* qg und name */
                       char *undName,

                       /* TS Dates and Times */
                       int *num_dates, long **dates, double **times,

                       int *num_factor,

                       /*	X1 TS	*/
                       double **sigma1, double **e1_1, double **e1_2,

                       /*	X2 TS	*/
                       double **sigma2, double **e2_1, double **e2_2,

                       /*	Correl TS	*/
                       double **rho);

#endif
