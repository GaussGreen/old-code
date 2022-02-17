/* ==============================================
FILENAME:  num_h_egarch.h

  PURPOSE:   Egarch functions
   ============================================== */
#ifndef NUM_H_EGARCH_H
#define NUM_H_EGARCH_H

#define DEFAULT_S2_INIT 0.01
#define DEFAULT_BETA_0 -0.01
#define DEFAULT_BETA_1 -0.05
#define DEFAULT_BETA_2 -0.05

typedef struct garchdata {
  Ddate date;
  double rate;
  double s1;
} GarchData;

Err srt_f_egarchmethod(
    Ddate date_chosen, /* date chosen to compute the GARCH */
    Ddate *dates,      /* input dates */
    double *rates,     /* input rates */
    double *s1,        /* input s1 */
    long num_data,     /* number of input data */
    long *fixed,       /* array describing if param are fixed or not */
    double *s2_init,   /* output & input Initialisation of S2 */
    double *beta_0,    /* output & input beta_0 in the Garch algorithm */
    double *beta_1,    /* output & input beta_0 in the Garch algorithm */
    double *beta_2,    /* output & input beta_0 in the Garch algorithm */
    double *s2_last,   /* output last value of S2 */
    double *zeta_last, /* output last value of zeta */
    long num_iter      /* number of iterations in the algorithm */
);

#endif
