/*--------------------------------------------------------------
        FILE: srt_h_range_accrual.h
        PURPOSE: Range accrual product implementation for GRFN
        AUTHOR: Dimitri Mayevski
        DATE: 12/06/2002
  --------------------------------------------------------------*/

#ifndef __SRT_H_RANGE_ACCRUAL_H__
#define __SRT_H_RANGE_ACCRUAL_H__

SrtProduct *
srt_f_init_range_accrual(char *name,        /* name given to the product */
                         char *yc_d,        /* domestic mkt yield curve */
                         char *yc_f,        /* index mkt yield curve */
                         char *volc_d,      /* domestic mkt vol curve */
                         char *volc_f,      /* index mkt vol curve */
                         char *refrate_d,   /* standard domestic ref
                                           rate */
                         char *refrate_f,   /* index ref rate */
                         int n_fxvol_dates, /*	\
                                             */
                         long *fxvol_dates, /*	 }	FX vol
                                           curve */
                         double *fxvol,     /*	/	(NULL if not quanto)
                                             */
                         double *qtocorrel, /* quanto correl term
                                           structure */
                         int n_periods, long *dates, /* n_periods + 1 : a start
                                                    + payment dates */
                         double *cpns,       /* coupons for each period */
                         char *basis,        /* day count basis */
                         double *upper_barr, /* upper barrier for
                                            each period */
                         double *lower_barr, /* lower barrier for
                                            each period */
                         char *recpay,       /* Payer - receiver */
                         double c_spread,    /* call spread margin */
                         int obs_freq,       /* observation frequency in
                                            days */
                         double rho_df);     /* domestic mkt - index
                                                correlation */

typedef struct _SrtRAObservStep SrtRAObservStep;
struct _SrtRAObservStep {
  long date;
  double K[4];           /* strikes corresponding to barriers */
                         /*	double		vol_atK[4];			/* volatilities at the
                          * barriers
                          */
  double betavol;        /* SABR betavol */
  double alpha;          /* SABR alpha */
  double beta;           /* SABR beta */
  double rho;            /* SABR rho */
  double correction;     /* quanto + DRS corrections */
  double spread;         /* libor spread */
  double cpn;            /* coupon * cvg */
  SrtRAObservStep *next; /* linked list */
};

typedef struct _SrtVolBumpData SrtVolBumpData;
struct _SrtVolBumpData {
  long date;   /* Exercise date */
  double bump; /* SABR SigmaBeta bump to calibrate at each exercise date */
  SrtVolBumpData *next;
};

typedef struct _SrtRangeAccrualSpec {
  int n_periods;
  long *dates; /* n_periods + 1 */

  int spotlag_f, tenor_f;
  SrtBasisCode basis_f;
  int quanto;

  SrtRAObservStep *
      *obs_steps;        /* pointers to first observation steps in periods */
  SrtVolBumpData *bumps; /* SigmaBeta bumps */

  /* Deterministic temporary values */
  long last_date;
  int period_idx;
  int ndfs_d, ndfs_f;
  long *dates_f;
  double *cvg_f;
  long *idx;

  SrtVolBumpData *bump_data;

} SrtRangeAccrualSpec;

#endif /* __SRT_H_RANGE_ACCRUAL_H__ */