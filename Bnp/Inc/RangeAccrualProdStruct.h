//----------------------------------------------------------------------
//---------- From D.M.'s srt_h_rangeaccrual.h --------------------------
//----------------------------------------------------------------------
#ifndef __RANGE_ACCRUAL_PROD_STRUCT_H
#define __RANGE_ACCRUAL_PROD_STRUCT_H

#define RA_NCPN 512

/*	Structures and functions for range accrual leg */
/*	-------------------------------------------------- */

typedef struct _RangeAccrualObservStep RangeAccrualObservStep;
struct _RangeAccrualObservStep {
  long date;
  double K[4];       /* strikes corresponding to barriers */
  double lognvol;    /* SABR lognvol */
  double normvol;    /* SABR normvol */
  double betavol;    /* SABR betavol */
  double alpha;      /* SABR alpha */
  double beta;       /* SABR beta */
  double rho;        /* SABR rho */
  double correction; /* quanto + DRS corrections */
  double spread;     /* libor spread */

  double cpn; /* coupon * cvg */

  double cvg;               /* cvg for floating coupon */
  double cpn_float_gearing; /* float coupon gearing */
  double cpn_float_spread;  /* float coupon basis spread */
  double cpn_float_past_fixing;
  double cpn_float_correction;
  double cpn_float_lower_barrier_adj;
  double cpn_float_upper_barrier_adj;

  RangeAccrualObservStep *next; /* linked list */
};

typedef struct _RAVolBumpData RAVolBumpData;
struct _RAVolBumpData {
  long date;   /* Exercise date */
  double bump; /* SABR SigmaBeta bump to calibrate at each exercise date */
  RAVolBumpData *next;
};

typedef struct _RangeAccrualStruct {
  double notional;
  int n_periods;
  long *dates; /* n_periods + 1 */

  int spotlag_f, tenor_f;
  SrtBasisCode basis_f;
  char tenor_f_char[256];
  int quanto;

  int typeVol; /* 0 : Beta Vol      , 1 : Normal Vol      , 2 : Lognormal Vol*/

  RangeAccrualObservStep *
      *obs_steps;       /* pointers to first observation steps in periods */
  RAVolBumpData *bumps; /* SigmaBeta bumps */

  /* Deterministic temporary values */
  long last_date;
  int period_idx;
  int ndfs_d, ndfs_f, ndfs_cpn_float;
  long *dates_f;
  double *cvg_f;
  long *idx;

  /* If RA floating coupon */
  int tenor_float_cpn;
  int float_cpn_is_dom_for;
  int cpn_type;
  SrtBasisCode basis_float_cpn;
  char tenor_float_cpn_char[256];
  long *float_cpn_dates;
  long nb_float_cpn_dates;
  double *cvg_float_cpn;

  RAVolBumpData *bump_data;

} RangeAccrualStruct;

Err ra_init_struct(
    long today, char *yc_d, char *yc_f, char *volc_d, char *volc_f,
    char *refrate_d, char *refrate_f, double notional, int ra_cpn_type,
    int n_fxvol_dates, long *fxvol_dates, double *fxvol, double *qtocorrel,
    int typeVol, int n_periods, long *dates, /* n_periods + 1 */
    double *cpns, char *basis, int *ra_nfixings, long **ra_fixingdates,
    double **ra_fixings, /*	Past coupon fixing if relevant */

    // RA floating coupon
    int ra_float_refrate_is_dom_for, char *ra_float_refrate,
    long *ra_float_startl, double *ra_float_past_fixings,
    double *ra_float_gearings,

    double *upper_barr, double *lower_barr, char *recpay, char *ra_refrate,
    double c_spread, int obs_freq, double rho_df,

    // Params for floating coupon
    double correl_start, double correl_end, int float_adj_strike,

    int eod_flag, /*	0: I      , 1: E */
    RangeAccrualStruct *RA);

Err ra_free(RangeAccrualStruct *RA);

Err RA_RequestDfDates(RangeAccrualStruct *RA, long date, int *pn_unds,
                      int **pn_dfs, long ***pdates);

Err RA_Payoff(RangeAccrualStruct *RA, long today, long date, double **dfs,
              double *payoff);

Err RA_FwdPV(char *yc_d, char *yc_f, RangeAccrualStruct *RA, long today,
             int nbDates, long *Dates, double *fwdPV);

typedef struct {
  long fixing_date;
  long start_date;
  long end_date;
  double K[4];       /* strikes corresponding to barriers */
  double betavol;    /* SABR betavol */
  double alpha;      /* SABR alpha */
  double beta;       /* SABR beta */
  double rho;        /* SABR rho */
  double correction; /* quanto + DRS corrections */
  double spread;     /* libor spread */
  double cpn;        /* coupon * cvg */

} ra_obs, *RA_OBS;

typedef struct {
  long start_date;   /*	Coupon start date */
  double start_time; /*	Coupon start time */
  long pay_date;     /*	Coupon pay date */
  double pay_time;   /*	Coupon pay time */
  double cvg;        /*	Cvg */
  double fixcpn;     /*	Fix Coupon */

  int num_fixings;   /*	number of fixings */
  long *fixingdates; /*	fixing dates */

  int num_obs;
  ra_obs *obs;

} ra_cpn, *RA_CPN;

typedef struct {
  double notional;
  int spotlag_f;
  int tenor_f;
  SrtBasisCode basis_f;
  char tenor_f_char[256];
  int num_cpn;
  /*	0..num_cpn-1 */
  ra_cpn *cpn;

} ra_leg, *RA_LEG;

#endif