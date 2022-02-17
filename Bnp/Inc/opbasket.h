#include "utallhdr.h"
#include <NUM_H_ALLHDR.H>

#ifndef OPBASKETCAP_H
#define OPBASKETCAP_H

Err srt_f_BasketCapPrice(int iNbRates, double dStartDate, double dStrike,
                         double *dForwards, double *dVols, double **dCorrels,
                         SrtGreekType SrtGreek, SrtCallPutType SrtCallPut,
                         SrtDiffusionType SrtVolType, double *dAnswer);

/*	Price X * max ( 0   , sum [ai * min (Yi   , Ki)] - Kx ) */
double
inflation_multiput(double T, int nidx, double *ai, double *fyi, double *ki,
                   double *si, double fx, double kx, double sx,
                   double **rho, /* index 0..nidx-1: xi  , index nidx: x */
                   int npth);

/*	Price X * max ( 0   , sum [ai * min (Yi   , Ki)] - Kx ) in a shifted-log
 * model */
double inflation_multiput_sl(
    double T, int nidx, double *ai, double *fyi, double *ki, double *shifti,
    double *voli, double fx, double kx, double sx,
    double **rho,              /* index 0..nidx-1: xi  , index nidx: x */
    int npth, int atmbumptype, /*	0: no bump  , 1: add  , 2: mult */
    int betabumptype,          /*	0: no bump  , 1: add  , 2: mult */
    double *atmbump, double *betabump);

/*	Find shifted-log parameters by calibration to 2 strikes */
/*	The model is:
. F = L - shift or L = F + shift where L is a lognormal martingale with a
certain vol . If the vol of L is negative in the case of a sub-normal skew  ,
then F = -L - shift or L = -F - shift where L is a lognormal martingale with vol
abs (vol) and led by a reversed Brownian Motion from F . The implied beta is
also calculated */
Err find_shifted_log_params(double T,   /*	Maturity in years */
                            double fwd, /*	Forward */
                            double k1,  /*	1st strike */
                            double v1,  /*	BS implied vol for 1st strike */
                            double k2,  /*	2nd strike */
                            double v2,  /*	BS implied vol for 2nd strike */
                            /*	Results */
                            double *shift, double *vol, double *beta);

#define MAX_EQD 50
typedef struct {
  double mat;
  double df;

  int num_eqd;
  double ref[MAX_EQD];
  double spot[MAX_EQD];
  double fwd[MAX_EQD];
  double strike1[MAX_EQD];
  double strike2[MAX_EQD];
  double vol1[MAX_EQD];
  double vol2[MAX_EQD];
  double shift[MAX_EQD];
  double used_shift[MAX_EQD];
  double vol[MAX_EQD];
  double beta[MAX_EQD];
  double qto[MAX_EQD];
  double used_fwd[MAX_EQD];

  double cpi_fwd;
  double used_cpi_fwd;
  double cpi_vol;

  double rho[MAX_EQD][MAX_EQD];
} i_stellar_info;

char *i_stellar_cpn(
    /*	The product */
    int prod_type,        /*	0: Stellar  , 1: I-Stellar */
    double fix_mat_years, /*	(Fixing date - max (Ref date  , Today)) / 365 */
    double vol_mat_years, /*	(Fixing date - Today) / 365 */
    int nidx,             /*	Number of eqd indices in play */
    double *
        idx_weights, /*	Weights of the indices in the payoff e.g. 0.45 */
    double
        *idx_ref,      /*	The base values for indices
                                                       historical fixing if Ref Date
                          <= Today      or interpolated from forward curve if Ref date >
                          Today */
    double *idx_str,   /*	Stellar strikes in % of ref values  , e.g. 1.065 */
    int *idx_qto,      /*	Wether eqd index is quantoed */
    double cpi_ref,    /*	The base value for CPI
                                                       historical fixing if Ref
                          Date <= Today    or interpolated from forward curve if Ref
                          date > Today */
    double global_str, /*	Global strike for the payoff */
    /*	The IR market */
    double pay_df, /*	DF to payment date times coverage */
    /*	The EQD market */
    double *idx_spots, /*	Spots of indices */
    double *idx_fwd,   /*	Forwards of indices for fixing date interpolated from
                                                       forward curves */
    char **idx_vol_name,  /*	Names of the vol curve of the indices */
    double (*get_eqd_vol)(/*	GetVol function for eqd indices */
                          double mat_years,   /*	maturity in years */
                          double strike_spot, /*	strike in % of spot  ,
                                                 e.g. 1.065 */
                          char *name),        /*	name of the vol curve */
    /*	Smile parameters
            0: ATMS
            1: ATMF
            2: ISTR
            3: Custom */
    int strike_1_type, int strike_2_type, double *strike_1_custom,
    double *strike_2_custom,
    /*	The CPI market */
    double cpi_fwd, /*	Forward CPI for fixing date interpolated from forward
                       curve */
    double
        cpi_vol, /*	CPI vol for fixing date interpolated from vol curve
                                                 WE NEED CUMULATIVE CPI VOL */
    /*	The Fx market */
    double *fx_vol,    /*	For each eqd index  , the volatility of the
                          corresponding Fx    interpolated from the implied BS
                          volatility curve of Fx    for the fixing date */
    double *fx_correl, /*	Correl between the eqd index and the fx
                                                       EXPRESSED IN EQD CCY /
                          TRADE(cpi) CCY interpolated from correlation curve by
                          tenor (not by fixed maturity) */
    /*	The basket correl as interpolated for the fixing date from the relevant
            term structures BY TENOR */
    double **rho, /*	0..nidx-1: eqd indices  , nidx: CPI */
    /*	Parameters */
    int integ_points, /*	10 to 15  , default 12 */
    /*	MAD: 0 everywhere */
    int atmf_bump_type, /*	0: none  , 1: add  , 2: mult */
    int beta_bump_type, double *atmf_bump, double *beta_bump,
    /*	Result */
    double *pv, int store_info, i_stellar_info *info);

/*	CMT Stellar */
/*	XXXXXXXXXXX	*/

/*	Price max ( ax * X   , sum [ai * min (Yi   , Ki)] - Kx ) */
double cmt_stellar(double T, int nidx, double *ai, double *fyi, double *ki,
                   double *si, double ax, double fx, double kx, double sx,
                   double **rho, /* index 0..nidx-1: xi  , index nidx: x */
                   int npth);

/*	Price max ( ax * X   , sum [ai * min (Yi   , Ki)] - Kx ) in a
 * shifted-log model */
double cmt_stellar_sl(
    double T, int nidx, double *ai, double *fyi, double *ki, double *shifti,
    double *voli, double ax, double fx, double kx, double sx,
    double **rho,              /*	index 0..nidx-1: xi  , index nidx: x */
    int npth, int atmbumptype, /*	0: no bump  , 1: add  , 2: mult */
    int betabumptype,          /*	0: no bump  , 1: add  , 2: mult */
    double *atmbump, double *betabump);

typedef struct {
  double mat;
  double df;

  int num_eqd;
  double ref[MAX_EQD];
  double spot[MAX_EQD];
  double fwd[MAX_EQD];
  double strike1[MAX_EQD];
  double strike2[MAX_EQD];
  double vol1[MAX_EQD];
  double vol2[MAX_EQD];
  double shift[MAX_EQD];
  double used_shift[MAX_EQD];
  double vol[MAX_EQD];
  double beta[MAX_EQD];
  double qto[MAX_EQD];
  double used_fwd[MAX_EQD];

  double cmt_fwd;
  double cmt_vol;

  double rho[MAX_EQD][MAX_EQD];
} cmt_stellar_info;

char *cmt_stellar_cpn(
    /*	The product */
    int prod_type,        /*	0: Stellar  , 1: CMT-Stellar */
    double fix_mat_years, /*	(Fixing date - max (Ref date  , Today)) / 365 */
    double vol_mat_years, /*	(Fixing date - Today) / 365 */
    int nidx,             /*	Number of eqd indices in play */
    double *
        idx_weights, /*	Weights of the indices in the payoff e.g. 0.45 */
    double
        *idx_ref,      /*	The base values for indices
                                                       historical fixing if Ref Date
                          <= Today      or interpolated from forward curve if Ref date >
                          Today */
    double *idx_str,   /*	Stellar strikes in % of ref values  , e.g. 1.065 */
    int *idx_qto,      /*	Wether eqd index is quantoed */
    double global_str, /*	Global strike for the payoff */
    /*	The IR market */
    double pay_df, /*	DF to payment date times coverage */
    /*	The EQD market */
    double *idx_spots, /*	Spots of indices */
    double *idx_fwd,   /*	Forwards of indices for fixing date interpolated from
                                                       forward curves */
    char **idx_vol_name,  /*	Names of the vol curve of the indices */
    double (*get_eqd_vol)(/*	GetVol function for eqd indices */
                          double mat_years,   /*	maturity in years */
                          double strike_spot, /*	strike in % of spot  ,
                                                 e.g. 1.065 */
                          char *name),        /*	name of the vol curve */
    /*	Smile parameters
            0: ATMS
            1: ATMF
            2: ISTR
            3: Custom */
    int strike_1_type, int strike_2_type, double *strike_1_custom,
    double *strike_2_custom,
    /*	The CMT market */
    double cmt_weight, double cmt_fwd, /*	Forward CMT for fixing date
                                          interpolated from forward curve */
    double
        cmt_vol, /*	CMT vol for fixing date interpolated from vol curve
                                                 WE NEED CUMULATIVE CMT VOL */
    /*	The Fx market */
    double *fx_vol,    /*	For each eqd index  , the volatility of the
                          corresponding Fx    interpolated from the implied BS
                          volatility curve of Fx    for the fixing date */
    double *fx_correl, /*	Correl between the eqd index and the fx
                                                       EXPRESSED IN EQD CCY /
                          TRADE(cpi) CCY interpolated from correlation curve by
                          tenor (not by fixed maturity) */
    /*	The basket correl as interpolated for the fixing date from the relevant
            term structures BY TENOR */
    double **rho, /*	0..nidx-1: eqd indices  , nidx: CMT */
    /*	Parameters */
    int integ_points, /*	10 to 15  , default 12 */
    /*	MAD: 0 everywhere */
    int atmf_bump_type, /*	0: none  , 1: add  , 2: mult */
    int beta_bump_type, double *atmf_bump, double *beta_bump,
    /*	Result */
    double *pv, int store_info, cmt_stellar_info *info);

#endif