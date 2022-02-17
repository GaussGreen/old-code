#ifndef __CPDVOL_H
#define __CPDVOL_H

/*===================================================

        SMILE MANAGEMENT FOR FX HYBRIDS

        supports:	SABR
                                BMM

===================================================*/

//=========================
// Structure for the vol market
//=========================
typedef struct {
  int smile_spec_type;

  int num_vols;

  double *forward;
  double *sigma;
  double *sigmabeta;
  double *alpha;
  double *beta;
  double *rho;
  double *pi;
  double *times;
  long *dates;
} smile_vol_market, *SMILE_VOL_MARKET;
//=========================
// Structure for the smile at one date
//=========================
typedef struct {
  int smile_spec_type;

  double forward;
  double sigma;
  double sigmabeta;
  double alpha;
  double beta;
  double rho;
  double pi;
  double times;
  long dates;

  // Specific to blksh beta
  int Use_BetaQuick_in_MC; // if 1  , uses SABR approx for blkschbeta in the
                           // forward pricing.
} smile_parameters, *SMILE_PARAMETERS;

Err cpd_alloc_smile_vol_market(int num_vols, SMILE_VOL_MARKET smile_mkt);

Err cpd_fill_smile_vol_market(
    long today,
    int smile_spec_type, // 0 SABR with ATMLOG  , 1 SABR with ATMBETA  , 2 BMM
                         // ...
    int num_vols, double *forward, double *sigma, double *alpha, double *beta,
    double *rho, double *pi, double *times, long *dates,
    SMILE_VOL_MARKET smile_mkt);

Err cpd_check_smile_vol_market(SMILE_VOL_MARKET smile_mkt);

Err cpd_free_smile_vol_market(SMILE_VOL_MARKET smile_mkt);

Err cpd_vol_get_vol(double forward, double fix_time, double strike,
                    SMILE_PARAMETERS smile_params, SrtDiffusionType output_vol,
                    double *smile_std);

Err cpd_vol_get_price(int type, double Forward, double Maturity, double Strike,
                      SMILE_PARAMETERS smile_params, double *Value);

Err cpd_vol_get_smile_params(int IsTime, // 0 will use date  , 1 will use time
                             long date, double time, SMILE_VOL_MARKET smile_mkt,
                             SMILE_PARAMETERS smile_params);

Err cpd_vol_get_linterpVol(int IsTime, // 0 will use date  , 1 will use time
                           long date, double time, SMILE_VOL_MARKET smile_mkt,
                           SMILE_PARAMETERS smile_params);

typedef struct _otcpd_precalc otcpd_precalc, *OTCPDPRECALC;

Err cpd_vol_get_smile_params_from_otc(int smile_spec_type, int row_idx,
                                      int col_idx, OTCPDPRECALC otc_precalc,
                                      SMILE_PARAMETERS smile_params);

#endif
