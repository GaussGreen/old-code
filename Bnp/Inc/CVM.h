/*	-------------------------------------------- */
/*	---- Continuous Volatility Model (CVM) ----- */
/*	-------------------------------------------- */

#ifndef __CVM_H
#define __CVM_H

typedef char char256[256];

typedef struct {
  // Market ID
  char ycid[256];
  char vcid[256];

  // Today      , SpotDate and SpotLag
  long today;
  long spot;
  long spotLag;

  // Swaption Details
  SrtCompounding srtFreq;
  SrtBasisCode srtBasis;
  char refRate[256];
  SrtCompounding floatsrtFreq;
  SrtBasisCode floatsrtBasis;

  // Shedules
  // Fix Schedule
  long iNFixPayDates;
  long iNFixDates;
  long *lFixPayDates;
  long *lFixStartDates;
  long *lFixEndDates;
  double *dFixCoverages;

  // Float Schedule
  long iNFloatPayDates;
  long iNFloatDates;
  long *lFloatFixingDates;
  long *lFloatPayDates;
  long *lFloatStartDates;
  long *lFloatEndDates;
  double *dFloatCoverages;
  double *dFloatSpreads;

  // Vol TS
  double
      shift; // -1 : Normal      , 0 : Lognormal      , others Shifted Lognormal
  int converging_sliding;  // 0 : converging      , 1: sliding
  int interpollation_mode; // 0 : piecewise constant      , 1 : continuous
                           // piecewise linear
  int nbOfDates;
  double *vol_dates; // vol_dates[DateIndex]
  double *vol_times; // vol_times[DateIndex]

  // Pillars
  int nbOfPillars;
  char256 *pillars_tenors;    // pillars_tenors[PillarIndex]
  double *pillars_maturities; // pillars_maturities[PillarIndex]
  double *pillars_times;      // pillars_maturities[PillarIndex]

  // Pillars Correlations and Vols
  double ***
      pillars_correls; // pillars_correlations[DateIndex][PillarIndex1][PillarIndex2]
  double ***pillars_chol; // pillars_chol[DateIndex][PillarIndex1][PillarIndex2]
  double **
      *pillars_angles; // pillars_angles[DateIndex][PillarIndex1][PillarIndex2]
  double ***pillars_vols;     // pillars_vols[DateIndex][PillarIndex][Component]
  double **pillars_vols_norm; // pillars_vols_norm[DateIndex][PillarIndex]

  // Parametric Sensitivities
  int sensi_type;          // 0 : use spline      , 1 : custom      , 2 : linear
  double **lambdas;        // lambdas[DateIndex][PillarIndex]
  double **sensitivities;  // pillars_sensitivities[PillarIndex][MaturityIndex]
  double **sensitivities2; // pillars_sensitivities[PillarIndex][MaturityIndex]

  // Custom Sensitivities
  int nbOfSensiTenors;
  char256 *sensi_tenors; // sensi_tenors[SensiTenorIndex]
  double *sensi_times;   // sensi_times[SensiTenorIndex][PillarIndex]

  // Numerical Parameter
  int nQuadLegendre;

  // Calibration Inputs
  int nbOfOptionTenors;
  char256 *option_tenors; // option_tenors[OptionIndex]
  int nbOfUnderlyingTenors;
  char256 *underlying_tenors; // underlying_tenors[UnderlyingIndex]
  double **market_vols;       // market_vols[OptionIndex][UnderlyingIndex]
  double **model_vols;        // model_vols[OptionIndex][UnderlyingIndex]
  int allow_negativeVols;     // 0 : Vols cannot be negative

} cvm_str, *CVM_STR;

Err cvm_fill_struct(
    // Market ID
    char *ycid, char *vcid,
    // Swaption Details
    SrtCompounding srtFreq, SrtBasisCode srtBasis, char *refRate,

    // TS Dates
    double shift, int converging_sliding, int interpollation_mode,
    int nbOfDates, double *vol_dates,

    // Pillars
    int nbOfPillars, char **pillars_tenors,

    // Pillars Correlations and Vols
    double ***pillars_correls, double **pillars_vols_norm,

    // Sensitivities
    int sensitype,

    // Parametric Sensitivities
    double **lambdas,

    // Custom Sensitivities
    int nbOfSensiTenors, char **sensi_tenors, double **sensitivities,

    // Numerical Parameter
    int nQuadLegendre,

    cvm_str *cvm);

Err cvm_set_pillars_vols(int DateIndex, int nPillars, double *pillarsVols,
                         cvm_str *cvm);

Err cvm_compute_fwdvol_at_voldate(int dateIndex, double maturity, int nPillars,
                                  double *fwdvol, cvm_str *cvm);

Err cvm_fwdvol(double date, double maturity, cvm_str *cvm);

Err cvm_convert_correls_to_angles(int nbOfPillars, int nbOfDates,
                                  double ***pillars_correls,
                                  double ***pillars_angles, cvm_str *cvm);

Err cvm_vol(double start, double end, double *vol, cvm_str *cvm);

Err cvm_partial_forward_vol(double date, double expiry, double start,
                            double end, double *vol, cvm_str *cvm);

Err cvm_partial_forward_cov(double date, double expiry, double start1,
                            double end1, double start2, double end2,
                            double *vol1, double *vol2, double *cov,
                            cvm_str *cvm);

Err cvm_vol_converging(double start, double end, double *vol, cvm_str *cvm);

Err cvm_vol_sliding(double start, double end, double *vol, cvm_str *cvm);

Err cvm_compute_sensitivities(double date, double maturity, int nPillars,
                              double *sensi, cvm_str *cvm);

Err sortspline(int nbasis, double *xbasis, double *ybasis, int n, double *x,
               double *y);

Err cvm_free_struct(cvm_str *cvm);

Err cvm_free_und_struct(SrtUndPtr pUndDesc);

Err cvm_get_struct_from_und(char *und, cvm_str **cvm);

Err cvm_copy_struct(cvm_str *cvm_target, cvm_str *cvm_source);

char *SrtInitCVMUnd(char *undName,
                    // Market ID
                    char *ycname, char *vcname,
                    // Swaption Details
                    SrtCompounding srtFreq, SrtBasisCode srtBasis,
                    char *refRate,

                    // TS Dates
                    double shift, int converging_sliding,
                    int interpollation_mode, int nbOfDates, double *vol_dates,

                    // Pillars
                    int nbOfPillars, char **pillars_tenors,

                    // Pillars Correlations and Vols
                    double ***pillars_correls, double **pillars_vols,

                    // Sensitivities
                    int sensitype, double **lambdas,

                    int nbOfSensiTenors, char **sensi_tenors,
                    double **sensitivities,

                    // Numerical Parameter
                    int nQuadLegendre);

Err cvm_bump_pillars_correls(int dateIndex, int nPillarsMinusOne,
                             double *correls, cvm_str *cvm);

Err cvm_set_pillars_correls(int dateIndex, int nPillars, double **correls,
                            cvm_str *cvm);
#endif