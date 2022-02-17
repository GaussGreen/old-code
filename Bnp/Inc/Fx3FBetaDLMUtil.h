
#ifndef FX3FBETADLMUTIL_H
#define FX3FBETADLMUTIL_H

#define TINY 1.0e-09

#include "srt_h_all.h"

typedef struct {
  long today;

  /* domestic LGM informations */
  char *yc_dom;

  double lam_dom;
  double tau_dom;
  int nb_dom;
  double *time_dom;
  double *sig_dom;

  /* foreign LGM informations */
  char *yc_for;

  double lam_for;
  double tau_for;
  int nb_for;
  double *time_for;
  double *sig_for;

  /* Fx informations */
  double spot_fx;

  double tstar;
  double b0;
  double c0;

  int nb_fx;
  double *time_fx;
  double *sig_fx;

  /* Correlations */
  int nb_3F_corr;
  double *time_3F_corr;
  double *dom_for_3F_corr;
  double *dom_fx_3F_corr;
  double *for_fx_3F_corr;

  int nb_DLM_corr;
  double *time_DLM_corr;
  double *dom_for_DLM_corr;
  double *dom_fx_DLM_corr;
  double *for_fx_DLM_corr;
  double *dom_fx_DLM_corr_down;
  double *for_fx_DLM_corr_down;

  /* Equi Fx Beta informations */
  double equi_beta;

  /* Alpha Fudge information */
  double alpha;
  double lambda;

} FxBetaDLM_str, *FXBETADLM_STR;

typedef struct {
  int iNbPoints;
  double dLowerBound;
  double dUpperBound;
  int iMethod;
  double dAlpha;
  double dMinTime;
  double dX0;

  /* For Smile Precalc */
  int iPrecalcSmile;
  int iPriceVol;
  int iNbX;
  double dNbStdX;
  int iNbFwd;
  double dNbStdFwd;
  int iUseCxtyAdj;
  int iEquiSpacedFwd;

  /* For The Correlation Mapping */
  double dMinTimeSpaceCorrel;
  double dMaxTimeCorrel;
  int iMappingMethod;
  double dFFxCorrelTolerance;

  /* For The Alpha Fudge */
  int iMidCorrelation;

  /* To Choose the input correlation */
  int iInitWith3FCorrel;

} FxBetaDLM_OptNumerParams, *FXBETADLM_OPTNUMERPARAMS;

void FxBetaDLM_free_str(FXBETADLM_STR str);

Err FxBetaDLM_get_C0_from_beta(double TStar, double beta, double *C0);

Err FxBetaDLM_get_beta_from_C0(double TStar, double C0, double *beta);

Err FxBetaDLM_get_C0_from_EquiBeta(double TStar, double B0, double beta,
                                   double C0[3]);

Err FxBetaDLM_get_EquiBeta_from_C0(double TStar, double B0, double C0,
                                   double *beta);

Err FxBetaDLM_get_limit_beta(double TStar, double B0, double *BetaLimit);

Err FxBetaDLM_free_und(SrtUndPtr pUndDesc);

Err FxBetaDLM_fill_str(
    char *undname_fx,  /* und name */
    char *undname_dom, /* domestic underlying name */
    char *undname_for, /* foreign underlying name */

    /*	FX Underlying	*/
    double tstar, double B0, double C0, double alpha, double lambda,
    double spot_fx, int nb_fx, long *date_fx, double *sigma_fx,

    /*	Correlations	*/
    /*	3 Factor correlations */
    int compute_3F_corr, int nb_3F_corr, long *date_3F_corr,
    double *dom_for_3F_corr, double *dom_fx_3F_corr, double *for_fx_3F_corr,

    /*	Model correlations */
    int compute_DLM_corr, double max_time_correl, double min_time_space_correl,

    int nb_DLM_corr, long *date_DLM_corr, double *dom_for_DLM_corr,
    double *dom_fx_DLM_corr, double *for_fx_DLM_corr,

    /* For Alpha case */
    double *dom_fx_DLM_corr_down, double *for_fx_DLM_corr_down,

    FXBETADLM_STR str, FxBetaDLM_OptNumerParams *NumParams);

/* Initialisation of the underlying in memory */
Err SrtInitFXBetaDLMUnd(
    char *undname_fx,  /* und name */
    char *undname_dom, /* domestic underlying name */
    char *undname_for, /* foreign underlying name */

    /*	FX DLM Underlying	*/
    double tstar, double B0, double C0, double alpha, double lambda,
    double spot_fx, int nb_fx, long *date_fx, double *sigma_fx,

    /*	Correlations	*/
    /*	3 Factor correlations */
    int compute_3F_corr, int nb_3F_corr, long *date_3F_corr,
    double *dom_for_3F_corr, double *dom_fx_3F_corr, double *for_fx_3F_corr,

    /*	Model correlations */
    int compute_DLM_corr, int nb_DLM_corr, long *date_DLM_corr,
    double *dom_for_DLM_corr, double *dom_fx_DLM_corr, double *for_fx_DLM_corr,
    double *dom_fx_DLM_corr_down, double *for_fx_DLM_corr_down,
    FxBetaDLM_OptNumerParams *NumParams);

Err FxBetaDLM_get_struct_from_und(char *und, FXBETADLM_STR *str);

Err FxBetaDLM_Get_TermStruct(char *und, long *today, double *tau_dom,
                             int *nb_dom, double **time_dom, double **sig_dom,
                             double *tau_for, int *nb_for, double **time_for,
                             double **sig_for, double *tstar, int *nb_fx,
                             double **time_fx, double **sig_fx, int *nb_3F_corr,
                             double **time_3F_corr, double **dom_for_3F_corr,
                             double **dom_fx_3F_corr, double **for_fx_3F_corr,
                             int *nb_DLM_corr, double **time_DLM_corr,
                             double **dom_for_DLM_corr,
                             double **dom_fx_DLM_corr, double **for_fx_DLM_corr,
                             double **dom_fx_DLM_corr_down,
                             double **for_fx_DLM_corr_down, double *equi_beta);

typedef struct {
  long lToday;

  /* discretisation */
  int iNbPWTime;
  double *dPWTime;

  /* domestic LGM informations */
  char *cYcDom;
  double dLambdaDom;
  double dTauDom;
  double *dSigmaDom;

  /* foreign LGM informations */
  char *cYcFor;
  double dLambdaFor;
  double dTauFor;
  double *dSigmaFor;

  /* Fx informations */
  double dSpotFx;
  double dCashFx;

  double dTStar;
  double dInitTStar;
  long lTStarDate;
  double dB0;
  double dC0;
  double dAlpha;
  double dLambda;
  double *dSigmaFx3F;
  double *dSigmaFx;
  double *dSigmaFxUp;
  double *dSigmaFxDown;

  /* Correlations */
  double *dCorrDomFor;
  double *dCorrDomFx;
  double *dCorrForFx;
  double *dCorrDomFxDown;
  double *dCorrForFxDown;
  double *dCorrDomFxInput;
  double *dCorrForFxInput;

  /* Equi Fx Beta informations */
  double dEquiBeta;

  /* Pointer on the structure if model is saved */
  FxBetaDLM_str *str;

} FxBetaDLM_model, *FXBETADLM_MODEL;

void free_FxBetaDLM_model(FxBetaDLM_model *model);

Err init_FxBetaDLM_model(double MinTimeCorrel, double MaxTimeCorrel, long today,
                         char *yc_dom, double tau_dom, int nb_dom,
                         double *time_dom, double *sig_dom, char *yc_for,
                         double tau_for, int nb_for, double *time_for,
                         double *sig_for, double spot_fx, double tstar,
                         double B0, double C0, double alpha, double lambda,
                         int nb_fx, double *time_fx, double *sig_fx3F,
                         double *sig_fx, int nb_corr, double *time_corr,
                         double *dom_for_corr, double *dom_fx_corr,
                         double *for_fx_corr, FxBetaDLM_model *model);

Err FxBetaDLM_Get_Model(char *und, FxBetaDLM_model *model);

void FxBetaDLM_SetDefaultOptNumerParams(FxBetaDLM_OptNumerParams *NumParams);

typedef struct {
  int iNbPoints;
  double *X;
  double *W;

} FxBetaDLM_Hermite, FXBETADLM_HERMITE;

void free_FxBetaDLM_Hermite(FxBetaDLM_Hermite *hermite);

Err initialise_FxBetaDLM_Hermite(FxBetaDLM_OptNumerParams *NumParams,
                                 FxBetaDLM_Hermite *hermite);

typedef struct {
  int iNbPWTime;
  int iNbDone;

  double *dVarRatesDom;
  double *dVarRatesFor;
  double *dVarRatesFor3D;
  double *dCovarRates;
  double *dCovarRates3D;

  double *dVarFFx;
  double *dCovarFFxDom;
  double *dCovarFFxFor;
  double *dCovarFFxFor3D;

  double *dVarX;
  double *dExpectFor;
  double *dIntegral;
  double *dIntegral2;
  double *dCovarRatesAdjust;

  double *dVarRatesFor_down;
  double *dVarRatesFor3D_down;
  double *dCovarRates_down;
  double *dCovarRates3D_down;

  double *dVarFFx_down;
  double *dCovarFFxDom_down;
  double *dCovarFFxFor_down;
  double *dCovarFFxFor3D_down;

  double *dVarX_down;
  double *dExpectFor_down;
  double *dIntegral_down;
  double *dIntegral2_down;
  double *dCovarRatesAdjust_down;

  double *dVarRatesFor_mid;
  double *dVarRatesFor3D_mid;
  double *dCovarRates_mid;
  double *dCovarRates3D_mid;

  double *dVarFFx_mid;
  double *dCovarFFxDom_mid;
  double *dCovarFFxFor_mid;
  double *dCovarFFxFor3D_mid;

  double *dVarX_mid;
  double *dExpectFor_mid;
  double *dIntegral_mid;
  double *dIntegral2_mid;
  double *dCovarRatesAdjust_mid;
} FxBetaDLM_ModelPrecalc, FXBETADLM_MODELPRECALC;

Err FxBetaDLM_Allocation_Precalculations(int iNbPWTime,
                                         FxBetaDLM_ModelPrecalc *Precalc);

void FxBetaDLM_Free_Precalculations(FxBetaDLM_ModelPrecalc *Precalc);

Err quadratic_bs(double constant_coef, double linear_coef,
                 double quadratic_coef, double mean, double variance,
                 double strike, SrtCallPutType callput, double *Premium);

typedef struct {
  double min_time;
  int numeraire; /*0: QBeta      , 1:QTStar */
  int do_pecs;
  int do_discount;

} FxBetaDLM_GRFNNumerParams, *FXBETADLM_GRFNNUMERPARAMS;

void FxBetaDLM_Init_GRFNNumerParams(FxBetaDLM_GRFNNumerParams *NumParams);

Err FxBetaDLM_correl_mapping(double volX, double B0, double C0, double Tstar,
                             double SigmaDom, double LambdaDom, double SigmaFor,
                             double LambdaFor, double RhoDF, double RhoDS,
                             double RhoFS, double previous_t, double next_t,
                             double VarFFx, double *RhoDX, double *RhoFX);

Err FxBetaDLM_correl_mapping_forward(
    double CorrelTolerance, double t, double previous_t, double B0, double C0,
    double varFFx, double TStar, double SigmaDom, double LambdaDom,
    double SigmaFor, double LambdaFor, double Sigma3F, double SigmaDLM,
    double RhoDF, double RhoDS, double RhoFS, double *RhoDX, double *RhoFX);

Err Get_Correl_FromModel(int index, double VarFFx, double VolX, double RhoXD,
                         double RhoXF, FxBetaDLM_model *model, double *RhoSD,
                         double *RhoSF);

Err FxBetaDLM_GetFirstGuessFromB0(double B0, double Tstar, double *exercise_opt,
                                  double *maturity_opt, long nbrOpt,
                                  double *maturity_rates_corr, long nbrMat,
                                  double *sig_curve_dom, double lda_dom,
                                  double *sig_curve_for, double lda_for,
                                  double *correl_dom_for, double *correl_dom_fx,
                                  double *correl_for_fx, double *cal_vol_3f,
                                  double *fx_vol_curve);

Err FxBetaDLM_GetModelFromStr(FXBETADLM_STR str, FxBetaDLM_model *model);

Err FxBetaDLM_get3FCorrelFromDLM(FXBETADLM_STR str);

Err FxBetaDLM_getDLMCorrelFrom3F(FXBETADLM_STR str,
                                 FxBetaDLM_OptNumerParams *NumParams);

Err FxBetaDLM_CheckOrMakePosMatrix(double RhoDF, double *RhoDX, double *RhoFX);

Err copy_FxBetaDLM_model(FxBetaDLM_model *model_source,
                         FxBetaDLM_model *model_dest);

#endif