#ifndef __COBRA_WRAP__

#ifndef LGMSVCALIBAPPROXH
#define LGMSVCALIBAPPROXH

#include "DiagCalibDLM.h"
#include "DiagcalibGen.h"
#include "LGMSVClosedformApprox.h"
#include "cpdcalib.h"
#include "lgmsvpde.h"
#include "srt_h_all.h"

Err LGMSVCalibApprox(
    /* Instrument informations	*/
    int    nex, /*	Total number of exercise dates */
    double ex_time[],
    double ex_lfwd[],
    double ex_llvl[],
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double shift[],
    double coef_vol[],
    double coef_meanrev[],

    /* First guess */
    double sig[],

    /* Model informations */
    double  dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,
    long*   lSigIndex,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Newton parameters */
    double Precision,
    int    NbIterMax,

    /* Output */
    double* dSigmaTS);

Err cpd_calibSV_approx(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    char* instr_freq, /*	Frequency and basis of instruments */
    char* instr_basis,
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates,
                                                      all supposed to be on or after today */
    long*  ex_date,   /*	Supposed to be sorted */
    char** end_tenor, /*	Tenors of the underlying instruments
                                                              or "DIAG" */
    long    end_date, /*	End date for diagonal */
    double* strike,   /*	Strikes
                                              0: ATM */
    /*	Model Parameters */
    double  dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Newton parameters */
    double Precision,
    int    NbIterMax,
    /*	Output */
    int*      numres, /*	Answer */
    double**  restime,
    double*** result,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data);

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
                with lambda calibration */
Err cpd_calib_diagonal_LGMSV_dlm(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /* Get Cash Vol ref */
    char* vol_ref_rate_name,

    /* Long Instruments */
    char* instr_long_freq, /*	Frequency and basis of instruments */
    char* instr_long_basis,
    char* long_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datesl, /*	Long Exercise dates */
    long*   ex_datel_,     /*	Supposed to be sorted */
    int*    cal_datel,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenorl,    /*	Tenors of the underlying instruments */
    long    end_datel,     /*	End date for diagonal */
    double* strikel_,      /*	Strikes */
    // double			*strikeS1l_,						/*	First strike for smile calibration
    // */ double			*strikeS2l_,						/*
    // Second strike for smile calibration */

    CPD_DIAG_CALIB_PARAM paraml,

    /* Short Instruments */
    char* instr_short_freq, /*	Frequency and basis of instruments */
    char* instr_short_basis,
    char* short_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datess, /*	Short Exercise dates */
    long*   ex_dates_,     /*	Supposed to be sorted */
    int*    cal_dates,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenors,    /*	Tenors of the underlying instruments */
    long    end_dates,     /*	End date for diagonal */
    double* strikes_,      /*	Strikes */

    CPD_DIAG_CALIB_PARAM params,

    /*	Model */
    int fix_lambda,
    // int				calib_smile,						/* 0: no calib, 1: total calib
    // */
    int     one2F,
    double* dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,
    double  lgm_alpha,
    double  lgm_gamma,
    double  lgm_rho,
    double* dRho2TS,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Newton parameters */
    double Precision,
    int    NbIterMax,

    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    double** alpha,
    double** lambdaeps,
    double** rho,
    double** rho2,

    /*	Parameters */
    DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

Err LGMSV_Get_Rho_From_RhoTarget(
    double  rho_target,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* rho1,
    double* rho2);

Err LGMSV_Get_RhoTarget_From_Rho(
    double  rho1,
    double  rho2,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* rho_target);

Err LGMSV_Get_RhoTargetSensitivity_From_Rho(
    double  rho1,
    double  rho2,
    double  rho_mat1,
    double  rho_mat2,
    double  LGM_alpha,
    double  LGM_gamma,
    double  LGM_rho,
    double* sensi1,
    double* sensi2);

typedef struct
{
    int fix_lambda;
    int use_lgm_lambda;

    int calib_flat_smile;

    int    use_sabr_calib;
    int    use_sabr_levenberg;
    double sabr_calib_min_time;
    double sabr_calib_max_lambda;
    double sabr_calib_default_lambda;
    double max_calib_alpha;

    int calib_rr_bt;
    int calib_alpha;
    int calib_lameps;
    int calib_rho;
    int calib_rho2;

    int calib_smile_on_prim;
    int novolcalib_on_smile_calib;
    int pricerho_on_alpha;

    int    onefac_rho;
    double rho_mat1;
    double rho_mat2;

    int use_lgm_first_guess;

    int recalib_at_end_lambda;
    int recalib_at_end_alpha;
    int recalib_at_end_lameps;
    int recalib_at_end_rho;

    double alpha_sv_shift;
    double lam_sv_shift;
    double rho_sv_shift;
    double rho2_sv_shift;

} LGMSV_CalibParams, *LGMSV_CALIBPARAMS;

void LGMSV_SetDefault_CalibParams(LGMSV_CALIBPARAMS);

void LGMSV_Copy_CalibParams(LGMSV_CALIBPARAMS source, LGMSV_CALIBPARAMS target);

typedef struct
{
    int     nb_lgm_vol;
    double* cum_var_sv;
    double* cum_var_lgm;
    double* lgm_vol;
    double* exp_fact;

} LGMSV_EquiLGM, *LGMSV_EQUILGM;

Err Initialise_EquiLGM(LGMSV_MODEL model, LGMSV_EQUILGM equi_lgm);

Err Update_EquiLGM(LGMSV_MODEL model, LGMSV_EQUILGM equi_lgm);

void Free_EquiLGM(LGMSV_EQUILGM equi_lgm);

typedef struct
{
    LGMSV_PRICINGCONST PricingConst;

    int                 iNbLongInst;
    void**              AllLongInst;
    long*               lSigLongIndex;
    LGMSV_HESTONINST    HestonLongInst;
    LGMSV_NUMERINST     NumerLongInst;
    CALIBCPNSCHEDULEDLM CalibCpnLongSchedule;
    CALIBEXESCHEDULEDLM CalibExeLongSchedule;

    int                 iNbShortInst;
    void**              AllShortInst;
    long*               lSigShortIndex;
    LGMSV_HESTONINST    HestonShortInst;
    LGMSV_NUMERINST     NumerShortInst;
    CALIBCPNSCHEDULEDLM CalibCpnShortSchedule;
    CALIBEXESCHEDULEDLM CalibExeShortSchedule;

    long* lSmileLongIndex;
    long* lSmileShortIndex;

    INSTPARAMS* LongInstParams;
    INSTPARAMS* ShortInstParams;

    CALIBGEN_PARAMS CalibParams;
    CALIBGEN_PARAMS CalibParamsLambda;
    CALIBGEN_PARAMS CalibParamsAlpha;
    CALIBGEN_PARAMS CalibParamsLamEps;
    CALIBGEN_PARAMS CalibParamsRho;
    CALIBGEN_PARAMS CalibParamsRho2;

    CalibFunctions CalibFunctionsForVol;
    CalibFunctions CalibFunctionsForLambda;
    CalibFunctions CalibFunctionsForAlpha;
    CalibFunctions CalibFunctionsForLamEps;
    CalibFunctions CalibFunctionsForRho;
    CalibFunctions CalibFunctionsForRho2;

    int    iIndexStrike;
    int    SkipLast;
    double MinFact;
    double MaxFact;

    CPD_DIAG_CALIB_PARAM long_param;
    CPD_DIAG_CALIB_PARAM short_param;

    LGMSV_CalibParams* calib_params;

    double sens_lambda;

    int index_vol;

    LGMSV_EQUILGM EquiLGM;

} LGMSV_AllParams, *LGMSV_ALLPARAMS;

Err Initialise_AllParams(
    LGMSV_MODEL model,

    long*                  lSigLongIndex,
    CALIBCPNSCHEDULEDLM    CalibLongCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibLongInst,
    CPD_DIAG_CALIB_PARAM   paraml,

    long*                  lSigShortIndex,
    CALIBCPNSCHEDULEDLM    CalibShortCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibShortInst,
    CPD_DIAG_CALIB_PARAM   params,

    LGMSV_CalibParams* CalibParams,

    LGMSV_NUMERPARAMS NumerParams,

    int NbLGMVol,

    LGMSV_ALLPARAMS AllParams);

void Free_AllParams(LGMSV_ALLPARAMS AllParams);

Err cpd_calib_diagonal_LGMSV_new_dlm(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /* Get Cash Vol ref */
    char* vol_ref_rate_name,

    /* Long Instruments */
    char* instr_long_freq, /*	Frequency and basis of instruments */
    char* instr_long_basis,
    char* long_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datesl, /*	Long Exercise dates */
    long*   ex_datel_,     /*	Supposed to be sorted */
    int*    cal_datel,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenorl,    /*	Tenors of the underlying instruments */
    long    end_datel,     /*	End date for diagonal */
    double* strikel_,      /*	Strikes */
    double* strikeS1l_,
    double* strikeS2l_,

    CPD_DIAG_CALIB_PARAM paraml,

    /* Short Instruments */
    char* instr_short_freq, /*	Frequency and basis of instruments */
    char* instr_short_basis,
    char* short_ref_rate_name, /*	Name of the reference rate */

    int     num_ex_datess, /*	Short Exercise dates */
    long*   ex_dates_,     /*	Supposed to be sorted */
    int*    cal_dates,     /*	1: use ex_date as calibration date, 0: don't */
    char**  end_tenors,    /*	Tenors of the underlying instruments */
    long    end_dates,     /*	End date for diagonal */
    double* strikes_,      /*	Strikes */
    double* strikeS1s_,    /*	Strike1 for smile calib */
    double* strikeS2s_,    /*	Strike2 for smile calib */
    double* weights_,      /*	Weights on secondary instruments */

    CPD_DIAG_CALIB_PARAM params,

    /*	Model */
    LGMSV_CALIBPARAMS calib_params,

    int     iOne2F,
    double* dLambdaX,
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dAlphaTS,
    double* dLambdaEpsTS,
    double* dRhoTS,
    double  dTStar,
    double  dLGMAlpha,
    double  dLGMGamma,
    double  dLGMRho,
    double* dRho2TS,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    double** alpha,
    double** lambdaeps,
    double** rho,
    double** rho2,

    /*	Parameters */
    DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

Err LGMSVCalibApprox_Given_Lambda_and_Smile(
    /* Instrument informations	*/
    void* AllCalibInst,

    /* Model informations */
    void* model,

    /* Parameter of grids */
    void* AllParams);

Err LGMSVBumpLambdaModel(void* model_, int iIndex, double dNewValue);

Err LGMSVPricingShortInstruments(/* Instrument informations	*/
                                 void* AllCalibInst_,

                                 /* Model informations */
                                 void* model_,

                                 /* Parameter of grids */
                                 void* AllParams_,

                                 int     iIndex,
                                 double* dPrice);

Err LGMSVCalibApprox_Given_Lambda(
    /* Instrument informations	*/
    void* AllCalibInst_,

    /* Model informations */
    void* model_,

    /* Parameter of grids */
    void* AllParams_);
#endif
#endif