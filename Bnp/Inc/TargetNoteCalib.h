#ifndef _TARGETNOTE_CALIB_H_
#define _TARGETNOTE_CALIB_H_

#include "CPDcalib.h"
#include "TargetNoteProdStruct.h"

typedef struct
{
    int     num_ex_dates;
    long    end_date;
    int*    cal_date;
    char**  end_tenor;
    double* strike;
    double* strikeP1;
    double* strikeM1;
} DiagCalibInstStruct;

void TARN_calib_freeInst(DiagCalibInstStruct* ptrInst);

typedef struct
{
    /* input lambda TS -> calibrated tau */
    int     nlam;
    double* lam;
    double* lam_time;

    /* 2F parameters */
    double alpha, gamma, rho;
    int    nFactors;

    /* SV output */
    double* calalpha;
    double* callambdaeps;
    double* calrho;
    double* calrho2;
    double  tstar;

    /* calibrated vol, tau */
    int     num_sig;
    double* sig_time;
    double* sig;

    /* Vol, tau, smile for init und */
    double** dmVol;
    double** dmTau;
    int      num_tau;
    int      num_smile, smile_col;
    double** smile_datas;

    /* output of SV get */
    long*   sigma_dates;
    long*   smile_dates;
    long*   tau_dates;
    double* tau;

} TARN_lgm;

typedef struct
{
    LGMSV_NumerParams    NumerParams;
    DiagCalibInstStruct  Primary;
    DiagCalibInstStruct  Secondary;
    diag_calib_lm_params LM_params;

} TARN_CALIB_AUX;

char* TARN_calc_KO_prob_2F(TARN_Struct* tarn, TARN_AUX* aux, TARN_CALIB_AUX* calib_aux);

char* TARN_calib(TARN_Struct* tarn, TARN_AUX* aux);

#endif