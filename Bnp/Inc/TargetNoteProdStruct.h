#ifndef __CTN_PROD_STRUCT_H
#define __CTN_PROD_STRUCT_H

#include "CCFProdStruct.h"
#include "CPDCalib.h"
#include "LGMSVCalib.h"
#include "LGMSVCalibApprox.h"
#include "LGMSVClosedForm.h"
#include "LGMSVClosedFormApprox.h"
#include "LGMSVGrfn.h"
#include "LGMSVMC.h"
#include "LGMSVPDE.h"
#include "LGMSVUtil.h"
#include "utError.h"
#include "utTypes.h"

/* ------------------------------------------------------------------------------------------------------------------
 */
/* ENUMs */

/* Different types of coupons */
typedef enum
{
    TN_ZC,
    TN_SWAP,
    TN_FINAL,
    TN_SWAPKO
} TARGET_NOTE_COUPON;

/* Different models */
typedef enum
{
    TN_2F,
    TN_1F_SV,
    TN_2F_SV
} TARGET_NOTE_MODEL;

/* ENUMs */
/* ------------------------------------------------------------------------------------------------------------------
 */

typedef struct
{
    /* market information */
    long   lToday;
    char*  szYieldCurve;
    char*  szFundYieldCurve;
    char*  szVolCurve;
    char*  szVolCurveRefRate;
    char*  szCcy;
    double fx_fund_dom;           /* Fx fund/dom, 2bd fwd */
    long   fx_fund_dom_spot_date; /* FX spot date */
                                  /* The underlyings (set to "CAL" if calibration required) */
    char* szLGM2FUnd;
    char* szLGM1FUnd;
    char* szLGM1FSVUnd;
    char* szLGMSVUnd;
    /* Market Functions */
    char* (*getCashVol)(
        char*   szVolCurve,
        double  dStartDate,
        double  dEndDate,
        double  cashStrike,
        int     iZero,
        char*   szRefRate,
        double* dVol,
        double* dPower);
    char* (*getDF)(char* szYieldCurve, double dStart, double dEnd, double* dDF);
    /* EOD flags */
    int eodFixFlag; /*	0: do not check for fixings today; 1: check for fixings today */
    int eodPayFlag; /*	0: ignore cash-flows with payment date today; 1: include cash-flows */

} TARN_Market_Struct;

typedef struct
{
    /* The target note details */
    double             dTarget;
    TARGET_NOTE_COUPON couponType;
    int                bRepayCoupon;
    /* The coupon schedule */
    double       dCpnNotional;
    int          nCouponDates;
    double*      dvCoupon;
    long*        lvCpnStartDates;
    long*        lvCpnPayDates;
    int*         ivIsFloored;
    double*      dvFloor;
    int*         ivIsCapped;
    double*      dvCap;
    SrtBasisCode basisCpn;
    double       dFixedCvg;
    /* The floater */
    double*      dvGearing;
    long*        lvFltrFixingDates;
    long*        lvFltrStartDates;
    long*        lvFltrEndDates;
    double*      dvFltrSpread;
    SrtBasisCode basisFltr;
    double*      dvFltrHistFixings;
    /* The funding */
    int          fund_ccy; /* 0: domestic,  1: foreign */
    double       dFundingNotional;
    int          nFundDates;
    long*        lvFundFixingDates;
    long*        lvFundStartDates;
    long*        lvFundEndDates;
    double*      dvFundMargin;
    double*      dvFundSpread;
    char*        szBasisFund;
    char**       szvBasisFund;
    SrtBasisCode basisFund;
    double*      dvFundHistFixings;
} TARN_Deal_Struct;

typedef struct
{
    /* LGM2F model parameters */
    double  alpha;
    double  gamma;
    double  rho;
    int     nlam;
    double* lam_time;
    double* lam;
    /* SV Model Parameters */
    int     nsmilepar;
    double* smilepartime;
    double* alphaepsts;
    double* ldaepsts;
    double* rhoepsts;
    double* rho2epsts;
    double  tstar;
    /* 1FSV Parameters */
    int     nsmilepar1F;
    double* smilepartime1F;
    double* alphaepsts1F;
    double* ldaepsts1F;
    double* rhoepsts1F;
} TARN_Model_Struct;

typedef struct
{
    cpd_diag_calib_param prim_param;
    cpd_diag_calib_param sec_param;
    /* Instruments calibration options */
    int skip_start;
    int num_prob_paths;

    /* Volatility calibration options */

    int    nb_iter_LM;
    double precision_LM;
    /* LGM2F parameter options */
    int fix_lambda;
    /* LGM2V calibration options */
    LGMSV_CalibParams lgmsv_calib_params;
    /* The calibration instrument parameters */
    char* szPrimRefRate;
    char* szPrimFreq;
    char* szPrimBasis;
    char* szPrimTenor;
    char* szSecRefRate;
    char* szSecFreq;
    char* szSecBasis;
    char* szSecTenor;

    /*	SV Numerical params */
    int    iNbX;
    double iNbSigmaXLeft;
    double iNbSigmaXRight;
    double dIntegParam;
    int    iIntegMethod;
    double dVolLimit;
    int    iCalibLGM;
    double dMinStd;
    double dMaxStd;
    double numer_tstar;
} TARN_Calibration_Struct;

typedef struct
{
    /* pricings required */
    int iPrice2FCV;
    int iPriceSV;
    /* MC */
    double dSpread;
    int    iSpreadType; /* -1:  underprice coupon, 0:  even spread, 1:  overprice coupon */
    int    nStepT;
    long   num_paths;
    int    do_pecs;
    int    iFeeAdjust; /* 1:  adjust price for fee,		0: do not adjust price for fee */
} TARN_Pricing_Struct;

typedef struct
{
    /* Output:  PV information */
    int      iIsKnockedOut;
    int      iIsExpired;
    double   dCumulCoupon;
    double   dAccretedCoupon;
    double** dmLGM2F_PV;
    double** dmLGM1F_PV;
    double** dmLGM1FSV_PV;
    double** dmLGMSV_PV;
    double*  dvFee;
    /* KO info */
    double** dmKO_prob;
    int      nProb;
    /* Output:  calibrated underlyings */
    char* szTARN_LGM2F_UND;
    char* szTARN_LGM1F_UND;
    char* szTARN_LGM1FSV_UND;
    char* szTARN_LGMSV_UND;
    /* Output:  calbiration instrument data */
    cpd_calib_inst_data inst_data_lgm2F;
    cpd_calib_inst_data inst_data_lgm1F;
    cpd_calib_inst_data inst_data_lgm1FSV;
    cpd_calib_inst_data inst_data_lgmSV;
} TARN_Output_Struct;

typedef struct
{
    TARN_Market_Struct      market;
    TARN_Deal_Struct        deal;
    TARN_Calibration_Struct calibration;
    TARN_Model_Struct       model;
    TARN_Pricing_Struct     pricing;
    TARN_Output_Struct      output;
} TARN_Struct;

/* ------------------------------------------------------------------------------------------------------------------
 */
/* deal definition structures */

typedef struct
{
    /* Coupon */
    double* dvCouponCvg;
    double* dvCumulCvg; /* the value used to calculate the cumulative coupon */
    double* dvCouponPayDF;
    /* The funding */
    int*    ivFundStartIndex;
    int*    ivFundEndIndex;
    double* dvFundCvg;
    double* dvFundForStartDF;
    double* dvFundDomStartDF;
    double* dvFundDomPayDF;
    /* This history */
    int    i1stCpn; /* the first coupon to be used in the MC */
    double dHistFundPV;
    double dHistCpnPV;
    double dHistCumulCpn;
    /* whether today is an event date or not */
    int iIsTodayFixing;
    /* the fees paid to the funding leg if the deal KOs */
    double* dvKnockOutFee;

} TARN_AUX;

char* init_TARN_AUX(TARN_Struct* tarn, TARN_AUX* out_TARN);

void free_TARN_AUX(TARN_AUX* out_TARN);

/* deal definition structures */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
 */
/* global variable for storing path dependency */
typedef struct
{
    double             m_dTarget;
    double             m_dAdjustedTarget;  // Target - Spread
    double             m_dSpread;
    double             m_dCumulCoupon;
    double             m_dHistCumulCoupon;
    double             m_dCouponNotionalRepaid;
    double             m_dCouponNotional;
    double             m_dFundingNotional;
    TARGET_NOTE_COUPON m_couponType;  // 0:  pay normal coupons		1:   payment at end
    double m_dFinalCoupon;            // 0 if target not reached, accreted at fltr rate if TN_FINAL
    double m_dPrevAccrual;            // 1 + Fltr * Cvg
    double m_dFundingNotionalExchange;
    double m_dCouponNotionalExchange;
    //	int m_UseNotionalExchange;
    TARN_Struct* tarn;
} TargetNote;

void initTargetNote(TARN_Struct* tarn);
void resetTargetNote();

/* global variable for storing path dependency */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
 */
/* variable used for storing and passing market information */

typedef struct
{
    /* time steps */
    int     nStps;
    int*    is_event;
    double* date;
    double* time;
    /* events */
    int     nEvents;
    double* evt_tms;
    double* evt_dts;
    void**  void_prm;
    /* Model information */
    char*             szUnd;
    char*             und_name;
    TARGET_NOTE_MODEL tnModel;
    void*             ptrModel;
    /* covariance structure */
    double*** covar;
} TARN_MC_AUX;

void  init_TARN_MC_AUX(TARN_MC_AUX* mc_aux, char* szUnd);
void  free_TARN_MC_AUX(TARN_MC_AUX* mc_aux);
char* alloc_TARN_MC_AUX(TARN_Struct* tarn, TARN_AUX* aux, TARN_MC_AUX* mc_aux);

/* variable used for storing and passing market information */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
 */
/* 2F specific information */

typedef struct
{
    /* LGM2F parameters */
    double  lam;
    double  lam2;
    double* ifr;
    double* fwd1;
    double* fwd2;
    double* exp1;
    double* exp2;
    double* phi1;
    double* phi2;
    double* phi12;
    double* gam1_fwd;
    double* gam2_fwd;
    double* bond_pay;
    double* gam1_pay;
    double* gam2_pay;
    /* numeraire quantities */
    long   pay_date;
    double pay_time;
} TARN_MC_AUX_2F;

char* init_TARN_MC_AUX_2F(TARN_MC_AUX_2F** mc_aux_2f);
void  free_TARN_MC_AUX_2F(TARN_MC_AUX_2F* mc_aux_2f);
char* alloc_TARN_MC_AUX_2F(TARN_Struct* tarn, TARN_AUX* aux, TARN_MC_AUX* aux_mc);

/* 2F specific information */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
 */
/* variable for each event date */
// LGM2F
typedef struct
{
    // Funding parameters
    int     num_fund_df;  // number of df required
    double* fund_start_gam1;
    double* fund_start_gam2;
    double* fund_start_gam12;
    double* fund_end_gam1;
    double* fund_end_gam2;
    double* fund_end_gam12;
    // Fixed Coupon
    double fixed_pay_gam1;
    double fixed_pay_gam2;
    double fixed_pay_gam12;
    // Floater details
    double fltr_start_gam1;
    double fltr_start_gam2;
    double fltr_start_gam12;
    double fltr_end_gam1;
    double fltr_end_gam2;
    double fltr_end_gam12;
    // Numeraire
    double numeraire_gam1;
    double numeraire_gam2;
    double numeraire_gam12;
} TARN_EVENT_2F;

void init_TARN_EVENT_2F(TARN_EVENT_2F* targetNotePtr);
void free_TARN_EVENT_2F(TARN_EVENT_2F* targetNotePtr);

typedef struct
{
    // Funding parameters
    int     num_fund_df;  // number of df required
    double* fund_start_gam1;
    double* fund_start_gam2;
    double* fund_start_gam1_2;
    double* fund_start_gam2_2;
    double* fund_start_gam12;
    double* fund_end_gam1;
    double* fund_end_gam2;
    double* fund_end_gam1_2;
    double* fund_end_gam2_2;
    double* fund_end_gam12;
    // Fixed Coupon
    double fixed_pay_gam1;
    double fixed_pay_gam2;
    double fixed_pay_gam1_2;
    double fixed_pay_gam2_2;
    double fixed_pay_gam12;
    // Floater details
    double fltr_start_gam1;
    double fltr_start_gam2;
    double fltr_start_gam1_2;
    double fltr_start_gam2_2;
    double fltr_start_gam12;
    double fltr_end_gam1;
    double fltr_end_gam2;
    double fltr_end_gam1_2;
    double fltr_end_gam2_2;
    double fltr_end_gam12;
} TARN_EVENT_SV;

void init_TARN_EVENT_SV(TARN_EVENT_SV* targetNotePtr);
void free_TARN_EVENT_SV(TARN_EVENT_SV* targetNotePtr);

typedef struct
{
    // Funding parameters
    int     num_fund_df;        // number of df required
    double* fund_start;         // pointer to the start evolved df
    double* fund_end;           // pointer to the end evolved df
    double* fund_start_tms;     // times df required
    double* fund_start_dts;     // dates df required
    double* fund_start_logdff;  // df[ RequiredDate ] / df[ EventDate ]
    double* fund_end_tms;       // times df required
    double* fund_end_dts;       // dates df required
    double* fund_end_logdff;    // df[ RequiredDate ] / df[ EventDate ]
    double* fund_cvgMargin;     // 1 - coverage times margin
                                // Fixed Coupon
    double fixed_pay;
    double fixed_pay_tms;     // times df required
    double fixed_pay_dts;     // dates df required
    double fixed_pay_logdff;  // log ( df[ RequiredDate ] / df[ EventDate ] )
    double fixed_cvg;
    double fixed_coupon;
    double cumul_cvg;
    // Floater details
    double fltr_start_tms;     // times df required
    double fltr_start_dts;     // dates df required
    double fltr_start_logdff;  // df[ RequiredDate ] / df[ EventDate ]
    double fltr_end_tms;       // times df required
    double fltr_end_dts;       // dates df required
    double fltr_end_logdff;    // df[ RequiredDate ] / df[ EventDate ]
    double fltr_cvg;
    double fltr_gearing;
    double fltr_start;
    double fltr_end;
    double fltr_spread;
    int    iIsFloored;
    double dFloor;
    int    iIsCapped;
    double dCap;
    // Model Parameters
    void*  ptrModel;
    int    iModelType;  // 0:  LGM2F, 1:  LGM1FSV, 2:  LGM2FSV
    double numeraire_logdff;
    // Control
    double dKnockOutFee;  // Fee paid to the funding leg if the deal KOs at this event.
    int    nEventDates;
    int    iEventDate;
    int    iIsFinal;  // 1 if final period

} TARN_EVENT;

void  init_TARN_EVENT(TARN_EVENT* event);
void  free_TARN_EVENT(TARN_EVENT* event);
char* alloc_TARN_EVENT(TARN_EVENT** event, TARN_MC_AUX* aux);

/* variable for each event date */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
 */
/* Different coupon PV calculators */

char* TargetNoteFunding(TARN_EVENT* targetNotePtr, double* res);

char* TargetNoteCoupon_ZC(
    TARN_EVENT* targetNotePtr,
    double      dCoupon,
    double      dNotionalRepay,
    double*     res,  // res[0] = fund, res[1] = coupon
    int*        stop_path);

char* TargetNoteCoupon_SWAP(
    TARN_EVENT* targetNotePtr,
    double      dCoupon,
    double      dNotionalRepay,
    double*     res,  // res[0] = fund, res[1] = coupon
    int*        stop_path);

char* TargetNoteCoupon_FINAL(
    TARN_EVENT* targetNotePtr,
    double      dCoupon,
    double      dFloater,
    double      dNotionalRepay,
    double*     res,  // res[0] = fund, res[1] = coupon
    int*        stop_path);

char* TargetNoteCoupon_SWAPKO(
    TARN_EVENT* event,
    double      dCoupon,
    double      dNotionalRepay,
    double*     res,  // res[0] = fund, res[1] = coupon
    int*        stop_path);

char* TargetNoteCoupon(
    TARN_EVENT* targetNotePtr,
    double*     res,  // res[0] = fund, res[1] = coupon
    int*        stop_path);

/* Different coupon PV calculators */
/* ------------------------------------------------------------------------------------------------------------------
 */

#endif