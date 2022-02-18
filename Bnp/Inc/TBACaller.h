#ifndef TBACALLER_H
#define TBACALLER_H
#include "MBSPPFuncs.h"
#include "MBSPTCalib.h"
#include "MBSPTProdStruct.h"
#include "MBSPrepay.h"
#include "TBAUtil.h"

char* mbsdf(
    char*   data_dir,
    long*   valuation,
    long*   sptdays,
    double* refisprd,
    double* tsycurve,
    double* swpcurve,
    double* volcurve,
    long*   volexp,
    double* meanrev,
    long*   choicen,
    long*   swapmat,
    long*   dealmat,
    double* coupon,
    double* inoas,
    long*   paydelays,
    double* pptweaks,
    double* inprices,
    double* outputs);

char* impp(
    TERM_DATA*  term_data,
    TREE_DATA*  tree_data,
    DEAL_DATA*  deal_data,
    HISTORY*    history,
    MBSDENSITY* density,
    TRIGGER*    trigger,
    DEAL*       deal,
    PREPAY*     prepay,
    double*     price,
    double*     HedgePrice);

char* Tweak(
    double       Option,      /* Original option price */
    TERM_DATA*   term_data,   /* Structure of term structure data */
    TREE_DATA*   tree_data,   /* Structure of tree data */
    DEAL_DATA*   deal_data,   /* Structure of deal data */
    PREPAY_DATA* prepay_data, /* Structure of prepayment data */
    HISTORY*     history,
    MBSDENSITY*  density,
    TRIGGER*     trigger,
    DEAL*        deal,
    PREPAY*      prepay,
    double*      price);

void PrintTweak(
    double  Option,    /* Option price */
    double* RateTweak, /* Interest rate tweak of the option */
    double* Duration); /* Duration of benchmark instruments */

void Main_Prepay(
    double** Amort,    /* Array of smms of various speed groups in the interest rate lattice */
    double*  ParYield, /* Par Yield in the lattice at the current period */
    int      CouponNb, /* Coupon Number: = deal_data->Term at the end of the tree */
    long     Month,    /* Current month of the year (1 to 12) */
    int      NbPer,    /* Current time period */
    int      N,
    MBS_Prepay_Engine* prepay_engine,
    TREE_DATA*         tree_data /* Structure of tree data */
);

int TweakSCurve(DEAL_DATA* deal_data, long today);  // 1 means succcessful

#endif