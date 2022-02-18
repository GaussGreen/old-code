/* ===================================================================================
   FILENAME:      InflationPricing.h / P.A

   PURPOSE:
   =================================================================================== */

#ifndef __INFLATIONPRICING_H
#define __INFLATIONPRICING_H

#include "math.h"
#include "srt_h_types.h"
#include "swp_h_df.h"

/* ------------------------- */
/* For df wrapping           */
/* ------------------------- */
typedef struct
{
    int iIsStripping; /* If we strip use the ptr of function if not use a yc name */
    double (*StripDf)(long lSettle, long ldate, void* BetaPtr); /* Ptr of function on bdf */
    char* CurveName;                                            /* Yield curve name */
    void* BetaPtr;
} WrapInfo;

/* ------------------------- */
/* Inflation Model structure */
/* ------------------------- */

typedef enum
{
    Infl3F
} InflationModelType;

typedef struct
{
    long               lToday;
    InflationModelType iModelType;

    /* Specific Model Parameter */
    void* ModelVolTs;
} InflationModel;

/* 3F Volatility TermStructure */
typedef struct
{
    /* Rates Time */
    long    lRateNbDate;
    double* dRateSigTime;

    /* Nominal TermStruct */
    double* dDomSig;
    double  dDomLambda;

    /* Real TermStruct */
    double* dForSig;
    double  dForLambda;

    /* CPI TermStruct */
    long    lFxNbDate;
    double* dFxSigtime;
    double* dFxSig;

    /* Correlation TermStruct */
    long    lCorrNbDate;
    double* dRhoTime;
    double* dRhoDomFor;
    double* dRhoDomFx;
    double* dRhoForFx;

} Inflation_3FVolTs;

/* ---------------- */
/* CPIMKT structure */
/* ---------------- */

/* Fixings term struct */
typedef struct
{
    long    lNbDate;
    long*   lFixingDate;
    double* dFixValue;

} Inflation_CPIFixing;

/* CPIMKT structure */
typedef struct
{
    Inflation_CPIFixing CPIFixing;
    long                lTodayDate;
    double              dTodayCPI;

} Inflation_CPIMkt;

/* -----------------*/
/* Coupon structure */
/* -----------------*/
typedef struct
{
    /* Description of one Coupon */
    /* at Tpay we receive if floored  Coupon*(max((alpha * CPI(Tf1) + beta * CPI(Tf2))/DRIRef +
     * gamma, Floor)) */
    /*					  if not     Coupon*((alpha * CPI(Tf1) + beta * CPI(Tf2))/DRIRef + gamma)
     */

    long   lTheoPayDate; /* Theoritical PayDate */
    long   lPayDate;     /* PayDate			 */
    double dCouponValue; /* CouponValue		 */

    double dDRIRef; /* DRI Ref */

    long lFixDateCPI1; /* fixing date of the CPI1 */
    long lFixDateCPI2; /* fixing date of the CPI2 */

    double dAlpha; /* DRI coef on CPI1 */
    double dBeta;  /* DRI coef on CPI2 */
    double dGamma; /* Constant */

    double dCvxtyAdj1; /* Convexty adj such that E under QTp(CPI(Tf)) = FwdCPI(t,Tf)* CvxtyAdj */
    double dCvxtyAdj2; /* Convexty adj such that E under QTp(CPI(Tf)) = FwdCPI(t,Tf)* CvxtyAdj */

    int    iIsFloored; /* Is floored coupon */
    double dFloor;     /* floor strike */
    double dVol1;      /* Cumulated vol of the CPI1 */
    double dVol2;      /* Cumulated vol of the CPI2 */
    double dCorr12;    /* Cumulated correlation between the CPI1 and the CPI2 */

    /* Temporary Pricing info used for assetswap */
    double dFwdPV; /* Nominal forward PV at PayDate */

    /* For infos */
    double dFwdCPI1;    /* Fwd CPI1 at Tf */
    double dFwdCPI2;    /* Fwd CPI1 at Tf */
    double dFwdPVFloor; /* forward PV of the floor at PayDate */
    double dDfPayDate;  /* Discount factor nominal at paydate */

} InflationCoupon;

/* ------------------- */
/* CashFlows structure */
/* ------------------- */
typedef struct
{
    int              iNbCoupon; /* Number of coupons */
    InflationCoupon* Coupon;    /* Coupon array */

    /* Temporary */
    double dTodayPV; /* Today's PV */

} InflationCashFlows;

/* ---------------*/
/* AssetSwap type */
/* ---------------*/
typedef enum
{
    FlatNotional,
    FlatIndexedNotional,
    IndexedNotional,
    FlatDirtyNotional,
    IndexedDirtyNotional
} InflationAssetSwapType;

/* Files */
/* ------------------------------------------------------------------------------------------------
        Inflation_GetFixing

        get the fixing from the array extrapolation and interpolation flat
        the array of fix date should be stored in decreasing date
  -------------------------------------------------------------------------------------------------
*/
double Inflation_GetFixing(long lValDate, Inflation_CPIFixing* CPIFixing);

/* ------------------------------------------------------------------------------------------------
        Inflation_FwdCPI

        calculate the classical FwdCPI(t,ValDate) = E(under QTval) (FwdCPI(ValDate,ValDate)
  -------------------------------------------------------------------------------------------------
*/
double Inflation_FwdCPI(
    long              lValDate,
    Inflation_CPIMkt* CPIMkt,
    const char*       NomCurve, /* Nominal Curve name*/
    WrapInfo*         RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */
);

/* ------------------------------------------------------------------------------------------------
        Inflation_FwdDRI_givenAdj

        calculate the E (under Qtp) (alpha*CPI(Tf1)+beta*CPI(Tf2))
        where the convexity adjustement are given
  -------------------------------------------------------------------------------------------------
*/
double Inflation_FwdDRI_givenAdj(/* Inputs */
                                 long              lFixDate1,
                                 long              lFixDate2,
                                 double            dAlpha,
                                 double            dBeta,
                                 double            dCvxtyAdj1,
                                 double            dCvxtyAdj2,
                                 Inflation_CPIMkt* CPIMkt,   /* CPI Market */
                                 const char*       NomCurve, /* Nominal Curve name*/
                                 WrapInfo* RealCurve, /* Real curve either a name if WestWrapper or
                                                         a c++ object ptr */

                                 /* Other Outputs */
                                 double* ptr_dFwdCPI1,
                                 double* ptr_dFwdCPI2);

/* -------------------------------------------------------------------------------------------------------
        Inflation_FillAdjInCashFlows

        Fill the convexity adj of all CPI and vols in the cash flows structure

        Assume that alpha, beta FixDate1, FixDate2, and paydate have been filled
  ---------------------------------------------------------------------------------------------------------
*/
Err Inflation_FillAdjInCashFlows(                       /* Input */
                                 InflationModel* Model, /*	Model used */
                                 /* Output */
                                 InflationCashFlows* CashFlows /* Description of the cash flows */);

/* -------------------------------------------------------------------------------------------------------
        Inflation_BondValue

        Calculate the fwd dirtyprice in nominal at settlement date of the bond
        BondCashFlows and AssetSwapCashFlows have been initialised

  ---------------------------------------------------------------------------------------------------------
*/
Err Inflation_BondValue(/* Inputs */
                        long                lSettleDate,
                        InflationCashFlows* BondCashFlows, /* Description of Bond cash flows */

                        InflationAssetSwapType iAssetSwapType,  /* Assetswap type */
                        InflationCashFlows* AssetSwapCashFlows, /* Description of Bond cash flows */

                        Inflation_CPIMkt* CPIMkt,   /* CPI Market */
                        const char*       NomCurve, /* Nominal Curve name*/
                        WrapInfo* RealCurve, /* Real curve either a name if WestWrapper or  a c++
                                                object ptr */

                        /* Outputs */
                        double* ptr_DirtyPrice, /* Fwd Dirty price at settle date */
                        double* ptr_DfSettle /* Disount factor nominal at settlement date */);

/* -------------------------------------------------------------------------------------------------------
        Inflation_ValueCashFlow

        calculate the nominal cash flows
        Fill the temporary fwdPV infos in the cash flows

        Description of the Cash flows
                 at Tpay we receive if floored  Coupon*(max(alpha * CPI(Tf1) + beta * CPI(Tf2) +
  gamma, Floor)) if not      Coupon*(alpha * CPI(Tf1) + beta * CPI(Tf2) + gamma)
  ---------------------------------------------------------------------------------------------------------
*/
double Inflation_ValueCashFlow(
    InflationCashFlows* CashFlows, /* Description of the cash flows */
    Inflation_CPIMkt*   CPIMkt,    /* CPI Market */
    const char*         NomCurve,  /* Nominal Curve name*/
    WrapInfo* RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */);

/* -------------------------------------------------------------------------------------------------------
        Inflation_InitNotionalAssetSwap

        Assume that
                CashFlows->dGamma contains the fundingcashspread * coverage
                Modify it to put the notional depending on the type of the AS
  ---------------------------------------------------------------------------------------------------------
*/
Err		Inflation_InitNotionalAssetSwap(/* Input */
										  InflationAssetSwapType	iAssetSwapType,			/* Assetswap type */	
										  InflationCashFlows		*BondCashFlows,			/* Description of Bond cash flows */
										  double					dTodayCBP,				/* dTodayCBP of the bond */	
										  double					dDf_Tsettle,			/* df nom at settle date */
										  /* OutPut */
										  InflationCashFlows		*AssetSwapCashFlows		/* Description of the assetswap cash flows	 */ );

#endif
