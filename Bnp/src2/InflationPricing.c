/* ===================================================================================
   FILENAME:      InflationPricing.cpp / P.A

   PURPOSE:
   =================================================================================== */
#include "InflationPricing.h"

#include "InflationUtils.h"
#include "RainbowOpt.h"
#include "opfnctns.h"

/* ------------------------------------	*/
/* Model independant pricing functions	*/
/* Low levels functions	Level 0			*/
/* ------------------------------------	*/

/* -----------------------------------------------------------------------------------------------
        Notations
                Cash Fx (t) = CPI (t)
   -----------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------
        Inflation_GetFixing

        get the fixing from the array extrapolation and interpolation flat
        the array of fix date should be stored in decreasing date
  -------------------------------------------------------------------------------------------------
*/
double Inflation_GetFixing(long lValDate, Inflation_CPIFixing* CPIFixing)
{
    int iNum;

    /* find the first date before lValDate */
    iNum = 0;
    while ((iNum < CPIFixing->lNbDate) && (lValDate < CPIFixing->lFixingDate[iNum]))
        iNum++;

    return CPIFixing->dFixValue[iNum];
}

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
)
{
    double dFwdCPI;
    long   lToday;

    lToday = CPIMkt->lTodayDate;

    if (lValDate < lToday)
    {
        /* Get the CPI from the HistFixing */
        dFwdCPI = Inflation_GetFixing(lValDate, &CPIMkt->CPIFixing);
    }
    else
    {
        /* Calculation of the FwdCPI */
        dFwdCPI = CPIMkt->dTodayCPI * dfReal(lToday, lValDate, RealCurve) /
                  swp_f_df(lToday, lValDate, NomCurve);
    }

    return dFwdCPI;
}

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
                                 double* ptr_dFwdCPI2)
{
    double dFwdCPI1, dFwdCPI2;
    double dFwdDRI;

    /* Calculation of the E (under Qtp) (S(Tf1)) */
    /* The adjustment have already been calculated */
    if (dAlpha != 0)
    {
        /* Calculation of the E (under Qtf1) (CPI(Tf1)) */
        dFwdCPI1 = Inflation_FwdCPI(
            lFixDate1,
            CPIMkt,
            NomCurve, /* Nominal Curve name*/
            RealCurve /* Real curve name */);

        /* Chgt of proba */
        dFwdCPI1 *= dCvxtyAdj1;
    }
    else
    {
        /* Default value if not used */
        dFwdCPI1 = 0.0;
    }

    /* Calculation of the E (under Qtp) (S(Tf2)) */
    /* The adjustment have already been calculated */
    if (dBeta != 0)
    {
        /* Calculation of the E (under Qtf2) (CPI(Tf2)) */
        dFwdCPI2 = Inflation_FwdCPI(
            lFixDate2,
            CPIMkt,
            NomCurve, /* Nominal Curve name*/
            RealCurve /* Real curve name */);

        /* Chgt of proba */
        dFwdCPI2 *= dCvxtyAdj2;
    }
    else
    {
        /* Default value if not used */
        dFwdCPI2 = 0.0;
    }

    /* Calculation of the FwdDRI */
    dFwdDRI = dAlpha * dFwdCPI1 + dBeta * dFwdCPI2;

    /* Return output */
    *ptr_dFwdCPI1 = dFwdCPI1;
    *ptr_dFwdCPI2 = dFwdCPI2;

    return dFwdDRI;
}

/* -------------------------------------------------------------------------------------------------------
        Inflation_FwdCouponValue

        calculate E(under Qtpay coupon)(Coupon)

        Description of the coupon
                 at Tpay we receive if floored  Coupon*(max(alpha * CPI(Tf1) + beta * CPI(Tf2) +
  gamma, Floor)) if not      Coupon*(alpha * CPI(Tf1) + beta * CPI(Tf2) + gamma)
  ---------------------------------------------------------------------------------------------------------
*/
double Inflation_FwdCouponValue(
    InflationCoupon*  Coupon,    /* Coupon */
    Inflation_CPIMkt* CPIMkt,    /* CPI Market */
    const char*       NomCurve,  /* Nominal Curve name*/
    WrapInfo*         RealCurve, /* Real curve either a name if WestWrapper or  a c++ object ptr */

    /* Other OutPuts */
    double* ptr_dFwdPVFloor, /* forward PV of the floor at PayDate */
    double* ptr_dFwdCPI1,    /* Fwd CPI1 at Tf */
    double* ptr_dFwdCPI2 /* Fwd CPI1 at Tf */)
{
    double dFwdDRI;
    double dTimeToExpiry, dFloor;
    double dPvIV, dPvFloor;
    long   lToday;
    double dFwdCPI1, dFwdCPI2;

    /* Init */
    dPvIV    = 0.0;
    dPvFloor = 0.0;
    dFloor   = 0;
    lToday   = CPIMkt->lTodayDate;

    /* Calculation of the cash nominal PV */

    /* Calculation of the fwd DRI (t, Tp, Tfix1, Tfix2)				*/
    /* = E (under QTp) (alpha*FWDCPI(Tf1,Tf1)+ beta*FWDCPI(Tf2,Tf2) */
    dFwdDRI = Inflation_FwdDRI_givenAdj(
        /* Inputs */
        Coupon->lFixDateCPI1,
        Coupon->lFixDateCPI2,
        Coupon->dAlpha,
        Coupon->dBeta,
        Coupon->dCvxtyAdj1,
        Coupon->dCvxtyAdj2,
        CPIMkt,    /* CPI Market */
        NomCurve,  /* Nominal Curve name*/
        RealCurve, /* Real curve name */

        /* Others OutPuts */
        &dFwdCPI1,
        &dFwdCPI2);

    /* Calculation of the floor on DRI */
    /*  E (under Qtpay) max(DRI(Tfix1,Tfix2)/DRIRef + gamma, Floor)		*/
    /* = E (under Qtpay) (  DRI(Tfix1,Tfix2)/DRIRef + gamma + max((Floor-gamma) -
     * DRI(Tfix1,Tfix2)/DRIRef , 0)*/

    /* Pricing of E (under Qtpay) max((Floor-gamma)*DRIRef - DRI(Tfix1,Tfix2), 0) */
    if (Coupon->iIsFloored)
    {
        dTimeToExpiry =
            ((double)(max(Coupon->lFixDateCPI1, Coupon->lFixDateCPI2) - lToday)) / DAYSINYEAR;

        if (dTimeToExpiry >= 0)
        {
            /* Spread Option pricing */
            /* E[max(K - (nx*Xt + ny*Yt), 0)] with Xt, Yt lognormal (using 1D numerical integration)
             */

            OptSpreadNew(
                dFwdCPI1,                          /* fwdx,	*/
                Coupon->dAlpha / Coupon->dDRIRef,  /* nx,		*/
                Coupon->dVol1,                     /* sigx,	*/
                dFwdCPI2,                          /* fwdy,	*/
                Coupon->dBeta / Coupon->dDRIRef,   /* ny,		*/
                Coupon->dVol2,                     /* sigy,	*/
                (Coupon->dFloor - Coupon->dGamma), /* K,		*/
                dTimeToExpiry,                     /* mat,		*/
                Coupon->dCorr12,                   /* rho,		*/
                SRT_PUT,
                &dFloor);
        }
        else
        {
            dFloor = max((Coupon->dFloor - Coupon->dGamma) - dFwdDRI / Coupon->dDRIRef, 0.0);
        }
    }

    /* Calculation of the Coupon Pv */
    dPvFloor = Coupon->dCouponValue * dFloor;
    dPvIV    = Coupon->dCouponValue * (dFwdDRI / Coupon->dDRIRef + Coupon->dGamma);

    *ptr_dFwdPVFloor = dPvFloor;
    *ptr_dFwdCPI1    = dFwdCPI1;
    *ptr_dFwdCPI2    = dFwdCPI2;

    /* Return results */
    return dPvFloor + dPvIV;
}

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
    WrapInfo*           RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */
)
{
    double dFwdCouponValue, dDf_Tpay;
    double dPv;
    int    iNum;
    long   lToday;
    double dFwdPVFloor, dFwdCPI1, dFwdCPI2;

    InflationCoupon* Coupon = NULL;

    /* Init */
    dPv    = 0.0;
    lToday = CPIMkt->lTodayDate;

    /* Calculation of the cash nominal PV */
    for (iNum = 0; iNum < CashFlows->iNbCoupon; iNum++)
    {
        /* Get the iNum coupon */
        Coupon = &(CashFlows->Coupon[iNum]);

        /* Get the Fwd PV of the coupon under QTpaycoupon */
        dFwdCouponValue = Inflation_FwdCouponValue(
            Coupon,       /* Coupon */
            CPIMkt,       /* CPI Market */
            NomCurve,     /* Nominal Curve name*/
            RealCurve,    /* Real curve either a name if WestWrapper or  a c++ object ptr */
            &dFwdPVFloor, /* FwdPVFloor */
            &dFwdCPI1,    /* FwdCPI1 */
            &dFwdCPI2 /* FwdCPI2 */);

        /* Calculation of the Coupon Pv */
        dDf_Tpay = swp_f_df(lToday, Coupon->lPayDate, NomCurve);
        dPv += dFwdCouponValue * dDf_Tpay;

        /* Save value in the structure */
        Coupon->dFwdPV      = dFwdCouponValue;
        Coupon->dFwdCPI1    = dFwdCPI1;
        Coupon->dFwdCPI2    = dFwdCPI2;
        Coupon->dFwdPVFloor = dFwdPVFloor;
        Coupon->dDfPayDate  = dDf_Tpay;
    }

    /* Put the today's Pv in the InflationCashFlows */
    CashFlows->dTodayPV = dPv;

    /* Return results */
    return dPv;
}

/* ------------------------------------	*/
/* Low levels functions	Level 1			*/
/* ------------------------------------	*/

/* -------------------------------------------------------------------------------------------------------
        Inflation_OneCouponFillAdjVols

        Fill the convexity adj of all CPI in the coupon and all the vols

        Assume that alpha, beta FixDate1, FixDate2, and paydate have been filled
  ---------------------------------------------------------------------------------------------------------
*/
Err Inflation_OneCouponFillAdjVols(                       /* Input */
                                   InflationModel* Model, /*	Model used */
                                   /* Output */
                                   InflationCoupon* Coupon /* Description of the cash flows */)
{
    Err err = NULL;

    /* Fill the convexty adjustment */

    /* Adjustement of FwdCPI 1 */
    if (Coupon->dAlpha != 0.0)
    {
        /* Calculation of the convexty adjustement */
        err = Inflation_GetAdjFwdCPI(
            Coupon->lFixDateCPI1, /*	Fixing date of the cash CPI */
            Coupon->lPayDate,     /*	Pay date  */
            Model,

            /* Output */
            &Coupon->dCvxtyAdj1);

        if (err)
            return err;
    }
    else
    {
        Coupon->dCvxtyAdj1 = 1.0;
    }

    /* Adjustement of FwdCPI 2 */
    if (Coupon->dBeta != 0.0)
    {
        /* Calculation of the convexty adjustement */
        err = Inflation_GetAdjFwdCPI(
            Coupon->lFixDateCPI2, /*	Fixing date of the cash CPI */
            Coupon->lPayDate,     /*	Pay date  */
            Model,

            /* Output */
            &Coupon->dCvxtyAdj2);

        if (err)
            return err;
    }
    else
    {
        Coupon->dCvxtyAdj2 = 1.0;
    }

    /* Fill the vols information if needed */
    if (Coupon->iIsFloored)
    {
        err = Inflation_GetCorrCPI1CPI2(               /* Inputs */
                                        Model->lToday, /* Start date */
                                        max(Coupon->lFixDateCPI1,
                                            Coupon->lFixDateCPI2), /* End date */

                                        Coupon->lFixDateCPI1, /* Val date1 */
                                        Coupon->lFixDateCPI2, /* Val date2 */

                                        Model, /*	Model used */

                                        /* Outputs */
                                        &(Coupon->dCorr12),
                                        &(Coupon->dVol1),
                                        &(Coupon->dVol2));
    }
    else
    {
        /* Default value */
        Coupon->dVol1   = 0.0;
        Coupon->dVol2   = 0.0;
        Coupon->dCorr12 = 0.0;
    }

    /* Return */
    return err;
}

/* -------------------------------------------------------------------------------------------------------
        Inflation_FillAdjInCashFlows

        Fill the convexity adj of all CPI and vols in the cash flows structure

        Assume that alpha, beta FixDate1, FixDate2, and paydate have been filled
  ---------------------------------------------------------------------------------------------------------
*/
Err Inflation_FillAdjInCashFlows(                       /* Input */
                                 InflationModel* Model, /*	Model used */
                                 /* Output */
                                 InflationCashFlows* CashFlows /* Description of the cash flows */)
{
    int iNum;
    Err err = NULL;

    for (iNum = 0; iNum < CashFlows->iNbCoupon; iNum++)
    {
        err = Inflation_OneCouponFillAdjVols(       /* Input */
                                             Model, /*	Model used */
                                             /* Output */
                                             &(CashFlows->Coupon
                                                   [iNum]) /* Description of the cash flows */);
        if (err)
            return err;
    }

    /* Return */
    return err;
}

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
										  InflationCashFlows		*AssetSwapCashFlows		/* Description of the assetswap cash flows	 */ )
{
    int iNum;
    Err err = NULL;

    double dBondLegTodayPv, dAssetSwapLegTodayPV, dSum, dDirtyPrice, dfTi_m1;

    /* Fill the Notional */
    switch (iAssetSwapType)
    {
    case FlatNotional:
        /* Put a flat notional = 1 */
        for (iNum = 0; iNum < AssetSwapCashFlows->iNbCoupon; iNum++)
        {
            AssetSwapCashFlows->Coupon[iNum].dCouponValue = 1.0;
        }

        break;

    case FlatIndexedNotional:
        /* Put a notional = CBP(today)/DIRRef */
        for (iNum = 0; iNum < AssetSwapCashFlows->iNbCoupon; iNum++)
        {
            AssetSwapCashFlows->Coupon[iNum].dCouponValue =
                dTodayCBP / BondCashFlows->Coupon[0].dDRIRef;
        }

        break;

    case IndexedNotional:
        /* Put a notional = E(DIR/DIRRef) * dDf_Tsettle  !!!!*/
        for (iNum = 0; iNum < AssetSwapCashFlows->iNbCoupon - 1; iNum++)
        {
            if (BondCashFlows->Coupon[iNum].dCouponValue == 0.0)
            {
                AssetSwapCashFlows->Coupon[iNum].dCouponValue = 0.0;
            }
            else
            {
                AssetSwapCashFlows->Coupon[iNum].dCouponValue =
                    dDf_Tsettle *
                    (BondCashFlows->Coupon[iNum].dFwdPV - BondCashFlows->Coupon[iNum].dFwdPVFloor) /
                    BondCashFlows->Coupon[iNum].dCouponValue;
            }
        }

        /* Special Case for the last coupon containing the redemption */
        iNum = AssetSwapCashFlows->iNbCoupon - 1;
        if (BondCashFlows->Coupon[iNum].dCouponValue +
                BondCashFlows->Coupon[iNum + 1].dCouponValue ==
            0.0)
        {
            AssetSwapCashFlows->Coupon[iNum].dCouponValue = 0.0;
        }
        else
        {
            AssetSwapCashFlows->Coupon[iNum].dCouponValue =
                dDf_Tsettle *
                (BondCashFlows->Coupon[iNum].dFwdPV + BondCashFlows->Coupon[iNum + 1].dFwdPV -
                 BondCashFlows->Coupon[iNum].dFwdPVFloor -
                 BondCashFlows->Coupon[iNum + 1].dFwdPVFloor) /
                (BondCashFlows->Coupon[iNum].dCouponValue +
                 BondCashFlows->Coupon[iNum + 1].dCouponValue);
        }

        break;

    case FlatDirtyNotional:

        /* Put a notional = dirty price */

        /* Calculation of the Dirty price if the notional of the asset swap is the dirty price */
        dBondLegTodayPv = BondCashFlows->dTodayPV;

        /* Calculation of the today asset swap PV if Notional = 1 */
        dSum = 0;
        for (iNum = 0; iNum < AssetSwapCashFlows->iNbCoupon; iNum++)
        {
            dSum +=
                AssetSwapCashFlows->Coupon[iNum].dGamma * BondCashFlows->Coupon[iNum].dDfPayDate;
        }

        dDirtyPrice = BondCashFlows->dTodayPV / (dSum + dDf_Tsettle);

        /* Put the notional  */
        for (iNum = 0; iNum < AssetSwapCashFlows->iNbCoupon; iNum++)
        {
            AssetSwapCashFlows->Coupon[iNum].dCouponValue = dDirtyPrice;
        }

        break;

    case IndexedDirtyNotional:

        /* Put notional(i) = fwd dirty price at Ti-1 */

        /* PV of the last coupon */
        dBondLegTodayPv = BondCashFlows->Coupon[AssetSwapCashFlows->iNbCoupon].dFwdPV *
                          BondCashFlows->Coupon[AssetSwapCashFlows->iNbCoupon].dDfPayDate;
        dAssetSwapLegTodayPV = 0.0;

        for (iNum = AssetSwapCashFlows->iNbCoupon - 1; iNum >= 0; iNum--)
        {
            /* Calculate the today's PV of Bond leg from TiNum */
            dBondLegTodayPv +=
                BondCashFlows->Coupon[iNum].dFwdPV * BondCashFlows->Coupon[iNum].dDfPayDate;

            /* Df at TiNum-1 */
            if (iNum == 0)
                dfTi_m1 = dDf_Tsettle;
            else
                dfTi_m1 = BondCashFlows->Coupon[iNum - 1].dDfPayDate;

            /* Calculate the dirty price at TiNum-1 */
            dDirtyPrice =
                (dBondLegTodayPv - dAssetSwapLegTodayPV) /
                (AssetSwapCashFlows->Coupon[iNum].dGamma * BondCashFlows->Coupon[iNum].dDfPayDate +
                 dfTi_m1);

            /* Update the notional of the asset swap */
            AssetSwapCashFlows->Coupon[iNum].dCouponValue = dDirtyPrice;

            /* Update the today pv of the asset swap leg */
            dAssetSwapLegTodayPV += AssetSwapCashFlows->Coupon[iNum].dGamma *
                                    AssetSwapCashFlows->Coupon[iNum].dCouponValue *
                                    BondCashFlows->Coupon[iNum].dDfPayDate;
        }

        break;

    default:
        err = "AssetSwapType not yet implemented ";
        return err;
    }

    /* Return result */
    return err;
}

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
                        double* ptr_DfSettle /* Disount factor nominal at settlement date */)
{
    double dPvBondLeg, dPvAssetSwapLeg;
    double dDf_Tsettle, dTodayCBP;
    Err    err = NULL;

    /* Calculation of the df at settle */
    dDf_Tsettle = swp_f_df(CPIMkt->lTodayDate, lSettleDate, NomCurve);

    /* Get dTodayCBP */
    dTodayCBP = Inflation_FwdCPI(
        CPIMkt->lTodayDate,
        CPIMkt,
        NomCurve, /* Nominal Curve name*/
        RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */
    );

    /* Pricing Bond leg cash flows seen as of today */
    /* All the intermediate PV seen as of today have been stored in the BondCashFlows */
    dPvBondLeg = Inflation_ValueCashFlow(
        BondCashFlows, /* Description of the cash flows */
        CPIMkt,        /* CPI Market */
        NomCurve,      /* Nominal Curve name*/
        RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */);
    /* Init notional of the assetswap */
    err = Inflation_InitNotionalAssetSwap(/* Input */
										  iAssetSwapType,			/* Assetswap type */
										  BondCashFlows,			/* Contain the fwd PV */
										  dTodayCBP,				/* dTodayCBP of the bond */	
										  dDf_Tsettle,
										  /* OutPut */
										  AssetSwapCashFlows		/* Description of the assetswap cash flows	 */ );
    if (err)
        return err;

    /* Pricing AssetSwap leg cash flows seen as of today */
    dPvAssetSwapLeg = Inflation_ValueCashFlow(
        AssetSwapCashFlows, /* Description of the cash flows */
        CPIMkt,             /* CPI Market */
        NomCurve,           /* Nominal Curve name*/
        RealCurve /* Real curve either a name if WestWrapper or  a c++ object ptr */);

    /* Calculation of the fwd dirty price at settle date */
    *ptr_DirtyPrice = (dPvBondLeg - dPvAssetSwapLeg) / dDf_Tsettle;
    *ptr_DfSettle   = dDf_Tsettle;

    /* Return results */
    return err;
}
