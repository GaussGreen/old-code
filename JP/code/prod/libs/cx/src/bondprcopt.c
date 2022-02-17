/*
*****************************************************************************
** Calculation implementation for bond price option pricer.
*****************************************************************************
*/

#include "bondprcopt.h"

#include <math.h>

#include <crxflow/include/crxutilio.h>
#include <crxflow/include/crxwrpio.h>
#include <crxflow/include/crxmacros.h>

#include <alib/bondcnst.h>
#include <alib/bondfwd.h>
#include <alib/bondsens.h>
#include <alib/lintrp.h>
#include <alib/normopt.h>
#include <alib/optprop.h>

#define YYYYMMDD(x) x>0?CrxTDate2DrDate(x):0

/*f
***************************************************************************
** Main calculator
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalc(
CrxTBondPriceOption*   deal,             /* (I) */
long                   distType,         /* (I) */
TDate                  today,            /* (I) */
CrxTBondPriceVolCurve* volCurve,         /* (I) */
CrxTBondPrice*         bondPrice,        /* (I) */
CrxTBondRepoCurve*     repoCurve,        /* (I) */
TCurve*                discCurve         /* (I) */
)
{
    static char routine[] = "CrxBondPriceOptionCalc";

    CrxTBondPriceOptionCalc *calc = NULL;

    TBond             *bond = NULL;
    TBondTradeData    *spotTrade = NULL;
    TBondTradeData    *fwdTrade = NULL;
    TBondRepoRateType *repoRateType = NULL;

    double             repoRate;
    double             spotPrice;
    double             fwdPrice;
    double             vol;
    TDate              settleDate;

    long               optionType;
    double             tExp;
    double             pvFactor;
    TOptionProperties  optResult;
    double             optionPrice;

    GTO_IF_LOGGING(CrxBondPriceOptionCalcLogInputs (1e6,
                                                    deal,
                                                    distType,
                                                    today,
                                                    volCurve,
                                                    bondPrice,
                                                    repoCurve,
                                                    discCurve));

    REQUIRE (deal != NULL);
    REQUIRE (repoCurve != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (distType == CRX_DIST_TYPE_LOGNORMAL || distType == CRX_DIST_TYPE_NORMAL);
    REQUIRE (today > 0);
    REQUIRE (volCurve != NULL);
    REQUIRE (bondPrice != NULL);
    REQUIRE (bondPrice->cleanPrice > 0.0);
    REQUIRE (deal->exerciseDate >= today);
    REQUIRE (deal->paymentDate >= deal->exerciseDate);
    REQUIRE (deal->paymentDate < deal->bond->maturityDate);
    REQUIRE (repoCurve->spotSettleDate > 0);

    settleDate = repoCurve->spotSettleDate;

    if (volCurve->today != 0)
    {
        REQUIRE (volCurve->today == today);
    }

    if (bondPrice->settleDate != 0)
    {
        REQUIRE (bondPrice->settleDate == settleDate);
    }

    switch (deal->optionType)
    {
    case GTO_OPTION_PUT: 
    case GTO_OPTION_CALL:
        optionType = deal->optionType;
        break;
    default:
        GtoErrMsg ("%s: Invalid deal option type %ld\n",
                   routine, deal->optionType);
        goto done; /* failure */
    }

    /* we copy the bond to maintain the const status of deal, since
       it is a property of the bond library that it changes the
       bond by filling in missing details */

    bond            = GtoCopyTBond (deal->bond);

    /* construct all the other objects required by the ALIB repo
       calculator */
    spotTrade       = GtoNewTBondTradeData (settleDate);
    fwdTrade        = GtoNewTBondTradeData (deal->paymentDate);
    repoRateType    = GtoNewTBondRepoRateType (repoCurve->dcc,
                                               GTO_BOND_REPO_SIMPLE);

    /*
     * interpolation
     */
    if (GtoLinInterpLongPoint1 (volCurve->dates,
                                sizeof(TDate),
                                volCurve->numDates,
                                volCurve->vols,
                                sizeof(double),
                                deal->exerciseDate,
                                NULL,
                                &vol) != SUCCESS)
    {
        GtoErrMsg ("%s: Could not interpolate bond price volatility\n",
                   routine);
        goto done; /* failure */
    }
    REQUIRE (vol >= 0.0);

    if (GtoLinInterpLongPoint1 (repoCurve->dates,
                                sizeof(TDate),
                                repoCurve->numDates,
                                repoCurve->rates,
                                sizeof(double),
                                deal->paymentDate,
                                NULL,
                                &repoRate) != SUCCESS)
    {
        GtoErrMsg ("%s: Could not interpolate repo rate for bond\n", routine);
        goto done; /* failure */
    }

    spotPrice = bondPrice->cleanPrice;
    
    if (GtoBondRepoCalc (bond,
                         spotTrade,
                         fwdTrade,
                         repoRateType,
                         GTO_REPO_CALC_FWD_PRICE,
                         &repoRate,
                         &spotPrice,
                         &fwdPrice) != SUCCESS)
    {
        GtoErrMsg ("%s: Could not compute forward price for bond\n", routine);
        goto done; /* failure */
    }
    
    /*
     * Price the option. 
     * Vol applies to the forward price
     */

    if (GtoDayCountFraction(today, deal->exerciseDate, GTO_ACT_365F,
                            &tExp) != SUCCESS)
        goto done; /* failure */

    if (GtoDiscountDateForward(today, deal->paymentDate, discCurve,
                               GTO_FLAT_FORWARDS,
                               &pvFactor) != SUCCESS)
        goto done; /* failure */

    switch (distType)
    {
    case CRX_DIST_TYPE_NORMAL:
    {
        TOptionProperties *tmp = GtoNormalOption (optionType,
                                                  fwdPrice,
                                                  deal->strikePrice,
                                                  tExp,
                                                  0.0, /*yearsToPayment*/
                                                  vol * fwdPrice,
                                                  0.0, /*discountRate*/
                                                  GTO_OPTION_PRICE);
        if (tmp == NULL)
            goto done; /* failure */
        optResult = *tmp;
        FREE (tmp);
        break;
    }

    case CRX_DIST_TYPE_LOGNORMAL:
         if(GtoOptionsAnalytics2(optionType,
                                 fwdPrice,
                                 deal->strikePrice,
                                 tExp,
                                 0.0, /*yearsToPayment*/
                                 vol,
                                 0.0, /*discountRate*/
                                 0.0, /*dividendRate*/
                                 0.0, /*growthRate*/
                                 GTO_OPTION_PRICE, 
                                 &optResult) != SUCCESS)
         {
             goto done; /* failure */
         }
         break;
    default:
        GtoErrMsg ("%s: Invalid distribution type %ld\n",
                   routine, distType);
        goto done; /* failure */
    }

    /* following line is a strange necessity to avoid annoying -0.0 prices */
    if (IS_EQUAL(optResult.fPrice, -0.0)) 
        optResult.fPrice = 0.0;

    optionPrice = optResult.fPrice * pvFactor;

    /* all done */

    calc = CrxBondPriceOptionCalcMake (optionPrice,
                                       repoRate,
                                       fwdPrice,
                                       optResult.fPrice);

 done:

    GtoFreeTBondRepoRateType (repoRateType);
    GtoFreeTBondTradeData (fwdTrade);
    GtoFreeTBondTradeData (spotTrade);
    GtoFreeTBond (bond);

    if (calc == NULL)
        GtoErrMsgFailure (routine);

    return calc;
}


void CrxBondPriceOptionCalcLogInputs (
    double                  notional,
    CrxTBondPriceOption*    deal,             
    long                    distType,         
    TDate                   today,            
    CrxTBondPriceVolCurve*  volCurve,
    CrxTBondPrice*          bondPrice,
    CrxTBondRepoCurve*      repoCurve,        
    TCurve*                 discCurve)
{
    FILE *fp = fopen("term.prn", "w");
    static char FREQ[] = "?AS?Q???????M";

    TBond  *bond = deal->bond;

    fprintf (fp, "\nFILE:crxbondprcopt_f.dat\n\n");
    fprintf (fp, "# option type (Call/Put)\n%c\n",  
             (char)deal->optionType);
    fprintf (fp, "# Notional\n%g\n",
             notional);
    fprintf (fp, "# expiration date\n%ld\n",
             YYYYMMDD(deal->exerciseDate));
    fprintf (fp, "# pay date\n%ld\n",
             YYYYMMDD(deal->paymentDate));
    fprintf (fp, "# Bond maturity date\n%ld\n",
             YYYYMMDD(bond->maturityDate));
    fprintf (fp, "# Bond issue date\n%ld\n",
             YYYYMMDD(bond->datedDate));
    fprintf (fp, "# Bond first coupon date\n%ld\n",
             YYYYMMDD(bond->firstCouponDate));
    fprintf (fp, "# Bond coupon interval (Q/S/A/M)\n%c\n",
             FREQ[bond->couponFrequency]);
    fprintf (fp, "# Bond coupon rate (in %%)\n%g\n",
             bond->couponRate * 1e2);
    fprintf (fp, "# Bond coupon day count conv\n%s\n",
             GtoFormatDayCountConv(bond->dayCountConv));
    fprintf (fp, "# strike price (par=100)\n%g\n",
             1e2 * deal->strikePrice);
    fprintf (fp, "# \"L\"ognormal, \"N\"ormal\n%c\n",
             (char)distType);

    fprintf (fp, "\nFILE:today.dat\n\n");
    fprintf (fp, "%ld\n", 
             YYYYMMDD(today));

    fprintf (fp, "\nFILE:bondprice_0.dat\n\n");
    CrxBondPriceDRWWrite (fp, bondPrice);
    
    fprintf (fp, "\nFILE:bondpricevolcurve_0.dat\n\n");
    CrxBondPriceVolCurveDRWWrite (fp, volCurve);
    
    fprintf (fp, "\nFILE:zero.dat\n\n");
    CrxWriteTCurve (fp, discCurve);

    fprintf (fp, "\nFILE:bondrepo_0.dat\n\n");
    CrxBondRepoCurveDRWWrite (fp, repoCurve);

    if (fp != NULL)
        fclose (fp);
}


