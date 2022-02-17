/*
*****************************************************************************
** Calculation implementation for bond spread option pricer.
**
** $Header$
*****************************************************************************
*/

#include "bondsprdopt.h"

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

#undef BASIS_POINT
#define BASIS_POINT 1e-4

/*f
***************************************************************************
** Main calculator
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalc(
CrxTBondSpreadOption*   deal,             /* (I) */
long                    distType,         /* (I) */
TDate                   today,            /* (I) */
CrxTBondSpreadVolCurve* volCurve,         /* (I) */
CrxTBondPrice*          bondPrice,        /* (I) */
CrxTBondPrice*          refBondPrice,     /* (I) */
CrxTBondRepoCurve*      repoCurve,        /* (I) */
CrxTBondRepoCurve*      refRepoCurve,     /* (I) */
TCurve*                 discCurve         /* (I) */
)
{
    static char routine[] = "CrxBondSpreadOptionCalc";

    CrxTBondSpreadOptionCalc *calc = NULL;

    TBond             *bond = NULL;
    TBond             *refBond = NULL;
    TBondTradeData    *spotTrade = NULL;
    TBondTradeData    *refSpotTrade = NULL;
    TBondTradeData    *fwdTrade = NULL;
    TBondTradeData    *refFwdTrade = NULL;
    TBondRepoRateType *repoRateType = NULL;
    TBondRepoRateType *refRepoRateType = NULL;

    TBond*             bonds[2];
    TBondTradeData*    spotTrades[2];
    TBondTradeData*    fwdTrades[2];
    TBondRepoRateType* repoRateTypes[2];
    CrxTBondRepoCurve*  repoCurves[2];
    double             repoRates[2];
    double             spotPrices[2];
    double             fwdPrices[2];
    double             fwdYields[2];

    double             vol;

    double             fwdPrice;
    double             refFwdPrice;
    double             fwdYield;
    double             refFwdYield;
    double             fwdSpread;
    double             yieldStrike;
    double             repoRate;
    double             refRepoRate;

    TBondSensOutputs   bondSens;

    long               yldOptionType;
    double             tExp;
    double             pvFactor;
    double             dpdy;
    TOptionProperties  optResult;
    double             optionPrice;

    char* names[] = {"Bond", "ReferenceBond"};

    int i;

    GTO_IF_LOGGING(CrxBondSpreadOptionCalcLogInputs (1e6,
                                                     deal,
                                                     distType,
                                                     today,
                                                     volCurve,
                                                     bondPrice,
                                                     refBondPrice,
                                                     repoCurve,
                                                     refRepoCurve,
                                                     discCurve));

    REQUIRE (deal != NULL);
    REQUIRE (repoCurve != NULL);
    REQUIRE (refRepoCurve != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (distType == CRX_DIST_TYPE_LOGNORMAL || distType == CRX_DIST_TYPE_NORMAL);
    REQUIRE (today > 0);
    REQUIRE (volCurve != NULL);
    REQUIRE (bondPrice != NULL);
    REQUIRE (bondPrice->cleanPrice > 0.0);
    REQUIRE (refBondPrice != NULL);
    REQUIRE (refBondPrice->cleanPrice > 0.0);
    REQUIRE (deal->exerciseDate >= today);
    REQUIRE (deal->paymentDate >= deal->exerciseDate);
    REQUIRE (deal->paymentDate < deal->bond->maturityDate);
    REQUIRE (deal->paymentDate < deal->refBond->maturityDate);

    if (volCurve->today != 0)
    {
        REQUIRE (volCurve->today == today);
    }

    if (bondPrice->settleDate != 0)
    {
        REQUIRE (bondPrice->settleDate == repoCurve->spotSettleDate);
    }

    if (refBondPrice->settleDate != 0)
    {
        REQUIRE (refBondPrice->settleDate == refRepoCurve->spotSettleDate);
    }

    switch (deal->optionType)
    {
    case GTO_OPTION_PUT: 
        yldOptionType = GTO_OPTION_CALL;
        break;
    case GTO_OPTION_CALL:
        yldOptionType = GTO_OPTION_PUT;
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
    refBond         = GtoCopyTBond (deal->refBond);

    /* construct all the other objects required by the ALIB repo
       calculator */
    spotTrade       = GtoNewTBondTradeData (repoCurve->spotSettleDate);
    refSpotTrade    = GtoNewTBondTradeData (refRepoCurve->spotSettleDate);
    fwdTrade        = GtoNewTBondTradeData (deal->paymentDate);
    refFwdTrade     = GtoNewTBondTradeData (deal->paymentDate);
    repoRateType    = GtoNewTBondRepoRateType (repoCurve->dcc,
                                               GTO_BOND_REPO_SIMPLE);
    refRepoRateType = GtoNewTBondRepoRateType (refRepoCurve->dcc,
                                               GTO_BOND_REPO_SIMPLE);

    /* set up little arrays of size 2 since we are going to do the
       same calculations for both the bond and the reference bond */

    bonds[0]         = bond;
    spotTrades[0]    = spotTrade;
    fwdTrades[0]     = fwdTrade;
    repoCurves[0]    = repoCurve;
    repoRateTypes[0] = repoRateType;
    spotPrices[0]    = bondPrice->cleanPrice;

    bonds[1]         = refBond;
    spotTrades[1]    = refSpotTrade;
    fwdTrades[1]     = refFwdTrade;
    repoCurves[1]    = refRepoCurve;
    repoRateTypes[1] = refRepoRateType;
    spotPrices[1]    = refBondPrice->cleanPrice;

    /*
     * vol interpolation
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

    /*
     * repo interpolation and forward price calculation
     */
    for (i = 0; i < 2; ++i)
    {
        if (GtoLinInterpLongPoint1 (repoCurves[i]->dates,
                                    sizeof(TDate),
                                    repoCurves[i]->numDates,
                                    repoCurves[i]->rates,
                                    sizeof(double),
                                    deal->paymentDate,
                                    NULL,
                                    &repoRates[i]) != SUCCESS)
        {
            GtoErrMsg ("%s: Could not interpolate repo rate for %s\n",
                       routine, names[i]);
            goto done; /* failure */
        }

        if (GtoBondRepoCalc (bonds[i],
                             spotTrades[i],
                             fwdTrades[i],
                             repoRateTypes[i],
                             GTO_REPO_CALC_FWD_PRICE,
                             &repoRates[i],
                             &spotPrices[i],
                             &fwdPrices[i]) != SUCCESS)
        {
            GtoErrMsg ("%s: Could not compute forward price for %s\n",
                       routine, names[i]);
            goto done; /* failure */
        }
        
        if (GtoBondYldToMaturity (bonds[i],
                                  fwdTrades[i],
                                  fwdPrices[i],
                                  FALSE,
                                  GTO_CALC_COMPOUND_ALWAYS,
                                  &fwdYields[i]) != SUCCESS)
        {
            GtoErrMsg ("%s: Could not compute forward compound yield for %s\n",
                       routine, names[i]);
            goto done; /* failure */
        }
    }

    repoRate    = repoRates[0];
    refRepoRate = repoRates[1];
    fwdPrice    = fwdPrices[0];
    refFwdPrice = fwdPrices[1];
    fwdYield    = fwdYields[0];
    refFwdYield = fwdYields[1];
    fwdSpread   = fwdYield - refFwdYield;
    yieldStrike = refFwdYield + deal->strikeSpread;

    if (GtoBondSensitivity (bond,
                            fwdTrade,
                            yieldStrike,
                            GTO_CALC_COMPOUND_ALWAYS,
                            &bondSens) != SUCCESS)
    {
        GtoErrMsg ("%s: Could not compute strike price for Bond\n", routine);
        goto done; /* failure */
    }

    if (fabs(yieldStrike - fwdYield) < BASIS_POINT)
    {
        dpdy = bondSens.pvbp * 1e4;
    }
    else
    {
        dpdy = (fwdPrice - bondSens.cleanPrice) / (yieldStrike - fwdYield);
    }

    /* Price option. 
     * Vol applies to the fwd spread, which is the difference of the two ytm's
     */

    if (GtoDayCountFraction(today, deal->exerciseDate, GTO_ACT_365F,
                            &tExp) != SUCCESS)
        goto done; /* failure */

    if (GtoDiscountDateForward(today, deal->paymentDate, discCurve, 
                               GTO_FLAT_FORWARDS, &pvFactor) != SUCCESS)
        goto done; /* failure */

    switch (distType)
    {
    case CRX_DIST_TYPE_NORMAL:
    {
        TOptionProperties *tmp = GtoNormalOption (yldOptionType,
                                                  fwdSpread,
                                                  deal->strikeSpread,
                                                  tExp,
                                                  0.0, /*yearsToPayment*/
                                                  vol * fabs(fwdSpread),
                                                  0.0, /*discountRate*/
                                                  GTO_OPTION_PRICE);
        if (tmp == NULL)
            goto done; /* failure */
        optResult = *tmp;
        FREE (tmp);
        break;
    }

    case CRX_DIST_TYPE_LOGNORMAL:
         if (fwdSpread <= 0.0)
         {
             GtoErrMsg("%s: forward spread is zero or negative", routine);
             goto done; /* failure */
         }

         if(GtoOptionsAnalytics2(yldOptionType,
                                 fwdSpread,
                                 deal->strikeSpread,
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

    optionPrice = optResult.fPrice * pvFactor * dpdy;

    /* all done */

    calc = CrxBondSpreadOptionCalcMake (optionPrice,
                                        repoRate,
                                        refRepoRate,
                                        fwdPrice,
                                        refFwdPrice,
                                        fwdSpread,
                                        dpdy * pvFactor,
                                        optResult.fPrice);

 done:

    GtoFreeTBondRepoRateType (refRepoRateType);
    GtoFreeTBondRepoRateType (repoRateType);
    GtoFreeTBondTradeData (refFwdTrade);
    GtoFreeTBondTradeData (fwdTrade);
    GtoFreeTBondTradeData (refSpotTrade);
    GtoFreeTBondTradeData (spotTrade);
    GtoFreeTBond (refBond);
    GtoFreeTBond (bond);

    if (calc == NULL)
        GtoErrMsgFailure (routine);

    return calc;
}

#define YYYYMMDD(x) x>0?CrxTDate2DrDate(x):0

void CrxBondSpreadOptionCalcLogInputs (
    double                   notional,
    CrxTBondSpreadOption*    deal,             
    long                     distType,         
    TDate                    today,            
    CrxTBondSpreadVolCurve*  volCurve,
    CrxTBondPrice*           bondPrice,
    CrxTBondPrice*           refBondPrice,
    CrxTBondRepoCurve*       repoCurve,        
    CrxTBondRepoCurve*       refRepoCurve,     
    TCurve*                  discCurve)
{
    FILE *fp = fopen("term.prn", "w");
    static char FREQ[] = "?AS?Q???????M";

    TBond               *bond      = deal->bond;
    TBond               *refBond   = deal->refBond;

    fprintf (fp, "\nFILE:crxbondsprdopt_f.dat\n\n");
    fprintf (fp, "# option type (Call/Put)\n%c\n",  
             (char)deal->optionType);
    fprintf (fp, "# Notional\n%g\n",
             notional);
    fprintf (fp, "# expiration date\n%ld\n",
             YYYYMMDD(deal->exerciseDate));
    fprintf (fp, "# pay date\n%ld\n",
             YYYYMMDD(deal->paymentDate));
    fprintf (fp, "# Reference maturity date\n%ld\n",
             YYYYMMDD(refBond->maturityDate));
    fprintf (fp, "# Reference issue date\n%ld\n",
             YYYYMMDD(refBond->datedDate));
    fprintf (fp, "# Reference first coupon date\n%ld\n",
             YYYYMMDD(refBond->firstCouponDate));
    fprintf (fp, "# Reference coupon interval (Q/S/A/M)\n%c\n",
             FREQ[refBond->couponFrequency]);
    fprintf (fp, "# Reference coupon rate (in %%)\n%g\n",
             refBond->couponRate * 1e2);
    fprintf (fp, "# Reference coupon day count conv\n%s\n",
             GtoFormatDayCountConv(refBond->dayCountConv));
    fprintf (fp, "# Underlying maturity date\n%ld\n",
             YYYYMMDD(bond->maturityDate));
    fprintf (fp, "# Underlying issue date\n%ld\n",
             YYYYMMDD(bond->datedDate));
    fprintf (fp, "# Underlying first coupon date\n%ld\n",
             YYYYMMDD(bond->firstCouponDate));
    fprintf (fp, "# Underlying coupon interval (Q/S/A/M)\n%c\n",
             FREQ[bond->couponFrequency]);
    fprintf (fp, "# Underlying coupon rate (in %%)\n%g\n",
             bond->couponRate * 1e2);
    fprintf (fp, "# Underlying coupon day count conv\n%s\n",
             GtoFormatDayCountConv(bond->dayCountConv));
    fprintf (fp, "# strike spread (in bps)\n%g\n",
             1e4 * deal->strikeSpread);
    fprintf (fp, "# \"L\"ognormal, \"N\"ormal\n%c\n",
             (char)distType);

    fprintf (fp, "\nFILE:today.dat\n\n");
    fprintf (fp, "%ld\n", 
             YYYYMMDD(today));

    fprintf (fp, "\nFILE:bondprice_0.dat\n\n");
    CrxBondPriceDRWWrite (fp, bondPrice);
    
    fprintf (fp, "\nFILE:bondprice_1.dat\n\n");
    CrxBondPriceDRWWrite (fp, refBondPrice);
    
    fprintf (fp, "\nFILE:bondspreadvolcurve_0.dat\n\n");
    CrxBondSpreadVolCurveDRWWrite (fp, volCurve);
    
    fprintf (fp, "\nFILE:zero.dat\n\n");
    CrxWriteTCurve (fp, discCurve);

    fprintf (fp, "\nFILE:bondrepo_0.dat\n\n");
    fprintf (fp, "# Repo curve for underlying bond\n");
    CrxBondRepoCurveDRWWrite (fp, repoCurve);

    fprintf (fp, "\nFILE:bondrepo_1.dat\n\n");
    fprintf (fp, "# Repo curve for reference bond\n");
    CrxBondRepoCurveDRWWrite (fp, refRepoCurve);

    if (fp != NULL)
        fclose (fp);
}

