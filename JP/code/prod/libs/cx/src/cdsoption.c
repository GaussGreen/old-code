/*
***************************************************************************
** FILENAME: cdsoption.c
**
** Calculation routine for single european option on CDS with adjusted
** forward and multi-Q
***************************************************************************
*/

#include "cdsoption.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <alib/dtivlo.h>
#include <common/include/drutils.h>
#include <crxmultiq/include/crmultiq.h>
#include <crxflow/include/crcrv.h>

#include <alib/lintrp.h>
#include <alib/rtbrent.h>

#include <cxutils/include/alibconv.h>
#include "cds.h"
#include "cdsbootstrap.h"
#include "recovery.h"

#include <cxutils/include/cxmacros.h>

#define INTERP_METHOD GTO_FLAT_FORWARDS

typedef struct _IMPLIED_VOL_PARAMS
{
    CrxTCdsOption*  deal;
    long            distType;
    TDate           today;
    double          price;
    CrxTQDist*      qdist;
    TCurve*         discCurve;
    CxTCreditCurve* sprdCurve;
    double          recoveryRate;
} IMPLIED_VOL_PARAMS;

static int ImpliedVolSolver (
    double              vol,
    IMPLIED_VOL_PARAMS *params,
    double             *priceDiff);

static int FwdDiscountFactor(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    double           *discountFactor);

static int FwdDefaultProbability(
    CrxTCdsOption    *deal,
    TDate             today,
    CxTCreditCurve   *sprdCurve,
    double           *defaultProb);

static int FwdSurvivalProbability(
    CrxTCdsOption    *deal,
    TDate             today,
    CxTCreditCurve   *sprdCurve,
    double           *survivalProb);

static int AdjustedFwdParSpread(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    CxTCreditCurve   *sprdCurve,
    double            recoveryRate,
    double           *fwdSpread,
    double           *fwdSpreadUnadj,
    double           *annuity
);

static int AdjustedFwdPrice(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    CxTCreditCurve   *sprdCurve,
    double            recoveryRate,
    double           *fwdPrice
);

static CrxTCdsOptionCalc* CdsOptionCalcSpreadStrike(
    CrxTCdsOption*  deal,             /* (I) */
    long            distType,         /* (I) */
    TDate           today,            /* (I) */
    double          vol,              /* (I) */
    CrxTQDist*      qdist,            /* (I) */
    TCurve*         discCurve,        /* (I) */
    CxTCreditCurve* sprdCurve,        /* (I) */
    double          recoveryRate      /* (I) */
    );

static CrxTCdsOptionCalc* CdsOptionCalcPriceStrike(
    CrxTCdsOption*  deal,             /* (I) */
    long            distType,         /* (I) */
    TDate           today,            /* (I) */
    double          vol,              /* (I) */
    CrxTQDist*      qdist,            /* (I) */
    TCurve*         discCurve,        /* (I) */
    CxTCreditCurve* sprdCurve,        /* (I) */
    double          recoveryRate      /* (I) */
    );


/*f
***************************************************************************
** Main calculator - price or spread strike
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalc(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          vol,              /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
)
{
    static char routine[] = "CrxCdsOptionCalc";

    CrxTCdsOptionCalc* calc = NULL;

    REQUIRE (deal != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (sprdCurve != NULL);

    REQUIRE (deal->exerciseDate >= today);
    REQUIRE (recoveryRate >= 0.0);
    REQUIRE (recoveryRate < 1.0);

    switch (deal->strikeType)
    {
    case CRX_CDS_OPTION_STRIKE_TYPE_SPREAD:
        calc = CdsOptionCalcSpreadStrike (deal,
                                          distType,
                                          today,
                                          vol,
                                          qdist,
                                          discCurve,
                                          sprdCurve,
                                          recoveryRate);
        break;
    case CRX_CDS_OPTION_STRIKE_TYPE_PRICE:
        calc = CdsOptionCalcPriceStrike (deal,
                                         distType,
                                         today,
                                         vol,
                                         qdist,
                                         discCurve,
                                         sprdCurve,
                                         recoveryRate);
        break;
    default:
        PROGRAM_BUG();
        goto done; /* failure */
    }
    
 done:

    if (calc == NULL)
        GtoErrMsgFailure(routine);

    return calc;
}


/*f
***************************************************************************
** Main calculator - spread strike only
***************************************************************************
*/
static CrxTCdsOptionCalc* CdsOptionCalcSpreadStrike(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          vol,              /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
)
{
    static char routine[] = "CdsOptionCalcSpreadStrike";
    int         status    = FAILURE;

    CrxTCdsOptionCalc *calc = NULL;

    MQDATA mq;
    double price;
    double t2exp;
    int q3OptionType;
    int i;
    int nQs;
    double *myQs = NULL;
    double *myDs = NULL;

    double fwdSpread;
    double fwdSpreadUnadj;
    double annuity;

    REQUIRE (deal->strikeType == CRX_CDS_OPTION_STRIKE_TYPE_SPREAD);

    switch (deal->optionType)
    {
    case GTO_OPTION_CALL: 
        q3OptionType = Q3_PUT;
        break;
    case GTO_OPTION_PUT: 
        q3OptionType = Q3_CALL;
        break;
    case CRX_OPTION_TYPE_STRADDLE: 
        /* put and straddle have the same q3OptionType, we
           will use this to determine if we do the No Ko adj*/
        q3OptionType = Q3_CALL; 
        break; 
    default:
        GtoErrMsg ("%s: Bad option type %ld\n", routine, deal->optionType);
        goto done;
    }

    switch (distType)
    {
    case CRX_DIST_TYPE_Q: 
        REQUIRE (qdist != NULL);
        nQs  = qdist->nQs;
        myQs = NEW_ARRAY(double, nQs);
        myDs = NEW_ARRAY(double, nQs);
        for (i = 0; i < nQs-1; ++i)
        {
            myQs[i] = 1.0 - qdist->Qs[i];
            myDs[i] = qdist->Ds[i];
        }
        myQs[nQs-1] = 1.0 - qdist->Qs[nQs-1];
        myDs[nQs-1] = 1.0; /* this value should be ignored */
        break;
    case CRX_DIST_TYPE_NORMAL:
    case CRX_DIST_TYPE_LOGNORMAL:
    {
        double q;
        
        q = (distType == CRX_DIST_TYPE_LOGNORMAL) ? 1.0 : 0.0;

        nQs = 6;
        myQs = NEW_ARRAY(double, nQs);
        myDs = NEW_ARRAY(double, nQs);
        for (i = 0; i < 6; ++i)
        {
            myQs[i] = q;
            myDs[i] = (double) (i+1) / (double)nQs;
        }
        break;
    }
    default: 
        GtoErrMsg("%s: Wrong distribution type %ld\n", routine, distType);
        goto done;
    }

    if (AdjustedFwdParSpread(deal,
                             today,
                             discCurve,
                             sprdCurve,
                             recoveryRate,
                             &fwdSpread,
                             &fwdSpreadUnadj,
                             &annuity) != SUCCESS)
    {
        GtoErrMsg("%s: Failed to compute adjusted fwd par spread\n", routine);
        goto done;
    };

    fwdSpread += deal->lossSoFar/annuity;  /* take loss history into account */

    t2exp = (double)(deal->exerciseDate - today) / 365.0;

    if ( Q3MQInitCR(
            &mq, 
            fwdSpread,
            vol, 
            q3OptionType, 
            t2exp, 
            nQs, 
            myQs, 
            nQs-1, 
            myDs) == FAILURE)
    {
        GtoErrMsg("%s: failed to initialise MQ structure\n", routine);
        goto done;
    }

    /* The following will calculate the price of a put if it is a straddle*/
    if (Q3MQPricerCR(
            &mq, 
            q3OptionType, 
            deal->strike, 
            &price) == FAILURE)
    {
        GtoErrMsg("%s: failed to price option under MQ dist.\n", routine);
        goto done;
    }

    /* If it is a straddle, we add the price of a call*/
    if (deal->optionType == CRX_OPTION_TYPE_STRADDLE)
    {
        double priceCall;
        if (Q3MQPricerCR(
            &mq, 
            Q3_PUT, 
            deal->strike, 
            &priceCall) == FAILURE)
        {
            GtoErrMsg("%s: failed to price option under MQ dist.\n", routine);
            goto done;
        }
        else price +=priceCall;
    }

    price = price * annuity;

    /* 
    ** We adjust for early default for single name put options which do not
    ** knock out on default.
    **
    ** In this case when we have the right to buy protection (a put), there
    ** is the chance that the name has defaulted already, in which case we have
    ** the right to get the protection payout at option expiry.
    **
    ** Note that index options behave differently - partial defaults are
    ** taken account of by the adjusted forward.
    */
    
    if ((q3OptionType == Q3_CALL) &&   /* put or straddle, no KO on default */
        (! deal->koOnDefault ) &&
        (! deal->isIndex))
    {    
        double discFactor;
        double defaultProb;

        if (FwdDiscountFactor(deal, today, discCurve, &discFactor) != SUCCESS)
            goto done; /* failure */
        
        if (FwdDefaultProbability(deal, today, sprdCurve,
                                  &defaultProb) != SUCCESS)
            goto done; /* failure */

        price += (1.0-recoveryRate) * discFactor * defaultProb;
    }

    calc = CrxCdsOptionCalcMake (price, 
                                 vol, 
                                 fwdSpread,
                                 fwdSpreadUnadj,
                                 annuity,
                                 t2exp,
                                 0.0); /* fwdPrice */
    if (calc == NULL) goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
        CrxCdsOptionCalcFree (calc);
        calc = NULL;
    }

    FREE (myQs);
    FREE (myDs);
    return calc;
}


/*f
***************************************************************************
** Main calculator - price strike only
***************************************************************************
*/
static CrxTCdsOptionCalc* CdsOptionCalcPriceStrike(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          vol,              /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
)
{
    static char routine[] = "CdsOptionCalcPriceStrike";
    int         status    = FAILURE;

    CrxTCdsOptionCalc *calc = NULL;

    MQDATA mq;
    double optionPrice;
    double discFactor;
    double survivalProb;
    double t2exp;
    int q3OptionType;
    int i;
    int nQs;
    double *myQs = NULL;
    double *myDs = NULL;

    double fwdPrice;
    
    REQUIRE (deal->strikeType == CRX_CDS_OPTION_STRIKE_TYPE_PRICE);
    REQUIRE (deal->coupon > 0.0);

    switch (deal->optionType)
    {
    case GTO_OPTION_CALL: q3OptionType = Q3_CALL; break;
    case GTO_OPTION_PUT:  q3OptionType = Q3_PUT;  break;
    case CRX_OPTION_TYPE_STRADDLE:  q3OptionType = Q3_PUT;  break;
    default:
        GtoErrMsg ("%s: Bad option type %ld\n", routine, deal->optionType);
        goto done;
    }

    switch (distType)
    {
    case CRX_DIST_TYPE_Q: 
        REQUIRE (qdist != NULL);
        nQs  = qdist->nQs;
        myQs = NEW_ARRAY(double, nQs);
        myDs = NEW_ARRAY(double, nQs);
        for (i = 0; i < nQs-1; ++i)
        {
            myQs[i] = 1.0 - qdist->Qs[i];
            myDs[i] = qdist->Ds[i];
        }
        myQs[nQs-1] = 1.0 - qdist->Qs[nQs-1];
        myDs[nQs-1] = 1.0; /* this value should be ignored */
        break;
    case CRX_DIST_TYPE_NORMAL:
    case CRX_DIST_TYPE_LOGNORMAL:
    {
        double q;
        
        q = (distType == CRX_DIST_TYPE_LOGNORMAL) ? 1.0 : 0.0;

        nQs = 6;
        myQs = NEW_ARRAY(double, nQs);
        myDs = NEW_ARRAY(double, nQs);
        for (i = 0; i < 6; ++i)
        {
            myQs[i] = q;
            myDs[i] = (double) (i+1) / (double)nQs;
        }
        break;
    }
    default: 
        GtoErrMsg("%s: Wrong distribution type %ld\n", routine, distType);
        goto done;
    }

    if (AdjustedFwdPrice(deal,
                         today,
                         discCurve,
                         sprdCurve,
                         recoveryRate,
                         &fwdPrice) != SUCCESS)
    {
        GtoErrMsg("%s: Failed to compute adjusted fwd price\n", routine);
        goto done;
    }

    fwdPrice += deal->lossSoFar;  /* take loss history into account */

    t2exp = (double)(deal->exerciseDate - today) / 365.0;

    if ( Q3MQInitCR(
            &mq, 
            fwdPrice,
            vol, 
            q3OptionType, 
            t2exp, 
            nQs, 
            myQs, 
            nQs-1, 
            myDs) == FAILURE)
    {
        GtoErrMsg("%s: failed to initialise MQ structure\n", routine);
        goto done;
    }


    if (Q3MQPricerCR(
            &mq, 
            q3OptionType, 
            deal->strike, 
            &optionPrice) == FAILURE)
    {
        GtoErrMsg("%s: failed to price option under MQ dist.\n", routine);
        goto done;
    }
    if (deal->optionType == CRX_OPTION_TYPE_STRADDLE)
    {
        double optionPrice_Call;
        if (Q3MQPricerCR(
            &mq, 
            Q3_CALL, 
            deal->strike, 
            &optionPrice_Call) == FAILURE)
        {
            GtoErrMsg("%s: failed to price option under MQ dist.\n", routine);
            goto done;
        }
        else optionPrice += optionPrice_Call;
    }
    
    /* Discount the price back to today */
    if (FwdDiscountFactor (deal, today, discCurve, &discFactor) != SUCCESS)
        goto done; /* failure */

    if (deal->isIndex)
    {
        /* forward price incorporates default probabilities of single names */
        /* so don't need to include the survival probability here */
        optionPrice = optionPrice * discFactor;
    }
    else
    {
        if (FwdSurvivalProbability (deal, today, sprdCurve, 
                                    &survivalProb) != SUCCESS)
            goto done; /* failure */
        
        optionPrice = optionPrice * discFactor * survivalProb;

        /* 
        ** We adjust for early default for single name put options which do
        ** not knock out on default.
        **
        ** In this case when we have the right to buy protection (a put),
        ** there is the chance that the name has defaulted already, in which
        ** case we have the right to get the protection payout at option
        ** expiry.
        */
        if (!deal->koOnDefault)
        {
            double adj = 0.0;
            double defaultProb = 1.0 - survivalProb;
            switch (deal->optionType)
            {
            case GTO_OPTION_CALL:
                adj = MAX(recoveryRate - deal->strike, 0.0);
                break;
            case GTO_OPTION_PUT:
                adj = MAX(deal->strike - recoveryRate, 0.0);
                break;
            case CRX_OPTION_TYPE_STRADDLE:
                adj = fabs(deal->strike - recoveryRate);
                break;
            }
            optionPrice += adj * defaultProb * discFactor;
        }
    }

    calc = CrxCdsOptionCalcMake (optionPrice, 
                                 vol, 
                                 0.0, /*fwdSpread*/
                                 0.0, /*fwdSpreadUnadj*/
                                 0.0, /*annuity*/
                                 t2exp,
                                 fwdPrice);
    if (calc == NULL) goto done;

    status = SUCCESS;

done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
        CrxCdsOptionCalcFree (calc);
        calc = NULL;
    }

    FREE (myQs);
    FREE (myDs);
    return calc;
}



/*f
***************************************************************************
** Implied volatility calculator
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionVolCalc(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          price,            /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
)
{
    static char routine[] = "CdsOptionVolCalc";

    CrxTCdsOptionCalc *calc = NULL;
    IMPLIED_VOL_PARAMS params;
    double vol;

    /* parameters for root finder */
    double boundLo       = 0.0;
    double boundHi       = 100.0;
    int    numIterations = 100;
    double guess         = 0.5;
    double initialXStep  = 0.01;
    double initialFDeriv = 0.0;
    double xacc          = 1e-10;
    double facc          = 1e-10;

    params.deal         = deal;
    params.distType     = distType;
    params.today        = today;
    params.price        = price;
    params.qdist        = qdist;
    params.discCurve    = discCurve;
    params.sprdCurve    = sprdCurve;
    params.recoveryRate = recoveryRate;

    if (GtoRootFindBrent ((TObjectFunc)ImpliedVolSolver,
                          &params,
                          boundLo,
                          boundHi,
                          numIterations,
                          guess,
                          initialXStep,
                          initialFDeriv,
                          xacc,
                          facc,
                          &vol) != SUCCESS)
        goto done; /* failure */

    calc = CrxCdsOptionCalc (deal,
                             distType,
                             today,
                             vol,
                             qdist,
                             discCurve,
                             sprdCurve,
                             recoveryRate);

 done:
    
    if (calc == NULL)
        GtoErrMsgFailure (routine);

    return calc;
}




static int AdjustedFwdParSpread(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    CxTCreditCurve   *sprdCurve,
    double            recoveryRate,
    double           *fwdSpread,
    double           *fwdSpreadUnadj,
    double           *annuity
)
{
    char routine[] = "AdjustedFwdParSpread";
    int  status = FAILURE;

    double rld, ryd, fwd, fwdann, adj;
    double fwdAnnK;

    CxTCreditCurve *strikeSpreadCurve = NULL;

    CxTRecoveryCurve *recoveryCurve = NULL;
    
    recoveryCurve = CxRecoveryCurveMakeFromRecoveryRate(recoveryRate);
    if (recoveryCurve == NULL)
        goto done; /* failure */

    /* TBD - we should doubtless discount these to today instead */

    GtoDiscountDate(deal->exerciseDate, discCurve, INTERP_METHOD, &rld);
    
    GtoDiscountDate(deal->exerciseDate, sprdCurve->tc, INTERP_METHOD, &ryd);

    if (CxCdsFeeLegPV (today,
                       today,
                       deal->exerciseDate,
                       deal->maturityDate,
                       deal->payAccOnDefault,
                       deal->feeInterval,
                       CX_SHORT_FRONT_STUB,
                       1.0, /* notional */
                       1.0, /* couponRate */
                       deal->dcc,
                       CX_BAD_DAY_NONE,
                       NULL,
                       discCurve,
                       sprdCurve,
                       FALSE,
                       FALSE,
                       &fwdann) != SUCCESS)
        goto done; /* failure */

    if (CxCdsParSpread (today,
                        today,
                        deal->exerciseDate,
                        deal->maturityDate,
                        0,   /* delay */
                        0.0, /* price */
                        deal->payAccOnDefault,
                        deal->feeInterval,
                        CX_SHORT_FRONT_STUB,
                        deal->dcc,
                        CX_BAD_DAY_NONE,
                        NULL,
                        discCurve,
                        sprdCurve,
                        recoveryCurve,
                        FALSE,
                        FALSE,
                        &fwd) != SUCCESS)
        goto done; /* failure */

    if (deal->isIndex)
    {
        adj = (1.0 - recoveryRate) * rld * (1. - ryd);
        switch (deal->optionPayoff)
        {
        case CRX_CDS_OPTION_PAYOFF_MARKET:
            break; /* do nothing */
        case CRX_CDS_OPTION_PAYOFF_STRIKE:
            /*
             * We need to calculate a credit curve as at the exercise date
             * as if the market spread was equal to the strike spread in
             * order to get the upfront charge required on exercising the
             * option. This then feeds into a second adjustment in the
             * forward spread (see Mehdi's note on the Adjusted(2)FwdSpread)
             */
            strikeSpreadCurve = CxCdsBootstrap(deal->exerciseDate,
                                               discCurve,
                                               deal->exerciseDate,
                                               deal->exerciseDate,
                                               1,
                                               &deal->maturityDate,
                                               &deal->strike,
                                               NULL,
                                               NULL,
                                               recoveryCurve,
                                               deal->payAccOnDefault,
                                               deal->feeInterval,
                                               deal->dcc,
                                               CX_SHORT_FRONT_STUB,
                                               CX_CURVE_TYPE_FLOW,
                                               NULL,
                                               NULL,
                                               FALSE, /* CMLib compatible */
                                               FALSE,
                                               0,
                                               CX_BAD_DAY_NONE,
                                               NULL);
            if (strikeSpreadCurve == NULL)
                goto done;
            if (CxCdsFeeLegPV (deal->exerciseDate,
                               deal->exerciseDate,
                               deal->exerciseDate,
                               deal->maturityDate,
                               deal->payAccOnDefault,
                               deal->feeInterval,
                               CX_SHORT_FRONT_STUB,
                               1.0,
                               1.0,
                               deal->dcc,
                               CX_BAD_DAY_NONE,
                               NULL,
                               discCurve,
                               strikeSpreadCurve,
                               FALSE, /* CMLib compatible */
                               FALSE,
                               &fwdAnnK) != SUCCESS)
                goto done;

            adj += (deal->strike - deal->coupon) * (fwdann - rld * fwdAnnK);
            break;
        default:
            GtoErrMsg ("%s: Unexpected value (%ld) for option payoff\n",
                       routine, (long)(deal->optionPayoff));
            goto done;
        }
    }
    else
    {
        adj = 0.0;
    }

    *fwdSpread      = fwd + adj/fwdann;
    *fwdSpreadUnadj = fwd;
    *annuity        = fwdann;

    status = SUCCESS;

 done:

    CxCreditCurveFree (strikeSpreadCurve);
    CxRecoveryCurveFree (recoveryCurve);

    if (status != SUCCESS)
        GtoErrMsgFailure(routine);

    return status;

}

static int AdjustedFwdPrice(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    CxTCreditCurve   *sprdCurve,
    double            recoveryRate,
    double           *fwdPrice
)
{
    char routine[] = "AdjustedFwdPrice";
    int  status = FAILURE;

    double upfrontCharge;
    TDate  riskStartDate;
    CxTRecoveryCurve *recoveryCurve = NULL;
    
    recoveryCurve = CxRecoveryCurveMakeFromRecoveryRate(recoveryRate);
    if (recoveryCurve == NULL)
        goto done; /* failure */

    /* upfront charge should be calculated on a knock-out on default basis
       for single names and no knock-out basis for index - hence
       the riskStartDate is different */

    if (deal->isIndex)
        riskStartDate = today;
    else
        riskStartDate = deal->exerciseDate;

    if (CxCdsPrice (riskStartDate,
                    deal->exerciseDate,
                    deal->exerciseDate,
                    deal->maturityDate,
                    0, /* delay */
                    deal->coupon,
                    deal->payAccOnDefault,
                    deal->feeInterval,
                    CX_SHORT_FRONT_STUB, /* stubType */
                    deal->dcc,
                    CX_BAD_DAY_NONE, /* badDayConv */
                    NULL,            /* calendar */
                    discCurve,
                    sprdCurve,
                    recoveryCurve,
                    FALSE, /* protectStart ??? */
                    FALSE, /* isPriceClean ??? */
                    &upfrontCharge) != SUCCESS)
        goto done; /* failure */

    if (deal->isIndex)
    {
        /* We need the probability that we knock out on default */
        /* This loss is then added to the upfront charge */

        double defaultProb;
        double adj;

        if (FwdDefaultProbability(deal,
                                  today,
                                  sprdCurve,
                                  &defaultProb) != SUCCESS)
            goto done; /* failure */

        adj = (1.0 - recoveryRate) * defaultProb;
        upfrontCharge += adj;
    }

    /* We quote the price with a par of 1.0 */
    *fwdPrice = 1.0 - upfrontCharge;

    status = SUCCESS;

 done:

    CxRecoveryCurveFree (recoveryCurve);

    if (status != SUCCESS)
        GtoErrMsgFailure(routine);
    
    return status;
}



static int ImpliedVolSolver (
    double              vol,
    IMPLIED_VOL_PARAMS *params,
    double             *priceDiff)
{
    int status = FAILURE;

    CrxTCdsOptionCalc* calc = CrxCdsOptionCalc (params->deal,
                                                params->distType,
                                                params->today,
                                                vol,
                                                params->qdist,
                                                params->discCurve,
                                                params->sprdCurve,
                                                params->recoveryRate);

    if (calc == NULL)
        goto done; /* failure */

    *priceDiff = params->price - calc->price;
    status = SUCCESS;

 done:

    CrxCdsOptionCalcFree (calc);
    return status;
}

/*
 * Computes the discount factor from the option payment date to today.
 */
static int FwdDiscountFactor(
    CrxTCdsOption    *deal,
    TDate             today,
    TCurve           *discCurve,
    double           *discountFactor)
{
    static char routine[] = "FwdDiscountFactor";

    if (GtoDiscountDateForward(today, 
                               deal->exerciseDate,
                               discCurve,
                               GTO_FLAT_FORWARDS, 
                               discountFactor) != SUCCESS)
        return GtoErrMsgFailure(routine);

    return SUCCESS;
}

/*
 * Computes the survival probability of the CDS at exercise, assuming
 * that survival probability at the end of today = 1.0
 */
static int FwdDefaultProbability(
    CrxTCdsOption    *deal,
    TDate             today,
    CxTCreditCurve   *sprdCurve,
    double           *defaultProb)
{
    static char routine[] = "FwdDefaultProbability";

    double survivalProb;

    if (FwdSurvivalProbability(deal, 
                               today, 
                               sprdCurve,
                               &survivalProb) != SUCCESS)
        return GtoErrMsgFailure(routine);

    *defaultProb = (1.0 - survivalProb);
    return SUCCESS;
}

/*
 * Computes the survival probability of the CDS at exercise, assuming
 * that survival probability at the end of today = 1.0
 */
static int FwdSurvivalProbability(
    CrxTCdsOption    *deal,
    TDate             today,
    CxTCreditCurve   *sprdCurve,
    double           *survivalProb)
{
    static char routine[] = "FwdSurvivalProbability";

    if (GtoDiscountDateForward(today, 
                               deal->exerciseDate, 
                               sprdCurve->tc, 
                               GTO_FLAT_FORWARDS,
                               survivalProb) != SUCCESS)
        return GtoErrMsgFailure(routine);

    return SUCCESS;
}
