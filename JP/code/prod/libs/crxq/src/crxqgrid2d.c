/******************************************************************************
 * Module:      CRXQ
 * Submodule:
 * File:        crxqgrid2d.c        
 * Function:    
 * Author:      Credit QRD
 * Revision:    $Header: $
 *****************************************************************************/
#include <math.h>
#include <malloc.h>
#include "crxq.h"
#include <crxmultiq/include/crmultiq.h>
#include <crxflow/include/crxerror.h>

#if !defined(NULL)
#define NULL 0
#define UNDEF_NULL
#endif

/*-----------------------------------------------------------------------------
 * Q3MCPricer
 *
 * Integrate payoff function using Monte Carlo method.
 *
 */
int Q3MCPricer(
    long      dimX,    /* (I) Number of variables.        */
    CRXQDATA  **smX,     /* (I) Measure data for variables. */ 
    double   *rhoX,    /* (I) Correlation of variables.   */
    CRXQPAYOFF   *pf,      /* (I) Payoff parameters.          */
    CRXQFPAYOFF  *payFunc, /* (I) Pointer to payoff function. */ 
    long      genTyp,  /* (I) Generator type.             */
    long      npts,    /* (I) Number of random samples.   */ 
    double   *result   /* (O) value, error estimate       */ 
    )
{

    static char routine[] = "Q3MCPricer";
    double      muX[CRXQ_DEV_DIM_MAX];
    double      sigX[CRXQ_DEV_DIM_MAX];
    double      yldX[CRXQ_DEV_DIM_MAX];
    double      fcnVal;
    double      fcnAvgVal;
    double      fcnAvgSqVal;
    double      errEst;
    double     *pts = NULL, *ppts = NULL;
    long        i, j;
    int         status = FAILURE;

    /* Check inputs. */
    if (dimX > CRXQ_DEV_DIM_MAX || dimX < 1 || npts < 1) goto RETURN;
    if (smX == NULL || payFunc == NULL) goto RETURN;
    if (genTyp != CRXQ_BV_RAN2 && genTyp != CRXQ_BV_SOBOL) goto RETURN;
    
    /* Allocate space for sample points. */
    pts = (double*) DR_Array(DOUBLE, 0, dimX * npts - 1);
    if (pts == NULL)
    {
        goto RETURN;
    }

    /* Collect variable mus and sigmas for call to sample generator. */
    for (i = 0; i < dimX; i++)
    {
        muX[i]  = (smX[i])->muMQ;
        sigX[i] = (smX[i])->sigMQ;
    }
     
    /* Generate independent n(0,1) sample points with dimensionality dimX. */
    if (BoxMuller(
        genTyp,
        dimX,
        npts,
        pts) == FAILURE) goto RETURN;

    /* Transform to n(mu,sigma) sample points with correlation rho. */
    if (Gauss(
        dimX,
        npts,
        sigX,
        muX,
        rhoX,
        pts) == FAILURE) goto RETURN;
    
    fcnAvgVal   = 0.;
    fcnAvgSqVal = 0;
    ppts        = pts;
    for (i =0; i < npts; i++)
    {
        /* Map sample points into yields. */
        for (j = 0; j < dimX; j++) 
        {
            if (CRXQMap(
                smX[j],
                ppts[j],
                &(yldX[j])) == FAILURE) goto RETURN;
        }

        /* Evaluate function value over sample points. */
        if ((*payFunc)(
            pf,
            yldX,
            1.0,
            &fcnVal) == FAILURE) goto RETURN;

        fcnAvgVal   += fcnVal;
        fcnAvgSqVal += fcnVal * fcnVal;

        /* Next sample point. */
        ppts += dimX;

    }

    fcnAvgVal   /= npts;
    fcnAvgSqVal /= npts;
    errEst = sqrt(fcnAvgSqVal - fcnAvgVal * fcnAvgVal);

    /* Return value. */
    result[0] = fcnAvgVal;

    status = SUCCESS;

 RETURN:

    /* Deallocate. */
    if (pts)
        Free_DR_Array(pts, DOUBLE, 0, dimX * npts - 1);

    return status;

}

/*f----------------------------------------------------------------------------
 * CRXQPay2D_FlrFlrOrCapCap
 *
 * expects params = (spread, strike, leverage)
 */
int CRXQPay2D_FlrFlrOrCapCap(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double spread   = (prm->params)[0];
    double strike   = (prm->params)[1];
    double leverage = (prm->params)[2];
    double cop      = prm->cop;
    smooth;

    /* apply leverage and spread */
    yield1 *= leverage;
    yield1 += spread;

    if (cop > 0.)
        *payoff = MAX(MAX(yield2, yield1), strike) - yield2;
    else
        *payoff = yield2 - MIN(MIN(yield2, yield1), strike);

    return SUCCESS;
}


/*f----------------------------------------------------------------------------
 * CRXQPay2D_FlrFlrOrCapCapEmbedFlt
 *
 * expects params = (spread, strike, leverage)
 */
int CRXQPay2D_FlrFlrOrCapCapEmbedFlt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double spread   = (prm->params)[0];
    double strike   = (prm->params)[1];
    double leverage = (prm->params)[2];
    double cop      = prm->cop;
    smooth;

    /* apply leverage and spread */
    yield1 *= leverage;
    yield1 += spread;

    if (cop > 0.)
        *payoff = MAX(MAX(yield2, yield1), strike);
    else
        *payoff = MIN(MIN(yield2, yield1), strike);

    return SUCCESS;
}


/*f----------------------------------------------------------------------------
 * CRXQPay2D_FlrCapOrCapFlr
 *
 * expects params = (spread, strike, leverage)
 */
int CRXQPay2D_FlrCapOrCapFlr(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double spread   = (prm->params)[0];
    double strike   = (prm->params)[1];
    double leverage = (prm->params)[2];
    double cop      = prm->cop;
    smooth;

    /* apply leverage and spread */
    yield1 *= leverage;
    yield1 += spread;

    if (cop > 0.)
        *payoff = MIN(MAX(yield2, yield1), strike) - yield2;
    else
        *payoff = yield2 - MAX(MIN(yield2, yield1), strike);

    return SUCCESS;
}


/*f----------------------------------------------------------------------------
 * CRXQPay2D_FlrCapOrCapFlrEmbedFlt
 *
 * expects params = (spread, strike, leverage)
 */
int CRXQPay2D_FlrCapOrCapFlrEmbedFlt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double spread   = (prm->params)[0];
    double strike   = (prm->params)[1];
    double leverage = (prm->params)[2];
    double cop      = prm->cop;
    smooth;

    /* apply leverage and spread */
    yield1 *= leverage;
    yield1 += spread;

    if (cop > 0.)
        *payoff = MIN(MAX(yield2, yield1), strike);
    else
        *payoff = MAX(MIN(yield2, yield1), strike);

    return SUCCESS;
}


/*f----------------------------------------------------------------------------
 * CRXQPay2D_Sum
 *
 * expects params = (strike, leverage)
 */
int CRXQPay2D_Sum(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double strike   = (prm->params)[0];
    double leverage = (prm->params)[1];
    double cop      = prm->cop;
    smooth;

    /* apply leverage */
    yield1 *= leverage;

    *payoff = MAX(cop * (yield2 + yield1 - strike), 0.);

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_Prod
 *
 * expects params = (strike)
 */
int CRXQPay2D_Prod(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double strike   = (prm->params)[0];
    double cop      = prm->cop;
    smooth;

    *payoff = MAX(cop * (yield2 * yield1 - strike), 0.);

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_Perc
 *
 * expects params = (strike)
 */
int CRXQPay2D_Perc(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double strike   = (prm->params)[0];
    double cop      = prm->cop;
    smooth;

    *payoff = MAX(cop * yield1 * (yield2 - strike), 0.);

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_PercWgt
 *
 * expects params = (weight1, weight2, strike)
 */
int CRXQPay2D_PercWgt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yld1   = pt[0];
    double yld2   = pt[1];
    double wgt1   = (prm->params)[0];
    double wgt2   = (prm->params)[1];
    double strike = (prm->params)[2];
    double cop    = prm->cop;
    smooth;

    *payoff = MAX(cop * yld1 * (wgt1 * yld1 + wgt2 * yld2 - strike), 0.);

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_YldYld
 *
 */
int CRXQPay2D_YldYld(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    smooth; prm;

    *payoff = yield2 * yield1;

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_YldNull_1D
 *
 * Test function. Prices yield on first variable, ignoring the second.
 */
int CRXQPay2D_YldNull(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    smooth; prm;
    
    *payoff = yield;
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_NullYld
 *
 * Test function. Prices yield on second variable, ignoring the first.
 */
int CRXQPay2D_NullYld(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[1];
    smooth; prm;
    
    *payoff = yield;
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_VnlNull
 *
 * Test function. Prices vanilla on first variable, ignoring the second.
 *
 * expects params = (strike)
 */
int CRXQPay2D_VnlNull(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double strike = (prm->params)[0];
    smooth;

    *payoff = MAX((yield - strike) * prm->cop,0.);
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay2D_NullVnl
 *
 * Test function. Prices vanilla on second variable, ignoring the first.
 *
 * expects params = (strike) 
 */
int CRXQPay2D_NullVnl(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[1];
    double strike = (prm->params)[0];
    smooth;

    *payoff = MAX((yield - strike) * prm->cop,0.);
    return SUCCESS;
}


/*f----------------------------------------------------------------------------
 * CRXQPay2D_MinMaxIn
 *
 * expects params = (LB, HB, leverage, spread, F, C, low epsilon, high epsilon)
 */
int CRXQPay2D_MinMaxIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double lb       = (prm->params)[0];
    double hb       = (prm->params)[1];
    double leverage = (prm->params)[2];
    double spread   = (prm->params)[3];
    double floor    = (prm->params)[4];
    double cap      = (prm->params)[5];
    double leps     = (prm->params)[6];
    double heps     = (prm->params)[7];

    double weight   = 0.;

    smooth;

    weight = Q3ProbInBarriers(
        yield2,
        lb,
        hb,
        leps,
        heps,
        smooth);

    *payoff = weight * MIN(MAX(leverage*yield1 + spread, floor), cap);

    return SUCCESS;

} /* CRXQPay2D_MinMaxIn */


/*f----------------------------------------------------------------------------
 * CRXQPay2D_MinMaxOut
 *
 * expects params = (LB, HB, leverage, spread, F, C, low epsilon, high epsilon)
 */
int CRXQPay2D_MinMaxOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double lb       = (prm->params)[0];
    double hb       = (prm->params)[1];
    double leverage = (prm->params)[2];
    double spread   = (prm->params)[3];
    double floor    = (prm->params)[4];
    double cap      = (prm->params)[5];
    double leps     = (prm->params)[6];
    double heps     = (prm->params)[7];

    double weight   = 0.;

    smooth;

    weight = Q3ProbOutBarriers(
        yield2,
        lb,
        hb,
        leps,
        heps,
        smooth);

    *payoff = weight * MIN(MAX(leverage*yield1 + spread, floor), cap);

    return SUCCESS;

} /* CRXQPay2D_MinMaxOut */


/*f----------------------------------------------------------------------------
 * CRXQPay2D_MinMaxSpdIn
 *
 * expects params = (LB, HB, leverage, spread, F, C, low epsilon, high epsilon)
 */
int CRXQPay2D_MinMaxSpdIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double lb       = (prm->params)[0] + yield1;
    double hb       = (prm->params)[1] + yield1;
    double leverage = (prm->params)[2];
    double spread   = (prm->params)[3];
    double floor    = (prm->params)[4];
    double cap      = (prm->params)[5];
    double leps     = (prm->params)[6];
    double heps     = (prm->params)[7];

    double weight   = 0.;

    smooth;

    weight = Q3ProbInBarriers(
        yield2,
        lb,
        hb,
        leps,
        heps,
        smooth);

    *payoff = weight * MIN(MAX(leverage*yield1 + spread, floor), cap);

    return SUCCESS;
} /* CRXQPay2D_MinMaxSpdIn */


/*f----------------------------------------------------------------------------
 * CRXQPay2D_MinMaxSpdOut
 *
 * expects params = (LB, HB, leverage, spread, F, C, low epsilon, high epsilon)
 */
int CRXQPay2D_MinMaxSpdOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1   = pt[0];
    double yield2   = pt[1];
    double lb       = (prm->params)[0] + yield1;
    double hb       = (prm->params)[1] + yield1;
    double leverage = (prm->params)[2];
    double spread   = (prm->params)[3];
    double floor    = (prm->params)[4];
    double cap      = (prm->params)[5];
    double leps     = (prm->params)[6];
    double heps     = (prm->params)[7];

    double weight   = 0.;

    smooth;

    weight = Q3ProbOutBarriers(
        yield2,
        lb,
        hb,
        leps,
        heps,
        smooth);

    *payoff = weight * MIN(MAX(leverage*yield1 + spread, floor), cap);

    return SUCCESS;
} /* CRXQPay2D_MinMaxSpdOut */


/*-----------------------------------------------------------------------------
 * CRXQPay2D_InAndIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 and LB2<yCnd<HB2 (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_InAndIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbInBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff = price;

    price = Q3ProbInBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff *= price;

    return SUCCESS;
} /* CRXQPay2D_InAndIn */

/*-----------------------------------------------------------------------------
 * CRXQPay2D_InOrIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 Or LB2<yCnd<HB2 (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_InOrIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbInBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff = price;

    price = Q3ProbInBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = *payoff + price - (*payoff) * price;

    return SUCCESS;
}/*CRXQPay2D_InOrIn*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_InAndOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) and LB2<yCnd<HB2 (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_InAndOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbInBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbOutBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff *= price;
    return SUCCESS;
}/*CRXQPay2D_InAndOut*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_InOrOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) or LB2<yCnd<HB2 (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_InOrOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbInBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbOutBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff = *payoff + price - (*payoff) * price;
    return SUCCESS;
}/*CRXQPay2D_InOrOut*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_OutAndIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 and (yCnd<LB2 or yCnd>HB2) (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_OutAndIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbOutBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbInBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff *= price;
    return SUCCESS;
}/*CRXQPay2D_OutAndIn*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_OutOrIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 Or (yCnd<LB2 or yCnd>HB2) (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_OutOrIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbOutBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbInBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff = *payoff + price - (*payoff) * price;
    return SUCCESS;
}/*CRXQPay2D_OutOrIn*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_OutAndOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) and (yCnd<LB2 or yCnd>HB2) (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_OutAndOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbOutBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbOutBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff *= price;
    return SUCCESS;
}/*CRXQPay2D_OutAndOut*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_OutOrOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) Or (yCnd<LB2 or yCnd>HB2) (with epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay2D_OutOrOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];

    double price;

    price = Q3ProbOutBarriers(
        yield2,
        lb2,
        hb2,
        leps2,
        heps2,
        smooth);
    *payoff = price;

    price = Q3ProbOutBarriers(
        yield1,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);
    *payoff = *payoff + price - (*payoff) * price;
    return SUCCESS;
}/*CRXQPay2D_OutOrOut*/

/*-----------------------------------------------------------------------------
 * CRXQPay2D_ComplexSPD
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * payoff = max(K0, a00*R1 + a10*R2 + 
 *                  b0 * max(K1, a01*R1 + a11*R2) +
 *                  b1 * max(K2, a02*R1 + a12*R2) +
 *                  b2 * max(K3, a03*R1 + a13*R2) +
 *                  b3 * max(K4, a04*R1 + a14*R2))
 */
int CRXQPay2D_ComplexSPD(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield1    = pt[0];
    double yield2    = pt[1];

    double K[5], b[4], a[2][5];
    int    i,j,idx;

    smooth;

    idx = 0;
    for(i=0;i<5;i++)
    {
        K[i] = prm->params[idx++];
    }

    for(i=0;i<5;i++)
    {
        for(j=0;j<2;j++)
        {
            a[j][i] = prm->params[idx++];
        }
    }

    for(i=0;i<4;i++)
    {
        b[i] = prm->params[idx++];
        if(fabs(b[i]) < TINY) /* b[i] == 0*/
        {
            K[i+1] = 0.0;
            a[0][i+1] = 0.0;
            a[1][i+1] = 0.0;
        }
    }
    *payoff = MAX(K[0], a[0][0] * yield1 + a[1][0] * yield2 +
                  b[0] * MAX( K[1], a[0][1] * yield1 + a[1][1] * yield2) +
                  b[1] * MAX( K[2], a[0][2] * yield1 + a[1][2] * yield2) +
                  b[2] * MAX( K[3], a[0][3] * yield1 + a[1][3] * yield2) +
                  b[3] * MAX( K[4], a[0][4] * yield1 + a[1][4] * yield2));
    return SUCCESS;
}/*CRXQPay2D_ComplexSPD*/




/*f----------------------------------------------------------------------------
 * DR_Array  
 *
 * Allocation of memory for an array of arbitrary type. Returns void*.
 */

void    *DR_Array ( int     type,   /* (I) Type         */
                    int     nl,     /* (I) Lower bound  */
                    int     nh)     /* (I) Higher bound */
{

    switch (type) 
    {
    case INT: 
        {
            int *v;

            v = (int *) calloc ((unsigned) (nh-nl+1), sizeof (int));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case INT_PTR: 
        {
            int **v;

            v = (int **) calloc ((unsigned) (nh-nl+1), sizeof (int *));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case INT_D_PTR: 
        {
            int ***v;

            v = (int ***) calloc ((unsigned) (nh-nl+1), sizeof (int **));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case LONG:
        {
            long *v;

            v = (long *) calloc ((unsigned) (nh-nl+1), sizeof (long));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case LONG_PTR: 
        {
            long **v;

            v = (long **) calloc ((unsigned) (nh-nl+1), sizeof (long *));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case LONG_D_PTR: 
        {
            long ***v;

            v = (long ***) calloc ((unsigned) (nh-nl+1), sizeof (long **));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case DOUBLE: 
        {
            double *v;

            v = (double *) calloc ((unsigned) (nh-nl+1), sizeof (double));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case DOUBLE_PTR: 
        {
            double **v;

            v = (double **) calloc ((unsigned) (nh-nl+1), 
                                    sizeof (double *));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }                          
    case DOUBLE_D_PTR: 
        {
            double ***v;

            v = (double ***) calloc ((unsigned) (nh-nl+1), 
                                     sizeof (double **));

            if (v == NULL) 
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));
        }
    case CHAR:
        {
            char *v;

            v = (char *) calloc ((unsigned) (nh-nl+1), sizeof (char));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        } 
    case CHAR_PTR:
        {
            char **v;

            v = (char **) calloc ((unsigned) (nh-nl+1), sizeof (char *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }
        

            return ((void *) (v-nl));
        }
    default: 
        {
            DR_Error ("Program bug: unsupported type in DR_Array!");
            return (NULL);
        }                          
    }  /* switch */

}  /* DR_Array */


/*f----------------------------------------------------------------------------
 * Free_DR_Array  
 *
 * Free DR Array
 */

int     Free_DR_Array ( void   *Array,  /* (I) Array        */
                        int     type,   /* (I) Type         */
                        int     nl,     /* (I) Lower bound  */
                        int     nh)     /* (I) Higher bound */
{
    /* To avoid warning message */
    nh += 0;

    switch (type) 
    {
    case INT: 
        {
            int *v = (int *) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case INT_PTR: 
        {
            int **v = (int **) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case INT_D_PTR: 
        {
            int ***v = (int ***) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case LONG: 
        {
            long *v = (long *) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case LONG_PTR: 
        {
            long **v = (long **) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case LONG_D_PTR: 
        {
            long ***v = (long ***) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case DOUBLE: 
        {
            double *v = (double *) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case DOUBLE_PTR: 
        {
            double **v = (double **) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case DOUBLE_D_PTR: 
        {
            double ***v = (double ***) Array;

            if (v != NULL) 
            {
                free ((char *) (v+nl));
            }

            break;
        }
    case CHAR:
        {
            char *v = (char *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
    case CHAR_PTR:
        {
            char **v = (char **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
    default: 
        {
            DR_Error("Program bug: unsupported type in Free_DR_Array!");
            return (FAILURE);
        }                          
    }  /* switch */

    Array = NULL;

    return (SUCCESS);

}  /* Free_DR_Array */


#ifdef UNDEF_NULL
#undef NULL
#endif


