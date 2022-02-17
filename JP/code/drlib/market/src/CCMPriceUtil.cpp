//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMPriceUtil.cpp
//
//   Description : CCM Price Utilities 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/Maths.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

#define TINY 1e-12

/* Tool to determine i such that Ti<= T <Ti+1  */
/* Uses an integer version of bisection        */
static int findIndex(
    double             T,     /* (I) Value of T to be localised in the array */ 
    const DoubleArray& t){    /* (I) Array of Loss events                    */ 
    int n = t.size(); // for ease
    double Tplus = T + 1e-14;
    if (Tplus < t[0])
        return -1;
    if (Tplus >= t[n-1])
        return n-1;
    int iUp  = n-1;
    int iLow = 0;
    int iMid = 0;
    while (iUp - iLow > 1)
    {
        iMid = (iLow + iUp)/2;    /* modulo division */
        if (Tplus >= t[iMid]){
            iLow = iMid;
        } else {
            iUp  = iMid;
        }
    }
    /* At this stage, data[iMid].date can be > or < T, and no one can tell */
    /* So here is a way around it */
    if (Tplus < t[iMid]){
        --iMid;
    }
    ASSERT(iMid >= 0 && iMid <n);
    ASSERT(t[iMid] <= Tplus && Tplus < t[iMid+1]);
    return iMid;
}

static double ccmZeroPriceF(
    const DoubleArray& tR,
    const DoubleArray& ZR,
    double             point){
    double r;
    int index = findIndex(point, tR);
    if (index == -1) /* point < tR[0] */
    {
        r = -1./tR[0]*log(ZR[0]);
        return exp(-r*point);
    }
    int nR = tR.size(); // for ease
    if (index == nR-1) /* no points in tR[] after point */
    {
        if (nR >= 2)
        {/* take last rate*/
            r = -1./(tR[nR-1]-tR[nR-2])*log(ZR[nR-1]/ZR[nR-2]);
            return ZR[nR-1]*exp(-r*(point-tR[nR-1]));
        }
        else /* only 1 point>0, and we return a ZC after it*/
        {
            r = -1./tR[0]*log(ZR[0]);
            return exp(-r*point);
        }
    }
    ASSERT(index+1<nR);
    r = -1./(tR[index+1]-tR[index])*log(ZR[index+1]/ZR[index]);
    return ZR[index]*exp(-r*(point-tR[index]));
}

static double ccmZeroPriceL(
    const DoubleArray& tR,
    const DoubleArray& ZR,
    double             point)
{
    int index = findIndex(point, tR);
    ASSERT(fabs(tR[0])<TINY); /* no default value for t=0 in Linear interp */
    if (index == -1){ /* point < tR[0] */
        throw ModelException("ccmZeroPriceL");
    }
    int nR = tR.size(); // for ease
    if (index == nR-1){ /* no points in tR[] after point */
        return ZR[nR-1]; /* convention: lambda=0 after last time point */
    }
    double r = -(ZR[index+1] - ZR[index]) / (tR[index+1]-tR[index]);
    return (ZR[index] - r*(point-tR[index]));
}

/* ---------------------------------------------------------------------
 * Payoff on contingent leg
 */
/* \int_a^b { exp(-rt) dt} */
static double IntegExp(double a, double b, double Za, double Zb)
{
    double r;
    if (fabs(b-a)<TINY) return 0.;
    r = -1./(b-a)*log(Zb/Za);
    if (fabs(r) < TINY) return (b-a)*Za;
    return (Za-Zb)/r;
}

/* \int_a^b { t exp(-rt) dt} */
static double IntegTExp(double a, double b, double Za, double Zb)
{
    double r;
    if (fabs(b-a)<TINY) return 0.;
    r = -1./(b-a)*log(Zb/Za);
    if (fabs(r) < TINY) return 0.5*(b-a)*(b-a)*Za;
    return ((r*a+1)*Za-(r*b+1)*Zb)/(r*r);
}

/** calculate \int_t0^t1 { Zr(t) dZlambda(t)} where Zlambda is flat fwd */
static double ccmCtgIntegF(
    double t0,    double t1,
    double ZR0,   double ZR1,
    double ZL0,   double ZL1)
{
    double a, integ;
    integ = IntegExp(t0,t1, ZR0*ZL0, ZR1*ZL1);
    a = fabs(integ) > TINY ? -1./(t1-t0)*log(ZL1/ZL0) : 0.;
    return a * integ;
}

/** calculate \int_t0^t1 { Zr(t) dZlambda(t)} where Zlambda is linear */
static double ccmCtgIntegL(
    double t0,    double t1,
    double ZR0,   double ZR1,
    double ZL0,   double ZL1)
{
    double a, integ;
    integ = IntegExp(t0,t1, ZR0, ZR1);
    a = fabs(integ) > TINY ? -(ZL1-ZL0)/(t1-t0) : 0.;
    return a*integ;
}

/** calculate \int_t0^t1 { t Zr(t) dZlambda(t)} where Zlambda is flat fwd */
static double ccmAccIntegF(
    double t0,    double t1,
    double ZR0,   double ZR1,
    double ZL0,   double ZL1)
{
    double a, integ;
    integ = IntegTExp(t0,t1, ZR0*ZL0, ZR1*ZL1);
    a = fabs(integ) > TINY ? -1./(t1-t0)*log(ZL1/ZL0) : 0.;
    return a*integ;
}

/** calculate \int_t0^t1 { t Zr(t) dZlambda(t)} where Zlambda is linear */
static double ccmAccIntegL(
    double t0,    double t1,
    double ZR0,   double ZR1,
    double ZL0,   double ZL1)
{
    double a, integ;
    integ = IntegTExp(t0,t1, ZR0, ZR1);
    a = fabs(integ) > TINY ? -(ZL1-ZL0)/(t1-t0) : 0.;
    return a*integ;
}

/**
 * Calculate the PV of a accrual on default leg with delay
 *     PV(T_s,T_e,d)  = \int_{T_s}^{T_e} {(a+bt) Z_r(t+d) dZ_\Lambda(t)}
 * With the following notations
 * Z_r(t)       & = & \exp(-\int_0^t{r_s ds})
 * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
 *
 * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
 * An interpolation method is chosen within I_i:
 * - riskless rate r_t are piecewise constant (flat forward)
 * - either \lambda_t is piecewise constant (FLATFWD expected loss)
 *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
 */
double CCMPriceUtil::contingentAndAccrPrice( 
    double             Ts,   /* (I) starting time                      */
    double             Te,   /* (I) end time                           */
    double             d,    /* (I) payment delay                      */
    const DoubleArray& tR,   /* (I) timeline for Zr                    */
    const DoubleArray& ZR,   /* (I) discount factor                    */
    const DoubleArray& tL,   /* (I) timeline for Zl                    */
    const DoubleArray& ZL,   /* (I) expected survival                  */
    double             a,    /* (I) ctg leg factor                     */
    double             b,    /* (I) acc payoff factor                  */
    ExpLossType        type) /* (I) FLATFWD or LINEAR                  */
{
    const char routine[] = "CCMPriceUtil::contingentAndAccrPrice";
    int i, j, n; 
    if (tR[0] < TINY) {
        throw ModelException(routine, "IR timeline must not start with a 0");
    }

    if (fabs(Te-Ts) < TINY) {
        return 0.0;
    }

    /* 2. calculate rates wrt tR
     * index choice: r[i] is the fwd btw tR[i-1] and tR[i], with the
     * convention that tR[-1] would be 0*/
    int nR = tR.size(); // for ease
    int nL = tL.size(); // for ease
    DoubleArray r(nR+1);
    r[0] = -1./tR[0]*log(ZR[0]);
    for (i=1; i<nR; ++i){
        r[i] = -1./(tR[i]-tR[i-1])*log(ZR[i]/ZR[i-1]);
    }
    /* timeline is assumed to be sorted and all points unique */
    r[nR] = r[nR-1];

    /* 3. calculate Zr(d) */
    for (i=0; i<nR; ++i) {
        if (tR[i]>d-TINY){
            break;
        }
    }

    DoubleArray Zrd(nR+1);
    Zrd[0] = (i>0) ? ZR[i-1] * exp(-r[i]*(d-tR[i-1]))
        : exp(-r[0]*d);

    /* 4. displaced timeline so that r(u+d) is piecewise constant */
    /* td[i+1] = max(tR[i]-d,0)
       index convention: r[i] represents r(u+d) for any u btw td[i] and td[i+1]
    */
    DoubleArray td(nR+1);
    td[0] = 0.;
    for (i=0; i<nR; ++i){
        td[i+1] = Maths::max(0.0, tR[i]-d);
    }
    /* 5. series of Zr(td+d), td[i] = max(tR[i]-d,0) */
    for (i=0; i<nR; ++i){
        Zrd[i+1] = Zrd[i] * exp(-r[i]*(td[i+1]-td[i]));
    }
    /* 6. merged timeline for ZL and Zrd */
    DoubleArray t(nR+nL+3);
    for (i=0; i<nR+1; ++i){
        t[i] = td[i];
    }
    for (i=0; i<nL; ++i){
        t[i+nR+1] = tL[i];
    }
    t[nR+nL+1] = Ts;
    t[nR+nL+2] = Te;

    sort(t.begin(), t.end());
    /* remove negative times */
    for (j=0; j<nR+nL+3; ++j){
        if (t[j]>-TINY){
            break;    
        }
    }
    n = j;

    /* remove duplicates */
    t[0] = t[n];
    j = 0;
    for (i=n+1; i<nR+nL+3; ++i)
    {
        if (t[i] > t[j]+TINY)
        {
            ++j;
            t[j] = t[i];
        }
    }
    n = j+1;  /* the last (nR+nL-n) elements are no longer used */ 
    t.resize(n); // to do: review
    int iStart = findIndex(Ts+TINY, t);
    int iEnd   = findIndex(Te+TINY, t);

    double tmp = 0.0;
    for (i=iStart; i<iEnd; ++i) {
        double pv1, pv2;
        double t0  = t[i]; 
        double t1  = t[i+1];
        double ZR0 = ccmZeroPriceF(td, Zrd, t0);
        double ZR1 = ccmZeroPriceF(td, Zrd, t1);
        if (type == FLATFWD) {
            double ZL0 = ccmZeroPriceF(tL, ZL, t0);
            double ZL1 = ccmZeroPriceF(tL, ZL, t1);
            pv1 = fabs(a)>TINY ? ccmCtgIntegF(t0, t1, ZR0, ZR1, ZL0, ZL1) : 0.;
            pv2 = fabs(b)>TINY ? ccmAccIntegF(t0, t1, ZR0, ZR1, ZL0, ZL1) : 0.;
            tmp += a*pv1 + b*pv2; 
        } else {
            double ZL0 = ccmZeroPriceL(tL, ZL, t0);
            double ZL1 = ccmZeroPriceL(tL, ZL, t1);
            pv1 = fabs(a)>TINY ? ccmCtgIntegL(t0, t1, ZR0, ZR1, ZL0, ZL1) : 0.;
            pv2 = fabs(b)>TINY ? ccmAccIntegL(t0, t1, ZR0, ZR1, ZL0, ZL1) : 0.;
            tmp += a*pv1 + b*pv2; 
        }
    }
    return tmp; 
}


double CCMPriceUtil::contingentPrice( 
    double             Ts,   /* (I) starting time                      */
    double             Te,   /* (I) end time                           */
    double             d,    /* (I) payment delay                      */
    const DoubleArray& tR,   /* (I) timeline for Zr                    */
    const DoubleArray& ZR,   /* (I) discount factor                    */
    const DoubleArray& tL,   /* (I) timeline for Zl                    */
    const DoubleArray& ZL,   /* (I) expected survival                  */
    const ExpLossType  type) /* (I) FLATFWD or LINEAR                  */
{
    return CCMPriceUtil::contingentAndAccrPrice (
        Ts,   /* (I) starting time                           */
        Te,   /* (I) end time                                */
        d,    /* (I) payment delay                           */
        tR,   /* (I) initial timeline for IR curve           */
        ZR,   /* (I) zero coupons                            */
        tL,   /* (I) initial timeline for Zl                 */
        ZL,   /* (I) credit zeros (survival)                 */
        1.,
        0.,
        type); /* (I) F or L for Flat of Linear interpolation */
}

double CCMPriceUtil::accrualOnDefault( 
    double            Ts,   /* (I) starting time                           */
    double            Te,   /* (I) end time                                */
    double            d,    /* (I) payment delay                           */
    const DoubleArray& tR,     /* (I) timeline for Zr                    */
    const DoubleArray& ZR,     /* (I) discount factor                    */
    const DoubleArray& tL,     /* (I) timeline for Zl                    */
    const DoubleArray& ZL,     /* (I) expected survival                  */
    ExpLossType        type)  /* (I) FLATFWD or LINEAR */
{
    return CCMPriceUtil::contingentAndAccrPrice (
        Ts,   /* (I) starting time                           */
        Te,   /* (I) end time                                */
        d,    /* (I) payment delay                           */
        tR,   /* (I) initial timeline for IR curve           */
        ZR,   /* (I) zero coupons                            */
        tL,   /* (I) initial timeline for Zl                 */
        ZL,   /* (I) credit zeros (survival)                 */
        -Ts,
        1.,
        type); /* (I) F or L for Flat of Linear interpolation */
}

/** calculate \int_t0^t1 { t Zlambda(t) dt} */
static double integAvgNtlF(double t0, double t1, double ZL0, double ZL1){
    return IntegTExp(t0,t1, ZL0, ZL1);
}

/** calculate \int_t0^t1 { t Zlambda(t) dt} */
static double integAvgNtlL(double t0, double t1, double ZL0, double ZL1) {
    /* a = -(ZL1-ZL0)/(t1-t0); */
    return ((t1-t0)*ZL0*(t1+t0)/2 + (ZL1-ZL0)*(t1*t1+t0*t1+t0*t0)/3);
}

/**
 * Calculate the expected average notional
 *     avgNtl(T_s,T_e)  = 
 \frac{1}{T_e-T_s} \int_{T_s}^{T_e} {t Z_\Lambda(t) dt}
 * With the following notations
 * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
 *
 * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
 * An interpolation method is chosen within I_i:
 * - either \lambda_t is piecewise constant (FLATFWD expected loss)
 *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
 */
double CCMPriceUtil::averageNotional( 
    double             Ts,     /* (I) ctg leg begin                     */
    double             Te,     /* (I) ctg leg end                       */
    const DoubleArray& tL,     /* (I) timeline for Zl                   */
    const DoubleArray& ZL,     /* (I) expected survival                 */
    ExpLossType        type){  /* (I) FLATFWD or LINEAR*/
    const char routine[] = "ccmAverageNotional";
    int i, j, n; 

    if (Te-Ts < TINY) {
        throw ModelException(routine, "Average notional requested for a "
                             "period of length 0");
    }

    /* merged timeline for ZL and Ts and Te */
    int nL = tL.size(); // for ease
    DoubleArray t(nL+2);
    for (i=0; i<nL; ++i){
        t[i] = tL[i];
    }
    t[nL]   = Ts;
    t[nL+1] = Te;

    sort(t.begin(), t.end());
    /* remove negative times */
    for (j=0; j<nL+2; ++j){
        if (t[j]>-TINY){
            break;    
        }
    }
    n = j;

    /* remove duplicates */
    t[0] = t[n];
    j = 0;
    for (i=n+1; i<nL+2; ++i) {
        if (t[i] > t[j]+TINY) {
            ++j;
            t[j] = t[i];
        }
    }
    n = j+1;  /* the last (nL-n) elements are no longer used */ 
    t.resize(n); // to do: review

    int iStart = findIndex(Ts+TINY, t);
    int iEnd   = findIndex(Te+TINY, t);

    double tmp=0.0;
    for (i=iStart; i<iEnd; ++i) {
        double bucket=0;
        double t0  = t[i]; 
        double t1  = t[i+1];
        if (type == FLATFWD) {
            double ZL0 = ccmZeroPriceF(tL, ZL, t0);
            double ZL1 = ccmZeroPriceF(tL, ZL, t1);
            bucket = integAvgNtlF(t0, t1, ZL0, ZL1);
            tmp += bucket; 
        } else {
            double ZL0 = ccmZeroPriceL(tL, ZL, t0);
            double ZL1 = ccmZeroPriceL(tL, ZL, t1);
            bucket = integAvgNtlL(t0, t1, ZL0, ZL1);
            tmp += bucket; 
        }
    }

    ASSERT(Te>Ts || fabs(tmp) < TINY);
    return (tmp / (Te-Ts)); 
}

#undef TINY

/**
 * Linear Loss interpolation
 *  */
double CCMPriceUtil::linearLossInterpolate(
    const DateTime&      maturity,
    const DateTimeArray& timeline,
    const DoubleArray&   expectedTrancheNotional) 
{
    const char routine[] = "linearInterpolate";
    int n = timeline.size();
    
    /* Intercept cases when maturity is in the timeline */
    for (int i = 0; i < n; ++i) {
        if (maturity == timeline[i]) {
            return expectedTrancheNotional[i];
        }
    }
    if (maturity > timeline[n-1]) {
        return expectedTrancheNotional[n-1];
    } else { /* means there is a point before and a point after */
        int i = 0; /* like FindIndex on longs, but knowing 
                      it is never overlapping */
        while (maturity>timeline[i] && i<n) {
            ++i;
        }
        --i;
        if (i < 0 || i >= n){
            throw ModelException(routine, "Internal error");
        }

        return (expectedTrancheNotional[i]+timeline[i].daysDiff(maturity)
                / ((double)(timeline[i].daysDiff(timeline[i+1])))
                * (expectedTrancheNotional[i+1]-expectedTrancheNotional[i]));
    }
}


double CCMPriceUtil::pastLossesInPeriod(
    const CashFlowArray& pastTrancheLosses,
    const DateTime&      periodStart,         // point to be returned
    const DateTime&      periodEnd,
    const BoolArray&     payPastTrancheLosses)
{
    double lossesInPeriod = 0.0;
    int numOfPastTrancheLosses = pastTrancheLosses.size();

    // Find the losses happening inside the protection period
    for (int i=0; i < numOfPastTrancheLosses; ++i) {
        if (pastTrancheLosses[i].date > periodStart) {
            // This loss happens after periodStart, so consider it

            if (pastTrancheLosses[i].date > periodEnd) {
                // This loss hapens after periodEnd, so there are no more losses
                // in the protection period
                break;
            }
            // This loss falls inside the protection period - take it into 
            // consideration if its payPastTrancheLosses is true
            if (payPastTrancheLosses[i]) {
                lossesInPeriod += (i == 0) ?
                    pastTrancheLosses[i].amount :
                    pastTrancheLosses[i].amount - pastTrancheLosses[i-1].amount;
            }
        }
    }
    return lossesInPeriod;
}


double CCMPriceUtil::expectedLossesInPeriod(
    const double                outstandingNotional, /* initialTrancheSize -
                                                      * pastTrancheLoss */
    const DateTime&             periodStart,
    const DateTime&             periodEnd,
    const IDiscountCurveRiskySP effectiveCurve)
{
    double lossesFromPotentialDefaults = 0.0;
    if (!Maths::isZero(outstandingNotional)) {
        const double loss1 =
            outstandingNotional * effectiveCurve->survivalProb(periodStart);

        const double loss2 =
            outstandingNotional * effectiveCurve->survivalProb(periodEnd);

        lossesFromPotentialDefaults = loss1 - loss2;
    }

    return lossesFromPotentialDefaults;
}


/* Returns the expected losses (in the fee leg) between today and maturity, 
   using the effectiveCurve */
double CCMPriceUtil::expectedFeeNotionalLoss(
    const double                outstandingNotional,
    const DateTime&             valueDate,           // reference date
    const DateTime&             maturity,            // point to be returned
    const IDiscountCurveRiskySP effectiveCurve)      // effective curve to
{
    double expectedNotionalLoss = 0.0;

    if ((maturity > valueDate) && !Maths::isZero(outstandingNotional)) {
        // Future dates involve a "real" calculation
        // Estimate the losses in the remaining "outstandingNotional"
        expectedNotionalLoss = outstandingNotional * 
            (1 - effectiveCurve->survivalProb(maturity));
    }
    return expectedNotionalLoss;
}


double CCMPriceUtil::outstandingNotionalAtMaturity(
    const double         initialNotional,
    const DateTime&      maturity,            // point to be returned
    const CashFlowArray& pastTrancheLosses)   // effective curve to
{
    int i = 0;
    for (; i < pastTrancheLosses.size(); ++i) {
        if (maturity < pastTrancheLosses[i].date) {
            break;
        }
    }

    return (i == 0? initialNotional : 
                    initialNotional - pastTrancheLosses[i-1].amount);
}


DRLIB_END_NAMESPACE
