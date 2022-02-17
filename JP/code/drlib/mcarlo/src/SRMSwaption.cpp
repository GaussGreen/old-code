//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMSwaption.cpp
//
//   Description : Helper for SRM - used for calibration against a swaption
//
//   Author      : Mark A Robson
//
//   Date        : 11 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMSwaption.hpp"
#include "edginc/SRMBlack.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/ModifiedNewtonRaphson.hpp"


DRLIB_BEGIN_NAMESPACE
#define TAU_MONTH 30
#define QSMALL    0.01
#define QCUTOFF       1E-4 
#define VOLPRECISION  0.0001

//// builds an empty SRMSwaption
SRMSwaption::SRMSwaption(): idx(0), prevIdx(0){}

/* Computes the deterministic squared mu factors (see section 2.2).
   From swapvol::MuFactor */
void SRMSwaption::muFactor(const SRMRatesHJMUtil&  ratesHJMUtil,
                           const DateTime&         currDate,
                           const DateTime&         swapStart,
                           const DateTime&         swapMat,
                           double&                 muFacSquare1, /* (O) */
                           double&                 muFacSquare2  /* (O) */) const{
    /* Zero coupon bond to swaption expiry */
    double ZeroToSt = ratesHJMUtil.pv(swapStart, true /* diffused curve */);
    /* Zero coupon bond to swap maturity */
    double ZeroToMat = ratesHJMUtil.pv(swapMat, true);
    /* factor weight */
    double alpha1 = ratesHJMUtil.getAlpha(0);
    vector<double> gfac(ratesHJMUtil.numFactors()); // reserve some space
    switch (ratesHJMUtil.numFactors()){
    case 1: {
        /* calculate mu1*/
        double mu1 = 0.0; /* mu value for each factor (lognormal assumption) */
        for (int j = 1; j <= NumCpns;j++) {
            double gfac = ratesHJMUtil.gFactor(currDate,CpnPayDates[j], 0);
            mu1   += gfac * weights[j] * alpha1 ;
        } /* for j*/
        double mu1_N = mu1 * ParYield; /* mu value for each factor 
                                          (normal assumption)  */
        double gfac =  ratesHJMUtil.gFactor(currDate, swapStart, 0);
        mu1   -=  ZeroToSt *  gfac * alpha1 / ( ZeroToSt - ZeroToMat);
        mu1_N -=  ZeroToSt *  gfac * alpha1 / Annuity;
        
        gfac = ratesHJMUtil.gFactor(currDate, swapMat, 0);
        mu1   +=  ZeroToMat * gfac * alpha1 / (ZeroToSt - ZeroToMat);
        mu1_N +=  ZeroToMat * gfac * alpha1 / Annuity;
        /* In the 1 factor case mu^2 is simply equal to mu1^2 */
        muFacSquare1 = Maths::square(mu1); // muFacSquare1 is scaled by alpha^2 (IR1F)
        muFacSquare2 = Maths::square(mu1_N); // muFacSquare2 is scaled by alpha^2 (IR1F)
        break;
    }
    case 2: {
        double alpha2 = ratesHJMUtil.getAlpha(1);
        /* calculate mu1 and mu2 respectively for the first and the
         * second factor*/
        double mu1 = 0.0; /* mu value for each factor (lognormal assumption) */
        double mu2 = 0.0; /* mu value for each factor (lognormal assumption) */
        for (int j=1;j <= NumCpns ;j++) {
            ratesHJMUtil.gFactor(currDate, CpnPayDates[j], gfac);
            mu1   += gfac[0] * weights[j] * alpha1;
            mu2   += gfac[1] * weights[j] * alpha2; 

        }/* for j*/
        double mu1_N = mu1 * ParYield;
        double mu2_N = mu2 * ParYield;
        ratesHJMUtil.gFactor(currDate, swapStart, gfac);
        mu1   -=  ZeroToSt *  gfac[0] * alpha1 / ( ZeroToSt - ZeroToMat);
        mu1_N -=  ZeroToSt *  gfac[0] * alpha1 / Annuity;
        
        mu2   -=  ZeroToSt *  gfac[1] * alpha2 / ( ZeroToSt - ZeroToMat);
        mu2_N -=  ZeroToSt *  gfac[1] * alpha2 / Annuity;

        ratesHJMUtil.gFactor(currDate, swapMat, gfac);
        mu1   +=  ZeroToMat * gfac[0] * alpha1 / ( ZeroToSt - ZeroToMat);
        mu1_N +=  ZeroToMat * gfac[0] * alpha1 / Annuity;

        mu2   +=  ZeroToMat * gfac[1] * alpha2 / ( ZeroToSt - ZeroToMat);
        mu2_N +=  ZeroToMat * gfac[1] * alpha2 / Annuity;

        /* mu^2 is so that mu^2 * dt = var (mu1 * dW1 + mu2 * dW2) */ 
        double rho1 = ratesHJMUtil.getRho(0);
        muFacSquare1 = Maths::square(mu1) + Maths::square(mu2) 
            + 2 * mu1   * mu2   * rho1;
        muFacSquare2 = Maths::square(mu1_N)+Maths::square(mu2_N) + 
            2 * mu1_N * mu2_N * rho1;
        break;
    }
    case 3: {

        /* factor weights */
        double alpha2 = ratesHJMUtil.getAlpha(1);
        double alpha3 = ratesHJMUtil.getAlpha(2);

        /* calculate mu1 and mu2 respectively for the first and the
         * second factor*/
        double mu1 = 0.0; /* mu value for each factor (lognormal assumption) */
        double mu2 = 0.0; /* mu value for each factor (lognormal assumption) */
        double mu3 = 0.0; /* mu value for each factor (lognormal assumption) */
        for (int j=1;j <= NumCpns;j++) {
            ratesHJMUtil.gFactor(currDate, CpnPayDates[j], gfac);
            mu1   += gfac[0] * weights[j] * alpha1;
            mu2   += gfac[1] * weights[j] * alpha2; 
            mu3   += gfac[2] * weights[j] * alpha3;

        }/* for j*/
        double mu1_N = mu1 * ParYield;
        double mu2_N = mu2 * ParYield;
        double mu3_N = mu3 * ParYield;

        ratesHJMUtil.gFactor(currDate, swapStart, gfac);
        mu1   -=  ZeroToSt *  gfac[0] * alpha1 / ( ZeroToSt - ZeroToMat);
        mu1_N -=  ZeroToSt *  gfac[0] * alpha1 / Annuity;

        mu2   -=  ZeroToSt *  gfac[1] * alpha2 / ( ZeroToSt - ZeroToMat);
        mu2_N -=  ZeroToSt *  gfac[1] * alpha2 / Annuity;

        mu3   -=  ZeroToSt *  gfac[2] * alpha3 / ( ZeroToSt - ZeroToMat);
        mu3_N -=  ZeroToSt *  gfac[2] * alpha3 / Annuity;

        ratesHJMUtil.gFactor(currDate, swapMat, gfac);
        mu1   +=  ZeroToMat * gfac[0] * alpha1 / ( ZeroToSt - ZeroToMat);
        mu1_N +=  ZeroToMat * gfac[0] * alpha1 / Annuity;

        mu2   +=  ZeroToMat * gfac[1] * alpha2 /( ZeroToSt - ZeroToMat);
        mu2_N +=  ZeroToMat * gfac[1] * alpha2 / Annuity;


        mu3   +=  ZeroToMat * gfac[2] * alpha3 / ( ZeroToSt - ZeroToMat);
        mu3_N +=  ZeroToMat * gfac[2] * alpha3 / Annuity;

        /* mu^2 is so that mu^2 * dt = var (mu1 * dW1 + mu2 * dW2 + mu3 * dW3)*/
        double rho1 = ratesHJMUtil.getRho(0);
        double rho2 = ratesHJMUtil.getRho(1);
        double rho3 = ratesHJMUtil.getRho(2);
 
        muFacSquare1 = Maths::square(mu1) + Maths::square(mu2) + 
            2 * mu1 * mu2 * rho1
            + Maths::square(mu3) + 2 * mu1 * mu3 * rho2 + 2 * mu2 * mu3 * rho3;
        muFacSquare2 = Maths::square(mu1_N) + Maths::square(mu2_N) + 
            2 * mu1_N * mu2_N * rho1
            + Maths::square(mu3_N) + 2 * mu1_N * mu3_N * rho2 + 
            2 * mu2_N * mu3_N * rho3 ;
        break;
    }
    }
}

/* Computes the deterministic squared mu factors (see section 2.2).
   From swapvol::FindTau */
DateTime SRMSwaption::findTau(const SRMRatesHJMUtil&  ratesHJMUtil,
                              const DateTime&         swapStart,
                              const DateTime&         swapMat,
                              const DateTime&         today) const{
    /* First Guess for tau : halfway between SwapSt and SwapMat*/
    DateTime Maturity = today.rollDate((int)(0.5 * (swapStart.daysDiff(today)
                                                    +swapMat.daysDiff(today))));

    /* Start Newton-Raphson*/
    int     iter = 0; /* iterations number in Newton-Raphson     */
    int     Num  = 0; /* number of consecutive iterations without change in mu*/
    double  temp = 0.0;         /* local variable for Xi */
    double  Xi;    /* xi factor      */
    do {
        double  derivative; /*  xi factor derivative    */
        double  maxXi;      /* if Xi << maxXi, the minimization succeeds  */
        ratesHJMUtil.xiFactor(Xi,          /* (O) : xi factor to cancel */
                         derivative,/*(O): derivative of xi 
                                      with respect to maturity */
                         maxXi,
                         weights,
                         CpnPayDates, NumCpns,
                         Maturity, today);

        if (fabs(derivative) < SRMConstants::SRM_TINY) {
            Maturity = Maturity.rollDate(TAU_MONTH);  /* avoid explosion in the
                                                     Newton-Raphson*/
        } else {
            /* General Newton-Raphson equation to cancel xi  */ 
            Maturity = Maturity.rollDate(-(int)(DateTime::DAYS_PER_YEAR *
                                                Xi / derivative)); 
        }

        Maturity = swapStart.max(Maturity); /* tau must always be > to the 
                                            swaption expiry	 */
        Maturity = swapMat.min(Maturity); /* tau must always be < to the 
                                             swaption maturity */
        iter++;  /* Number of N-R iterations */

        if (Maths::areEqualWithinTol(temp,Xi, SRMConstants::SRM_TINY)) {
            Num++;  /* cumulative number of iterations with no change for Xi */
        } else {
            Num = 0;
        }
        temp = Xi;
    }
    while ((fabs(Xi) > SRMConstants::SRM_TINY)  &&  (iter < SRMConstants::MAX2QITER)  &&
            (Num < 5));  /*condition for the Newton-Raphson*/

    /* Transfer values to output */
    return Maturity;
}

/* Computes the Zero Coupon weights, the annuity, the ParYield, 
 * the Coupon payment dates and the DCCFrac for a given swaption.
 From swapvol::SwapDetails */ 
void SRMSwaption::swapDetails(const SRMRatesHJMUtil&  ratesHJMUtil,
                              const DateTime&         swapStart,
                              const DateTime&         swapMat,
                              bool                    stubAtEnd,
                              bool                    useDiffusedCurve) {
    static const string method("SRMSwaption::swapDetails");
    MaturityPeriodSP fixedPayInterval(ratesHJMUtil.getSwapFrequency());
    string payInterval;
    int    payFrequency;
    // split eg 3M into a '3' and a 'M'
    fixedPayInterval->decompose(payFrequency, payInterval);
    // This set of dates includes the swap start
    DateTimeArraySP dates(SwapTool::dateArray(swapStart, swapMat,
                                              payFrequency, payInterval,
                                              stubAtEnd));
    CpnPayDates = *dates;
    DayCountConventionSP dcc(ratesHJMUtil.getSwapDCC());
    NumCpns =  CpnPayDates.size()-1; /* Number of Coupon payment dates 
                                        (-1 because of swap start date) */
    weights.resize(NumCpns+1); // weight[0] not used I think
    double AnnuityL = 0.0; /* Initialization      */
    /* Get coupon payment dates and DCC Frac */
    int j;
    for (j = 1; j < CpnPayDates.size(); j++) {
        double DCCFrac = dcc->years(CpnPayDates[j-1], CpnPayDates[j]);
        /* Zero to current coupon payment date */
        double ZerotoCpn = ratesHJMUtil.pv(CpnPayDates[j], useDiffusedCurve);
        weights[j] = DCCFrac * ZerotoCpn;
        AnnuityL  += DCCFrac * ZerotoCpn;
    }/* for j */

    /* Compute weights : divide by the annuity*/
    for (j = 1; j < CpnPayDates.size(); j++) {
        weights[j] /= AnnuityL;
    }/*for j*/

    /* Zero to Swap expiry */
    double ZeroToSt = ratesHJMUtil.pv(swapStart, useDiffusedCurve);
    /* Zero to Swap maturity */
    double ZeroToMat = ratesHJMUtil.pv(swapMat, useDiffusedCurve);
    /* Transfer values to output */
    Annuity  = AnnuityL; /* annuity at time 0               */
    ParYield = (ZeroToSt - ZeroToMat)/ AnnuityL; /* par yield at time 0 */
}

/*  Populate the new structure created to avoid recalculating the factors 
   (mu, tau, variance, v(t)) that do not depend on the spot vol at each 
   step of the Newton-Raphson algorithm in the main bootstrapping routine
   see Appendix C of report. From swapvol::PopulateSwpn */
void SRMSwaption::populate(const SRMRatesHJMUtil&  ratesHJMUtil,
                           const DateTime&         swaptionExpiry,
                           const DateTime&         swapStart,
                           const DateTime&         swapMat,
                           int                     prevIdx,
                           const DateTimeArray&    dates){
    const DateTime& today = ratesHJMUtil.getBaseDate();
    this->swapStart = swapStart;
    this->swapMat = swapMat;
    this->swapExp = swaptionExpiry;
    this->prevIdx = prevIdx; // prevIdx is idx of the previous swaption
    idx = swaptionExpiry.findLower(dates);
    expiry = today.yearFrac(swapStart);
    if (Maths::isZero(expiry)){
        throw ModelException("SRMSwaption::populate",
                             "Swap start must be > today");
    }
    // calculate NumCpns, Annuity, CpnPayDates, ParYield, weights
    swapDetails(ratesHJMUtil, swapStart, swapMat,
                false /* stub at front */, true /* use diffused curve */);
    /*  Find tau : the maturity of the vol rate driver */
    tau = findTau(ratesHJMUtil, swapStart, swapMat, today); // independent of alpha (IR1F)
    /*  Populate the variance array and the exp_ integrals arrays */ 
    variance.resize(idx);
    rbart.resize(idx);
    del_t.resize(idx);
    muFacSquare1.resize(idx);
    muFacSquare2.resize(idx);
    mr.resize(idx);
    for (int i = 0; i< idx; i++) {
        double delt_1 = dates[i].yearFrac(tau);
        double delt_2 = dates[i+1].yearFrac(tau);
        /* value of rbar at time t */
        rbart[i] = ratesHJMUtil.rBar(dates[i]);
        /* time step for the integration*/
        del_t[i] =  dates[i].yearFrac(dates[i+1]);        
        variance[i] = ratesHJMUtil.modelVariance(delt_1, delt_2); // scaled by alpha^2 (IR1F)
        ratesHJMUtil.modelMeanReversionIntegral(delt_1, delt_2, mr[i]); // scaled by alpha (IR1F)
        /* compute the Squared mu factors */
        muFactor(ratesHJMUtil, dates[i], swapStart, swapMat,
                 muFacSquare1[i], muFacSquare2[i]); // mu is scaled by alpha (IR1F)
    }
}

/* Compute the expectation of (F_t)^2 for the 2Q case see section 3.2 
   From swapvol::MeanFSquare */
double SRMSwaption::meanFSquare(
    double t,    /* (I) time t : cumulative variance v(t) (see report)*/
    double qR,   /* (I) QRight                                        */
    double qL){  /* (I) QLeft                                         */
    /*  To avoid explosion if t is too small */
    if (fabs(t) < SRMConstants::SRM_TINY) {
        return 1.0;
    }
    /*first compute the expectation when F_t is > 0 */
    /* new parameters as introduced by Karatzas & Shreve */
    double theta1 = - qL / 2.0;
    double theta0 = - qR / 2.0;
    /* time function for the integration */
    double phi = exp( 1 / 4.0 * (theta1 * theta1 - theta0 * theta0) *
                      t + 1 /8.0 * (Maths::square(theta0 - theta1)) * t
                      - 1 / 2.0 * theta1 * theta1 * t);
    /* Equivalent theta */
    double theta = - (theta0 - theta1) / 2.0;
    /* New q for the integration */
    double q = 2 * qR + theta0 - theta;
    /* expression of the expectation of f^2 when F>0 */
    double Expec = phi * ((1 + theta / q) * exp(q * theta * t + q * 
                                                q * t /2.0) * 
							SRMBlack::NormalH((theta + q) * sqrt(t))
                          - theta / q * SRMBlack::NormalH(theta * sqrt(t)) );
    /* add expectation when qL prevails (symmetry properties of brownians)*/
    /* symmetry properties : see report */
    theta0 = - theta1;
    theta1 = -theta0;
    theta = - (theta0 - theta1) / 2.0;

    /* time function for the integration */
    phi = exp( 1 / 4.0 * (theta1 * theta1 - theta0 * theta0) * t + 
               1 /8.0 * (Maths::square(theta0 - theta1)) * t
               - 1 / 2.0 * theta1 * theta1 * t );
    /* New q for the integration */
    q = -2 * qL + theta0 - theta;
    /* expression of the expectation of f^2 when F<0 */
    Expec += phi * ((1 + theta / q) * exp(q * theta * t + q *q * t /2.0) *
                    SRMBlack::NormalH((theta + q) * sqrt(t))
                    - theta / q * SRMBlack::NormalH(theta * sqrt(t)) );	
    return (Expec);
}

/* Compute the correlation between mu factor square and F_t square
   (see section 3.3). From swapvol::CorrelMuF */
void SRMSwaption::correlMuF(
    const SRMRatesHJMUtil&        ratesHJMUtil,
    const DateTimeArray&          mergedList, // zero dates + vol dates
    const DateTime&               swapStart,
    const DateTime&               swapMat,
    const vector<double>&         extendedSpotVol, /* : descretized spotVol */
    vector<double>&               correl) const{ //(O)correlation between mu^2 and F^2 

    double qL = ratesHJMUtil.getQLeft();
    double qR = ratesHJMUtil.getQRight();

    if (fabs(qR) < QSMALL  && fabs(qL) < QSMALL) {
        for (int i = 0; i< idx; i++) {
            correl[i] = 1.0;
        }
        return;
    }

    /* semi-empirical formula   */
    double correlFac = fabs(qL) * fabs(qL)/
        (fabs(qL) + fabs(qR)) + fabs(qR) * fabs(qR)/(fabs(qL) + fabs(qR));
    double q         = correlFac;

    /* Zero to expiry */
    double ZeroToSt = ratesHJMUtil.pv(swapStart, true);
    /* Zero to maturity */
    double ZeroToMat = ratesHJMUtil.pv(swapMat,true);
    double zeroDiff = ZeroToSt - ZeroToMat;
    double zeroDiffSq = Maths::square(zeroDiff);
    double sum_cov = 0.0; /* cumulative covariance between F_t and mu_t */
    int numFactors = ratesHJMUtil.numFactors();
    switch (numFactors) {
    case 1:{
        /* factor weight */
        double alpha = ratesHJMUtil.getAlpha(0);
        /* calculate the correlation at each time t of the time line
         * MERGEDLIST */
        for (int i = 0; i < idx; i++) {
            /* note: role of alpha in IR1F
                idx=0 (first swaption) => correl=1 (no dependency on alpha)
                idx>0 (2nd to nth swaption) => see comments below */
            /* main expression for the correlation term   */
            correl[i] =  exp(correlFac * sum_cov);
            /* calculate the cumulative covariance */
            double gfac1 = ratesHJMUtil.gFactor(mergedList[i], swapStart, 0);
            double gfac2 = ratesHJMUtil.gFactor(mergedList[i], swapMat, 0);
            double temp = -Maths::square(alpha) * ZeroToSt * ZeroToMat * 
                Maths::square(gfac1 - gfac2) / 
                (zeroDiffSq *  sqrt(muFacSquare1[i]) );
            gfac1 = ratesHJMUtil.gFactor(mergedList[i],tau, 0);
            for (int j = 1; j <= NumCpns; j++) {
                double gfac2 = ratesHJMUtil.gFactor(mergedList[i], 
                                               CpnPayDates[j], 0);
                temp += Maths::square(alpha) * weights[j] *
                    (gfac1 - gfac2) * gfac2 / sqrt(muFacSquare1[i]) ;
            }
            temp *=  extendedSpotVol[i] * rbart[i];
            // sigma_mu is scaled by alpha and extendedSpotVol is scaled by 1/alpha => cancellation
            // sigma_F (via mr) is scaled by alpha and extendedSpotVol is scaled by 1/alpha => cancellation
            // thus: sum_cov does not depend on alpha in IR1F
            sum_cov += temp * q * extendedSpotVol[i] * mr[i][0] ;
        }
        break;
    }
    case 2:{
        // reserve some space
        vector<double> alpha(numFactors);
        vector<double> gfac1(numFactors);
        vector<double> gfac2(numFactors);
        double derivative[2][2];
        /* factor weights */
        ratesHJMUtil.getAlpha(alpha);

        /* correlation */
        double rho1 = ratesHJMUtil.getRho(0);

        /* calculate the correlation at each time t of the time line 
           MERGEDLIST */
        for (int i = 0; i < idx; i++) {
            /* main expression for the correlation term   */
            correl[i] =  exp(correlFac * sum_cov);

            /* we calculate the mu factor and the derivative of mu[p] with 
               respect to the kth factor : derivative[p][k] */
            ratesHJMUtil.gFactor(mergedList[i], swapStart,gfac1);
            ratesHJMUtil.gFactor(mergedList[i], swapMat,gfac2);

            double mu1 = alpha[0] * (- gfac1[0] * ZeroToSt + gfac2[0] *
                                     ZeroToMat) / zeroDiff;
            double mu2 = alpha[1] * (- gfac1[1] * ZeroToSt + gfac2[1] *
                                     ZeroToMat) / zeroDiff;
            for (int k = 0; k < 2; k++) {
                for (int p = 0; p < 2; p++) {
                    /* derivative of mu_p with respect to the kth factor */
                    derivative[p][k] = alpha[p] * alpha[k] *( 
                        ( gfac1[p] * gfac1[k] * ZeroToSt - 
                          gfac2[p] * gfac2[k] * ZeroToMat) / zeroDiff
                        + gfac1[p] * ZeroToSt * (gfac2[k] * ZeroToMat -
                                                 gfac1[k] * ZeroToSt) / 
                        zeroDiffSq 
                        - gfac2[p] * ZeroToMat* (gfac2[k] * ZeroToMat -
                                                 gfac1[k] * ZeroToSt) /
                        zeroDiffSq);
                }/* for p*/
            }/* for k*/

            ratesHJMUtil.gFactor(mergedList[i], tau, gfac1);
            for (int j = 1; j <= NumCpns; j++) {
                ratesHJMUtil.gFactor(mergedList[i], CpnPayDates[j], gfac2);
                mu1 += gfac2[0] * weights[j] * alpha[0];
                mu2 += gfac2[1] * weights[j] * alpha[1];
                
                for (int k= 0; k < 2; k++) {
                    for (int p = 0; p < 2; p++) {
                        derivative[p][k] += weights[j] * gfac2[p] *
                            (gfac1[k] - gfac2[k]) * alpha[p] * alpha[k];
                    }/* for p*/
                }/* for k*/

            }/* for j*/
            const vector<double> mr = this->mr[i]; // for ease
            /* we calculate the covariance term for each factor */
            double cov1 = derivative[0][0] * 
                mr[0] + rho1 * derivative[0][0] * mr[1] +
                derivative[0][1] * mr[1] + rho1 * derivative[0][1]*mr[0];

            double cov2 = derivative[1][0] * mr[0] + 
                rho1 * derivative[1][0] * 
                mr[1] + derivative[1][1] * mr[1] + 
                rho1 * derivative[1][1]*mr[0];

            double temp = mu1 * cov1 + mu2 * cov2 + rho1 * mu1 * cov2 + 
                rho1 * mu2 * cov1;
            temp *=  extendedSpotVol[i] * rbart[i] / muFacSquare1[i];
            sum_cov += temp * q * extendedSpotVol[i]  ; 

        }/* for i*/
        break;
    }
    case 3:{
        // reserve some space
        vector<double> alpha(numFactors);
        vector<double> gfac1(numFactors);
        vector<double> gfac2(numFactors);
        double derivative[3][3];
        /* factor weights */
        ratesHJMUtil.getAlpha(alpha);

        /* correlation */
        double rho1 = ratesHJMUtil.getRho(0);
        double rho2 = ratesHJMUtil.getRho(1);
        double rho3 = ratesHJMUtil.getRho(2);

        /* factor weights */
        ratesHJMUtil.getAlpha(alpha);
        
        for (int i=0; i < idx; i++) {
            /* main expression for the correlation term  */
            correl[i] =  exp(correlFac *sum_cov);

            /* we calculate the mu factors and the derivative of mu[p]
               with respect to the kth factor : derivative[p][k] */ 
            ratesHJMUtil.gFactor(mergedList[i], swapStart,gfac1);
            ratesHJMUtil.gFactor(mergedList[i], swapMat,gfac2);
            double mu1 = alpha[0] * ( - gfac1[0] * ZeroToSt + 
                                      gfac2[0] * ZeroToMat) / zeroDiff ;
            double mu2 = alpha[1] * ( - gfac1[1] * ZeroToSt + 
                                      gfac2[1] * ZeroToMat) / zeroDiff;
            double mu3 = alpha[2] * ( - gfac1[2] * ZeroToSt + 
                                      gfac2[2] * ZeroToMat) / zeroDiff;
            for (int k = 0; k < 3; k++) {
                for (int p = 0; p < 3; p++) {
                    /* derivative of mu_p with respect to the kth factor */
                    derivative[p][k] = alpha[p] * alpha[k] * (
                        ( gfac1[p] * gfac1[k] * ZeroToSt - 
                          gfac2[p] * gfac2[k] * ZeroToMat) / zeroDiff
                        + gfac1[p] * ZeroToSt * (gfac2[k] * ZeroToMat - 
                                                 gfac1[k] * ZeroToSt) / 
                        zeroDiffSq 
                        - gfac2[p] * ZeroToMat * (gfac2[k] * ZeroToMat - 
                                                  gfac1[k] * ZeroToSt) /
                        zeroDiffSq);
                }/* for p */
            }/* for k*/
            
            ratesHJMUtil.gFactor(mergedList[i], tau, gfac1);
            for (int j = 1; j <= NumCpns; j++) {
                ratesHJMUtil.gFactor(mergedList[i], CpnPayDates[j], gfac2);
                mu1 += gfac2[0] * weights[j] * alpha[0];
                mu2 += gfac2[1] * weights[j] * alpha[1];
                mu3 += gfac2[2] * weights[j] * alpha[2];
                for (int k = 0; k < 3; k++) {
                    for (int p = 0; p < 3; p++) {
                        derivative[p][k] += weights[j] * gfac2[p] *
                            (gfac1[k] - gfac2[k]) * alpha[p] * alpha[k];
                    }/* for p*/
                }/* for k*/
            }/* for j*/
            const vector<double> mr = this->mr[i]; // for ease
            /* we calculate the covariance term for each factor */ 
            double cov1 = derivative[0][0] * (mr[0] + rho1 * mr[1] + 
                                              rho2 * mr[2])
                + derivative[0][1] * (rho1 * mr[0] + mr[1] + 
                                      rho3 * mr[2]) 
                + derivative[0][2] * (rho2 * mr[0] + rho3 * mr[1] + 
                                      mr[2]);

            double cov2 =  derivative[1][0] * (mr[0] + rho1 * mr[1] +
                                        rho2 * mr[2])
                + derivative[1][1] * (rho1 * mr[0] + mr[1] + 
                                      rho3 * mr[2]) 
                + derivative[1][2] * (rho2 * mr[0] + rho3 * mr[1] + 
                                      mr[2]);

            double cov3 =  derivative[2][0] * (mr[0] + rho1 * mr[1] + 
                                               rho2 * mr[2])
                  + derivative[2][1] * (rho1 * mr[0] + mr[1] + 
                                        rho3 * mr[2]) 
                  + derivative[2][2] * (rho2 * mr[0] + rho3 * mr[1] + 
                                        mr[2]);
            
            double temp =  mu1 * cov1 + mu2 * cov2 + mu3 * cov3 
                + mu2 * rho1 * cov1 + mu1 * rho1 * cov2
                + mu1 * rho2 * cov3 + mu3 * rho2 * cov1
                + mu2 * rho3 * cov3 + mu3 * rho3 * cov2;
            
            temp *=  extendedSpotVol[i] * rbart[i] / muFacSquare1[i];
            sum_cov += temp * q * extendedSpotVol[i];
        }/* for i*/
        break;
    }}
}

/** Swaption implied volatility  see Appendix C. populate must have
    been called before this function */
double SRMSwaption::impVol2Q( // from SRM3::SwpnVol2Q
    const SRMRatesHJMUtil&        ratesHJMUtil,
    const DateTimeArray&          mergedList, // zero dates + vol dates
    double                        sigma,      /* spotvol prevailing between
                                              expiry(i-1) and expiry(i) */
    vector<double>&               extendedSpotVol) const{ /* : descretized spotVol */
    static const                  string method("SRMSwaption::impVol2Q");
    try {

        /* for convenience */
        double qR = ratesHJMUtil.getQRight();
        double qL = ratesHJMUtil.getQLeft();
        double q = 0.5 * fabs(qL) + 0.5 * fabs(qR);
        
        /* if only one q is 0, to avoid explosion in the expression of E(F_t^2)  */
        if (fabs(qL) < QSMALL && fabs(qR) > QSMALL){
            qL = QSMALL;
        }
        if (fabs(qR) < QSMALL && fabs(qL) > QSMALL){
            qR = QSMALL;
        }
        double volL = 0.0; /* Lognormal and nomal approx. vol  */  
        double valN = 0.0; /* value assuming normal distributed yield */
        double var1 = 0.0; /* cumulative variance for 1q case  */  
        double var2 = 0.0; /* cumulative variance */
        /* prevIdx is the index of the previous swaption:
           the new sigma prevails between the expiry of the previous
           swaption and the expiry of the swaption we are calibrating */ 
        int currIdx = prevIdx;
        int i;
        for (i = currIdx; i< idx; i++) {
            extendedSpotVol[i] = sigma; 
        }/* for i*/
        /* get correlations between F^2 and mu^2 */
        vector<double> correl(idx);
        correlMuF(ratesHJMUtil, mergedList, swapStart, swapMat, extendedSpotVol, correl);
        /* calculate the implied vol as  an integral  */
        double MeanSigmaRSquare;
        for (i = 0; i< idx; i++) {
            if (Maths::areEqualWithinTol(qR, 0.0, SRMConstants::SRM_TINY) && 
                Maths::areEqualWithinTol(qL,0.0, SRMConstants::SRM_TINY)){  /*  Normal case*/
                MeanSigmaRSquare = Maths::square(extendedSpotVol[i] * rbart[i]);
            } else {
                /* 1q case treated separately */
                if (Maths::areEqualWithinTol(qR, qL, SRMConstants::SRM_TINY)){
                    MeanSigmaRSquare = Maths::square(
                        extendedSpotVol[i] * rbart[i])* exp(var1); 
                } else { /* 2q case */
                    MeanSigmaRSquare = Maths::square(extendedSpotVol[i] * rbart[i]) 
                        * meanFSquare(var2,qR,qL);
                }
            }
            /* cumulative variance */
            volL    += muFacSquare1[i] * correl[i]*   MeanSigmaRSquare  *  del_t[i];
            valN    += muFacSquare2[i] *  MeanSigmaRSquare  *  del_t[i];
            if (!Maths::finite(volL) || !Maths::finite(valN)){
                throw ModelException("SRMSwaption::impVol2Q",
                                     "Variance has become infinite at "
                                     "index "+Format::toString(i));
            }
            /* single q variance */
            var1 += Maths::square(qR * extendedSpotVol[i]) * variance[i];
            /* 2q variance */
            var2 += Maths::square(extendedSpotVol[i]) *  variance[i];
        }/* for i*/
        volL = volL/expiry;
        volL = sqrt(volL);
        /* Swaption price when normal approximation */
        valN  /=  (2.0 * Maths::PI);
        valN   =  sqrt(valN);
        double volN = SRMBlack::ImpVol_BS2Q(ParYield, ParYield, expiry , 
                                   valN, 'C' , 1.0,1.0,0.0, 0.5);
        /* If ImpVol_BS2Q fails, then volN = volL */
        /* note that ImpVol_BS2Q returns -10e5 if failed */
        if (volN < 0.0) {
            volN = volL;
        }
        /* note: both, volN and volL are not dependent on alpha (modulo numerical diffs), 
            since extendedSpotVol is scaled by 1/alpha and thus the dependency on alpha 
            via eg muFacSquare or variance cancels out */
        return (q < 1.0? 
                (1.0 - q) * volN + q * volL: /*Interpolation between lognormal and
                                               normal vol */
                volL);
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/* Given the values of tau at each vo date, populate the array ir->tau
   whose size is the number of zero dates. From swapvol::ExtendTau */
DateTimeArray SRMSwaption::extendTau(
    const SRMRatesHJMUtil&       ratesHJMUtil,
    const DateTimeArray& tau){ /* (I) : tau value at each vol dates  */
    const DateTimeArray& ZDate = ratesHJMUtil.getExtendedTimeLine();
    const DateTimeArray& VolDate = ratesHJMUtil.getSwaptionExpiryDates();
    DateTimeArray extendedTau(ZDate.size());
    int NumSwap = 0;       /* Vol date  index */ 
    for (int i = 0; i < ZDate.size(); i++) {
        /* if zero date > swaption expiry, then increase the vol date index */
        /* except when this happens for last swaption                       */
        if(ZDate[i] >= VolDate[NumSwap]) {
            if (NumSwap == VolDate.size()-1){   /* last swaption */
                extendedTau[i] = tau[NumSwap];
            } else {
                NumSwap ++; /* increment the swaption index */
                extendedTau[i] = tau[NumSwap]; /* fill in ir->tau */
            }
        } else {
            /* zero date is < than the vol date and there is no
               need to increase the vol date index */
            extendedTau[i] = tau[NumSwap];
        }
    }
    return extendedTau;
}

/** Constructor using explicit data */
SRMSwaption::SRMSwaption(const SRMRatesHJMUtil&    ratesHJMUtil,
                         const DateTime&           swaptionExpiry,
                         const DateTime&           swapStart,
                         const DateTime&           swapMat,
                         const DateTime&           prevSwaptionExpiry,
                         const DateTimeArray&      dates)
{
    int prevIdx = prevSwaptionExpiry.findLower(dates);
    // to do: how much of the work done in populate is actually needed
    // for pricing/recalibrating against the swaption
    populate(ratesHJMUtil, swaptionExpiry, swapStart, swapMat, prevIdx, dates);
}
   
/** Implied vol of a call(?) swaption. Returns < 0 if implied vol cannot
    be found */
double SRMSwaption::impliedVol(double swaptionPrice) const{
    // can surely do something better for the initial guess (eg what the
    // current vol is)
    return SRMBlack::ImpVol_BS2Q (ParYield, ParYield, expiry , swaptionPrice/ Annuity , 
                        'C' ,1.0,1.0,0.0, 0.50 /* initial guess */);
}

/* Calibrates the Swaption to the target implied vol target  
 * From swapvol::CalibSwaption */
void SRMSwaption::calibrate(const SRMRatesHJMUtil&  ratesHJMUtil,
                            double          target,
                            vector<double>& SpotVol /* (M) */ ) const{

    double sigma_1 = SpotVol[prevIdx];
    double sigma   = sigma_1;
    const DateTimeArray& extendedTimeLine = ratesHJMUtil.getExtendedTimeLine();

    SRMSwaptionFunctor swapFunc(*this, ratesHJMUtil, extendedTimeLine, SpotVol, target);
    NewtonRootFinder<SRMSwaptionFunctor> rootFinder(true);

    double low = 0.0;
    double high = 10.0; //Root finder will extend the right hand side until f becomes positive.

    try
    {
        rootFinder.findRoot(low, high, sigma_1, VOLPRECISION, swapFunc, sigma_1, true);
        QLIB_VERIFY(sigma_1 >= 0.0, "Negative sigma_1"); //should never happen: root solver will already throw an exception

        // if this is the last swaption (irrespective of last sim date)
        // then extend the spot vol. Currently, I can't see the point in doing
        // this as the SpotVol is not used beyond the last simulation date
        // but perhaps SpotVol will be used somewhere else .... (which I 
        // haven't got too it)
        if (ratesHJMUtil.getLastExpiry().equals(swapExp)){
            double sigma = SpotVol[prevIdx];
            for (unsigned int i = idx; i < SpotVol.size(); i++){
                SpotVol[i] = sigma; 
            }
        }
    }
    catch (exception &/*e*/)
    {
        /* if number of iterations exceeds 20 , the NR has failed:
        rollback to the previous spotvol */
        for (int i= prevIdx; i < idx; i++) {
            SpotVol[i] = sigma; 
        }/* for i*/
    }
}

/** Construct a Pricer which can be used to price paths for this swaption */
SRMSwaptionPricer* SRMSwaption::createPricer(const SRMRatesHJMUtil& ratesHJMUtil) const{
    IYieldCurveConstSP diffCurve(ratesHJMUtil.getDiffYC());
    IYieldCurveConstSP stochDiscCurve(ratesHJMUtil.getDiscYC());
    // these casts are a bit weak - should be IYieldCurve 
    YieldCurveConstSP discCurve(YieldCurveConstSP::dynamicCast(stochDiscCurve));
    DayCountConventionSP dcc(ratesHJMUtil.getSwapDCC());
    // note: CpnPayDates includes swap start date
    SVGenIRSwapSP swapGen(new SVGenIRSwap(YieldCurveConstSP::dynamicCast(diffCurve),
                                    discCurve,
                                    swapStart,
                                    DateTimeArray(CpnPayDates.begin()+1, 
                                                  CpnPayDates.end()),
                                    dcc));
    SVGenDiscFactorSP dfGen(new SVGenDiscFactor(ratesHJMUtil.getBaseDate(),
                                          discCurve,
                                          swapStart));
    return new SRMSwaptionPricer(ParYield, swapGen, dfGen);
}


/* Main calibration routine. From swapvol::CalibVol2Q */
DateTimeArray SRMSwaption::calibVol2Q(const SRMRatesHJMUtil&  ratesHJMUtil,
                                      bool                    skipFlag,
                                      DoubleArray&            SpotVol /* (O) */ )
{
    static const string method("SRMHJMCalib::calibVol2Q");
    try {
        const DateTimeArray& extendedTimeLine = ratesHJMUtil.getExtendedTimeLine();
        const DateTimeArray& VolDate = ratesHJMUtil.getSwaptionExpiryDates();
        const DoubleArray& Vol = ratesHJMUtil.getSwaptionVols();
        /* Compute the mergedList : zero dates + vol dates */
        DateTimeArray MergedList(DateTime::merge(extendedTimeLine, VolDate));

        /*  Initialise the swaption structure  and allocate memory */
        SRMSwaption SwpnStruc;
        vector<double> ExtendedSpotVol(MergedList.size());

        // get hold of the details of underlying swap of the swaption vols
        const DateTimeArray& SwapSt = ratesHJMUtil.getSwapStartDates();
        const DateTimeArray& SwapMat = ratesHJMUtil.getSwapMatDates();
        DateTimeArray tau(SwapSt.size());
        SpotVol.resize(SwapSt.size());
        
        //first check for zero vols.....
        bool allZero = true;
        for (int p = 0;p < SwapSt.size(); p++) {
            if (Vol[p] != 0.0) {
                allZero = false;
            }
        }
        if (allZero) {
            for (int p = 0;p < SwapSt.size(); p++) {

                /* populate SwpnStruc */ 
                SwpnStruc.populate(ratesHJMUtil, VolDate[p], SwapSt[p], SwapMat[p],
                    SwpnStruc.idx, MergedList);

                /* store the value of tau in the array tau */
                tau[p]           = SwpnStruc.tau;

                SpotVol[p] = 0.0;
            }
        }
        else
        {

            //Not all vols are zero: proceed to usual calibration
            
            /* first guess for the spot volatility of the first swaption */
            double sigma_1 = 1.0;
            /* main calibration loop */
            for (int p = 0; p < SwapSt.size(); p++) {
                /* i is the number of iterations in the Newton-Raphson */
                int iter = 0;

                /* populate SwpnStruc */ 
                SwpnStruc.populate(ratesHJMUtil, VolDate[p], SwapSt[p], SwapMat[p],
                    SwpnStruc.idx, MergedList);

                /* store the value of tau in the array tau */
                tau[p]           = SwpnStruc.tau;
                /* start Newton-Raphson */

                double low = 0.0;
                double high = 10.0; //Root finder will extend the right hand side until f becomes positive 

                SRMSwaptionFunctor swapFunc(SwpnStruc, ratesHJMUtil, MergedList, ExtendedSpotVol, Vol[p]);
                NewtonRootFinder<SRMSwaptionFunctor> rootFinder;            
                try
                {
                    rootFinder.findRoot(low, high, sigma_1, VOLPRECISION, swapFunc, sigma_1, true);
                    SpotVol[p] = sigma_1;
                }
                catch (exception &e)
                {
                    if (!skipFlag) {
                        throw ModelException(e, method, "Unable to bootstrap vol "
                                                        "(Newton-Raphson has failed)!");
                    }
                    /* MAR: This is change to native SRM3. Revert back to original
                    guess as otherwise you end up in some 'random' platform-
                    dependent place */
                    if (p == 0)
                    {
                        if (Vol[p])
                            SpotVol[p] = Vol[p];
                        else
                            SpotVol[p] = SpotVol[p-1]; 
                    }
                }

                /* if the bootstrapped spot vol is negative, the NR had failed */
                if (sigma_1 < 0.0 && !skipFlag) {
                    throw ModelException("negative vol failed at (" + 
                        Format::toString(p) + ")");
                }
            }
        }

        /* finally, populate extended tau array with the values of taus
        contained in the array tau */
        return (extendTau(ratesHJMUtil, tau));
    } catch (exception& e){
        throw ModelException(e, method);
    }
} /* CalibVol2Q */
        
/*****  From swapvol.c::SpotVol   ****************************************/
/* MAR: This is the old calibration method - used currently for quanto 
   adjustment to CDS Par Spreads */
/*                                                                           */
/*      Calculate the spot vol curve from a series of base/swaption vols     */
/*                                                                           */
/*      The curve has the same points as the diffused zero curve which is a  */
/*      superset of the timeline of the simulation.                          */
/*                                                                           */
/*      The vols stored in ir->SpotVol[i] are unadjusted for the forward shift*/
/*      effect. i.e. they are applied to an effective rate which Vladimirises*/
/*      to rbar.                                                             */
/*      It is necessary to store these unadjusted vols so that the spot vols */
/*      are compatible with other calibration functions e.g. fx spot vol     */
/*                                                                           */
DoubleArray SRMSwaption::spotVolVF(const SRMRatesHJMUtil&  ratesHJMUtil,
                                   bool                    skipFlag) 
{
    static const string method("SRMSHJMCalib::spotVolVF");
    /* only used for quanto adjustment for cds par spread curves - so only
       bothered porting the single factor case */
    if (ratesHJMUtil.numFactors() > 1){
        throw ModelException(method, "Only single factor case supported");
    }
    double aw = ratesHJMUtil.getAlpha(0); // Historical aweights
    const DoubleArray& swaptionVols = ratesHJMUtil.getSwaptionVols();
    const DateTimeArray& swaptionExpiries = ratesHJMUtil.getSwaptionExpiryDates();
    const DateTimeArray& swapStartDates = ratesHJMUtil.getSwapStartDates();
    const DateTimeArray& swapMatDates = ratesHJMUtil.getSwapMatDates();
    MaturityPeriodSP swapFrequency = ratesHJMUtil.getSwapFrequency();
    DayCountConventionSP swapDCC = ratesHJMUtil.getSwapDCC();
    int NbVol = swaptionVols.size();
    DoubleArray   Aweight(NbVol); // return value
    vector<vector<double> > B(NbVol); // B in Christian memo
    vector<double>          VolT(NbVol); //  Expiries in years
    vector<double>          pY(NbVol); //  Fwd par yield
    double beta = ratesHJMUtil.getBeta(0);
    double qLeft = ratesHJMUtil.getQLeft();
    double qRight = ratesHJMUtil.getQRight();
    double fwdShift = ratesHJMUtil.getFwdShift();
    const DateTime& baseDate = ratesHJMUtil.getBaseDate();
    int i; // MSVC broken
    for (i = 0; i < NbVol; i++) {
        /* Time to expiry in years */
        VolT[i] = baseDate.yearFrac(swaptionExpiries[i]);
        B[i] = ratesHJMUtil.bFactor(i); // B Factor this swaption with index i
        B[i][0] *= exp(-beta * swaptionExpiries[i].yearFrac(swapStartDates[i]));
        pY[i] = ratesHJMUtil.getDiffYC()->couponRate(swapStartDates[i],
                                                swapMatDates[i], 
                                                *swapFrequency,
                                                false /* stub At front */,
                                                swapDCC.get());
    }
    double L = 0.0; // Integrals of lambda factors
    int LastCalibIdx = -1;              /* no calibrated points so far */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */

    for (i = 0; i < NbVol; i++) {       
        double T = VolT[i]; //  Time to current expiry
        // calc time between two consecutive expiries
        double t = ((LastCalibIdx == -1) ? 
                    VolT[i] : (VolT[i] - VolT[LastCalibIdx]));
        
        double D = aw * B[i][0]; //  aweight*B
        /* Lognormal to X-space vol adjustment. The single and */
        /* two q cases are treated differently for consistency */
        /* with old model.                                     */
        double y; // Vector for lambda system;
        if ((fabs(qLeft - qRight) < SRMConstants::SRM_TINY) && 
            (fabs(fwdShift) < SRMConstants::SRM_TINY)) {
            if (fabs(qLeft) > QCUTOFF){
				y = SRMBlack::Normal_InvH(
					  qLeft * (SRMBlack::NormalH(0.5*sqrt(T)*swaptionVols[i])-
                                         0.5)+0.5) / (0.5*qLeft);
            } else {
                y = (2. * SRMBlack::NormalH(0.5*sqrt(T)*swaptionVols[i]) - 1.) *
                    sqrt (2. * Maths::PI);
            }
            y = Maths::square(y);
        } else {
            double atmPr = 
                SRMBlack::Option_BS2Q(pY[i], pY[i], T, swaptionVols[i], 'C', 1., 1., 0.);
            y = SRMBlack::ImpVol_BS2Q (pY[i], pY[i], T, atmPr, 'C',
                             qLeft, qRight, fwdShift, swaptionVols[i]);
            y = T * Maths::square(y);
        }

        y -= D * D * L * exp(-2. * beta * t);
        // Matrix for lambda system
        double M = D * D * ratesHJMUtil.expDecay (0, 2. * t) * t;
        /* *  Solve using Gauss-Jordan method */
        if (fabs(M) < SRMConstants::SRM_TINY) {
            throw ModelException(method, "Problem in bootstrapping "+
                                 swaptionExpiries[i].toString()+
                                 " volatility (negative variance)");
        }

        y /= M;
        
        /* * Check if solution is positive */
        if (y < SRMConstants::SRM_TINY) {
            if (!skipFlag) {
                throw ModelException(method, "problem in bootstrapping "+
                                     swaptionExpiries[i].toString()+
                                     "volatility (negative value)");
            }   
            continue;
        }
        
        lambda = sqrt(y);

        /* Relative weights of factors */
        for (int k = LastCalibIdx+1; k <= i; k++) {
            Aweight[k] = lambda * aw;
        }

        /* * Reset LastCalibIdx to current index */
        LastCalibIdx = i;

        /* Update integrals */
        L *= exp (-2. * beta * t);
        L += lambda * lambda * ratesHJMUtil.expDecay(0, 2. * t) * t;
    }  /* for i */


    for (int k = LastCalibIdx+1; k < NbVol; k++) {
        Aweight[k] = lambda * aw;
    }  /* for k */
    
    if (LastCalibIdx == -1) {
        throw ModelException(method, "none of the points calibrated!");
    }
    return Aweight;
}

DRLIB_END_NAMESPACE
