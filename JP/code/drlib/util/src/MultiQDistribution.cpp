//----------------------------------------------------------------------------
//
//   Group       : Credit QRD
//
//   Filename    : MultiQDistribution.cpp
//
//   Description : Helper class for calibrating Multi-Q Distributions
//                 And calculating option prices and greeks
//
//   Author      : Charles Morcom
//
//   Date        : August 5, 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MULTIQDISTRIBUTION_CPP
#include "edginc/MultiQDistribution.hpp"
#include "edginc/Addin.hpp"  
#include "edginc/mathlib.hpp"
#include "edginc/Black.hpp"
#include "edginc/RootFinderND.hpp"
#include "edginc/DiscreteDistribution.hpp"

DRLIB_BEGIN_NAMESPACE

#define MULTIQ_TINY 1E-7
#define MULTIQ_MINQ 1E-8
#define MULTIQ_MAX_STDEV 10.

/* simple function to calculate (exp(qi(x-xi))-1)/qi */
double qxq(double q, double x, double x0) {
  if(fabs(q) < MULTIQ_MINQ) {
    return x-x0;
  } else {
    return (exp(q*(x-x0))-1.0)/q;
  }
}

/*=============================================================================
 * FUNCTION CLASS THAT MAPS A SET OF Qs, muQ, sigQ 
 * TO A FORWARD AND A SET OF OPTION PRICES - USED FOR ROOT FINDER DURING
 * CALIBRATION
 *===========================================================================*/
class MapQsToPriceDiffs : public FunctionNDWithJacobian {
public:

    virtual void operator() (const CDoubleArray&  x,
                            CDoubleArray&        f) const {

        const static string method("MapQsToPriceDiffs::operator()");

        if (x.size()!=(mq.qL.size()+mq.qR.size()+2)) {
            char buf[100];
            sprintf(buf,"Wrong no of vars: nX=%d, nQL=%d, nQR=%d\n", x.size(),mq.qL.size(),mq.qR.size());
            throw ModelException(method, buf);
        }
        // make sure the output vector is cleared.
        int i= 0;
        f.clear();

        /*===================================================================*/
        /* Load up mu, sigma, and qs from argument vector                    */
        /*===================================================================*/
        mq.muQ = x[0];
        mq.sigQ = x[1];
        for (i=0; i<mq.qL.size(); i++) {
            mq.qL[i] = x[i+2];
        }
        for (i=0; i<mq.qR.size(); i++) {
            mq.qR[i] = x[i+2+mq.qL.size()];
        }
        mq.ensureNonZeroQs();

        /*===================================================================*/
        /* Calibrate ds and xs - you already know the ks                     */
        /*===================================================================*/
        mq.xL.clear();
        mq.dL.clear();
        mq.xL.push_back(0.0);
        mq.dL.push_back(1.0);
        for (i=0; i<mq.qL.size()-1; i++) {
            double expQDeltaX = 1 + (mq.kL[i+1]-mq.kL[i])*mq.qL[i]/mq.dL[i];
            mq.dL.push_back(mq.dL[i] * expQDeltaX);
            mq.xL.push_back(mq.xL[i] + log(expQDeltaX)/mq.qL[i]);
        }
        mq.xL.push_back(mq.muQ - mq.sigQ * MULTIQ_MAX_STDEV); // cut off beyond this
        
        /* right side */
        mq.xR.clear();
        mq.dR.clear();
        mq.xR.push_back(0.0);
        mq.dR.push_back(1.0);
        for (i=0; i<mq.qR.size()-1; i++) {
            double expQDeltaX = 1 + (mq.kR[i+1]-mq.kR[i])*mq.qR[i]/mq.dR[i];
            mq.dR.push_back(mq.dR[i] * expQDeltaX);
            mq.xR.push_back(mq.xR[i] + log(expQDeltaX)/mq.qR[i]);
        }
        mq.xR.push_back(mq.muQ + mq.sigQ * MULTIQ_MAX_STDEV); // cut off beyond this

        /*===================================================================*/
        /* Calculate the expected value of the Y variable                    */
        /* The first element of the returned vector is the difference
         * between the calculated forward and the value you are trying to match */
        /*===================================================================*/
        double expectedH = mq.calculateExpectedH();
        double actualForward = mq.f*expectedH;
        f.push_back(actualForward-mq.f);

        /*===================================================================*/
        /* Price all the options - puts to the left, calls to the right, and
         * add the value differences to the result vector                    */
        /*===================================================================*/
        // LEFT SIDE
        for (i=0; i<=mq.qL.size(); i++) {
            double optPrice = mq.vanillaOptionPrice(false /*PUT*/, mq.f*mq.kL[i]);
            f.push_back(optPrice - bsPriceL[i]); 
        }
        // RIGHT SIDE (EXCLUDING THE ATM OPTION - DONE IN LEFT SIDE)
        for (i=1; i<=mq.qR.size(); i++) {
            double optPrice = mq.vanillaOptionPrice(true /*CALL*/, mq.f*mq.kR[i]);
            f.push_back(optPrice - bsPriceR[i]); 
        }

        return;

    };

    MapQsToPriceDiffs(const MultiQDistribution& uncalibratedMQ, const RangeArray& ranges
        ) : FunctionNDWithJacobian(uncalibratedMQ.inputStrikeBoundaryParameters.size()+2, ranges), 
        mq(uncalibratedMQ), bsPriceL(), bsPriceR() {

        // calculate dollar prices for the options
        // ATM price first
        bsPriceL.push_back(mq.atmPx);
        bsPriceR.push_back(mq.atmPx);
        int i;
        int nQL = mq.qL.size();
        for (i=0; i<nQL; i++) {
            double vol = (mq.inputVolParameters)[nQL-i-1];
            bsPriceL.push_back(Black::price(false, mq.f, mq.f*mq.kL[i+1], 1.0, 
                mq.atmT*vol*vol));
        }
        int nQR = mq.qR.size();
        for (i=0; i<nQR; i++) {
            double vol = (mq.inputVolParameters)[nQL+i];
            bsPriceR.push_back(Black::price(true, mq.f, mq.f*mq.kR[i+1], 1.0, 
                mq.atmT*vol*vol));
        }
    };
    
    const MultiQDistribution& mq;
    DoubleArray bsPriceL;
    DoubleArray bsPriceR;

};


double MultiQDistribution::BSQIntegral(double a, double b, double mu, double sigma, double q) {

    if (sigma<=0)
        throw ModelException("MultiQDistribution::BQSIntegral","sigma must be positive!");

    double A = q*sigma;
    double B = q*mu;
    double x1 = (a-mu)/sigma;
    double x2 = (b-mu)/sigma;

    if (fabs(A)>10 && fabs(x1-A)>1 && fabs(x2-A)>1) {
        /* FOR LARGE q sigma, exp(B+0.5*A*A)(N(x2-A)-N(x1-A)) MUST BE CALCULATED 
           IN ONE PIECE, ELSE MULTIPLYING VERY SMALL BY VERY LARGE WILL GIVE
           BAD PRECISION */
        /* N(x2-A)-N(x1-A) = 0.5*(erf(z2)-erf(z1)) this method uses a continued
           fraction representation of erf(x) */
        double z1 = (x1-A)/sqrt(2.0);
        double tz1 = 2.0*z1;
        double z2 = (x2-A)/sqrt(2.0);
        double tz2 = 2.0*z2;
        double cf1 
            = z1+1/(tz1+2/(z1+3/(tz1+4/(z1+5/(tz1+6/(z1+7/(tz1+8/(z1+9/(tz1+10/z1)))))))));
        double cf2 
            = z2+1/(tz2+2/(z2+3/(tz2+4/(z2+5/(tz2+6/(z2+7/(tz2+8/(z2+9/(tz2+10/z2)))))))));

        double w = B + 0.5*A*A;
        double res = 
            (exp(w-z1*z1)/cf1 - exp(w-z2*z2)/cf2)/(2*sqrt(3.14159265358979323846264338328));

        return res;

    } else {
        /* OK TO EVALUATE x1 AND x2 PIECES SEPARATELY */

        double erfc1, erfc2; 
        double erfc = 0; 
        int    sign = 0;

        erfc1 = ExpCErrFcn (A,B,x1); 
        erfc2 = ExpCErrFcn (A,B,x2);

        /* process x1 ********************************************************/
        if (x1-A > 0) {
            erfc += erfc1; 
            sign -= 1;
        } else {
            erfc -= erfc1;
        }
        /* process x2 ********************************************************/
        if (x2-A > 0) {
            erfc -= erfc2;
            sign += 1;
        } else {
            erfc += erfc2;
        }

        /* NormCum has an extra coefficient of 1/2 */
        erfc *= 0.5;
        /* put back value at infinity as in NormCum */
        if (sign) erfc += sign * exp (B + 0.5 * A * A);
 
        return erfc;
    }
}

/*=============================================================================
 * DESTRUCTOR
 *===========================================================================*/
MultiQDistribution::~MultiQDistribution() {};

/*=============================================================================
 * CONVERT INPUT PARAMETERS INTO INTERNAL CALIBRATED REPRESENTATION
 *===========================================================================*/
void MultiQDistribution::calibrateFromInputs() const {

    static const string method("MultiQDistribution::calibrateFromInputs");

    if (calibrated) return; // only do this once

    calibrationFailed = false;
    dontCalibrate = true; // stop recursive recalibration while this method calls others

    /*=========================================================================
     * COMMON CHECKS AND SET-UP
     *=======================================================================*/
    /* VOL AND TIME MUST BE POSITIVE *****************************************/
    if (inputTimeToExpiry<0) {   
        calibrationFailed = true;
        throw ModelException(method, 
                         "Negative time to option expiry");
    } else if (inputATMVolatility<0) {
        calibrationFailed = true;
        throw ModelException(method, 
                         "Negative ATM BS volatility");
    }
    /* FORWARD MUST BE POSITIVE **********************************************/
    if (inputForward<MULTIQ_TINY) {
        calibrationFailed = true;
        throw ModelException(method, 
            "inputForward must be strictly positive.");
    }
    double sigTotal  = inputATMVolatility * sqrt(inputTimeToExpiry);
    /* CHECK TOTAL VOL *******************************************************/
    if (sigTotal<MULTIQ_TINY) {
        calibrationFailed = true;
        throw ModelException(method, 
                         "Total ATM vol (sigma*sqrt(T)) is too small!");
    }
    /* BASIC DISTRIBUTION PARAMETERS *****************************************/
    f = inputForward;
    atmVol = inputATMVolatility;
    atmT = inputTimeToExpiry;
    atmPx = Black::price(true, f, f, 1.0, sigTotal*sigTotal);
    /* DEFAULT CALIBRATION SETTINGS ******************************************/
    C = 1.0;
    K = 1.0;
    muQ = 0.0;
    sigQ = sigTotal;
    /* ENSURE THAT ARRAYS ARE EMPTY ******************************************/
    qL.clear();
    qR.clear();
    dL.clear();
    dR.clear();
    kL.clear();
    kR.clear();
    xL.clear();
    xR.clear();
    /*=======================================================================*/

    if (inputCalibrationType=="explicitQ") {
        calibrateExplicitQAndDeltaInputs();
    } else if (inputCalibrationType=="bsImpliedVol") {
        calibrateBSImpliedVolInputs();
    } else {
        calibrationFailed = true;
        throw ModelException(method, "inputCalibrationType must have value explicitQ or bsImpliedVol!");
    }

    calibrated = true; // if suceeded, then don't do it again
    dontCalibrate = false;
    return;
}

/*=============================================================================
 * CONSTRUCTOR - DIRECTLY FROM ATM OPTION, QS AND DELTAS
 *===========================================================================*/
MultiQDistribution::MultiQDistribution(
        double              optionTimeToExpiry,
        double              forwardRate,
        double              atmBSVolatility,
        const DoubleArray&  inputQs,
        const DoubleArray&  inputDeltas
        ) : CObject(TYPE), inputTimeToExpiry(optionTimeToExpiry), inputForward(forwardRate),
            inputATMVolatility(atmBSVolatility), inputCalibrationType("explicitQ"), 
            calibrated(false), calibrationFailed(false), dontCalibrate(false), inputVolParameters(inputQs),
            inputStrikeBoundaryParameters(inputDeltas) {

    calibrateFromInputs();
}

MultiQDistribution::MultiQDistribution(
        double              optionTimeToExpiry,
        double              forwardRate,
        double              atmBSVolatility,
        const DoubleArray&  inputQs,
        const DoubleArray&  inputDeltas, /**<Must be ordered.*/
        const string&       inputType
        )
    : CObject(TYPE), inputTimeToExpiry(optionTimeToExpiry), inputForward(forwardRate),
      inputATMVolatility(atmBSVolatility), inputCalibrationType(inputType), 
      calibrated(false), calibrationFailed(false), dontCalibrate(false), inputVolParameters(inputQs),
      inputStrikeBoundaryParameters(inputDeltas) {

    calibrateFromInputs();
}

 //Create a log normal multiQ distribution
MultiQDistributionSP MultiQDistribution::createLogNormalMq(double              optionTimeToExpiry,
							   double              forward,
							   double              atmBSVolatility)

{
  DoubleArray logNormalQs(6,0.0);
  DoubleArray logNormalDs(5);
  logNormalDs[0]=0.2;
  logNormalDs[1] = 0.4;
  logNormalDs[2]=0.5;
  logNormalDs[3]=0.6;
  logNormalDs[4]=0.8;
  return MultiQDistributionSP(new MultiQDistribution(optionTimeToExpiry,
						     forward,
						     atmBSVolatility,
						     logNormalQs,
						     logNormalDs));

}
/*=============================================================================
 * CALLED BY CONSTRUCTOR TO CALIBRATE FROM EXPLICIT Qs AND STRIKE DELTAS
 *===========================================================================*/
void MultiQDistribution::calibrateExplicitQAndDeltaInputs() const {
    
    const static string method("MultiQDistribution::calibrateExplicitQAndDeltaInputs");
    int i; 

    /* CHECK THAT THERE ARE THE RIGHT NUMBER OF QS AND DELTAS */
    if (inputStrikeBoundaryParameters.size()!=(inputVolParameters.size()-1)) {
        calibrationFailed = true;
        throw ModelException(method, 
            "If there are N q values (inputVolParameters), there must be exactly (N-1) strike deltas (inputStrikeBoundaryParameters).");
    }

    /* CHECK ORDERING AND VALUES OF DELTAS, AND FIND DELTA=0.5 MIDPOINT */
    int halfDeltaIndex = -1;
    double nextDelta = -1;
    for (i=0; i < inputStrikeBoundaryParameters.size()-1; i++) {
        double delta = inputStrikeBoundaryParameters[i];
        double nextDelta = inputStrikeBoundaryParameters[i+1];
        if (fabs(delta-0.5)<MULTIQ_TINY) halfDeltaIndex=i;
        double deltaDiff = 
            nextDelta-delta;
        if (deltaDiff<0) {
            calibrationFailed = true;
            throw ModelException(method, 
                         "strike deltas in inputStrikeBoundaryParameters are not in ascending order.");
        } else if (deltaDiff < MULTIQ_TINY) {
            calibrationFailed = true;
           throw ModelException(method, 
                         "strike deltas in inputStrikeBoundaryParameters are not distinct.");
        }
    }
    int nQL;
    int nQR;
    if (halfDeltaIndex>=0) {
        nQL = halfDeltaIndex+1;
        nQR = inputStrikeBoundaryParameters.size() - halfDeltaIndex;
    } else if (fabs(nextDelta-0.5)<MULTIQ_TINY) {
        nQL = inputStrikeBoundaryParameters.size();
        nQR = 1;
    } else if (inputStrikeBoundaryParameters.size()==1 
        && fabs(inputStrikeBoundaryParameters[0]-0.5)<MULTIQ_TINY) {
        nQL = 1;
        nQR = 1;
    } else {
        calibrationFailed = true;
        throw ModelException(method, 
                         "One of the strike deltas must be 0.5 (corresponds to ATM forward strike).");
    }
    /* NOTE THAT ENSURING DELTA=0.5 IS ONE OF THE BREAKS MEANS THAT xL[0]=xR[0] = 0 */
    
    /* LEFT SIDE */
    for(i = 0; i < nQL; i++) {
        qL.push_back(inputVolParameters[nQL - i - 1]);
        xL.push_back(sigQ * N1InverseBetter(inputStrikeBoundaryParameters[nQL - i - 1]));
    }
    xL.push_back(-sigQ * MULTIQ_MAX_STDEV); // cut off beyond this
    /* Build k and d coefficients using freedom to set first values 
     * so that H(0)=1 and H'(0)=1 and then
     * by requiring continuity and differentiability */
    kL.push_back(1.0);
    dL.push_back(1.0);
    for(i = 1; i<nQL; i++) {
        dL.push_back(dL[i-1]*exp(qL[i-1] * (xL[i]-xL[i-1])));
        kL.push_back(kL[i-1] + dL[i-1] * qxq(qL[i-1], xL[i], xL[i-1]));
    }  

    /* RIGHT SIDE */
    for(i = 0; i < nQR; i++) {
        qR.push_back(inputVolParameters[nQR + i]);
        xR.push_back(sigQ * N1InverseBetter(inputStrikeBoundaryParameters[nQR + i - 1]));
    }
    xR.push_back(sigQ * MULTIQ_MAX_STDEV); // cut off beyond this
    /* Build k and d coefficients using freedom to set first values 
     * so that H(0)=1 and H'(0)=1 and then
     * by requiring continuity and differentiability */
    kR.push_back(1.0);
    dR.push_back(1.0);
    for(i = 1; i<nQR; i++) {
        dR.push_back(dR[i-1]*exp(qR[i-1] * (xR[i]-xR[i-1])));
        kR.push_back(kR[i-1] + dR[i-1] * qxq(qR[i-1], xR[i], xR[i-1]));
    }

    /* Make sure qs are bounded away from zero */
    ensureNonZeroQs();

    /* NOW COMPUTE C AND K SO THAT f IS THE EXPECTATION, AND SO THAT THE
       ATM OPTION PRICE MATCHES BLACK-SCHOLES */
    calibrateAffineAdjustments();
 
};

void MultiQDistribution::calibrateBSImpliedVolInputs() const {

    const static string method("MultiQDistribution::calibrateBSImpliedVolInputs");
    int i; 

    /* CHECK THAT THERE ARE THE RIGHT NUMBER OF VOLS AND STRIKES */
    if (inputStrikeBoundaryParameters.size()!=inputVolParameters.size()) {
        calibrationFailed = true;
        throw ModelException(method, 
            "If there are N implied vols (inputVolParameters), there must be exactly N strike ratios (inputStrikeBoundaryParameters).");
    }
    
    /* CHECK ORDERING AND SIZE OF INPUTS */
    int switchIdx = -1;
    double lastVol = inputVolParameters[0];
    double lastStrike = inputStrikeBoundaryParameters[0];
    if (lastStrike>1.0-MULTIQ_TINY) {
        calibrationFailed = true;
        throw ModelException(method,
            "The lowest strike ratio must be below ATM (1.0)");
    }
    for (i=1; i<inputVolParameters.size(); i++) {
        double vol = inputVolParameters[i];
        double strike = inputStrikeBoundaryParameters[i];

        if (vol<MULTIQ_TINY) {
            calibrationFailed = true;
            throw ModelException(method,"BS implied vols must all be strictly positive.");
        } else if (strike<lastStrike) {
            calibrationFailed = true;
            throw ModelException(method,"Strike ratios must be strictly increasing.");
        } else if (strike<lastStrike+MULTIQ_TINY) {
            calibrationFailed = true;
            throw ModelException(method,"Strike ratios must be distinct.");
        } else if (fabs(strike-1.0)<MULTIQ_TINY) {
            calibrationFailed = true;
            throw ModelException(method,
                "Input impled vols (inputVolParameters) and strike ratios "
                "(inputStrikeBoundaryParameters) should not include the ATM strike (strike ratio 1.0).");
        }
        if (switchIdx<0 && strike>1.0) switchIdx = i;

        lastVol = vol;
        lastStrike = strike;
    }
    if (lastStrike<1+MULTIQ_TINY) {
        calibrationFailed = true;
        throw ModelException(method,
            "The highest strike ratio must be above ATM (1.0)");
    }

    int nQL = switchIdx;
    int nQR = inputVolParameters.size()-nQL;

    /* LEFT SIDE */
    kL.push_back(1.0);
    for(i = 0; i < nQL; i++) {
        qL.push_back(1.0); // need dummy values to start!
        kL.push_back(inputStrikeBoundaryParameters[nQL - i - 1]);
    }

    /* RIGHT SIDE */
    kR.push_back(1.0);
    for(i = 0; i < nQR; i++) {
        qR.push_back(1.0); // dummy values to start
        kR.push_back(inputStrikeBoundaryParameters[nQL + i]);
    }
    
    /* NOW CREATE THE MAPPING FROM QS TO OPTION PRICES */
    Infinity plusInfinity(Infinity::Minus);
    ClosedBoundary zeroBoundary(0);
    RangeArray ranges(2+nQL+nQR);
    ranges.push_back(InfiniteRangeSP()); // muQ
    ranges.push_back(RangeSP(new Range(zeroBoundary, plusInfinity))); //sigQ
    for (i=0; i<nQL+nQR; i++) {
        ranges.push_back(InfiniteRangeSP());
    }
    MapQsToPriceDiffs mqf(*this, ranges);
    RootFinderNDNewtonSafe rootFinder(
        50 /* max iterations */, 0.00001 /* tolx */, 
        0.00001 /* tolf */, 100);

    /*=========================================================================
     * SOLVE FOR THE Qs
     *=======================================================================*/
    DoubleArray initialGuess;
    initialGuess.push_back(muQ); // muQ = 0
    initialGuess.push_back(sigQ); // muQ = totalVol
    for (i=0; i<inputVolParameters.size(); i++) {
        initialGuess.push_back(1.0); // initial guesses are lognormal for all qs.
    }
    rootFinder.solve(mqf, initialGuess);

    if (29<rootFinder.getInfo()->iter) {
        calibrationFailed = true;
        throw ModelException(method, "MultiQDistribution failed to converge after 30 iterations.");
    }
    /* NOTE THAT THE Qs ARE ALSO CHANGED DIRECTLY IN mqf, SO NO NEED TO
       COPY SOLUTIONS TO this */
    
    return;
};

/*=============================================================================
 * qMap()
 *===========================================================================*/
double MultiQDistribution::qMap(double xVal) const {

    static const char* method = "MultiQDistribution::qMap";

    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    /* for zero volatility, return forward no matter what the point is */
    if (sigQ < MULTIQ_TINY) {
        return f;
    }   
   
    int   nq; /* number of intervals      */
    bool isLeft;
    int   idx;
    double hVal;
    locateInGaussianSpace(xVal, true, &isLeft, &nq, &idx, &hVal);

    /* add effect of affine adjustments */
	return f*(K*(C * hVal - 1.0) + 1.0);
	
}
double MultiQDistribution::getNormalVol() const {

    static const char* method = "MultiQDistribution::getNormalVol";

    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();
return sigQ;
}
double MultiQDistribution::calculateExpectedH() const {

    const static string method("MultiQDistribution::calculateExpectedH");
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    double expectedH = 0.0;
    
    if (sigQ<MULTIQ_TINY) {
        expectedH = qMap(muQ);
    } else {
        /* FIND THE EXPECTATION OF THE MAPPING FUNCTION, H(x) */
        /* Left side first */
        int i;
        for (i=0; i<qL.size(); i++) {
            double qTemp = (fabs(qL[i])<MULTIQ_MINQ ? MULTIQ_MINQ : qL[i]);
            double d2q = dL[i]/qTemp;
            expectedH += (kL[i] - d2q) * fabs(BSQIntegral(xL[i], xL[i+1], muQ, sigQ , 0.0));
            expectedH += d2q*exp(-qTemp*xL[i]) 
                * fabs(BSQIntegral(xL[i], xL[i+1], muQ, sigQ , qTemp));
        }
        /* Now right side */
        for (i=0; i<qR.size(); i++) {
            double qTemp = (fabs(qR[i])<MULTIQ_MINQ ? MULTIQ_MINQ : qR[i]);
            double d2q = dR[i]/qTemp;
            expectedH += (kR[i] - d2q) * fabs(BSQIntegral(xR[i], xR[i+1], muQ, sigQ , 0.0));
            expectedH += d2q*exp(-qTemp*xR[i]) 
                * fabs(BSQIntegral(xR[i], xR[i+1], muQ, sigQ , qTemp));
        }
    }
    return expectedH;
}

/*=============================================================================
 * calibrateAffineAdjustments()
 *===========================================================================*/
void MultiQDistribution::calibrateAffineAdjustments() const {
    
    static const char* method = "MultiQDistribution::calibrateAffineAdjustments";
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    K = 1.0;
    C = 1.0;

    double expectedH = calculateExpectedH();
    if (fabs(expectedH)<MULTIQ_TINY) {
        throw ModelException(method, 
          "E[H(x)] is effectively zero - unable to calibrate affine parameters.");
    }
    C = 1/expectedH;

    /* Now, since E[Y] = fwdRate, (K (C expecH-1) + 1) = 1. Also, if fwdX is the point
       in gaussian X-space corresonding to the forward rate,
       fwdRate = fwdRate (K (C H(fwdX)-1) + 1)
       fwdX = H_inv(expecH)
       Since the price of the call is E[(Y-fwdRate)+] = fwdRate K C E[(H(x) - H(fwdX)+]
       we can price the call assuming that volC=1 and then rescale to get the correct call price
       */
    //calibrated = true; // avoid infinite loop!
    double k1CallPrice = this->vanillaOptionPrice(true, f);
    if (k1CallPrice<MULTIQ_TINY) {
        throw ModelException(method, 
          "ATM Multi-Q Call Price with K=1 is too small!");
    }
    //calibrated = false; // set back to false
 
    /* Now set K so that the ATM option price matches, and the expectation is the
       forward rate */
    K = atmPx/k1CallPrice;
    return;
}

/**Ensure all qs are bounded away from zero */
void MultiQDistribution::ensureNonZeroQs() const {
    int i;

    /* LEFT SIDE */
    for (i=0; i<qL.size(); i++) {
        if (fabs(qL[i])<MULTIQ_MINQ) {
            qL[i] = (qL[i]<0 ? -MULTIQ_MINQ : MULTIQ_MINQ);
        }
    }

    /* RIGHT SIDE */
    for (i=0; i<qR.size(); i++) {
        if (fabs(qR[i])<MULTIQ_MINQ) {
            qR[i] = (qR[i]<0 ? -MULTIQ_MINQ : MULTIQ_MINQ);
        }
    }
}

/*=============================================================================
 * forwardPrice()
 *===========================================================================*/
double MultiQDistribution::forward() const {
    static const string method("MultiQDistribution::forward");
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    return this->f;
}

/**Find a point's location in gaussian x-space.*/
void MultiQDistribution::locateInGaussianSpace(
    double val, /**<Value to find*/
    bool isX, /**<Is val a forward space value, or an x-space value?*/
    bool* isLeft, /**<Returns true if left-side, else false*/
    int* nQ, /**<Number of qs on side*/
    int* idx, /**<(*x)[idx] is just closer to 0 than val */
    double* hVal /**<value of H(x) corresponding to x or Y=F(K(CH(x)-1)+1)*/
    ) const {

    static const char* method = "CMutiQDistribution::locateInGaussianSpace";
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    /* If val is in y-space, then locate based on H(x), else use x directly */
    double strkMultMQ;
    if (!isX) {
        double strkMult = val/f;
        strkMultMQ = (K-1)/(K*C) + strkMult/ (K*C);
    }

    /* pick left or right side, depending on where point is */
    int sign;
    DoubleArray* q;
    DoubleArray* k;
    DoubleArray* d;
    DoubleArray* x;
    if (isX && val<0 || !isX && strkMultMQ < 1.) {
        *isLeft = true;
        *nQ = qL.size();
        q = &qL;
        k = &kL;
        d = &dL;
        x = &xL;
        sign = -1;
    
    } else {
        *isLeft=false;
        *nQ=qR.size();
        q = &qR;
        k = &kR;
        d = &dR;
        x = &xR;
        sign = 1;
    }
    DoubleArray* searchVec = (isX ? x : k);
      
    /* find strike location */
    *idx = 0;
    while ((*idx)<(*nQ)-1 && sign*((isX ? val : strkMultMQ)-(*searchVec)[(*idx)+1])>= 0) 
        (*idx)++;

    *hVal = (isX ? (*k)[*idx] + (*d)[*idx] * qxq((*q)[*idx], val, (*x)[*idx]) : strkMultMQ);

    return;
    
}

/*=============================================================================
 * vanillaOptionPrice()
 *===========================================================================*/
double MultiQDistribution::vanillaOptionPrice (
    bool   isCall, /**<true if call, else put*/
    double strike  /**<Option strike price */
    ) const {

    static const char* method = "CMutiQDistribution::vanillaOptionPrice";
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();
    double copSign = (isCall ? 1.0 : -1.0);

    /* for zero volatility, price intrinsic value */
    if (sigQ < MULTIQ_TINY) {
        return max(copSign * (f - strike), 0.0);
    }

    /* pick left or right integration, depending on where strike is */
    DoubleArray* q;
    DoubleArray* k;
    DoubleArray* d;
    DoubleArray* x;
    int nQ;
    int idx;
    double hVal;
    bool isLeft;
    locateInGaussianSpace(strike,false,&isLeft,&nQ,&idx,&hVal);
    if (isLeft) {
        q = &qL;
        k = &kL;
        d = &dL;
        x = &xL;
    } else {
        q = &qR;
        k = &kR;
        d = &dR;
        x = &xR;
    }
    
    /* find normal-space point corresponding to strike */
    double d2q = (*d)[idx]/(*q)[idx];
    double temp = 1 + (hVal - (*k)[idx])/d2q;
    double xtemp = (temp > 0 ? (*x)[idx] + log(temp)/(*q)[idx] : (*x)[nQ]);

    /* key step(!): price put for low strikes, call for high strikes */

    /* last interval (closest to ATM) */
    temp = exp(-(*q)[idx] * (*x)[idx]);
    double price = 
        hVal      * BSQIntegral((*x)[nQ], xtemp, muQ, sigQ, 0) -
        ((*k)[idx]-d2q) * BSQIntegral((*x)[idx+1], xtemp, muQ, sigQ, 0) -
        d2q*temp        * BSQIntegral((*x)[idx+1], xtemp, muQ, sigQ, (*q)[idx]);

    /* other intervals further away from x=0 */
    for (int i=idx+1; i<nQ; i++) 
    {
        d2q = (*d)[i]/(*q)[i];
        temp = exp(-(*q)[i] * (*x)[i]);
        price -= 
            ((*k)[i] - d2q) * BSQIntegral((*x)[i+1], (*x)[i], muQ, sigQ, 0) +
             d2q * temp     * BSQIntegral((*x)[i+1], (*x)[i], muQ, sigQ, (*q)[i]);
    }

    /* apply affine adjustment and rescale by forward */
    price *= K * C * f;
    /* Call/Put parity to convert option type */
    if (isCall && strike<f*C*K) {
        price += f - strike;
    } else if (!isCall && strike>f*C*K) {
        price += strike - f;
    }
    return price;

}

//--------------------------------------------------------------------
//   IDistribution1D methods
//--------------------------------------------------------------------
/**Cumulative density function of the mapped value. P(Y<=y) = mq.cdf(y)
  */
double MultiQDistribution::cdf(double y) const {

    static const char* method = "CMutiQDistribution::cdf";
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    /* for zero volatility, cdf is 1.0 if above forward, else zero */
    if (sigQ < MULTIQ_TINY) {
        return (y>=f ? 1.0 : 0.0);
    }

    /* pick left or right integration, depending on where strike is */
    DoubleArray* q;
    DoubleArray* k;
    DoubleArray* d;
    DoubleArray* x;
    int nQ;
    int idx;
    double hVal;
    bool isLeft;
    locateInGaussianSpace(y,false,&isLeft,&nQ,&idx,&hVal);
    if (isLeft) {
        q = &qL;
        k = &kL;
        d = &dL;
        x = &xL;
    } else {
        q = &qR;
        k = &kR;
        d = &dR;
        x = &xR;
    }

    /* find normal point corresponding to strike */
    double d2q = (*d)[idx]/(*q)[idx];
    double temp = 1 + (hVal - (*k)[idx])/d2q;
    double xtemp = (temp > 0 ? (*x)[idx] + log(temp)/(*q)[idx] : (*x)[nQ]);

    return N1((xtemp - muQ) / sigQ);
    
}

/**Probability density function of the mapped value: P(y<=Y<=y+dy) = mq.pdf(y) dy.*/
double MultiQDistribution::pdf(double y) const {

    static const char* method = "CMutiQDistribution::pdf";
    if (calibrationFailed) throw ModelException(method, "MultiQDistribution object cannot be calibrated.");
    if (!calibrated && !dontCalibrate) calibrateFromInputs();

    /* for zero volatility, pdf not defined*/
    if (sigQ < MULTIQ_TINY) {
        throw ModelException(method, 
            "The pdf of a multi-q distribution is not well defined for a very small total volatility.");
    }

    /* pick left or right integration, depending on where strike is */
    DoubleArray* q;
    DoubleArray* k;
    DoubleArray* d;
    DoubleArray* x;
    int nQ;
    int idx;
    double hVal;
    bool isLeft;
    locateInGaussianSpace(y,false,&isLeft,&nQ,&idx,&hVal);
    if (isLeft) {
        q = &qL;
        k = &kL;
        d = &dL;
        x = &xL;
    } else {
        q = &qR;
        k = &kR;
        d = &dR;
        x = &xR;
    }

    /* find normal point corresponding to strike */
    double d2q = (*d)[idx]/(*q)[idx];
    double temp = 1 + (hVal - (*k)[idx])/d2q;
    double xtemp = (temp > 0 ? (*x)[idx] + log(temp)/(*q)[idx] : (*x)[nQ]);

   return N1Density((xtemp - muQ) / sigQ)
     / sigQ
     / ( (*d)[idx] * temp)
     / (f*K*C);

}

double MultiQDistribution::expectation() const {
    return forward();
}

double MultiQDistribution::expectation(const Function1DDouble& payoff) const {
    throw ModelException("MultiQDistribution::expactation(payoff)", "not implemented yet");
}

double MultiQDistribution::mgf(double z, int n) const {
    throw ModelException("MultiQDistribution::mgf", "not implemented yet");
    return 0;
}

double MultiQDistribution::variance() const {
    throw ModelException("MultiQDistribution::variance", "not implemented yet");
    return 0;
}

DiscreteDistributionConstSP MultiQDistribution::discretise() const {
    throw ModelException("MultiQDistribution::discretise", "not implemented yet");
}    


/*=============================================================================
 * binaryOptionPrice()
 *===========================================================================*/
double MultiQDistribution::binaryOptionPrice (
    bool   isCall, /**<true if call, else put*/
    double strike  /**<Option strike price */
    ) const {

    double probabilityBelow = cdf(strike);
    return (isCall ? 1-probabilityBelow : probabilityBelow);
}

/*=============================================================================
 * Reflection static load method
 *===========================================================================*/
void MultiQDistribution::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(MultiQDistribution, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultMultiQDistribution);

    /* INTERFACE FIELDS */
    FIELD(inputForward, "Forward rate/price at option expiry.");
    FIELD(inputTimeToExpiry, "Time to option expiry in years.");
    FIELD(inputATMVolatility, "Black-Scholes implied vol of ATM option.");
    FIELD(inputCalibrationType, "How should the input parameters be interpreted? explicitQ=>(q, delta); bsImpliedVol=>(vols, strike ratios).");
    FIELD(inputVolParameters, "Volatility calibration parameters.");
    FIELD(inputStrikeBoundaryParameters, "Strike-boundary calibration parameters.");

};

/*=============================================================================
 * Register class loader
 *===========================================================================*/
CClassConstSP const MultiQDistribution::TYPE = 
    CClass::registerClassLoadMethod("MultiQDistribution", typeid(MultiQDistribution), load);

/*=============================================================================
 * Default Constructor (needed for reflection)
 *===========================================================================*/
MultiQDistribution::MultiQDistribution() : CObject(TYPE), calibrated(false), 
    calibrationFailed(false), dontCalibrate(false) {};

/*=============================================================================
 * Reflection callback for default constructor.
 *===========================================================================*/
IObject* MultiQDistribution::defaultMultiQDistribution() {return new MultiQDistribution();};


/*=============================================================================
 * Class to provide add-in functions
 *===========================================================================*/
class MultiQDistributionAddin : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    MultiQDistributionSP mqDist;
    string fieldName;
    int arrayIndex;
 
    double getValue() {
        if (mqDist->calibrationFailed) throw ModelException("getValue", "MultiQDistribution object cannot be calibrated.");
        if (!mqDist->calibrated) mqDist->calibrateFromInputs();
        if (fieldName=="forward") {
            return mqDist->f;
        } else if (fieldName=="muQ") {
            return mqDist->muQ;
        } else if (fieldName=="sigQ") {
            return mqDist->sigQ;
        } else if (fieldName=="atmVol") {
            return mqDist->atmVol;
        } else if (fieldName=="atmPx") {
            return mqDist->atmPx;
        } else if (fieldName=="atmT") {
            return mqDist->atmT;
        } else if (fieldName=="C") {
            return mqDist->C;
        } else if (fieldName=="K") {
            return mqDist->K;
        } else if (fieldName=="nQL") {
            return (mqDist->qL).size();
        } else if (fieldName=="nQR") {
            return (mqDist->qR).size();
        } else {
            throw ModelException("MultiQDistributionAddin::getValue",
                fieldName+" is not a legitimate field name.");
        }
    }

    double getArrayValue() {
        DoubleArray* ap;
        if (mqDist->calibrationFailed) throw ModelException("getArrayValue", "MultiQDistribution object cannot be calibrated.");
        if (!mqDist->calibrated) mqDist->calibrateFromInputs();
        if (fieldName=="qL") {
            ap= &(mqDist->qL);
        } else if (fieldName=="qR") {
            ap= &(mqDist->qR);
        } else if (fieldName=="dL") {
            ap= &(mqDist->dL);
        } else if (fieldName=="dR") {
            ap= &(mqDist->dR);
        } else if (fieldName=="kL") {
            ap= &(mqDist->kL);
        } else if (fieldName=="kR") {
            ap= &(mqDist->kR);
        } else if (fieldName=="xL") {
            ap= &(mqDist->xL);
        } else if (fieldName=="xR") {
            ap= &(mqDist->xR);
        } else {
            throw ModelException("MultiQDistributionAddin::getArray",
                fieldName+" is not a legitimate field name.");
        }
        if (arrayIndex<0 || arrayIndex>=ap->size()) {
            return -99999.0;
        } else {
            return (*ap)[arrayIndex];
        }
    }



    
    MultiQDistributionAddin() : CObject(TYPE), fieldName(""), arrayIndex(0) {};

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(MultiQDistributionAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMultiQDistributionAddin);
        FIELD(mqDist, "MultiQDistribution object.");
        FIELD(fieldName, "Name of field whose value you want.");
        FIELD(arrayIndex, "Index into array.");
        FIELD_MAKE_OPTIONAL(arrayIndex);

        Addin::registerDoubleMethod("MULTIQ_GET_VALUE", 
            Addin::UTILITIES, 
            "Returns the value of a field of a multi-q distribution.", 
            &MultiQDistributionAddin::getValue);
        
        Addin::registerDoubleMethod("MULTIQ_GET_ARRAY_VALUE",
            Addin::UTILITIES, 
            "Returns a value array",
            &MultiQDistributionAddin::getArrayValue);
    }

    static IObject* defaultMultiQDistributionAddin() {
        return new MultiQDistributionAddin();
    }
};

CClassConstSP const MultiQDistributionAddin::TYPE = CClass::registerClassLoadMethod(
    "MultiQDistributionAddin", typeid(MultiQDistributionAddin), MultiQDistributionAddin::load);


/*=============================================================================
 * Class to provide add-in functions - values and option prices
 *===========================================================================*/
class MultiQDistributionAddin2 : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    MultiQDistributionSP mqDist;
    double value;
    string callOrPut;

    double qMap() {
        return mqDist->qMap(value);
    }

    double pdf() {
        return mqDist->pdf(value);
    }

    double cdf() {
        return mqDist->cdf(value);
    }

    double vanillaPrice() {
        bool isCall = false;
        if (callOrPut.size()<1) 
            throw ModelException("MultiQDistributionAddin::vanillaPrice",
            "callOrPut may not be blank and must begin with 'C' for call or 'P' for put.");
        char cp = toupper(callOrPut[0]);
        if (cp!='C' && cp!='P')
            throw ModelException("MultiQDistributionAddin::vanillaPrice",
            "callOrPut must begin with 'C' for call or 'P' for put.");

        if (cp=='C') isCall = true;

        return mqDist->vanillaOptionPrice(isCall, value);
    }

    double binaryPrice() {
        bool isCall = false;
        if (callOrPut.size()<1) 
            throw ModelException("MultiQDistributionAddin::binaryPrice",
            "callOrPut may not be blank and must begin with 'C' for call or 'P' for put.");
        char cp = toupper(callOrPut[0]);
        if (cp!='C' && cp!='P')
            throw ModelException("MultiQDistributionAddin::binaryPrice",
            "callOrPut must begin with 'C' for call or 'P' for put.");

        if (cp=='C') isCall = true;

        return mqDist->binaryOptionPrice(isCall, value);
    }
    
    MultiQDistributionAddin2() : CObject(TYPE), value(0), callOrPut("") {};

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(MultiQDistributionAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMultiQDistributionAddin2);
        FIELD(mqDist, "MultiQDistribution object.");
        FIELD(value,"Input value or strike for addin functions");
        FIELD(callOrPut,"Should begin with 'C' for call or 'P' for put.");
        FIELD_MAKE_OPTIONAL(callOrPut);

        Addin::registerDoubleMethod("MULTIQ_QMAP",
            Addin::UTILITIES,
            "Maps a gaussian driver X value to F(K(CH(X)-1)+1), where H(X) is the q-mapping.",
            &MultiQDistributionAddin2::qMap);

        Addin::registerDoubleMethod("MULTIQ_VANILLA_PRICE",
            Addin::UTILITIES,
            "Prices a vanilla call or put with given strike.",
            &MultiQDistributionAddin2::vanillaPrice);

        Addin::registerDoubleMethod("MULTIQ_BINARY_PRICE",
            Addin::UTILITIES,
            "Prices a binary call or put with given strike.",
            &MultiQDistributionAddin2::binaryPrice);

        Addin::registerDoubleMethod("MULTIQ_CDF",
            Addin::UTILITIES,
            "Returns the cumulative probability density of a given value.",
            &MultiQDistributionAddin2::cdf);

        Addin::registerDoubleMethod("MULTIQ_PDF",
            Addin::UTILITIES,
            "Returns the probability density function at a given value.",
            &MultiQDistributionAddin2::pdf);
    }

    static IObject* defaultMultiQDistributionAddin2() {
        return new MultiQDistributionAddin2();
    }
};

CClassConstSP const MultiQDistributionAddin2::TYPE = CClass::registerClassLoadMethod(
    "MultiQDistributionAddin2", typeid(MultiQDistributionAddin2), MultiQDistributionAddin2::load);

bool MultiQDistributionLinkIn() {
    return MultiQDistribution::TYPE!=0 && MultiQDistributionAddin::TYPE!=0 && MultiQDistributionAddin2::TYPE!=0;
}



DRLIB_END_NAMESPACE
