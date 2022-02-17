//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Black.cpp
//
//   Description : Black-Scholes forumla
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Black.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

/** price2Q
    Price of a call or put using 2Q version of Black&Scholes.
    Note: vol is defined by Y = Y/(1+fwdSh) * F((1+fwdSh)/(1+q*fwdSh)*vol)
    From SRMSwaption.cpp:Option_BS2Q */

static const double ccDelta = 0.0001;   /* needed by double Black::convexity2Q */
static const double ccEqErr = 0.000001; /* needed by double Black::convexity2Q */
static const double qCutoff = 1E-4;     /* Normal model for |q|<qCutoff, needed by double Black::*2Q functions */
static const double tinyDouble = 1E-10;

double Black::price2Q(
    bool   isCall,
    double fwd,      /* Forward yield */
    double strike,
    double pv,
    double variance,
    double qLeft,
    double qRight,
    double fwdSh     /* Forward shift */
) {
    /* lognormal case */
    if (Maths::areEqualWithinTol(qLeft, 1.0, tinyDouble) &&
        Maths::areEqualWithinTol(qRight, 1.0, tinyDouble) &&
        Maths::areEqualWithinTol(fwdSh, 0.0, tinyDouble)
    ) {
        return price(isCall, fwd, strike, pv, variance);
    }
    if (fwd == 0.) return (0.);

    double
        C, P,     /* Call, put price               */
        d,        /* d in N(d) in Black & Scholes  */
        fwd1, strike1;   /* Modified fwd and strike       */

    /* Q for fwd shift correction */
    double QM = (qLeft + qRight) / 2; 

    /* sqrt(variance) * Q */
    double St = sqrt(variance) * (1 + fwdSh) / (1 + QM * fwdSh); 

    /* Yield correction */
    double YCorr = fwd / (1. + fwdSh);

    /* Calibration const */
    double CC = convexity2Q (St, qLeft, qRight, fwdSh);

    /* Q used in part of density */
    double Q = (strike > YCorr) ? qRight : qLeft;
    
    if (fabs(Q) > qCutoff) {
        d   = (CC - log ((strike / YCorr - 1.) * Q + 1.) / Q) / St;

        fwd1  = exp (CC * Q + 0.5 * St * St * Q * Q);
        fwd1 *= YCorr / Q;
        strike1  = strike - YCorr * (1. - 1. / Q);

        if (strike > YCorr) {
            C = + fwd1 * N1 (+ d + St * Q) - strike1 * N1 (+ d);
            P = C - fwd + strike;
        } else {
            P = - fwd1 * N1 (- d - St * Q) + strike1 * N1 (- d);
            C = fwd - strike + P;
        }
    } else {
        d   = (CC - (strike / YCorr - 1.)) / St;

        strike1  = strike - YCorr * (1. + CC);

        if (strike > YCorr) {
            C = + YCorr * St * N1Density(d) - strike1 * N1 (+ d);  
            P = C - fwd + strike;
        } else {
            P = + YCorr * St * N1Density(d) + strike1 * N1 (- d);  
            C = fwd - strike + P;
        }
    }
    return pv * (isCall ? C : P);
}


/** convexity2Q
    Calibrate constant in 2q formulas. From SRMSwaption.cpp:ConvexityC_BS2Q */

double Black::convexity2Q(
    double vol,    /* Annualized volatility */
    double qLeft,  /* Q left                */
    double qRight, /* Q right               */
    double fwdSh   /* Fwd shift             */
) {
    int count = 0;

    double 
        QM,     
        CC0,
        CCEq,
        CCEq2;

    QM = (qLeft + qRight) * 0.5;
    if (fabs (QM) > qCutoff) {
        CC0 = log (fwdSh * QM + 1.) / QM - 0.5 * QM * vol * vol;
    } else {
        CC0 = fwdSh;
    }

    bool found = false;
    do {
        CCEq  = equation2Q(vol, qLeft, qRight, fwdSh, CC0);
        CCEq2 = equation2Q(vol, qLeft, qRight, fwdSh, CC0 + ccDelta);
        
        if (fabs(CCEq) > ccEqErr) {
            CC0  -= CCEq * ccDelta / (CCEq2 - CCEq);
        } else {
            found = true;
        }
        count++ ;
    }
    while (!found && (count < 10));

    if (count < 10) {
        return CC0;
    } else {
        throw ModelException("ConvexityC_BS2Q", "Failed to calibrate constant in 2q formulas");
    }
}

/** equation2Q
    Evalute calibration equation for calibration constant.
    From SRMSwaption.cpp:CCEquation_BS2Q */

double Black::equation2Q(
    double vol,    /* Annualized volatility */
    double qLeft,  /* Q left                */
    double qRight, /* Q right               */
    double fwdSh,  /* Fwd shift             */
    double CC      /* Constant              */
) {
    double SQR  = vol * qRight;
    double SQL  = vol * qLeft;
    double CCIS = CC / vol;
    double CCEq = fwdSh;

    if (fabs (qRight) > qCutoff) {
        CCEq -= (exp (CCIS * SQR + 0.5 * SQR * SQR) * N1 (+ CCIS + SQR)
                 - N1 (+ CCIS) ) / qRight;
    } else {
        CCEq -= (CCIS * N1(+ CCIS) + N1Density(CCIS)) * vol;  
    }

    if (fabs (qLeft) > qCutoff) {
        CCEq -= (exp (CCIS * SQL + 0.5 * SQL * SQL) * N1 (- CCIS - SQL)
                 - N1 (- CCIS) ) / qLeft;
    } else {
        CCEq -= (CCIS * N1(- CCIS) - N1Density(CCIS)) * vol;  
    }

    return CCEq;
}

/*********************************************************************/

double Black::price(
    bool   isCall, 
    double fwd,
    double strike,
    double pv,
    double variance
) {
    // this is bad
    if (!Maths::isPositive(fwd)) {
        return (isCall ? 0.0 : pv * strike);
    }
    if (!Maths::isPositive(strike)) {
        return (isCall ? pv * (fwd - strike) : 0.0);
    }
    if (!Maths::isPositive(variance)) {
        return pv * (isCall ? Maths::max(fwd-strike, 0.0) : 
                              Maths::max(strike-fwd, 0.0) );
    }
    double sqrtVar = sqrt(variance);
    double d1 = (log(fwd/strike) + (0.5 * variance)) / sqrtVar;
    double d2 = d1 - sqrtVar;
    return pv * (isCall ? (fwd * N1(d1) - strike * N1(d2)) :
                          (strike * N1(-d2) - fwd * N1(-d1)) );
}
                        
double Black::delta(
    bool   isCall, 
    double spot,
    double fwd,
    double strike,
    double pv,
    double variance
){
    // this is bad
    if (!Maths::isPositive(fwd) 
        || !Maths::isPositive(spot)) {
        return 0.0;
    }
    if (!Maths::isPositive(strike)) {
        if (isCall) {
            return (pv * fwd / spot);
        }
        return 0.0;
    }
    if (!Maths::isPositive(variance)) {
        if (isCall) {
            return (Maths::isPositive(fwd - strike) ? pv * fwd / spot : 0.0);
        }
        return (Maths::isPositive(strike - fwd) ? - pv * fwd / spot : 0.0);
    }
    double sqrtVar = sqrt(variance);            
    double d1 = (log(fwd / strike) + (0.5 * variance)) / sqrtVar;
    if (isCall) {
        return (pv * fwd / spot * N1(d1));
    }
    return (- pv * fwd / spot * N1(-d1));
}

double Black::gamma(
    bool   isCall, 
    double spot,
    double fwd,
    double strike,
    double pv,
    double variance
){
    // this is bad
    if (!Maths::isPositive(fwd) 
        || !Maths::isPositive(spot)) {
        return 0.0;
    }
    if (!Maths::isPositive(strike)) {
        if (isCall) {
            return (pv * fwd / spot);
        }
        return 0.0;
    }
    if (!Maths::isPositive(variance)) {
        return 0.0;
    }
    double sqrtVar = sqrt(variance);            
    double d1 = (log(fwd / strike) + (0.5 * variance)) / sqrtVar;

	return N1Density(d1)/(spot * sqrtVar);
}

double Black::vega(
    bool   isCall,
    double fwd,
    double strike,
    double pv,
    double yrsToMat,
    double vol
) {
    static string method = "Black::vega";

    try{
        /* Vega is zero in some limiting cases... */
        if (!Maths::isPositive(fwd)       || 
            !Maths::isPositive(strike)    ||
            !Maths::isPositive(vol)       ||
            !Maths::isPositive(yrsToMat)) {
            return 0.0;
        }

        /* Vega - main case. */
        double d1 = log(fwd / strike) + (0.5 * vol * vol * yrsToMat); 
        d1 /= vol * sqrt(yrsToMat);
    
        double vega  = fwd * N1Density(d1);        
        vega *= pv * sqrt(yrsToMat);
        return vega;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/* Implied var calculation using combined Newton/bisection method.
 * Should be fairly full proof. */
double Black::impliedVariance(
    bool   isCall,
    double fwd,
    double strike,
    double pv,
    double guessVar,     // initial guess
    double targetPrice,
    double epsilon       // i.e., var accuracy
){        
    double impliedVar;
    impliedVariance(isCall, fwd, strike, pv, guessVar, targetPrice,
                    epsilon, true, impliedVar);
    return impliedVar;
}

/** Same as above but does not throw exception if implied variance
    cannot be found for numerical reasons. Instead returns status */
bool Black::impliedVariance(
    bool    isCall,
    double  fwd,
    double  strike,
    double  pv,
    double  guessVar,    /* initial guess */
    double  targetPrice,
    double  epsilon,  /* i.e., var accuracy */
    double& impliedVar
){
    return impliedVariance(isCall, fwd, strike, pv, guessVar, targetPrice,
                           epsilon, false, impliedVar);
}

// used by impliedVariance
struct Black::VegaDiffFunc{
    double weight;
    bool   isCall;
    double fwd;
    double strike;
    double pv;
    double targetPrice;

    double operator()(double var) const{
        double f = Black::price(isCall, fwd, strike, pv, var);

        f -= targetPrice;
        f *= weight;
        return f;
    }

    void operator()(
        double  var,
        double& f,
        double& df
    ) const{
        /* check that variance is always non-negative */
        if (Maths::isNegative(var)) {
            throw ModelException("Balck::vegaDiffWithDeriv", 
                                 "Variance is negative.");
        }

        f = this->operator()(var);

        double stdDeviation = sqrt(var);
        df = Black::vega(isCall,
                         fwd,
                         strike,
                         pv,
                         1.0,                           // years to mat 
                         stdDeviation);         // vega wants the std dev 

        /* we want deriv wrt variance, not std; hence, the "/ (2.0 * stdDeviation)" */
        df *= weight / (2.0 * stdDeviation);
    }
};

bool Black::impliedVariance(
    bool    isCall,
    double  fwd,
    double  strike,
    double  pv,
    double  guessVar,    /* initial guess */
    double  targetPrice,
    double  epsilon,  /* i.e., var accuracy */
    bool    throwOnError,
    double& impliedVar
){
    static string method = "Black::impliedVariance";
    
    try{
        VegaDiffFunc vegaDiff;
        vegaDiff.isCall       = isCall;
        vegaDiff.fwd          = fwd;
        vegaDiff.strike       = strike;
        vegaDiff.pv           = pv;
        vegaDiff.targetPrice  = targetPrice;

        /* check that guessVar is positive */
        if (!Maths::isPositive(guessVar)) {
            if (!throwOnError){
                return false; // failure
            }
            throw ModelException(method, "guessVar is non-positive.");
        }

        /* Computes vega at initial var guess. */
        vegaDiff.weight = vega(isCall,
                               fwd,
                               strike,
                               pv,
                               1.0,         /* assuming 1 year to maturity */
                               sqrt(guessVar));    /* vega wants the std dev */

        double cp = (isCall ? 1.0 : -1.0);
        double intrValue = pv * Maths::max(0.0, cp * (fwd - strike));

        /* Throw an error if option has no time value. */
        if (!Maths::isPositive(targetPrice - intrValue)) {
            if (!throwOnError){
                return false; // failure
            }
            throw ModelException(method, "Option has no time value.");
        }

        /* If option's vega is zero at initial guess */
        if (!Maths::isPositive(vegaDiff.weight)) {
            /* Note: we should rarely pass here, as this case is likely to be
               cleared by case above */
            vegaDiff.weight = 1.0;
        }
        /* otherwise, define weight as 1/vega.  I.e., we minimize/zero
           price diff where diff is expressed in vegas */
        else {
            vegaDiff.weight = 1.0 / vegaDiff.weight;
        }

        /* Need to bracket the root first before we solve for it. */
        double lowerVar = 0.0;         // lower bound for variance
        double upperVar = guessVar;    // upper bound for variance
        if (!ZBracPositive_bracket(vegaDiff,
                                   lowerVar,
                                   upperVar,
                                   throwOnError) &&  !throwOnError){
            return false; // failure
        }

        /* Root is now bracketed. Can look for it. */
        if (!RtSafe_solve(vegaDiff,
                          lowerVar,
                          upperVar,
                          epsilon,
                          throwOnError,
                          impliedVar) && !throwOnError){
            return false; // failure
        }

        if (!Maths::isPositive(impliedVar)) {
            if (!throwOnError){
                return false; // failure
            }

            throw ModelException(method, 
                                 "RtSafe_solve returned negative impliedVar.");
        }

        return true; // success
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double Black::density(
    double fwd,
    double strike,
    double tradYear,
    double vol,
    double dvol_dK,
    double d2vol_dK2
){
    double sqrtTradYear = sqrt(tradYear);
    double sqrtVar = vol * sqrtTradYear;
    double var = sqrtVar * sqrtVar;
    double d1 = (log(fwd/strike) + 0.5 * var) / sqrtVar;
    double d2 = d1 - sqrtVar;
    double densityBS = N1Density(d2) / (strike * sqrtVar);
    double adjustment = 1.0;
    adjustment += 2.0 * strike * sqrtTradYear * d1 * dvol_dK;
    adjustment += strike * strike * tradYear * d1 * d2 * dvol_dK * dvol_dK;
    adjustment += vol * strike * strike * tradYear * d2vol_dK2;
    return (adjustment * densityBS);
}

double Black::probability(
    double  fwd,
    double  strike,
    double  tradYear,
    double  vol,
    double  dvol_dK
){
    double sqrtTradYear = sqrt(tradYear);
    double sqrtVar = vol * sqrtTradYear;
    double var = sqrtVar * sqrtVar;
    double d2 = (log(fwd/strike) - 0.5 * var) / sqrtVar;
    return (N1(-d2) + sqrtTradYear * strike * N1Density(d2) * dvol_dK);
}

/** Black price ADDIN method */
class BlackPriceAddin: public CObject{
    static CClassConstSP const TYPE;

    bool   isCall;
    double fwd;
    double strike;
    double pv;
    double variance;

    static IObjectSP go(BlackPriceAddin* params){
        static const string routine = "BlackPriceAddin::go";
        try {
            double result = Black::price(params->isCall,
                                         params->fwd,
                                         params->strike,
                                         params->pv,
                                         params->variance);            
            CDoubleSP res(CDouble::create(result));
            return IObjectSP(res.clone());
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    BlackPriceAddin():  
    CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BlackPriceAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBlackPriceAddin);
        FIELD(isCall, "isCall");
        FIELD(fwd, "fwd");
        FIELD(strike, "strike");
        FIELD(pv, "pv");
        FIELD(variance, "variance");

        Addin::registerClassObjectMethod("BLACK_PRICE",
                                         Addin::RISK,
                                         "Computes the Black price",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)go);

    }

    static IObject* defaultBlackPriceAddin(){
        return new BlackPriceAddin();
    }
};

CClassConstSP const BlackPriceAddin::TYPE = CClass::registerClassLoadMethod(
    "BlackPriceAddin", typeid(BlackPriceAddin), load);


/** implied variance ADDIN method */
class ImpliedVarianceAddin: public CObject{
    static CClassConstSP const TYPE;

    bool   isCall;
    double fwd;
    double strike;
    double pv;
    double guessVar;
    double targetPrice;
    double epsilon;

    static IObjectSP computeImpliedVariance(ImpliedVarianceAddin* params){
        static const string routine = "ImpliedVarianceAddin::computeImpliedVariance";
        try {
            double result = Black::impliedVariance(params->isCall,
                                                   params->fwd,
                                                   params->strike,
                                                   params->pv,
                                                   params->guessVar,
                                                   params->targetPrice,
                                                   params->epsilon);
            
            CDoubleSP res(CDouble::create(result));
            return IObjectSP(res.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    ImpliedVarianceAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ImpliedVarianceAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultImpliedVarianceAddin);
        FIELD(isCall, "isCall");
        FIELD(fwd, "fwd");
        FIELD(strike, "strike");
        FIELD(pv, "pv");
        FIELD(guessVar, "guessVar");
        FIELD(targetPrice, "targetPrice");
        FIELD(epsilon, "epsilon");

        Addin::registerClassObjectMethod("IMPLIED_VARIANCE",
                                         Addin::RISK,
                                         "Computes the BS Implied Variance",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)computeImpliedVariance);

    }

    static IObject* defaultImpliedVarianceAddin(){
        return new ImpliedVarianceAddin();
    }
};

CClassConstSP const ImpliedVarianceAddin::TYPE = CClass::registerClassLoadMethod(
    "ImpliedVarianceAddin", typeid(ImpliedVarianceAddin), load);

DRLIB_END_NAMESPACE
