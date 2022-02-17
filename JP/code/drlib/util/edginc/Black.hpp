//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Black.hpp
//
//   Description : Black-Scholes forumla
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef BLACK_HPP
#define BLACK_HPP

DRLIB_BEGIN_NAMESPACE

/** Set of methods used for calculating option prices using Black-Scholes
    methodology */
class UTIL_DLL Black {
public:

    /** Price of a call or put using 2Q version of Black&Scholes.
        Note: vol is defined by Y = Y/(1+fwdSh) * F((1+fwdSh)/(1+q*fwdSh)*vol) */
	static double price2Q(
        bool   isCall,
        double fwd,      /* Forward yield */
        double strike,
        double pv,
        double variance,
        double qLeft,
        double qRight,
        double fwdSh     /* Forward shift */
    );


    /** Calculates the price for a spot starting european vanilla
        option, given the forward at maturity, the strike, the
        discount factor between today and maturity, and the variance
        between today and maturity */
    static double price(
        bool   isCall, 
        double fwd,
        double strike,
        double pv,
        double variance
    );

    /** Calculates the delta for a spot starting european vanilla
        option, given the forward at maturity, the strike, the
        discount factor between today and maturity, and the variance
        between today and maturity */
    static double delta(
        bool   isCall, 
        double spot,
        double fwd,
        double strike,
        double pv,
        double variance
    );

	/** Calculates the gamma for a spot starting european vanilla
        option, given the forward at maturity, the strike, the
        discount factor between today and maturity, and the variance
        between today and maturity */
    static double gamma(
        bool   isCall, 
        double spot,
        double fwd,
        double strike,
        double pv,
        double variance
    );

    /** Calculates the vega for a spot starting european vanilla
        option, given the forward at maturity, the strike, the
        discount factor between today and maturity, and the variance
        between today and maturity */
    static double vega(
        bool   isCall,
        double fwd,
        double strike,
        double pv,
        double yrsToMat,
        double vol
    );

    /** Implied var calculation using combined Newton/bisection method.
        Should be fairly full proof. */
    static double impliedVariance(
        bool   isCall,
        double fwd,
        double strike,
        double pv,
        double guessVar,    /* initial guess */
        double targetPrice,
        double epsilon      /* i.e., var accuracy */
    );

    /** Same as above but does not throw exception if implied variance
        cannot be found for numerical reasons. Instead returns status 
        true = success */
    static bool impliedVariance(
        bool    isCall,
        double  fwd,
        double  strike,
        double  pv,
        double  guessVar,    /* initial guess */
        double  targetPrice,
        double  epsilon,     /* i.e., var accuracy */
        double& impliedVar
    );  

    /** Given an implied volatility and its first and second strike derivatives
        at a given strike, calculates the density of the spot price at that strike */
    static double density(
        double fwd,
        double strike,
        double tradYear,
        double vol,
        double dvol_dK,
        double d2vol_dK2
    );

    /** Given an implied volatility and its first derivative at a given strike, 
        calculates the probability mass of the spot price up to that strike */
    static double probability(
        double fwd,
        double strike,
        double tradYear,
        double vol,
        double dvol_dK
    );

private:

    /** begin of price2Q private functions */
    static double convexity2Q(
        double vol,    /* Annualized volatility */
        double qLeft,  /* Q left                */
        double qRight, /* Q right               */
        double fwdSh   /* Fwd shift             */
    );
    static double equation2Q(
        double vol,    /* Annualized volatility */
        double qLeft,  /* Q left                */
        double qRight, /* Q right               */
        double fwdSh,  /* Fwd shift             */
        double CC      /* Constant              */
    );
    /** end of price2Q private functions */

    static bool impliedVariance(
        bool    isCall,
        double  fwd,
        double  strike,
        double  pv,
        double  guessVar,    /* initial guess */
        double  targetPrice,
        double  epsilon,     /* i.e., var accuracy */
        bool    throwOnError,
        double& impliedVar
    );  

    // used by impliedVariance
    struct VegaDiffFunc;
};

DRLIB_END_NAMESPACE

#endif
