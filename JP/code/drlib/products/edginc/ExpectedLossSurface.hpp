//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ExpectedLossSurface.hpp
//
//   Description : Class to store a surface of base expected losses
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/AtomicArray.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DECLARE.hpp"

#ifndef QLIB_EXPECTED_LOSS_SURFACE
#define QLIB_EXPECTED_LOSS_SURFACE

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL ExpectedLossSurface : public CObject 
{
     
    //FIELDS
    static const double TINY; // A small number.

        
    /** date grid points */
    DateTimeArraySP  dates; 
        
    /** Strike grid points */
    DoubleArraySP strikes;

	/** Double Array for expected loss at grid points [strike][time] = [cols][rows] */
    DoubleArrayArraySP losses; 
        
    /** recovery rate */
    double RR;
        

    // METHODS
    /** linear interpolation of exp loss by strike 
		no validation, inline for performance 
		assumes homogeneous strike grid */
    double strikeInterp(double strike, int timeIndex) const;
        
protected:

    ExpectedLossSurface();

public:
        
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultExpectedLossSurface();

    /** Called immediately after object constructed */
    void validatePop2Object();    

    /** CONSTRUCTOR */
    ExpectedLossSurface(DateTimeArraySP dates,
                        DoubleArraySP strikeGrid,
                        DoubleArrayArraySP expLosses,
                        double recoveryRate);

    /** Returns expected loss for base tranche with given upper strike and date (uses linear interpolation) */
    double getBaseEL(double strike,const DateTime & date) const; 

    /** return effective curve (of contingent leg) for given strikes  at timeline points*/
    DoubleArraySP getEffectiveCurve(double lowStrike, double highStrike) const;

    /** return effective curve of fee leg for given strikes.
		This can be different from contingent leg effective curve for super senior teranches.
		We assume constant RR in the calculation */
    DoubleArraySP getFeeLegEffectiveCurve(double lowStrike, double highStrike) const;

    
	/**	Returns effective curve value of tranche [low,high] for feeLeg in index tranche convention 
		assumes constant recovery rate 
		input expected loss is that of the contingent leg 
		baseELLow is expected loss for tranche [0,A] where A = (1-RR)/RR*(1-low) 
		baseELHigh is expected loss for tranche [0,B] where B = (1-RR)/RR*(1-high) */
    static double feeLegEffCurve(
        double expectedLoss,  // unscaled (cont leg) expected loss of tranche [low,high]
        double trancheSize,
        double baseELLow,
        double baseELHigh,
        double RR
        );
        

    /** get the loss surface */
    DoubleArrayArraySP getLosses() const;
        
    /** get the timeline */
    DateTimeArrayConstSP getDates() const;

	/** get the strikes */
	DoubleArrayConstSP getStrikes() const;

	/** get the recovery rate */
	double getRecovery() const;
};

DECLARE(ExpectedLossSurface);


DRLIB_END_NAMESPACE


#endif


