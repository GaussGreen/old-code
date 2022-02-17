//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RiskyDurationCalculator.hpp
//
//   Description : Holds the information required to calculate the risky
//                 duration when computing the index basis
//
//   Author      : Jose Hilera
//
//   Date        : 30 August 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_RISKYDURATIONCALCULATOR_HPP
#define QLIB_RISKYDURATIONCALCULATOR_HPP

#include "edginc/Object.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICDSParSpreads);
FORWARD_DECLARE(Expiry);
FORWARD_DECLARE(CreditIndexBasis);


/** Class used to calculate the risky durations of the associated CDS Par
 * Spreads (an object of this class is returned by
 * ICDSParSpreads::getRiskyDurationCalculator, to optimize the calculation
 * for that specific type of CDS Par Spreads). */
class MARKET_DLL RiskyDurationCalculator : public virtual VirtualDestructorBase {
public:
    RiskyDurationCalculator();
    virtual ~RiskyDurationCalculator();

    // Method to get the ICDSParSpreadsSP back
    virtual ICDSParSpreadsConstSP getParSpreads() const = 0;

    // Calculate the risky duration for the associated cdsParSpreads,
    // in the calculatingTimeIndex-th expiry
    virtual double riskyDuration(CreditIndexBasisConstSP indexBasis,
                                 int calculatingTimeIndex) const = 0;

    // Validate the durations
    virtual void validateDurations() const = 0;

    // Get name of the curve associated to this RDCalculator
    virtual string getName() const = 0;
};

typedef smartPtr<RiskyDurationCalculator> RiskyDurationCalculatorSP;
typedef vector<RiskyDurationCalculatorSP> RiskyDurationCalculatorArray;
typedef refCountPtr<RiskyDurationCalculatorArray> RiskyDurationCalculatorArraySP;
typedef refCountPtr<const RiskyDurationCalculatorArray> RiskyDurationCalculatorArrayConstSP;

DRLIB_END_NAMESPACE

#endif
