//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : Computes the conditional survival probability at a given 
//               timepoint. Essentially multiplies the "inner names" 
//               conditional survival probabilities - And then integrates
//               them across market factors, to compute the fee leg price.
//               This class is very FtD specific.
//
// Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONDITIONALFTDLOSSINTEGRAND_HPP
#define QLIB_CONDITIONALFTDLOSSINTEGRAND_HPP

#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IDiscountCurve);

/** Computes the  conditional survival probability at a given timepoint.
    Essentially calls the "convolution algorithm" on "inner names" 
    conditional survival probabilities. */
class CONVOLUTION_DLL ConditionalFtDLossIntegrand : public virtual FunctionNDDouble {

public:
    // destructor
    virtual ~ConditionalFtDLossIntegrand();
    
    // constructor
    ConditionalFtDLossIntegrand(
        const vector<ICondLossDistributionsGenKeyArraySP> keysByDate,
        DateTimeArraySP timeLine,
        DateTimeArrayConstSP integrationDates,
        CDoubleArrayConstSP namesLoss,
        int mfDim,
        const DateTime& valueDate,
        const double recoveryDelay,
        IDiscountCurveSP discount);

    // Function
    virtual double operator()(
        const CDoubleArray&  marketFactor) const;

private:
    // Fields
    DateTimeArrayConstSP integrationDates;

    // Array of keys (one per date and per name)
    const vector<ICondLossDistributionsGenKeyArraySP> keysByDate;

    DateTimeArraySP timeLine;

    // Array of "loss given default" (one per name)
    CDoubleArrayConstSP namesLoss;

    const DateTime valueDate;
    const double recoveryDelay;
    IDiscountCurveSP discount;
};

DRLIB_END_NAMESPACE

#endif
