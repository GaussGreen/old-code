//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Abstract base class for per name engine specific 
//                 credit parameters.
//                 It is identical to CreditEngineParameters but provides a 
//                 different (sub)type. All "new-style" credit engine
//                 parameters derive from it - therefore allowing us to specify
//                 that, eg, only new-style parameters are allowed in certain
//                 places.
//
//   Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RATIONALISEDCREDITENGINEPARAMETERS_HPP
#define QLIB_RATIONALISEDCREDITENGINEPARAMETERS_HPP

#include "edginc/CreditEngineParameters.hpp"

DRLIB_BEGIN_NAMESPACE


class RationalisedCreditEngineParameters;
typedef smartPtr<RationalisedCreditEngineParameters> RationalisedCreditEngineParametersSP;
typedef smartConstPtr<RationalisedCreditEngineParameters> RationalisedCreditEngineParametersConstSP;

// support for wrapper class
typedef MarketWrapper<RationalisedCreditEngineParameters> RationalisedCreditEngineParametersWrapper;
typedef smartPtr<RationalisedCreditEngineParametersWrapper> RationalisedCreditEngineParametersWrapperSP;


/** Abstract base class for per name engine specific credit parameters*/
class MARKET_DLL RationalisedCreditEngineParameters: public CreditEngineParameters {

public:
    static CClassConstSP const TYPE;

    virtual ~RationalisedCreditEngineParameters();
    
protected:
    RationalisedCreditEngineParameters(const CClassConstSP& clazz);

private:
    RationalisedCreditEngineParameters(
        const RationalisedCreditEngineParameters& rhs); // don't use
    RationalisedCreditEngineParameters& operator=(
        const RationalisedCreditEngineParameters& rhs); // don't use
    static void load(CClassSP& clazz);
};

// support for arrays
typedef array<RationalisedCreditEngineParametersSP, 
              RationalisedCreditEngineParameters> RationalisedCreditEngineParametersArray;

DRLIB_END_NAMESPACE

#endif
