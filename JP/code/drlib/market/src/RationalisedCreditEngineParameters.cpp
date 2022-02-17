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

#include "edginc/config.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"

DRLIB_BEGIN_NAMESPACE


RationalisedCreditEngineParameters::~RationalisedCreditEngineParameters() 
{}

RationalisedCreditEngineParameters::RationalisedCreditEngineParameters(
        const CClassConstSP& clazz) :
    CreditEngineParameters(clazz)
{}


/** Invoked when Class is 'loaded' */
void RationalisedCreditEngineParameters::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RationalisedCreditEngineParameters, clazz);
    SUPERCLASS(CreditEngineParameters);
}

CClassConstSP const RationalisedCreditEngineParameters::TYPE = 
    CClass::registerClassLoadMethod("RationalisedCreditEngineParameters", 
                                    typeid(RationalisedCreditEngineParameters), 
                                    RationalisedCreditEngineParameters::load);

DEFINE_TEMPLATE_TYPE(RationalisedCreditEngineParametersWrapper);

// initialise type for array of RationalisedCreditEngineParameters
DEFINE_TEMPLATE_TYPE(RationalisedCreditEngineParametersArray);

DRLIB_END_NAMESPACE
