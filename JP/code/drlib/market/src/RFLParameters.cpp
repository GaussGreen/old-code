//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : RFLParameters.cpp
//
//   Description : RFLParameters contains 
//                 parameters for RFL and optional CCMParameters.
//                 Represents 'model oriented' parameters for a given name.
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated: Use Cm(Ccm)RflParameters instead
//      This class is now deprecated: Use Cm(Ccm)RflParameters instead
//      All non-essential methods in this class have been removed or throw
//      exceptions (note that instances of this class are automatically
//      converted to the new-style parameters in 
//      CreditEngineParameters::convertToNewParamStyle().
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/RFLParameters.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/PiecewiseFlatMappingFunction.hpp"
#include "edginc/PiecewiseFlatIncrementalMappingFunction.hpp"
#include "edginc/PiecewiseLinearMappingFunction.hpp"
#include "edginc/Nrfns.hpp" //for zbrentUseful root finder
#include "edginc/PortfolioName.hpp" // for the betaTweak function

/* ---------------------------------------------------------------------------
 * Integ Method
 */
#define CCM_INTEG_LOW       -7.5
#define CCM_INTEG_UP         7.5
#define CCM_INTEG_NB         101

DRLIB_BEGIN_NAMESPACE


/** Destructor */
RFLParameters::~RFLParameters() {}

/** 
 * Returns the name of this object
 * */
string RFLParameters::getName() const {
    return name;
}


/** Populates the object with the market data that this object needs */
void RFLParameters::getMarket(const IModel* model, const MarketData* market) {
    if (!ccmParameters.isEmpty()) {
        ccmParameters.getData(model, market);    
    }        
}

/** Returns the CCMParameters of this object */
CCMParametersConstSP RFLParameters::getCCMEngineParams() const {
    return ccmParameters.getSP();
}

/** Invoked when Class is 'loaded' */
void RFLParameters::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RFLParameters, clazz);
    SUPERCLASS(CompositeCreditEngineParameters);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "Identifier");
    FIELD_MAKE_OPTIONAL(name); /* will need to be mandatory when it goes in the
                                  market data cache */

    FIELD(correctMean, "Flag to correct or not the mean");
    FIELD(correctVariance, "Flag to correct or not the variance");
    FIELD(isTweak, "Flag set to true if the beta curve is relative a "
                 "tweak from the historical beta");
    FIELD_MAKE_OPTIONAL(isTweak);
    FIELD(useClosedForm, "Boolean set to true if closed form formulae "
                 "should be used when available (e.g. with piecewise flat "
                 "beta curve)");
    FIELD_MAKE_OPTIONAL(useClosedForm);
    FIELD(   betaCurve,   "Mapping function for beta as a function of market factor");
    FIELD(ccmParameters, "Parameters for CCM model");
    FIELD_MAKE_OPTIONAL(ccmParameters);
}

/** Only build instances of that class using reflection */
RFLParameters::RFLParameters() : CompositeCreditEngineParameters(TYPE), 
                                 isTweak(false), 
                                 useClosedForm(true),
                                 ccmParameters(0)
{}
    
/** Default constructor */
IObject* RFLParameters::defaultConstructor() {
    return new RFLParameters();
}

/** TYPE (for reflection) */        
CClassConstSP const RFLParameters::TYPE =
    CClass::registerClassLoadMethod("RFLParameters",
                                    typeid(RFLParameters),
                                    RFLParameters::load);


DRLIB_END_NAMESPACE
