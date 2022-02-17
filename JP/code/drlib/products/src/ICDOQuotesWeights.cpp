//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOQuotesWeights.hpp
//
//   Description : ICDOQuotesWeights is an interface for the definition of weights
//                  of each quote in the calibration objective function
//                  (see trancheIndexLEastSquareFit.[c|h]pp)
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICDOQuotesWeights.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(ICDOQuotesWeightsWrapper);

/** Destructor */
ICDOQuotesWeights::~ICDOQuotesWeights() {}

/** Invoked when Class is 'loaded' */
void ICDOQuotesWeights::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ICDOQuotesWeights, clazz);
    EXTENDS(IObject);
}

/** TYPE for CDOQuotesWeights */
CClassConstSP const ICDOQuotesWeights::TYPE = CClass::registerInterfaceLoadMethod(
    "ICDOQuotesWeights", typeid(ICDOQuotesWeights), ICDOQuotesWeights::load);
    
DRLIB_END_NAMESPACE
