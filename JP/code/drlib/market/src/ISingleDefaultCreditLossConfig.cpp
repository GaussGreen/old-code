//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : An interface for ICreditLossConfigs that have a single
//                 default date and all the notional goes at the time of
//                 that one default, e.g., single names or NtDs
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ISingleDefaultCreditLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

ISingleDefaultCreditLossConfig::~ISingleDefaultCreditLossConfig() 
{}

ISingleDefaultCreditLossConfig::ISingleDefaultCreditLossConfig() 
{}

static void ISingleDefaultCreditLossConfigLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(ISingleDefaultCreditLossConfig, clazz);
    EXTENDS(ICreditLossConfig);
    clazz->setPublic(); // make visible to EAS/spreadsheet, and allow creation of arrays
}


CClassConstSP const ISingleDefaultCreditLossConfig::TYPE = 
    CClass::registerInterfaceLoadMethod("ISingleDefaultCreditLossConfig", 
                                        typeid(ISingleDefaultCreditLossConfig), 
                                        ISingleDefaultCreditLossConfigLoad);

DEFINE_TEMPLATE_TYPE(ISingleDefaultCreditLossConfigWrapper);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(ISingleDefaultCreditLossConfigArray);

DRLIB_END_NAMESPACE
