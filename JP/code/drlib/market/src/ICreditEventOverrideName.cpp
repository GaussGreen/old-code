//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : ICreditEventOverrideName.cpp
//
//   Description : Interface used by credit instruments in order to override
//                 (at instrument level) a name's default-related parameters.
//                 Due to design changes, this interface now derives from
//                 ITrancheCreditEventOverride and defines just one new 
//                 method, so conceptually it is almost identical to 
//                 ITrancheCreditEventOverride.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/IBadDayAdjuster.hpp"

DRLIB_BEGIN_NAMESPACE

ICreditEventOverrideName::~ICreditEventOverrideName()
{}


void ICreditEventOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(ICreditEventOverrideName, clazz);
    EXTENDS(ITrancheCreditEventOverride);
}


CClassConstSP const ICreditEventOverrideName::TYPE = 
    CClass::registerInterfaceLoadMethod("ICreditEventOverrideName", 
                                        typeid(ICreditEventOverrideName), 
                                        load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(ICreditEventOverrideNameArray);

DRLIB_END_NAMESPACE
