//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICreditLegConvention.cpp
//
//   Description : Interface for "credit legs" conventions
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICreditLegConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Virtual destructor */
ICreditLegConvention::~ICreditLegConvention() {}

/** Invoked when Class is 'loaded' */
static void ICreditLegConventionLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ICreditLegConvention, clazz);
    EXTENDS(IGetMarket);
}

/** TYPE for ICreditLegConvention */
CClassConstSP const ICreditLegConvention::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ICreditLegConvention",
        typeid(ICreditLegConvention),
        ICreditLegConventionLoad);

DRLIB_END_NAMESPACE

