//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IBootstrapper.cpp
//
//   Description : Generic interface for objects capable to produce a "state"
//                 for each step of a loop
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IBootstrapper.hpp"

DRLIB_BEGIN_NAMESPACE


/** Virtual destructor */
IBootstrapper::~IBootstrapper() {}

/** Invoked when Class is 'loaded' */
static void IBootstrapperLoad(CClassSP& clazz){
    REGISTER_INTERFACE(IBootstrapper, clazz);
    EXTENDS(IObject);
}

/** TYPE for IBootstrapper */
CClassConstSP const IBootstrapper::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IBootstrapper",
        typeid(IBootstrapper),
        IBootstrapperLoad);

/** TYPE for IBootstrapperArray */
DEFINE_TEMPLATE_TYPE(IBootstrapperArray);

DRLIB_END_NAMESPACE
