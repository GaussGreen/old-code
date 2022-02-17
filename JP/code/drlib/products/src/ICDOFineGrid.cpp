//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOFineGrid.cpp
//
//   Description : Interface for "CDO fine grids" responsible for
//                 "extending" market quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ICDOFineGrid.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

/** Invoked when Class is 'loaded' */
static void ICDOFineGridLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ICDOFineGrid, clazz);
    EXTENDS(IObject);
}

/** TYPE for ICDOFineGrid */
CClassConstSP const ICDOFineGrid::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ICDOFineGrid",
        typeid(ICDOFineGrid),
        ICDOFineGridLoad);

DRLIB_END_NAMESPACE

