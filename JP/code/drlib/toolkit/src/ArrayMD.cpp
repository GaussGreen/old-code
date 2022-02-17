//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ArrayMD.cpp
//
//   Description : Multidimensional Array
//
//   Date        : Sep 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_ARRAYMD_CPP
#include "edginc/ArrayMD.hpp"
#ifdef DEBUG
// so we instantiate the template
#include "edginc/Lattice.hpp"
#endif
DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL ArrayMD<double>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DoubleArrayMD>);


template<> CClassConstSP const DoubleArrayMD::TYPE = CClass::registerClassLoadMethod(
    "DoubleArrayMD", typeid(DoubleArrayMD), load);

DRLIB_END_NAMESPACE

