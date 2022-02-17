//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AsMultiFactors.hpp
//
//   Description : AsMultiFactors interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASMULTIFACTORS_HPP
#define EDR_ASMULTIFACTORS_HPP
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
class IMultiFactors;
/** Interface for classes (essentially collections of assets) that
    give a view of themselves as an IMultiFactors object */
class MARKET_DLL IAsMultiFactors {
public:
    IAsMultiFactors(); // in MultiFactors.cpp
    virtual ~IAsMultiFactors();  // in MultiFactors.cpp

    virtual IMultiFactors* asMultiFactors() const = 0;
};

DRLIB_END_NAMESPACE
#endif
