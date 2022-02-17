//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MappingFunction.hpp
//
//   Description : Mapping function interface
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_MAPPING_FUNCTION_HPP
#define QLIB_MAPPING_FUNCTION_HPP

#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** Mapping function interface */
class MARKET_DLL IMappingFunction: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** The mapping function */
    virtual double map(double x) const = 0;
    
    /** Set a name for this function (useful for tweaking) */
    virtual void setName(string name) = 0;
};

typedef smartPtr<IMappingFunction> IMappingFunctionSP;

DRLIB_END_NAMESPACE

#endif

