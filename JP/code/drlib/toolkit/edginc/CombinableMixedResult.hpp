//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CombinableResult.hpp
//
//   Description : Interface for results that can be added together with
//                 objects of a different type
//
//   Author      : Mark A Robson
//
//   Date        : 15 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_COMBINABLE_MIXED_RESULT_HPP
#define EDR_COMBINABLE_MIXED_RESULT_HPP

#include "edginc/CombinableResult.hpp"

DRLIB_BEGIN_NAMESPACE
    
/* Interface for results that can be added together with other types of
   results */
class TOOLKIT_DLL CombinableMixedResult: public CombinableResult{
public:
    static CClassConstSP const TYPE; // defined in CombinableResult.cpp

    /** create a new object by adding an object (scaled by
        scaleFactor) to this result. Implementations should cope with
        the case when x is of different type to this */
    virtual IObject* addResult(const IObject& x, 
                               double scaleFactor) const = 0;
protected:
    CombinableMixedResult(); // defined in CombinableResult.cpp
};


DRLIB_END_NAMESPACE
#endif
