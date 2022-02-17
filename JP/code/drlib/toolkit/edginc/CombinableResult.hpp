//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CombinableResult.hpp
//
//   Description : Interface for results that can be added together
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef COMBINABLERESULT_HPP
#define COMBINABLERESULT_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for results that can be scaled and added together */
class TOOLKIT_DLL CombinableResult: virtual public IObject {
public:
    virtual ~CombinableResult();

    static CClassConstSP const TYPE;
    /** scale by factor x */
    virtual void scale(double x) = 0;

    /** add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor) = 0;

protected:
    CombinableResult();
};

DRLIB_END_NAMESPACE
#endif
