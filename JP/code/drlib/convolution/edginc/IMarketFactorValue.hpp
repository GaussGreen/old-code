//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 03-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IMARKETFACTORVALUE_HPP
#define QLIB_IMARKETFACTORVALUE_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * IMarketFactorValue is a marker interface to represent
 * some "realisation" of the market factor.
 * For example, it could be a path of a process (MC) or just
 * a "point" in a multi-dimensional space.
 * */
class CONVOLUTION_DLL IMarketFactorValue: public virtual VirtualDestructorBase
{
public:
    virtual ~IMarketFactorValue();

protected:
    IMarketFactorValue();

private:    
    IMarketFactorValue(const IMarketFactorValue& rhs); // don't use
    IMarketFactorValue& operator=(const IMarketFactorValue& rhs); // don't use
};

DECLARE(IMarketFactorValue);

DRLIB_END_NAMESPACE

#endif /*QLIB_IMARKETFACTORVALUE_HPP*/
