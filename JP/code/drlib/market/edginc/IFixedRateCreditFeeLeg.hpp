//----------------------------------------------------------------------------
//
//   Filename    : IFixedRateCreditFeeLeg.hpp
//
//   Description : Interface to describe a credit fee leg which is
//                 constrained to have a constant fixed fee
//
//----------------------------------------------------------------------------

#ifndef QR_IFIXEDRATECREDITFEELEG_HPP
#define QR_IFIXEDRATECREDITFEELEG_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IFixedRateCreditFeeLeg : public virtual IObject
{
public:

    /**Get the fee rate of the fee leg.*/
    virtual double getRate() const = 0;
    
    /**Change the fee rate of the fee leg. Note that, if the rate is currently zero,
       this will throw an exception.*/
    virtual void setRate(double newRate) = 0;

    static CClassConstSP const TYPE;
    IFixedRateCreditFeeLeg();
    virtual ~IFixedRateCreditFeeLeg();
private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
