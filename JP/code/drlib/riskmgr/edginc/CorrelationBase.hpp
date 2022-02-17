//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationBase.hpp 
//
//   Description : CorrelationBase  - allows models to manipulate correlations
//
//   Author      : Mark A Robson
//
//   Date        : 19 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CORRELATIONBASE_HPP
#define EDR_CORRELATIONBASE_HPP
#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE
class IPerNameSensitivity;

/** Interface for correlations - the main idea being that correlations between
    different factors can often reuse the same structures often it is required
    that they produce different greeks */
class RISKMGR_DLL CorrelationBase: public MarketObject {
public:
    static CClassConstSP const TYPE;

    virtual ~CorrelationBase();

    /** Configure this correlation object so that under tweaking it behaves
        properly when it is a correlation between object of type clazz1 and
        an object of type clazz2 */
    virtual void configureForSensitivities(CClassConstSP clazz1,
                                           CClassConstSP clazz2) = 0;

    /** Returns true if this correlation is [really] sensitive to the
        supplied sensitivity */
    virtual bool isSensitiveTo(const IPerNameSensitivity* sens) const = 0;
protected:
    CorrelationBase(CClassConstSP clazz);
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<CorrelationBase> CorrelationBaseSP;
#ifndef QLIB_CORRELATIONBASE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CorrelationBase>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CorrelationBase>);
#endif

DRLIB_END_NAMESPACE
#endif
