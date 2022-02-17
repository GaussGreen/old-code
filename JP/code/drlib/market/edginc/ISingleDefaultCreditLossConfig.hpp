//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : An interface for ICreditLossConfigs that have a single
//                 default date and all the notional goes at the time of
//                 that one default, e.g., single names or NtDs
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ISINGLEDEFAULTCREDITLOSSCONFIG_HPP
#define QLIB_ISINGLEDEFAULTCREDITLOSSCONFIG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/ICreditLossConfig.hpp"


DRLIB_BEGIN_NAMESPACE


class MARKET_DLL ISingleDefaultCreditLossConfig : virtual public ICreditLossConfig
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    virtual ~ISingleDefaultCreditLossConfig();

    /** Returns true if this loss config has defaulted */
    virtual bool hasDefaulted() const = 0;

    /** Returns the default date, or an empty DateTime if the loss config
        has not defaulted */
    virtual DateTime getDefaultDate() const = 0;

    /** Returns the timeline used in this loss config, ie, the timepoints 
        defining the regions at which it is safe to assume that the (flat 
        forwards) default rate is constant - If this concept does not make 
        sense for this loss config, a null SP should be returned */
    virtual DateTimeArraySP getTimeLine() const = 0;


protected:
    ISingleDefaultCreditLossConfig();


private:
    /** Default constructor (for reflection) */
    static IObject* defaultConstructor();

	/** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    ISingleDefaultCreditLossConfig(const ISingleDefaultCreditLossConfig&); // do not use
    ISingleDefaultCreditLossConfig& operator=(const ISingleDefaultCreditLossConfig&); // do not use
};

DECLARE (ISingleDefaultCreditLossConfig);
typedef MarketWrapper<ISingleDefaultCreditLossConfig> ISingleDefaultCreditLossConfigWrapper;


DRLIB_END_NAMESPACE

#endif
