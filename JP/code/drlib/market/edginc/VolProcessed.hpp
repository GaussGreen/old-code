//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessed.hpp
//
//   Description : Abstract Processed Vol Interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLPROCESSED_HPP
#define VOLPROCESSED_HPP

#include "edginc/DateTime.hpp"
#include "edginc/TimeMetric.hpp"

DRLIB_BEGIN_NAMESPACE

/** defines base object used for holding the result of combining instrument
    specific data with the volatility market data */

class MARKET_DLL IVolProcessed: public virtual IObject {
public:
    static CClassConstSP const TYPE;
    /** identifies the market data name of the volatility */
    virtual string getName() const = 0;
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const = 0;
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const = 0;

    virtual ~IVolProcessed();
protected:
    IVolProcessed();
private:
};

typedef smartConstPtr<IVolProcessed> IVolProcessedConstSP;
typedef smartPtr<IVolProcessed> IVolProcessedSP;
#ifndef QLIB_VOLPROCESSED_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IVolProcessed>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IVolProcessed>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IVolProcessed>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IVolProcessed>);
#endif

// for backwards compatibility
typedef IVolProcessed  CVolProcessed;
typedef smartConstPtr<IVolProcessed> CVolProcessedConstSP;
typedef smartPtr<IVolProcessed> CVolProcessedSP;

DRLIB_END_NAMESPACE

#endif
