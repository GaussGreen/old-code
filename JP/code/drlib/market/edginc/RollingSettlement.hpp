//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RollingSettlement.hpp
//
//   Description : Prototype rolling settlement representation
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef ROLLINGSETTLEMENT_HPP
#define ROLLINGSETTLEMENT_HPP

#include "edginc/Settlement.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/smartPtr.hpp"

DRLIB_BEGIN_NAMESPACE
class MARKET_DLL RollingSettlement : public Settlement {
public:
    static CClassConstSP const TYPE;
    friend class RollingSettlementHelper;

    // Create a RollingSettlement
    // param period Number of business days to settlement (i.e. T+n)
    // param hols Holiday calendar - identifies non-business days

    RollingSettlement(int period, const HolidayConstSP& hols);

    // Creates a new, instantaneous RollingSettlement
    RollingSettlement();

    virtual ~RollingSettlement();

    // Returns the settlement date corresponding to a trade date
    virtual DateTime settles(const DateTime& tradeDate) const;

    // Returns the trade date corresponding to a settlement date
    DateTime tradeDate(const DateTime& settleDate) const;

    /** returns smart pointer to market holiday object */
    virtual HolidayConstSP getMarketHolidays() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

     int getPeriod() const { return period; }

protected:
    RollingSettlement(const RollingSettlement &rhs);
    RollingSettlement& operator=(const RollingSettlement& rhs);
    int            period;
    HolidayWrapper hols;
};

typedef smartConstPtr<RollingSettlement> RollingSettlementConstSP;
typedef smartPtr<RollingSettlement> RollingSettlementSP;

DRLIB_END_NAMESPACE

#endif
