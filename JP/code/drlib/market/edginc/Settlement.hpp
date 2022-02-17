//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Settlement.hpp
//
//   Description : Abstract settlement interface
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef SETTLEMENT_HPP
#define SETTLEMENT_HPP

#include "edginc/DateTime.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/GetMarket.hpp"

DRLIB_BEGIN_NAMESPACE
class YieldCurve;    // forward declaration - avoids #include

// defines an interface for settlement
class MARKET_DLL Settlement: public CObject,
                  virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    virtual ~Settlement();
    
    // return the settlement date for a given date
    virtual DateTime settles(const DateTime& date) const = 0;

    /** returns smart pointer to market holiday object */
    virtual HolidayConstSP getMarketHolidays() const = 0;

    // discount factor between a date & its settlement date
    // default concrete implementation
    virtual double pv(const DateTime& date, YieldCurve *yc) const;
protected:
    Settlement(CClassConstSP clazz);
private:
    Settlement(const Settlement &rhs);
    Settlement& operator=(const Settlement& rhs);
};

typedef smartConstPtr<Settlement> SettlementConstSP;
typedef smartPtr<Settlement> SettlementSP;
#ifndef QLIB_SETTLEMENT_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<Settlement>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<Settlement>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<SettlementSP>(SettlementSP* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<SettlementSP>(SettlementSP* t,
                                                    IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<Settlement>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<Settlement>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<SettlementSP>(SettlementSP* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<SettlementSP>(SettlementSP* t,
                                                         IObjectSP o));
#endif

DRLIB_END_NAMESPACE

#endif
