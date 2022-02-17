//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashSettlePeriod.hpp
//
//   Description : Cash settlement on a rolling basis
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef CASHSETTLEPERIOD_HPP
#define CASHSETTLEPERIOD_HPP

#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/HolidayCollector.hpp"

DRLIB_BEGIN_NAMESPACE

/** Cash settlement on a rolling basis */

class MARKET_DLL CashSettlePeriod : public InstrumentSettlement {
public:
    static CClassConstSP const TYPE;

    CashSettlePeriod(int period, const Holiday* hols);
    CashSettlePeriod(int period);

    virtual ~CashSettlePeriod();
    
    /** given a trade date, when does this settle ? */
    virtual DateTime settles(const DateTime& date,
                             const CAsset*   asset) const;
   
    /** is this a physical settlement ? */
    virtual bool isPhysical() const;

    /** is this a margin option ? */
    virtual bool isMargin() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

private:
    /** Holiday collector - rather than a get method
        cos not all settlement types have holidays */
    static void acceptHoliday(const CashSettlePeriod* settle,
                       HolidayCollector*   collector);
    
    CashSettlePeriod();  // for reflection
    friend class CashSettlePeriodHelper;

    int            period;
    HolidayWrapper hols;    
};
typedef smartPtr<CashSettlePeriod> CashSettlePeriodSP;
DRLIB_END_NAMESPACE

#endif
