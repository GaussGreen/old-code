//   Filename    : MixedSettlement.hpp
//
//   Description : Prototype mixed settlement representation
//
//   Author      : Stephen Hope
//
//   Date        : 25 Jan 2001
//
//
//   $Log:
//----------------------------------------------------------------------------

#ifndef MIXEDSETTLEMENT_HPP
#define MIXEDSETTLEMENT_HPP

#include "edginc/Settlement.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

/** A mixed settlement is made up of two other settlements which
    can be any combination of 'Fixed', 'Rolling' or 'Mixed' again */
class MARKET_DLL MixedSettlement : public Settlement {
public:
    static CClassConstSP const TYPE;
    friend class MixedSettlementHelper;
    
    /** Constructor : settle1 and settle2 can be any combination of
        'Fixed',  'Rolling' or 'Mixed' again */
    MixedSettlement(Settlement*       settle1,
                    Settlement*       settle2,
                    const DateTime&   switchDate);        
    
    virtual ~MixedSettlement();

    /** Given a date, return the corresponding settlement date.
        For dates less than the switch date settle1 is used.
        For dates >= the switch date, settle2 is used. */
    virtual DateTime settles(const DateTime& date) const;

    virtual HolidayConstSP getMarketHolidays() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

private:
    MixedSettlement();
    MixedSettlement(const MixedSettlement &rhs);
    MixedSettlement& operator=(const MixedSettlement& rhs);
    SettlementSP  settle1;
    SettlementSP  settle2;
    DateTime      switchDate;
};

DRLIB_END_NAMESPACE

#endif
