//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FixedSettlement.hpp
//
//   Description : Prototype fixed settlement representation
//
//   Author      : Andrew J Swain
//
//   Date        : 16 November 2000
//
//

//
//----------------------------------------------------------------------------

#ifndef FIXEDSETTLEMENT_HPP
#define FIXEDSETTLEMENT_HPP

#include "edginc/Settlement.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE
class MARKET_DLL FixedSettlement : public Settlement {
public:
    static CClassConstSP const TYPE;
    friend class FixedSettlementHelper;

    /** Constructor */
    FixedSettlement(const DateTimeArray& firstTradingDates,
                    const DateTimeArray& settlementDates);        
    
    /** Validation checks */
    virtual void validatePop2Object();

    /** Destructor */
    virtual ~FixedSettlement();

    /** Given a date, return the corresponding settlement date */
    virtual DateTime settles(const DateTime& date) const;

    /** returns smart pointer to market holiday object */
    virtual HolidayConstSP getMarketHolidays() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

     const DateTimeArray& getFirstTradingDates() const { return firstTradingDates; }
     const DateTimeArray& getSettlementDates() const { return settlementDates; }

private:
    FixedSettlement();
    FixedSettlement(const FixedSettlement &rhs);
    FixedSettlement& operator=(const FixedSettlement& rhs);
    // fields
    DateTimeArray  firstTradingDates;
    DateTimeArray  settlementDates;
    HolidayWrapper marketHols;
    // transient fields
    SettlementSP  longDatedSettlement; /* what to use after the settlement runs
                                          out of dates - see comments in
                                          src file for why we need this */
};

DRLIB_END_NAMESPACE

#endif
