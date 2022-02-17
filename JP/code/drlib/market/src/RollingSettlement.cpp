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

#include "edginc/config.hpp"
#include "edginc/RollingSettlement.hpp"

DRLIB_BEGIN_NAMESPACE


// Create a RollingSettlement
// param period Number of business days to settlement (i.e. T+n)
// param hols Holiday calendar - identifies non-business days

RollingSettlement::RollingSettlement(int period, const HolidayConstSP& hols) :
    Settlement(TYPE),
    period(period), hols(hols.clone()) {
    // empty
}

// Creates a new, instantaneous RollingSettlement 
RollingSettlement::RollingSettlement() : 
    Settlement(TYPE),
    period(0), hols(Holiday::noHolidays()) {
    // empty
}

RollingSettlement::~RollingSettlement(){
}

// Returns the settlement date corresponding to a trade date
DateTime RollingSettlement::settles(const DateTime& tradeDate) const {
    return hols->addBusinessDays(tradeDate, period);
}

// Returns the trade date corresponding to a settlement date
DateTime RollingSettlement::tradeDate(const DateTime& settleDate) const {
    return hols->addBusinessDays(settleDate, -period);
}

/** returns smart pointer to market holiday object */
HolidayConstSP RollingSettlement::getMarketHolidays() const {
    return hols.getSP();
}

/** populate from market cache */
void RollingSettlement::getMarket(
    const IModel* model, const MarketData* market) {
    hols.getData(model, market);
}

class RollingSettlementHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RollingSettlement, clazz);
        SUPERCLASS(Settlement);
        EMPTY_SHELL_METHOD(defaultRollingSettlement);
        FIELD(period, "Offset in days");
        FIELD(hols, "Settlement hols");
    }

    static IObject* defaultRollingSettlement(){
        return new RollingSettlement();
    }
};

CClassConstSP const RollingSettlement::TYPE = CClass::registerClassLoadMethod(
    "RollingSettlement", typeid(RollingSettlement), 
    RollingSettlementHelper::load);

DRLIB_END_NAMESPACE
