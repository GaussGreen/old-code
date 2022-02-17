//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MixedSettlement.cpp
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

#include "edginc/config.hpp"
#include "edginc/MixedSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor : settle1 and settle2 can be any combination of
    'Fixed',  'Rolling' or 'Mixed' again */
MixedSettlement::MixedSettlement(Settlement*       settle1,
                                 Settlement*       settle2,
                                 const DateTime&   switchDate)
    :Settlement(TYPE), switchDate(switchDate)
{
    // destructor is not called if we throw an exception, so copy components
    // before in case of failure 
    SettlementSP   settle1Copy(copy(settle1));
    SettlementSP   settle2Copy(copy(settle2));
    this->settle1 = settle1Copy;
    this->settle2 = settle2Copy; 
}

MixedSettlement::~MixedSettlement() 
{
    // empty
}

/** For dates less than the switch date settle1 is used.
    For dates >= the switch date, settle2 is used. */
DateTime MixedSettlement::settles(const DateTime& date) const 
{
    // select which settlement to use based on the switch date
    if (date.isLess(switchDate))
    {
        return settle1->settles(date);
    }
    else
    {
        return settle2->settles(date);
    }
    
}

HolidayConstSP MixedSettlement::getMarketHolidays() const {
   return settle1->getMarketHolidays();
}

/** populate from market cache */
void MixedSettlement::getMarket(
    const IModel* model, const MarketData* market) {
    settle1->getMarket(model, market);
    settle2->getMarket(model, market);
}

// for reflection
MixedSettlement::MixedSettlement(): Settlement(TYPE){}

class MixedSettlementHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MixedSettlement, clazz);
        SUPERCLASS(Settlement);
        EMPTY_SHELL_METHOD(defaultMixedSettlement);
        FIELD(settle1, "settlement before switch date");
        FIELD(settle2, "settlement on or after switch date");
        FIELD(switchDate, "switch over date");
    }

    static IObject* defaultMixedSettlement(){
        return new MixedSettlement();
    }
};

CClassConstSP const MixedSettlement::TYPE = CClass::registerClassLoadMethod(
    "MixedSettlement", typeid(MixedSettlement), MixedSettlementHelper::load);

DRLIB_END_NAMESPACE
