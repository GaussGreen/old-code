//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KOStubRule.hpp
//
//   Description   KO stub rule for coupon stream
//
//
//   $Log: KOStubRule.hpp,v $
//----------------------------------------------------------------------------
#ifndef EDR_KO_STUB_RULE_HPP
#define EDR_KO_STUB_RULE_HPP

#include "edginc/config.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Class.hpp"
#include "edginc/InstrumentUtil.hpp"

DRLIB_BEGIN_NAMESPACE

typedef enum{AsSchedule, AsInstSettle} TpayTimingRule;
class KOSettle;
FORWARD_DECLARE(KOStubRule)

class PRODUCTS_DLL KOStubRule: public CObject {
public:
    static CClassConstSP const TYPE;
    friend class KOSettle;

    KOStubRule();

    KOStubRule(const string koStubRule,
        const bool isAccrueUpToSettle,
        const string payTimingRule);

    void validatePop2Object();

    KOSettle* makeKOSettle(KOStubRuleSP koRule,
                           const DateTimeArray& payDates,
                           const DateTimeArray& accrueDates);
    
    KOSettle* makeKOSettle(KOStubRuleSP koRule,
                           const DateTimeArray& payDates,
                           const DateTimeArray& accrueDates,
                           InstrumentSettlementSP settle);
    //
    bool isInstSettle();

    string getKOStubRule();    
    bool   getIsAccrueUpToSettle();
    string getpayTimingRule();
    
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    static IObject* defaultKOStubRule();

    string koStubRule;
    bool   isAccrueUpToSettle;
    string payTimingRule;

    TpayTimingRule payTimeRule; // $unregistered
};

typedef smartPtr<KOStubRule> KOStubRuleSP;

//////////////////////////////
// KO Settle Class          //
//////////////////////////////
class PRODUCTS_DLL KOSettle {
public:
    
    // constructor
    KOSettle(KOStubRuleSP            koRule,
             const InstrumentSettlementSP  settle,
             const DateTimeArray           payDates,
             const DateTimeArray           accrueDates,
             const bool                    isRollingSettle);

	KOSettle(KOStubRuleSP        koRule,
             const DateTimeArray           payDates,
             const DateTimeArray           accrueDates);


    KOSettle(KOStubRuleSP            koRule,
             const InstrumentSettlementSP  settle,
             const DateTimeArray           payDates,
             const DateTimeArray           accrueDates);

    bool isEmpty() const;

    // return the settleDate.  Those data should be set up by product class.
    bool getSettleDate(const DateTime hitDate, DateTime& settleDate);

    bool getAccrueUpToSettle();
    
    DateTimeArray getAccrueDates();

    DateTime getLastAccrueDate();

    DateTime getFirstAccrueDate();

    string	 getKOStubRule();

    // is instrument settlement given?
	bool	 hasSettle();

private:
    KOStubRuleSP            koRule;
    InstrumentSettlementSP  settle;                 //settle object (Not used when isRollingSettle = false)
    bool                    isRollingSettle;        
    
    DateTimeArray           payDates;               //payment Schedule.  
    DateTimeArray           accrueDates;        //accrue end dates.

};

typedef refCountPtr<KOSettle> KOSettleSP;
                

DRLIB_END_NAMESPACE

#endif
