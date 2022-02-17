#ifndef _INDEXSPECIR_HPP
#define _INDEXSPECIR_HPP

#include "edginc/IndexSpec.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

const string NamedZeroBankTypeValues();

class MARKET_DLL IndexSpecIR : public IndexSpec
{
public:
    static CClassConstSP const TYPE;

    struct NamedZeroBankType {
        enum Enum { STANDARD, CALLTURFLOWS_T, SWAPTION_T,
            CALLFXRIB_T, EXPLICIT, CLAIMBANK, RIBOBS };
    };
    typedef BoxedEnum<NamedZeroBankType::Enum> NamedZeroBankTypeBoxedEnum;

    /************************ exported fields ************************/
    YieldCurveWrapper factor; // market factor this spec is based on

    MaturityPeriodSP     tenor;
    MaturityPeriodSP     frequency;
    DayCountConventionSP dcc;
    AssetHistoryWrapper  history;// historical values
    // optional fields to specify holidays etc.
    BadDayConventionSP   accrualBadDayConv;
    BadDayConventionSP   paymentBadDayConv;
    ExpirySP             fwdRateOffset;  // allows both period and fixed date
    MaturityPeriodSP     paymentOffset;
    HolidayWrapper       holidays;

    NamedZeroBankType::Enum  zeroBankMethod;

    // these fields allow one to explicity supply zero dates (used when converting
    // existing wrappers and synching critical dates and results
    DateTimeArray        explicitZeroDates;
    int                  nbZerosAtOneTime;

    /*************************** methods ***************************/
    /* IProdCreator:: */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    virtual IMarketFactorConstSP getFactor() const { return factor.getSP(); }
    virtual IMarketFactorSP getFactor() { return factor.getSP(); }

    virtual void validatePop2Object(void);

    IndexSpecIR(const string &name, const YieldCurveWrapper &factor, 
                const AssetHistoryWrapper &history = AssetHistoryWrapper() ) 
        : IndexSpec(TYPE), factor(factor), history(history), zeroBankMethod(NamedZeroBankType::STANDARD), 
                nbZerosAtOneTime(-1) {}

protected:
    IndexSpecIR(const CClassConstSP& type = TYPE) : IndexSpec(type), 
        zeroBankMethod(NamedZeroBankType::STANDARD), nbZerosAtOneTime(-1) {}
    virtual AssetHistoryConstSP getAssetHistory(const string &source) const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IndexSpecIR(); }
};

DECLARE(IndexSpecIR); // declares smart pointers and arrays

DRLIB_END_NAMESPACE

#endif
