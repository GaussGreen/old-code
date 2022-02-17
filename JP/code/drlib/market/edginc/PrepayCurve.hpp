//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : PrepayCurve.hpp
//
//   Description : A prepayment curve for ABS
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#ifndef PREPAYCURVE_HPP
#define PREPAYCURVE_HPP

#include "edginc/DecretionCurve.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/PrepayParallelTP.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(Settlement);

class MARKET_DLL PrepayCurve: 
    public DecretionCurve,
    virtual public ITweakableWithRespectTo<PrepayParallelTP> {
public:
    static CClassConstSP const TYPE;
    
    virtual ~PrepayCurve();

    /** checks parameters immediately after object is constructed */
    void validatePop2Object();

    //------------------------------------------
    // IDecretionCurve methods
    //------------------------------------------
    /** pv, which actually is balance here, can be relative to a start date */
    virtual double pv(const DateTime& startDate,
                      const DateTime& endDate) const;
    
    /** pv (balance) relative to initial balance */
    virtual double pv(const DateTime& endDate) const;

    /** return decretion speed on a date */
    virtual double getDecretionSpeed(const DateTime& date) const;

    /**/
    virtual DateTimeArraySP getStepDates() const;
    /** Subtly different from the above - see the code */
    virtual DateTimeArraySP getAlternativeStepDates() const;

    /** return if balances are stepwise or continuous */ 
    virtual bool isStepBalances() const;

    /** Returns a reference to the value date */
    virtual double getFactor(const DateTime& date) const;

    /** Returns settlement rule for prepay */
    virtual SettlementConstSP getSettlement() const;

    //------------------------------------------
    // MarketObject methods
    //------------------------------------------
    bool equalTo(const IObject *obj) const;
    void getMarket(const IModel* model, const MarketData *market);
    int hashCode() const;
    IObject *clone() const;

    //------------------------------------------
    // Tweak methods
    //------------------------------------------
    virtual string sensName(const PrepayParallelTP*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<PrepayParallelTP>&);

protected:
    PrepayCurve(CClassConstSP clazz = TYPE);

    DateTime        today;
    bool            isFlatBalances;

    DateTimeArraySP factorDates;
    DoubleArraySP   factors;

    DateTimeArraySP dates;
    mutable DateTimeArraySP stepDates; // dates on which step wise decretion changes $unregistered
    DoubleArraySP   rates;
    DoubleArraySP   balances;
    DoubleArraySP   speeds;
    SettlementSP    prepaySettle;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultPrepayCurve();
};

typedef smartConstPtr<PrepayCurve> PrepayCurveConstSP;
typedef smartPtr<PrepayCurve>      PrepayCurveSP;
typedef MarketWrapper<PrepayCurve> PrepayCurveWrapper;

DRLIB_END_NAMESPACE
#endif
