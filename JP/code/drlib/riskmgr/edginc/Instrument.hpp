
#ifndef EDG_INSTRUMENT_H
#define EDG_INSTRUMENT_H

#include "edginc/DateTime.hpp"
#include "edginc/SensControl_forward.hpp"
#include "edginc/MarketData_forward.hpp"
#include "edginc/Instrument_forward.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE
class IModel;
//class CTree1fProduct;
//class CTree1f;

/** Interface that instruments implement. Note that, in general, most of
    the infrastucture has not been changed to use this type and instead still
    uses CInstrument */
class RISKMGR_DLL IInstrument: public virtual IObject{
public:
    static CClassConstSP const TYPE;

    ~IInstrument();
    IInstrument();

    /** instrument GetMarket to initiate market data selection.
        all it needs to do usually is to specify data tags for the market data
        this instrument depends on and then call data select in model.
        it allows instrument to implement specific get market actions */
    virtual void GetMarket(const IModel*, const CMarketDataSP) = 0 ;

    /** Called once before the initial pricing */
    virtual void Validate() = 0;

    /** override a control shift (eg for delta on trees) - may return
        null to use original. */
    virtual CSensControl* AlterControl(
        const IModel*          modelParams,
        const CSensControl*    sensControl) const = 0;

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const = 0;

    /** price a dead instrument until settlement - exercised, expired,
        knocked out etc.  returns true if it is dead (and priced), false
        if it is not dead */
    virtual bool priceDeadInstrument(CControl* control,
                                     CResults* results) const = 0;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const = 0;

private:
    static void load(CClassSP& clazz);
};


/** A base class instruments. Actual instruments will, in addition, need
    to implement an interface specified by one or more CModel classes */
class RISKMGR_DLL CInstrument: public CObject,
                   public virtual IInstrument{
public:
    friend class InstrumentHelper;
    static CClassConstSP const TYPE;

    virtual ~CInstrument();

    /** override a control shift (eg for delta on trees) - may return
        null to use original. Default implementation returns null. */
    virtual CSensControl* AlterControl(
        const IModel*          modelParams,
        const CSensControl*    sensControl) const;

    /** price a dead instrument until settlement - exercised, expired,
        knocked out etc.  returns true if it is dead (and priced), false
        if it is not dead */
    virtual bool priceDeadInstrument(CControl* control,
                                     CResults* results) const;

protected:
    CInstrument(CClassConstSP clazz);
};

DRLIB_END_NAMESPACE
#endif
