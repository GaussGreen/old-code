//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSOption.hpp
//
//   Description : Option on CDS (index or single name). 
//                 Only models implemented so far are European
//                 BS-implied or multi-Q smile. 
//
//   Author      : Charles Morcom
//
//   Date        : February 16, 2006
//
//----------------------------------------------------------------------------

#ifndef QR_CDSOPTION_HPP
#define QR_CDSOPTION_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormBSImpliedSmile.hpp"
#include "edginc/ClosedFormMultiQSmile.hpp"
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CDSOption)
FORWARD_DECLARE(Settlement)
FORWARD_DECLARE(Schedule)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(Control)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE_WRAPPER(YieldCurve)

/**Option on a CDS with fixed rate fee payments. 
 *
 * CALL/PUT DEFINITION
 * A call is the right to buy risk (i.e. sell protection); a put is the right to sell risk 
 * (buy protection).
 *
 * UNDERLYING
 * How the field this->underlying is interpreted depends on the option type (see below)
 * The notional of the underlying is ALWAYS ignored - this->notional should be set as
 * desired and the handling of the notional underlying is controlled by the option type.
 *
 * FIXED/RELATIVE MATURITY UNDERLYING
 * Exercise is into a fixed underlying unless ulMaturityDefinedByOption is true. 
 * This would be the case for, e.g. a single-name american option with a spread-strike,
 * where the maturity was defined as 5Y relative to the exercise date.
 *
 * INDEX VS SINGLE_NAME UNDERLYINGS
 * Index options are always on fixed underlyings, whether spread or price struck.
 * Note the field this->isIndexOption: ideally this should be detected from the
 * underlying's curve; however, currently, index underlyings are not always
 * of an "index" type (e.g. IndexCDSParSpreads) in Pyramid, so this is not
 * a safe thing to try and do.
 * If the strike is a spread, the underlying exercised into is a CDS with the
 * strike spread and a start date the same as the exercise settlement. i.e. there
 * is never any accrued interest.
 *
 * KNOCK-OUT ON DEFAULT
 * If this->koOnDefault is true, the option vanishes on default. This makes a difference
 * to the price of put options, but not usually to call options. Index options may not
 * knock-out on default, and the instrument will throw an exception if you try this.
 *
 * EXERCISE SCHEDULES
 * The schedule must be in strictly ascending date order.
 * For now, only European options are supported, and the instrument will throw an exception
 * if you try and create an option with multiple exercise dates, or with a type not
 * equal to EXERCISE_TYPE_EUROPEAN
 *
 * CASH/PHYSICAL SETTLEMENT
 * Spread-struck index options may not be physically settled: the instrument will throw
 * an exception if you try this. This is because the exercise value is determined by the
 * Bloomberg CDSW flat-curve unwind convention, and this is only meaningful as a cash value.
 *
 * PAY AT DEFAULT
 * If this->payAtDefault is true, then a default triggers an immediate recovery payment
 * (unless koOnDefault!). If the underlying is determined only at the time of exercise
 * (e.g. if this->ulMaturityDefinedByOption or if it is a single-name spread option), the
 * default value of the underlying is assumed to be just (1-R), even though the underlying
 * has not started yet.
 *
 * BINARY/DIGITAL UNDERLYINGS
 * The closed-form pricing models only work for a market recovery underlying CDS.
 *
 * TODO:
 * -Default handling for index options is not implemented: the model throws an exception if
 *  you try to price an index option with defaults. Implementation is pending full default
 *  implementation at the curve level.
 * -Payment delays not implemented in CDSParSpreadCurve, yet
 * -Fix issue in CredDefSwap::generateCDS which causes weird crashes. For now, not used.
 * -CDSW calculations for index forward adjustment need doing.
 */
class PRODUCTS_DLL CDSOption :    public CInstrument, 
                     public virtual ClosedFormMultiQSmile::IIntoProduct,
                     public virtual ClosedFormBSImpliedSmile::IIntoProduct,
                     public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;

    static const string EXERCISE_TYPE_EUROPEAN;
    static const string EXERCISE_TYPE_AMERICAN;
    static const string EXERCISE_TYPE_BERMUDAN;
    static const string STRIKE_TYPE_SPREAD;
    static const string STRIKE_TYPE_CLEAN_PRICE;
    static const string STRIKE_TYPE_DIRTY_PRICE;
    
    virtual ~CDSOption();
    /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    virtual void GetMarket(const IModel* model, const CMarketDataSP);
    virtual void Validate();
    virtual DateTime getValueDate() const;
    virtual string discountYieldCurveName() const;

    /*=========================================================================
     * Pricing model interface implementations
     *=======================================================================*/
    virtual ClosedFormMultiQSmile::IProduct* createProduct(ClosedFormMultiQSmile* model) const;
    virtual ClosedFormBSImpliedSmile::IProduct* createProduct(ClosedFormBSImpliedSmile* model) const;

    /*=========================================================================
     * Tweak interface methods
     *=======================================================================*/
    bool sensShift(Theta* theta);

protected:
    CDSOption();

private:
    CDSOption(const CDSOption& rhs);
    CDSOption& operator=(const CDSOption& rhs);

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Trade value date*/
    DateTime        valueDate;
    /**When is the premium paid, given the trade date? Default is T+3 */
    SettlementSP    premiumSettlement;
    /**Notional of the option/underlying*/
    double          notional;
    /**Start date of option*/
    DateTime        issueDate;
    /**True if knock-out on default*/
    bool            koOnDefault;
    /**True if call (right to buy risk/sell protection, else put.*/
    bool            isCall;
    /**True if the option is on an index underlying name; affects how the forward is adjusted.*/
    bool            isIndexOption;
    /**Exercise type of option: european, american, or bermudan/multi-european.
       Must have one of the following values "EXERCISE_TYPE_EUROPEAN", 
       "EXERCISE_TYPE_BERMUDAN", or "EXERCISE_TYPE_AMERICAN"*/
    string          exerciseType;
    /**Exercise schedule and interpolation type. Note that interpolation type 
       matters only if the option is American, and then it must be INTERP_STAIRS or
       INTERP_LINEAR. A european option with more than one exercise date will be
       treated as if the earliest remaining date is the only one. */
    ScheduleSP      exerciseSchedule;
    /**If true, payments happen at moment of default. 
       If not, then default payment waits until the next exercise after default. */
    bool            payAtDefault;
    /**If true, spread strike, else Price (100%-upfront fee for protection buyer)*/
    string          strikeType;
    /**Underlying CDS.*/
    ICDSSP          underlying; 
    /**Defines underlying settlement relative to exercise. Default is immediate settlement */
    SettlementSP    exerciseSettlement;
    /**If true, the underlying maturity is determined relative to the exercise date. If
       false, the underlying maturity is fixed in underlying.*/
    bool            ulMaturityDefinedByOption;
    /**Underlying maturity (used if ulMaturityDefinedByOption).*/
    ExpirySP        ulExpiry;
    /*=========================================================================
     * IF EXERCISED
     *=======================================================================*/
    /**If true, exercise is cash-settled; if false, physical settlement.*/
    bool            isCashSettled;
    /**True if the option has been exercised.*/
    bool            isExercised; 
    /**Date that option was exercised if (isExercised).*/
    DateTime        exerciseDate;
    /**Market value of the U/L at exercise. If this is a spread strike, this will
       be the 'CDSW' value at the strike when the exercise happened. Note that this
       is _not_ the spread. Note also that it is the dirty PV - be careful with accrued. */
    double          exerciseValue; 
    /*=========================================================================
     * NEED TO ADD FIELDS TO CAPTURE DEFAULTS - DATES, PAYMENTS, STRIKE
     * ADJUSTMENTS, ETC.
     *=======================================================================*/
    
    // these all come from the underlying discount and credit curves and
    // are always overwritten by those in the underlying.
    mutable ICDSParSpreadsSP        crv; 
    mutable YieldCurveSP            rfCrv;
    ICDSParSpreadsWrapper           cdsParSpreads;    
    YieldCurveWrapper               discount;        

    /*=========================================================================
     * FOR REFLECTION
     *=======================================================================*/
    static void load(CClassSP& clazz);
    static IObject* defaultCDSOption();

    /*=========================================================================
     * OTHER METHODS
     *=======================================================================*/
    /**Calculates the next exercise date and strike - returns false if expired, else true.*/
    bool findNextExercise(const DateTime* dateOverride, DateTime* exDt, double* strk) const;
    void priceClosedForm(
        CResults* results, 
        Control* control, 
        const CModel* model) const;

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    friend class CDSOptionHelper;
    friend class CDSOptionClosedFormBS;
    friend class CDSOptionClosedFormMultiQ;
};

DRLIB_END_NAMESPACE
#endif




