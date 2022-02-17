//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : QuantoCDSParSpreads.hpp
//
//   Description : Holds the par spreads in a different currency to the name's
//                 native currency (eg IBM par spreads in euros)
//
//   Author      : Mark A Robson
//
//   Date        : November 26, 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_QUANTO_CDSPARSPREADS_HPP
#define EDR_QUANTO_CDSPARSPREADS_HPP

#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/Duration.hpp"
#include "edginc/VolProcessedBSIR.hpp"

DRLIB_BEGIN_NAMESPACE

/** Holds the par spreads in a different currency to the name's
    native currency (eg IBM par spreads in euros) */
class MARKET_DLL QuantoCDSParSpreads: public MarketObject,
                           public virtual ICDSParSpreads,
                           public virtual Duration::IParHandlerWithClosedForm,
                           public virtual ParSpreadMaxTweakSize::IShift {
public:
    /** Models must be able to return an instant of this interface in order
        for the quanto'd curve to be bootstrapped */
    class MARKET_DLL IAlgorithm : public virtual VirtualDestructorBase {

        public:

        IAlgorithm();
        virtual ~IAlgorithm();

        /** Return clean spreads in the new currency
            (ie fx base/dom currency) */
        virtual DefaultRatesConstSP currencyAdjust(
            IYieldCurveConstSP    domYC,   // set to Projection curve
            IYieldCurveConstSP    forYC,   // set to Projection curve
            ICDSParSpreadsConstSP           parSpreads,
            FXAssetConstSP                  fx,
            CorrelationConstSP              corrFxCDS, // fx and cds par spreads
            CorrelationConstSP              corrFxDom, // fx and dom IR
            CorrelationConstSP              corrFxFor, // fx and for IR
            CorrelationConstSP              corrCDSDom,// cds and domestic iR
            CorrelationConstSP              corrCDSFor,// cds and foreign IR
            CorrelationConstSP              corrDomFor) const = 0;// dom and for

        //-------------------------------------
        // IR Grid point cache related methods
        //-------------------------------------

        /** Cache IR grid points */
        virtual void cacheGridPoints(const IYieldCurveConstSP domYC,
                                     const IYieldCurveConstSP forYC) const = 0;

        //------------------------------
        // Cache-related methods
        //------------------------------

        /** Hash code function - the cache needs improved performance compared
         * to the default "hashCode" function in CObjet: only the required
         * components are hashed */
        virtual int hashCodeOpt() const = 0;

        /** Comparison function - the cache needs improved performance compared
         * to the default "equalTo" function in CObjet: only the required
         * components are compared */
        virtual bool equalToOpt(const IAlgorithm* algorithm) const = 0;
    };
    typedef smartPtr<IAlgorithm> IAlgorithmSP;
    typedef smartConstPtr<IAlgorithm> IAlgorithmConstSP;

    /** Interface that models must implement if they want to support quanto'd
        CDS Par Spread Curves. Model implementations can depend upon the
        behaviour that the QuantoCDSParSpreads will not cache the
        IAlgorithmBuilder but rather the IAlgorithm instead. This IAlgorithm
        will be used across tweaks etc */
    class MARKET_DLL IAlgorithmBuilder: public virtual IObject{
    public:
        //// interface is typed so we can use reflection information to see
        //// which models support this feature
        static CClassConstSP const TYPE;

        virtual ~IAlgorithmBuilder();

        /** Returns an instance of IAlgorithm. Typically the quantoCDSParSpreads
            parameter would be ignored but is there in case you want to switch
            the algorithm dependent upon some property of the quanto'd curve */
        virtual IAlgorithmSP cdsQuantoAlgorithm(
            const QuantoCDSParSpreads* quantoCDSParSpreads) const = 0;

        /** sets debug state to specified value. If true, then any calls to
            cdsQuantoAlgorithm will return an object that will cache debug data.
            This can be retrieved via getDebugInfo() */
        virtual void selectDebugState(bool switchOn) = 0;
        /** Returns debug info - may be null eg if the algorithm has not been
            used */
        virtual IObjectSP getDebugInfo() const = 0;
    private:
        static void load(CClassSP& clazz);
    };

    static CClassConstSP const TYPE;
    ~QuantoCDSParSpreads();

    /** Override CObject default to preserve the IAlgorithm and cached data */
    virtual IObject* clone() const;

    /** Override CObject::fieldsUpdated() to invalidate cached data */
    virtual void fieldsUpdated(const CFieldArray& fields);

    /** Returns the name of this object.*/
    virtual string getName() const;

    virtual string getParSpreadsName() const;

    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Quanto'd par spreads, computed from defaultRates() */
    ParSpreadCurveConstSP getParSpreadCurve() const;

    /** Quanto'd par spreads, computed from defaultRates() */
    DoubleArrayConstSP getParSpreads() const;

    /** Allows a shift to multiple points on the curve simultaneously */
    // tenorShifts is expected to be the same size as the expiries/spreads
    // shifts are additive
    // modifies the object
    // This method is outside the usual sensitivity framework
    virtual void bucketShift(const DoubleArray& tenorShifts);

    /** Expiries (relative: eg. 1M) */
    ExpiryArrayConstSP getParSpreadsExpiries() const;

    //// Returns the clean spreads in the new currency
    virtual DefaultRatesSP defaultRates() const;

    /** return the bootstrapped dates */
    virtual DateTimeArray zeroDates() const;

    // accessor methods for logOfDiscFactorKey
    virtual IDiscountCurve::IKey* getDiscountKey() const;
    virtual DefaultRates::IKey* getRiskyKey() const;

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const;

    /** Returns the recovery rate */
    virtual double getRecovery() const;

    /** Returns true if the name has defaulted */
    virtual bool defaulted() const;

    /** Returns the date that this name defaulted. If a default has not
        happened an 'empty' DateTime is returned */
    virtual const DateTime& getDefaultDate() const;

    /** return CDS Par Spreads currency [isocode] */
    virtual string getCcy() const;

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const;

    /** This currently fails as it's not clear what the right behaviour is.
        One possibility is to just use the vol of the curve that was quanto'd */

    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;

    /** Returns the name of the curve that was quanto'd - not ideal but then
        a) we currently haven't a way of specifying such a correlation and
        b) you probably want the model to consider the separate factors
        anyway */
    virtual string getCorrelationName() const;

    /** Returns the date when to stop tweaking after lastDate */
    virtual const DateTime stopTweaking(const DateTime& lastDate) const;

    /** Quanto'd spread at a particular date, interpolated if necessary
        (part of ICreditCurve interface) */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate) const;

    /** Quanto'd spread at a particular date, interpolated if necessary
        (part of ICreditCurve interface) */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate,
                                    const BadDayConvention* bdc,
                                    const Holiday* hols) const;

    /** Calculates the maximum tweak sizes for this ICDSParSpreads.
        Currently just passes this down to original (ie unquantoed) curve */
    virtual DoubleArrayConstSP calculateTweakSizes() const;

    /** Returns the day count convention used for par CDS */
    virtual DayCountConvention* dayCountConv() const;

    /** Returns holidays used to adjust cash flow dates */
    virtual HolidayConstSP getHolidays() const;

    /** Returns bad day convention used to adjust cash flow dates */
    virtual BadDayConventionConstSP getBadDayConvention() const;

    /** Returns the effective date */
    virtual DateTime spotDate(const DateTime& valueDate) const;

    /** Returns TRUE if paying fee accrued */
    bool isFeeAccrued() const;

    int getSwapFrequency() const;

#ifdef CDS_BACKWARD_COMPATIBILITY
    /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc);

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays);

    /** ParSpreadMaxTweakSize support */
    virtual string sensName(ParSpreadMaxTweakSize* shift) const;
    virtual DoubleArrayConstSP sensMaxTweakSize(ParSpreadMaxTweakSize* shift) const;

#endif

    /** builds a QuantoCDSParSpreads object with the names of the relevant
        market data filled in. Need to call getMarket after this in order
        for the market data to be populated */
    QuantoCDSParSpreads(const string& parSpreadsName,
                        const string& baseCcyName);



    //------------------------------
    // Cache-related methods
    //------------------------------

    /** Returns a DefaultRates object from the cache.
     * This method collects from the cache the default rates associated to the
     * ICDSParSpreads passed as parameter (calculating them if not there in
     * the first place).
     * It is typically invoked from the "defaultRates" method if required (there
     * can be another layer of caching there) and "external classes" are
     * strongly DISCOURAGED from using this method directly: the "defaultRates"
     * method should be used instead */
    virtual const DefaultRatesSP getCachedDefaultRates(
         const ICDSParSpreads* cds,
         const TypeOfEntry     entryType) const;

    /** Provides a facility for previously constructed DefaultRates
     *  objects to be added to the cache */
    virtual const bool cacheDefaultRates(const ICDSParSpreads* cds,
                                         const TypeOfEntry     entryType,
                                         const DefaultRatesSP  entry);
    /** Returns a DefaultRates object.
     * This method does NOT use the (potentially availabe) caching mechanism
     * which avoids the slow process of calculating default rates everytime.
     * It is invoked from inside the cache to calculate the default rates
     * for the first time and its direct use from any external classes is
     * strongly DISCOURAGED: the "defaultRates" method should be used instead */
    virtual const DefaultRatesSP computeDefaultRatesForCache() const;

    /** Resets the pointer to the cache - The cache will not be used in this
     * object from this point onwards */
    virtual void resetCache() const;


    /** Hash code function - the cache needs improved performance compared
     * to the default "hashCode" function in CObjet: only the required
     * components are hashed (see comments in equalToOpt below) */
    virtual int hashCodeOpt() const;

    /** Comparison function - the cache needs improved performance compared
     * to the default "equalTo" function in CObjet: only the required
     * components are compared. This means that the objects being compared
     * are not required to be identical: only equal as far as the default
     * rates calculation is concerned (e.g, if the volType is not used in
     * the calculation, the equalToOpt method should not fail if the volTypes
     * are different) */
    virtual bool equalToOpt(const ICDSParSpreads* cdsBoot) const;

    //------------------------------
    // Duration::IParHandler methods
    //------------------------------
    //return benchmarks for par instruments
    //these will determine the durations to calculate unless specified
    virtual ExpiryArrayConstSP getParBenchmarks() const;

    //return closed form solution
    virtual ExpiryResultArraySP getDuration(const Duration* duration) const;

    //------------------------------------------
    //  IBadDayAdjuster method
    //------------------------------------------
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

    //
    // Implementation of IDiscountCurve
    //

    /** Compute discount factor between two dates. */
    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const ;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. */
    virtual double pv(const DateTime& date) const ;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve) */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const ;

    //
    // Implementation of IDiscountCurveRisky
    //

    /**Returns the probability that there are no default events from d1 to d2, 
       conditional that there are no default events up to d1*/
    virtual double survivalProb(
        const DateTime& d1, 
        const DateTime& d2) 
        const ;
    
    /**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
    virtual double survivalProb(
        const DateTime& dt
        ) const ;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
     * in case of default between startDate and endDate, and zero otherwise */
    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        double              recoveryDelay=0
        ) const ;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
     * in case of default between startDate and endDate, and zero otherwise */
    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        const DateTime&     recoveryDate
	) const ;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a sequence of payments, with simple linear accrued-interest
     * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
     * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for 
     * protectionPV. */
    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
        RecoveryType            accruedRecTyp,
        double                  recoveryTiming=0,
	DateTime                accrueStartDate = DateTime()
        ) const ;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a sequence of payments, with simple linear accrued-interest
     * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
     * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for 
     * protectionPV. */
    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
	RecoveryType            accruedRecTyp,
        const DateTime&         recoveryTiming,
	DateTime                accrueStartDate = DateTime()
        ) const ;

    //------------------------------------------------------
    // Additional Riskless PV methods for flexibility
    //------------------------------------------------------

    /** Compute discount factor between two dates. This is the price
     * at time date1 of a zero-coupon (discount) bond that pays 1
     * at time date2 if it still exists then. 
     * Note that, if this bond is risky, this method
     * in IDiscountCurveRisky means the PV of such a bond which knocks
     * out on default with no recovery at all, and the price is 
     * contingent on no default before date1.
     * @param date1 payment settlement date (conditional on existence at date1)
     * @param date2 payment/zero-coupon bond maturity date
     * @return Discount factor between date1 & date2
     */
    virtual double risklessPV(const DateTime& date1, 
                              const DateTime& date2) const;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. Note settlement is to TODAY, not to
     * the curve spot date. This is because some curves
     * may have ambiguous spot-dates - for example should a combined
     * credit and rates curve have spot date T+1 or T+2?
     * @param date To get discount factor/PV for
     * @return Discount factor between today & given date
     */
    virtual double risklessPV(const DateTime& date) const;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve)
        Cash-flows on or before baseDate are ignored. No
        ordering of the cashflows is assumed. */
    virtual double risklessPV(const CashFlowArray& cashFlows,
                              const DateTime&      baseDate) const;
private:
    QuantoCDSParSpreads();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //utility for IDiscountCurveRisky methods
    YieldCurveConstSP getYieldCurve() const;
    DateTime getValueDate() const;

    ///// fields /////
    string                   name;
    IYieldCurveSP            domYC;              // set to Projection curve
    IYieldCurveSP            forYC;              // set to Projection curve
    ICDSParSpreadsWrapper    parSpreads;         // in original currency
    FXAssetWrapper           fx;                 // from orig ccy to new ccy
    CorrelationSP            corrDomFor;         // dom and foreign IR
    CorrelationSP            corrFxDom;          // fx and dom IR
    CorrelationSP            corrFxFor;          // fx and for IR
    CorrelationSP            corrFxCDS;          // fx and cds par spreads
    CorrelationSP            corrCDSDom;         // cds and domestic iR
    CorrelationSP            corrCDSFor;         // cds and foreign IR
    IAlgorithmConstSP        algorithm;          // how to compute quanto adj $unregistered
    mutable DefaultRatesSP   theDefaultRates;    // buffer for defaultRates() $unregistered
    mutable ParSpreadCurveSP quantodSpreadCurve; // buffer for getParSpreads() $unregistered

    struct MARKET_DLL hashCodeCache {
        int hc;      // Contains the pre-computed hashCode
        bool valid;  // Indicates whether the hashCode is valid or not
    };
    mutable hashCodeCache    theHashCode;        // buffer for hashCodeOpt() $unregistered
};

typedef smartPtr<QuantoCDSParSpreads> QuantoCDSParSpreadsSP;
typedef smartConstPtr<QuantoCDSParSpreads> QuantoCDSParSpreadsConstSP;

DRLIB_END_NAMESPACE

#endif
