//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : AdjustedCDSParSpreads.hpp
//
//   Description : CDSParSpreads with Legal Basis
//   
//   Author      : Antoine Gregoire 
//
//   Date        : Jan 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_ADJUSTED_CDSPARSPREADS_HPP
#define QLIB_ADJUSTED_CDSPARSPREADS_HPP

#include "edginc/CDSParSpreadsAdjustment.hpp"
#include "edginc/Duration.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * Implementation of ICDSBootstrappable used
 * for Legal Basis adjustment.
 * AdjustedCDSParSpreads has a reference to both the ICDSBootstrappable
 * to adjust and the CDSParSpreadsLegalBasis.
 * */
class MARKET_DLL AdjustedCDSParSpreads :
    public BootstrappableCDSParSpreads,
    public virtual Duration::IParHandlerWithClosedForm,
    public virtual ParSpreadMaxTweakSize::IShift {
    
    friend class Dummy; // suppress compiler warning
public:
    static CClassConstSP const TYPE;

    virtual ~AdjustedCDSParSpreads();
    
    // ---------------------
    // MarketObject methods:
    // ---------------------   
    
    virtual string getName() const;

    virtual string getParSpreadsName() const;

    /** Get cdsToAdjust and legalBasisAdjustment from the market data cache */
    virtual void getMarket(const IModel* model, const MarketData *market);
    
    /** Clears cached clean spreads */
    virtual void fieldsUpdated(const CFieldArray& fields);

    IObject* clone() const;

    // -----------------------
    // ICDSParSpreads methods:
    // -----------------------
    
    /** Returns a DefaultRates object which gives access to 
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object */
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
    
    /** @return Yield curve's currency [isocode] */
    virtual string getCcy() const;

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const;
    
    /** Like getName() but returns the name that should be used when retrieving
        correlations against this ICDSParSpreads (eg with Legal Basis you
        don't want a different correlation for each adjusted cds curve) */
    virtual string getCorrelationName() const;

    /** Returns the date when to stop tweaking after lastDate */
    virtual const DateTime stopTweaking(const DateTime& lastDate) const;
    
    /** Returns the day count convention used for par CDS */
    virtual DayCountConvention* dayCountConv() const;
    
    /** Returns the effective date */
    virtual DateTime spotDate(const DateTime& valueDate) const;

    int getSwapFrequency() const;
    
   /** Returns an processed vol - which combines the vol market data with the
       instrument data in the volRequest. A model must have used an 
       appropriate MarketDataFetcher in order for the volatility data to be
       present */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

#ifdef CDS_BACKWARD_COMPATIBILITY
    /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc);

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays);
#endif

    /** ICreditCurve implementation */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate,
                                    const BadDayConvention* bdc,
                                    const Holiday* hols) const;

    /** ICreditCurve implementation */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate) const;
    
    /** ParSpreadMaxTweakSize support */
    virtual string sensName(ParSpreadMaxTweakSize* shift) const;
    virtual DoubleArrayConstSP sensMaxTweakSize(ParSpreadMaxTweakSize* shift) const;
    
    // ---------------------------
    // ICDSBootstrappable methods:
    // ---------------------------
    
    /** Returns holidays used to adjust cash flow dates */
    virtual HolidayConstSP getHolidays() const;
    
    /** Returns bad day convention used to adjust cash flow dates */
    virtual BadDayConventionConstSP getBadDayConvention() const;
    
    /** Returns the corresponding discount curve */
    virtual YieldCurveConstSP getYieldCurve() const;
    
    /** Returns the value date */
    virtual const DateTime& getValueDate() const;
    
    /** Returns the 'effective date' = date used as reference to
     * compute cash flow dates */
    virtual const DateTime getEffDate() const;
    
    /** Returns the expiries dates (absolute: eg. 27/02/2007) */
    virtual const DateTimeArray getExpiryDates() const;
    
    /** Returns the expiries (relative: eg. 1M) */
    virtual ExpiryArrayConstSP getParSpreadsExpiries() const;
    
    /** Returns cash flow dates corresponding to expiry = index
     * (before bad days and holidays adjustment) */
    virtual const DateTimeArray getCashFlowDates(int index) const;
    
    /** Returns the protection end date corresponding to expiry = index */
    virtual const DateTime getProtectionEndDate(int index) const;

    /** Returns TRUE if paying fee accrued */
    virtual bool isFeeAccrued() const;
    
    /** Returns (legal basis adjusted) CDS par spreads */
    virtual DoubleArrayConstSP getParSpreads() const;
 
    /** Allows a shift to multiple points on the curve simultaneously */
    // tenorShifts is expected to be the same size as the expiries/spreads
    // shifts are additive
    // modifies the object
    // This method is outside the usual sensitivity framework
    virtual void bucketShift(const DoubleArray& tenorShifts);

    //------------------------------
    // Cache-related methods
    //------------------------------

    //convenience method for resetting the clean spreads component
    //otherwise only available via fieldsUpdated
    void clearLocalCache() const;

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

    //-----------------------------
    //AdjustedCDSParSpreads methods
    //-----------------------------

    /** Allow construction from "unwrapped" components */
    // NB the input SP's are simply assigned within, no clone
    // is made
    AdjustedCDSParSpreads(BootstrappableCDSParSpreadsSP cds,
                          CDSParSpreadsAdjustmentSP adjustment);

    //------------------------------------------
    //  IBadDayAdjuster method
    //------------------------------------------
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

    //----------------------------
    // IDiscountCurveRisky methods
    //----------------------------

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

    /**Returns the market recovery rate on defaulted assets given default at 
       time defaultDate. This allows the possibility of a term-structure of recovery rates. */
    virtual double getRecovery(
        const DateTime& defaultDate
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
        double                  recoveryDelay=0
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
        const DateTime&         recoveryDate
        ) const ;

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
    /** Constructor */
    AdjustedCDSParSpreads();
    
    /** Invoked when Class is 'loaded', used for reflection */
    static void load(CClassSP& clazz);

    /** Used for reflection */
    static IObject* defaultConstructor();
    
    /** Copy constructor : don't use */  
    AdjustedCDSParSpreads(
        const AdjustedCDSParSpreads& adjustedCDSParSpreads);

    /** Override '=' operator : don't use */
    AdjustedCDSParSpreads& operator=(
        const AdjustedCDSParSpreads& adjustedCDSParSpreads);
    
// FIELDS
    /** Name of the adjusted CDS par spreads */
    string name;
    
    /** Reference to the CDS par spreads to adjust */
    ICDSBootstrappableWrapper cdsToAdjust;
    
    /** arbitrary adjustment */
    //Note that the field name is in the public domain
    //and relates to when this class explitly supported legal basis adjustments.
    //By opening up the type of adjustments, the field name is not meaningful
    //but there is no appetite to change it at the moment.
    //Since (at the time of writing) only index basis adjustments are the only
    //other type of adjustment, and these adjusted types are created on-the-fly
    //and not by the clients
    CDSParSpreadsAdjustmentWrapper legalBasisAdjustment;   

    mutable DefaultRatesSP cleanSpreads;  // "local" cache - not registered, "copy on write" $unregistered
};

typedef smartConstPtr<AdjustedCDSParSpreads> AdjustedCDSParSpreadsConstSP;
typedef smartPtr<AdjustedCDSParSpreads> AdjustedCDSParSpreadsSP;
typedef array<AdjustedCDSParSpreadsSP, AdjustedCDSParSpreads> AdjustedCDSParSpreadsArray;
typedef smartPtr<AdjustedCDSParSpreadsArray> AdjustedCDSParSpreadsArraySP;
typedef smartConstPtr<AdjustedCDSParSpreadsArray> AdjustedCDSParSpreadsArrayConstSP;

// support for wrapper class
typedef MarketWrapper<AdjustedCDSParSpreads> AdjustedCDSParSpreadsWrapper;


DRLIB_END_NAMESPACE

#endif
