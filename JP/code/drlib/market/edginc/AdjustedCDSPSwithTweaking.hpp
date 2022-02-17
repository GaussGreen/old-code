//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : AdjustedCDSPSwithTweaking.hpp
//
//   Description : CDSParSpreads with Legal/Index Basis plus tweaking
//   
//   Author      : Gordon Stephens 
//
//   Date        : Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_ADJUSTED_CDSPS_WITH_TWEAKING_HPP
#define QLIB_ADJUSTED_CDSPS_WITH_TWEAKING_HPP

#include "edginc/CDSParSpreadsAdjustment.hpp"
#include "edginc/Duration.hpp"
#include "edginc/CreditIndexBasisParallel.hpp"
#include "edginc/CreditIndexBasisPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * Sibling class to AdjustedCDSParSpreads that supports
 * par spread tweaking after adjustments have been applied
 * to the underlying curve
 **/
class MARKET_DLL AdjustedCDSPSwithTweaking :
    public BootstrappableCDSParSpreads,
    virtual public ITweakableWithRespectTo<ParSpreadParallel>,
    virtual public ITweakableWithRespectTo<ParSpreadParallelRelative>,
    virtual public ITweakableWithRespectTo<ParSpreadPointwise>,
    virtual public ITweakableWithRespectTo<CreditIndexBasisParallel>,
    virtual public ITweakableWithRespectTo<CreditIndexBasisPointwise>,
    virtual public ParSpreadParallelShift::IShift,
    virtual public ParSpreadPropShift::IShift,
    virtual public Duration::IParHandlerWithClosedForm,
    virtual public ParSpreadMaxTweakSize::IShift {
    
    friend class Dummy; // suppress compiler warning
public:
    static CClassConstSP const TYPE;

    virtual ~AdjustedCDSPSwithTweaking();
    
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
    
    /** Returns CDS par spreads 
      * Embellished in this class to apply tweaks on the fly
      */
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
    AdjustedCDSPSwithTweaking(BootstrappableCDSParSpreadsSP cds,
                          CDSParSpreadsAdjustmentSP adjustment);

    //-----------------
    // Tweaking support
    //-----------------

    // PAR SPREAD TWEAKS

    // Parallel tweak (sensitivity) support
    virtual string sensName(const ParSpreadParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadParallel>&);

    // Parallel tweak relative (sensitivity) support
    virtual string sensName(const ParSpreadParallelRelative*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadParallelRelative>&);

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const ParSpreadPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadPointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const ParSpreadPointwise*) const;

    // Proportional shift support
    virtual string sensName(ParSpreadPropShift* shift) const;
    virtual bool sensShift(ParSpreadPropShift* shift);

    // Parallel shift (scenario) support
    virtual string sensName(ParSpreadParallelShift* shift) const;
    virtual bool sensShift(ParSpreadParallelShift* shift);

    // INDEX BASIS TWEAKS

    // Parallel tweak (sensitivity) support
    virtual string sensName(const CreditIndexBasisParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexBasisParallel>&);

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const CreditIndexBasisPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexBasisPointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const CreditIndexBasisPointwise*) const;

    //------------------------------------------
    //  IBadDayAdjuster method
    //------------------------------------------
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

private:

    /** Constructor */
    AdjustedCDSPSwithTweaking();
    
    /** Invoked when Class is 'loaded', used for reflection */
    static void load(CClassSP& clazz);

    /** Used for reflection */
    static IObject* defaultConstructor();
    
    /** Copy constructor : don't use */  
    AdjustedCDSPSwithTweaking(
        const AdjustedCDSPSwithTweaking& adjustedCDSParSpreads);

    /** Override '=' operator : don't use */
    AdjustedCDSPSwithTweaking& operator=(
        const AdjustedCDSPSwithTweaking& adjustedCDSParSpreads);
    
// FIELDS
    /** Name of the adjusted CDS par spreads */
    string name;
    
    /** Reference to the CDS par spreads to adjust */
    ICDSBootstrappableWrapper cdsToAdjust;
    
    /** arbitrary adjustment */
    CDSParSpreadsAdjustmentWrapper adjustment;   

// TRANSIENT TWEAKING FIELDS
    bool   isTweaking;   // true if we are, and expect the state variables to be set
    bool   isAdditive;   // should tweakSize be added or multiplied onto the underlying curve
    double tweakSize;    // state variable for the size of the tweak currently being calculated
    bool   isPointwise;  // true  = apply tweakSize to all curve points
                         // false = apply to tweakPoint only
    int    tweakPoint;   // the point being tweaked, if a pointwise tweak is occuring

    mutable DefaultRatesSP cleanSpreads;  // "local" cache - not registered, "copy on write" $unregistered
};

typedef smartConstPtr<AdjustedCDSPSwithTweaking> AdjustedCDSPSwithTweakingConstSP;
typedef smartPtr<AdjustedCDSPSwithTweaking> AdjustedCDSPSwithTweakingSP;
typedef array<AdjustedCDSPSwithTweakingSP, AdjustedCDSPSwithTweaking> AdjustedCDSPSwithTweakingArray;
typedef smartPtr<AdjustedCDSPSwithTweakingArray> AdjustedCDSPSwithTweakingArraySP;

// support for wrapper class
typedef MarketWrapper<AdjustedCDSPSwithTweaking> AdjustedCDSPSwithTweakingWrapper;


DRLIB_END_NAMESPACE

#endif
