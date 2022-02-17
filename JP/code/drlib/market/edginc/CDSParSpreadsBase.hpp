//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CDSParSpreadsBase.hpp
//
//   Description : Base class for CDS par spreads
//
//   Date        : August 30, 2001
//
//----------------------------------------------------------------------------

#ifndef CDSPARSPREADSBASE_HPP
#define CDSPARSPREADSBASE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Theta.hpp"
#include "edginc/RecoveryShiftBase.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/CreditSpreadRhoPointwise.hpp"
#include "edginc/BootstrappableCDSParSpreads.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/ICDSVol.hpp"
#include "edginc/Duration.hpp"
#include "edginc/RiskyDurationCalculator.hpp"
#include "edginc/DecretionCurve.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

// Some hacking to ensure backward compatibility
// in CDSParSpreads related test cases.
// The following fields in CDSParSpreadsBase were NOT
// mandatory in the previous versions so we have to
// populate them if necessary at the code level :
// - discount
// - bad day convention
// - holidays
// THIS SHOULD BE REMOVED once the test cases are updated
#define CDS_BACKWARD_COMPATIBILITY


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(CVolRequest);
FORWARD_DECLARE(DayCountConvention);
FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(CDSPricer);
FORWARD_DECLARE(ICDSParSpreadsCache);
FORWARD_DECLARE(PrepayCurve);
FORWARD_DECLARE_REF_COUNT(DefaultRates);
FORWARD_DECLARE_WRAPPER(Holiday);
FORWARD_DECLARE_WRAPPER(YieldCurve);
FORWARD_DECLARE_WRAPPER(ParSpreadCurve);
FORWARD_DECLARE_WRAPPER(IDecretionCurve);
FORWARD_DECLARE_WRAPPER(DecretionCurve);
FORWARD_DECLARE(CreditIndexBasis);




/*******************************************************************************
 *                         CDSParSpreadsBase declarations
 *******************************************************************************/

/** Holds the current par spreads for CDSs */
class MARKET_DLL CDSParSpreadsBase: public BootstrappableCDSParSpreads,
                         public virtual Duration::IParHandlerWithClosedForm,
                         public virtual Duration::IParHandlerWithoutClosedForm,
                         public virtual Theta::IShift,
                         public virtual RecoveryShiftBase::IRestorableShift,
                         public virtual CreditSpreadRhoParallel::RestorableShift,
                         public virtual CreditSpreadRhoPointwise::IRestorableShift,
                         public virtual ParSpreadMaxTweakSize::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class CDSParSpreadsBaseHelper;
    friend class LiquiditySpreadCurve;
    friend class CreditCurveDDE;

    virtual ~CDSParSpreadsBase();
    //// checks parameters immediately after object is constructed
    void validatePop2Object();

    /** Clears cached clean spreads */
    virtual void fieldsUpdated(const CFieldArray& fields);

    virtual IObject* clone() const;

    virtual string getName() const;

    virtual string getParSpreadsName() const;

    /** @return [discount] yield curve's currency [isocode] */
    virtual string getCcy() const;

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const;

    /** Returns a DefaultRates object which gives access to
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object
        The DefaultRate returned is immutable and therefore it will not be
        changed - which means that there is no need to clone it */
    virtual DefaultRatesSP defaultRates() const;

    /** return the bootstrapped dates */
    virtual DateTimeArray zeroDates() const;

    int                getSpotOffset() const;
    DoubleArrayConstSP getUpfrontFees() const;

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

   /** Returns a processed vol - which combines the vol market data with the
       instrument data in the volRequest. A model must have used an
       appropriate MarketDataFetcher in order for the volatility data to be
       present */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;

    /** Returns the same as getName() */
    virtual string getCorrelationName() const;

    /** Returns the day count convention used for par CDS */
    virtual DayCountConvention* dayCountConv() const;

    /** Returns the effective date */
    virtual DateTime spotDate(const DateTime& valueDate) const;

    /** Returns TRUE if paying fee accrued */
    virtual bool isFeeAccrued() const;

    /** Returns CDS par spreads */
    virtual DoubleArrayConstSP getParSpreads() const;

    /** Swap frequency (= CDSParSpreadsBase::parSwapFreq) */
    int getSwapFrequency() const;

    DateTimeArray getExpiryDates(DateTime today,
                                 const BadDayConvention* bdc,
                                 const Holiday* hols) const;

    ExpiryArrayConstSP getExpiries() const;

    int numSpreads();


    /** Allows a shift to multiple points on the curve simultaneously */
    // tenorShifts is expected to be the same size as the expiries/spreads
    // shifts are additive
    // modifies the object
    // This method is outside the usual sensitivity framework
    virtual void bucketShift(const DoubleArray& tenorShifts);

    /** Shifts the object using given shift (see Theta::Shift)*/
    virtual bool sensShift(Theta* shift);

    /** Returns name identifying yield curve for rho parallel */
    virtual string sensName(IRecoveryPerturb* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(IRecoveryPerturb* shift);
    /** Restores the object to its original form */
    virtual void sensRestore(IRecoveryPerturb* shift);

    /** SpreadRhoParallel support */
    virtual string sensName(CreditSpreadRhoParallel* shift) const;
    virtual bool sensShift(CreditSpreadRhoParallel* shift);
    virtual void sensRestore(CreditSpreadRhoParallel* shift);

    /** SpreadRhoPointwise support */
    virtual string sensName(CreditSpreadRhoPointwise* shift) const;
    virtual ExpiryArrayConstSP sensExpiries(CreditSpreadRhoPointwise* shift) const;
    virtual bool sensShift(CreditSpreadRhoPointwise* shift);
    virtual void sensRestore(CreditSpreadRhoPointwise* shift);

    /** ParSpreadMaxTweakSize support */
    virtual string sensName(ParSpreadMaxTweakSize* shift) const;
    virtual DoubleArrayConstSP sensMaxTweakSize(ParSpreadMaxTweakSize* shift) const;

    /** Get the par spreads curve from the market data cache */
    virtual void getMarket(const IModel* model, const MarketData *market);

    /**Returns the market recovery rate on defaulted assets given default at 
    time defaultDate. This allows the possibility of a term-structure of recovery rates. */
    virtual double getRecovery(const DateTime& defaultDate) const;

    double getRecovery() const;

    bool getAccrueFee() const;

    /** to be removed - should use bdc and hols from this object. Until
        then bdc and hols are optional, if null then values from this are
        used */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate,
                                    const BadDayConvention* bdc,
                                    const Holiday* hols) const;

    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate) const;

    /** Returns FALSE, indicating the name has NOT defaulted.
     * Derived classes can override this method to keep track of defaults */
    virtual bool defaulted() const;

    /** Returns the date that this name defaulted. Since the name has
     * not happened, an 'empty' DateTime is returned.
     * Derived classes can override this method to keep track of defaults */
    virtual const DateTime& getDefaultDate() const;

    /** Returns the date when to stop tweaking after lastDate */
    virtual const DateTime stopTweaking(const DateTime& lastDate) const;

    /** Returns holidays used to adjust cash flow dates */
    virtual HolidayConstSP getHolidays() const;

    /** Returns bad day convention used to adjust cash flow dates */
    virtual BadDayConventionConstSP getBadDayConvention() const;

#ifdef CDS_BACKWARD_COMPATIBILITY
    /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc);

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays);
#endif
    
    // ---------------------------
    // ICDSBootstrappable methods:
    // ---------------------------

//    /** Returns day count convention used to compute cash flow amounts */
//    virtual DayCountConventionConstSP getDCC() const;

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

    /** Returns the protection end date corresponding to expiry = index */
    virtual const DateTime getProtectionEndDate(int index) const;

    /** Return prepay curve */
    virtual IDecretionCurveConstSP getPrepayCurve() const;
    virtual DecretionCurveConstSP getPrepayCurveObj() const;  //return object version

    /** Returns cash flows corresponding to expiry = index.
     ** Reuse BootstrappableCDSParSpreads::getCashFlowArray(),
     ** and add upfront fees to the first element corresponding to
     ** start date
     */
    virtual CashFlowArraySP getCashFlowArray(int index) const;
    virtual CashFlowArraySP getCashFlowArray(CreditIndexBasisConstSP indexBasis,
                                             int index) const;

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

    //return a means of tweaking the MarketObject in a pointwise manner
    //assumption is that the results of this is a VectorResult.....
    virtual VectorRiskPropertySensitivityConstSP getPointwiseTweaker() const;

    //return a par instrument for the specified maturity
    virtual InstrumentSP getParInstrument(const ExpirySP maturity) const;

    //return benchmarks for par instruments
    //these will determine the durations to calculate unless specified
    virtual ExpiryArrayConstSP getParBenchmarks() const;

    //return closed form solution
    virtual ExpiryResultArraySP getDuration(const Duration* duration) const;

    //static function to return closed form solution - used by several
    //different subclasses of CDSParSpreadsBase
    static ExpiryResultArraySP CDSGetDuration(const Duration* duration,
                                              ExpiryArrayConstSP altBenchmarks,
                                              string altYCName,
                                              string altName);

    //--------------------------------------
    //Support for creating risky zero curves
    //--------------------------------------
    IYieldCurveSP makeRiskyCurve(
	    const IYieldCurve& yc,
	    const DateTime*    maturityDate) const;


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


protected:
    CDSParSpreadsBase(CClassConstSP clazz = TYPE);

    /** A constructor for convenience to be called by other QLib classs */
    CDSParSpreadsBase(string                   name,
                      DateTime                 valueDate,
                      int                      spotOffset,
                      double                   parRecovery,
                      int                      parSwapFreq,
                      ExpiryArraySP            expiries,
                      CDoubleArraySP           spreads,
                      CDoubleArraySP           upfronts,
                      bool                     parAccrueFee,
                      string                   parDCC,       // day count convention for par CDS
                      string                   parBDC,       // bad day convention for par CDS
                      HolidayConstSP           parHols,      // holidays for par CDS
                      YieldCurveConstSP        discount,     // corresponding discount curve
                      DecretionCurveConstSP    prepay        // prepay curve
                      );

    /** Allow setting the par spreads. Used, e.g., to allow interpolating
     * the original par spreads - note this method is protected */
    void setParSpreadCurveWrapper(ParSpreadCurveWrapper curveWrapper);

//     static void acceptWrapperNameCollector(CDSParSpreadsBase* parSpreads,
//                                            WrapperNameCollector* collector);

    /// non-static fields ////////
    string                 name;
    DateTime               valueDate;
    double                 parRecovery;
    mutable DefaultRatesSP cleanSpreads; // "local" cache - not registered, "copy on write" $unregistered

private:
    //--------------------------------------
    //Support for creating risky zero curves
    //--------------------------------------
    IObjectSP getRiskyComponent() const;

    class RDCalculator;

    /** Checks that recovery is in [0,1] */
    bool isRecoveryValid() const;

    /** Create a creditSpreadCurve from the CDSCurve */
    void buildCreditSpreadCurve() const;

    /// non-static fields ////////
    ParSpreadCurveWrapper  parSpreads;
    int                    parSwapFreq;
    string                 parDCC; // day count convention for par CDS
    string                 parBDC; // bad day convention for par CDS
    HolidayWrapper         parHols;// holidays for par CDS
    bool                   parAccrueFee;
    int                    spotOffset;
    YieldCurveWrapper      discount; // corresponding discount curve
    ICDSVolWrapper         cdsVol;   /* holds the cds vol, optional */
    DecretionCurveWrapper  prepay;   //prepay curve

    // Transient
    DayCountConventionSP   parDayCountConv;
    BadDayConventionSP     parBadDayConv;

    // Transient, to support credit spread sensitivities
    mutable bool                   priceViaCreditSpreadCurve;
    mutable CreditSpreadCurveSP    creditSpreadCurve;
    mutable CreditSpreadCurveSP    csRestore; //to restore after tweaking

    mutable ICDSParSpreadsCacheSP cache;     // <ICDSParSpreads, cleanSpreads> cache $unregistered
};

typedef smartConstPtr<CDSParSpreadsBase> CDSParSpreadsBaseConstSP;
typedef smartPtr<CDSParSpreadsBase>      CDSParSpreadsBaseSP;
typedef MarketWrapper<CDSParSpreadsBase> CDSParSpreadsBaseWrapper;

DRLIB_END_NAMESPACE

#endif
