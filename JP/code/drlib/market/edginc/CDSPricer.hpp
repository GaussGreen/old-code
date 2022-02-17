//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CDSPricer.hpp
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CDS_PRICER_HPP
#define QLIB_CDS_PRICER_HPP

#include "edginc/CashFlow.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/NullDecretionCurve.hpp"
#include "edginc/CDSHelper.hpp" //for CParSpreadsDefaultRates

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICDSParSpreads);
FORWARD_DECLARE(IDecretionCurve);
FORWARD_DECLARE(IDiscountCurve);
FORWARD_DECLARE(DayCountConvention);
FORWARD_DECLARE(Expiry);
FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(ICDSBootstrappable);
FORWARD_DECLARE(Holiday);
FORWARD_DECLARE(CreditIndexBasis);
FORWARD_DECLARE_REF_COUNT(DefaultRates);


/**
 * CDSPricer is a "light" pricer for Credit Default Swap.
 * It doesn't live inside CredDefSwap class (product directory)
 * since the pricer can be used to bootstrap par spread curves (market
 * directory).
 * */
class MARKET_DLL CDSPricer : public VirtualDestructorBase {
public:

    enum TimeDiscountType {
        /** No discounting at all for time value of value - users of calculateDefaultPayment()
            need to discount by themself - this is usefull for delay to a specific date */
        NO_TIME_DISCOUNTING,
        /** The delay is given by a period - a default at date t is discounted for t + delay */
        DELAYED_DISCOUNTING,
        /** No delay at all - a default at date t is discounted for t */
        NO_DELAY_DISCOUNTING};


    /*
    Little class to help do the integration. It basically supports computing
    two different integrals along the timeline. It is optimised. Actually
    optimised assuming that calculating log of the default probability is
    quicker than calculating it's exponential. This class could be moved to
    the appropriate DefaultRate class if needed.
    */
    class MARKET_DLL Integrator {
    public:
        //// simple constructor
        Integrator(const DateTime&           today,
                   DateTimeArrayConstSP      dates,
                   IDiscountCurveConstSP     discount,
                   DefaultRatesSP            defaultRates,
                   IDecretionCurveConstSP    prepaySP,
                   DayCountConventionConstSP swpAccrualDCC,
                   TimeDiscountType          timeDiscountType,
                   double                    recoveryDelay);

        //// Returns true if there are no dates to integrate over
        bool empty() const;
        
        //// Returns the first integration date (and not the current interval start
        //// date). Can be called at any time
        const DateTime& firstIntegrationDate() const;

        //// Returns the last integration date (and not the current interval start
        //// date). Can be called at any time
        const DateTime& lastIntegrationDate() const;

        //// integrate over the next interval, returns false if the end has already
        //// been reached
        bool nextStep();

        //// Returns the date of the start of the current integration interval
        //// Only valid once nextStep() has been called
        const DateTime& intervalStartDate() const;

        //// Returns the date of the end of the current integration interval
        //// Only valid once nextStep() has been called
        const DateTime& intervalEndDate() const;

        //// Returns the integral of (Prob of 1st default) across interval
        double integralOfProbDensity() const;

        //// Returns the integral of (PV(t) * Prob of 1st default) across interval
        double integralOfPVTimesProbDensity() const; 

        //// Returns the integral of (Prob of 1st default * t) across interval
        double integralOfProbDensityTimesT() const;

        //// Returns the integral of (PV(t) * Prob of 1st default * t) across interval
        double integralOfPVTimesProbDensityTimesT() const;

    private:
        DateTimeArrayConstSP           dates;          // for integration - includes start
                                                       // and end dates
        auto_ptr<IDiscountCurve::IKey> discFactorKey;
        auto_ptr<DefaultRates::IKey>   defaultPVKey;
        int                            currentIdx;     // where we are in the integration
        double                         logOfDF;        // forward rate of yield curve *
                                                       // -(t2-t1)
        double                         logOfDefaultPV; // forward rate of default
                                                       // curve * -(t2-t1)
        double                         logOfPrepayPV;  // log of pv2/pv1
        IDecretionCurveConstSP         prepaySP;
        double                         riskyPV;        // riskyPV across this interval
        double                         pv0;
        DateTime                       today;          // Defines base
                                                       // date for computing year fracs.
                                                       // It's important to use the same base
                                                       // date consistently throughout the calc */
        TimeDiscountType               timeDiscountType;
        double                         recoveryDelay; 
    #ifndef FIXED_DCC_BUG
        DayCountConventionConstSP      swpAccrualDCC; // bogus complexity

    #endif
    };


    /** Little class to help do the integration required to compute the 
        contingent leg of an FtD. It basically supports computing
        two different integrals along the timeline. It is optimised. Actually
        optimised assuming that calculating log of the default probability is
        quicker than calculating it's exponential. */
    class MARKET_DLL FtDIntegrator {
    public:
        //// simple constructor
        FtDIntegrator(const DateTime&           today,
                      DateTimeArrayConstSP      dates,
                      IDiscountCurveConstSP     discount,
                      DefaultRatesSP            nameDefaultRates,
                      DefaultRatesSP            jointDefaultRates,
                      IDecretionCurveConstSP    prepaySP,
                      TimeDiscountType          timeDiscountType,
                      double                    recoveryDelay);

        /** Returns the integral of (PV(t) * Prob of 1st default) across interval */
        double integralOfFtDPVTimesProbDensity() const; 

        /** Integrate over the next interval, returns false if the end has already
            been reached */
        bool nextStep();

        /** Returns true if there are no dates to integrate over */
        bool empty() const;

    private:
        /** Returns the date of the start of the current integration interval
            Only valid once nextStep() has been called */
        const DateTime& intervalStartDate() const;

        /* Returns the date of the end of the current integration interval
           Only valid once nextStep() has been called */
        const DateTime& intervalEndDate() const;

        DateTimeArrayConstSP           dates;          // for integration - includes start
                                                       // and end dates
        auto_ptr<IDiscountCurve::IKey> discFactorKey;
        auto_ptr<DefaultRates::IKey>   defaultPVKeyName;
        auto_ptr<DefaultRates::IKey>   defaultPVKeyJoint;
        IDecretionCurveConstSP         prepaySP;
        int              currentIdx;     // Where we are in the integration
        double           logOfDF;        // forward rate of yield curve *
                                         // -(t2-t1)
        double           logOfDefaultPVName; // Forward rate of default
                                             // curve * -(t2-t1)
        double           logOfDefaultPVJoint; // Forward rate of default
                                              // curve * -(t2-t1)
        double           logOfPrepayPV;  // Log of pv2/pv1
        double           riskyPV;        // RiskyPV across this interval
        double           pv0;
        DateTime         today;          // Defines base date
        TimeDiscountType timeDiscountType;
        double           recoveryDelay; 
    };


    /** Constructor with optional prepay curve, and fee payments */
    CDSPricer(CashFlowArraySP feePayments,
              DefaultRatesSP defaultRates,
              double notional,
              double recovery,
              bool payAccruedFee,
              DayCountConventionConstSP swpAccrualDCC,
              const DateTime& valueDate,
              const DateTime& protectionStartDate,
              const DateTime& protectionEndDate,
              const DateTime& accruedEffectiveDate,
              IDiscountCurveConstSP discount,
              IDecretionCurveConstSP prepay = IDecretionCurveConstSP()
              );
              
    
    /** Constructor with risky AND riskless fee payments */
    CDSPricer(CashFlowArraySP feePayments,
              CashFlowArraySP risklessFeePayments,
              DefaultRatesSP defaultRates,
              double notional,
              double recovery,
              bool payAccruedFee,
              DayCountConventionConstSP swpAccrualDCC,
              const DateTime& valueDate,
              const DateTime& protectionStartDate,
              const DateTime& protectionEndDate,
              const DateTime& accruedEffectiveDate,
              IDiscountCurveConstSP discount,
              IDecretionCurveConstSP prepay = IDecretionCurveConstSP()
              );

    /** Constructor with fee payments but the defaultRates are
     * produced locally (since MutableDefaultRates is private) */
    CDSPricer(CashFlowArraySP feePayments,
              const DateTimeArray& defaultDates,
              double notional,
              double recovery,
              bool payAccruedFee,
              DayCountConventionConstSP swpAccrualDCC,
              const DateTime& valueDate,
              const DateTime& protectionStartDate,
              const DateTime& protectionEndDate,
              const DateTime& accruedEffectiveDate,
              IDiscountCurveConstSP discount,
              IDecretionCurveConstSP prepay = IDecretionCurveConstSP()
              );
       
    /** Constructor with coupon rate and payment dates */
    CDSPricer(double couponRate,
              const DateTimeArray& paymentDates,
              DefaultRatesSP defaultRates,
              double notional,
              double recovery,
              bool payAccruedFee,
              DayCountConventionConstSP swpAccrualDCC,
              const DateTime& valueDate,
              const DateTime& protectionStartDate,
              const DateTime& protectionEndDate,
              const DateTime& accruedStartDate,
              IDiscountCurveConstSP discount,
              IDecretionCurveConstSP prepay = IDecretionCurveConstSP()
              );

    /** Destructor */
    ~CDSPricer();

    /** Computes the price of that CDS */        
    virtual double price() const;
    
    /** Compute implied par spread, even if CDS is not at par */
    double impliedParSpread(double price) const;

    /** Computes implied par spread and risky duration */
    void impliedParSpreadAndDuration(double& impliedParSpread,       // (O)
                                     double& riskyDuration) const;   // (O)
    
    /** Calculates effect of default. ie returns the sum of minus the
        expected default payment and expected accrual. Unless
        accrualFeeOnly is true in which case just the expected accrual
        is returned. Expected accrual is zero if CDS does not pay
        accrued fees */
    double calculateDefaultPayment(bool accrualFeeOnly) const;

    /** The following is needed for adding simple delay functionality to CDSPricer
        Deliberately this has not been added to CDSPricer's constructor and just to 
        the Integrator sub class and calculateDefaultPayment() methods. This way 
        consumers of CDSPricer know that the effects of delay is not included in all 
        methods. Delay functionality it only available via IDiscountCurveRisky interface and
        CDSParSpreads.
    
        DELAYED_DISCOUNTING simply delays the discounting by a number of calendar days.
        Delay doesn't consider business day adjustments. Although the corresponding input to 
        calculateDefaultPayment() is a double currently only int delay are possible 
        (implicit cast from double to int).

        NO_TIME_DISCOUNTING doesn't discount at all for time value of money - or it just discounts 
        with factor of 1.0.

        NO_DELAY_DISCOUNTING is the usual discounting of time value of money - equivalent to the previous 
        implementation.

        The constructors/method of CDSPricer, Integrator, calculateDefaultPayment() have been modified in a 
        way that it is transparent for existing code and no modifications outside are necessary.

        Search for _DELAY_ in comments CDSPricer.cpp to see all modification

    */
        
    /** Calculates effect of default. ie returns separately the the expected
        default payment and the expected accrual. If accrualFeeOnly is true then
        value of expected default payment will be zero. If the CDS does not pay
        accrued fees then the value of the expected accrual will be zero. */
    void calculateDefaultPayment(bool    accrualFeeOnly,
                                 double& protectionValue,     // (O)
                                 double& accrualValue,        // (O)  
                                 TimeDiscountType timeDiscountType = NO_DELAY_DISCOUNTING,
                                 double recoveryDelay = .0) const; 
    
    /** Calculates the fee payment */
    double calculateFeePayment() const;
    
    /** Bootstraps clean spreads from CDS par spreads */      
    static DefaultRatesSP bootstrap(ICDSBootstrappableConstSP cdsBootstrappable);

    /** Bootstraps clean spreads from CDS par spreads + index basis and
     * returns the "lastProtectionValueLeg / coupon", ie, the risky duration
     * of the cds (one double).
     * Note it is not a static function */
    double bootstrap(ICDSBootstrappableConstSP cdsBoots,
                     ExpiryConstSP             expiryFrom,
                     ExpiryConstSP             expiryTo,
                     CreditIndexBasisConstSP   indexBasis); 

    /** Bootstraps clean spreads from CDS par spreads */      
    static void bootstrapAllowNeg(ICDSBootstrappableConstSP cdsBootstrappable,
                                  const DateTime            maturity,      //(I) allow truncation of the bootstrapping
                                  const DateTime            prevMaturity,  //(I) where we've already bootstrapped
                                  DefaultRatesSP&           defaultRates); //(M) the bootstrapping so far

    /** Computes the risky durations for a given set of expiries */
    static ExpiryResultArraySP riskyDurations(
        const ICDSParSpreadsConstSP cdsParSpreads,
        const IDiscountCurveConstSP discount,
        const ExpiryArrayConstSP expiries);
        
    /** Computes risky durations and implied spreads for the given dates */
    static void impliedParSpreadsAndDurations(
        const ICDSParSpreadsConstSP cdsParSpreads,
        const IDiscountCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads, /* (Output) */
        DoubleArray& durations); /* (Output) */

    /** Computes (forward) risky durations and implied spreads for the given dates */
    static void impliedParSpreadsAndDurations(
        const DateTimeConstSP accrualStartDate,
        const DateTimeConstSP protectionStartDate,
        const ICDSParSpreadsConstSP cdsParSpreads,
        const IDiscountCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads, /* (Output) */
        DoubleArray& durations); /* (Output) */
        
    /** Validate the default rates between the two expiries
     * Will throw an exception if there are issues */
    void validateDefRates() const;

    // Methods to set several attributes in the cds pricer on the fly
    void setDefaultRates(DefaultRatesSP newDefaultRates);
    void setFeePayments(CashFlowArraySP newFeePayments);
    void setProtectionEndDate(const DateTime& newProtectionEndDate);

    /** Set state if in bootstrapping */
    void setBootstrapFlag(bool isBootstrapping);

private:
    static const double DEF_RATE_MIN;
    static const double DEF_RATE_MAX;

    /** Computes the integration dates */
    DateTimeArrayConstSP getIntegrationDates(TimeDiscountType timeDiscountType = NO_DELAY_DISCOUNTING,
                                             double recoveryDelay = 0) const;
    
    /** Creates fee payments from dates with given coupon rate */
    CashFlowArraySP createFeePayments(const DateTimeArray& paymentDates, 
                                      double couponRate) const;
      
    /** Used for bootstrapping only : set the index of the spread
     * currently computed */
    void setCurrentRateIndex(int index);
    
    /** Used for bootstrapping only : set the value of the spread
     * currently computed */
    void setCurrentRate(double rate); 
    
    /** Reset locally cached data */
    void clearCache(); 
    
    void validate();

    /** 
     * Used for bootstrapping only :
     * CallBack used by the zBrentUseful routine to find the default rates 
     * from par spreads
     * */
    static double priceCDSwithGivenRate(double rate, void* inCds);
    
    /** Fee payments */
    CashFlowArraySP feePayments;
    CashFlowArraySP risklessFeePayments;

    /** Notional */                
    double notional;
    
    /** Recovery */
    double recovery;
    
    /** Do we pay accrued fee */
    bool payAccruedFee;
    
    /** Day count convention for accrual periods */
    DayCountConventionConstSP swpAccrualDCC;
    
    /** Value date */
    DateTime valueDate;
    
    /** Protection start date */
    // = "cdsParSpreads->spotDate(valueDate)" in CDS pricing context
    // = "cdsBootstrappable->getEffDate()" in CDS bootstrapping context
    DateTime protectionStartDate;

    /** Protection end date */
    // unadjusted
    DateTime protectionEndDate;
    
    /** */
    //TODO
    // = "swapEffectiveDate" in CDS pricing context
    // = "effDate" in CDS bootstrapping context
    DateTime accruedStartDate;
    //// when the accrual ends. Currently defaults to protectionEndDate.
    //// Should be an explicit input (at least in one constructor)
    DateTime accruedEndDate;
    
    /** Discount curve */
    IDiscountCurveConstSP discount;

    /** Prepay curve */
    IDecretionCurveConstSP prepay;
    
    /** Used for bootstrapping only : index of the spread currently computed */
    int currentRateIndex;

    /** Default rates */
    DefaultRatesSP defaultRates;

    /** Is it in bootstrapping context? */
    bool isBootstrapping;
    
    // -------
    // CACHING
    // -------
    
    /** Cache to avoid computing integration dates many times */
    mutable DateTimeArraySP integrationDatesCache;
    
    // Cache the lastDefaultPayment - when we finish bootstrapping it
    // simplifies calculating the risky duration (since the CDS will
    // be at par, risky duration = lastDefaultPayment / coupon)
    mutable double lastProtectionValueLeg;

    class MutableDefaultRates;
};


typedef smartPtr<CDSPricer> CDSPricerSP;
typedef array<CDSPricerSP, CDSPricer> CDSPricerArray;
typedef smartPtr<CDSPricerArray> CDSPricerArraySP;

DRLIB_END_NAMESPACE

#endif //QLIB_CDS_PRICER_HPP
