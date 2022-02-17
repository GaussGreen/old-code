//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ICDSParSpreads.hpp
//
//   Description : Interface to the CDS par spreads
//
//   Date        : 18 August 2005
//
//----------------------------------------------------------------------------

#ifndef ICDSPARSPREADS_HPP
#define ICDSPARSPREADS_HPP

#include "edginc/DefaultRates.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ParSpreadCurve.hpp"
#include "edginc/RecoveryShiftBase.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CreditDefaultSensBase.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/IABCDSDecretion.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

// Some hacking to ensure backward compatibility
// in CDSParSpreads related test cases.
// The following fields in CDSParSpreads were NOT
// mandatory in the previous versions so we have to 
// populate them if necessary at the code level :
// - discount
// - bad day convention
// - holidays
// THIS SHOULD BE REMOVED from ICDSParSpreads and CDSParSpreads 
// once the test cases are updated
#define CDS_BACKWARD_COMPATIBILITY

FORWARD_DECLARE_WRAPPER(ICDSParSpreads);
#ifndef QLIB_ICDSPARSPREADS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ICDSParSpreads>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ICDSParSpreads>);
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSParSpreads>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<ICDSParSpreadsWrapper>(
                    ICDSParSpreadsWrapper* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetInLine<ICDSParSpreadsWrapper>(
                    ICDSParSpreadsWrapper* t, IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ICDSParSpreads>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ICDSParSpreads>);
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSParSpreads>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<ICDSParSpreadsWrapper>(
                         ICDSParSpreadsWrapper* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetInLine<ICDSParSpreadsWrapper>(
                         ICDSParSpreadsWrapper* t, IObjectSP o));
#endif
FORWARD_DECLARE(CreditIndexBasis);
class RiskyDurationCalculator;
typedef smartPtr<RiskyDurationCalculator> RiskyDurationCalculatorSP;

/* Type of entry - used by ICDSParSpreadsCache to determine the best way
 * to cache the clean spreads associated to this ICDSParSpreads.
 * Note that this is kind of type information but we do not need such fine 
 * granularity (actually at the moment ICDSParSpreadsCache only cares whether
 * the curve is quanto or not) */
typedef enum {
    QUANTO_CDS_PAR_SPREADS,
    ADJUSTED_CDS_PAR_SPREADS,
    BASE_CDS_PAR_SPREADS,
    UNTWEAKABLE_CDS_PAR_SPREAD
} TypeOfEntry;


/** Interface for objects representing CDSParSpreads - allows us to implement
    different representations. This is still a little in development */
class MARKET_DLL ICDSParSpreads: public virtual IABCDSDecretion,
                                 public virtual IMarketFactor,
                                 public virtual ICreditCurve,
                                 public virtual IDiscountCurveRisky,
                                 public virtual IBadDayAdjuster
{
public:
    static CClassConstSP const TYPE;
    virtual ~ICDSParSpreads();

	/**
	 * Create risky yield curve.
	 */
	virtual IYieldCurveSP makeRiskyCurve(
		const IYieldCurve& yc, 
		const DateTime*    maturityDate=NULL) const;

    CreditSpreadCurveSP makeCreditSpreadCurve(
        const DateTime&    date, 
        const IYieldCurve& yc) const;

    /** Questions:
        1. What other methods should go in here?
        Answers:
        1. Actual getParSpreads() and getParSpreadsExpiries(). [William]
    */

    /** Returns CDS par spreads (may be adjusted) */
    virtual DoubleArrayConstSP getParSpreads() const = 0;

    /** Returns the RiskyDurationCalculator associated to this ICDSParSpreads */
    virtual RiskyDurationCalculatorSP getRiskyDurationCalculator() const;

    /** Allows a shift to multiple points on the curve simultaneously */
    // tenorShifts is expected to be the same size as the expiries/spreads
    // shifts are additive
    // modifies the object
    // This method is outside the usual sensitivity framework
    virtual void bucketShift(const DoubleArray& tenorShifts) = 0;


    /** Returns the expiries (relative: eg. 1M) */
    virtual ExpiryArrayConstSP getParSpreadsExpiries() const = 0;
    
    /** Returns a DefaultRates object which gives access to 
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object.
        The DefaultRate returned is immutable and therefore it will not be
        changed - which means that there is no need to clone it */
    virtual DefaultRatesSP defaultRates() const = 0;

    //default implementation of IDiscountCurveRiskyMethod
    //just calls defaultRates()
    virtual DefaultRatesSP getDefaultRates() const;

    // to do - add methods from DefaultRates here. Consider retiring above
    // although might be too much of a culture shock.

    //NEW METHODS:
    
    /** Returns market name of the object (same as MarketObject::getName())*/ 
    virtual string getName() const = 0;

    /** Returns the market name of the underlying par spreads */
    virtual string getParSpreadsName() const = 0;

    /** Calculates the maximum tweak sizes for this ICDSParSpreads */
    /** Returns an array of the maximum tweak sizes to apply to the CDS par spreads
        for ParSpreadRhoPointwise - same algorithm as used by CMLib.
        The algorithm is as follows :
        1) Define the vector S( l , r ) to be the benchmark par spreads as a function of clean spreads  l, and interest rates r
        2) Find the Gradient of S with respect to l denoted by GradS -- it is a matrix
        3) dS = GradS * dl or equivalently dl = GradS-1 * dS
        4) Define bounds for dl such that the curve remains valid.
        5) Given a tweak dS if dl goes beyond its prescribed bounds then scale down the tweak size dS to conform to the bound. */
    virtual DoubleArrayConstSP calculateTweakSizes() const = 0;
        
    //TEMPORARY METHODS:
    
    /** Returns the day count convention used for par CDS */
    virtual DayCountConvention* dayCountConv() const = 0;
    
    /** Returns holidays used to adjust cash flow dates */
    virtual HolidayConstSP getHolidays() const = 0;
    
    /** Returns bad day convention used to adjust cash flow dates */
    virtual BadDayConventionConstSP getBadDayConvention() const = 0;
    
    /** Returns the effective date */
    virtual DateTime spotDate(const DateTime& valueDate) const = 0;

    /** Returns TRUE if paying fee accrued */
    virtual bool isFeeAccrued() const = 0;

    /** Swap frequency (= CDSParSpreads::parSwapFreq) */
    virtual int getSwapFrequency() const = 0;

#ifdef CDS_BACKWARD_COMPATIBILITY
     /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc) = 0;

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays) = 0;
#endif

    //END OF MODIFICATIONS.
    
    /** Returns the recovery rate */
    virtual double getRecovery() const = 0;

    /** Returns true if the name has defaulted */
    virtual bool defaulted() const = 0;
    
    /** Returns the date that this name defaulted. If a default has not
        happened an 'empty' DateTime is returned */
    virtual const DateTime& getDefaultDate() const = 0;

    /** return CDS Par Spreads currency [isocode] */
    //// fails if no yield curve present
    virtual string getCcy() const = 0;

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const = 0;

    /** Returns the date when to stop tweaking after lastDate */
    virtual const DateTime stopTweaking(const DateTime& lastDate) const = 0;

   /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest. Note that by default the 
        ICDSParSpreads object will not have retrieved any volatilty data from
        the market data cache and in such circumstances this method will fail.
        To avoid this, a model must use an appropriate MarketDataFetcher */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const = 0;

    /** Like getName() but returns the name that should be used when retrieving
        correlations against this ICDSParSpreads (eg with Legal Basis you
        don't want a different correlation for each adjusted cds curve) */
    virtual string getCorrelationName() const = 0;


    //////////////////// cache-related methods //////////////////////

    // For reference see the "Caching of default rates (clean spreads)" 
    // design document in the QLib DB 

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
        const TypeOfEntry     entryType) const = 0;

    /** Provides a facility for previously constructed DefaultRates
     *  objects to be added to the cache */
    virtual const bool cacheDefaultRates(const ICDSParSpreads* cds,
                                         const TypeOfEntry     entryType,
                                         const DefaultRatesSP  entry) = 0;
    /** Returns a DefaultRates object.
     * This method does NOT use the (potentially availabe) caching mechanism
     * which avoids the slow process of calculating default rates everytime.
     * It is invoked from inside the cache to calculate the default rates
     * for the first time and its direct use from any external classes is
     * strongly DISCOURAGED: the "defaultRates" method should be used instead */
    virtual const DefaultRatesSP computeDefaultRatesForCache() const = 0;

    /** Resets the pointer to the cache - The cache will not be used in this
     * object from this point onwards */
    virtual void resetCache() const = 0;


    /** Hash code function - the cache needs improved performance compared
     * to the default "hashCode" function in CObjet: only the required
     * components are hashed (see comments in equalToOpt below) */
    virtual int hashCodeOpt() const = 0;

    /** Comparison function - the cache needs improved performance compared
     * to the default "equalTo" function in CObjet: only the required
     * components are compared. This means that the objects being compared
     * are not required to be identical: only equal as far as the default
     * rates calculation is concerned (e.g, if the volType is not used in
     * the calculation, the equalToOpt method should not fail if the volTypes
     * are different) */
    virtual bool equalToOpt(const ICDSParSpreads* cdsBoot) const = 0;


// methods coming from IDiscountCurveRisky
// should be overwritten by derived classes
// in the long term these methods should be become pure virtual,
// but they have been added here temporarilt so that not all derived
// classes need to implement them for now
    virtual double survivalProb(
        const DateTime& /*d1*/, 
        const DateTime& /*d2*/) 
        const  { throw ModelException("ICDSParSpreads::survivalProb(d1,d2)",
                                      "not implemented yet"); }
 
    virtual double survivalProb(
        const DateTime& /*dt*/
        ) const { throw ModelException("ICDSParSpreads::survivalProb(dt)",
                                       "not implemented yet"); }


    virtual double getRecovery(
        const DateTime& /*defaultDate*/
        ) const { throw ModelException("ICDSParSpreads::getRecovery(dt)",
                                       "not implemented yet"); }
    

    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        double              recoveryDelay=0
        ) const { throw ModelException("ICDSParSpreads::protectionPV(5.1)",
                                       "not implemented yet"); }

    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        const DateTime&     recoveryDate
        ) const { throw ModelException("ICDSParSpreads::protectionPV(5.2)",
                                       "not implemented yet"); }

    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
        RecoveryType            accruedRecTyp,
        double                  recoveryDelay=0,
	DateTime                accrueStartDate = DateTime()
        ) const { throw ModelException("ICDSParSpreads::annuityPV(5.1)",
                                       "not implemented yet"); }

    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
        RecoveryType            accruedRecTyp,
        const DateTime&         recoveryDate,
	DateTime                accrueStartDate = DateTime()
        ) const { throw ModelException("ICDSParSpreads::annuityPV(5.2)",
                                       "not implemented yet"); }

    virtual double risklessPV(
        const DateTime& date1, 
        const DateTime& date2
        ) const { throw ModelException("ICDSParSpreads::risklessPV(dt,dt)",
                                       "not implemented yet"); }
    
    virtual double risklessPV(
        const DateTime& date
        ) const { throw ModelException("ICDSParSpreads::risklessPV(dt)",
                                       "not implemented yet"); }
    
    virtual double risklessPV(
        const CashFlowArray& cashFlows,
        const DateTime&      baseDate
        ) const { throw ModelException("ICDSParSpreads::risklessPV(cfs,dt)",
                                       "not implemented yet"); }

    double pv(
        const DateTime& date1, 
        const DateTime& date2
        ) const { throw ModelException("ICDSParSpreads::pv(dt,dt)",
                                       "not implemented yet"); }
    
    double pv(
        const DateTime& date
        ) const { throw ModelException("ICDSParSpreads::pv(dt)",
                                       "not implemented yet"); }
    
    double pv(
        const CashFlowArray& cashFlows,
        const DateTime&      baseDate
        ) const { throw ModelException("ICDSParSpreads::pv(cfs,dt)",
                                       "not implemented yet"); }

    //////////////////// static utility methods //////////////////////

    /** Static method to create a cash flow array using the suplied arguments */
    static CashFlowArraySP getCashFlowArray(const DateTime&           maturity, 
                                            const DateTime&           start,
                                            const BadDayConvention*   bdc,
                                            const DayCountConvention* dcc,
                                            const Holiday*            hols,
                                            const int                 swapFreq,
                                            double                    spread);

    /** If supplied CDSParSpreads wrapper is using the market data cache, then
        retrieves appropriate CDSParSpreads object, and if necessary, applies
        quanto adjustment. The model must have been informed of the 
        'domestic' currency (see Model::getDomesticYCName()) */
    static void getMarketData(
        const IModel*            model, 
        const MarketData*        market,
        ICDSParSpreadsWrapper&   parSpreadsWrapper);

    /** Same as above but takes domestic yield curve name in directly */
    static void getMarketData(
        const IModel*            model, 
        const MarketData*        market,
        const string&            domesticYCName,
        ICDSParSpreadsWrapper&   parSpreadsWrapper);

    /** Retrieve, and make if necessary, a ICDSParSpreads object which has the
        specified currency dictated by the discount curve */
    static MarketObjectSP buildUsingCache(
        const IModel*     model, 
        const MarketData* market,
        const string&     domesticYC,
        const string&     cdsParSpreadName);

    //------------------------------------------
    //  IABCDSDecretion methods
    //------------------------------------------
    virtual IDecretionCurveConstSP getPrepayCurve() const;
    virtual DecretionCurveConstSP getPrepayCurveObj() const;

private:
    class CleanSpreadsBuilder;
    class ImpliedSpreadAndDurationBuilder;
};


/** Specialisation of arrayObjectCast (needed as the array is not an array
 *  of pointers) */
template <> class MARKET_DLL arrayObjectCast<ICDSParSpreadsWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const ICDSParSpreadsWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(ICDSParSpreadsWrapper& value);

    /** Turns the IObjectSP into a ICDSParSpreadsWrapper */
    static const ICDSParSpreadsWrapper& fromIObject(IObjectSP& value);
};

typedef array<ICDSParSpreadsWrapper, ICDSParSpreadsWrapper> ICDSParSpreadsWrapperArray;
typedef smartPtr<ICDSParSpreadsWrapperArray>                ICDSParSpreadsWrapperArraySP;
typedef array<ICDSParSpreadsSP, ICDSParSpreads> ICDSParSpreadsArray;


DRLIB_END_NAMESPACE

#endif
