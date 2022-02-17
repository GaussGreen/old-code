//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CreditIndex.hpp
//
//   Description : Holds the information required for a credit index
//
//   Author      : Jose Hilera
//
//   Date        : August 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITINDEX_HPP
#define QLIB_CREDITINDEX_HPP

#include "edginc/DateTime.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/CDSIndexParSpreads.hpp" //required for wrapper array somehow
#include "edginc/CreditIndexBase.hpp"
#include "edginc/CreditIndexMap.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/AdjustedCDSParSpreads.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CreditIndexSpreadParallel.hpp"
#include "edginc/CreditIndexSpreadPointwise.hpp"
//#include "edginc/CreditIndexSpreadBypass.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(AdjustedCDSPSwithTweaking);
FORWARD_DECLARE_WRAPPER(CreditIndexBasis);
FORWARD_DECLARE_WRAPPER(CDSIndexParSpreads);
FORWARD_DECLARE_WRAPPER(ICDSParSpreads);
FORWARD_DECLARE(ICreditEventOverride);
typedef array<ICDSParSpreadsWrapper, ICDSParSpreadsWrapper> ICDSParSpreadsWrapperArray;
typedef smartPtr<ICDSParSpreadsWrapperArray>                ICDSParSpreadsWrapperArraySP;


/* The CreditIndex class holds all the information required to turn the spreads 
 * of a portfolio of names into the spreads of the associated index.
 * To do this, it accepts an optional input parameter with the index basis 
 * adjustment. If not provided, this class computes it (it is its main purpose
 * really) - the computation is performed such that:
 *   If spreads_i(t) are the spreads of name i within the portfolio
 *  (indexPortfolio), we define:
 *
 *   Unadjusted duration weighted average:
 *   S(t) =  Sum_i{spreads_i(t) * duration_i(t)} / Sum_i{duration_i(t)}
 *
 *   Basis-adjusted spreads:
 *   adjSpreads_i(t) = spreads_i(t) * [1 + basis(t)/S(t)]
 *
 *   With those basis-adjusted spreads, the corresponding adjusted durations
 *   adjDuration_i(t) are calculated.
 *
 *   Finally, the theoretical index spreads are computed:
 *   indexSpreads(t) = 
 *      Sum_i{adjSpreads_i(t) * adjDuration_i(t)} / Sum_i{adjDuration(t)}
 *
 * We solve for the basis so that indexSpreads(t) match the market provided
 * spreads of the index (indexParSpreads).
 * 
 * To solve for the basis we progress one maturity at a time, using zbrac to 
 * bracket the range where the basis lives, and zbrent to actually find it. 
 * When zbrent (or zbrac) want to see if a certain basis is "the correct one", 
 * they invoke a method in the BasisSolver to compute the difference between 
 * the theoretical indexSpreads and the indexParSpreads.
 * Computing the theoretical indexSpreads requires computing the durations 
 * associated to the adjusted spreads and this implies computing the clean 
 * spreads of the adjusted curves and therefore solving for the default rates 
 * (using zbrent again, inside the zbrac/zbrent solving for the basis).
 * The BasisSolver holds an array of Risky Duration Calculators in order to
 * compute the risky durations as quickly as possible.
 */
class MARKET_DLL CreditIndex : 
    public CreditIndexBase,
    public virtual ICDSParSpreads,
    public virtual ICreditIndexMap,
    public virtual IBadDayAdjuster,
    public virtual ITweakableWithRespectTo<CreditIndexSpreadParallel>,
    public virtual ITweakableWithRespectTo<CreditIndexSpreadPointwise>
{
public:
    static CClassConstSP const TYPE;

    virtual ~CreditIndex();

    /** Returns the name of this object */
    virtual string getName() const;

    /** Get the credit index data from the market data cache, calculating or
     * fetching the index basis as required. */
    virtual void getMarket(const IModel* model, const MarketData *market);

    /** Return an array of AdjustedCDSPSwithTweaking, with one element for each
     * of the underlying names in the index. The adjustment is the index basis
     * calculated/fetched in this classs.
     * CAUTION the index basis is calculated/fetched within the getMarket
     * method, so any attemps to use it for adjustments before that will return,
     * at best, garbage. */
    AdjustedCDSPSwithTweakingArrayConstSP getAllBasisAdjustedCDSPSwithTweaking() const;

    /** Return an array of AdjustedCDSParSpreads, with one element for each
     * of the underlying names in the index. The adjustment is the index basis
     * calculated/fetched in this classs.
     * CAUTION the index basis is calculated/fetched within the getMarket
     * method, so any attemps to use it for adjustments before that will return,
     * at best, garbage. */
    AdjustedCDSParSpreadsArrayConstSP getAllBasisAdjustedCDSParSpreads() const;

    /** Returns the date when to stop tweaking after lastDate */
    const DateTime stopTweaking(const DateTime& lastDate) const;

    /** Validate inputs */
    void validatePop2Object();

    /** Return the index (=names') currency */
    virtual string getCcy() const;

    /** Return the singleName's weight in the index. 
     * Returns 0 if the name is NOT present in the index. */
    virtual double getNamesWeight(const string& singleName) const;

    /** Returns the CDSIndexParSpreads of the index, interpolating
     * the spreads (using the index basis) if required. 
     * Note this method returns a clone of the original index curve 
     * (interpolated if required), so changing the returned curve 
     * will have no impact on the original curve */
    virtual CDSIndexParSpreadsSP getIndexCDSCurve(bool interpolate) const;

    /** Returns the local creditEventOverride, if available.
     * This can be used to keep credit event information in the market data
     * if it is common accross trades, rather than having to enter it on every
     * instrument */
    ICreditEventOverrideSP getCreditEventOverride() const;

    //-------------------------------------------
    // Supporting methods for index sensitivities
    //-------------------------------------------

    // Pointwise tweak (sensitivity) support to the index curve
    virtual string sensName(const CreditIndexSpreadParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexSpreadParallel>& shift);

    // Pointwise tweak (sensitivity) support to the index curve
    virtual string sensName(const CreditIndexSpreadPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CreditIndexSpreadPointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const CreditIndexSpreadPointwise*) const;

    //get the risk mapping matrix that supports
    //CreditIndexSpreadParallel & CreditIndexSpreadPointwise
    // ie we do not tweak the index par spreads directly
    // since that would require re-calculating the basis
    // but rather shift the (per index expiry) mapped points on the
    // single name curves and regress the results from there
    // parallel result is simply the sum of the pointwise
    RiskMappingMatrixConstSP getRiskMappingMatrix() const;

    //supporting method for CreditIndexSpreadPointwise and the
    //generation of its risk mapping matrix
    void bucketShiftSingleNames(const DoubleArrayArray& singleNameTenorShifts);

    //we use a bypass sensitivity where the qualifier
    //holds the bucket shifts for the scenario
    //TweakOutcome CreditIndex::sensShift(const PropertyTweak<CreditIndexSpreadBypass>& tweak);

    //supporting method for CreditIndexSpreadPointwise and the
    //generation of its risk mapping matrix
    //returns one implied spread per input expiry
    DoubleArraySP indexImpliedSpreads(ExpiryArrayConstSP indexExpiries) const;

    //------------------------
    // CreditIndexBase methods
    //------------------------

    /** Return the index basis once calculated/fetched from the market.
     * CAUTION the index basis is calculated/fetched within the getMarket
     * method, so any attemps to get the index basis before that will return,
     * at best, garbage. */
    virtual CreditIndexBasisConstSP getIndexBasis() const;

    //-------------------------------------------------------
    // ICDSParSpreads methods
    // We allow the index to behave as a single name curve
    // which really means we delegate this functionality down
    // the index par spread curve
    //-------------------------------------------------------

    virtual DateTimeArray zeroDates() const;
    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const;
    virtual IDiscountCurve::IKey* getDiscountKey() const;
    virtual DefaultRates::IKey* getRiskyKey() const;

    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate) const;

    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate,
                                    const BadDayConvention* bdc,
                                    const Holiday* hols) const;

    /** Returns CDS par spreads (may be adjusted) */
    virtual DoubleArrayConstSP getParSpreads() const;

    /** Allows a shift to multiple points on the curve simultaneously */
    // tenorShifts is expected to be the same size as the expiries/spreads
    // shifts are additive
    // modifies the object
    // This method is outside the usual sensitivity framework
    virtual void bucketShift(const DoubleArray& tenorShifts);

    /** Returns the expiries (relative: eg. 1M) */
    virtual ExpiryArrayConstSP getParSpreadsExpiries() const;

    /** Returns a DefaultRates object which gives access to 
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object.
        The DefaultRate returned is immutable and therefore it will not be
        changed - which means that there is no need to clone it */
    virtual DefaultRatesSP defaultRates() const;

    /** Returns the market name of the underlying par spreads */
    virtual string getParSpreadsName() const;

    /** Calculates the maximum tweak sizes for this ICDSParSpreads */
    /** Returns an array of the maximum tweak sizes to apply to the CDS par spreads */
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
    virtual bool isFeeAccrued() const;

    /** Swap frequency (= CDSParSpreads::parSwapFreq) */
    virtual int getSwapFrequency() const;

#ifdef CDS_BACKWARD_COMPATIBILITY
     /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc);

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays);
#endif

    /** Returns the recovery rate */
    virtual double getRecovery() const;

    /** Returns true if the name has defaulted */
    virtual bool defaulted() const;
    
    /** Returns the date that this name defaulted. If a default has not
        happened an 'empty' DateTime is returned */
    virtual const DateTime& getDefaultDate() const;

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const;

   /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest. Note that by default the 
        ICDSParSpreads object will not have retrieved any volatilty data from
        the market data cache and in such circumstances this method will fail.
        To avoid this, a model must use an appropriate MarketDataFetcher */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;

    /** Like getName() but returns the name that should be used when retrieving
        correlations against this ICDSParSpreads (eg with Legal Basis you
        don't want a different correlation for each adjusted cds curve) */
    virtual string getCorrelationName() const;

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

    virtual double survivalProb(const DateTime& d1, 
                                const DateTime& d2) const;
 
    virtual double survivalProb(const DateTime& dt) const;

    virtual double getRecovery(const DateTime& defaultDate) const;

    virtual double protectionPV(const DateTime&     paymentDate, 
                                const DateTime&     startDt, 
                                const DateTime&     endDt,
                                RecoveryType        recTyp,
                                double              recoveryDelay=0) const;

    virtual double protectionPV(const DateTime&     paymentDate, 
                                const DateTime&     startDt, 
                                const DateTime&     endDt,
                                RecoveryType        recTyp,
                                const DateTime&     recoveryDate) const;

    virtual double annuityPV(const CashFlowArray&    payments,
                             const DateTime&         paymentDate,
                             RecoveryType            accruedRecTyp,
                             double                  recoveryDelay=0,
	                         DateTime                accrueStartDate = DateTime()) const;

    virtual double annuityPV(const CashFlowArray&    payments,
                             const DateTime&         paymentDate,
                             RecoveryType            accruedRecTyp,
                             const DateTime&         recoveryDate,
	                         DateTime                accrueStartDate = DateTime()) const;

    virtual double risklessPV(const DateTime& date1, 
                              const DateTime& date2) const;
    
    virtual double risklessPV(const DateTime& date) const;
    
    virtual double risklessPV(const CashFlowArray& cashFlows,
                              const DateTime&      baseDate) const;

    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const;
    
    virtual double pv(const DateTime& date) const;
    
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const;

    //------------------------
    // ICreditIndexMap methods
    //------------------------

    /** Return the credit index associated to the single name
     * as defined by this map */
    virtual CreditIndexBaseConstSP getIndex(const string& singleName) const;

    /** Return the single credit index defined by the map */
    virtual CreditIndexBaseConstSP getIndex() const;

    //------------------------------------------
    //  IBadDayAdjuster methods
    //------------------------------------------
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

    /** Calculate the basis so that the index curve calculated using the 
     * basis-adjusted single-name curves matches the market index par spreads */
    // was previously private - now exposed for IndexBasisCache functionality
    CreditIndexBasisSP computeIndexBasis(string name);

    /** Calculates the duration weighted averages of the portfolio, using yc as
        the discount curve */
    static ExpiryResultArraySP calcUnadjustedDurationWeightedAverage(
        ICDSParSpreadsWrapperArraySP indexPortfolio,
        DoubleArraySP                indexWeights,
        YieldCurveConstSP            yc);

private:
    class BasisSolver;
    typedef refCountPtr<BasisSolver> BasisSolverSP;

    // For reflection
    CreditIndex();
    static void load (CClassSP& clazz);
    static IObject* defaultCreditIndex();
    
    // Private methods, only used internally

    /** Utility routine to determine if we are in a fully
     *  defaulted state */
    bool allNamesDefaulted() const;

    /** Computes the index spread at "expiry" by computing the duration
     * weighted average of the names' spreads. Asumes that the index basis
     * has been computed already */
    virtual double interpolateIndexSpread(ExpirySP expiry) const;

    /** Creates a ParSpreadCurveSP by interpolating the index in the names'
     * tenors if required. */
    virtual ParSpreadCurveSP getInterpolatedIndexSpreads() const;

    /** Verify that the names in the portfolio have the right maturities */
    void validatePortfolio() const;

    /** Verify that the index basis and the names in the portfolio are consistent */
    void validateBasis() const;

    /** Factor out the code around solving */
    void doSingleSolve(bool                isPriming,     //drives the way the basis is applied
                       int                 t,             //tenor
                       ExpiryArrayConstSP  indexExpiries,
                       BasisSolver&        basisSolver,
                       double              bestGuess);

    /** provide a best guess to initiate the solving of the basis */
    DoubleArraySP solverBestGuesses(const DateTime&    idxSpotDate,
                                    const DoubleArray& upfronts, //pre-calculated
                                    const IntArray&    matches); //pre-calculated


    // Fields
    string name;                                  // Name of this market object
    ICDSParSpreadsWrapperArraySP indexPortfolio;  // Porfolio of names within the index
    DoubleArraySP                indexWeights;    // Per name weight in the index (sum to 1.0)
    CDSIndexParSpreadsWrapper    indexParSpreads; // The market Par Spreads for the index
    CreditIndexBasisWrapper      indexBasis;      // The basis required to transform the
                                                  // porfolio of names into the index
    CBoolArraySP                 indexSwapsProvideProtection; // Whether this credit index offers
                                                              // protection on names' defaults in swaps 
    CBoolArraySP                 indexIncludesInBasis; // Whether each name should be included in
                                                       // the index basis computation
    ICreditEventOverrideSP creditEventOverride;   // Default parameters overrides

};

typedef smartConstPtr<CreditIndex>   CreditIndexConstSP;
typedef smartPtr<CreditIndex>        CreditIndexSP;
typedef MarketWrapper<CreditIndex>   CreditIndexWrapper;
typedef smartPtr<CreditIndexWrapper> CreditIndexWrapperSP;
typedef array<CreditIndexWrapperSP,CreditIndexWrapper> CreditIndexWrapperArray;

DRLIB_END_NAMESPACE

#endif
