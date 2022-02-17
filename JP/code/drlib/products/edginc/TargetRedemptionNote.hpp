//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TargetRedemptionNote.hpp
//
//   Description : 
//
//   Date        : Mar 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/Format.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/Events.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

/** TargetRedemptionNote product - option on doubly rainbowed perfs - across time and asset */
class PRODUCTS_DLL TargetRedemptionNote: public GenericNFBase, 
                            virtual public LegalTerms::Shift,
                            virtual public TargetRedemption::IEventHandler,
                            virtual public IMCIntoProduct{
protected:
    TargetRedemptionNote(CClassConstSP clazz);

    /// fields ////////

    // before aggregating the assetBasket, apply the assetPerf
    IDoubleArrayModifierMakerSP  assetPerf;

    // after aggregating the assetBasket, apply the performance
    IDoubleArrayModifierMakerSP  performance;   

    IAggregateMakerSP            assetBasket;
    bool                         avgFromStart;
    DateTimeArray                averageOutDates;
    DateTimeArray                couponDates;
    ScheduleSP                   frontCoupons;
    
    bool                         isCliquetStyle;
    bool                         isRefPrevAvgOut;
    
    bool                         isMakeWhole;
    double                       targetLevel;
    double                       ecoTargetLevel;
    
    bool                         hasFloater;
    LiborLegSP                   floater;
    
    bool                         floorWithPreviousCoup;
    bool                         receiveOvershoot;
    
    // if true then swap is cancelled when target breached
    // if false then swap is KI'n when target is breached
    // only applies for swap form:
    bool                         isKO; 

    // if true, investor recieves matPerformance if the target is not
    // breached at maturity
    bool                         optionAtMat; 
    
    // at IMS, expose interface to the subclass xxxxxtodo
    IDoubleArrayModifierMakerSP  matPerformance; 
    
    DoubleArray                  bonusCoupons;

    // pays all the cahsflow at the earlier of redemption/last coupon date if true
    bool                         onePayment; 
        
public:
    static CClassConstSP const TYPE;
    friend class TargetRedemptionNoteMC;
    friend class TargetRedemptionNoteSVMC;
    friend class TargetRedemptionNoteHelper;
    friend class AllWeatherTARNMC;
    
    // validation
    void validatePop2Object();

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market); 

    /** Implementation of MonteCarlo::IntoProduct interface - the
    implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; 
    
    bool sensShift(Theta* shift);    

    /** Satisfy LegalTerms::Shift interface */
    bool sensShift(LegalTerms* shift);

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const;

    // TargetRedemption::IEventHandler interface
    void getEvents(const TargetRedemption* target, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

private:
    TargetRedemptionNote(): GenericNFBase(TYPE), assetPerf(0), performance(0), isCliquetStyle(false), 
                            isRefPrevAvgOut(false), isMakeWhole(false), ecoTargetLevel(0.0), hasFloater(false), 
                            floater(0), floorWithPreviousCoup(false), receiveOvershoot(false), isKO(true), 
                            optionAtMat(false), onePayment(false) {} // for reflection
    TargetRedemptionNote(const TargetRedemptionNote& rhs); // not implemented
    TargetRedemptionNote& operator=(const TargetRedemptionNote& rhs); // not implemented
};

/* MC product class for super TargetRedemptionNote */
class PRODUCTS_DLL TargetRedemptionNoteMC: public IMCProduct,
                              virtual public IMCProductLN,
                              virtual public IHandlePaymentEvents {

protected:      
    const TargetRedemptionNote*   inst;     // reference to original instrument
    int                      nbAssets;      // nicer
    IntArray                 couponMap;     // [nbAvgDates] - convenient way to track coupon dates
    SimpleDoubleArray        coupons;  
    
    // Operational aggregation and performance calcs
    IAggregateSP             assetBasket;
    IDoubleArrayModifierSP   assetPerf;

    // record of TargetRedemption event
    bool hasTargetRedeemed;
    struct RedemptionDetails {
        int couponIdx;
        double bonus;
        double finalCoupon;
        double totalCoupon;
        double target;
        double redemption;
    } cachedRedemptionDetails;

private:
    DoubleArray              sum;       // [nbAssets], saves alloc 
    DoubleArray              refLevel;  // [nbAssets], saves alloc 
    SimpleDoubleArray        assetComps;    // asset performances to be processed and aggregated
    
    // historical values
    DoubleArray              sumSoFar;      // [nbAssets]
    DoubleArray              refLevelSoFar; // [nbAssets]
    SimpleDoubleArray        couponsSoFar;  // [nbCoupons]
    int                      iCouponSoFar;
    
    // Operational aggregation and performance calcs
    IDoubleArrayModifierSP   performance;
    
    IntArray                 nbAvgOutPerCouponDate; // [nbCoupons] 
    DoubleArray              couponsFV;
    DoubleArray              liborParticipation;
    
    // for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS
    DateTimeArray            instPaymentDates;          
    OutputRequestUtil::KnownCashFlows knownCFs;
    int                      iFirstUnKnownFloater;  //first unknown floater cash flow.

    // Returns a double array of same size as the couponDates array.
    // The i'th value represents the maturity settlement value of the unpaid funding payments
    // participated in, given that the target is breached at coupon i.
    // Note: if the target is breached at coupon i, then the last funding payment participated in (paid or unpaid)
    // is the first one on or after the coupon date when the breach occured.
    virtual DoubleArraySP computeLiborParticipation(const DateTime& matSettlementDate);

    virtual double computeTargetRedemptionPayoff(
                double aggregateAtMaturity,
                const IPathGenerator*  pathGen);
                    
    virtual void cacheKnownCFs(
                double aggregateAtMaturity,
                const IPathGenerator*  pathGen);
    
    // Returns the coupon idx when redemption occurs
    virtual int processCoupons(
                double aggregateAtMaturity, 
                DoubleArray &newCoupons, /* M */
                bool &earlyRedeemed,
                struct RedemptionDetails& details,
                const IPathGenerator*  pathGen);

    // Modify coupons according to any special (e.g. all-weather) conditions
    // Returns the coupon idx corresponding to the redemption date (ignoring target)
    virtual int couponsOverride(const IPathGenerator*  pathGen,
                                DoubleArray& newCoupons,
                                bool& earlyRedeemed);

    // This function sets 'newCoupons' to the actual cash amounts of the payoff, based on the performances in 'coupons'.
    // It also sets 'couponIdxAtBreach' to the index of the coupon at the point of knockout.  If there is no knockout, then
    // couponIdxAtBreach is set to the coupon.size() - 1.
    virtual void applyTarget(
                double aggregateAtMaturity,
                int &earlyRedeemCpn, // earlier of maturity and early redemption 
                bool &earlyRedeemed,        /* M */
                struct RedemptionDetails& details,
                DoubleArray& newCoupons,
                const double myTarget);   
   
    // Store dates against which to report barriers (if applicable)
    virtual void barrierDates(const int earlyRedeemCpn);

    // Record barriers (if applicable)
    virtual void recordBarriers(Control* control, Results* results) const;

public:

    /** equivalent to InstIntoMCProduct */
    TargetRedemptionNoteMC(const TargetRedemptionNote*         inst,
                           const SimSeriesSP&             simSeries);
    
    /** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
    virtual void initialiseLN(const  IMCPathGenerator*  pathGen) const
    {} // empty
     
    // control variate is done here
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;
    
    virtual void payoff(const IPathGenerator*  pathGen, IMCPrices& prices);
    
    // Satisfy IHandlePaymentEvents interface
    // for Tarn specific record of events (pay dates, cashflows etc)
    virtual void recordEvents(Control* control, Results* results);
            
    // for the LogNormal path generator
    virtual CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;

    // Override MCProduct::retrieveEvents.
    // Mechanism to support TargetRedemption event reporting.
    virtual void retrieveEvents(EventResults* events) const;
};

typedef smartConstPtr<TargetRedemptionNote>  TargetRedemptionNoteConstSP;
typedef smartPtr<TargetRedemptionNote>       TargetRedemptionNoteSP;


/* State Var version */
class PRODUCTS_DLL TargetRedemptionNoteSVMC: public MCProductClient,
                                public IMCStatelessProductClient,
                                virtual public IMCProductLN,
                                virtual public IHandlePaymentEvents {

protected:      
    const TargetRedemptionNote*   inst;     // reference to original instrument
    int                      nbAssets;      // nicer
    IntArray                 couponMap;     // [nbAvgDates] - convenient way to track coupon dates
    SimpleDoubleArray        coupons;  
    
    // Operational aggregation and performance calcs
    IAggregateSP             assetBasket;
    IDoubleArrayModifierSP   assetPerf;

    // record of TargetRedemption event
    bool hasTargetRedeemed;
    struct RedemptionDetails {
        int couponIdx;
        double bonus;
        double finalCoupon;
        double totalCoupon;
        double target;
        double redemption;
    } cachedRedemptionDetails;

private:
    DoubleArray              sum;       // [nbAssets], saves alloc 
    DoubleArray              refLevel;  // [nbAssets], saves alloc 
    SimpleDoubleArray        assetComps;    // asset performances to be processed and aggregated
    
    // historical values
    DoubleArray              sumSoFar;      // [nbAssets]
    DoubleArray              refLevelSoFar; // [nbAssets]
    SimpleDoubleArray        couponsSoFar;  // [nbCoupons]
    int                      iCouponSoFar;
    int                      lastPastDateIdx; // For stateless payoff

    // Operational aggregation and performance calcs
    IDoubleArrayModifierSP   performance;
    
    IntArray                 nbAvgOutPerCouponDate; // [nbCoupons] 
    CashFlowArraySP          liborCashFlows; // working area shared with SV
    
    // for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS
    DateTimeArray            instPaymentDates;          
    OutputRequestUtil::KnownCashFlows knownCFs;
    int                      iFirstUnKnownFloater;  //first unknown floater cash flow.

    // State var stuff
    SVGenSpotSP                  spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    SVGenDiscFactorSP            matDfGen;     //!< Generator for discount factors
    SVGenDiscFactorSP            couponDfGen;  //!< Generator for discount factors from coupon dates
    LiborLeg::LiborLegSVGenSP liborLegGen;  //!< Generator for Libor flows
    SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVDiscFactorSP matDfSV;      //!< Df state variable
    LiborLeg::LiborLegSVSP    liborLegSV;   //!< Libor flows state variable
    SVDiscFactorSP couponDfSV;   //!< Df state variable

    // Not relevant for stoch rates - this always throws an exception
    // The work is now done in computeTargetRedemptionPayoff.
    virtual DoubleArraySP computeLiborParticipation(const DateTime& matSettlementDate);

    virtual double computeTargetRedemptionPayoff(
                double aggregateAtMaturity);
                    
    virtual void cacheKnownCFs(
                double aggregateAtMaturity);
    
    // Returns the coupon idx when redemption occurs
    virtual int processCoupons(
                double       aggregateAtMaturity, 
                DoubleArray& newCoupons, /* M */
                bool&        earlyRedeemed,
                struct RedemptionDetails& details);

    // Modify coupons according to any special (e.g. all-weather) conditions
    // Returns the coupon idx corresponding to the redemption date (ignoring target)
    virtual int couponsOverride(DoubleArray&           newCoupons,
                                bool&                  earlyRedeemed);

    // This function sets 'newCoupons' to the actual cash amounts of the payoff, based on the performances in 'coupons'.
    // It also sets 'couponIdxAtBreach' to the index of the coupon at the point of knockout.  If there is no knockout, then
    // couponIdxAtBreach is set to the coupon.size() - 1.
    virtual void applyTarget(
                double aggregateAtMaturity,
                int &earlyRedeemCpn, // earlier of maturity and early redemption 
                bool &earlyRedeemed,        /* M */
                struct RedemptionDetails& details,
                DoubleArray& newCoupons,
                const double myTarget);   
   
    // Store dates against which to report barriers (if applicable)
    virtual void barrierDates(const int earlyRedeemCpn);

    // Record barriers (if applicable)
    virtual void recordBarriers(Control* control, Results* results) const;

    // IMCStatelessProductClient
    IHistoricalContextSP createHistoricalContext();

    IHistoricalContextSP getInitialHistoricalContext();

    virtual DateTimeArray getPastDates();

    vector<int> finalize( const DateTimeArray& simDates );

    void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices );

protected:      

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

public:

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** equivalent to InstIntoMCProduct */
    TargetRedemptionNoteSVMC(const TargetRedemptionNote*         inst,
                             const SimSeriesSP&             simSeries);
    
    /** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
    virtual void initialiseLN(const  IMCPathGenerator*  pathGen) const
    {} // empty
     
    // control variate is done here
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;
    
    virtual void payoff(const IPathGenerator*  pathGen, IMCPrices& prices);
    
    // Satisfy IHandlePaymentEvents interface
    // for Tarn specific record of events (pay dates, cashflows etc)
    virtual void recordEvents(Control* control, Results* results);
            
    // for the LogNormal path generator
    virtual CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;

    // Override MCProduct::retrieveEvents.
    // Mechanism to support TargetRedemption event reporting.
    virtual void retrieveEvents(EventResults* events) const;
};

DRLIB_END_NAMESPACE
