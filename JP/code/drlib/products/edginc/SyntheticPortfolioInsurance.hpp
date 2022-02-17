//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SyntheticPortfolioInsurance.hpp
//
//   Description : Synthetic Portfolio Insurance (aka SPI)
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_HPP
#define EDR_SPI_HPP

#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/DeltaTPlusN.hpp"
#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/Events.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/SPIRunTime.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/SPIBondFloor.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/CombinableResult.hpp"

DRLIB_BEGIN_NAMESPACE

// Don't worry about this yet
#define FWD_START_VOL_INTERP 0

#define SPI_ASSET_TYPE_SPI      "SPI Asset"
#define SPI_ASSET_TYPE_ALPHA    "Alpha SPI Asset"
#define SPI_ASSET_TYPE_NON_SPI  "Non SPI Asset"

class SyntheticPortfolioInsurance: public GenericNFBase, 
                                   virtual public IMCIntoProduct,
                                   virtual public DeltaTPlusN::IDeltaTPlusNImnt,
                                   virtual public Theta::Shift, 
                                   virtual public TargetRedemption::IEventHandler,
                                   virtual public SPIFixing::IEventHandler {
protected:
    /// fields ////////
    SPIDynamicBasketSP      basket;
    DateTimeArray           averageOutDates; // dates on which dynamic basket is averaged
    ILockInSPIInterfaceSP   lockIn;
    double                  dayCountBasis;     // day count basis
    ICutoffSPIInterfaceSP   cutoff;
    bool                    isCall;
    double                  strike;
    // optional extra cash flow at maturity, simply added to call/put without 
    // validation as to sign/size. In units of notional (i.e. it's %)
    double                  matCashFlow;
    HolidayWrapper          ccyHols;
    bool                    excludeAssetHols;

    StringArray             gapRiskBucketOffsets; // would be best as ExpiryArray but IMS would barf
    bool                    isRainbowSPI;
    
    bool                    debugOn;
//    bool                    weightedBasket;
//    DoubleArray             tempWeights;
public:
    static CClassConstSP const TYPE;
    friend class SyntheticPortfolioInsuranceMC;
    friend class SyntheticPortfolioInsuranceMCSV;
    friend class RainbowSPIMC;
    friend class RainbowSPIMCSV;
    friend class SPIRunTime;
    friend class SPIRunTimeSV;

    /** Validate instrument having aquired market data */
    void Validate();

    /** Validate dates in instrument having aquired market data 
        in case dates need to be built from date builders whicb
        need market data*/
    void validateDates();

    // validate the basic instrument
    void validatePop2Object();
   
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // for DeltaTPlusN we need to move forward lag many rebalance periods
    virtual int deltaOffset(const Holiday* hols) const;

    //// roll through time (setting historic values)
    virtual bool sensShift(Theta* theta);

    // used by the EDR_GET_BOND_FLOOR_HISTORY addin function to retrospectively
    // fill in the bond floor history given a set of correct yield curves
    // for historic dates
    CashFlowArraySP getBondFloorHistory(MonteCarlo* mc,
                                        DateTimeArray baseDates,
                                        YieldCurveArray ycArray);

    // we validate the SPI before we call getBondFloorHistory it will fail if it detects zeros in historic BF
    // we should make sure any dates we're calculating have something non-zero
    void setRequiredBFDatesToValid(DateTimeArray baseDates);

    // implementation of TargetRedemption::IEventHandler interface
    void getEvents(const TargetRedemption* tarn, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

    // implementation of SPIFixing::IEventHandler interface
    void getEvents(const SPIFixing* fixing, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const;

protected:
        SyntheticPortfolioInsurance(CClassConstSP clazz); // for reflection

private:
    SyntheticPortfolioInsurance(); // for reflection
    SyntheticPortfolioInsurance(const SyntheticPortfolioInsurance& rhs); // not implemented
    SyntheticPortfolioInsurance& operator=(const SyntheticPortfolioInsurance& rhs); // not implemented

    static IObject* defaultSyntheticPortfolioInsurance();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartPtr<SyntheticPortfolioInsurance> SyntheticPortfolioInsuranceSP;

/***************************************************************************/
/** class for holding a couple of results around a gap event occuring
    - gapRiskProbability is the probability that the final basket value
      is below the strike (for a call) or the lock in
    - expectedConditionalValue is the average value of the final basket level
      conditional one of these events occuring */
class SPIGapEventStatistics : public CObject,
                    virtual public CombinableResult {
public:
     static CClassConstSP const TYPE;
     int    numGapEvents;
     double gapEventProbability;
     double expectedConditionalValue;

     SPIGapEventStatistics(); // for reflection

     void evaluateFinalLevel(bool isCall, double finalBasket, 
                             double finalLockIn, double strike);

     void finalise(int numIters, double discFactor);

     void reset();

     //The CombinableResult interface
    /** scale by factor x */
    void scale(double x);

    /** add an object to this result */
    void add(const CombinableResult& x, double scaleFactor);

private:

     SPIGapEventStatistics(const SPIGapEventStatistics& rhs); // not implemented
     SPIGapEventStatistics& operator=(const SPIGapEventStatistics& rhs); // not implemented

     static IObject* defaultStatistics();
        
     /** Invoked when Class is 'loaded' */
     static void load(CClassSP& clazz);
};
typedef smartPtr<SPIGapEventStatistics> SPIGapEventStatisticsSP;

/*****************************************************************************/
// Provides gap risk too
// and expect cutoff info
class PricesSPI: public MCPricesSimple {
public:
    PricesSPI(int                    NbIter,
              int                    NbSubSamples,
              const DateTimeArray&   gapRiskBucketDates);

    /** Clears out SumSubPrices and resets iSubSample */
    virtual void reset();

    virtual int storagePerPath(IMCProduct* product) const;

    virtual void configureCache(const IntArray& changedAssets);

    void addGapRisk(const DoubleArray& gapRiskByBucket);

    void addCutoffData(bool isCutoff,
                       const DateTime& when);

    void addGapEventStats(bool isCall, double finalBasket,
                          double lockIn, double strike);
    
    double getGapRisk() const;

    CashFlowListSP getGapRiskBucketted() const;
    
    double getExpectedCutoff() const;

    SPIGapEventStatisticsSP getGapEventStats(double discFactor) const;

    DateTimeSP getExpectedCutoffDate() const;

    virtual ~PricesSPI();

protected:
    IMCPrices* emptyConstructor() const;

    // fields
    int                     sumCutoffDate;
    int                     numCutoff;
    DoubleArray             gapRiskBuckets;
    DateTimeArray           gapRiskBucketDates;
    SPIGapEventStatisticsSP gapEventStats;
};

class SPIReport : public CObject {
public:
     static CClassConstSP const TYPE;
     int                                   historySize; // how many points (helps limit mem alloc)
     DateTimeArray                         historySampleDate;
     DoubleArray                           historyDynBasket;
     DoubleArray                           historyBond;
     DoubleArray                           historyBF;
     DoubleArray                           historyUE;
     DoubleArray                           historyUC;
     DoubleArray                           historySE;
     DoubleArray                           historyTE;
     DoubleArray                           historySC;
     DoubleArray                           historyBL;
     BoolArray                             historyIsRebal;
     BoolArray                             historyIsCutoff;
     DoubleArray                           historynZ;
     DoubleArray                           historynE0;
     DoubleArray                           historynE1;
     DoubleArray                           historyE0;
     DoubleArray                           historyE1;
     DoubleArray                           historyA0;
     DoubleArray                           historyA1;
     DoubleArray                           historyLB;
     DoubleArray                           historyAccumulatedFee;
     DoubleArray                           historyFeeAmount;
     DoubleArray                           historyCouponAmt;
     DoubleArray                           historyCumDynBasket;
     DoubleArray                           historyCumSE;
     DoubleArray                           historyCumTE;
     DoubleArray                           historyPayoff;
     DoubleArray                           historyOutPerfAsset0;
     DoubleArray                           historyOutPerfAsset1;
     DoubleArray                           historyRC0;
     DoubleArray                           historyRC1;

     SPIReport(int                numDates);

     void downsize();

     void dumpReportStep(int iStep, bool redeemedEarly, SPIRunTime* rtSPI);

private:
     SPIReport(); // for reflection

     SPIReport(const SPIReport& rhs); // not implemented
     SPIReport& operator=(const SPIReport& rhs); // not implemented

     static IObject* defaultReport();
        
     /** Invoked when Class is 'loaded' */
     static void load(CClassSP& clazz);
};
typedef smartPtr<SPIReport> SPIReportSP;
typedef array<SPIReportSP, SPIReport> SPIReportArray;
typedef smartPtr<SPIReportArray> SPIReportArraySP;

/* MC product class for 1 or 2 equity case */
class SyntheticPortfolioInsuranceMC : public IMCProduct,
                                      virtual public IMCProductLN,
                                      virtual public IHandlePaymentEvents {
public:
    double                                unsettledCash;
private:
    // for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS
    DateTimeArray                         instPaymentDates;

    // extra outputs. At the moment some for debug too
    bool                                  today_haveValues; // not defined if forward starting
    DateTime                              today_date; // key for these values
    double                                today_dynBasket;
    double                                today_Z;
    double                                today_nZ;
    DoubleArray                           today_E;
    DoubleArray                           today_nE;
    double                                today_UE;
    double                                today_SE;
    double                                today_TE;
    double                                today_BF;
    double                                today_BL;

    static const string                   DEBUG_ITE; // reports initial target exposure for information
    bool                                  haveDebugITE;
    double                                debugITE;


protected:
    // the following are protected so RainbowSPIMC can see them
    const SyntheticPortfolioInsurance*    inst;         // reference to original instrument
    SPIRunTimeSP                          dynBaskRT;
    int                                   numSimAssets; // as opposed to number of risky assets in the dynamic allocation
    DateTimeArray                         feeNotifDates;
    ILockInSPI*                           lockIn;

public:
    SPIReportSP report;

//public:
//    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    SyntheticPortfolioInsuranceMC(const SyntheticPortfolioInsurance*   inst,
                                  const SimSeriesSP&                   simSeries);

    // override hasFuture so that early termination can be captured in the past
    virtual bool hasFuture() const;

    // Satisfy IHandlePaymentEvents interface
    // for SPI specific record of events (pay dates, cashflows etc)
    void recordEvents(Control* control,
                      Results* results);

    IMCPrices* createOrigPrices(int  nbIter,
                                        int  nbSubSamples,
                                        int  mode);

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  MCPathGenerator*  pathGen)const;

    /** Called within the simulation loop */
    virtual void payoff(const MCPathGenerator*  pathGen,
                        IMCPrices&              prices);

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const MCPathGenerator* pathGen,
                                    int                     iAsset) const;

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const;

    // collects events which should be sitting there
    // called after past is run
    virtual void retrieveEvents(EventResults* events) const;
};

/* MC product class for 1 or 2 equity case */
class SyntheticPortfolioInsuranceMCSV : public MCProductClient,
                                        virtual public IMCProductLN,
                                        virtual public IHandlePaymentEvents {
public:
    double                                unsettledCash;
private:
    // for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS
    DateTimeArray                         instPaymentDates;

    // extra outputs. At the moment some for debug too
    bool                                  today_haveValues; // not defined if forward starting
    DateTime                              today_date; // key for these values
    double                                today_dynBasket;
    double                                today_Z;
    double                                today_nZ;
    DoubleArray                           today_E;
    DoubleArray                           today_nE;
    double                                today_UE;
    double                                today_SE;
    double                                today_TE;
    double                                today_BF;
    double                                today_BL;

    static const string                   DEBUG_ITE; // reports initial target exposure for information
    bool                                  haveDebugITE;
    double                                debugITE;

protected:
    // the following are protected so RainbowSPIMC can see them
    // State variables and generators
    SVGenSpotSP               spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    SVGenDiscFactorSP         dfGen;        //!< Generator for discount factors
    SVGenSpot::IStateVarSP    spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVDiscFactorSP            dfSV;         //!< Df state variable

    const SyntheticPortfolioInsurance*    inst;         // reference to original instrument
    SPIRunTimeSP                          dynBaskRT;
    int                                   numSimAssets; // as opposed to number of risky assets in the dynamic allocation
    DateTimeArray                         feeNotifDates ;
    ILockInSPI*                           lockIn;

public:
    SPIReportSP report;

//public:
//    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    SyntheticPortfolioInsuranceMCSV(const SyntheticPortfolioInsurance*   inst,
				    const SimSeriesSP&                   simSeries);

    // override hasFuture so that early termination can be captured in the past
    virtual bool hasFuture() const;

    // Satisfy IHandlePaymentEvents interface
    // for SPI specific record of events (pay dates, cashflows etc)
    void recordEvents(Control* control,
                      Results* results);


    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    IMCPrices* createOrigPrices(int  nbIter,
                                        int  nbSubSamples,
                                        int  mode);

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  MCPathGenerator*  pathGen)const;

    /** Called within the simulation loop */
    virtual void payoff(const MCPathGenerator*  pathGen,
                        IMCPrices&              prices);

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const MCPathGenerator* pathGen,
                                    int                     iAsset) const;

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const;

    // collects events which should be sitting there
    // called after past is run
    virtual void retrieveEvents(EventResults* events) const;
};

// for calculating bond floor history
class SPIAddin: public CObject{
    static CClassConstSP const TYPE;

    SyntheticPortfolioInsuranceSP spi;
    YieldCurveArray ycArray;
    DateTimeArray baseDates;
    IModelSP model;
    MarketDataSP market;
    
    static IObjectSP getBondFloorHistory(SPIAddin* params);

    /** for reflection */
    SPIAddin();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultSPIAddin();
    
};

DRLIB_END_NAMESPACE

#endif
