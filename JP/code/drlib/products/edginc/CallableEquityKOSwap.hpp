//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallableEquityKOSwap.hpp
//
//   Description : Callable Equity KO Swap 
//
//----------------------------------------------------------------------------

#ifndef EDG_CKS_HPP
#define EDG_CKS_HPP

#include "edginc/Generic1Factor.hpp"

#include "edginc/LastSensDate.hpp"
#include "edginc/LegalTerms.hpp"

#include "edginc/LatticeProdEDR.hpp"

#include "edginc/KOStubRule.hpp"
#include "edginc/EqLinkCashFlow.hpp"
#include "edginc/KOLiborLeg.hpp"
#include "edginc/KOFixedLeg.hpp"
#include "edginc/IFinalPerf.hpp"

#include "edginc/EventResults.hpp"
//#include "edginc/RateTree.hpp"


DRLIB_BEGIN_NAMESPACE

typedef enum{NO_SMOOTHING=0, INCREASE=1, DECREASE=-1} TSmoothingType;

class CallableEquityKOSwap: public Generic1Factor,
                          virtual public FDModel::IIntoProduct,
                          virtual public LastSensDate,
                          virtual public LegalTerms::Shift,
                          virtual public BarrierBreach::IEventHandler,
                          virtual public KnownCashflows::IEventHandler 
{
public:
    static CClassConstSP const TYPE;  

    typedef enum{AsSchedule, AsInstSettle} TpayTimingRule;
    typedef enum{Alive, UpHit, LowHit, Called} TDeadType;

    virtual void Validate();

    //------------------------------------------//
    //  Class for Callable                 //
    //------------------------------------------//
    class CallSchedule: public CObject {
    public:
        static CClassConstSP const TYPE;  
        
        bool    isPuttable;       
        ScheduleSP callLevels;
        //DateTimeArray notification;
        //CashFlowArray schedule;

        CallSchedule(); 

        void validatePop2Object() ;

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);
        
        static IObject* defaultCallSchedule();
    };

    typedef smartPtr<CallSchedule> CallScheduleSP;

    //------------------------------------------//
    //  Class for Barrier (Trigger)        //
    //------------------------------------------//
    class Trigger: public CObject {
    public:
        static CClassConstSP const TYPE;  
        
        bool    isUp;       
        ScheduleSP barrier;
        ScheduleSP ecoBarrier;
        bool        intraDayMonitor;
        bool        hasRebate;
        ScheduleSP  rebate;
        DateTimeArray paymentDates;         // use specific payment date rather than instSettle
        string      smoothingType;
        double      peakDelta;
        
        Trigger();

        Trigger(const bool isUp,const ScheduleSP barrier,
            const ScheduleSP ecoBarrier, const bool intraDayMonitor,
            const bool hasRebate, const ScheduleSP  rebate,
            const string smoothingType, const double peakDelta, DateTimeArray paymentDates);

        void validatePop2Object() ;

        // convert input schedule to the array along to the time line
        void setStepRebatePVed(const DateTimeArray timeLine,       // (I) Date Array to be based on.
                          const vector<bool>& isKO,            // (I) is the step KO.
                          const YieldCurveConstSP discount,
                          const InstrumentSettlementSP  instSettle,
                          const CAssetWrapper asset,
                          vector<double>& reb) const;          // (O) the rebete value

        // convert barrier schedule to absolute level along with the time line
        void setStepBarrier(const DateTimeArray timeLine,      // (I) Date Array to be based on.
                            const bool          isUP,          // (I) is Up Barrier or not.  
                            const double        scale,         // (I) reference Level (convert to absolute Lvl).                             
                            const vector<bool>& isKO,          // (I) is the step KO.                            
                            vector<double>& bar) const;        // (O) the rebete value

        // check is already touched before.
        bool isTriggered(const CashFlowArray samples,const DateTime valDate, 
                         const double spot, const double refLevel, DateTime* hitDate) const;

        TSmoothingType getSmoothType() const;

        bool useEcoBar();

        // a barrier event maker 
        void makeBarEvnt(DateTime eventDate, DateTime breachDate, EventResults* events, 
                         CAssetWrapper asset, double initialSpot) const;
    
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);
        
        static IObject* defaultTrigger();
    };

    typedef smartPtr<Trigger> TriggerSP;

    // return average out closed form price.
//    double avgValue(const AvgOutPerfSP aop, double const spot, double SpotRef, const DateTime& when) const ;
    
    // return average out price, only after all simulation date.
//    double avgValue(const AvgOutPerfSP aop, double SpotRef) const ;

    // return intrinsic value when avgOut dates is one day sample.
//    double avgIntrinsic(const AvgOutPerfSP aop,double spot,double SpotRef, DateTime when) const;

    // return average out start date for tree.
//    DateTime getAvgOutStart() const;

    DateTime getStartDate() const;

    bool isMultiState() const;

    CallableEquityKOSwap(//Generic1Fcator
                        const DateTime                valueDate,
                        const bool                    oneContract,
                        const double                  notional,
                        const double                  initialSpot,
                        const string                  ccyTreatment,
                        const InstrumentSettlementSP  instSettle,
                        const InstrumentSettlementSP  premiumSettle,
                        const CAssetWrapper           asset,
                        const YieldCurveWrapper       discount,
                        // for CKS
                        const bool hasEqLeg,
                        EqCpnKOMakerSP koELC,
                        IFinalPerfMakerSP  finalPerfMaker,         
                        const double       redemption,    
                        const bool hasFixedLeg,
                        const bool hasFloater,
                        KOFixedLegMaker  koFixedLeg,
                        KOLiborLegMaker  koFloater,
                        const bool hasUpTrigger,
                        const bool hasLowTrigger,
                        Trigger     upTrigger,
                        Trigger     lowTrigger,
                        const bool hasCallSchedule,
                        CallSchedule callSchedule,    
                        const bool          isAvgIn,
                        CashFlowArray avgIn,
                        const bool  useIniSpotAsRefLvl,
                        CashFlowArray samples,
                        int numOfAvgInSV,
                        InstrumentSettlementSP instSettleForFP
#ifdef  TREE_THETA_CAP
                        ,
                        const double thetaSmoothThrehold,
                        const bool isPositiveThetaSmooth
#endif
                        );


public:
    friend class CallableEquityKOSwapHelper;
    friend class CallableEquityKOSwapProd;
    friend class CallableEquityKOSwapFDProd;
   
    CallableEquityKOSwap();


    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const Model *model);

    // make a LN vol request - not in real use as prefer LV tree
    CVolRequestLN* makeVolRequest() const ;
    
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                        const Model*      model);

    DateTime endDate(const Sensitivity* sensControl) const;

    bool sensShift(Theta* shift) ;

    /** for LegalTerms::Shift */
    virtual bool sensShift(LegalTerms* shift);
    
    /** extra output requests */
    void recordOutputRequests(Control* control, 
                              Results* results) const ;  

   
    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;
 
    /** Get the asset and discount market data */
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market) ;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const;

    bool isAlreadyHit(DateTime* hitDate, TDeadType &deadType) const;

    CashFlowArraySP getKnownCashFlow(const DateTime hitDate, const TDeadType deadType) const;

    DateTime getLastSimDate() const;

    DateTimeArray getCritDates() const;

    //------------------------------------------//    
    // for event handling
    //------------------------------------------//    
    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const; 

    // KnownCashflows::IEventHandler interface
    void getEvents(const KnownCashflows* flows,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const ;

private:
    // for access from product class.
    // internal member, as those variable type are not good in IMS.
    SampleListSP avgInSL;
    
    //AvgOutPerfSP aop;
    
    //------------
    // field
    // Equity Linked Coupon Parts
    bool hasEqLeg;
    EqCpnKOMakerSP koELC;

    // Final Option
    IFinalPerfMakerSP  finalPerfMaker;     //AvgOutPerfMakerSP  avgOutPerf;         
    double           redemption;    
    InstrumentSettlementSP  instSettleForFP;  // instrument settle for Final Perf & redemption @ final.

    // Swap Part
    bool hasFixedLeg;
    bool hasFloater;
    KOFixedLegMaker  koFixedLeg;
    KOLiborLegMaker  koFloater;

    // Barrier Part
    bool hasUpTrigger;
    bool hasLowTrigger;
    Trigger     upTrigger;
    Trigger     lowTrigger;
    
    // Callable Part
    bool hasCallSchedule;
    CallSchedule callSchedule;    

    // Average In
    bool          isAvgIn;
    CashFlowArray avgIn;    
    bool          useIniSpotAsRefLvl; // $unregistered
    CashFlowArray    samples;         
    int           numOfAvgInSV;
    //------------

    mutable CashFlowArraySP knownCFL; // $unregistered

#ifdef  TREE_THETA_CAP
    // for theta smoothing
    double  thetaSmoothThrehold;    
    bool    isPositiveThetaSmooth;
#endif

};
typedef smartPtr<CallableEquityKOSwap>       CallableEquityKOSwapSP;

/****************************************************************************************************************/
/////////////////////////////////////////////////////////
//     private class for all tree/FD product
//     new state variable interface for any num of factors
/////////////////////////////////////////////////////////
class CallableEquityKOSwapFDProd : public LatticeProdEDRIns{
protected:
    const CallableEquityKOSwap*     cks;
    int                 numOfPrice;

public:    
    CallableEquityKOSwapFDProd(const CallableEquityKOSwap* cks, FDModel* m);

    ~CallableEquityKOSwapFDProd();
    
    class StateInfo 
    {
    public:           
        double                    refLevel;
        IFinalPerfSP              finalPerf;
        vector<double>            stepUpBarrier;
        vector<double>            stepDnBarrier;

        // barrier at a tree slice after discrete monitor adjustment
        double                    adjUpBarrier; 
        double                    adjDnBarrier; 

        // for keeping insert node level;
        bool foundUpLvl;
        bool foundDnLvl;
        IntArray iLstDnNode;
        IntArray iFstUpNode;
        double prePriceAtUpBar;
        double prePriceAtDnBar;

        StateInfo(){};
    };

    // Set Boolean Array along time steps, which indicate the time step is 
    // corresponding to the event date of schedule.
    // Jus using SetStepExercise is not good 
    // because maturity date automatically true or last event date will be missed when it's not = lastNode. 
    void SetStepEvent(vector<bool>& stepIsEvent, ScheduleConstSP sched) ;

    void SetPerfIndex(const DateTimeArray timeLine, 
                      const DateTimeArray observDates);

    void CalcPayoff(int step, const TreeSlice & spot, const vector< TreeSliceSP > & price); 

    void CalcPayoff_oper(int step, const TreeSlice & spot, const vector< TreeSliceSP > & price); 


    /** CallableEquityKOSwap, quanto or struck */
    virtual string getCcyTreatment() const 
    {
        return cks->ccyTreatment;
    }
    
    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd();

    void         InitProdFD1D();

    void         InitFD1D(Control*    control) const;


    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, 
                        FDProduct::UpdateType type);

    virtual void preCalc(int step);

    void postCalc(int step, const TreeSlice & spot, const vector< TreeSliceSP > & price); 

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, 
                              YieldCurveConstSP disc, 
                              Results* results);

    virtual DateTime getStartDate() const
    {
        return cks->fwdStarting ? cks->startDate : cks->valueDate;
    }

    // !!! this is to work with current tree which uses this to set up timeline only */
    virtual CVolRequestConstSP GetLNRequest() const;

    bool    priceDeadInst(Control * control,Results * results);

private:
    EqCpnKOMakerSP koELC;

    KOFixedLegSP koFix;
    KOLiborLegSP koFlt;

    vector<bool>              stepIsUpKO;     
    vector<bool>              stepIsDnKO;     
    vector<double>            stepUpRebate;
    vector<double>            stepDnRebate;
    vector<double>            stepCallLevel;
    vector<bool>              stepIsCallable;

    vector<double>            stepFixedFee;
    vector<double>            stepFloatFee;
    vector<double>            stepKOCpns;
    vector<double>            stepEqCpnFactor;

    vector<StateInfo>         State;
    int                       numState;
    
    int                  baseNP;

    DateTime                  matDate;
    DateTime                  avgOutStart;
    int                       iStepAvgOutStart;
    int                       iStepAvgInEnd;
    DoubleArray               pdf_AvgIn;
    
    EqCpnKOLeg eqCpnKOLeg;        // plain without any KOStubRule
    // Operational aggregation and performance calcs
    SimpleDoubleArray         eqPerfs;    // equity performances to be processed and aggregated for eqLinkCpn

    IDoubleArrayModifierSP    performance;
    vector<int> stepPerfIndex;

    TSmoothingType upSmoothType;
    TSmoothingType lowSmoothType;

    // used at RefinePrice...    
    double                      notPaidValues;       //just not paid yet, but never cancelled.

    DoubleArray db_UpBarSLow;     // barrier smoothing lower level.
    DoubleArray db_UpBarSHi;      // barrier smoothing higher level.
    DoubleArray db_UpBarPLow;     // barrier smoothing lower price.
    DoubleArray db_UpBarPHi;      // barrier smoothing higher price.
    DoubleArray db_DnBarSLow;     // barrier smoothing lower level.
    DoubleArray db_DnBarSHi;      // barrier smoothing higher level.
    DoubleArray db_DnBarPLow;     // barrier smoothing lower price.
    DoubleArray db_DnBarPHi;      // barrier smoothing higher price.
    CashFlowArraySP db_UpBarDelta;     // barrier smoothing lower level.
    CashFlowArraySP db_DnBarDelta;     // barrier smoothing lower level.

    // to allow switching between original and using slice operators update
    typedef void ( CallableEquityKOSwapFDProd::*CalcPayoff_FUNC )(
        int step, const TreeSlice & s, const vector< TreeSliceSP > & price );
    CalcPayoff_FUNC calcPayoff;

//    ZeroBondProdSP   zeroProd;
    
};

/// ------------  end of new FD product -----------

///////////////////////////////////////////////////////
////   for Barrier Events                 //////////////
///////////////////////////////////////////////////////
// a helper that makes a barrier event once one is detected, for one factor Barrier.
class Barrier1F 
{
public:
    // DATA constructor
    Barrier1F(ScheduleSP barrier,
                bool  isUp,
                bool isCont,
                DateTime breachDate,                                  
                bool isBreached,
                CAssetWrapper asset,
                double initialSpot);
                        
    void makeBarrierEvent(EventResults* events,
                        const DateTime& eventDate,
                        string barrDesc,
                        string barrierType);
private:
    ScheduleSP barrier;
    bool  isUp;
    bool isCont;
    DateTime breachDate;
    bool isBreached;
    CAssetWrapper asset;
    double initialSpot;
};



DRLIB_END_NAMESPACE
#endif
