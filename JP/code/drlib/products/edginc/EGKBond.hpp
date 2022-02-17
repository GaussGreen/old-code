//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EGKBond.hpp
//
//   Description : As CallableEquityKOSwap but with limited input for EGK trades:
//
//   Date        : May 2005
//
//
//----------------------------------------------------------------------------

DRLIB_BEGIN_NAMESPACE
#ifndef EDG_EQK_BOND_HPP
#define EDG_EQK_BOND_HPP

/*****************************************************************************/

/// EGKBond product - as CKS but with limited interface
class PRODUCTS_DLL EGKBond: public Generic1Factor, 
              virtual public FDModel::IIntoProduct,
              virtual public LastSensDate {

public:
    //---------------------------------------------//
    //  Class for Simple Equity Linked CashFlow   //
    //---------------------------------------------//
    class PRODUCTS_DLL     SimpleEqLinkCashFlow: public CObject {
    public:
        static CClassConstSP const TYPE;  

        void validatePop2Object(){
            static const string method("CallableEquityKOSwap::EqLinkCashFlow::validatePop2Object");
            try {
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        };
        
        // data
        IDoubleArrayModifierMakerSP  eqLinkPerf;        // performace of each observation d
        DateTimeArray                observDates;      // eq link coupons observation Dates
        //DateTimeArray                paymentDates;      // eq link coupons payment Dates.
        //int                          offset;            // offset for maturity of eqCpn.

        SimpleEqLinkCashFlow():CObject(TYPE){}; 
   
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("EGK Bond's Equity Linked CashFlows");
            REGISTER(EGKBond::SimpleEqLinkCashFlow, clazz);
            SUPERCLASS(CObject);
            FIELD(eqLinkPerf,       "performace of each observation date");
            FIELD(observDates,        "each eq link coupons observation dates");
            EMPTY_SHELL_METHOD(defaultSimpleEqLinkCashFlow);
        };
        
        static IObject* defaultSimpleEqLinkCashFlow(){
            return new EGKBond::SimpleEqLinkCashFlow();
        };
    };

    //---------------------------------------------//
    //  Class for Trigger  (Simple as no Rebate)   //
    //---------------------------------------------//
    class PRODUCTS_DLL     SimpleTrigger: public CObject {
    public:
        static CClassConstSP const TYPE;  
        
        ScheduleSP barrier;
        ScheduleSP ecoBarrier;

        // for user support...
        DoubleArray bonusCoupons;
        string      smoothingType;
        double      peakDelta;
       
        SimpleTrigger():CObject(TYPE) {        
            bonusCoupons = DoubleArray(0);
            smoothingType = "NO_SMOOTHING";
            peakDelta = 0.0;
        };

        void validatePop2Object(){
            static const string method("EGKBond::SimpleTrigger::validatePop2Object");
            try {
                if (bonusCoupons.size() > 0){
                    if (bonusCoupons.size() != barrier->length())
                        throw ModelException(method, "bonusCoupons should have same size array to barDates.");
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        };
        
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("EGK Bond Simple Trigger");
            REGISTER(EGKBond::SimpleTrigger, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultSimpleTrigger);
            FIELD(barrier, "Trigger Schedule");
            FIELD(ecoBarrier, "Trigger Schedule in contract (Not used for pricing)");
            FIELD(bonusCoupons, "Add this cash amount to redemption, when KO occurs");
            FIELD_MAKE_OPTIONAL(bonusCoupons);
            FIELD(smoothingType, "smoothing type of trigger, NONE, UPPER or LOWER");
            FIELD_MAKE_OPTIONAL(smoothingType);
            FIELD(peakDelta, "smoothing size.  model will smooth so as to reduce the peak delta by this number");
            FIELD_MAKE_OPTIONAL(peakDelta);
        };
        
        static IObject* defaultSimpleTrigger(){
            return new EGKBond::SimpleTrigger();
        };
    };

    //---------------------------------------------//
    //  Data Class for Average In                  //
    //---------------------------------------------//
    class PRODUCTS_DLL     AverageIn: public CObject {
    public:
        static CClassConstSP const TYPE;  
        
        bool          useIniSpotAsRefLvl;
        CashFlowArray avgIn;
        
        AverageIn():CObject(TYPE) {};

        void validatePop2Object(){
            static const string method("EGKBond::AverageIn::validatePop2Object");
            try {
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        };
        
        int numAvgDates(){
            return avgIn.size();
        }

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("EGK Bond Avgerate In");
            REGISTER(EGKBond::AverageIn, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultAverageIn);
            FIELD(useIniSpotAsRefLvl, "true : use SpotAtStart value rather than sampled average value for reference Level");
            FIELD(avgIn,    "average-in schedule");
        };
        
        static IObject* defaultAverageIn(){
            return new EGKBond::AverageIn();
        };
    };


private: 
    /// fields ////////
    SimpleEqLinkCashFlow eqCpnSimple;
    AvgOutPerfMakerSP  avgOutPerf;         
    double           redemption;    
    FixedLegSP  fixedLeg;
    SimpleTrigger  trigger;

    // Average In
    bool          isAvgIn;
    AverageIn     averageIn;
    int           numOfAvgInSV; // $unregistered

    CashFlowArray samples;
    bool          isSkipCheckSamples;

    CallableEquityKOSwapSP       originalCKS;    
    IFinalPerfMakerSP finalPerfMaker;

    InstrumentSettlementSP myFinalSettle; // $unregistered

#ifdef  TREE_THETA_CAP
    // for theta smoothing
    double  thetaSmoothThrehold;    
    bool    isPositiveThetaSmooth;
#endif

public:
    static CClassConstSP const TYPE;
    
    // validation
    void validatePop2Object();

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market); 

    // Pass-through to internal CallableEquityKOSwap.
    // The EGK Bond itself has not acquired market data.
    void Validate();

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const ;

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const;

    bool priceDeadInstrument(CControl* control, CResults* results) const;
    
private: 
    friend class EGKBondHelper;
    friend class EGKBondProd;

    EGKBond():Generic1Factor(TYPE), numOfAvgInSV(1), isSkipCheckSamples(false)
#ifdef  TREE_THETA_CAP
                ,thetaSmoothThrehold(0.0),isPositiveThetaSmooth(false)
#endif
                {};// for reflection

    EGKBond(const EGKBond& rhs); // not implemented
    EGKBond& operator=(const EGKBond& rhs); // not implemented    
};



DRLIB_END_NAMESPACE
#endif

