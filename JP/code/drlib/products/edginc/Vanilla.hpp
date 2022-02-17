//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Vanilla.hpp
//
//   Description : Vanilla instrument
//
//   Author      : Andre X Segger
//
//   Date        : 23 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VANILLA_HPP
#define EDR_VANILLA_HPP
#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/FD1FGeneric.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/VanillaCreditSupport.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/VegaMatrixLite.hpp"

DRLIB_BEGIN_NAMESPACE

/** Vanilla instrument */
class PRODUCTS_DLL CVanilla: public CInstrument, public CClosedFormLN::IIntoProduct, 
                public ISupportVegaMatrixLite,
				public FDModel::IIntoProduct,
                public Theta::IShift,
                public ISensitiveStrikes,
                public LastSensDate,
                public IMCIntoProduct,
                public FD1F::IIntoProduct,
                public FD1FGeneric::IIntoProduct,
                public FourierEngine::IIntoProduct,
                public CreditSupport::Interface,
                public ITaxableInst::Basic {
public:
    static CClassConstSP const TYPE;

    /** Class that creates a DeltaImpledStrike calculator */
    class PRODUCTS_DLL DeltaImpliedStrikeMakerLN: public CObject, 
                                     virtual public IDeltaToStrikeMaker {
    public:
        static CClassConstSP const TYPE;

        /** Factory method */
        virtual const IDeltaToStrike* make(const DateTime&             valueDate,
                                           const DateTime&             maturityDate,
                                           bool                        isCall,
                                           const CAsset*               asset,
                                           const YieldCurve*           discount,
                                           const InstrumentSettlement* settle,
                                           double                      deltaShiftSize,
                                           double                      tgtDelta,
                                           const string&               volType,
                                           bool                        allowNegativeFwdVar);
    
    private:
        IModelSP model; // $unregistered
        
        /** Default constructor */
        DeltaImpliedStrikeMakerLN();

        /** Empty shell */
        static IObject* defaultDeltaImpliedStrikeMakerLN();
        
        /** Load method */
        static void load(CClassSP& clazz);
    };
    
    /** Class that lets you imply a strike from a delta level
        for a simple European vanilla option */
    class PRODUCTS_DLL DeltaImpliedStrike: virtual public IDeltaToStrikeMaker::IDeltaToStrike {
    public:
        struct PRODUCTS_DLL Helper{
            static void getMarket(const MarketData*        market,
                                  const IModel*             model,
                                  CAssetWrapper&           asset,
                                  YieldCurveWrapper&       discount,
                                  InstrumentSettlement*    settle);
        };

        DeltaImpliedStrike(const IModel&                model,
                           const DateTime&             valueDate,
                           const DateTime&             maturityDate,
                           bool                        isCall,
                           const CAsset*               asset,
                           const YieldCurve*           discount,
                           const InstrumentSettlement* settle,
                           double                      deltaShiftSize,
                           double                      tgtDelta);

        /** Given lower and upper bounds, and given a bracketer and
            tolerance number, calculate strike implied by delta level */
        double calcStrike(double             lowerStrike,
                          double             upperStrike,
                          double             strikeAbsAcc) const;
    
        ~DeltaImpliedStrike();
        class Imp;
    private:
        DeltaImpliedStrike(const DeltaImpliedStrike& rhs);
        DeltaImpliedStrike& operator=(const DeltaImpliedStrike& rhs);
        auto_ptr<Imp> me;
    };
    friend class DeltaImpliedStrike::Imp;
        
    virtual void avoidVegaMatrixLite(const IModel* model);
    
    static double priceSpread(const DateTime& valueDate,
                              const DateTime& startDate,
                              const DateTime& matDate,
                              bool isCall,
                              bool fwdStarting,
                              bool oneContract,
                              double notional,
                              double initialSpot,
                              double lowStrike,
                              double highStrike,
                              const InstrumentSettlement* instSettle,
                              const Asset* asset,
                              const YieldCurve* discount);

    static double priceBS(const DateTime& valueDate,
                          const DateTime& startDate,
                          const DateTime& matDate,
                          bool isCall,
                          bool fwdStarting,
                          bool oneContract,
                          double notional,
                          double initialSpot,
                          double strike,
                          const InstrumentSettlement* instSettle,
                          const Asset* asset,
                          const YieldCurve* discount);

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** retrieve market data needed by Vanilla - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Implementation of MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** Implementation of FourierEngine::IntoProduct interface */
    virtual FourierProduct* createProduct(const FourierEngine* model) const;

    /** create a fd payoff product */
    virtual FD1F::IProduct* createProduct(FD1F* model) const;

    /** create a fd payoff product for DDE */
    virtual FD1FGeneric::IProduct* createProduct(FD1FGeneric* model) const;

	//////////////////////// **************** /////////////////
    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual CSensControl* AlterControl(const IModel*       modelParams,
                                       const CSensControl* sensControl) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);
  
    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    /** make a simple started vanilla ready for pricing */
    static CVanilla* make(
        const DateTime&             valueDate,
        bool                        isCall,
        bool                        american,
        const Schedule*             exerciseSchedule,
        const CAsset*               asset,
        const YieldCurve*           discount,
        const InstrumentSettlement* settle,
        int                         noExerciseWindow);

    /** make a simple started vanilla ready for pricing from CAssetWrapper */
    static CVanilla* make(
        const DateTime&             valueDate,
        bool                        isCall,
        bool                        american,
        bool                        oneContract,
        double                      notional,
        double                      initialSpot,
        const Schedule*             exerciseSchedule,
        const CAssetWrapper&        asset,
        const YieldCurveWrapper&    discount,
        const InstrumentSettlement* settle,
        int                         noExerciseWindow);

    /** make a vanilla ready for pricing */
    static CVanilla* make(
        const DateTime&             valueDate,
        bool                        isCall,
        bool                        american,
        const Schedule*             exerciseSchedule,
        bool                        fwdStarting,
        const DateTime&             startDate,
        bool                        oneContract,
        double                      notional,
        double                      initialSpot,
        const CAsset*               asset,
        const string&               ccyTreatment, 
        const YieldCurve*           discount,
        const InstrumentSettlement* settle);


    virtual CreditSupportSP createCreditSupport(CMarketDataSP market){
                return CreditSupportSP(new VanillaCreditSupport(this, market));}

    /** for ITaxableInst::Basic */
    const DateTime getFinalPaymentDate() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    friend class CVanillaHelper;
    friend class CVanillaClosedForm;
    friend class Vanilla1fProd;
    friend class Vanilla1fProdDDE;
    friend class VanillaMC;
    friend class VanillaMCSV;
    friend class VanillaFP;
    friend class VanillaCreditSupport;

    friend class FastQuoteEnv;

    friend class VanillaTSO;
    friend class VanillaFDProd;

    CVanilla();
    CVanilla(const CVanilla& rhs);
    CVanilla& operator=(const CVanilla& rhs);
    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

protected:
    CVanilla(CClassConstSP clazz);

    DateTime                valueDate;
    DateTime                startDate;
    InstrumentSettlementSP  premiumSettle;

    bool                    isCall;
    bool                    fwdStarting;
    bool                    canExerciseEarly;
    ScheduleSP              exerciseSchedule;

    bool                    oneContract;
    double                  notional;
    double                  initialSpot;
    double                  spotAtMaturity;

    DateTime                dateExercised;
    bool                    isExercised;

    CAssetWrapper           asset;
    string                  ccyTreatment;
    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;

    int                     noExerciseWindow;
};

typedef CVanilla Vanilla;
typedef smartPtr<CVanilla> CVanillaSP;

class PRODUCTS_DLL VanillaISAP: public CObject, 
                   public FourierEngine::ISAP {
public:
    static CClassConstSP const TYPE;
    friend class VanillaISAPHelper;
    friend class VanillaFP;
    
    void validatePop2Object(){
        static const string method = "VanillaISAP::validatePop2Object";    
        if(!Maths::isPositive(payoffToProcFreqBoundWeight) || !Maths::isPositive(1.0 - payoffToProcFreqBoundWeight)){
            throw ModelException(method,
                                 "payoffToProcFreqBoundWeight must be strictly between 0.0 and 1.0.");
        }
    }

protected:
    VanillaISAP():
    CObject(TYPE),
    useOneIntegral(true),
    payoffToProcFreqBoundWeight(0.9){}

private:
    bool useOneIntegral;
    double payoffToProcFreqBoundWeight;
};

typedef smartPtr<VanillaISAP> VanillaISAPSP;
typedef smartConstPtr<VanillaISAP> VanillaISAPConstSP;


DRLIB_END_NAMESPACE
#endif
