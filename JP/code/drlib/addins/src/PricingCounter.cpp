//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PricingCounter.cpp
//
//   Description : EDRAction to count pricings for an instrumetn given a ctrl
//
//   Author      : Andrew J Swain
//
//   Date        : 6 November 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

class PricingCounter: public CInstrument,
                      public ISensitiveStrikes,
                      public LastSensDate,
                      virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;  
    
    virtual DateTime getValueDate() const {
        return inst->getValueDate();
    }

    /** Get the asset and discount market data */
    virtual void GetMarket(const IModel*       model, 
                           const CMarketDataSP market) {
        (const_cast<IModel*>(model))->getInstrumentAndModelMarket(market.get(), inst.get());
    }
    
    virtual void Validate() {
        inst->Validate();
    }


    /** Returns the name of the instrument's discount currency */
    string discountYieldCurveName() const {
        return inst->discountYieldCurveName();
    }


private:
    class BogusModel: public CModel {
    public:
        static CClassConstSP const TYPE;
        
        CIntArraySP pricings; // $unregistered
        IModelConstSP model; // $unregistered

        /** Simple constructor */
        BogusModel(IModelConstSP model):
            Model(TYPE),
            pricings(new IntArray(1)),
            model(model)
        {}

        /** calculate single price and store result in CResult */
        virtual void Price(CInstrument*  instrument, 
                           CControl*     control, 
                           CResults*     results) {
            (*pricings)[0]++;
            results->storePrice(0.0, "");
        }        

        /** Whether to enable RiskMapping when computing sensitivities for
         * instruments priced using this model */
        virtual IModel::WantsRiskMapping wantsRiskMapping() const {
            return model->wantsRiskMapping();
        }

        /** Override clone method to reference pricings */
        IObject* clone() const {
            BogusModel* m = new BogusModel(model);
            m->pricings = this->pricings;
            return m;
        }    
    };

    friend class PricingCounterHelper;
    friend class BogusModelHelper;

    PricingCounter(): CInstrument(TYPE) {};

    PricingCounter(const PricingCounter& rhs);
    PricingCounter& operator=(const PricingCounter& rhs);

    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model) {
        ISensitiveStrikes* vm = dynamic_cast<ISensitiveStrikes*>(inst.get());
        if (vm) {
            return (vm->avoidVegaMatrix(this->model.get()));
        }
        return true;
    }

    /** Returns all strikes the PricingCounter is sensitve to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        ISensitiveStrikes* vm = dynamic_cast<ISensitiveStrikes*>(inst.get());
        if (vm) {
            return (vm->getSensitiveStrikes(outputName, this->model.get()));
        }
        DoubleArraySP strikes(new DoubleArray(0));
        strikes->push_back(0);
        return (strikes);
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const {
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(inst.get());
        if (lsd) {
            return (lsd->endDate(sensControl));
        }
        // ugh !
        MaturityPeriod ages("50Y");
        return ages.toDate(inst->getValueDate());
    }
           
    // EdrAction 
    IObjectSP run() {
        BogusModel bogus(model);

        CResultsSP results(bogus.go(CInstrumentSP::attachToRef(this),
                                    ScenarioSP(),
                                    ctrl,
                                    market));

        return IObjectSP(CInt::create((*bogus.pricings)[0]));
    }

private:
    IModelSP      model;
    CInstrumentSP inst;
    CControlSP    ctrl;
    CMarketDataSP market;  
};

class PricingCounterHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PricingCounter, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate); 
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultPricingCounter);
        FIELD(model, "model");
        FIELD(inst, "inst");
        FIELD(ctrl, "ctrl");
        FIELD(market, "market");
    }

    static IObject* defaultPricingCounter(){
        return new PricingCounter();
    }
};

CClassConstSP const PricingCounter::TYPE = CClass::registerClassLoadMethod(
    "PricingCounter", typeid(PricingCounter), PricingCounterHelper::load);


class BogusModelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PricingCounter::BogusModel, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultBogusModel);
        // no fields
    }

    static IObject* defaultBogusModel(){
        return new PricingCounter::BogusModel(IModelSP());
    }
};

CClassConstSP const PricingCounter::BogusModel::TYPE = CClass::registerClassLoadMethod(
    "BogusModel", typeid(PricingCounter::BogusModel), BogusModelHelper::load);

/* for class loading */
bool PricingCounterLoad() {
    return (PricingCounter::TYPE != 0);
}


DRLIB_END_NAMESPACE
