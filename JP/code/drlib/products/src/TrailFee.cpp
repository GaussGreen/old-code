//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TrailFee.cpp
//
//   Description : trail fee - difference in price of underlying instrument
//                 with and without additional borrowing costs corresponding
//                 to fee levels of each fund
//
//   Author      : Andrew J Swain
//
//   Date        : 12 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNAssetAndImnt.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/BorrowParallelShift.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/ContangoCommodity.hpp"
#include "edginc/ContangoRhoParallel.hpp"

DRLIB_BEGIN_NAMESPACE

class TFModel: public CModel {
public:
    static CClassConstSP const TYPE;
    friend class TFModelHelper;

    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const {
        return model->GetMarket(market, name, type);
    }

    MarketDataFetcherSP createMDF() const {
        return model->getMDF();
    }
    
    virtual void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments) {
        model->getMarket(market, instruments);
    }

    /** the class that the product must be able to create */
    class IProduct{
    public:
        virtual void price(TFModel* model,
                           Control* control, 
                           Results* results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        friend class TFModelHelper;
        virtual IProduct* createProduct(TFModel* model) const = 0;
    };

    virtual void Price(CInstrument*  instrument, 
                          CControl*     control, 
                       CResults*     results){
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException("TFModel::Price", "Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support TFModel::IntoProduct");
        }
        IProduct*  product = 0;
        try{
            // cast to TFModel::IIntoProduct
            IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
            // create the product
            product = intoProd.createProduct(this);
            // and the invoke the pricing
            product->price(this, control, results);
        } catch (exception& e){
            delete product;
            throw ModelException(e, "TFModel::Price");
        }
        delete product;
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const {
        return model->wantsRiskMapping();
    }

    // fields
    IModelSP model;
protected:
    TFModel(CClassConstSP clazz);
    
private:
    /* for reflection */
    TFModel();
};

TFModel::TFModel(CClassConstSP clazz): CModel(clazz) {}

// for reflection
TFModel::TFModel():CModel(TYPE) {}

class TFModelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TFModel, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultTFModel);
        FIELD(model, "model");
    }

    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(TFModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultTFModel(){
        return new TFModel();
    }
};

CClassConstSP const TFModel::TYPE = CClass::registerClassLoadMethod(
    "TFModel", typeid(TFModel), TFModelHelper::load);


CClassConstSP const TFModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("TFModel::IIntoProduct",
                                    typeid(TFModel::IIntoProduct), 
                                    TFModelHelper::loadIntoProduct);


/** A MonteCarlo that can be captured in Pyramid using 
     current IMS as a model inside a model */
class MonteCarloImpliedTF: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloImpliedTF():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloImpliedTF, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloImpliedTF);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloImpliedTF(){
        return new MonteCarloImpliedTF();
    }
};

CClassConstSP const MonteCarloImpliedTF::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloImpliedTF", typeid(MonteCarloImpliedTF), load);

// ditto for the top-level TrailFee model itself
class TrailFeeMonteCarlo: public TFModel {
public:
    static CClassConstSP const TYPE;
    
private:
    TrailFeeMonteCarlo():TFModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(TrailFeeMonteCarlo, clazz);
        SUPERCLASS(TFModel);
        EMPTY_SHELL_METHOD(defaultTrailFeeMonteCarlo);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultTrailFeeMonteCarlo(){
        return new TrailFeeMonteCarlo();
    }
};

CClassConstSP const TrailFeeMonteCarlo::TYPE = 
CClass::registerClassLoadMethod(
    "TrailFeeMonteCarlo", typeid(TrailFeeMonteCarlo), load);

// LN flavour
class TrailFeeMonteCarloLN: public TFModel {
public:
    static CClassConstSP const TYPE;
    
private:
    TrailFeeMonteCarloLN():TFModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(TrailFeeMonteCarloLN, clazz);
        SUPERCLASS(TFModel);
        EMPTY_SHELL_METHOD(defaultTrailFeeMonteCarloLN);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultTrailFeeMonteCarloLN(){
        return new TrailFeeMonteCarloLN();
    }
};

CClassConstSP const TrailFeeMonteCarloLN::TYPE = 
CClass::registerClassLoadMethod(
    "TrailFeeMonteCarloLN", typeid(TrailFeeMonteCarloLN), load);

// the Trail Fee proper
class TrailFee: public GenericNAssetAndImnt,
                virtual public CClosedFormLN::IIntoProduct,
                virtual public TFModel::IIntoProduct,
                virtual public ISensitiveStrikes,
                virtual public LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void Validate() {
        static const string method = "TrailFee::Validate";
        try {
            if (fees.size() != assets.size()) {
                throw ModelException(method, 
                                     "fees must be same length as assets");
            }            
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model) {
        ISensitiveStrikes* vegaMatrixImnt = 
            dynamic_cast<ISensitiveStrikes*>(imnt.get());

        if (vegaMatrixImnt) {
            return vegaMatrixImnt->avoidVegaMatrix(model);
        }
        else {
            return true;  // can't do it
        }
    }
        

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
        if (avoidVegaMatrix(model)) {
            throw ModelException("TrailFee::getSensitiveStrikes", 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        ISensitiveStrikes& strikes = dynamic_cast<ISensitiveStrikes&>(*imnt);

        return strikes.getSensitiveStrikes(outputName, model);
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TrailFee, clazz);
        SUPERCLASS(GenericNAssetAndImnt);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(TFModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultTrailFee);
        FIELD(fees, "fees");
        FIELD(useFeeModel, "price using Trail fee model");
    }

    static IObject* defaultTrailFee(){
        return new TrailFee();
    }
    
private:
    friend class TrailFeePricer;
    friend class TFModel;

    TrailFee():GenericNAssetAndImnt(TYPE), useFeeModel(false) {}; 
    TrailFee(const TrailFee& rhs);
    TrailFee& operator=(const TrailFee& rhs);

    /** Get the asset and discount market data:
        overrides GenericNAssetAndImnt::GetMarket because we need to branch on
        useFeeModel */
    void GetMarket(const IModel*       feeModel, const CMarketDataSP market) {
        try {
            market->GetReferenceDate(valueDate);
            discount.getData(feeModel, market);
            instSettle->getMarket(feeModel, market.get());

            for (int i = 0; i < assets.size(); i++) {
                CAsset::getAssetMarketData(feeModel,
                                           market.get(), 
                                           CAsset::CCY_TREATMENT_NONE,
                                           discount.getName(), 
                                           assets[i]);
            }

            (useFeeModel ? const_cast<IModel*>(feeModel) : imntModel.get())->
                getInstrumentAndModelMarket(market.get(), imnt.get());
        }
        catch (exception& e) {
            throw ModelException(e, "GenericNAssetAndImnt::GetMarket");
        }
    }

    void price(IModel* feeModel,Control* control,CResults* results) const{
        static const string method = "TrailFee::price";

        try  {
            CControlSP ctrl(Control::makeFromFlags("", 0.0));
            IModel*    model = useFeeModel ? feeModel : imntModel.get();
            // for extra caution we'll clone the model
            IModelSP modelCopy1(copy(model));
            double value = 0.0;
            // 'Run' should route through control as appropriate
            ResultsSP baseResults(modelCopy1->Run(imnt.get(), ctrl.get()));
            double basePrice = baseResults->retrievePrice();

            // shift borrowing cost in each asset by fee amount
            IScenarioShiftArraySP shifts(
                new IScenarioShiftArray(fees.size()));
            int i;

            for (i = 0; i < fees.size(); i++) {
                // We can't "just" get the asset name as this may not match 
                // the name given to the borrow/contango curve, so pull it out via
                // SensMgr as that's guaranteed to give the right thing.
                // Need paranoia checks in case we've got a composite
                // asset (which component do we shift ?) or somehow get
                // a multi-name key for shifting
                CAssetSP asset(copy(assets[i].get()));
                SensMgr sensMgr(asset);
                OutputNameArrayConstSP names;

                // for commodity assets we shift the contango by -fee
                // as commodities have no borrow curve
                if (ContangoCommodity::TYPE->isInstance(assets[i].get())) {
                    smartPtr<ScalarShift> scalarShift(new ContangoRhoParallel(-fees[i]));
                    names = sensMgr.allNames(scalarShift.get());
                } else {
                    smartPtr<ScalarPerturbation> perturbation(new BorrowParallelShift(fees[i]));
                    names = sensMgr.allNames(perturbation.get());
                }

                if (names->size() > 1) {
                    throw ModelException(method,
                                         "asset " + assets[i]->getName() +
                                         " (" + 
                                         assets[i]->getClass()->getName() + 
                                         ") has more than one borrow/contango "
                                         "sensitive component");
                }
                 
                // if no name two cases
                // 1) fees with the asset: fail
                // 2) no fees: skip the asset
                if (names->empty()) {
                    if (!Maths::isZero(fees[i])) {
                        throw ModelException(method,
                                             "asset " + assets[i]->getName() +
                                             " (" + 
                                             assets[i]->getClass()->getName() + 
                                             ") has no borrow/contango "
                                             "sensitive component (it is a " + asset->getClass()->getName()+")");
                    }
                    else {
                        // skip the asset
                    }
                }
                else {
                    // see if it's gone horribly wrong
                    if ((*names)[0]->idCount() != 1) {
                        throw ModelException(method,
                                             "asset " + assets[i]->getName() +
                                             " has more than 1 name ??");
                    }
                    
                    smartPtr<IPerturbation> shift;
                    // for commodity assets we shift the contango by -fee
                    // as commodities have no borrow curve
                    if (ContangoCommodity::TYPE->isInstance(assets[i].get())) {
                        shift = smartPtr<IPerturbation>(new ContangoRhoParallel(-fees[i]));
                    } else {
                        shift = smartPtr<IPerturbation>(new BorrowParallelShift(fees[i]));
                    }

                    IScenarioShiftSP scenShift(
                        new ScenarioShift(shift, (*names)[0]->idGet(0)));
                    
                    (*shifts)[i] = scenShift;
                }
            }

            // the imnt AND the model may BOTH contain (potentially the same)
            // market data, so need to shift both
            CInstrumentSP tweaked(copy(imnt.get()));
            IModelSP modelCopy2(copy(model)); // extra safety

            // shift all borrowing costs before pricing
            for (i = 0; i < shifts->size(); i++) {
                // shift only if necessary
                if (!Maths::isZero(fees[i])) {
                    (*shifts)[i]->applyScenario(tweaked);
                    (*shifts)[i]->applyScenario(modelCopy2);
                }
            }
            // 'Run' should route through control as appropriate
            ResultsSP tweakResults(modelCopy2->Run(tweaked.get(), ctrl.get()));
            double tweakPrice = tweakResults->retrievePrice();
                         
            value = basePrice - tweakPrice;

            results->storePrice(value, discount->getCcy());

            if (control && control->isPricing()) {
                // In Pyramid, all flavours of SPI are set up with compute index 75
                // to trigger a certain overnight Merlin configuration. To ensure
                // that SPIs are set correctly, hardwire inside QLib rather than
                // rely on some external process.
                // Ditto any Trail Fee on an SPI
                OutputRequest*request = control->requestsOutput(OutputRequest::COMPUTE_INDEX);
                if (request) {
                    CClassConstSP spiType = CClass::forName("SyntheticPortfolioInsurance");
                    if (spiType->isInstance(imnt.get())) {
                        results->storeRequestResult(request, 75);
                    }
                }
            }                
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime endDate(const Sensitivity* sensitivity) const {
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(imnt.get());

        if (lsd) {
            // always use the imnt's model to figure out the end date
            // rationale - the instrument knows its end date when it prices
            // and this avoid problems where the instrument end date only
            // works with its own model (i.e. GenericNFBase, I thank you)
            SensitivityConstSP shift(sensitivity->spawn(imntModel.get()));

            return lsd->endDate(shift.get());
        }

        throw ModelException("TrailFee::endDate",
                             "underlying imnt (" + 
                             imnt->getClass()->getName() + 
                             ") does not implement LastSensDate");
    }

    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
    TFModel::IProduct* createProduct(TFModel* model) const;

private:
    DoubleArray fees;
    bool        useFeeModel;
};

CClassConstSP const TrailFee::TYPE = CClass::registerClassLoadMethod(
    "TrailFee", typeid(TrailFee), TrailFee::load);


/** private class */
class TrailFeePricer: virtual public CClosedFormLN::IProduct,
                      virtual public TFModel::IProduct {
private:
    const TrailFee* opt; // a reference

public:
    TrailFeePricer(const TrailFee* opt): opt(opt){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        opt->price(model, control, results);
    }
    void price(TFModel*  model,
               Control*  control, 
               CResults* results) const{
        opt->price(model->model.get(), control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* TrailFee::createProduct(CClosedFormLN* model) const
{
    return new TrailFeePricer(this);
}

TFModel::IProduct* TrailFee::createProduct(TFModel* model) const
{
    return new TrailFeePricer(this);
}


// for class loading 
bool TrailFeeLoad() {
    return (TrailFee::TYPE != 0);
}

DRLIB_END_NAMESPACE
