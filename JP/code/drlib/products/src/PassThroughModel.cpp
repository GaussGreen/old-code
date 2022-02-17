//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PassThroughModel.cpp
//
//   Description : Model that wraps another model - 
//                 similar to ClosedForm in that asks instrument to drive the
//                 process, but routes pricing etc through wrapped model
//
//   Author      : Andrew J Swain
//
//   Date        : 19 August 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/PassThroughModel.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

PassThroughModel::PassThroughModel(CClassConstSP clazz): Model(clazz){}

PassThroughModel::PassThroughModel(): Model(TYPE){}

// price using wrapped model - gateway for instrument to get access
void PassThroughModel::calculate(CInstrument* instrument, 
                                 Control*     control, 
                                 Results*     results) {
    model->Price(instrument, control, results);
}


/** calculate single price and store result in Results */
void PassThroughModel::Price(Instrument*  instrument, 
                             Control*     control, 
                             Results*     results){
    static const string method = "PassThroughModel::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support PassThroughModel::IntoProduct");
    }
    IProduct* product = 0;
    try{
        product = intoProd->createProduct(this);
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}
  

MarketObjectSP PassThroughModel::GetMarket(const MarketData*    market,
                                           const string&        name,
                                           const CClassConstSP& type) const {
    return model->GetMarket(market, name, type);
}

void PassThroughModel::getMarket(const MarketData* market, IInstrumentCollectionSP instruments){
    model->getMarket( market, instruments );
}

IModel::WantsRiskMapping PassThroughModel::wantsRiskMapping() const {
    return model->wantsRiskMapping();
}

class PassThroughModelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PassThroughModel, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultPassThroughModel);
        FIELD(model, "model");
    }
    // for PassThroughModel::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(PassThroughModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultPassThroughModel(){
        return new PassThroughModel();
    }
};

CClassConstSP const PassThroughModel::TYPE = CClass::registerClassLoadMethod(
    "PassThroughModel", typeid(PassThroughModel), PassThroughModelHelper::load);

CClassConstSP const PassThroughModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("PassThroughModel::IIntoProduct",
                                    typeid(PassThroughModel::IIntoProduct),
                                    PassThroughModelHelper::loadIntoProduct);

// a series of bogus flavours of PassThroughModel for IMS as it has no
// polymorphic support for the "model" element

// this is the ClosedFormLN version
class PassThroughCFLN: public PassThroughModel{
public:
    static CClassConstSP const TYPE;
private:
    PassThroughCFLN(): PassThroughModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PassThroughCFLN, clazz);
        SUPERCLASS(PassThroughModel);
        EMPTY_SHELL_METHOD(defaultPassThroughCFLN);
    }

    static IObject* defaultPassThroughCFLN(){
        return new PassThroughCFLN();
    }
};
CClassConstSP const PassThroughCFLN::TYPE = CClass::registerClassLoadMethod(
    "PassThroughCFLN", typeid(PassThroughCFLN), load);

// this is the Tree1fLN version
class PassThroughTR1F: public PassThroughModel{
public:
    static CClassConstSP const TYPE;
private:
    PassThroughTR1F(): PassThroughModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PassThroughTR1F, clazz);
        SUPERCLASS(PassThroughModel);
        EMPTY_SHELL_METHOD(defaultPassThroughTR1F);
    }

    static IObject* defaultPassThroughTR1F(){
        return new PassThroughTR1F();
    }
};
CClassConstSP const PassThroughTR1F::TYPE = CClass::registerClassLoadMethod(
    "PassThroughTR1F", typeid(PassThroughTR1F), load);

// MC Implied
class PassThroughMCIM: public PassThroughModel{
public:
    static CClassConstSP const TYPE;
private:
    PassThroughMCIM(): PassThroughModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PassThroughMCIM, clazz);
        SUPERCLASS(PassThroughModel);
        EMPTY_SHELL_METHOD(defaultPassThroughMCIM);
    }

    static IObject* defaultPassThroughMCIM(){
        return new PassThroughMCIM();
    }
};
CClassConstSP const PassThroughMCIM::TYPE = CClass::registerClassLoadMethod(
    "PassThroughMCIM", typeid(PassThroughMCIM), load);

// MC Lognormal
class PassThroughMCLN: public PassThroughModel{
public:
    static CClassConstSP const TYPE;
private:
    PassThroughMCLN(): PassThroughModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PassThroughMCLN, clazz);
        SUPERCLASS(PassThroughModel);
        EMPTY_SHELL_METHOD(defaultPassThroughMCLN);
    }

    static IObject* defaultPassThroughMCLN(){
        return new PassThroughMCLN();
    }
};
CClassConstSP const PassThroughMCLN::TYPE = CClass::registerClassLoadMethod(
    "PassThroughMCLN", typeid(PassThroughMCLN), load);

DRLIB_END_NAMESPACE
