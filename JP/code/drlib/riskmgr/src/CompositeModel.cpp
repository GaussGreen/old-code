//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompositeModel.cpp
//
//   Description : Class that wraps one or more CModels. To be used when the 
//                 pricing is delegated to the instrument and possibly requires 
//                 more than one model
//
//   Author      : Regis Guichard
//
//   Date        : 31 March 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/CompositeModel.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor takes type of vol to use */
CompositeModel::CompositeModel():
CModel(TYPE){}

/** throws exception if expected nb of models differ from actual */
void CompositeModel::checkNbModels(int    nbModelsExpected,
                                   string className) const{
    if (nbModelsExpected != models.size()){
        throw ModelException("CompositeModel::checkNbModels",
                             className + " expects " + Format::toString(nbModelsExpected)
                             + " model(s). Yet the actual nb of models is " 
                             + Format::toString(models.size()));
    }
}

/** returns the nber of models */
int CompositeModel::getNbModels() const{
    return models.size();
}

/** returns the i-th model */
IModelSP CompositeModel::getModel(int i) const{
    return models[i];
}

/** calculate single price and store result in CResult */
void CompositeModel::Price(CInstrument*  instrument, 
                           CControl*     control, 
                           CResults*     results){
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("CompositeModel::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support CompositeModel::IntoProduct");
    }
    IProduct*  product = 0;
    try{
        // cast to CompositeModel::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        product = intoProd.createProduct(this);
        // and the invoke the pricing
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(&e, "CompositeModel::Price");
    }
    delete product;
}

MarketObjectSP CompositeModel::GetMarket(const MarketData*    market,
                                         const string&        name,
                                         const CClassConstSP& type) const{
    throw ModelException("CompositeModel::GetMarket",
                         "Not supported");
}

IModel::WantsRiskMapping CompositeModel::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class CompositeModelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CompositeModel, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultCompositeModel);
        FIELD(models, "models");
    }

    // for CompositeModel::IIntoProduct
    static void IntoProduct_load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(Model::IModelIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultCompositeModel(){
        return new CompositeModel();
    }
};

CClassConstSP const CompositeModel::TYPE = CClass::registerClassLoadMethod(
    "CompositeModel", typeid(CompositeModel), CompositeModelHelper::load);

CClassConstSP const CompositeModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("CompositeModel::IIntoProduct",
                                    typeid(CompositeModel::IIntoProduct), 
                                    CompositeModelHelper::IntoProduct_load);


DRLIB_END_NAMESPACE
