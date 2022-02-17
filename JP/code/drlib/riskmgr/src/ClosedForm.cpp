//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedForm.cpp
//
//   Description : Closed Form Algorithm (asks instrument to do it)
//
//   Author      : Andrew J Swain
//
//   Date        : 16 February 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE
ClosedForm::ClosedForm(): Model(TYPE){}

static void purifyFix(ClosedForm::IProduct*     product,
                      ClosedForm*               closedForm,
                      CControl*                 control,
                      CResults*                 results){
    product->price(closedForm, control, results);
}

/** calculate single price and store result in CResult */
void ClosedForm::Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results){
    static const string method = "ClosedForm::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedForm::IntoProduct");
    }
    IProduct*     product = 0;
    try{
        product = intoProd->createProduct(this);
        //product->price(this, control, results);
        purifyFix(product, this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

ClosedForm::IProduct::~IProduct()
{}

IModel::WantsRiskMapping ClosedForm::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class ClosedFormHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedForm, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultClosedForm);
        // no fields
    }
    // for ClosedForm::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(ClosedForm::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultClosedForm(){
        return new ClosedForm();
    }
};

CClassConstSP const ClosedForm::TYPE = CClass::registerClassLoadMethod(
    "ClosedForm", typeid(ClosedForm), ClosedFormHelper::load);

CClassConstSP const ClosedForm::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedForm::IIntoProduct",
                                    typeid(ClosedForm::IIntoProduct),
                                    ClosedFormHelper::loadIntoProduct);


DRLIB_END_NAMESPACE
