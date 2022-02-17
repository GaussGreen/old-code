//----------------------------------------------------------------------------
//
//   Group       : GCCG DR
//
//   Filename    : ClosedFormEnergy.cpp
//
//   Description : Input info for Energy Closed form model
//
//   Author      : Sean Chen
//
//   Date        : August 25, 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormEnergy.hpp"
#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE

static void purifyFix(ClosedFormEnergy::IProduct*     product,
                      ClosedFormEnergy*               closedForm,
                      CControl*                 control,
                      CResults*                 results){
    product->price(closedForm, control, results);
}

ClosedFormEnergy::ClosedFormEnergy(): Model(TYPE) {}

/** calculate single price and store result in CResult */
void ClosedFormEnergy::Price(CInstrument*  instrument, 
                         CControl*     control, 
                         CResults*     results){
    static const string method = "ClosedFormEnergy::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormEnergy::IntoProduct");
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

IModel::WantsRiskMapping ClosedFormEnergy::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class ClosedFormEnergyHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormEnergy, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultClosedFormEnergy);

        FIELD(smileOff, "Not to use smile vols");
        FIELD_MAKE_OPTIONAL(smileOff);
    }

    static IObject* defaultClosedFormEnergy(){
        return new ClosedFormEnergy();
    }
};

CClassConstSP const ClosedFormEnergy::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormEnergy", typeid(ClosedFormEnergy), ClosedFormEnergyHelper::load);

typedef ClosedFormEnergy::IIntoProduct ClosedFormEnergy_IIntoProduct;
CClassConstSP const ClosedFormEnergy_IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormEnergy::IIntoProduct",
                                    typeid(ClosedFormEnergy_IIntoProduct), 0);

bool  ClosedFormEnergyLoad() {
    return ClosedFormEnergy::TYPE != 0 &&
           ClosedFormEnergy::IIntoProduct::TYPE != 0;
}

DRLIB_END_NAMESPACE
