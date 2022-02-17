//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedFormLN.hpp
//
//   Description : Closed Form Algorithm (asks instrument to do it)
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"

DRLIB_BEGIN_NAMESPACE

CClosedFormLN::IIntoProduct::IIntoProduct(){}
CClosedFormLN::IIntoProduct::~IIntoProduct(){}
CClosedFormLN::IProduct::IProduct(){}
CClosedFormLN::IProduct::~IProduct(){}


/** Constructor takes type of vol to use */
CClosedFormLN::CClosedFormLN(const string& volType,
                             bool          allowNegativeFwdVar):
CModelLN(TYPE, volType, allowNegativeFwdVar){}

/** calculate single price and store result in CResult */
void CClosedFormLN::Price(CInstrument*  instrument, 
                          CControl*     control, 
                          CResults*     results){
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("CClosedFormLN::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support CClosedFormLN::IntoProduct");
    }
    IProduct*  product = 0;
    try{
        // cast to CClosedFormLN::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        product = intoProd.createProduct(this);
        // and the invoke the pricing
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(&e, "CClosedFormLN::Price");
    }
    delete product;
}

/** Essentially relies on instrument implementing ISensitiveIRVolPoints.
    If not returns null. */
IRGridPointAbsArraySP CClosedFormLN::getSensitiveIRVolPoints(
    OutputNameConstSP outputName,
    const Instrument* inst) const{
    const ISensitiveIRVolPoints* vpImnt =
        dynamic_cast<const ISensitiveIRVolPoints*>(inst);
    return vpImnt? vpImnt->getSensitiveIRVolPoints(outputName, this): 
        IRGridPointAbsArraySP();
}

CClosedFormLN::CClosedFormLN():
CModelLN(TYPE){};

class CClosedFormLNHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CClosedFormLN, clazz);
        SUPERCLASS(CModelLN);
        EMPTY_SHELL_METHOD(defaultCClosedFormLN);
        // no fields
    }

    // for ClosedFormLN::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(CClosedFormLN::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultCClosedFormLN(){
        return new CClosedFormLN();
    }
};

CClassConstSP const CClosedFormLN::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormLN", typeid(CClosedFormLN), CClosedFormLNHelper::load);

CClassConstSP const CClosedFormLN::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormLN::IIntoProduct",
                                    typeid(CClosedFormLN::IIntoProduct), 
                                    CClosedFormLNHelper::loadIntoProduct);

bool  CClosedFormLNLoad() {
    return (CClosedFormLN::TYPE != 0);
   }



DRLIB_END_NAMESPACE
