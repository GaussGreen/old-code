//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CFDGridPass.cpp
//
//   Description : FD associate a product to a model
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CFDGridPass.hpp"
#include "edginc/Instrument.hpp"


DRLIB_BEGIN_NAMESPACE

CFDGridPass::CFDGridPass(): CModel(TYPE){}


/** calculate single price and store result in CResult */
void CFDGridPass::Price(CInstrument*  instrument, 
                          CControl*     control, 
                          CResults*     results){
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("CFDGridPass::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support CFDGridPass::IntoProduct");
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
        throw ModelException(&e, "CFDGridPass::Price");
    }
    delete product;
}

IModel::WantsRiskMapping CFDGridPass::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class CFDGridPassHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CFDGridPass, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultCFDGridPass);
		FIELD(volType, "Type of vol to use");
		FIELD(iMax, "number of steps for S");
		FIELD(nMax, "number of steps for T");
		FIELD(divType, "0 - as yield, 1 - dollar, 2 - mix");
		FIELD_MAKE_OPTIONAL(divType);

    }

    static IObject* defaultCFDGridPass(){
        return new CFDGridPass();
    }
};

CClassConstSP const CFDGridPass::TYPE = CClass::registerClassLoadMethod(
    "FDGrid", typeid(CFDGridPass), CFDGridPassHelper::load);
bool   CFDGridPassLoad () {
    return (CFDGridPass::TYPE != 0);
   }


CClassConstSP const CFDGridPass::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("CFDGridPass::IIntoProduct",
                                    typeid(CFDGridPass::IIntoProduct), 0);


DRLIB_END_NAMESPACE
