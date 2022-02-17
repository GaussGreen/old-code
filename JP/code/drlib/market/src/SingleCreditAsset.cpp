//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Empty class that just derives from CreditAsset. See 
//                 CreditAsset.hpp for details
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SingleCreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE

/** Destructor */
SingleCreditAsset::~SingleCreditAsset() 
{}

SingleCreditAsset::SingleCreditAsset() : CreditAsset(TYPE)
{}

IObject* SingleCreditAsset::defaultConstructor(){
    return new SingleCreditAsset();
}

void SingleCreditAsset::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SingleCreditAsset, clazz);
    SUPERCLASS(CreditAsset);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const SingleCreditAsset::TYPE = CClass::registerClassLoadMethod(
    "SingleCreditAsset", typeid(SingleCreditAsset), load);

DRLIB_END_NAMESPACE

