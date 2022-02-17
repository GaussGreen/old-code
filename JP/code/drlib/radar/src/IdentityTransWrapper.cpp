#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/IdentityTransWrapper.hpp"
#include "edginc/IFittingVarTransform.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IFittingVarTransformSP IdentityTransWrapper::getTransform(void) const {
    return IFittingVarTransformSP(new IdentityTransform());
}

void IdentityTransWrapper::validatePop2Object() {
}

void IdentityTransWrapper::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("PolyBasis object");
    REGISTER(IdentityTransWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IFittingVarTransWrapper);

    EMPTY_SHELL_METHOD(defaultConstructor);
    Addin::registerConstructor(Addin::UTILITIES, IdentityTransWrapper::TYPE);
}

CClassConstSP const IdentityTransWrapper::TYPE = CClass::registerClassLoadMethod(
    "IdentityTransform", typeid(IdentityTransWrapper), IdentityTransWrapper::load);

/******************************/
// for type linking
bool IdentityTransWrapperLoad(void){
    return (IdentityTransWrapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
