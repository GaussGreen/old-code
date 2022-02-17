#include "edginc/config.hpp"
#include "edginc/IFittingVarTransWrapper.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IFittingVarTransWrapper::TYPE = CClass::registerInterfaceLoadMethod(
    "IFittingVarTransWrapper", typeid(IFittingVarTransWrapper), IFittingVarTransWrapper::load);

void IFittingVarTransWrapper::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IFittingVarTransWrapper, clazz);
    EXTENDS(IObject);
}

bool IFittingVarTransWrapperLoad() {
    return IFittingVarTransWrapper::TYPE != 0;
}
DRLIB_END_NAMESPACE
