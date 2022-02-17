#include "edginc/config.hpp"
#include "edginc/IFuncBasisWrapper.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IFuncBasisWrapper::TYPE = CClass::registerInterfaceLoadMethod(
    "IFuncBasisWrapper", typeid(IFuncBasisWrapper), IFuncBasisWrapper::load);

void IFuncBasisWrapper::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IFuncBasisWrapper, clazz);
    EXTENDS(IObject);
}

bool IFuncBasisWrapperLoad() {
    return IFuncBasisWrapper::TYPE != 0;
}
DRLIB_END_NAMESPACE