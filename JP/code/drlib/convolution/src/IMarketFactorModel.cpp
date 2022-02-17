#include "edginc/config.hpp"
#include "edginc/IMarketFactorModel.hpp"

DRLIB_BEGIN_NAMESPACE

static void IMarketFactorModelLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IMarketFactorModel, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IMarketFactorModel::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IMarketFactorModel", 
        typeid(IMarketFactorModel), 
        IMarketFactorModelLoad);

DEFINE_TEMPLATE_TYPE(IMarketFactorModelArray);

/* external symbol to allow class to be forced to be linked in */
bool IMarketFactorModelLoad(){
    return (IMarketFactorModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
