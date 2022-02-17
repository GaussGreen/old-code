#include "edginc/config.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

static void IConditionalDefaultsModelLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IConditionalDefaultsModel, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IConditionalDefaultsModel::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IConditionalDefaultsModel", 
        typeid(IConditionalDefaultsModel), 
        IConditionalDefaultsModelLoad);

DEFINE_TEMPLATE_TYPE(IConditionalDefaultsModelArray);

/* external symbol to allow class to be forced to be linked in */
bool IConditionalDefaultsModelLoad(){
    return (IConditionalDefaultsModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
