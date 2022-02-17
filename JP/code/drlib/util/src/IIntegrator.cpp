#include "edginc/config.hpp"
#include "edginc/IIntegrator.hpp"

DRLIB_BEGIN_NAMESPACE

static void IIntegratorLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IIntegrator, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IIntegrator::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IIntegrator", 
        typeid(IIntegrator), 
        IIntegratorLoad);

DEFINE_TEMPLATE_TYPE(IIntegratorArray);

/* external symbol to allow class to be forced to be linked in */
bool IIntegratorLoad(){
    return (IIntegrator::TYPE != 0);
}

DRLIB_END_NAMESPACE
