#include "edginc/config.hpp"
#include "edginc/IConvolutor.hpp"
DRLIB_BEGIN_NAMESPACE

static void IConvolutorLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IConvolutor, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IConvolutor::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IConvolutor", 
        typeid(IConvolutor), 
        IConvolutorLoad);

DEFINE_TEMPLATE_TYPE(IConvolutorArray);

/* external symbol to allow class to be forced to be linked in */
bool IConvolutorLoad(){
    return (IConvolutor::TYPE != 0);
}

DRLIB_END_NAMESPACE
