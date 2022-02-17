//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Description : generic indicator function.
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IIndicCreator.hpp"

DRLIB_BEGIN_NAMESPACE

void IIndicCreator::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IIndicCreator, clazz);
    EXTENDS(CModel::IProdCreator);
}

CClassConstSP const IIndicCreator::TYPE = 
    CClass::registerInterfaceLoadMethod("IIndicCreator", typeid(IIndicCreator), load);
DEFINE_TEMPLATE_TYPE(IIndicCreatorArray);

DRLIB_END_NAMESPACE

