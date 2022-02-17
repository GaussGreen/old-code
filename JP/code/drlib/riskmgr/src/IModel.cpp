//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IModel.cpp
//
//   Description : Interface class for Models
//
//   Author      : Jose Hilera
//
//   Date        : 28 June 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IModel.hpp"


DRLIB_BEGIN_NAMESPACE

IModel::IModel(){}
IModel::~IModel() {}

/** Invoked when Class is 'loaded' */
void IModel::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IModel, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IModel::TYPE = CClass::registerInterfaceLoadMethod(
    "IModel", typeid(IModel), load);


DRLIB_END_NAMESPACE
