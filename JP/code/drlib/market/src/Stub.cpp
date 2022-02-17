//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Stub.cpp
//
//   Description : Stub payment interface
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

Stub::~Stub() {
    // empty
}

Stub::Stub(CClassConstSP clazz): CObject(clazz){}

/** Invoked when Class is 'loaded' */
static void loadStub(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Stub, clazz);
    SUPERCLASS(CObject);
    // no fields
}

CClassConstSP const Stub::TYPE = CClass::registerClassLoadMethod(
    "Stub", typeid(Stub), loadStub);

DRLIB_END_NAMESPACE



