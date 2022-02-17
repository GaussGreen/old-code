//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestDVF.cpp
//
//   Description : Deterministic Vol Func request
//
//   Author      : JNJ
//
//   Date        : 06 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestDVF.hpp"

DRLIB_BEGIN_NAMESPACE

CVolRequestDVF::~CVolRequestDVF(){}

CVolRequestDVF::CVolRequestDVF(const CClassConstSP& clazz): CVolRequest(clazz){}

static void CVolRequestDVFLoad(CClassSP& clazz){
    REGISTER(CVolRequestDVF, clazz);
    SUPERCLASS(CVolRequest);
}

CClassConstSP const CVolRequestDVF::TYPE = CClass::registerClassLoadMethod(
    "VolRequestDVF", typeid(CVolRequestDVF), CVolRequestDVFLoad);

DRLIB_END_NAMESPACE


