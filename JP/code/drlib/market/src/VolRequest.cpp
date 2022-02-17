//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequest.cpp
//
//   Description : Abstract vol request interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLREQUEST_CPP
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

CVolRequest::~CVolRequest(){}

CVolRequest::CVolRequest(const CClassConstSP& clazz): CObject(clazz){}

static void myLoad(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CVolRequest, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const CVolRequest::TYPE = CClass::registerClassLoadMethod(
    "VolRequest", typeid(CVolRequest), myLoad);

// initialise type for array of CVolRequest
DEFINE_TEMPLATE_TYPE_WITH_NAME("VolRequestArray", CVolRequestArray);

DRLIB_END_NAMESPACE


