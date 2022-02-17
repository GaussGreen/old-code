//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFRequest.cpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFRequest.hpp"

DRLIB_BEGIN_NAMESPACE

static void myLoad(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(PDFRequest, clazz);
    SUPERCLASS(CObject);
}


PDFRequest::~PDFRequest(){}

PDFRequest::PDFRequest(const CClassConstSP& clazz):CObject(clazz) {}


CClassConstSP const PDFRequest::TYPE = CClass::registerClassLoadMethod(
    "PDFRequest", typeid(PDFRequest), myLoad);


DRLIB_END_NAMESPACE

