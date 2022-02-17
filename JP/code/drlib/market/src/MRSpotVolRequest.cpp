//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : MRSpotVolRequest.cpp
//
//   Description : A way of asking for 'spot' vols. 
//                 Here 'spot' means in the traditional IR sense
//
//   Author      : Mark A Robson
//
//   Date        : 14 Dec 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

//// initially this is going to be used for a flat cds vol curve so we
//// can defer for now what interpolation data we need
MRSpotVolRequest::MRSpotVolRequest(): CVolRequest(TYPE){}

IObject* MRSpotVolRequest::defaultConstructor(){
    return new MRSpotVolRequest();
}

CClassConstSP const MRSpotVolRequest::TYPE = CClass::registerClassLoadMethod(
    "MRSpotVolRequest", typeid(MRSpotVolRequest), load);

void MRSpotVolRequest::load(CClassSP& clazz){
    // class deliberately private - no need for clients to know it exists
    REGISTER(MRSpotVolRequest, clazz);
    SUPERCLASS(CVolRequest);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

MRSpotVolProcessed::~MRSpotVolProcessed(){}

MRSpotVolProcessed::MRSpotVolProcessed(CClassConstSP clazz):
    CObject(TYPE){}

void MRSpotVolProcessed::load(CClassSP& clazz){
    // class deliberately private - no need for clients to know it exists
    REGISTER(MRSpotVolProcessed, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
}    

CClassConstSP const MRSpotVolProcessed::TYPE = CClass::registerClassLoadMethod(
    "MRSpotVolProcessed", typeid(MRSpotVolProcessed), load);

DRLIB_END_NAMESPACE
