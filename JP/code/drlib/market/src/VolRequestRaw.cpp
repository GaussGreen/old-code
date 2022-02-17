//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestRaw.cpp
//
//   Description : Vol request for accessing the 'raw' vol for cases where
//                 you know the type of vol that you are expecting and the
//                 vol implements IVolProcessed
//
//   Author      : Mark A Robson
//
//   Date        : 16 May 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

VolRequestRaw::~VolRequestRaw(){}

VolRequestRaw::VolRequestRaw(): CVolRequest(TYPE){}

static CVolBaseSP copyVol(CVolProcessedSP volProc){
    static const string method = "VolRequestRaw::copyVol";
    if (!CVolBase::TYPE->isInstance(volProc)){
        throw ModelException(method, "Incorrect type returned for "
                             "VolRequestRaw");
    }
    CVolBaseSP tempVol(CVolBaseSP::dynamicCast(volProc));
    return CVolBaseSP(tempVol.clone());
}

CVolBaseSP VolRequestRaw::copyVolBase(const IMultiFactors& mAsset,
                                      int                  iAsset){
    static const string method = "VolRequestRaw::copyVolBase";
    VolRequestRaw volRequest;
    CVolProcessedSP volProc(mAsset.factorGetProcessedVol(iAsset, &volRequest));
    return copyVol(volProc);
}

CVolBaseSP VolRequestRaw::copyVolBase(const Asset& asset){
    VolRequestRaw volRequest;
    CVolProcessedSP volProc(asset.getProcessedVol(&volRequest));
    return copyVol(volProc);
}

void VolRequestRaw::load(CClassSP& clazz){
    REGISTER(VolRequestRaw, clazz);
    SUPERCLASS(CVolRequest);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* VolRequestRaw::defaultConstructor(){
    return new VolRequestRaw();
}

CClassConstSP const VolRequestRaw::TYPE = 
CClass::registerClassLoadMethod("VolRequestRaw", typeid(VolRequestRaw), load);

DRLIB_END_NAMESPACE

