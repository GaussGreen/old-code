//----------------------------------------------------------------------------
//
//   Filename    : SamplingConvention.cpp
//
//   Description : Classes for handling holiday adjustment 
//                 for use with MarketObservable and ObservableHistory
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares   
//
//   Date        : February 2 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SamplingConvention.hpp"

DRLIB_BEGIN_NAMESPACE
    
CClassConstSP const SamplingConvention::TYPE = CClass::registerInterfaceLoadMethod(
    "SamplingConvention", typeid(SamplingConvention), load);

void SamplingConvention::load(CClassSP& clazz){
    REGISTER_INTERFACE(SamplingConvention, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

// Returns false if obs is to be omitted
bool UnadjustedConvention::observationDate(const DateTime& sampleDate, 
                                           const Holiday*  hols,
                                           DateTime*       obsDate) const {
    *obsDate = sampleDate;
    return true;
}

bool UnadjustedConvention::rollAssetsTogether() const {
    return false;
}

bool UnadjustedConvention::isOmit() const {
    return false;
}

bool UnadjustedConvention::scheduleMoves() const {
    return false;
}

bool UnadjustedConvention::isUnadjusted() const {
    return true;
}

UnadjustedConvention::UnadjustedConvention() : CObject(TYPE) {}

void UnadjustedConvention::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(UnadjustedConvention, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(SamplingConvention);
    EMPTY_SHELL_METHOD(defaultUnadjustedConvention);
}

IObject* UnadjustedConvention::defaultUnadjustedConvention() {
    return new UnadjustedConvention();
}

CClassConstSP const UnadjustedConvention::TYPE = CClass::registerClassLoadMethod(
    "UnadjustedConvention", typeid(UnadjustedConvention), load);

DRLIB_END_NAMESPACE





