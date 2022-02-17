//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : VolTypeSensitiveStrikes.cpp
//
//   Description : For Vega Matrix we call getSensitiveStrikes on the asset
//                 In some instances (e.g. ProxyVol) we also need to call
//                 the actual vol for the sensitive strikes
//
//   Author      : Ian S Stares
//
//   Date        : 09 October 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/VolTypeSensitiveStrikes.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IVolTypeSensitiveStrikes::TYPE = 
CClass::registerInterfaceLoadMethod("IVolTypeSensitiveStrikes",
                                    typeid(IVolTypeSensitiveStrikes), 0);

IVolTypeSensitiveStrikes::IVolTypeSensitiveStrikes(){}

IVolTypeSensitiveStrikes::~IVolTypeSensitiveStrikes(){}

// forced linking
bool VolTypeSensitiveStrikesLoad() {
    return (IVolTypeSensitiveStrikes::TYPE != NULL);
}

DRLIB_END_NAMESPACE
