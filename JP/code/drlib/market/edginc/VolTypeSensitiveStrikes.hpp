//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : VolTypeSensitiveStrikes.hpp
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

#ifndef VOL_TYPE_SENSITIVE_STRIKES_HPP
#define VOL_TYPE_SENSITIVE_STRIKES_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Asset.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IVolTypeSensitiveStrikes: virtual public IObject{
public:
    IVolTypeSensitiveStrikes();
    virtual ~IVolTypeSensitiveStrikes();

    /** returns sensitive strikes for a given vol request and asset */
    virtual void getSensitiveStrikes(const CAsset* asset,
                             const CVolRequest* volRequest,
                             OutputNameConstSP outputName,
                             const SensitiveStrikeDescriptor& sensStrikeDesc,
                             DoubleArraySP sensitiveStrikes) const = 0;

    static CClassConstSP const TYPE;
};

DRLIB_END_NAMESPACE
#endif
