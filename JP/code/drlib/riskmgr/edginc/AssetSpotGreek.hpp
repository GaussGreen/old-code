//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetSpotGreek.hpp
//
//   Description : Marker interface for sensitivities that are calculated 
//                 solely via altering the spot of an asset (eg Delta, 
//                 FXDelta, DeltaProxy, CrossGamma etc)
//
//   Author      : Mark A Robson
//
//   Date        : 11 July 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASSETSPOTGREEK_HPP
#define EDR_ASSETSPOTGREEK_HPP


DRLIB_BEGIN_NAMESPACE

/** Marker interface for sensitivities that are calculated solely via
    altering the spot of an asset (eg Delta, FXDelta, DeltaProxy,
    CrossGamma etc) */
class RISKMGR_DLL IAssetSpotGreek {
public:
    virtual ~IAssetSpotGreek(); // empty

protected:
    IAssetSpotGreek(); // empty
};


DRLIB_END_NAMESPACE
#endif
