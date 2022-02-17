//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TwoSidedDeriv.cpp
//
//   Description : Marker interface for sensitivities that are calculated via
//                 altering the spot of an asset (eg Delta, FXDelta, 
//                 DeltaProxy, CrossGamma etc)
//
//   Author      : Mark A Robson
//
//   Date        : 11 July 2002
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

/** Marker interface for sensitivities that are calculated via
    altering the spot of an asset (eg Delta, FXDelta, DeltaProxy,
    CrossGamma etc) */

IAssetSpotGreek::~IAssetSpotGreek(){}
IAssetSpotGreek::IAssetSpotGreek() {}

DRLIB_END_NAMESPACE
