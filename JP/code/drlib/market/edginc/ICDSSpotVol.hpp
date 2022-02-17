//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CDSSpotVol.hpp
//
//   Description : Interface defining what CDS Spot Vols can do (here 'spot'
//                 is in the C&R sense - the vol used in the diffusion equation
//                 and not the Black vols)
//
//   Author      : Mark A Robson
//
//   Date        : November 26, 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ICDSSPOTVOL_HPP
#define EDR_ICDSSPOTVOL_HPP

#include "edginc/ICDSVol.hpp"

DRLIB_BEGIN_NAMESPACE

/** Essentially a marker interface for CDS 'Spot' Vols can do (here 'spot' is in
    the C&R sense - the vol used in the diffusion equation
    and not the Black vols) */
class MARKET_DLL ICDSSpotVol: public virtual ICDSVol{
public:
    static CClassConstSP const TYPE; // in CDSParSpreads.cpp

    ~ICDSSpotVol();
private:
    static void load(CClassSP& clazz);  // in CDSParSpreads.cpp
};


DRLIB_END_NAMESPACE

#endif




