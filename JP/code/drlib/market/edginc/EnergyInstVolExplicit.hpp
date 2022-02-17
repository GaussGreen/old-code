//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolExplicit.hpp
//
//   Description : Energy 2 factor vol surface. Based on drcommodityvolsurfacei.cpp
//                    and drcommodityInstaneousvolsurface.cpp in FXLIB.
//                 
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_INST_VOL_EXPLICIT_HPP
#define ENERGY_INST_VOL_EXPLICIT_HPP

#include "edginc/EnergyInstVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyInstVolExplicit : public EnergyInstVolBase
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyInstVolExplicitHelper;

    friend class GetEnergyATMInstVolExplicitAddin;

    ~EnergyInstVolExplicit();

    void validatePop2Object();

    int calibrate();

private:
    
    EnergyInstVolExplicit();

    int deriveRatios();

    DoubleArray sigma1s;
    DoubleArray sigma2s;
    DoubleArray sigma1Bars;
    DoubleArray sigma2Bars;
};

typedef smartConstPtr<EnergyInstVolExplicit> EnergyInstVolExplicitConstSP;
typedef smartPtr<EnergyInstVolExplicit> EnergyInstVolExplicitSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolExplicit> EnergyInstVolExplicitWrapper;
    
DRLIB_END_NAMESPACE

#endif
