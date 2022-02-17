//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolExplicitTier2.hpp
//
//   Description : Model volatility (not market vol) class for energy spreads
//
//   Author      : Lawrence Siu
//
//   Date        : May 01, 2006
//
//----------------------------------------------------------------------------
#ifndef ENERGY_INST_VOL_SPREAD_HPP
#define ENERGY_INST_VOL_SPREAD_HPP

#include "edginc/EnergyInstVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

class EnergyInstVolExplicitTier2 : public EnergyInstVolBase
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyInstVolExplicitTier2Helper;

    ~EnergyInstVolExplicitTier2();

    void validatePop2Object();

    virtual int calibrate() { return 0; } // not doing anything

private:
    
    EnergyInstVolExplicitTier2();

    //int deriveRatios();

    DoubleArray sigmas;

};

typedef smartConstPtr<EnergyInstVolExplicitTier2> EnergyInstVolExplicitTier2ConstSP;
typedef smartPtr<EnergyInstVolExplicitTier2> EnergyInstVolExplicitTier2SP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolExplicitTier2> EnergyInstVolExplicitTier2Wrapper;
    
DRLIB_END_NAMESPACE

#endif
