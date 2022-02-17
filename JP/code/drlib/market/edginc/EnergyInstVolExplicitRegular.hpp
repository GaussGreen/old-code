//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolExplicitRegular.hpp
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
#ifndef ENERGY_INST_VOL_EXPLICIT_REGULAR_HPP
#define ENERGY_INST_VOL_EXPLICIT_REGULAR_HPP

#include "edginc/EnergyInstVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyInstVolExplicitRegular : public EnergyInstVolBase
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyInstVolExplicitRegularHelper;

    ~EnergyInstVolExplicitRegular();

    void validatePop2Object();

    int calibrate() {return 0;}

    double getX() const {return x;}
    double getY() const {return y;}
    double getZ() const {return z;}

private:
    
    EnergyInstVolExplicitRegular();

    //int deriveRatios();

    //DoubleArray sigma1s;
    //DoubleArray sigma2s;
    //DoubleArray sigma1Bars;
    //DoubleArray sigma2Bars;

    // QSRM Oil model calibrated params:
    double x, y, z;

    // Sampras Oil model calibrated params:
    double w1, w2;
    double rho;
    DoubleArray sigma1s;

    // The following functions are used for converting from Sampras Oil
    // model calibrated params to QSRM Oil model params:
    void convertInputs();
    double computeZ(); 
    double computeX();
};

typedef smartConstPtr<EnergyInstVolExplicitRegular> EnergyInstVolExplicitRegularConstSP;
typedef smartPtr<EnergyInstVolExplicitRegular> EnergyInstVolExplicitRegularSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolExplicitRegular> EnergyInstVolExplicitRegularWrapper;
    
DRLIB_END_NAMESPACE

#endif
