//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolRegular.hpp
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
#ifndef ENERGY_INST_VOL_REGULAR_HPP
#define ENERGY_INST_VOL_REGULAR_HPP

#include "edginc/EnergyInstVolCalibrated.hpp"
//#include "edginc/EnergyAlpha.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyInstVolRegular: public EnergyInstVolCalibrated
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyInstVolRegularHelper;
    friend class EnergyInstVolRegularWrappper;

    virtual ~EnergyInstVolRegular();

    virtual void validatePop2Object();

    virtual int deriveSigmas(const int numBenchmarks,
			     const DoubleArray& vols);

	DoubleMatrix calibSigmas(const DateTimeArray & futureMaturities,
							 const DateTimeArray & optionExpiries,
							 const DoubleArray & atmVols,
							 double x,
							 double z) 
							 const;

	DoubleMatrix deriveSigmas(const DateTimeArray & fSigmasDates, // future maturities
							  const DateTimeArray & fOptionExpiries, // option maturities
							  const DoubleArray & vols, // atm vols
							  double x,
							  double z) 
							  const;

	double calibX() const;
	double calibY() const {return 1.0; }
	double calibZ(double x) const;

private:

    EnergyInstVolRegular();

};

typedef smartConstPtr<EnergyInstVolRegular> EnergyInstVolRegularConstSP;
typedef smartPtr<EnergyInstVolRegular> EnergyInstVolRegularSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolRegular> EnergyInstVolRegularWrapper;
    
DRLIB_END_NAMESPACE

#endif
