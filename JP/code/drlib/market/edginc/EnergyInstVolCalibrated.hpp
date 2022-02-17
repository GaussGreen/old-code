//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolCalibrated.hpp
//
//   Description : Energy 2 factor vol surface. Based on drcommodityvolsurfacei.cpp
//             and drcommodityInstaneousvolsurface.cpp in FXLIB.
//          
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_INST_VOL_CALIBRATED_HPP
#define ENERGY_INST_VOL_CALIBRATED_HPP

#include "edginc/EnergyInstVolBase.hpp"
//#include "edginc/EnergyAlpha.hpp"


DRLIB_BEGIN_NAMESPACE

class EnergyImpliedVolSurface;

class MARKET_DLL EnergyInstVolCalibrated: public EnergyInstVolBase          
{

public:    

    static CClassConstSP const TYPE;

    friend class EnergyInstVolCalibratedHelper;
    friend class EnergyInstVolCalibratedWrappper;

    virtual ~EnergyInstVolCalibrated();

    virtual void validatePop2Object();

    double getVolRatio() const;
    double getInstCorr() const;
    int getFirstContract() const;
    int getSecondContract() const;
	const DoubleArray & getVols() const;

    virtual int getTweakedSigmas(int numContracts,
                                 int contractTweak,
                                 double tweakSize);

    virtual void reconfigure(const EnergyImpliedVolSurface& volSurface,
                             double delta);

    virtual int calibrate();

    virtual int deriveSigmas(const int numBenchmarks, 
                             const DoubleArray& vols);

protected:    

    EnergyInstVolCalibrated(const CClassConstSP&);

    DoubleArray     fVols;
    int             fFirstContract; // $unregistered
    int             fSecondContract; // $unregistered
    double          fInstVolRatio; // $unregistered
    double          fInstCorr; // $unregistered

    // utility function
    double getNormalizedSigma1(const double alpha,
                               const int firstContract, 
                               const int secondContract,
                               const double InstVolRatio,
                               const double InstCorr) const;
    double getNormalizedSigma2bar(const double alpha,
                                  const double NormSigma1,
                                  const int firstContract, 
                                  const int secondContract,
                                  const double InstVolRatio) const;
};

typedef smartConstPtr<EnergyInstVolCalibrated> EnergyInstVolCalibratedConstSP;
typedef smartPtr<EnergyInstVolCalibrated> EnergyInstVolCalibratedSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolCalibrated> EnergyInstVolCalibratedWrapper;
    
DRLIB_END_NAMESPACE

#endif
