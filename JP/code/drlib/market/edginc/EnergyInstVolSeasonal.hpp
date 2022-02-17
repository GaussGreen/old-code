//----------------------------------------------------------------------------
//
//   Group    : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolSeasonal.hpp
//
//   Description : Energy 2 factor vol surface. Based on drcommodityvolsurfacei.cpp
//              and drcommodityInstaneousvolsurface.cpp in FXLIB.
//           
//
//   Author      : Sean Chen
//
//   Date     : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_INST_VOL_SEASONAL_HPP
#define ENERGY_INST_VOL_SEASONAL_HPP

#include "edginc/EnergyInstVolCalibrated.hpp"
//#include "edginc/EnergyAlpha.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyInstVolSeasonal: public EnergyInstVolCalibrated
{

public:    

    static CClassConstSP const TYPE;

    enum ECalibrationMethod {kSigma1Bar_X, kSeasonal};

    friend class EnergyInstVolSeasonalHelper;
    friend class EnergyInstVolSeasonalWrappper;

    virtual ~EnergyInstVolSeasonal();

    virtual void validatePop2Object();

    virtual int deriveSigmas(const int numBenchmarks,
                             const DoubleArray& vols);

    virtual void calculateLongTermMarketVol(double& marketVol, 
                                     int& contractID,
                                     int& numCalibratableBenchmark) const;

    virtual void calculateLongTermModelVol(double& modelVol,  
                                     int& contractID, 
                                     int& numCalibratableBenchmark) const;

private:

    EnergyInstVolSeasonal();

    // Calibration method: kSigma1Bar_X, compute X(T) from sigma1Bar
    int deriveSigmasSigma1Bar_X(const int numBenchmarks,
                                const DoubleArray& vols);

    int deriveSigmasFromSigma1BarGuess_X(const int numBenchmarks,
                                         const DoubleArray& vols,
                                         const int groupStart,
                                         const int groupSize,
                                         const double sigma1BarGuess);

    // Calibration method: kSeasonal, compute sigmas from X(T) and y(t)
    int deriveSigmasSeasonal(const int numBenchmarks,
                             const DoubleArray& vols);

    // Works when y(t) is given
    int calculateMaxCalibratableSigma1Bar(DoubleArray& maxCalibratableSigma1Bar,
                                          const int numBenchmarks,
                                          const DoubleArray& vols,
                                          const int groupStart,
                                          const int groupSize) const;

    double getInstVolRatio(const double alpha,
                           const double beta,
                           const double normV1_1,
                           const double normV1_2,
                           const double normV2Bar,
                           const double normV2,
                           const double alphaFactor1,
                           const double betaFactor1,
                           const double alphaFactor2,
                           const double betaFactor2) const;


     double          beta;
     DoubleArray     fNSigma1s; //x(T)
     DoubleArray     fNSigma2s; //y(t)
     double          fNSigma2Bar; //z
     double          fNSigma1Bar; //sigma1bar
     double          fSigma1Bar; //calculated sigma1bar

     ECalibrationMethod fCalibrationMethod; // $unregistered

    
};

typedef smartConstPtr<EnergyInstVolSeasonal> EnergyInstVolSeasonalConstSP;
typedef smartPtr<EnergyInstVolSeasonal> EnergyInstVolSeasonalSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolSeasonal> EnergyInstVolSeasonalWrapper;
    
DRLIB_END_NAMESPACE

#endif
