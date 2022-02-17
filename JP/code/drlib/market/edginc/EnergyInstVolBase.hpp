//----------------------------------------------------------------------------
//
//   Group    : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolBase.hpp
//
//   Description : Energy 2 factor vol model. Based on drcommodityvolmodel.cpp
//              and drcommodityvolmodeltwofactor....cpp in FXLIB.
//           
//
//   Author      : Sean Chen
//
//   Date     : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_INST_VOL_BASE_HPP
#define ENERGY_INST_VOL_BASE_HPP

#include "edginc/EnergyVolBase.hpp"
#include "edginc/DoubleMatrix.hpp"
//#include "edginc/EnergyAlpha.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyInstVolBase: public EnergyVolBase               
{

public:    

     static CClassConstSP const TYPE;

     friend class EnergyInstVolBaseHelper;
     friend class EnergyInstVolBaseWrappper;
 
     virtual ~EnergyInstVolBase();
 
    virtual void validatePop2Object();
     
     // These two functions are for capital T model only
     virtual void calculateLongTermMarketVol(double& marketVol, int& contractID, int& numCalibratableBenchmark) const;
     virtual void calculateLongTermModelVol(double& marketVol, int& contractID, int& numCalibratableBenchmark) const;

     // These four implement parents'
     virtual double getATMVol(const DateTime& expiry, 
                              const DateTime& futuresExpiry = DateTime()) const;

    virtual double getSmileVolByStrike(const DateTime& expiry, double strike) const;
    virtual double getSmileVolByDelta(const DateTime& expiry, double delta) const;

     // Major one.... that all children use
     virtual void analyticApproximation(
                      double* fwdPositive,
                      double* fwdNegative,
                      double* volPositive,
                      double* volNegative,
                      double* correlation,
                      const DateTime& expiry,
                      const DateTimeArray& fixingDates,
                      const DoubleArray& weights,
                      const DateTimeArray& maturities,
                      const DoubleArray& prices) const;

     const DoubleArray& getAlphas() const; 
     const DateTime& getMaxPricingDate() const;
     const DateTimeArray& getFutureMaturityDates() const;
     const DateTimeArray& getOptionExpiryDates() const;
     const DoubleArray& getRatios() const;
     const DoubleMatrix& getSigmas() const;

     virtual int calibrate() = 0;

protected:

     EnergyInstVolBase();
     EnergyInstVolBase(const CClassConstSP& clazz);
     
     DoubleArray         fAlphas;
     DateTime            fMaxPricingDate;
     DoubleArray         fRatios;
     DateTimeArray       fFutureMaturityDates; //fSigmasDates_;
     DateTimeArray       fOptionExpiryDates;
     DoubleMatrix        fSigmas;     

     //All inst vols models have alpha
     double  alpha; 

     // utility function
     virtual double getInstVolRatio(    
         const double alpha,
         const double beta,
         const double v1Bar,
         const double v1,
         const double v2Bar,
         const double v2,
         const double alphaFactor1,
         const double betaFactor1,
         const double alphaFactor2,
         const double betaFactor2) const;
     
private:

     void AnalyticVariance(
         double* fwd,
         double* var,
         const DoubleArray& x_T,
         const DoubleArray& fixingTimes,
         const DoubleArray& fixingPrices,
         const DoubleArray& fixingVols,
         const DoubleArray& expAlpha_T,
         const DoubleArray& expBeta_T,
         const DoubleArray& expAlpha_t,
         const DoubleArray& expBeta_t,
         const IntArray& intervalEndpoints,
         const DoubleArray& boundaryTimes,
         const DoubleArray& expAlpha_S,
         const DoubleArray& expBeta_S) const;

     void AnalyticCorrelation(
         double* corr,
         double fwd1,
         double vol1,
         const DoubleArray& x_T1,
         const DoubleArray& fixingTimes1,
         const DoubleArray& fixingPrices1,
         const DoubleArray& expAlpha_T1,
         const DoubleArray& expBeta_T1,
         const DoubleArray& expAlpha_t1,
         const DoubleArray& expBeta_t1,
         const IntArray& intervalEndpoints1,
         double fwd2,
         double vol2,
         const DoubleArray& x_T2,
         const DoubleArray& fixingTimes2,
         const DoubleArray& fixingPrices2,
         const DoubleArray& expAlpha_T2,
         const DoubleArray& expBeta_T2,
         const DoubleArray& expAlpha_t2,
         const DoubleArray& expBeta_t2,
         const IntArray& intervalEndpoints2,
         const DoubleArray& boundaryTimes,
         const DoubleArray& expAlpha_S,
         const DoubleArray& expBeta_S) const;

     void AnalyticExpectations(
         DoubleArray* expectations,
         int fixingIndex,
         double t,
         double expAlpha_t,
         double expBeta_t,
         const IntArray& intervalEndpoints,
         const DoubleArray& boundaryTimes,
         const DoubleArray& expAlpha_S,
         const DoubleArray& expBeta_S) const;


};

typedef smartConstPtr<EnergyInstVolBase> EnergyInstVolBaseConstSP;
typedef smartPtr<EnergyInstVolBase> EnergyInstVolBaseSP;

// support for wrapper class
typedef MarketWrapper<EnergyInstVolBase> EnergyInstVolBaseWrapper;
    
DRLIB_END_NAMESPACE

#endif

