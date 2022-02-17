//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyImpliedVolSurface.hpp
//
//   Description : 
//                 Energy Implied vol surface. Based on drcommodityvolsurfacei.cpp
//                 and drcommodityImpliedvolsurface.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------
#ifndef ENERGY_IMPLIED_VOL_SURFACE_HPP
#define ENERGY_IMPLIED_VOL_SURFACE_HPP

#include "edginc/EnergyVolBase.hpp"

#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyVolCurve.hpp"
#include "edginc/DoubleMatrix.hpp"

#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/VolParallel.hpp"

#include <string>
#include <map>
using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyImpliedVolSurface: public EnergyVolBase,
                               virtual public ITweakableWithRespectTo<VolParallel>
{
public:

    static CClassConstSP const TYPE;

    friend class EnergyImpliedVolSurfaceHelper;
    friend class EnergyImpliedVolSurfaceWrappper;

    virtual ~EnergyImpliedVolSurface();

    void validatePop2Object();

    double getATMFwd(const DateTime& expiry) const;
    double getATMVol(const DateTime& expiry) const;
    double getATMVol(const DateTime& expiry, const DateTime& futureExpiry) const;
    double getSmileVolByStrike(const DateTime& expiry, double strike) const;
    double getSmileVolByDelta(const DateTime& expiry, double delta) const;
    
    // tweaking vol surface
    virtual string sensName(const VolParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>&);

    void getMarket(const IModel* model, const MarketData* market);
    void popVolSurface();

private:

    EnergyImpliedVolSurface();

    EnergyUnderlyerWrapper        energyUnderlyerWrapper;
    EnergyFuturesCurveWrapper    curveWrapper;

    CStringArray expiryLabels;
    CDoubleArray atmVols;
    CDoubleMatrix buckets;
    CDoubleMatrix smileVols;
    
    map<string, EnergyVolCurveSP> volSurface; // $unregistered
    typedef map<string, EnergyVolCurveSP>::const_iterator VolSurfaceIter;
};

typedef smartConstPtr<EnergyImpliedVolSurface> EnergyImpliedVolSurfaceConstSP;
typedef smartPtr<EnergyImpliedVolSurface> EnergyImpliedVolSurfaceSP;

// support for wrapper class
typedef MarketWrapper<EnergyImpliedVolSurface> EnergyImpliedVolSurfaceWrapper;
    
DRLIB_END_NAMESPACE

#endif
