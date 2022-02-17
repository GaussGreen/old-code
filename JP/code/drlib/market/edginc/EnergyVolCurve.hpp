//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyVolCurve.hpp
//
//   Description : Created after DRCommodityVolCurve.hpp
//
//   Author      : Sean Chen
//
//   Date        : July 27, 2005
//
//----------------------------------------------------------------------------

#ifndef _energyvolcurve_
#define _energyvolcurve_

#include "edginc/MarketObject.hpp"
#include "edginc/Spline.hpp"
#include "edginc/Interpolator.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyVolCurve : public MarketObject
{
	friend class EnergyVolCurveHelper;

    public:

        static CClassConstSP const TYPE;

		 EnergyVolCurve();
        ~EnergyVolCurve();

        enum EBucketFormat { kStrike, kDelta };
        
        void build(
            EBucketFormat format,
            double time,
            double fwd,
            const DoubleArray& buckets,
            const DoubleArray& vols);

        double getATMFwd() const;
        double getATMVol() const;
        double getSmileVolByStrike(double strike) const;
        double getSmileVolByDelta(double delta) const; 
	
	    string getName() const { return "";}
        Interpolator::InterpolantConstSP getInterpolant() const { return interpolantSP;}

	
    private:

        EBucketFormat                bucketFormat; // $unregistered
        double                       expiryTime; // $unregistered
        double                       aTMFwd; // $unregistered
        Interpolator::InterpolantConstSP     interpolantSP; // $unregistered
};

typedef smartConstPtr<EnergyVolCurve> EnergyVolCurveSConstSP;
typedef smartPtr<EnergyVolCurve> EnergyVolCurveSP;

DRLIB_END_NAMESPACE

#endif
