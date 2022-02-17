//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : DDEParams.hpp
//
//   Description : parameters for DDE calculations: asset vol, sprd mean reversion etc
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#ifndef EDG_DDE_PARAMS_H
#define EDG_DDE_PARAMS_H

#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Random.hpp"

DRLIB_BEGIN_NAMESPACE


//////////////////////////////////////////////
//											//
//			DDEParams						//
//											//
//  implement asset vol and mean reversion	//
//  as linear function of log(clean sprd)	//
//											//
//////////////////////////////////////////////

class MARKET_DLL DDEParams: public CObject
{
public:
// not done yet
    static CClassConstSP const TYPE;
	static const double DEF_E2C_ASSET_VOL;
	static const double DEF_E2C_BARRIER_VOL;
	static const double DEF_E2C_HORIZON;
   	static const double DEF_MEAN_REVERSION;
	static const double DEF_SPREAD_Q;
	static const double DEF_RISKYNESS;
	static const double DEF_CORR_BASE_SPRD;

    friend class DDEParamsHelper;

	double getAssetVol(double spread) const;
	double getMeanReversion(double spread) const;
	double spreadVol(double spread) const;
	double spreadEqCorr(double spread) const;

	// used to calibrate asset/debt ratio
	static double calibFuncE2C(double lnd, void *calibDataMem);

	virtual void validatePop2Object();

private:
    DDEParams();

public:
	double			assetVol100bp;
	double			volSlope;			// asset vol slope per 100bp change in clean spread
	double			horizonE2C;			// E2C params for spread vol calc
	double			spreadQ;			// for spread dynamics. 1/0 is normal/lognormal
	double			mr100bp;
	double			mrSlope;			// mean reversion slope per 100bp change in clean spread
	double			corrBaseSprd;		// characteristic spread for eq-sprd corr
	bool			isRelatedVol;		// if merge vol and eq-sprd corr
};

typedef smartConstPtr<DDEParams> DDEParamsConstSP;
typedef smartPtr<DDEParams>      DDEParamsSP;

DRLIB_END_NAMESPACE
#endif
