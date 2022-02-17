//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : DDEParams.cpp
//
//   Description : parameters for DDE calculations: asset vol, sprd mean reversion etc
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DDEParams.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

const double DDEParams::DEF_E2C_ASSET_VOL	= 0.17;
const double DDEParams::DEF_MEAN_REVERSION	= 0.1;
const double DDEParams::DEF_E2C_BARRIER_VOL = 0.3;	//E2C barrier vol
const double DDEParams::DEF_E2C_HORIZON		= 2.0;
const double DDEParams::DEF_SPREAD_Q		= 0;		// lognormal spread
const double DDEParams::DEF_CORR_BASE_SPRD	= 0.05;		// 500 bp

typedef struct _CalibE2CData
{
	double delta;
	double lnpv0;
} CalibE2CData;


//////////////////////////////////////////////
//											//
//			DDEParams						//
//											//
//////////////////////////////////////////////

DDEParams::DDEParams() 
: CObject(TYPE),
	assetVol100bp(DEF_E2C_ASSET_VOL),
	volSlope(0),
	horizonE2C(DEF_E2C_HORIZON),
	spreadQ(DEF_SPREAD_Q),
	mr100bp(DEF_MEAN_REVERSION),
	mrSlope(0),
	corrBaseSprd(DEF_CORR_BASE_SPRD),
	isRelatedVol(true)
{}

// for root finder for lnd
double DDEParams::calibFuncE2C(double lnd, void *calibDataMem)
{
	double delta = ((CalibE2CData *)calibDataMem)->delta;
	double pv = N1(-0.5*delta + lnd/delta) - exp(lnd) * N1(-0.5*delta - lnd/delta);
	return log(pv) -((CalibE2CData *)calibDataMem)->lnpv0;
}

double DDEParams::getAssetVol(double spread) const
{	
	return assetVol100bp + volSlope * (spread*100 - 1);  
}

double DDEParams::getMeanReversion(double spread) const
{	
	return mr100bp + mrSlope * (spread*100 - 1);  
}

double DDEParams::spreadVol(double spread) const
{
	try {

	static double MIN_LNPV = 100*DBL_EPSILON;

	double T = horizonE2C; 
	double assetVol = getAssetVol(spread);
	double barrVar = DEF_E2C_BARRIER_VOL * DEF_E2C_BARRIER_VOL; // barrier vol^2
	double delta = sqrt(assetVol * assetVol * T + barrVar);

	// root searching for d factor, ie. asset/debt ratio, to calibrate spread
	CalibE2CData calibData;
	calibData.delta = delta;
	calibData.lnpv0 = - spread * T;

	// avoid unnecessary root finding error
	if( calibData.lnpv0 > -MIN_LNPV ) calibData.lnpv0 = -MIN_LNPV;

	double lnd = zbrentUseful(
                  &calibFuncE2C,				/* (I) function to find the root of			*/
                  &calibData,					/* (I) only 1 other parameter to function	*/
                  barrVar,						/* (I) low value for lnd					*/ 
                  barrVar + 3,                  /* (I) high value for lnd					*/
                  barrVar*5e-4);				/* (I) tolerance. about 1bp for vol		*/


	// factor in parenthesis
	double c = (1 - lnd * (1-barrVar/delta/delta) / (exp(lnd - barrVar) - 1)) / delta;
	double factor = 2 * c * N1Density(-0.5*delta+lnd/delta) - exp(lnd) * N1(-0.5*delta - lnd/delta);

	// mean reversion factor, using the default MR to isolate spread vol from MR spec
	double x = DEF_MEAN_REVERSION * T;
	double fx = (x<1e-6?1:((1-exp(-x))/x));
	return assetVol * factor * exp(-calibData.lnpv0)/(-calibData.lnpv0)/fx;

	} catch (exception &e) {
		throw ModelException(e, "DDEParams::spreadVol", "Failed");
	}
}

double DDEParams::spreadEqCorr(double spread) const
{
	return tanh(spread/corrBaseSprd);
}


void DDEParams::validatePop2Object()
{
    static const string method = "DDEParams::validatePop2Object";

    if( Maths::isZero( assetVol100bp + 1 ) ) assetVol100bp = DEF_E2C_ASSET_VOL;
    if( Maths::isZero( volSlope + 1 ) ) volSlope = 0;
	if( Maths::isZero(horizonE2C + 1) ) horizonE2C = DEF_E2C_HORIZON;
	if( Maths::isZero(corrBaseSprd + 1) ) corrBaseSprd = DEF_CORR_BASE_SPRD;
    if( Maths::isZero( mr100bp + 1 ) ) mr100bp = DEF_MEAN_REVERSION;
    if( Maths::isZero( mrSlope + 1 ) ) mrSlope = 0;

    if( volSlope < 0 )
        throw ModelException("DDEAssetVol", "Asset vol slope must be non-negative or is -1 (ie. use default)");
    if( assetVol100bp < volSlope )
        throw ModelException("DDEAssetVol", "Asset vol slope must not results in negative asset vol at spread ~0");
	if( !Maths::isPositive(horizonE2C) )
		throw ModelException(method, "horizonE2C must be >= 0 yr");
    if( mrSlope < 0 )
        throw ModelException("DDEMeanReversion", "Mean reversion slope must be non-negative or is -1 (ie. use default)");
    if( mr100bp < mrSlope )
        throw ModelException("DDEMeanReversion", "Mean reversion slope must not results in negative mr at spread ~0");
	if( !Maths::isZero(spreadQ) )
		throw ModelException(method, "Only support lognormal credit spread. spreadQ must be 0");
	if( !Maths::isPositive(corrBaseSprd) )
		throw ModelException(method, "corrBaseSprd must be positive");

}

class DDEParamsHelper{
public:
	
	/** Invoked when Class is 'loaded' */

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DDEParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDEParams);
		FIELD(assetVol100bp, 	"asset vol for name with clean spread 100bp");
		FIELD_MAKE_OPTIONAL(assetVol100bp);
		FIELD(volSlope,			"asset vol slope per 100bp change in clean spread");
		FIELD_MAKE_OPTIONAL(volSlope);
        FIELD(horizonE2C,		"E2C horizon");
        FIELD_MAKE_OPTIONAL(horizonE2C);
		FIELD(spreadQ,			"spread Q. 1/0 for normal/lognormal");
		FIELD_MAKE_OPTIONAL(spreadQ);
		FIELD(mr100bp,		"mean reversion for name with clean spread 100bp");
		FIELD_MAKE_OPTIONAL(mr100bp);
		FIELD(mrSlope,		"mean reversion slope per 100bp change in clean spread");
		FIELD_MAKE_OPTIONAL(mrSlope);
		FIELD(corrBaseSprd,		"characteristic spread for eq-sprd corr");
        FIELD_MAKE_OPTIONAL(corrBaseSprd);
		FIELD(isRelatedVol,		"if merge spread vol with eq spread corr");
        FIELD_MAKE_OPTIONAL(isRelatedVol);
    }

    static IObject* defaultDDEParams(){
        return new DDEParams();
    }
};


CClassConstSP const DDEParams::TYPE = CClass::registerClassLoadMethod(
              "DDEParams", typeid(DDEParams), DDEParamsHelper::load);

DRLIB_END_NAMESPACE
