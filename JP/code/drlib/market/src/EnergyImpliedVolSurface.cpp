//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyImpliedVolSurface.hpp
//
//   Description : Energy Implied vol surface. Based on drcommodityvolsurfacei.cpp
//                 and drcommodityImpliedvolsurface.cpp in FXLIB.
//                 
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EnergyImpliedVolSurface.hpp"
#include "edginc/EnergyVolCurve.hpp"
#include "edginc/EnergyVolSurfaceUtils.hpp"


DRLIB_BEGIN_NAMESPACE

EnergyImpliedVolSurface::EnergyImpliedVolSurface():EnergyVolBase(TYPE)
{
}

EnergyImpliedVolSurface::~EnergyImpliedVolSurface()
{
}

void EnergyImpliedVolSurface::validatePop2Object()
{
    
    static const string method = "EnergyImpliedVolSurface::validatePop2Object";

    int numExpiries = expiryLabels.getLength();

	if (numExpiries != atmVols.size())
        throw ModelException(method,"Inconsistent number of expiries for atm vols");
    if (numExpiries != buckets.numRows())
        throw ModelException(method,"Inconsistent number of expiries for buckets");
    if (buckets.numRows()!=smileVols.numRows() || buckets.numCols()!=smileVols.numCols())
        throw ModelException(method,"Inconsistent number of buckets and smileVols");
}

void EnergyImpliedVolSurface::popVolSurface()
{
    static const string method = "EnergyImpliedVolSurface::popVolSurface()";

    int i,j;
    int numLevels = buckets.numCols();
    DoubleArray bucketLevels, volLevels;
    double aFwdRate, time;
	double atmBucket, bucket;
	bool notSaved;
   
    int numExpiries = expiryLabels.getLength();
    for (i=0; i<numExpiries; ++i)
    {
		notSaved = true;
        if (baseDate <= energyUnderlyerWrapper->expiryDate(expiryLabels[i]) )
        {
            bucketLevels.clear();
            volLevels.clear();

			time = energyUnderlyerWrapper->expiryDate(expiryLabels[i],1).daysDiff(baseDate)/365.0;
			if (Maths::isNegative(time)) time = 0.0;

            if ( (aFwdRate = curveWrapper->getFwdRateForLabel(expiryLabels[i])) < 0 )
                throw ModelException(method,"No future price for"+expiryLabels[i]);

			atmBucket = EnergyVolSurfaceUtils::ConvertStrikeToDelta(1,time,1,atmVols[i]);
            for (j=0; j<numLevels; ++j)
            {
                // deltas are ascending
			
				bucket = buckets[j][i];
				if(atmBucket<bucket && notSaved)
				{
				    bucketLevels.push_back(atmBucket);
					volLevels.push_back(atmVols[i]);
					notSaved = false;
				}

				bucketLevels.push_back(bucket);
                volLevels.push_back(smileVols[j][i]);
            }
	

            EnergyVolCurveSP aVolCurve(new EnergyVolCurve());
            
			// use option rule to get option expiry date
            aVolCurve->build(EnergyVolCurve::kDelta,
                             time, aFwdRate,bucketLevels, volLevels);
            volSurface[expiryLabels[i]] = aVolCurve;
        }
    }
}

double EnergyImpliedVolSurface::getATMFwd(const DateTime& expiry) const
{
    static const string method = "EnergyImpliedVolSurface::getATMFwd";
    
    string aLabel = energyUnderlyerWrapper->calculateContractLabel(expiry,EnergyUnderlyer::INDEX).asString();
    VolSurfaceIter iter = volSurface.find(aLabel);
	if ( iter == volSurface.end() )
        throw ModelException(method,"No vol curve for"+aLabel);

    return dynamic_cast<EnergyVolCurve*>(iter->second.get())->getATMFwd();
}


double EnergyImpliedVolSurface::getATMVol(const DateTime& expiry) const
{
    static const string method = "EnergyImpliedVolSurface::getATMVol";
    
    string aLabel = energyUnderlyerWrapper->calculateContractLabel(expiry,EnergyUnderlyer::INDEX).asString();
    VolSurfaceIter iter = volSurface.find(aLabel);
	if (iter == volSurface.end())
        throw ModelException(method, "No vol curve for " + aLabel);

    return dynamic_cast<EnergyVolCurve*>(iter->second.get())->getATMVol();
}

double EnergyImpliedVolSurface::getATMVol(const DateTime& expiry, const DateTime& futureExpiry) const
{
    return getATMVol(expiry);
}

double EnergyImpliedVolSurface::getSmileVolByStrike(const DateTime& expiry, double strike) const
{
    static const string method = "EnergyImpliedVolSurface::getSmileVolByStrike";

    string aLabel = energyUnderlyerWrapper->calculateContractLabel(expiry,EnergyUnderlyer::INDEX).asString();
    VolSurfaceIter iter = volSurface.find(aLabel);
	if (iter == volSurface.end())
        throw ModelException(method, "No vol curve for " + aLabel);

    return dynamic_cast<EnergyVolCurve*>(iter->second.get())->getSmileVolByStrike(strike);
}

double EnergyImpliedVolSurface::getSmileVolByDelta(const DateTime& expiry, double delta) const
{
    static const string method = "EnergyImpliedVolSurface::getSmileVolByDelta";

    string aLabel = energyUnderlyerWrapper->calculateContractLabel(expiry,EnergyUnderlyer::INDEX).asString();
    VolSurfaceIter iter = volSurface.find(aLabel);
	if (iter == volSurface.end())
        throw ModelException(method, "No vol curve for " + aLabel);
    return dynamic_cast<EnergyVolCurve*>(iter->second.get())->getSmileVolByDelta(delta);
}
    
// This is VolParalell tweakable

string EnergyImpliedVolSurface::sensName(const VolParallel*) const
{
    return getName();
}

TweakOutcome EnergyImpliedVolSurface::sensShift(const PropertyTweak<VolParallel>& tweak)
{
    static const string method = "EnergyImpliedVolSurface::sensShift()";
    
    double shiftSize = tweak.coefficient;
    
    int i,j;
    int numLevels = buckets.numCols();
    DoubleArray bucketLevels, volLevels;
    double aFwdRate, time;
	double atmBucket, bucket;
	bool notSaved;
   
    int numExpiries = expiryLabels.getLength();
    for (i=0; i<numExpiries; ++i)
    {
		notSaved = true;
        if (baseDate <= energyUnderlyerWrapper->expiryDate(expiryLabels[i]) )
        {
            bucketLevels.clear();
            volLevels.clear();

			time = energyUnderlyerWrapper->expiryDate(expiryLabels[i],1).daysDiff(baseDate)/365.0;
			if (Maths::isNegative(time)) time = 0.0;

            if ( (aFwdRate = curveWrapper->getFwdRateForLabel(expiryLabels[i])) < 0 )
                throw ModelException(method,"No future price for"+expiryLabels[i]);

			atmBucket = EnergyVolSurfaceUtils::ConvertStrikeToDelta(1,time,1,atmVols[i]);
            for (j=0; j<numLevels; ++j)
            {
                // deltas are ascending
			
				bucket = buckets[j][i];
				if(atmBucket<bucket && notSaved)
				{
				    bucketLevels.push_back(atmBucket);
					volLevels.push_back(atmVols[i]*(1.0+shiftSize));
					notSaved = false;
				}

				bucketLevels.push_back(bucket);
                volLevels.push_back(smileVols[j][i]*(1.0+shiftSize));
            }
	

            EnergyVolCurveSP aVolCurve(new EnergyVolCurve());
            
			// use option rule to get option expiry date
            aVolCurve->build(EnergyVolCurve::kDelta,
                             time, aFwdRate,bucketLevels, volLevels);
            volSurface[expiryLabels[i]] = aVolCurve;
        }
    }

    return TweakOutcome(tweak.coefficient, false);
}

/***
void EnergyImpliedVolSurface::sensRestore(EnergyVega* shift)
{
}
*******/

void EnergyImpliedVolSurface::getMarket(const IModel* model, const MarketData* market)
{
    energyUnderlyerWrapper.getData(model, market);
    curveWrapper.getData(model, market);
    popVolSurface();
}

class EnergyImpliedVolSurfaceHelper
{
public:

       static IObject* defaultEnergyImpliedVolSurface()
       {
         return new EnergyImpliedVolSurface();
       }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyImpliedVolSurface, clazz);
        SUPERCLASS(EnergyVolBase);
        IMPLEMENTS(ITweakableWithRespectTo<VolParallel>);
        EMPTY_SHELL_METHOD(defaultEnergyImpliedVolSurface);
        FIELD(energyUnderlyerWrapper,  "Energy Underlyer Wrapper");
        FIELD(curveWrapper,  "Energy Futures Curve Wrapper");
        FIELD(expiryLabels,     "Contract Expiry Labels");
		FIELD(atmVols,"ATM Vol Array");
        FIELD(buckets,  "Buckets Level Matrix");
        FIELD(smileVols,"Vol Matrix");
    }
};

CClassConstSP const EnergyImpliedVolSurface::TYPE = CClass::registerClassLoadMethod(
        "EnergyImpliedVolSurface", typeid(EnergyImpliedVolSurface), EnergyImpliedVolSurfaceHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyImpliedVolSurfaceWrapper);

// Addins....

class GetEnergySmileVolByDeltaAddin : public CObject 
{

public:

    static CClassConstSP const TYPE;

    // addin parameters
    EnergyImpliedVolSurfaceSP theVolObj; // a handle 
    DateTime expiryDate;
    double delta;


    static IObjectSP getSmileVolDelta(GetEnergySmileVolByDeltaAddin *params)
        {

        static const string method = "GetEnergySmileVolByDeltaAddin::getSmileVolDelta";
        

        if (EnergyImpliedVolSurface::TYPE->isInstance(params->theVolObj))
        {
            EnergyImpliedVolSurface& aVolObj = dynamic_cast<EnergyImpliedVolSurface&>(*params->theVolObj);

            double output = aVolObj.getSmileVolByDelta(params->expiryDate, params->delta);
            CDoubleSP res(CDouble::create(output));
            return IObjectSP(res.clone());
        }
        else
        {
            throw ModelException("Objecy is not an Energy Implied Vol Surface");
            
        }
    
    }
    
    GetEnergySmileVolByDeltaAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        REGISTER(GetEnergySmileVolByDeltaAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetEnergySmileVolByDeltaAddin);
        // order of registration effects order of parameters in addin function
        FIELD(theVolObj, "Handle to EnergyImpliedVolSurface object");
        FIELD(expiryDate, "Expiry");
        FIELD(delta, "Delta");

        Addin::registerClassObjectMethod("GET_ENERGY_SMILE_VOL_BY_DELTA",
                                         Addin::MARKET,
                                         "Get energy Smile vol by Delta",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getSmileVolDelta);
    }
    
    static IObject* defaultGetEnergySmileVolByDeltaAddin()    
    {
        return new GetEnergySmileVolByDeltaAddin();
    }
};

CClassConstSP const GetEnergySmileVolByDeltaAddin::TYPE = CClass::registerClassLoadMethod(
                        "GetEnergySmileVolByDeltaAddin", typeid(GetEnergySmileVolByDeltaAddin),
                         GetEnergySmileVolByDeltaAddin::load);

class GetEnergySmileVolByStrikeAddin : public CObject 
{

public:
    static CClassConstSP const TYPE;

    // addin parameters
    EnergyImpliedVolSurfaceSP theVolObj; // a handle 
    DateTime expiryDate;
    double strike;


    static IObjectSP getSmileVolStrike(GetEnergySmileVolByStrikeAddin *params)
        {

        static const string method = "GetEnergySmileVolByStrikeAddin::getSmileVolStrike";
        

        if (EnergyImpliedVolSurface::TYPE->isInstance(params->theVolObj))
        {
            EnergyImpliedVolSurface& aVolObj = dynamic_cast<EnergyImpliedVolSurface&>(*params->theVolObj);

            double output = aVolObj.getSmileVolByStrike(params->expiryDate, params->strike);
            CDoubleSP res(CDouble::create(output));
            return IObjectSP(res.clone());
        }
        else
        {
            throw ModelException("Objecy is not an Energy Implied Vol Surface");
            
        }
    
    }
    
    GetEnergySmileVolByStrikeAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        REGISTER(GetEnergySmileVolByStrikeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetEnergySmileVolByStrikeAddin);
        // order of registration effects order of parameters in addin function
        FIELD(theVolObj, "Handle to EnergyImpliedVolSurface object");
        FIELD(expiryDate, "Expiry");
        FIELD(strike, "Strike");

        Addin::registerClassObjectMethod("GET_ENERGY_SMILE_VOL_BY_STRIKE",
                                         Addin::MARKET,
                                         "Get energy Smile vol by Strike",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getSmileVolStrike);
    }
    
    static IObject* defaultGetEnergySmileVolByStrikeAddin()    
    {
        return new GetEnergySmileVolByStrikeAddin();
    }
};

CClassConstSP const GetEnergySmileVolByStrikeAddin::TYPE = CClass::registerClassLoadMethod(
                        "GetEnergySmileVolByStrikeAddin", typeid(GetEnergySmileVolByStrikeAddin),
                         GetEnergySmileVolByStrikeAddin::load);




DRLIB_END_NAMESPACE
