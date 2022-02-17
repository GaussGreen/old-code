//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedDVF.cpp
//
//   Description : interface for processed local vol object
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/Asset.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MarketDataFetcherLN.hpp"

DRLIB_BEGIN_NAMESPACE

static void volProcessedDVFLoad(CClassSP& clazz){
    REGISTER(CVolProcessedDVF, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
}

CVolProcessedDVF::CVolProcessedDVF(const CClassConstSP& clazz): CObject(clazz){}

CClassConstSP const CVolProcessedDVF::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedDVF", typeid(CVolProcessedDVF), volProcessedDVFLoad);

class DVFInterpVolAddin: public CObject{
private:

public:
    static CClassConstSP const TYPE;

    // addin parameters
    DateTime            startDate; 
    DateTimeArray       maturities;
    CDoubleArray        strikes;
    CAssetWrapper       asset;
    string              volType;
    MarketDataSP        market;

    bool                useUnscaledTimeTweak;
    bool                useTweakingForTimeDerivs;
    bool                useMidPoint;
    double              tweakStrikeUnscaled;
    double              tweakTimeUnscaled;
    double              probDensRatioMin;

    bool                isIntraDayInterp;

    // return the vols from the given date to each future date
    static IObjectSP calcVolsDVF(DVFInterpVolAddin* params){
        static const string method = "DVFInterpVolAddin::calcVolsDVF";
        // if market supplied use it
        if (params->market.get()){
            IModelSP npm(new NonPricingModel());
            MarketDataFetcherLNSP mdf(new MarketDataFetcherLN(params->volType));
            npm = IModelSP(new NonPricingModel(mdf));
            params->asset.getData(npm.get(), params->market.get());
        }

        // populate lattice of strikes
        int iMat;
        int iStrike;
        int nbMat = params->maturities.size();
        int nbStrikes = params->strikes.size();

        /*
        vector<int> sizes(nbMat, nbStrikes);
        CLatticeDoubleSP latticeStrikes(new CLatticeDouble(sizes));
        CLatticeDoubleSP locVolLattice(new CLatticeDouble(sizes));
        for(iMat = 0; iMat < nbMat; iMat++ ) {
            for(iStrike = 0; iStrike < nbStrikes; iStrike++ ) {
                (*latticeStrikes)[iMat][iStrike] = params->strikes[iStrike];
            }
        }
        */

        // create a simple vol request
        LocVolRequestSP volRequest(new LocVolRequest(
                    params->startDate,                  // startDate
                    false,                              // fwdStarting
                    params->useUnscaledTimeTweak,       // useUnscaledTimeTweak
                    params->useTweakingForTimeDerivs,   // useNextStepDerivs
                    params->useMidPoint,                // useMidPoint
		            params->tweakStrikeUnscaled,        // strikeTweakUnscaled
                    params->tweakTimeUnscaled,          // timeTweakUnscaled
                    params->probDensRatioMin));           // probDensRatioMin

        // process the vol
        CVolProcessedSP procVol(params->asset->getProcessedVol(volRequest.get()));
        CVolProcessedDVF& procVolDVF = dynamic_cast<CVolProcessedDVF&>(*procVol);

        // then interpolate
        CDoubleMatrixSP locVolMatrix(new CDoubleMatrix(nbStrikes,nbMat));
        for(iMat = 0; iMat < nbMat; iMat++ ) {
            for(iStrike = 0; iStrike < nbStrikes; iStrike++ ) {
                (*locVolMatrix)[iStrike][iMat]
                    = procVolDVF.CalcLocVol( params->maturities[iMat], params->strikes[iStrike], params->isIntraDayInterp);
            }
            // procVolDVF.CalcLocVol( latticeStrikes.get(), params->maturities, locVolLattice.get());
        }
        
        /*
        // lattice to matrix
        CDoubleMatrixSP locVolMatrix(new CDoubleMatrix(nbStrikes,nbMat));
        for(iMat = 0; iMat < nbMat; iMat++ ) {
            for(iStrike = 0; iStrike < nbStrikes; iStrike++ ) {
                (*locVolMatrix)[iStrike][iMat] = (*locVolLattice)[iMat][iStrike];
            }
        }
        */

        return IObjectSP(locVolMatrix.clone());
    }
    
    DVFInterpVolAddin(): CObject(TYPE), 
                         useUnscaledTimeTweak(true),
                         useTweakingForTimeDerivs(false),
                         useMidPoint(true),
                         tweakStrikeUnscaled(0.005),
                         tweakTimeUnscaled(0.001),
                         probDensRatioMin(0.01),
                         isIntraDayInterp(true) {}
                                    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DVFInterpVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLNInterpVolAddin);
        FIELD(startDate, "startDate");
        FIELD(maturities, "maturities");
        FIELD(strikes, "strikes");
        FIELD(asset, "asset");
        FIELD(volType, "volType");
        FIELD(market, "market data cache");

        FIELD(useUnscaledTimeTweak, "useUnscaledTimeTweak (true)");
        FIELD_MAKE_OPTIONAL(useUnscaledTimeTweak);
        FIELD(useTweakingForTimeDerivs, "useTweakingForTimeDerivs (false)");
        FIELD_MAKE_OPTIONAL(useTweakingForTimeDerivs);
        FIELD(useMidPoint, "useMidPoint (true)");
        FIELD_MAKE_OPTIONAL(useMidPoint);
        FIELD(tweakStrikeUnscaled, "tweakStrikeUnscaled (0.005)");
        FIELD_MAKE_OPTIONAL(tweakStrikeUnscaled);
        FIELD(tweakTimeUnscaled, "tweakTimeUnscaled (0.001)");
        FIELD_MAKE_OPTIONAL(tweakTimeUnscaled);
        FIELD(probDensRatioMin, "probDensRatioMin (0.01)");
        FIELD_MAKE_OPTIONAL(probDensRatioMin);
        FIELD(isIntraDayInterp, "True: linearize fwd intraday, False: step fwd");
        FIELD_MAKE_OPTIONAL(isIntraDayInterp);

        Addin::registerClassObjectMethod("LOCAL_VOL",
                                         Addin::MARKET,
                                         "Calculates local vol",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)calcVolsDVF);
    }
    
    static IObject* defaultLNInterpVolAddin(){
        return new DVFInterpVolAddin();
    }
 
};

CClassConstSP const DVFInterpVolAddin::TYPE = CClass::registerClassLoadMethod(
    "DVFInterpVolAddin", typeid(DVFInterpVolAddin), load);

DRLIB_END_NAMESPACE
