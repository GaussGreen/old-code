//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SubReprice.cpp
//
//   Description : Captures whether paths in monte carlo simulation can be 
//                 skipped when calculating sensitivities
//
//   Author      : Mark A Robson
//
//   Date        : 29 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SubReprice.hpp"
//#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MCProduct.hpp"

DRLIB_BEGIN_NAMESPACE

ISubReprice::~ISubReprice(){}
ISubReprice::IOneSidedGreekParamsSet::~IOneSidedGreekParamsSet(){}
ISubReprice::ITwoSidedGreekParamsSet::~ITwoSidedGreekParamsSet(){}
ISubReprice::IXGammaParamsSet::~IXGammaParamsSet(){}

/** Constructor. 'store' should be set to true if the object needs
    to record values */
ISubReprice::Rainbow::Rainbow(int mode, int numPaths): 
    isStoringMode((mode & IMCProduct::QUICK_X_GAMMA_BIT) != 0){
    if (isStoringMode){
        perfDiff = DoubleArrayArraySP(new DoubleArrayArray(0));
        perfDiff->reserve(numPaths);
    }
}
/** Returns a deep copy of this object. */
ISubReprice* ISubReprice::Rainbow::clone() const{
    Rainbow* copy = new Rainbow(0, 0);
    copy->isStoringMode = isStoringMode;
    copy->perfDiff = perfDiff; // shallow copy
    copy->asset1 = asset1;
    copy->asset2 = asset2;
    copy->perf1Tol = perf1Tol;
    copy->perf2Tol = perf2Tol;
    return copy;
}

/** Sets object in mode ready to pricing for a [1st order] greek
    which is done by a one sided taylor expansion. Use NULL for params */
void ISubReprice::Rainbow::setForOneSidedGreek(
    IOneSidedGreekParamsSet* params){
    isStoringMode = false;
}

/** Sets object in mode ready to pricing for a [1st or 2nd order] greek
    which is done by a two sided taylor expansion. Use NULL for params */
void ISubReprice::Rainbow::setForTwoSidedGreek(
    ITwoSidedGreekParamsSet* params){
    isStoringMode = false;
}

ISubReprice::Rainbow::XGammaParamsSet::XGammaParamsSet(
    int             asset1,
    double          perf1Tol,
    int             asset2,
    double          perf2Tol): asset1(asset1), asset2(asset2),
                               perf1Tol(perf1Tol), perf2Tol(perf2Tol){}

ISubReprice::Rainbow::XGammaParamsSet::XGammaParamsSet(
            const DoubleArray&               maxCmptScalingFactor,
            const MCPathGenerator*           futurePathGen,
            const IMultiFactors*             assets, 
            const ScalarShiftArray&          sens){
    int                   assetIDs[2];    // which assets are being tweaked
    double                tolerances[2]; // tolerances for these assets
    crossGammaTolerances(maxCmptScalingFactor, futurePathGen, assets,
                         sens, assetIDs, tolerances);
    this->asset1 = assetIDs[0];
    this->perf1Tol = tolerances[0];
    this->asset2 = assetIDs[1];
    this->perf2Tol = tolerances[1];
}

ISubReprice::Rainbow::XGammaParamsSet::~XGammaParamsSet(){}

/** Sets object in mode ready to pricing for cross gamma. Use
    XGammaParamsSet for an instance of IXGammaParamsSet */
void ISubReprice::Rainbow::setForXGamma(IXGammaParamsSet* params){
    XGammaParamsSet& myParams = dynamic_cast<XGammaParamsSet&>(*params);
    isStoringMode = false;
    this->asset1 = myParams.asset1;
    this->perf1Tol = myParams.perf1Tol;
    this->asset2 = myParams.asset2;
    this->perf2Tol = myParams.perf2Tol;
}

/** Stores the necessary data for the current path on initial
    pricing run and thus allows methods below to be
    implemented. The method must be invoked for each simulated
    path and called in order. The IndexedPerfList must have already
    been sorted */
void ISubReprice::Rainbow::store(
    const IndexedPerfList&           indexPerfs,
    const IMCProduct*                 product,
    const MCPathGenerator* pathGen){
    if (isStoringMode){
        // get some memory
        perfDiff->push_back(DoubleArray());
        perfDiff->back().resize(indexPerfs.size());
        DoubleArray& perfDiffs = perfDiff->back();
        indexPerfs.storeNearestPerfDiff(perfDiffs);
#if 0
        /* removed for now as divideByPathWiseMaxDriftProduct is not
           correct for implied path generator */
        /* our condition for skipping cross gamma paths is if
           perf Diff(path) > tolerance(path)
           where tolerance(path) = 
           delta S * scalingFactor * maxDriftProduct(path)
           So divide through by maxDriftProduct(path) to make RHS independent
           of path */
        divideByPathWiseMaxDriftProduct(perfDiffs, pathGen, 0);
        // (hard coded to zero'th path)
#endif
    }
}

/** For best or worst of payoffs. Uses the supplied information to
    store data for the current path on initial pricing run and this
    allows methods below to be implemented. The method must be invoked
    for each simulated path and called in order. The performances must not
    be sorted */
void ISubReprice::Rainbow::store(
    const DoubleArray&               performances,
    int                              pickedAsset,
    const IMCProduct*                 product,
    const MCPathGenerator* pathGen){
    if (isStoringMode){
        // get some memory
        perfDiff->push_back(DoubleArray());
        perfDiff->back().resize(performances.size());
        DoubleArray& perfDiffs = perfDiff->back();
        // need to store difference to nearest performer with non zero weight
        double pickedPerf = fabs(performances[pickedAsset]);
        bool   smallestDifSet = false;
        double smallestDif = 0.0;
        for (int i = 0; i < performances.size(); i++){
            if (i != pickedAsset){
                perfDiffs[i] = fabs(performances[i]-pickedPerf);
                // assets with equal performances are of equal 
                // rank (so rank is arbitrary)
                if ((!smallestDifSet || perfDiffs[i] < smallestDif) &&
                    !Maths::isZero(perfDiffs[i])){
                    smallestDif = perfDiffs[i];
                    smallestDifSet = true;
                }
            }
        }
        perfDiffs[pickedAsset] = smallestDif;
#if 0
        /* removed for now as divideByPathWiseMaxDriftProduct is not
           correct for implied path generator */
        // see explanation in above store() method
        divideByPathWiseMaxDriftProduct(perfDiffs, pathGen, 0);
        // (hard coded to zero'th path)
#endif
    }
}

/** Returns false */
bool ISubReprice::Rainbow::firstDerivZero(int path) const{
    return false;
}

/** Returns false */
bool ISubReprice::Rainbow::firstNumericalDerivZero(int path) const{
    return false;
}

/** Returns true if order of assets in rainbow will not change under
    tweak */
bool ISubReprice::Rainbow::crossNumericalDerivZero(int path) const{
    const DoubleArray& thisPerfDiff = (*perfDiff)[path];
    if (Maths::equals(thisPerfDiff[asset1], thisPerfDiff[asset2])){
        double perfDiff = thisPerfDiff[asset1];
        if (Maths::isZero(perfDiff)){
            // all assets equally ranked - so can viewed as being in any order
            // already
            return true;
        }
        /* Can skip code below if we assume that shifts are of the
           same sign (ie both assets up or both assets down) */
#if 0
        // chances are then the two assets are next to each other in
        // terms of ranking by performance
        if (perfDiff > perf1Tol + perf2Tol){
            // both assets move - but the combined effect is not enough
            return true;
        } else {
            return false;
        }
#endif
    }
    if (thisPerfDiff[asset1] > perf1Tol && thisPerfDiff[asset2] > perf2Tol){
        //  order of assets won't change.
        return true;
    }
    return false;
}

/** Returns the largest absolute value in a double array */
double ISubReprice::Rainbow::largestAbsWeight(const DoubleArray& dbles){
    double maxDb = 0.0;
    for (int i = 0; i < dbles.size(); i++){
        double db = fabs(dbles[i]);
        if (db > maxDb){
            maxDb = db;
        }
    }
    return maxDb;
}

/** Returns null SP */
refCountPtr<ISubReprice> ISubReprice::Rainbow::getRepriceForCmpt(){
    return ISubRepriceSP();
}

/** Calculates tolerances needed for Rainbow under cross gamma. It returns
    the id's of the two assets being shifted together with the maximum size
    that they can move (scaled by maxCmptScalingFactor eg 'participation' 
    but note rainbow weights should not be included in maxCmptScalingFactor) */
void ISubReprice::Rainbow::crossGammaTolerances(
    const DoubleArray&      maxCmptScalingFactor,
    const MCPathGenerator*  futurePathGen,
    const IMultiFactors*    assets,        /* should change to IMultiFactors */
    const ScalarShiftArray& sens,
    int                   assetIDs[2],    // which assets are being tweaked
    double                tolerances[2]){ // tolerances for these assets
    static const string method("ISubReprice::Rainbow::crossGammaTolerances");
    try{
        if (sens.size() != 2){
            throw ModelException(method, "Internal error - num shifts != 2");
        }
        // need to know how much each of the two 'forwards' might move by
        for (int i = 0; i < 2; i++){
            // find which asset is sensitive to this market data
            IntArray sensAssets(assets->getSensitiveAssets(
                sens[i].get(), false /* excude phi */));
            if (sensAssets.size() != 1){
                // could happen if fx cross gamma is set up wrongly
                throw ModelException(method,
                                     "Internal error - num sens assets != 1");
            }
            int assetID = sensAssets.front();
            assetIDs[i] = assetID;
            tolerances[i] = assets->assetGetSpot(assetID) * 
                fabs(maxCmptScalingFactor[assetID])*
                sens[i]->getShiftSize() 
#if 1
                /* included until we make pathWiseMaxDriftProduct work for
                   implied */
                * futurePathGen->maxDriftProduct(assetID)
#endif
                ;
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

#if 0
    /* Removed for now as we need a path generator specific method to
       do this calculation. The implementation here is ok for log normal
       but not implied */
/** Calculates the largest 'derivative' of the current simulated path
    wrt to spot price (of specified asset) across the different
    simulation dates. If it returns a value of lambda then the value
    of spot on a simulated path cannot move more than lambda * delta S
    when the initial spot S is moved by delta S. */
double ISubReprice::Rainbow::pathWiseMaxDriftProduct(
    const MCPathGenerator* pathGen,
    int                              iAsset,
    int                              iPath){
    int begin = pathGen->begin(iAsset);
    int end = pathGen->end(iAsset);
    const double* path = pathGen->Path(iAsset, iPath);
    // this is wrong - it should be comparing to spot
    double refLevel = pathGen->refLevel(iAsset, iPath);
    double maxPerf = DBL_EPSILON; // just in case begin == end
    for (int i = begin; i < end; i++){
        double perf = path[i]/refLevel;
        if (perf > maxPerf){
            maxPerf = perf;
        }
    }
    return maxPerf;
}
     
void ISubReprice::Rainbow::divideByPathWiseMaxDriftProduct(
    DoubleArray&                     toScale,
    const MCPathGenerator* pathGen,
    int                              iPath){
    for (int i = 0; i < toScale.size(); i++){
        toScale[i] /= pathWiseMaxDriftProduct(pathGen, i, iPath);
    }
}
#endif

DRLIB_END_NAMESPACE
