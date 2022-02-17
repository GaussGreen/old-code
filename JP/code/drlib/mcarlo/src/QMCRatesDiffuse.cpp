//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCRatesDiffuse.cpp (code previously in MCPathConfigSRM.cpp)
//
//   Description : base class for SRM IR diffusion, takes care of dates
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/QMCFXBaseDiffuse.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"


#include <cassert>
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

// make sure all elementary types are initialized to some unreal values
QMCRatesDiffuse::QMCRatesDiffuse() :
    randomIndex(-999),
    fx(NULL),
    saveDomLnMoney(true),
    saveSigmaR(true),
    NbSigmasMax(-999),
    NbSigmasMin(-999),
    todayIndex(-1), //
    calibrateAgainstSwaptionVols(true),
    todayIdx(-999),
    initialized(false),
    assetDatesFixed(false),
//    timeLogic(new QMCHelperCachingTimeLogic),
    diffBound(NULL) // no fx

{}

QMCRatesDiffuse::~QMCRatesDiffuse()
{
/*    delete timeLogic; //FIXME move to SP
    timeLogic = NULL;*/
}
    
const vector<double>* QMCRatesDiffuse::getSigmaFX() 
{ 
    return fx.get() ?  &(fx->getSigmaFX()) : NULL;
}

const vector<double>* QMCRatesDiffuse::getSpotFXFullPath()
{ 
    return fx.get() ?  &(fx->requestFullSpotPath()) : NULL; 
}

void QMCRatesDiffuse::setFXBase(QMCFXBaseDiffuseSP _fx)
{
    fx=_fx;
    getDiffusionBound()->add(
            dynamic_cast<QMCHelperBoundedDiffusion *>(fx->getDiffusionBound())); // notify that we're given FX
}

double* QMCRatesDiffuse::getInternalDFArray()
{
    double * result= (df.empty() ? NULL : & df[0]);
    return result;
} // this asset stores DFs

size_t QMCRatesDiffuse::registerYCFlavor(IYieldCurveConstSP yc)
{
    static const string method("QMCRatesDiffuse::registerYCFlavor");
    if (yc.get() == NULL)
        throw ModelException(method, "Error: YieldCurveConstSP points to NULL");

    for(size_t i=0; i < ycForwardsDB.size(); ++i)
	{
        if (ycForwardsDB[i].first.get() == yc.get())
            return i;
    }

	if (initialized) // it's too late for adding new YC flavor at this point
        throw ModelException(method, "Logical error: out of order function call.");

    ycForwardsDB.push_back(make_pair(yc, vector<double>()));
    return ycForwardsDB.size()-1;
}

SpotIdx QMCRatesDiffuse::getFwdIdx2EdfIdx(FwdIdx idx) // translates FwdIdx into RequestedIdx
{
    return getTimeLogic()->getReqEDFIdx(idx);
}


DRLIB_END_NAMESPACE


