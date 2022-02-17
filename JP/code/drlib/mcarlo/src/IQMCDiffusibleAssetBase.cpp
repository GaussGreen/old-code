//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : IQMCDiffusibleAssetBase.cpp
//
//   Description : Non-virtual and default methods common to all diffusible assets
//
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
//#include "edginc/IQMCRNGManager.hpp"

#include "edginc/QMCHelperCachingTimeLogic.hpp"
//#include "edginc/QMCHelperDateTimeCache.hpp"
#include "edginc/QMCHelperDateTimeCacheNew.hpp"

#include "edginc/IQMCDiffusibleAssetBase.hpp"

DRLIB_BEGIN_NAMESPACE

IQMCDiffusibleAsset::IQMCDiffusibleAsset() 
    : timeLogic( new QMCHelperCachingTimeLogic(IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCacheNew),
                                                IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCacheNew),
                                                IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCacheNew))
                   )
{}

/** Sophisticated clients may override the default (fast/memory inefficient) timeLogic impl. */
IQMCHelperTimeLogicSP IQMCDiffusibleAsset::getTimeLogic() const
{
    return timeLogic;
}
/** Default implementation via timeLogic */
DateTimeArray IQMCDiffusibleAsset::getAssetDates()
{
    return getTimeLogic() ->getAssetDates();
}


/// addSpotDates is a simple redirect to addAggregatedDates; no need to be  virtual
void IQMCDiffusibleAsset::addSpotDates( const DateTimeArray& measurementDates )
{
    addAggregatedDates( measurementDates, DateTimeArray(), DateTimeArray() );
}

/// addForwardDates  is a simple redirect to addAggregatedDates; no need to be  virtual
void IQMCDiffusibleAsset::addForwardDates( const DateTime& measurementDate,
        const DateTimeArray& futureDates )
{
    addAggregatedDates( DateTimeArray(), DateTimeArray( 1, measurementDate ), futureDates );
}

/// addAggregatedDates passes all dates to the timeLogic class.
/// Derived classes may redefine this function if they need to inform other assets about their dates (see CR example)
void IQMCDiffusibleAsset::addAggregatedDates(
    const DateTimeArray& spot,
    const DateTimeArray& fwd,
    const DateTimeArray& fwdfwd )
{
    getTimeLogic() ->addAggregatedDates( spot, fwd, fwdfwd );
    getTimeLogic() ->addAggregatedDates( DateTimeArray(), DateTimeArray(), fwd ); ///< include ForwardDates into ForwardForwardDates
    updateMaxMaturity(spot, fwd, fwdfwd);
}

/// Same function using SPs.: simply redirect to non-SP version
/// Default implementation redirects to non-SP version and updates maxMaturity
void IQMCDiffusibleAsset::addAggregatedDates(
    DateTimeArrayConstSP spot,
    DateTimeArrayConstSP fwd,
    DateTimeArrayConstSP fwdfwd,
    const DateTime& maxDiffDate,
    const DateTime& maxCurveMaturity)
{
    addAggregatedDates(maxDiffDate.getPastDates(*spot), maxDiffDate.getPastDates(*fwd), maxCurveMaturity.getPastDates(*fwdfwd));
    
    getDiffusionBound()->updateMaxDiffDate(maxDiffDate);
    getDiffusionBound()->updateCurveMat(maxCurveMaturity);
}

/// Simple way to update maxMaturity information by looking what is passed to the asset.
void IQMCDiffusibleAsset::updateMaxMaturity(
    const DateTimeArray& spot,
    const DateTimeArray& fwd,
    const DateTimeArray& fwdfwd )
{
    if (! spot.empty())
        getDiffusionBound()->updateMaxDiffDate(spot.back());
    if (!fwd.empty())
        getDiffusionBound()->updateMaxDiffDate(fwd.back());
    if (!fwdfwd.empty())
        getDiffusionBound()->updateCurveMat(fwdfwd.back());
}

SpotIdx IQMCDiffusibleAsset::getSpotIndex( const DateTime& measurementDate ) const
{
    return getTimeLogic() ->getDFIdx( measurementDate );
}

SpotIdx IQMCDiffusibleAsset::getForwardIndex( const DateTime& forwardDate ) const
{
    return getTimeLogic() ->getReqEDFIdx( forwardDate );
}

FwdIdx IQMCDiffusibleAsset::getForwardForwardIndex( const DateTime& forwardDate ) const
{
    return getTimeLogic() ->getFwdEDFIdx( forwardDate );
}

/// Return requested Spot dates trimmed by maxMaturity calculated from SVGens posted for this and dependent assets
DateTimeArray IQMCDiffusibleAsset::getSpotDates()
{
    DateTime maxDate = getDiffusionBound()->getMaxDiffDate();
    DateTimeArray dates = getTimeLogic()->getDFDates();
    return maxDate.getPastDates( dates );
}

/// Return requested ForwardDates trimmed by maxMaturity
DateTimeArray IQMCDiffusibleAsset::getForwardDates()
{
    DateTime maxDate = getDiffusionBound()->getMaxDiffDate();
    DateTimeArray dates = getTimeLogic()->getReqEDFDates();
    return maxDate.getPastDates( dates );
}

/// Reqturn requested ForwardForward dates trimmed by maxCurveMaturity
DateTimeArray IQMCDiffusibleAsset::getForwardForwardDates()
{
    DateTime maxDate = getDiffusionBound()->getMaxCurveMat();
    DateTimeArray dates = getTimeLogic()->getFwdEDFDates();
    return maxDate.getPastDates( dates );
}

size_t IQMCDiffusibleAsset::getNumESDFDates()
{
    return getForwardDates().size();
} // classes will override for eff.

DRLIB_END_NAMESPACE
