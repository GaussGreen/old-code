//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCDiffusibleAssetBase.hpp
//
//   Description : An interface of internal hooks for asset diffusion and
//                 diffused data retrieval
//

//
//----------------------------------------------------------------------------

#ifndef EDR_IQMCDIFFUSIBLEASSETBASE_HPP
#define EDR_IQMCDIFFUSIBLEASSETBASE_HPP

#include "edginc/VirtualDestructorBase.hpp"
//#include "edginc/DoubleMatrix.hpp"
#include "edginc/DateTime.hpp"

#include "edginc/DECLARE.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

#include "edginc/TemplateIdx.hpp"
#include "edginc/IQMCHelperTimeLogic.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/QMCStrata.hpp"

DRLIB_BEGIN_NAMESPACE


class IQMCHelperBoundedDiffusion;
FORWARD_DECLARE_REF_COUNT(IQMCRNGManager);


/************************************************************************
 * These classes form the interface for all the diffusible assets       *
 * they are convenient for a few generic methods acting on all the      *
 * different assets uniformly and for model-independent access to the   *
 * data the models produce                                              *
 ************************************************************************/

/** Abstract base of any generator of paths for any asset type  */
class MCARLO_DLL IQMCDiffusibleAssetBase : public virtual VirtualDestructorBase
{
public:
    virtual ~IQMCDiffusibleAssetBase() {}

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;

    /** get the relative probability of the given strata. Default is 1.0 */
    virtual double getStrataProbability(QMCStrataConstSP /*strata*/)
    { return 1.0; }

    /** get all the asset-specific dates */
    virtual DateTimeArray getAssetDates() = 0;

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;

    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;

    /** getting information about asset specific diffusion bounds */
    /* Note that we require to return a concrete class here to enforce consistency of max. diffusion propagation across all assets */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound()  = 0;
};
DECLARE(IQMCDiffusibleAssetBase);


/** this class provides partial implementation for the most common functionality --
    keeping the list of dates for which the diffused values were requested */
class MCARLO_DLL IQMCDiffusibleAsset : public IQMCDiffusibleAssetBase
{

public:
    IQMCDiffusibleAsset();

    /** Sophisticated clients may override the default (fast/memory inefficient) timeLogic impl. */
    virtual IQMCHelperTimeLogicSP   getTimeLogic()  const;
    /** Default implementation via timeLogic */
    virtual DateTimeArray    getAssetDates();


    /// addSpotDates is a simple redirect to addAggregatedDates; no need to be  virtual
    void    addSpotDates(const DateTimeArray& measurementDates);

    /// addForwardDates  is a simple redirect to addAggregatedDates; no need to be  virtual
    void    addForwardDates(const DateTime&      measurementDate,
                            const DateTimeArray& futureDates);

    /// addAggregatedDates passes all dates to the timeLogic class.
    /// Derived classes may redefine this function if they need to inform other assets about their dates (see CR example)
    virtual void    addAggregatedDates(
                const DateTimeArray& spot,
                const DateTimeArray& fwd,
                const DateTimeArray& fwdfwd);

    /// Same function using SPs. Only AggregatedSV will use this function and will override corresponding classes
    virtual void    addAggregatedDates(
                DateTimeArrayConstSP spot,
                DateTimeArrayConstSP fwd,
                DateTimeArrayConstSP fwdfwd,
                const DateTime& maxDiffDate,
                const DateTime& maxCurveMaturity);
    
    /// Simple way to update maxMaturity information by looking what is passed to the asset.
    /// No need to make it virtual
    void updateMaxMaturity(
                const DateTimeArray& spot,
                const DateTimeArray& fwd,
                const DateTimeArray& fwdfwd );

    /// Return requested Spot dates trimmed by maxMaturity calculated from SVGens posted for this and dependent assets
    DateTimeArray getSpotDates();

    SpotIdx getSpotIndex( const DateTime& measurementDate ) const; // returns index inside getSpotDates[]

    /// Return requested ForwardDates trimmed by maxMaturity
    DateTimeArray getForwardDates();

    SpotIdx getForwardIndex( const DateTime& forwardDate ) const;

    /// Reqturn requested ForwardForward dates trimmed by maxCurveMaturity
    DateTimeArray getForwardForwardDates();

    FwdIdx getForwardForwardIndex( const DateTime& forwardDate ) const;

    virtual size_t  getNumESDFDates(); // classes will override for eff.

private:
    IQMCHelperTimeLogicSP timeLogic; ///< default impl. of timeLogic; sophisticated clients should override getTimeLogic()
};

/** a hookup for across-the-assets diffusion collector */
DECLARE(IQMCDiffusibleAsset);


DRLIB_END_NAMESPACE
#endif //EDR_IQMCDIFFUSIBLEASSETBASE_HPP

