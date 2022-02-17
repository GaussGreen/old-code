//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCPureJumps.hpp
//
//   Description : An interface of internal hooks for pure jumps asset and
//                 diffused data retrieval
//
//
//----------------------------------------------------------------------------

#ifndef IQMCPureJumps_HPP
#define IQMCPureJumps_HPP

#include "edginc/IQMCDiffusibleAssetBase.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class creates interface for a generatog of several poisson processes 
    and declares methods to return jump times from each of them. 

    Design note: 
    
    even though the jump processes can be independent, it makes sense to have 
    a single multi-poisson generator for the following reasons:
        stratified sampling, e.g. "at least 1 jump happened"
        efficiency of storing the generated jump dates
  */

class MCARLO_DLL IQMCPureJumps : public virtual IQMCDiffusibleAssetBase
{
public:

    /** Default implementation for pure jump process -- no special dates */
    virtual DateTimeArray    getAssetDates() {
        return DateTimeArray();
    }

    // Common hookups for the engine

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;

    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;

    /** getting information about asset specific diffusion bounds */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound() = 0;


    /** The simulated jump time retrieval in a generic format */

    /** after the path was generated, accessing the generated list of jump dates for all the sources */
    virtual size_t      getTotalNumberOfJumps() = 0;
    virtual size_t      getNumberOfJumps(size_t sourceIdx) = 0;
    virtual DateTime    getJumpDateByIdx(size_t sourceIdx, size_t jumpIdx) = 0;
    virtual double      getJumpDateByIdxAsDouble(size_t sourceIdx, size_t jumpIdx) =0;
    virtual double      getJumpIntensity(size_t sourceIdx) = 0;

    // There are so far no individual SV associated with this Asset type,
    // so the constructor of SV is bypassed for now.

};

/** a hookup for across-the-assets diffusion collector */
DECLARE(IQMCPureJumps);


DRLIB_END_NAMESPACE
#endif //IQMCPureJumps_HPP

