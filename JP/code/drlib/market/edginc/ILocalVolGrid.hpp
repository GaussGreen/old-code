//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ILocalVolGrid.hpp
//
//   Description : local vol grid support interpolations
//
//   Author      : Ning Shen
//
//   Date        : 6 Jan 2003
//
//
//----------------------------------------------------------------------------

#ifndef ILOCALVOLGRID_HPP
#define ILOCALVOLGRID_HPP

#include "edginc/VolProcessedDVF.hpp"

DRLIB_BEGIN_NAMESPACE

/** ILocalVolGrid interface */
class MARKET_DLL ILocalVolGrid: virtual public CVolProcessedDVF::IVolCalculator {
public:
    /** Linear interpolation of local vol from grid */ 
    virtual double interpLocVar(int step, double logSpot) const = 0;
    
    virtual double maxDrift() const = 0;
    
    /** Rolls LocalVolGrid forward */
    virtual void rollTime(const DateTime& valueDate) = 0;
};

DECLARE_REF_COUNT(ILocalVolGrid);

DRLIB_END_NAMESPACE

#endif
