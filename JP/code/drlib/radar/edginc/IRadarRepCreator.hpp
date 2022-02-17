#ifndef _IRADAR_REP_CREATOR_HPP
#define _IRADAR_REP_CREATOR_HPP

#include "edginc/config.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IRadarRep.hpp"

DRLIB_BEGIN_NAMESPACE

class  RADAR_DLL IRadarRepCreator {
public:    
    virtual void addSample(FittingArray x, double val) = 0;
    //virtual IRadarRepSP getRadarRep() const = 0;
    virtual IRadarRepSP getRadarRep() = 0;
    virtual void clearSamples() = 0;
    virtual ~IRadarRepCreator() {}
};

DECLARE_REF_COUNT(IRadarRepCreator);

DRLIB_END_NAMESPACE

#endif
