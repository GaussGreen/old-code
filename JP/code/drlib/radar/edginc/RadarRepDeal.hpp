// Author: Kranthi Gade
// Purpose: Evaluate RadarRepresentation at given time && fitting variables
// Method: pick the correct function object (using time) and apply it to the fitting variables

#ifndef RADAR_REP_DEAL_H
#define RADAR_REP_DEAL_H

#include "edginc/DECLARE.hpp"
#include "edginc/IRadarRep.hpp"
#include "edginc/RadarRepUtil.hpp"
#include "edginc/DateTime.hpp"

#include <map>
DRLIB_BEGIN_NAMESPACE

// a class that calculates the value of a radar deal, essentially a
// vector of RadarRepAtTimeTs, one for each time
class  RADAR_DLL RadarRepDeal {
public:
    RadarRepDeal(const map<DateTime, IRadarRepSP>& radar) : m_radar(radar) {};
    double getRadarRepValue(const DateTime& t, const FittingArray& fvars);

private:
    std::map<DateTime, IRadarRepSP> m_radar;
};

DECLARE_REF_COUNT( RadarRepDeal);

DRLIB_END_NAMESPACE

#endif
