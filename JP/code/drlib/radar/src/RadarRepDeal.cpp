#include "edginc/config.hpp"
#include "edginc/RadarRepDeal.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/*double RadarRepDeal::getRadarRepValue(DateTime t, vector<double>& mktObs ) {
    double result = (* m_radar[t]) ( evaluteFittingVar(t, mktObs ));
    return result;
}*/

double RadarRepDeal::getRadarRepValue(const DateTime& t, const FittingArray& fvars) {
    if (m_radar.find(t) == m_radar.end()) {
        throw ModelException("RadarRep object for date:"+t.toString()+" not found");
    }
    IRadarRep * radarTimeT = m_radar[t].get();
    if (radarTimeT == NULL) {
        throw ModelException("RadarRep object for date:"+t.toString()+" found, but is null");
    }
    double result = (* radarTimeT) ( fvars);
    return result;
}

DRLIB_END_NAMESPACE
