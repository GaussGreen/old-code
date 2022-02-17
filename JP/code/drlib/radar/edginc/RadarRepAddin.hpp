#ifndef RADAR_REP_WRAPPER_HELPER_H
#define RADAR_REP_WRAPPER_HELPER_H

#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/RadarRepDealWrapper.hpp"

DRLIB_BEGIN_NAMESPACE

class  RADAR_DLL RadarRepAddin : public CObject 
{
public: 
    static CClassConstSP const TYPE;
    RadarRepAddin() : CObject(TYPE){}
    RadarRepAddin(CClassConstSP const &type) : CObject(type) {}
    RadarRepAddin(RadarRepDealWrapperSP p);

    virtual void validatePop2Object(void);
    virtual double getValue();

private:
    RadarRepDealWrapperSP radarDealW;
    DoubleArraySP fittingVals;
    DateTime radarDate;

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new RadarRepAddin(TYPE); }
};

DECLARE(RadarRepAddin);

DRLIB_END_NAMESPACE

#endif
