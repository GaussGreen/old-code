#ifndef RADAR_REP_DEAL_WRAPPER_H
#define RADAR_REP_DEAL_WRAPPER_H

#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/RadarRepUtil.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/RadarRepDeal.hpp"

DRLIB_BEGIN_NAMESPACE

class RADAR_DLL RadarRepDealWrapper : public CObject 
{
public: 
    static CClassConstSP const TYPE;
    RadarRepDealWrapper() : CObject(TYPE){}
    RadarRepDealWrapper(CClassConstSP const &type) : CObject(type) {}
    RadarRepDealWrapper(RadarRepDealSP newRadarRep);

    virtual void validatePop2Object(void);
    virtual double getValue(DateTime date, DoubleArraySP dArr);

    virtual IObject* clone() const;

private:
    RadarRepDealSP underlyingRadarRep;
    
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new RadarRepDealWrapper(TYPE); }
};

DECLARE(RadarRepDealWrapper);

DRLIB_END_NAMESPACE

#endif
