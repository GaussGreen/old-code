#ifndef RADAR_REP_WRAPPER_H
#define RADAR_REP_WRAPPER_H

#include "edginc/Object.hpp"
#include "edginc/IFuncBasisWrapper.hpp"
#include "edginc/IFittingVarTransWrapper.hpp"
#include "edginc/RadarRepUtil.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IRadarRep.hpp"

DRLIB_BEGIN_NAMESPACE

class  RADAR_DLL RadarRepWrapper : public CObject 
{
public: 
    static CClassConstSP const TYPE;
    RadarRepWrapper() : CObject(TYPE){}
    RadarRepWrapper(CClassConstSP const &type) : CObject(type) {}
    //RadarRepWrapper::RadarRepWrapper(DoubleArraySP coef, IFuncBasisWrapperSP funcBasis, IFittingVarTransWrapperSP fittingTransform);
    RadarRepWrapper(RadarRepSP radar);

    virtual void validatePop2Object(void);
    //double getValue(DoubleArraySP dArr) const;
    virtual double getValue(DoubleArraySP dArr);

    virtual IObject* clone() const;

private:
    /*DoubleArraySP coefArr;
    IFuncBasisWrapperSP funcBasis;
    IFittingVarTransWrapperSP fittingTransform;*/

    RadarRepSP underlyingRadarRep;

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new RadarRepWrapper(TYPE); }
};

DECLARE( RadarRepWrapper);

DRLIB_END_NAMESPACE

#endif
