//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StochasticYieldCurve.hpp
//
//   Description : Basically a traditional YieldCurve with Vol
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STOCHASTICYIELDCURVE_HPP
#define EDR_STOCHASTICYIELDCURVE_HPP
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
class CVolRequest;
class CAsset;

/** Basically a traditional YieldCurve with Vol */
class MARKET_DLL IStochasticYieldCurve: public virtual IYieldCurve{
public:
    static CClassConstSP const TYPE; 

    virtual ~IStochasticYieldCurve();

    // perhaps a more specialized verion of getProcessedVol?
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IStochasticYieldCurve> IStochasticYieldCurveSP;
typedef smartConstPtr<IStochasticYieldCurve> IStochasticYieldCurveConstSP;

DRLIB_END_NAMESPACE
#endif
