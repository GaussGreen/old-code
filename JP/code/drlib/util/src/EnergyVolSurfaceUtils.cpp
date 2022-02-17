// ----------------------------------------------------------------------------
//
//   Group        : GCCG Derivatives Research
//
//   Filename     : EnergyVolSurfaceUtils.cpp
//
//   Description  : Modified after drcommodityvolsurfaceutils.cpp
//
//   Author       : Sean Chen
//
//   Date         : July 28, 2005
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/EnergyVolSurfaceUtils.hpp"
#include "edginc/Interpolator.hpp"
//#include "edginc/DRProbability.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

#include <math.h>

double EnergyVolSurfaceUtils::ConvertDeltaToStrike(
    double delta, 
    double time, 
    double fwd, 
    double vol)
{
    
	return fwd * exp(0.5*vol*vol*time-vol*sqrt(time)*N1Inverse(delta));
}

double EnergyVolSurfaceUtils::ConvertInvertedDeltaToStrike(
    double invertedDelta, 
    double time, 
    double fwd, 
    double vol)
{
    return fwd * exp(0.5*vol*vol*time-vol*sqrt(time)*invertedDelta);
}

double EnergyVolSurfaceUtils::ConvertStrikeToDelta(
    double strike, 
    double time, 
    double fwd, 
    double vol)
{
    //return DRProbability::CN((log(fwd/strike)+0.5*vol*vol*time)/(vol*sqrt(time)));
	return N1((log(fwd/strike)+0.5*vol*vol*time)/(vol*sqrt(time)));
}

double EnergyVolSurfaceUtils::ConvertDeltaToStrike(
    double delta, 
    double time, 
    double fwd, 
    const Interpolator::InterpolantConstSP& volCurveByStrike)
{

    const int maxIterations = 50;
    double strikeBound1 = 0.0;
    double strikeBound2 = 5.0*fwd;
    double strike = (strikeBound2+strikeBound1)/2.0;
    int iteration;

    double error = delta-ConvertStrikeToDelta(
             strike, time, fwd, volCurveByStrike->value(strike));
    if (error<0.0)
        strikeBound1 = strike;
    else
        strikeBound2 = strike;
    for (iteration=0; iteration<maxIterations && fabs(error)>1e-6; ++iteration)
    {
        strike = (strikeBound2+strikeBound1)/2.0;
        error = delta-ConvertStrikeToDelta(
            strike, time, fwd, volCurveByStrike->value(strike));
        if (error<0.0)
            strikeBound1 = strike;
        else
            strikeBound2 = strike;
    }
    
    if (iteration==maxIterations)
        throw ModelException("EnergyVolSurfaceUtils::ConvertDeltaToStrike",
                             "Cannot convert delta to strike");

    return strike;
}

double EnergyVolSurfaceUtils::ConvertStrikeToDelta(
    double strike, 
    double time, 
    double fwd, 
    const Interpolator::InterpolantConstSP& volCurveByDelta)
{
    const int maxIterations = 50;
	const double maxError = fwd*1e-6;
    double deltaBound1 = 1.0;
    double deltaBound2 = 0.0;
    double delta = (deltaBound1+deltaBound2)/2.0;

    double error = strike-ConvertDeltaToStrike(
        delta, time, fwd, volCurveByDelta->value(delta));
    int iteration;

    if (error<0.0)
        deltaBound2 = delta;
    else
        deltaBound1 = delta;
    for (iteration=0; iteration<maxIterations && deltaBound1>0.001 && deltaBound2<0.999 && fabs(error)>maxError; ++iteration)
    {
        delta = (deltaBound1+deltaBound2)/2.0;
        error = strike-ConvertDeltaToStrike(delta, time, fwd,
            volCurveByDelta->value(delta));
        if (error<0.0)
            deltaBound2 = delta;
        else
            deltaBound1 = delta;
    }
    
    if (iteration==maxIterations)
        throw ModelException("EnergyVolSurfaceUtils::ConvertDeltaToStrike",
                             "Cannot convert strike to delta");

    return delta;
}


DRLIB_END_NAMESPACE
