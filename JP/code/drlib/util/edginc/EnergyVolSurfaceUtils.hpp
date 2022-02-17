// ----------------------------------------------------------------------------
//
//   Group        : GCCG Derivatives Research
//
//   Filename     : EnergyVolSurfaceUtils.hpp
//
//   Description  : Modified after drcommodityvolsurfaceutils.h
//
//   Author       : Sean Chen
//
//   Date         : July 28, 2005
//----------------------------------------------------------------------------

#ifndef _energyvolsurfaceutils_
#define _energyvolsurfaceutils_

#include "edginc/Interpolator.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL EnergyVolSurfaceUtils
{
	public:

		static double ConvertDeltaToStrike(
			double delta, 
			double time, 
			double fwd, 
			double vol);

		static double ConvertInvertedDeltaToStrike(
			double invertedDelta, 
			double time, 
			double fwd, 
			double vol);

		static double ConvertStrikeToDelta(
			double strike, 
			double time, 
			double fwd, 
			double vol);

		static double ConvertDeltaToStrike(
			double delta, 
			double time, 
			double fwd, 
			const Interpolator::InterpolantConstSP& volCurveByStrike);

		static double ConvertStrikeToDelta(
			double strike, 
			double time, 
			double fwd, 
			const Interpolator::InterpolantConstSP& volCurveByDelta);
		

	private:
};

DRLIB_END_NAMESPACE

#endif
