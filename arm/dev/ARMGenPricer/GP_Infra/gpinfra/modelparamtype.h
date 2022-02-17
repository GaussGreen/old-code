/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelparamtype.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou, JM Prié, O. Croissant
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_MODELPARAMTYPE_H
#define _INGPINFRA_MODELPARAMTYPE_H

#include "gpbase/port.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )


struct ARM_ModelParamType
{
	/// this list is used to associate a number to a model param object 
	/// this enables us to avoid a seach for the model param.. 
	/// in the vector of model params, the model param is precisely
	/// the element with nb ParamNb
    enum ParamNb
	{
		Volatility = 0,
		MeanReversion,
		VolatilityRatio,
		MeanReversionSpread,
        Shift,
        Alpha,
		Beta,
        Correlation,
        VolOfVol,
		VolDrift,
		Dividend,
		Multiplier,
        Skew,
        InitialVol,
		LongTermVol,
		VolMeanReversion,
		JumpProba,
		JumpSize,
		JumpVol,
		QParameter,
		LNVol,
		NVol,
		QVol,
		Drift,
		Smile,
		BrownianCorrelation,
		ForwardAdjustment,
		StrikeAdjustment,
		Hump,
		BetaCorrelation,
		CrossFactor,
		ReCorrelation,
		Sigma,
		ScalingVol,
		CompoundVol,

		/// last parameter (do not move it... it is crucial to be at the last position)!
		Unknown			/// just to see that an ARM_ModelParam is unknown
	};


	enum DataNb
	{
		BreakPointTimes,
		Values,
		Tenors,
		Strikes
	};

	/// conversion from the nb to a string
	static string GetTypeString( ParamNb nb );
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


