/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curve.cpp
 *  \brief file for the definition of templated curves
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date May 2004
 */

#include "gpbase/curve.h"
#include "gpbase/curvefactory.h"
#include "gpbase/curvetypedef.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )


ARM_RootObject* ARM_CurveFactory::CreateGenericCurve(
	const vector<double>& abscisses,
	const vector<double>& ordinates,
	long rowsNb,
	long colsNb,
	bool sortAbscisses,
	const string& interpolatorName,
	bool alwaysMulti)
{
	ARM_RootObject* genericCurve = NULL;

	bool isMultiCurve = colsNb != 1 || alwaysMulti;

	if(!isMultiCurve)
	{

		if( ordinates.size() != abscisses.size() )
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, " abscisses and ordinates must have same number of rows!" );
		}

		if( interpolatorName == "LINEAR" )
		{
			genericCurve = new ARM_Curve( ARM_GP_Vector(abscisses),ARM_GP_Vector(ordinates),new ARM_LinInterpCstExtrapolDble, sortAbscisses);
		}
		else if( interpolatorName == "STEPUPRIGHT" )
		{
			genericCurve = new ARM_Curve(ARM_GP_Vector(abscisses),ARM_GP_Vector(ordinates),new ARM_StepUpRightOpenCstExtrapolDble,sortAbscisses);
		}
		else if( interpolatorName == "STEPUPLEFT" )
		{
			genericCurve = new ARM_Curve(ARM_GP_Vector(abscisses),ARM_GP_Vector(ordinates),new ARM_StepUpLeftOpenCstExtrapolDble, sortAbscisses);
		}
		else if( interpolatorName == "LINEARMID" )
		{
			genericCurve = new ARM_Curve( ARM_GP_Vector(abscisses),ARM_GP_Vector(ordinates),new ARM_LinInterpMidCstExtrapolDble(ARM_GP_Vector(abscisses)), sortAbscisses);
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, " unknown interpolator type: permitted are LINEAR, STEPUPRIGHT, STEPUPLEFT!");
		}

	}
	else
	{
		if( rowsNb != abscisses.size() )
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, " abscisses and ordinates must have same number of rows!");
		}
	
		ARM_GP_T_Vector< ARM_GP_Vector > newOrdinates( rowsNb );

		for( size_t i=0; i<rowsNb; ++i )
		{
			newOrdinates[i] = ARM_GP_Vector( colsNb);
			for(size_t j=0; j<colsNb; ++j )
				newOrdinates[i][j] = ordinates[i*colsNb+j];
		}

		if( interpolatorName == "LINEAR" )
		{
			genericCurve = new ARM_MultiCurve(ARM_GP_Vector(abscisses),newOrdinates,new ARM_LinInterpCstExtrapolVec, sortAbscisses);
		}
		else if( interpolatorName == "STEPUPRIGHT" )
		{
			genericCurve = new ARM_MultiCurve(ARM_GP_Vector(abscisses),newOrdinates,new ARM_StepUpRightOpenCstExtrapolVec, sortAbscisses);
		}
		else if( interpolatorName == "STEPUPLEFT" )
		{
			genericCurve = new ARM_MultiCurve(ARM_GP_Vector(abscisses),newOrdinates,new ARM_StepUpLeftOpenCstExtrapolVec, sortAbscisses);
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, " unknown interpolator type: permitted are LINEAR, STEPUPRIGHT, STEPUPLEFT!");
		}
	}

	return genericCurve;
}



CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

