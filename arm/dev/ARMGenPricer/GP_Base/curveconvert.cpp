/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *! \file curveconvert.cpp
 *
 *  \brief file to convert curve into a refvalue!
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2005
 */

#include "gpbase/curveconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"

CC_BEGIN_NAMESPACE( ARM )

// Function to convert a curve into a reference value
ARM_ReferenceValue* CurveToRefValue( const ARM_Curve& curve, double asOfDate, double base )
{
	ARM_Vector* dates =  To_pARM_Vector(const_cast<ARM_GP_Vector*>(&(curve.GetAbscisses())));
	ARM_Vector* values = To_pARM_Vector(const_cast<ARM_GP_Vector*>(&(curve.GetOrdinates())));

	(*dates) += asOfDate;
	(*values) *= base;

	int calcMethod;

	if (dynamic_cast<ARM_StepUpRightOpenCstExtrapolDble*>(curve.GetInterpolator()))
	{
		calcMethod = K_STEPUP_RIGHT;
	}
	else if (dynamic_cast<ARM_StepUpLeftOpenCstExtrapolDble*>(curve.GetInterpolator()))
	{
		calcMethod = K_STEPUP_LEFT;
	}
	else
	{
		calcMethod = K_LINEAR;
	}

	ARM_ReferenceValue* refVal =  new ARM_ReferenceValue(dates, values);
	refVal->SetCalcMethod(calcMethod);

	return refVal;
}

// Function to convert a curve into a reference value
ARM_ReferenceValue GPCurveToRefValue( const ARM_Curve& curve, double asOfDate, double base )
{
	ARM_Vector* dates =  To_pARM_Vector(const_cast<ARM_GP_Vector*>(&(curve.GetAbscisses())));
	ARM_Vector* values = To_pARM_Vector(const_cast<ARM_GP_Vector*>(&(curve.GetOrdinates())));

	(*dates) += asOfDate;
	(*values) *= base;

	int calcMethod;

	if (dynamic_cast<ARM_StepUpRightOpenCstExtrapolDble*>(curve.GetInterpolator()))
	{
		calcMethod = K_STEPUP_RIGHT;
	}
	else if (dynamic_cast<ARM_StepUpLeftOpenCstExtrapolDble*>(curve.GetInterpolator()))
	{
		calcMethod = K_STEPUP_LEFT;
	}
	else
	{
		calcMethod = K_LINEAR;
	}
	int valueType = 1;
	int conversion = 0;

	return ARM_ReferenceValue(dates, values,valueType,conversion,calcMethod);;
}

// Function to convert a reference value into a curve
ARM_Curve* RefValueToCurve( const ARM_ReferenceValue& refVal, double asOfDate, double base )
{
	ARM_GP_Vector gpDates;
	if (refVal.GetDiscreteDates())
	{
		gpDates = To_ARM_GP_Vector(refVal.GetDiscreteDates());
		// We replace the julian date by a lag from the asof dates
		gpDates -= asOfDate;
	}
	else
		gpDates = ARM_GP_Vector(1,0.0);
	ARM_GP_Vector gpValues = To_ARM_GP_Vector(refVal.GetDiscreteValues());

	gpValues /= base;

	ARM_Curve*     curve    = NULL;

	ARM_Interpolator<double,double>* interpolator = NULL;

	switch (refVal.GetCalcMethod())
	{
	case K_STEPUP_RIGHT:
			interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
		break;
	case K_STEPUP_LEFT:
			interpolator = new ARM::ARM_StepUpLeftOpenCstExtrapolDble();
		break;
	default:
		interpolator = new ARM::ARM_LinInterpCstExtrapolDble();
	}

	return new ARM_Curve(gpDates,gpValues, interpolator);	
}

// Function to convert a reference value into a curve
ARM_Curve RefValueToCurve( ARM_ReferenceValue* refVal, 
						  double asOfDate,
						  double base, 
						  int calcMethod)
{
	if(!refVal) 
		return ARM_Curve();

	ARM_GP_Vector gpDates;
	if (refVal->GetDiscreteDates())
	{
		gpDates = To_ARM_GP_Vector(refVal->GetDiscreteDates());
		// We replace the julian date by a lag from the asof dates
		gpDates -= asOfDate;
	}
	else
		gpDates = ARM_GP_Vector(1,0.0);
	ARM_GP_Vector gpValues = To_ARM_GP_Vector(refVal->GetDiscreteValues());

	gpValues /= base;

	ARM_Curve*     curve    = NULL;

	ARM_Interpolator<double,double>* interpolator = NULL;

	if (calcMethod == -1111) calcMethod = refVal->GetCalcMethod();

	switch (calcMethod)
	{
	case K_STEPUP_RIGHT:
			interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
		break;
	case K_STEPUP_LEFT:
			interpolator = new ARM::ARM_StepUpLeftOpenCstExtrapolDble();
		break;
	case K_LINEAR:
			interpolator = new ARM::ARM_LinInterpCstExtrapolDble();
		break;
	default:
		interpolator = new ARM::ARM_LinInterpCstExtrapolDble();
	}

	return  ARM_Curve(gpDates,gpValues, interpolator);	
}
CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

