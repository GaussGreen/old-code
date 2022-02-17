/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curvetypedef.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_CURVETYPEDEF_H
#define _INGPBASE_CURVETYPEDEF_H

/// this headers has to come first
/// this header comes firts as it includes some preprocessor constants!
#include "removeidentifiedwarning.h"

#include "port.h"
#include "gplinalgtypedef.h"


CC_BEGIN_NAMESPACE( ARM )


/// forward declaration
template <typename T = double, typename U = T>							class ARM_T_Curve; 
template <typename T = double, typename U = T>							class ARM_T_FlatCurve; 
template <typename T, typename U=T>										struct ARM_Interpolator;
template <typename T = double, typename U=T>							struct ARM_LinearInterpolatorCstExtrapol;
template <typename T = double, typename U=T>							struct ARM_StepUpRightOpenCstExtrapol;
template <typename T = double, typename U=T>							struct ARM_StepUpLeftOpenCstExtrapol;
template <typename T = double, typename U=T>							struct ARM_LinearInterpolatorMidCstExtrapol;

/// simple curves
typedef ARM_T_Curve<double,double>								ARM_Curve;
typedef ARM_T_FlatCurve<double,double>							ARM_FlatCurve;
typedef ARM_Interpolator<double, double >						ARM_CurveInterpolator;
typedef ARM_LinearInterpolatorCstExtrapol<double,double>		ARM_LinInterpCstExtrapolDble;
typedef ARM_StepUpRightOpenCstExtrapol<double,double>			ARM_StepUpRightOpenCstExtrapolDble;
typedef ARM_StepUpLeftOpenCstExtrapol<double,double>			ARM_StepUpLeftOpenCstExtrapolDble;
typedef ARM_LinearInterpolatorMidCstExtrapol<double,double>		ARM_LinInterpMidCstExtrapolDble;
	
/// multi curves
typedef ARM_T_Curve<double, ARM_GP_T_Vector<double> >							ARM_MultiCurve;
typedef ARM_Interpolator<double, ARM_GP_T_Vector<double> >						ARM_MultiCurveInterpolator;
typedef ARM_LinearInterpolatorCstExtrapol<double, ARM_GP_T_Vector<double> >		ARM_LinInterpCstExtrapolVec;
typedef ARM_StepUpRightOpenCstExtrapol<double, ARM_GP_T_Vector<double> >		ARM_StepUpRightOpenCstExtrapolVec;
typedef ARM_StepUpLeftOpenCstExtrapol<double, ARM_GP_T_Vector<double> >			ARM_StepUpLeftOpenCstExtrapolVec;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
