/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file interpolatorvector.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_INTERPOLATORVECTOR_H
#define _INGPBASE_INTERPOLATORVECTOR_H

#include "port.h"
#include "env.h"
#include "gpvector.h"
#include "gplinalgtypedef.h"
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )


/// function for point interpolation
inline double GetValueAtPointDerivative(const ARM_GP_Vector& times, const ARM_GP_Vector& values, size_t i)
{ return (values[i+1]-values[i])/(times[i+1]-times[i]); }
inline double GetValueAtPoint(const ARM_GP_Vector& times, const ARM_GP_Vector& values, size_t i)
{ return values[i]; }

inline double ExtrapolateLeftDerivative( const ARM_GP_Vector& times, const ARM_GP_Vector& values )
//{ return (values[1]-values[0])/(times[1]-times[0]); }
{ return 0.0; }
inline double ExtrapolateLeft( const ARM_GP_Vector& times, const ARM_GP_Vector& values )
{ return values[0]; }

inline double ExtrapolateRightDerivative( const ARM_GP_Vector& times, const ARM_GP_Vector& values )
{ return 0.0; }
inline double ExtrapolateRight( const ARM_GP_Vector& times, const ARM_GP_Vector& values )
{ return values[values.size()-1]; }

inline double InterpolateBetweenPointsDerivative( const ARM_GP_Vector& times, const ARM_GP_Vector& values, size_t i, double time )
{ return (values[i]-values[i-1])/(times[i]-times[i-1]); }
inline double InterpolateBetweenPoints( const ARM_GP_Vector& times, const ARM_GP_Vector& values, size_t i, double time )
{ return values[i-1]+(time-times[i-1])*(values[i]-values[i-1])/(times[i]-times[i-1]); }

/// interpolation value and derivative
/// interpolation at point
void VectorValueAndDerivativeLinearAtPoint( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps, ARM_GP_Vector& func, ARM_GP_Vector& der, size_t i, size_t j );

/// vectorial interpolation but with normal points
void VectorValuesAndDerivativesLinear( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func, ARM_GP_Vector& der);

/// vectorial interpolation but at mid points
void VectorValuesAndDerivativesLinearMidPoints( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func, ARM_GP_Vector& der);

/// interpolation value
/// interpolation at point
void VectorValueLinearAtPoint( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps, ARM_GP_Vector& func, size_t i, size_t j );

/// vectorial interpolation 
void VectorValuesLinear( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func );

/// vectorial interpolation 
void VectorValuesLinearMidPoints( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func );


double FunctionSpecialInterpolation(double t, const ARM_GP_Vector& times, const ARM_GP_Vector& values);
double DerivativeSpecialInterpolation(double t, const ARM_GP_Vector& times, const ARM_GP_Vector& values);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
