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
#include "expt.h"


CC_BEGIN_NAMESPACE( ARM )


/// function for point interpolation
inline double GetValueAtPointDerivative(const std::vector<double>& times, const std::vector<double>& values, size_t i)
{ return (values[i+1]-values[i])/(times[i+1]-times[i]); }
inline double GetValueAtPoint(const std::vector<double>& times, const std::vector<double>& values, size_t i)
{ return values[i]; }

inline double ExtrapolateLeftDerivative( const std::vector<double>& times, const std::vector<double>& values )
//{ return (values[1]-values[0])/(times[1]-times[0]); }
{ return 0.0; }
inline double ExtrapolateLeft( const std::vector<double>& times, const std::vector<double>& values )
{ return values[0]; }

inline double ExtrapolateRightDerivative( const std::vector<double>& times, const std::vector<double>& values )
{ return 0.0; }
inline double ExtrapolateRight( const std::vector<double>& times, const std::vector<double>& values )
{ return values[values.size()-1]; }

inline double InterpolateBetweenPointsDerivative( const std::vector<double>& times, const std::vector<double>& values, size_t i, double time )
{ return (values[i]-values[i-1])/(times[i]-times[i-1]); }
inline double InterpolateBetweenPoints( const std::vector<double>& times, const std::vector<double>& values, size_t i, double time )
{ return values[i-1]+(time-times[i-1])*(values[i]-values[i-1])/(times[i]-times[i-1]); }

/// interpolation value and derivative
/// interpolation at point
void VectorValueAndDerivativeLinearAtPoint( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps, std::vector<double>& func, std::vector<double>& der, size_t i, size_t j );

/// vectorial interpolation but with normal points
void VectorValuesAndDerivativesLinear( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps,	std::vector<double>& func, std::vector<double>& der);

/// vectorial interpolation but at mid points
void VectorValuesAndDerivativesLinearMidPoints( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps,	std::vector<double>& func, std::vector<double>& der);

/// interpolation value
/// interpolation at point
void VectorValueLinearAtPoint( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps, std::vector<double>& func, size_t i, size_t j );

/// vectorial interpolation 
void VectorValuesLinear( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps,	std::vector<double>& func );

/// vectorial interpolation 
void VectorValuesLinearMidPoints( const std::vector<double>& times, const std::vector<double>& values,
	const std::vector<double>& timeSteps,	std::vector<double>& func );


double FunctionSpecialInterpolation(double t, const std::vector<double>& times, const std::vector<double>& values);
double DerivativeSpecialInterpolation(double t, const std::vector<double>& times, const std::vector<double>& values);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
