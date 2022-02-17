/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file interpolatorvector.cpp
 *  \brief file for vectorial interpolation
 *	\author  E.Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpbase/interpolatorvector.h"

// Kernel
#include <glob/expt.h>
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

#define K_STD_DOUBLE_TOL 1.0e-12

///////////////////////////
/// interpolation at point
///////////////////////////
void VectorValueLinearAtPoint( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps, ARM_GP_Vector& func, size_t i, size_t j )
{	
	if( i < times.size() && fabs( times[i] - timeSteps[0] ) < K_STD_DOUBLE_TOL )
		func[j] = GetValueAtPoint(times,values,i);
	else
	{
		if( i == 0 )
		{
			func[j] = ExtrapolateLeft(times,values);
		}
		else if( i == times.size() )
		{
			func[j] = ExtrapolateRight(times,values);(times,values);
		}
		else 
		{
			func[j] = InterpolateBetweenPoints(times,values,i,timeSteps[j]);
		}
	}
}


///////////////////////////
/// interpolation at point
///////////////////////////
void VectorValueAndDerivativeLinearAtPoint( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps, ARM_GP_Vector& func, ARM_GP_Vector& der, size_t i, size_t j )
{	
	if( i < times.size() && fabs( times[i] - timeSteps[0] ) < K_STD_DOUBLE_TOL )
	{
		func[j] = GetValueAtPoint(times,values,i);
		der[j]  = GetValueAtPointDerivative(times,values,i);
	}
	else
	{
		if( i == 0 )
		{
			func[j] = ExtrapolateLeft(times,values);
			der[j]  = ExtrapolateLeftDerivative(times,values);
		}
		else if( i == times.size() )
		{
			func[j] = ExtrapolateRight(times,values);
			der[j]  = ExtrapolateRightDerivative(times,values);
		}
		else 
		{
			func[j] = InterpolateBetweenPoints(times,values,i,timeSteps[j]);
			der[j]  = InterpolateBetweenPointsDerivative(times,values,i,timeSteps[j]);
		}
	}
}


////////////////////////////////////
/// vectorial interpolation but with normal points
////////////////////////////////////

void VectorValuesAndDerivativesLinear( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func, ARM_GP_Vector& der)
{
#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

	if( times.size() == 1 )
	{
		func = ARM_GP_Vector(timeSteps.size(),values[0]);
		double derivativeValue = 0;
		if( timeSteps[0] )
			derivativeValue = values[0]/timeSteps[0];
		der  = ARM_GP_Vector(timeSteps.size(),derivativeValue);
	}
	else
	{
		/// find the first point
		int i = std::lower_bound(times.begin(),times.end(), timeSteps[0] ) - times.begin();
		func.resize(timeSteps.size());
		der.resize(timeSteps.size());

		VectorValueAndDerivativeLinearAtPoint(times,values,timeSteps,func,der,i,0);
		for( size_t j=1;j<timeSteps.size();++j)
		{
			while( i<times.size() && times[i] <= timeSteps[j] )
				++i;
			VectorValueAndDerivativeLinearAtPoint(times,values,timeSteps,func,der,i,j);
		}
	}
}


////////////////////////////////////
/// vectorial interpolation 
////////////////////////////////////
void VectorValuesAndDerivativesLinearMidPoints( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func, ARM_GP_Vector& der)
{
#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

	ARM_GP_Vector midPoints;
	if( times.size() == 1 )
		midPoints = times;
	else
	{
        size_t j;
        if(times[0] > 0)
        {
		    midPoints.resize(times.size());
            midPoints[0] = 0.5*times[0];
            j=1;
        }
        else
        {
            midPoints.resize(times.size()-1);
            j=0;
        }
		size_t i;
		for( i=0; i<times.size()-1; ++i,++j )
			midPoints[j] = (times[i+1]+times[i])*0.5;
	}
	VectorValuesAndDerivativesLinear(midPoints,values, timeSteps, func, der );
}


////////////////////////////////////
/// vectorial interpolation 
////////////////////////////////////
void VectorValuesLinear( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func )
{
#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

	if( times.size() == 1 )
		func = ARM_GP_Vector(timeSteps.size(),values[0]);
	else
	{
		/// find the first point
		int i = std::lower_bound(times.begin(),times.end(), timeSteps[0] ) - times.begin();
		func.resize(timeSteps.size());

		VectorValueLinearAtPoint(times,values,timeSteps,func,i,0);
		for( size_t j=1;j<timeSteps.size();++j)
		{
			while( i<times.size() && times[i] <= timeSteps[j] )
				++i;
			VectorValueLinearAtPoint(times,values,timeSteps,func,i,j);
		}
	}
}

////////////////////////////////////
/// vectorial interpolation 
////////////////////////////////////
void VectorValuesLinearMidPoints( const ARM_GP_Vector& times, const ARM_GP_Vector& values,
	const ARM_GP_Vector& timeSteps,	ARM_GP_Vector& func )
{
#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

	ARM_GP_Vector midPoints;
	if( times.size() == 1 )
		midPoints = times;
	else
	{
		midPoints.resize(times.size()-1);
		size_t i;
		for( i=0; i<times.size()-1; ++i )
			midPoints[i] = (times[i+1]+times[i])*0.5;
	}
	VectorValuesLinear(midPoints,values, timeSteps, func );
}


////////////////////////////////////
/// Point by point interpolation used by Tree3F code...
////////////////////////////////////
double DerivativeSpecialInterpolation(double t, const ARM_GP_Vector& times, const ARM_GP_Vector& values)
{
    int iFlag=0.;
    size_t uiK=0;

#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

    if(times[0]!=0.)
		ARM_THROW(ERR_INVALID_ARGUMENT, " times must begin at time asOfDate" );

    if(times.size()<3)
		ARM_THROW(ERR_INVALID_ARGUMENT, " times must contain at least 3 values" );


    if(t >= times[times.size()-2] + 0.5*(times[times.size()-1]-times[times.size()-2]) )
        return 0.;

    // No special interpolation possible
    if(times.size()==3) return values[0];

    // Initial derivative
    double dX10=0.5*(times[1]-times[0]);
    double dX20=0.5*(times[2]-times[1]);
    if(t<times[0]+dX10)
    {
        double dXX1=times[0]+dX10;
        double dXX2=times[1]+dX20;
        double dYY1=values[0];
        double dYY2=values[1];
        return (dYY1-dYY2)/(dXX1-dXX2);
    }


    while(iFlag==0 && uiK<times.size()-2)
    {
        double dX1=0.5*(times[uiK+1]-times[uiK]);
        double dX2=0.5*(times[uiK+2]-times[uiK+1]);
        if(t>=times[uiK]+dX1 && t<times[uiK+1]+dX2)
        {
            iFlag=1;
            double dXX1=times[uiK]+dX1;
            double dXX2=times[uiK+1]+dX2;
            double dYY1=values[uiK];
            double dYY2=values[uiK+1];
            return (dYY1-dYY2)/(dXX1-dXX2);
        }
        ++uiK;
    }

    ARM_THROW(ERR_INVALID_ARGUMENT, " out of bounds in derivative interpolation" );
}

double FunctionSpecialInterpolation(double t, const ARM_GP_Vector& times, const ARM_GP_Vector& values)
{

#ifdef __GP_STRICT_VALIDATION
	if( times.empty() )
		ARM_THROW(ERR_INVALID_ARGUMENT, " times cannot be empty!" );
#endif

    if(t>times[times.size()-1]) return values[values.size()-1];

    if(times[0]!=0.)
		ARM_THROW(ERR_INVALID_ARGUMENT, " times must begin at time asOfDate" );

    if(times.size()<3)
		ARM_THROW(ERR_INVALID_ARGUMENT, " times must contain at least 3 values" );


    if(values[values.size()-1]!=values[values.size()-2])
		ARM_THROW(ERR_INVALID_ARGUMENT, " two last values must be equal to force constant extrapolation" );

    // No special interpolation possible
    if(times.size()==3) return values[0];

    double a=0.;
    double b=0.;

    // Left edge
    if( t>=times[0] && t < times[1] + 0.5*(times[2]-times[1]) )
    {
        double dT01=0.5*(times[1]);
        a=DerivativeSpecialInterpolation(t,times,values);
        b=values[0]-a*dT01;
    }

    // Middle
    size_t uiK=1;
    int iFlag=0;
    while(iFlag==0 && uiK<times.size()-2)
    {
        double dX1=0.5*(times[uiK+1]-times[uiK]);
        double dX2=0.5*(times[uiK+2]-times[uiK+1]);
        if(t>=times[uiK]+dX1&&t<times[uiK+1]+dX2)
        {
            iFlag=1;
            double dT01=times[uiK]+dX1;
            a=DerivativeSpecialInterpolation(t,times,values);
            b=values[uiK]-a*dT01;
        }
        ++uiK;
    }

    // Right edge
    if( t >= times[times.size()-2] + 0.5*(times[times.size()-1]-times[times.size()-2]) )
    {
        a=0.0;
        b=values[values.size()-1];
    }

    return a*t+b;
}

#undef K_STD_DOUBLE_TOL


CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
