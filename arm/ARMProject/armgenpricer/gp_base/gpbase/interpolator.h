/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file interpolator.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_INTERPOLATOR_H
#define _INGPBASE_INTERPOLATOR_H

#include "port.h"
#include "env.h"
#include "checkarg.h"
#include "expt.h"
#include "gpvector.h"
#include <algorithm>


#define K_STD_DOUBLE_TOL 1.0e-12


#if defined( __GP_STRICT_VALIDATION )
	#define CHECK_NOT_EMPTY(x,y) CC_NS(ARM_Check,CheckNotEmpty)(x,y)
#else
	#define CHECK_NOT_EMPTY(x,y) 
#endif

CC_BEGIN_NAMESPACE( ARM )

/// default class for interpolation
template <typename T, typename U=T>
	struct ARM_Interpolator
{
    /// default constructor
    ARM_Interpolator<T,U>() {}
    /// destructor
    virtual ~ARM_Interpolator() {}

	/// on x=ai returns bi
	virtual U GetValueAtPoint( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i ) const
	{
#if defined( __GP_STRICT_VALIDATION )
			CC_NS(ARM_Check,CheckRange)( abscisses,i,"abscisses" );
#endif
		return ordinates[i];
	}

	/// on [a0,...an] with value [b0,...,bn]
	/// if x<= a0 y=b0
	/// if x>= an y=bn
	virtual U ExtrapolateRight( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates ) const
	        {	CHECK_NOT_EMPTY(ordinates, "ordinates"); return ordinates[ordinates.size()-1];  }
	virtual U ExtrapolateLeft( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates ) const
	        {	CHECK_NOT_EMPTY(ordinates,"ordinates"); return ordinates[0];  }
	
	/// standard interpolation
	virtual U Interpolate( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, T val ) const
	{
		size_t i = CC_NS(std,lower_bound)(abscisses.begin(),abscisses.end(), val ) - abscisses.begin();
		if( i < abscisses.size() && fabs( abscisses[i] - val ) < K_STD_DOUBLE_TOL )
			return GetValueAtPoint(abscisses,ordinates,i);
		else
		{
			if( i == 0 )
				return ExtrapolateLeft(abscisses,ordinates);
			else if( i == abscisses.size() )
				return ExtrapolateRight(abscisses,ordinates);
			else 
				return InterpolateBetweenPoints(abscisses,ordinates,i,val);
		}
	}

	virtual U InterpolateBetweenPoints( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val ) const=0;
	virtual string toString() const	= 0;
    virtual ARM_Interpolator<T,U>* Clone() const = 0;
};


/// structure to check interpolation arguments
template <typename T, typename U=T>
	struct ARM_InterpolCheck
{
	static void CheckInterpolArg(const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val )
	{
#if defined( __GP_STRICT_VALIDATION )
		CC_NS(ARM_Check,CheckSameArgSize)( abscisses, ordinates, "abscisses", "ordinates" );
		CC_NS(ARM_Check,CheckRange)( abscisses,i,"abscisses" );
		if( i== 0 )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": lower bound cannot be equal to 0!" );
		if( abscisses[i] == val )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": you should use GetValueAtPoint!" );
#endif
	}
};


/// on [ai-1,ai] with value [bi-1,bi]
/// y=bi-1+(bi-bi-1)/(ai-ai-1)*(x-ai-1)
template <typename T = double, typename U=T>
	struct ARM_LinearInterpolatorCstExtrapol : public ARM_Interpolator<T,U>
{
        
    /// default constructor for template instantiation
	ARM_LinearInterpolatorCstExtrapol<T,U>() {}
    /// destructor
    virtual ~ARM_LinearInterpolatorCstExtrapol() {}

	virtual U InterpolateBetweenPoints( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val ) const
	{	
#if defined( __GP_STRICT_VALIDATION )
		ARM_InterpolCheck<T,U>::CheckInterpolArg( abscisses, ordinates, i, val );
		if( abscisses[i] == abscisses[i-1] )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": division by zero!" );	

#endif
		return ordinates[i-1]+(ordinates[i]-ordinates[i-1])* ( (val-abscisses[i-1])/(abscisses[i]-abscisses[i-1]) );
	}

	virtual string toString() const { return "LINEAR"; }
	virtual ARM_Interpolator<T,U>* Clone() const { return new ARM_LinearInterpolatorCstExtrapol<T,U>(*this); }
};


/// on [ai-1,ai] with value [bi-1,bi]
/// y=1{ai-1=<x}*bi-1+1{x==ai}*bi ... basically [ai-1,ai[
template <typename T = double, typename U=T>
	struct ARM_StepUpLeftOpenCstExtrapol :  public ARM_Interpolator<T,U>
{
    
    /// default constructor for template instantiation
	ARM_StepUpLeftOpenCstExtrapol<T,U>() {}
    /// destructor
    virtual ~ARM_StepUpLeftOpenCstExtrapol() {}

	virtual U InterpolateBetweenPoints( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val ) const
	{
		ARM_InterpolCheck<T,U>::CheckInterpolArg( abscisses, ordinates, i, val );
		return ordinates[i-1];
	}
    
	virtual string toString() const { return "STEPUPLEFT"; }
	virtual ARM_Interpolator<T,U>* Clone() const { return new ARM_StepUpLeftOpenCstExtrapol<T,U>(*this); }
};


/// on [ai-1,ai] with value [bi-1,bi]
/// y=1{ai-1==x}*bi-1+1{x<=ai}*bi ... basically ]ai-1,ai]
template <typename T = double, typename U=T>
	struct ARM_StepUpRightOpenCstExtrapol :  public ARM_Interpolator<T,U>
{

     	
    /// default constructor for template instantiation
	ARM_StepUpRightOpenCstExtrapol<T,U>() {}
    /// destructor
    virtual ~ARM_StepUpRightOpenCstExtrapol() {}

	virtual U InterpolateBetweenPoints( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val ) const
	{
		ARM_InterpolCheck<T,U>::CheckInterpolArg( abscisses, ordinates, i, val );
		return ordinates[i];
	}
    
	virtual string toString() const { return "STEPUPRIGHT"; }
	virtual ARM_Interpolator<T,U>* Clone() const { return new ARM_StepUpRightOpenCstExtrapol<T,U>(*this); }
};


template <typename T = double, typename U=T>
	struct ARM_LinearInterpolatorMidCstExtrapol : public ARM_Interpolator<T,U>
{
    /// default constructor for template instantiation
	ARM_LinearInterpolatorMidCstExtrapol<T,U>(const ARM_GP_T_Vector<T>& abscisses) : itsMidPoints(abscisses.size()-1)
	{
		if( abscisses.size()<2)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": size<2!" );	
		for( size_t i=0; i<abscisses.size()-1; ++i )
			itsMidPoints[i] = (abscisses[i+1]+abscisses[i])*0.5;
	}

    /// destructor
    virtual ~ARM_LinearInterpolatorMidCstExtrapol() {}

	virtual string toString() const { return "LINEAR MID"; }
	virtual ARM_Interpolator<T,U>* Clone() const { return new ARM_LinearInterpolatorMidCstExtrapol<T,U>(*this); }

	virtual U Interpolate( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, T val ) const
	{
		/// if the abscisses have been modified, then reinit the midpoint
		if( itsMidPoints.size() != abscisses.size()-1 )
		{
			itsMidPoints.resize(abscisses.size()-1);
			for( size_t i=0; i<abscisses.size()-1; ++i )
				itsMidPoints[i] = (abscisses[i+1]+abscisses[i])*0.5;
		}

		size_t i = CC_NS(std,lower_bound)(itsMidPoints.begin(),itsMidPoints.end(), val ) - itsMidPoints.begin();
		if( i < itsMidPoints.size() && fabs( itsMidPoints[i] - val ) <K_STD_DOUBLE_TOL )
			return GetValueAtPoint(itsMidPoints,ordinates,i);
		else
		{
			if( i == 0 )
				return ExtrapolateLeft(itsMidPoints,ordinates);
			else if( i == itsMidPoints.size() )
				return ExtrapolateRight(itsMidPoints,ordinates);
			else 
				/// linear interpolation
				return InterpolateBetweenPoints(itsMidPoints,ordinates,i,val);
		}
	}
	virtual U InterpolateBetweenPoints( const ARM_GP_T_Vector<T>& abscisses, const ARM_GP_T_Vector<U>& ordinates, size_t i, T val ) const
	{	
		return ordinates[i-1]+(ordinates[i]-ordinates[i-1])* ( (val-abscisses[i-1])/(abscisses[i]-abscisses[i-1]) );
	};

private:
	mutable ARM_GP_T_Vector<T> itsMidPoints;
};

#undef K_STD_DOUBLE_TOL

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
