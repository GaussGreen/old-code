/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: interpolator2D.h,v $
 * Revision 1.1  2004/05/02 07:52:17  ebenhamou
 * Initial revision
 *
 */

/*! \file interpolator2D.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date Ocotber 2004
 */

#ifndef _INGPBASE_INTERPOLATOR2D_H
#define _INGPBASE_INTERPOLATOR2D_H

#include "port.h"
#include "env.h"
#include <glob/expt.h>
#include "gpvector.h"
#include "gpmatrix.h"

CC_BEGIN_NAMESPACE( ARM )


#if defined( __GP_STRICT_VALIDATION )
	#define CHECK_NOT_EMPTY(x,y) CC_NS(ARM_Check,CheckNotEmpty)(x,y)
#else
	#define CHECK_NOT_EMPTY(x,y) 
#endif

/// the x1Pos and x2Pos are lowerbound -1 indexes!
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >
	class ARM_T_Interporlator2D
{
public:
	/// 2D stuff (9 functions)
	/// Interpolate X1 Interpolate X2
	virtual T3 InterpolateX1InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const = 0;

	/// Extrapolate X1 and Interpolate X2
	virtual T3 ExtrapolateLeftX1InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return InterpolateX2( x1, x2, x3, 0, x2Pos, x1[0], x2Value );	}
	
	virtual T3 ExtrapolateRightX1InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return InterpolateX2( x1, x2, x3, x1.size()-1, x2Pos, x1[x1.size()-1], x2Value ); }


	/// Interpolate X1 and Extrapolate X2
	virtual T3 InterpolateX1ExtrapolateLeftX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return InterpolateX1( x1, x2, x3, x1Pos, 0, x1Value, x2[0] );	}

	virtual T3 InterpolateX1ExtrapolateRightX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return InterpolateX1( x1, x2, x3, x1Pos, x2.size()-1, x1Value, x2[x2.size()-1]); }

	/// Extrapolate X1 Extrapolate X2
	virtual T3 ExtrapolateLeftX1ExtrapolateLeftX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return x3( x1Pos,x2Pos); }
	
	virtual T3 ExtrapolateRightX1ExtrapolateLeftX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return x3( x1Pos,x2Pos); }
	
	virtual T3 ExtrapolateLeftX1ExtrapolateRightX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return x3( x1Pos,x2Pos); }
	
	virtual T3 ExtrapolateRightX1ExtrapolateRightX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return x3( x1Pos,x2Pos); }

	/// 1D stuff (2 functions)
	virtual T3 InterpolateX1( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
        size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const = 0;
    virtual T3 InterpolateX1WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t Lx1Pos, size_t Ux1Pos,size_t x2Pos,T1 x1Value, T2 x2Value ) const =0;

	virtual T3 InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const = 0;

    virtual T3 InterpolateX2WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t Lx1Pos, size_t Ux1Pos,size_t x2Pos,T1 x1Value, T2 x2Value ) const=0;

	virtual T3 ExtrapolateLeftX1( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{	return GetValueAtPoint( x1, x2, x3, 0, x2Pos, x1[0], x2Value ); 	}

	virtual T3 ExtrapolateRightX1( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const 
	{	return GetValueAtPoint( x1, x2, x3, x1.size()-1, x2Pos, x1[x1.size()-1], x2Value ); }

	virtual T3 ExtrapolateLeftX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{	return GetValueAtPoint( x1, x2, x3, x1Pos, 0, x1Value, x2[0] ); }

	virtual T3 ExtrapolateRightX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{	return GetValueAtPoint( x1, x2, x3, x1Pos, x2.size()-1, x1Value, x2[x2.size()-1] ); }

	/// GetValueAtPoint
	virtual T3 GetValueAtPoint( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{	return x3( x1Pos,x2Pos); }

	virtual string toString() const	= 0;
    virtual ARM_T_Interporlator2D<T1,T2,T3>* Clone() const = 0;
};


template <typename T1 = double, typename T2 = T1, typename T3 = T1 >
	class ARM_T_Interporlator2DLin1CstExtrapol : public ARM_T_Interporlator2D<T1,T2,T3>
{
public:
	/// 2D stuff
	virtual T3 InterpolateX1InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
		T3 x1lowInterpolx3Value  = InterpolateX2( x1, x2, x3, x1Pos,   x2Pos, x1[x1Pos],   x2Value );
		T3 x1highInterpolx3Value = InterpolateX2( x1, x2, x3, x1Pos+1, x2Pos, x1[x1Pos+1], x2Value );
		return x1lowInterpolx3Value + (x1Value-x1[x1Pos]) * (x1highInterpolx3Value-x1lowInterpolx3Value)/(x1[x1Pos+1]-x1[x1Pos]);
	}

	/// 1D stuff (2 functions)
	/// interpolation on X1 assuming we are at an x2 Point
	virtual T3 InterpolateX1( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x1[x1Pos] - x1Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x1 points, this is not the case!" );

		if( fabs( x2[x2Pos] - x2Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x2 point, this is not the case!" );
#endif
		return x3(x1Pos,x2Pos) + (x1Value-x1[x1Pos]) * (x3(x1Pos+1,x2Pos)-x3(x1Pos,x2Pos))/(x1[x1Pos+1]-x1[x1Pos]);
	}

	/// interpolation on X2 assuming we are at an x1 Point
	virtual T3 InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x2[x2Pos] - x2Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x2 points, this is not the case!" );

		if( fabs( x1[x1Pos] - x1Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x1 point, this is not the case!" );
#endif
		return x3(x1Pos,x2Pos) + (x2Value-x2[x2Pos]) * (x3(x1Pos,x2Pos+1)-x3(x1Pos,x2Pos))/(x2[x2Pos+1]-x2[x2Pos]);
	}
    virtual T3 InterpolateX1WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t Lx1Pos, size_t Ux1Pos,size_t x2Pos,T1 x1Value, T2 x2Value ) const
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": InterpolateX1WithBounds:: we don't need to define this function!" );
	}

    /// interpolation on X2 assuming we are at an x1 Point with upper and lower bound
	virtual T3 InterpolateX2WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t Lx2Pos,size_t Ux2Pos,T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x2[Lx2Pos] - x2Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x2 points, this is not the case!" );

		if( fabs( x1[x1Pos] - x1Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x1 point, this is not the case!" );
#endif
			return x3(x1Pos,Lx2Pos) + (x2Value-x2[Lx2Pos]) * (x3(x1Pos,Ux2Pos)-x3(x1Pos,Lx2Pos))/(x2[Ux2Pos]-x2[Lx2Pos]);
	}

	virtual string toString() const	{ return "Interporlator2DLin2CstExtrapol"; }
	virtual ARM_T_Interporlator2D<T1,T2,T3>* Clone() const { return new ARM_T_Interporlator2DLin1CstExtrapol<T1,T2,T3>(*this); }
};


template <typename T1 = double, typename T2 = T1, typename T3 = T1 >
	class ARM_T_Interporlator2DLin2CstExtrapol : public ARM_T_Interporlator2D<T1,T2,T3>
{
public:
	/// 2D stuff
	virtual T3 InterpolateX1InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
		T3 x2lowInterpolx3Value  = InterpolateX1( x1, x2, x3, x1Pos,   x2Pos, x1Value, x2[x2Pos]);
		T3 x2highInterpolx3Value = InterpolateX1( x1, x2, x3, x1Pos, x2Pos+1, x1Value, x2[x2Pos+1]);
		return x2lowInterpolx3Value + (x2Value-x2[x2Pos]) * (x2highInterpolx3Value-x2lowInterpolx3Value)/(x2[x2Pos+1]-x2[x2Pos]);
	}

	/// 1D stuff (2 functions)
	/// interpolation on X1 assuming we are at an x2 Point
	virtual T3 InterpolateX1( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x1[x1Pos]-x1Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x1 points, this is not the case!" );

		if( fabs( x2[x2Pos] - x2Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x2 point, this is not the case!" );
#endif

		return x3(x1Pos,x2Pos) + (x1Value-x1[x1Pos]) * (x3(x1Pos+1,x2Pos)-x3(x1Pos,x2Pos))/(x1[x1Pos+1]-x1[x1Pos]);
	}

    virtual T3 InterpolateX2( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t x1Pos, size_t x2Pos, T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x2[x2Pos]-x2Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x2 points, this is not the case!" );

		if( fabs( x1[x1Pos] - x1Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x1 point, this is not the case!" );
#endif

		return x3(x1Pos,x2Pos) + (x2Value-x2[x2Pos]) * (x3(x1Pos,x2Pos+1)-x3(x1Pos,x2Pos))/(x2[x2Pos+1]-x2[x2Pos]);

	}
    /// interpolation on X2 assuming we are at an x1 Point with upper and lower bound
	virtual T3 InterpolateX1WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t Lx1Pos, size_t Ux1Pos,size_t x2Pos,T1 x1Value, T2 x2Value ) const
	{
#if defined( __GP_STRICT_VALIDATION )
		if( fabs( x1[Lx1Pos] - x1Value ) < K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be between two x1 points, this is not the case!" );

		if( fabs( x2[x2Pos] - x2Value ) > K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": we expect to be at an x2 point, this is not the case!" );
#endif
			return x3(Lx1Pos,x2Pos) + (x1Value-x1[Lx1Pos]) * (x3(Ux1Pos,x2Pos)-x3(Lx1Pos,x2Pos))/(x1[Ux1Pos]-x1[Lx1Pos]);
	}
    virtual T3 InterpolateX2WithBounds( const ARM_GP_T_Vector<T1>& x1, const ARM_GP_T_Vector<T2>& x2, const ARM_GP_T_Matrix<T3>& x3,
		size_t Lx1Pos, size_t Ux1Pos,size_t x2Pos,T1 x1Value, T2 x2Value ) const
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": InterpolateX2WithBounds::we don't need to define this function!" );
	}

	virtual string toString() const	{ return "Interporlator2DLin1CstExtrapol"; }
	virtual ARM_T_Interporlator2D<T1,T2,T3>* Clone() const { return new ARM_T_Interporlator2DLin2CstExtrapol<T1,T2,T3>(*this); }
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
