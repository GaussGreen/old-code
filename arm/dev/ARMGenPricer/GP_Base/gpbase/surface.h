/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file surface.h
 *  \brief file for the definition of templated surfaces
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */

#ifndef _INGPBASE_SURFACE_H
#define _INGPBASE_SURFACE_H

#include "port.h"
#include "env.h"			/// to have strict validation in debug mode
#include "interpolator2D.h"
#include "rootobject.h"
#include "comparisonfunctor.h"
#include "numericconstant.h"
#include "checkinputs.h"
#include "assignop.h"
#include "typedef.h"

#include <cstdio>			/// for fclose
#include <algorithm>
#include <iomanip>

#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T1 = double, typename T2 = T1, typename T3 = T2 > class ARM_T_SurfaceWithInterpol;

////////////////////////////////////////////
/// \brief Surface BASE CLASS
///		ARM_T_SurfaceBase is a template base class for all the 
///			surface objects. The key operator is 
///			virtual T3 operator( T1 x1, T2 x2 ) const =0 defined as
///			pure virtual
////////////////////////////////////////////

template <typename T1 = double, typename T2 = T1, typename T3 = T2 >
	class ARM_T_SurfaceBase : public ARM_RootObject
{
private:
        T3 itsDefaultValue;
public:
	inline T3 operator()( T1 x1, T2 x2 ) const { return Interpolate(x1,x2); }
    inline T3 GetDefaultValue() const { return itsDefaultValue;}
	virtual T3 Interpolate(T1 x1, T2 x2 ) const = 0;
	virtual T3 InterpolateAtPoint(size_t x1Pos, size_t x2Pos ) const = 0;
	virtual void insert(T1 x1, T2 x2, T3 value ) = 0;
	virtual void insertAtPoint(size_t x1Pos, size_t x2Pos, T3 value ) = 0;
	virtual void insert( const ARM_T_SurfaceBase<T1,T2,T3>* rhs ) = 0;
	virtual void reserve(size_t x1Size, size_t x2Size ) = 0;
	virtual size_t size() const = 0;
	ARM_T_SurfaceBase(T3 defaultValue = ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER) : itsDefaultValue(defaultValue),ARM_RootObject() {};
	ARM_T_SurfaceBase(const ARM_T_SurfaceBase<T1,T2,T3>& rhs ) : ARM_RootObject(rhs), itsDefaultValue(rhs.itsDefaultValue) {};
	ARM_T_SurfaceBase& operator=(const ARM_T_SurfaceBase<T1,T2,T3>& rhs )
	{
		if( this != &rhs )
        {
			ARM_RootObject::operator=(rhs);
            itsDefaultValue=rhs.itsDefaultValue;
        }
		return *this;
	}

	/// downcast to a discretized surface
	virtual ARM_T_SurfaceWithInterpol<T1,T2,T3>* AsDiscreteSurface()
	{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": makes no sense to downcast to a discretized surface!" ); }

	/// standard ARM_Object support
	virtual ~ARM_T_SurfaceBase() {}
	virtual string ExportShortName() const { return "LSURF";}

		/// accessors
	virtual const ARM_GP_T_Vector<T1>& GetX1() const=0 ;
	virtual const ARM_GP_T_Vector<T2>& GetX2() const=0 ;
	virtual const ARM_GP_T_Matrix<T3>& GetX3() const=0 ;

	virtual bool operator ==(const T3& val) const = 0;

};


///////////////////////////////////
///// Flat surface
///////////////////////////////////
template <typename T1 = double, typename T2 = T1, typename T3 = T2 >
	class ARM_T_FlatSurface : public ARM_T_SurfaceBase<T1,T2,T3>
{
public:
	/// interpolation, insert, reserve ..
	virtual T3 Interpolate(T1 x1, T2 x2 ) const { return itsX3; }
	virtual size_t size() const { return 1; }
	virtual T3 InterpolateAtPoint(size_t x1Pos, size_t x2Pos ) const { return itsX3; }
	virtual void insert(T1 x1, T2 x2, T3 value ) 
	{ 
		itsX3 = value;
		itsX2 = 0.0;
		itsX1 = 0.0; 
	}
	double GetValue() const { return itsX3.Elt(0,0); }
	virtual void insert( const ARM_T_SurfaceBase<T1,T2,T3>* surface )
	{
		if( const ARM_T_FlatSurface<T1,T2,T3>* flatSurface = dynamic_cast< const ARM_T_FlatSurface<T1,T2,T3>* >(surface) )
		{
			itsX3 = flatSurface->GetX3().Elt(0,0);
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only flat surface can be inserted!" );
	}

	virtual void insertAtPoint(size_t x1Pos, size_t x2Pos, T3 value ) 
	{
		itsX3 = value;
		itsX2 = 0.0;
		itsX1 = 0.0; 
	}
	virtual void reserve(size_t x1Size, size_t x2Size ) {};

	/// constructor, copy constructor, assignment operator and destructor
	virtual ~ARM_T_FlatSurface() {}
	ARM_T_FlatSurface( T3 value,
        T3 defaultValue = ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER )
	: ARM_T_SurfaceBase<T1,T2,T3>(defaultValue),
	itsX3(value),
	itsX2(0.0),
	itsX1(0.0) {}
	ARM_T_FlatSurface( const ARM_T_FlatSurface<T1,T2,T3>& rhs )
	:	ARM_T_SurfaceBase<T1,T2,T3>(rhs),
	itsX3(rhs.itsX3),
	itsX2(rhs.itsX2),
	itsX1(rhs.itsX1)
	{};
	
		/// assignment operator
	// ASSIGN_OPERATOR(ARM_T_SurfaceBase)
	ASSIGN_OPERATOR(ARM_T_FlatSurface)
	
	virtual ARM_Object* Clone() const { return new ARM_T_FlatSurface<T1,T2,T3>(*this); }
	virtual string toString( const string& indent="", const string& nextIndent="" ) const { return "ARM_T_FlatSurface"; }
	virtual const ARM_GP_T_Matrix<T3>& GetX3() const {return *new ARM_GP_T_Matrix<T3>(1,1,itsX3);};
	virtual const ARM_GP_T_Vector<T2>& GetX2() const {return *new ARM_GP_T_Vector<T2>(1,itsX2);};
	virtual const ARM_GP_T_Vector<T1>& GetX1() const {return *new ARM_GP_T_Vector<T1>(1,itsX1);}
	virtual bool operator ==(const T3& val) const {return itsX3 == val;}


private:
	T1 itsX1;
	T2 itsX2;
	T3 itsX3;
};

///////////////////////////////////
///// Std surface with an interpolator
///////////////////////////////////
template <typename T1 = double, typename T2 = T1, typename T3 = T1>
	class ARM_T_SurfaceWithInterpol : public ARM_T_SurfaceBase<T1,T2,T3>
{
public:
	/// constructor, copy constructor, assignment operator and destructor
	ARM_T_SurfaceWithInterpol( 
		const ARM_GP_T_Vector<T1>& x1,
		const ARM_GP_T_Vector<T2>& x2,
		const ARM_GP_T_Matrix<T3>& x3,
		//ARM_T_Interporlator2D<T1,T2,T3>*  interpolator,
		const ARM_InterpolType& type,
        T3 defaultValue = ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER)
	:	ARM_T_SurfaceBase<T1,T2,T3>(defaultValue), itsX1(x1), itsX2(x2), itsX3(x3), itsInterpolType(type), itsInterpolator( NULL)
	{
		switch(type)
		{
		case ARM_InterpolationType::linear_column_extrapoleCst:
		case ARM_InterpolationType::linear_column_row_extrapoleCst_column_row:
			itsInterpolator  = new ARM_2DLin2Interpol;
			break;
		case ARM_InterpolationType::linear_row_extrapoleCst:
		case ARM_InterpolationType::linear_row_column_extrapoleCst_row_column:
			itsInterpolator  = new ARM_2DLin1Interpol;
			break;
		default:
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,	"Interpolation Type::unknown type / unvalaible interpolation!" );
		}

	}
	
	ARM_T_SurfaceWithInterpol(const ARM_T_SurfaceWithInterpol<T1,T2,T3>& rhs )
	:	ARM_T_SurfaceBase<T1,T2,T3>(rhs), itsX1(rhs.itsX1), itsX2(rhs.itsX2), itsX3( rhs.itsX3), itsInterpolType(rhs.itsInterpolType)
	{
        itsInterpolator = rhs.itsInterpolator ? (ARM_T_Interporlator2D<T1,T2,T3>* ) rhs.itsInterpolator->Clone(): NULL;
	};

	ARM_T_SurfaceWithInterpol& operator=(const ARM_T_SurfaceWithInterpol<T1,T2,T3>& rhs )
	{
		if( this != &rhs )
		{
			this->~ARM_T_SurfaceWithInterpol();
			new (this) ARM_T_SurfaceWithInterpol<T1,T2,T3> (rhs);
		}
		return *this;
	}
	virtual ~ARM_T_SurfaceWithInterpol()
	{
		delete itsInterpolator;
		itsInterpolator = NULL;
	}

	inline ARM_T_Interporlator2D<T1,T2,T3>* GetInterpolator() const { return itsInterpolator; }
	inline ARM_InterpolType GetInterpolType() const { return itsInterpolType; }
	virtual ARM_T_SurfaceWithInterpol<T1,T2,T3>* AsDiscreteSurface() {	return this; }

	/// accessors
	const ARM_GP_T_Vector<T1>& GetX1() const { return itsX1; }
	const ARM_GP_T_Vector<T2>& GetX2() const { return itsX2; }
	const ARM_GP_T_Matrix<T3>& GetX3() const { return itsX3; }

	/// size
	virtual size_t size() const { return itsX3.size(); }
	virtual T3 InterpolateAtPoint(size_t x1Pos, size_t x2Pos ) const { return itsX3(x1Pos,x2Pos); }
	virtual void insertAtPoint(size_t x1Pos, size_t x2Pos, T3 value ) { itsX3(x1Pos,x2Pos)= value; }
	virtual void reserve(size_t x1Size, size_t x2Size )
	{
		itsX1.reserve(x1Size);
		itsX2.reserve(x2Size);
		itsX3.reserve(x1Size,x2Size);
	};

	inline void insertOneRow( int x1Pos, T1 x1 )
	{
		itsX1.insert( itsX1.begin()+x1Pos, x1);
		ARM_GP_T_Vector<T1>::iterator end = itsX3.end();
		itsX3.resize( itsX3.rows()+1,itsX3.cols());
		int i,j;

		for( i=itsX3.rows()-1; i>x1Pos; --i)
			for( j=itsX3.cols()-1; j>=0; --j )
				*(itsX3.begin()+i*itsX3.cols()+j) = *(itsX3.begin()+(i-1)*itsX3.cols()+j);
		for( j=itsX3.cols()-1; j>=0; --j )
			*(itsX3.begin()+i*itsX3.cols()+j) = GetDefaultValue();
	}

	inline void insertLastRow(int x1Pos, T1 x1 )
	{
		itsX1.push_back(x1);
		itsX3.resize( itsX3.rows()+1,itsX3.cols());
		CC_NS(std,fill)(itsX3.begin()+(x1Pos+1)*itsX3.cols(),itsX3.end(),GetDefaultValue());
	}


	inline void insertOneCol( int x2Pos, T2 x2 )
	{
		itsX2.insert( itsX2.begin()+x2Pos,x2);
		itsX3.resize( itsX3.rows(),itsX3.cols()+1);
		int i,j;

		/// the various lines
		for( i=itsX3.rows()-1; i>0; --i)
		{
			for( j=itsX3.cols()-1; j>x2Pos; --j )
				*(itsX3.begin()+i*itsX3.cols()+j) = *(itsX3.begin()+i*itsX3.cols()+j-i-1);
			*(itsX3.begin()+i*itsX3.cols()+x2Pos) = GetDefaultValue();
			for( j=x2Pos-1; j>=0; --j )
				*(itsX3.begin()+i*itsX3.cols()+j) = *(itsX3.begin()+i*itsX3.cols()+j-i);
		}
		/// last line
		for( j=itsX3.cols()-1; j>x2Pos; --j )
			*(itsX3.begin()+i*itsX3.cols()+j) = *(itsX3.begin()+i*itsX3.cols()+j-1);
		*(itsX3.begin()+x2Pos) = GetDefaultValue();
	}


	inline void insertLastCol( int x2Pos, T2 x2 )
	{
		itsX2.push_back( x2);
		itsX3.resize( itsX3.rows(),itsX3.cols()+1);
		int i,j;

		/// the various lines
		for( i=itsX3.rows()-1; i>0; --i)
		{
			*(itsX3.begin()+i*itsX3.cols()+x2Pos+1) = GetDefaultValue();
			for( j=x2Pos; j>=0; --j )
				*(itsX3.begin()+i*itsX3.cols()+j) = *(itsX3.begin()+i*itsX3.cols()+j-i);
		}
		/// last line
		*(itsX3.begin()+i*itsX3.cols()+x2Pos+1) = GetDefaultValue();
	}

	inline void insertInX1( int& x1Pos, T1 x1 )
	{
		if( x1Pos<0 )
			insertOneRow(x1Pos=0,x1);	
		else  if( x1Pos == itsX1.size()-1 && itsX1[x1Pos]+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE<x1)
			insertLastRow(x1Pos++,x1);
		else if( fabs(x1-itsX1[x1Pos])>ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
			insertOneRow(++x1Pos,x1);
	}

	inline void insertInX2( int& x2Pos, T2 x2 )
	{		
		if( x2Pos < 0 )
			insertOneCol(x2Pos=0,x2);
		else if( x2Pos == itsX2.size()-1 && itsX2[x2Pos]+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE<x2)
			insertLastCol(x2Pos++,x2);
		else if( fabs(x2-itsX2[x2Pos])>ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
			insertOneCol(++x2Pos,x2);
	}

	virtual void insert(T1 x1, T2 x2, T3 value ) 
	{
		/// range check before
		RangeCheck();
		
		int x1Pos = lower_boundPosWithPrecision(itsX1,x1)-1;
		int x2Pos = lower_boundPosWithPrecision(itsX2,x2)-1;
		insertInX1(x1Pos,x1);
		insertInX2(x2Pos,x2);
		itsX3(x1Pos,x2Pos) = value;
		
		/// range check after
		RangeCheck();
	}

	virtual void insert( const ARM_T_SurfaceBase<T1,T2,T3>* surface )
	{
		ARM_T_SurfaceWithInterpol<T1,T2,T3>* linSurface = dynamic_cast< ARM_T_SurfaceWithInterpol<T1,T2,T3>* >( const_cast<ARM_T_SurfaceBase<T1,T2,T3>*>(surface) );
		if(	linSurface != NULL &&  typeid(itsInterpolator ) == typeid(linSurface->GetInterpolator()) )
		{
			for( size_t i=0; i<linSurface->itsX1.size(); ++ i)
				for( size_t j=0; j<linSurface->itsX2.size(); ++ j )
				{
					if(		fabs(itsX1[i]-linSurface->itsX1[i]) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE 
						&&  fabs(itsX2[j]-linSurface->itsX2[j]) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
						itsX3(i,j) = linSurface->itsX3(i,j);
					else
						insert( linSurface->itsX1[i], linSurface->itsX2[j], linSurface->itsX3(i,j) );
				}
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only flat surface can be inserted!" );
	}

	/// Interpolation part
	virtual T3 Interpolate(T1 x1, T2 x2 ) const
	{
		RangeCheck();

		if ( (itsX1.size() == 0) || (itsX2.size() == 0) )
			return 0.0;
		
        double value;
		int x1Pos = lower_boundPosWithPrecision(itsX1,x1)-1;
		int x2Pos = lower_boundPosWithPrecision(itsX2,x2)-1;

        if(dynamic_cast<ARM_2DLin1Interpol*>(itsInterpolator))
        {
		    /// 1.1) left out of bound for x1
		    if(x1Pos < 0 )
			    value = interpolateOutOfBoundLeftX1( x1, x2, x1Pos, x2Pos );
		    ///	1.2) check whether we are exactly at an x1 point
		    else if( fabs(x1-itsX1[x1Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			    value = valueAtx1Point(x1,x2,x1Pos,x2Pos); 
		    /// 1.3) right out of bound for x1
		    else if( x1Pos == itsX1.size()-1)
			    value = interpolateOutOfBoundRightX1(x1, x2, x1Pos, x2Pos );
		    /// 1.4) interpolate in x1
		    else 
                value = interpolateX1(x1,x2,x1Pos,x2Pos);
        }
        else if(dynamic_cast<ARM_2DLin2Interpol*>(itsInterpolator))
        {
		    /// 1.1) left out of bound for x2
		    if(x2Pos < 0 )
			    value = interpolateOutOfBoundLeftX2( x1, x2, x1Pos, x2Pos );
		    ///	1.2) check whether we are exactly at an x2 point
		    else if( fabs(x2-itsX2[x2Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			    value = valueAtx2Point(x1,x2,x1Pos,x2Pos); 
		    /// 1.3) right out of bound for x2
		    else if( x2Pos == itsX2.size()-1)
			    value = interpolateOutOfBoundRightX2(x1, x2, x1Pos, x2Pos );
		    /// 1.4) interpolate in x2
		    else 
                value = interpolateX2(x1,x2,x1Pos,x2Pos);
        }

        return value;
	}

	inline T3 interpolateX1( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x2
		if( x2Pos < 0 )
			return itsInterpolator->InterpolateX1ExtrapolateLeftX2( itsX1, itsX2, itsX3, x1Pos, 0, x1, x2 );
		///	1.2) check whether we are exactly at an x2 point
		else if( fabs(x2-itsX2[x2Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->InterpolateX1( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 1.3) right out of bound for x2
		else if( x2Pos == itsX2.size()-1 )
			return itsInterpolator->InterpolateX1ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos, itsX2.size()-1, x1, x2 );
		/// 1.4) interpolate
		else return itsInterpolator->InterpolateX1InterpolateX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );;
	}


	inline T3 valueAtx1Point( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1 left out of bound for x2
		if( x2Pos < 0 )
			return itsInterpolator->ExtrapolateLeftX2( itsX1, itsX2, itsX3, x1Pos, 0, x1, x2 );
		///	1.2) check whether we are exactly at an x2 point
		else if( fabs(x2-itsX2[x2Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
		{
			ARM_GP_Vector* vect = itsX3.GetRow(x1Pos);
			size_t size = vect->size();
            int Ux2Pos = upper_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
			int Lx2Pos = lower_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
            delete vect;
			///First validation to check no full colomun with default value
			if(Ux2Pos == size && Lx2Pos == -1)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only linSurface without flat infinity colomun is valid!" );
			else if(Lx2Pos == x2Pos)
				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
			
			else if(Ux2Pos == size){
				if(itsInterpolType == ARM_InterpolationType::linear_row_column_extrapoleCst_row_column){
					ARM_GP_Vector* vect = itsX3.GetColumn(x2Pos);
					int Ux1Pos = upper_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
					int Lx1Pos = lower_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
					if(Ux1Pos != vect->size() && Lx1Pos != -1)
					{
						ARM_2DLin2Interpol Interpolator;
						return Interpolator.InterpolateX1WithBounds( itsX1, itsX2, itsX3, Lx1Pos,Ux1Pos, x2Pos, x1, x2 );
					}
					delete vect;
				}
				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, x1Pos, Lx2Pos, x1, x2 );
			}
			else if(Lx2Pos == -1){
				if(itsInterpolType == ARM_InterpolationType::linear_row_column_extrapoleCst_row_column){
					ARM_GP_Vector* vect = itsX3.GetColumn(x2Pos);
					int Ux1Pos = upper_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
					int Lx1Pos = lower_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
					if(Ux1Pos != vect->size() && Lx1Pos != -1)
					{
						ARM_2DLin2Interpol Interpolator;
						return Interpolator.InterpolateX1WithBounds( itsX1, itsX2, itsX3, Lx1Pos,Ux1Pos, x2Pos, x1, x2 );
					}
					delete vect;
				}
				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, x1Pos, Ux2Pos, x1, x2 );
			}
			else 
				return itsInterpolator->InterpolateX2WithBounds( itsX1, itsX2, itsX3, x1Pos, Lx2Pos, Ux2Pos, x1, x2 );
		}
		/// 1.3) right out of bound for x2
		else if( x2Pos == itsX2.size()-1 )
			return itsInterpolator->ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos, itsX2.size()-1, x1, x2 );
		/// 1.4) interpolate
		else return itsInterpolator->InterpolateX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
	}

	inline T3 interpolateOutOfBoundLeftX1( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x2
		if( x2Pos<0)
			return itsInterpolator->ExtrapolateLeftX1ExtrapolateLeftX2( itsX1, itsX2, itsX3, 0, 0, x1, x2 );
        ///	1.2) check whether we are exactly at an x2 point
		else if( fabs(x2-itsX2[x2Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->ExtrapolateLeftX1( itsX1, itsX2, itsX3, 0, x2Pos, x1, x2 );
		/// 1.3) right out of bound for x2
		else if( x2Pos == itsX2.size()-1 )
			return itsInterpolator->ExtrapolateLeftX1ExtrapolateRightX2( itsX1, itsX2, itsX3, 0, x2Pos, x1, x2 );
		/// 1.4) simple interpolation on X2 
		else 
			return itsInterpolator->ExtrapolateLeftX1InterpolateX2( itsX1, itsX2, itsX3, 0, x2Pos, x1, x2 );
	}

	inline T3 interpolateOutOfBoundRightX1( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x2
		if( x2Pos<0)
			return itsInterpolator->ExtrapolateRightX1ExtrapolateLeftX2( itsX1, itsX2, itsX3, itsX1.size()-1, 0, x1, x2 );
		/// 1.2) right out of bound for x2
		else if( fabs(x2-itsX2[x2Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->ExtrapolateRightX1( itsX1, itsX2, itsX3, itsX1.size()-1, x2Pos, x1, x2 );
		/// 1.3) exactly in an x2 point?
		else if( x2Pos == itsX2.size()-1 )
			return itsInterpolator->ExtrapolateRightX1ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 4) simple interpolation on X2 
		else 
			return itsInterpolator->ExtrapolateRightX1InterpolateX2( itsX1, itsX2, itsX3, itsX1.size()-1, x2Pos, x1, x2 );
	}




    inline T3 interpolateX2( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x1
		if( x1Pos < 0 )
			return itsInterpolator->ExtrapolateLeftX1InterpolateX2( itsX1, itsX2, itsX3, 0, x2Pos, x1, x2 );
		///	1.2) check whether we are exactly at an x1 point
		else if( fabs(x1-itsX1[x1Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->InterpolateX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 1.3) right out of bound for x1
		else if( x1Pos == itsX1.size()-1 )
			return itsInterpolator->ExtrapolateRightX1InterpolateX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 1.4) interpolate
		else return itsInterpolator->InterpolateX1InterpolateX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );;
	}


	inline T3 valueAtx2Point( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1 left out of bound for x1
		if( x1Pos < 0 )
			return itsInterpolator->ExtrapolateLeftX1( itsX1, itsX2, itsX3, 0, x2Pos, x1, x2 );
		///	1.2) check whether we are exactly at an x1 point
		else if( fabs(x1-itsX1[x1Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
		{
			ARM_GP_Vector* vect = itsX3.GetColumn(x2Pos);
			size_t size = vect->size();
            int Ux1Pos = upper_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
			int Lx1Pos = lower_boundPosWithDefault(*vect,x1Pos,ARM_NumericConstants::ARM_INFINITY);
            delete vect;
			///First validation to check no full colomun with default value
			if(Ux1Pos == size && Lx1Pos == -1)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only linSurface without flat infinity row is valid!" );
			else if(Lx1Pos == x1Pos)
				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
			else if(Ux1Pos == size){
				if(itsInterpolType == ARM_InterpolationType::linear_column_row_extrapoleCst_column_row){
					ARM_GP_Vector* vect = itsX3.GetRow(x1Pos);
					int Ux2Pos = upper_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
					int Lx2Pos = lower_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
					if(Ux2Pos != vect->size() && Lx2Pos != -1)
					{
						ARM_2DLin1Interpol Interpolator;
						return Interpolator.InterpolateX2WithBounds( itsX1, itsX2, itsX3, x1Pos,Lx2Pos, Ux2Pos, x1, x2 );

					}
					delete vect;
				}

				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, Lx1Pos, x2Pos, x1, x2 );
			}
			else if(Lx1Pos == -1){
				if(itsInterpolType == ARM_InterpolationType::linear_column_row_extrapoleCst_column_row){
					ARM_GP_Vector* vect = itsX3.GetRow(x1Pos);
					int Ux2Pos = upper_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
					int Lx2Pos = lower_boundPosWithDefault(*vect,x2Pos,ARM_NumericConstants::ARM_INFINITY);
					if(Ux2Pos != vect->size() && Lx2Pos != -1)
					{
						ARM_2DLin1Interpol Interpolator;
						return Interpolator.InterpolateX2WithBounds( itsX1, itsX2, itsX3, x1Pos,Lx2Pos, Ux2Pos, x1, x2 );

					}
					delete vect;
				}

				return itsInterpolator->GetValueAtPoint( itsX1, itsX2, itsX3, Ux1Pos, x2Pos, x1, x2 );
			}
			else 
				return itsInterpolator->InterpolateX1WithBounds( itsX1, itsX2, itsX3, Lx1Pos, Ux1Pos, x2Pos, x1, x2 );
		}
		/// 1.3) right out of bound for x1
		else if( x1Pos == itsX1.size()-1 )
			return itsInterpolator->ExtrapolateRightX1( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 1.4) interpolate
		else return itsInterpolator->InterpolateX1( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
	}

	inline T3 interpolateOutOfBoundLeftX2( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x1
		if( x1Pos<0)
			return itsInterpolator->ExtrapolateLeftX1ExtrapolateLeftX2( itsX1, itsX2, itsX3, 0, 0, x1, x2 );
        ///	1.2) check whether we are exactly at an x2 point
		else if( fabs(x1-itsX1[x1Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->ExtrapolateLeftX2( itsX1, itsX2, itsX3, x1Pos, 0, x1, x2 );
		/// 1.3) right out of bound for x1
		else if( x1Pos == itsX1.size()-1 )
			return itsInterpolator->ExtrapolateRightX1ExtrapolateLeftX2( itsX1, itsX2, itsX3, x1Pos, 0, x1, x2 );
		/// 1.4) simple interpolation on x1
		else 
			return itsInterpolator->ExtrapolateLeftX1InterpolateX2( itsX1, itsX2, itsX3, x1Pos, 0, x1, x2 );
	}

	inline T3 interpolateOutOfBoundRightX2( T1 x1, T2 x2, int x1Pos, int x2Pos ) const
	{
		/// 1.1) left out of bound for x1
		if( x1Pos<0)
			return itsInterpolator->ExtrapolateLeftX1ExtrapolateRightX2( itsX1, itsX2, itsX3, 0,itsX2.size()-1, x1, x2 );
		/// 1.2) right out of bound for x1
		else if( fabs(x1-itsX1[x1Pos])<ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
			return itsInterpolator->ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos,itsX2.size()-1, x1, x2 );
		/// 1.3) exactly in an x1 point?
		else if( x1Pos == itsX1.size()-1 )
			return itsInterpolator->ExtrapolateRightX1ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos, x2Pos, x1, x2 );
		/// 4) simple interpolation on x1 
		else 
			return itsInterpolator->InterpolateX1ExtrapolateRightX2( itsX1, itsX2, itsX3, x1Pos,itsX2.size()-1, x1, x2 );
	}


	/// standard ARM_Object support
	virtual ARM_Object* Clone() const { return new ARM_T_SurfaceWithInterpol<T1,T2,T3>(*this); }
	virtual string toString( const string& indent="", const string& nextIndent="" ) const
	{ 
		CC_Ostringstream os;
		os << " Interpolated Surface with rows = " << itsX1.size() << " cols = " << itsX2.size() << "\n";
		os << " Interpolator : " << itsInterpolator->toString() << "\n\n";
		os << CC_NS(std,setw)(13)<<CC_NS(std,left)<<"X1/X2";
        for( size_t j=0; j<itsX2.size(); ++j )
			os << CC_NS(std,setw)(13)<< CC_NS(std,left)<< itsX2[j];
		os << "\n";

		for( size_t i=0; i<itsX1.size(); ++i )
		{
            os << CC_NS(std,setw)(13)<<  CC_NS(std,fixed) << CC_NS(std,setprecision)(1) << itsX1[i];;
			for( size_t j=0; j<itsX2.size(); ++j )
            {
                if(fabs(itsX3(i,j)) >ARM_NumericConstants::ARM_INFINITY)
                    os << CC_NS(std,setw)(13)<< CC_NS(std,scientific) << CC_NS(std,setprecision)(1) << itsX3(i,j);
                else
                    os << CC_NS(std,setw)(13)<< CC_NS(std,fixed) << CC_NS(std,setprecision)(5) << itsX3(i,j);
            }
			os << "\n";
		}
		return os.str();
	}

	virtual bool operator ==(const T3& val) const
	{
		for( size_t i=0; i<itsX1.size(); ++i )
		{
			for( size_t j=0; j<itsX2.size(); ++j )
            {
				if(itsX3(i,j) != val)
					return false;
			}
		}
		return true;
	}

private:
	ARM_GP_T_Vector<T1> itsX1;
	ARM_GP_T_Vector<T2> itsX2;
	ARM_GP_T_Matrix<T3> itsX3;
	ARM_T_Interporlator2D<T1,T2,T3>* itsInterpolator;
	ARM_InterpolType itsInterpolType;

	void RangeCheck() const
	{
#if defined( __GP_STRICT_VALIDATION )
		CheckVectorStrictlyIncreasing(itsX1,"X1", "insert",__LINE__,__FILE__);
		CheckVectorStrictlyIncreasing(itsX2,"X2", "insert",__LINE__,__FILE__);
#endif
	}
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
