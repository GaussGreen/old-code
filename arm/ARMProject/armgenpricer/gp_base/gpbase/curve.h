/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curve.h
 *  \brief file for the definition of templated curves
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_CURVE_H
#define _INGPBASE_CURVE_H

#include "port.h"
#include "env.h"			/// to have strict validation in debug mode
#include "interpolator.h"
#include "utilityport.h"
#include "gpvector.h"
#include "checkinputs.h"
#include "rootobject.h"
#include "specialsort.h"

#include <string>
CC_USING_NS(std,string)

#include <cmath>

/// kernel headers

#include "expt.h"



CC_BEGIN_NAMESPACE( ARM )



#define K_STD_DOUBLE_TOL 1.0e-12

/// MACRO SECTION

///		-CHECK MACROS
#if defined( __GP_STRICT_VALIDATION )
	#define CHECK_RANGE(i) CheckRange(i)
	#define CHECK_SORTED_ABSCISSES() CheckSortedAbscisses()
	#define CHECK_SAME_SIZE() CheckSameSize()
#else
	#define CHECK_RANGE(i) 
	#define CHECK_SORTED_ABSCISSES() 
	#define CHECK_SAME_SIZE() 
#endif

///		-UNARYOP
#define UNARYOP_ONEVAL(op)															\
	ARM_GP_T_Vector<T>::iterator iter;												\
	for( iter=itsOrdinates.begin(); iter!=itsOrdinates.end(); ++iter )				\
		*iter op val;																\
	return *this;																	
																					
/*#define UNARYOP_CURVE(op)															\
	ARM_GP_T_Vector<T> abscisses;													\
	abscisses.reserve(itsAbscisses.size()+rhs.GetAbscisses().size());				\
	ARM_GP_T_Vector<U> ordinates;													\
	ordinates.reserve(itsOrdinates.size()+rhs.GetOrdinates().size());				\
	size_t i,j;																		\
																					\
	for( i=0,j=0; i<itsAbscisses.size(); ++i )										\
	{																				\
		while( rhs.GetAbscisses()[j] < itsAbscisses[i] )							\
		{																			\
			abscisses.push_back( rhs.GetAbscisses()[j] );							\
			ordinates.push_back( rhs.GetOrdinates()[j] );							\
			++j;																	\
		}																			\
		if( fabs( rhs.GetAbscisses()[j] - itsAbscisses[i] ) < K_STD_DOUBLE_TOL )	\
		{																			\
			abscisses.push_back( itsAbscisses[i] );									\
			ordinates.push_back( itsOrdinates[i] op	rhs.GetOrdinates()[j]);			\
			++j;																	\
		}																			\
		else																		\
		{																			\
			abscisses.push_back( itsAbscisses[i] );									\
			ordinates.push_back( itsOrdinates[i] );									\
		}																			\
	}																				\
																					\
	for( ;j<rhs.GetAbscisses().size(); ++j )										\
	{																				\
		abscisses.push_back( rhs.GetAbscisses()[j] );								\
		ordinates.push_back( rhs.GetOrdinates()[j] );								\
	}																				\
																					\
	itsAbscisses.swap(abscisses);													\
	itsOrdinates.swap(ordinates);													\
	return *this;*/

#define UNARYOP_CURVE(op)															\
	ARM_T_Curve<T,U>* newCurve = Interpolate(rhs.GetAbscisses());					\
																					\
	itsOrdinates op rhs.GetOrdinates();												\
	return *this;


/// ARM_T_Curve is a template for curve
template <typename T = double, typename U = T>
	class ARM_T_Curve : public ARM_RootObject
{
public:
	/// constructor, destructor, copy constructor, assignment operator
	explicit ARM_T_Curve( const ARM_GP_T_Vector<T>& abscisses = ARM_GP_T_Vector<T>(),
        const ARM_GP_T_Vector<U>& ordinates  = ARM_GP_T_Vector<U>(), 
        ARM_Interpolator<T,U>* interpolator = NULL, 
        bool sortAbscisses = false )
	:	
        ARM_RootObject(), 
        itsAbscisses(abscisses),
        itsOrdinates(ordinates),
        itsInterpolator( interpolator )
	{
		if(sortAbscisses)
			ARM_T_Sort<T,U>::sortTwoVectorsWithSameSize(itsAbscisses,itsOrdinates);

		CHECK_SAME_SIZE();
		CHECK_SORTED_ABSCISSES();
		
		CC_ARM_SETNAME( ARM_GENERIC_CURVE );
	};

	
	ARM_T_Curve( const ARM_T_Curve<T,U>& rhs )
	: 
        ARM_RootObject(rhs), 
        itsAbscisses(rhs.itsAbscisses), 
        itsOrdinates(rhs.itsOrdinates)
	{
        itsInterpolator = rhs.itsInterpolator ? (ARM_Interpolator<T,U>* ) rhs.itsInterpolator->Clone(): NULL;
    };
	
	ARM_T_Curve<T,U>& operator=( const ARM_T_Curve<T,U>& rhs )
	{
		if( this != &rhs )
		{
			ARM_RootObject::operator =(rhs);
			itsAbscisses	= rhs.itsAbscisses;
			itsOrdinates	= rhs.itsOrdinates;
			delete itsInterpolator;
			itsInterpolator = rhs.itsInterpolator ? (ARM_Interpolator<T,U>* ) rhs.itsInterpolator->Clone(): NULL;
		}
		return *this;
	}

	~ARM_T_Curve()
	{
		delete itsInterpolator;
		itsInterpolator = NULL;
	}

	/// accessor (all ARM_GP_T_Vector and ith element)
	inline const ARM_GP_T_Vector<T>& GetAbscisses() const { return itsAbscisses; }
	inline ARM_GP_T_Vector<T>& GetAbscisses() { return itsAbscisses; }
    inline void SetAbscisses(const ARM_GP_T_Vector<T>& abscisses) { itsAbscisses = abscisses; }
	inline const ARM_GP_T_Vector<U>& GetOrdinates() const { return itsOrdinates; }
	inline ARM_GP_T_Vector<U>& GetOrdinates() { return itsOrdinates; }

	inline const T& GetAbscisse( size_t i ) const { CHECK_RANGE(i); return itsAbscisses[i]; }
	inline T& GetAbscisse( size_t i ) { CHECK_RANGE(i); return itsAbscisses[i]; }
	inline const U& GetOrdinate( size_t i ) const { CHECK_RANGE(i); return itsOrdinates[i]; }
	inline U& GetOrdinate( size_t i ) { CHECK_RANGE(i); return itsOrdinates[i]; }

    inline void SetOrdinates(const ARM_GP_T_Vector<U>& ordinates) { itsOrdinates = ordinates; }
    inline void SetInterpolator( ARM_Interpolator<T,U>* interpolator ){ itsInterpolator = interpolator;}
	inline ARM_Interpolator<T,U>* GetInterpolator() const { return itsInterpolator; }

	/// operator +=,-=,*=,/= on T
	ARM_T_Curve<T,U>& operator += (T val){ UNARYOP_ONEVAL(+=); }
	ARM_T_Curve<T,U>& operator -= (T val){ UNARYOP_ONEVAL(-=); }
	ARM_T_Curve<T,U>& operator *= (T val){ UNARYOP_ONEVAL(*=); }
	ARM_T_Curve<T,U>& operator /= (T val){ UNARYOP_ONEVAL(/=); }

	/// operator +=,-=,*=,/= on ARM_T_Curve<T,U>
	ARM_T_Curve<T,U>& operator += (const ARM_T_Curve<T,U>& rhs){ UNARYOP_CURVE(+=); }
	ARM_T_Curve<T,U>& operator -= (const ARM_T_Curve<T,U>& rhs){ UNARYOP_CURVE(-=); }
	ARM_T_Curve<T,U>& operator *= (const ARM_T_Curve<T,U>& rhs){ UNARYOP_CURVE(*=); }
	ARM_T_Curve<T,U>& operator /= (const ARM_T_Curve<T,U>& rhs){ UNARYOP_CURVE(/=); }
	
	/// unary minus
	ARM_T_Curve<T,U>& operator-(){	return operator *=(-1);	}
	ARM_T_Curve<T,U> operator-() const{ return -ARM_T_Curve<T,U>(*this); }

	/// Function and interpolate function
	U operator()( T val ) const { return Interpolate(val); }
	U Interpolate( T val ) const
	{	
		return itsInterpolator->Interpolate(itsAbscisses,itsOrdinates,val); 
	}

	/// ARM_GP_T_Vector like function
	size_t lower_bound( T abscisse ) const
	{ return CC_NS(std,lower_bound)(itsAbscisses.begin(), itsAbscisses.end(), abscisse) - itsAbscisses.begin();}

	inline U& operator[](size_t i) { CHECK_RANGE(i); return itsOrdinates[i]; }
	inline U operator[](size_t i) const{ CHECK_RANGE(i); return itsOrdinates[i]; }
	inline size_t size() const { return itsAbscisses.size(); }
	inline bool empty() const { return itsAbscisses.empty() || itsOrdinates.empty();}
    inline void clear() {itsAbscisses.clear(); itsOrdinates.clear();}
	size_t insert( T abscisse, U ordinate );
	void insert( const ARM_T_Curve<T,U>* rhs )
	{
	    for( size_t i=0; i<rhs->size(); ++i)
			insert( rhs->GetAbscisse(i), rhs->GetOrdinate(i) );
	}

	CC_NS(std,pair)<bool,size_t> contains( T abscisse ) const;

    double ExactSearch(T abscisse ) const;

	void push_back( T lastAbscisse, U lastOrdinate )
	{
		itsAbscisses.push_back(lastAbscisse );
		itsOrdinates.push_back(lastOrdinate );
	}

	ARM_T_Curve<T,U>* CptCurve( const ARM_GP_T_Vector<T>& abscisses ) const
	{
		ARM_GP_T_Vector<U> ordinates( abscisses.size() );
		for(size_t i=0;i<abscisses.size(); ++i)
			ordinates[i] = Interpolate(abscisses[i]);
		return new ARM_T_Curve<T,U>( abscisses, ordinates, itsInterpolator->Clone() );
	}

	ARM_T_Curve<T,U>* Interpolate( const ARM_GP_T_Vector<T>& abscisses ) const
	{
		return CptCurve(abscisses);
	}

	/// standard ARM_Object support
	virtual ARM_Object* Clone() const { return new ARM_T_Curve<T,U>(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	void sort() { ARM_T_Sort<T>::sortTwoVectorsWithSameSize(itsAbscisses,itsOrdinates); }

// FIXMEFRED: mig.vc8 (21/05/2007 10:54:04): typedef typename
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<T>::iterator t_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<T>::const_iterator t_const_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<U>::iterator u_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<U>::const_iterator u_const_iterator;
#else					// Visual C++ 6
	typedef CC_NS(ARM,ARM_GP_T_Vector)<T>::iterator t_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<T>::const_iterator t_const_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<U>::iterator u_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<U>::const_iterator u_const_iterator;

#endif 

private:
	ARM_GP_T_Vector<T> itsAbscisses;
	ARM_GP_T_Vector<U> itsOrdinates;
	ARM_Interpolator<T,U>* itsInterpolator;

	void CheckRange(size_t i) const;
	void CheckSameSize() const;
	void CheckSortedAbscisses() const;
};

// ARM_T_Curve is a template for curve
template <typename T = double, typename U = T>
	class ARM_T_FlatCurve : public ARM_T_Curve<T,U>
{
public:
     ///Constructor, destructor,copy constructor, assignmet operator
	explicit ARM_T_FlatCurve( const U value)
	:	ARM_T_Curve<T,U>(ARM_GP_T_Vector<T>(1,0.0),ARM_GP_T_Vector<U>(1,value),new ARM_StepUpRightOpenCstExtrapol<T,U>)
	{}
	ARM_T_FlatCurve ( const ARM_T_FlatCurve<T,U>& rhs)
	:ARM_T_Curve<T,U>(rhs)
	{}

	ARM_T_FlatCurve<T,U>& operator = ( const ARM_T_FlatCurve<T,U>& rhs)
	{
		if(this != &rhs)
			ARM_T_Curve ::operator = (rhs);
		return *this;
	}
	~ARM_T_FlatCurve() {}

	/// standard ARM_Object support
	virtual ARM_Object* Clone() const { return new ARM_T_FlatCurve<T,U>(*this); }
	virtual string toString(const string& indent="", const string& nextIndent="") const {return "ARM_T_FlatCurve";}
};

/// symmetric operator
/// operator +
template <typename T, typename U> inline ARM_T_Curve<T,U> operator +(const ARM_T_Curve<T,U>& lhs, T rhs )
{	return ARM_T_Curve<T,U>(lhs) += rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator +(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	
	return ARM_T_Curve<T,U>( lhs.GetAbscisses(), lhs.GetOrdinates() ) += rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator +(T lhs, const ARM_T_Curve<T,U>& rhs )
{	return ARM_T_Curve<T,U>( rhs.GetAbscisses(), ARM_GP_T_Vector<T>(rhs.size(), lhs ) ) +=rhs; }


/// operator -
template <typename T, typename U> inline ARM_T_Curve<T,U> operator -(const ARM_T_Curve<T,U>& lhs, T rhs )
{	return ARM_T_Curve<T,U>(lhs) -= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator -(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	
	return ARM_T_Curve<T,U>( lhs.GetAbscisses(), lhs.GetOrdinates() ) -= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator -(T lhs, const ARM_T_Curve<T,U>& rhs )
{	return ARM_T_Curve<T,U>( rhs.GetAbscisses(), ARM_GP_T_Vector<T>(rhs.size(), lhs ) ) -=rhs; }


/// operator *
template <typename T, typename U> inline ARM_T_Curve<T,U> operator *(const ARM_T_Curve<T,U>& lhs, T rhs )
{	return ARM_T_Curve<T,U>(lhs) *= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator *(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	
	return ARM_T_Curve<T,U>( lhs.GetAbscisses(), lhs.GetOrdinates() ) *= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator *(T lhs, const ARM_T_Curve<T,U>& rhs )
{	return ARM_T_Curve<T,U>( rhs.GetAbscisses(), ARM_GP_T_Vector<T>(rhs.size(), lhs ) ) *=rhs; }


/// operator /
template <typename T, typename U> inline ARM_T_Curve<T,U> operator /(const ARM_T_Curve<T,U>& lhs, T rhs )
{	return ARM_T_Curve<T,U>(lhs) /= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator /(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	
	return ARM_T_Curve<T,U>( lhs.GetAbscisses(), lhs.GetOrdinates() ) /= rhs; }

template <typename T, typename U> inline ARM_T_Curve<T,U> operator /(T lhs, const ARM_T_Curve<T,U>& rhs )
{	return ARM_T_Curve<T,U>( rhs.GetAbscisses(), ARM_GP_T_Vector<T>(rhs.size(), lhs ) ) /=rhs; }


/// operator ==
template <typename T, typename U> inline bool operator ==(const ARM_T_Curve<T,U>& lhs, U val )
{
    for( ARM_T_Curve<T,U>::u_const_iterator iter = lhs.GetOrdinates().begin(); iter != lhs.GetOrdinates().end(); ++iter )
		if( *iter != val )
			return false;
	return true;	
}

template <typename T, typename U> inline bool operator ==(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	
	return lhs.GetOrdinates() == rhs.GetOrdinates() &&
		lhs.GetAbscisses() == rhs.GetAbscisses();
}

template <typename T, typename U> inline bool operator ==(T val, const ARM_T_Curve<T,U>& rhs )
{	return rhs==val; }

/// operator !=
template <typename T, typename U> inline bool operator !=(const ARM_T_Curve<T,U>& lhs, T val )
{	return !(lhs==val); }

template <typename T, typename U> inline bool operator !=(const ARM_T_Curve<T,U>& lhs, const ARM_T_Curve<T,U>& rhs )
{	return !(lhs==rhs); }

template <typename T, typename U> inline bool operator !=(T val, const ARM_T_Curve<T,U>& rhs )
{	return !(val==rhs); }

/// operator <
template <typename T, typename U> inline bool operator <(T val, const ARM_T_Curve<T,U>& rhs )
{
	std::vector<double>::const_iterator ordIter;
	for(ordIter = rhs->GetOrdinates().begin(); ordIter != rhs->GetOrdinates().end(); ++ordIter)
		if (val >= *ordIter)
			return false;
	return true;
}
template <typename T, typename U> inline bool operator <(const ARM_T_Curve<T,U>& lhs, T val)
{
	std::vector<double>::const_iterator ordIter;
	for(ordIter = lhs.GetOrdinates().begin(); ordIter != lhs.GetOrdinates().end(); ++ordIter)
		if (*ordIter >= val)
			return false;
	return true;
}

/// operator <=
template <typename T, typename U> inline bool operator <=(T val, const ARM_T_Curve<T,U>& rhs )
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 10:30:29): it's not double but U
	std::vector<double>::const_iterator ordIter;
	for(ordIter = rhs.GetOrdinates().begin(); ordIter != rhs.GetOrdinates().end(); ++ordIter)
		if (val > *ordIter)
			return false;
	return true;
}
template <typename T, typename U> inline bool operator <=(const ARM_T_Curve<T,U>& lhs, T val)
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 10:30:29): it's not double but T
	std::vector<double>::const_iterator ordIter;
	for(ordIter = lhs.GetOrdinates().begin(); ordIter != lhs.GetOrdinates().end(); ++ordIter)
		if (*ordIter > val)
			return false;
	return true;
}


/// operator >
template <typename T, typename U> inline bool operator >(T val, const ARM_T_Curve<T,U>& rhs )
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 10:30:29): it's not double but U
	std::vector<double>::const_iterator ordIter;
	for(ordIter = rhs.GetOrdinates().begin(); ordIter != rhs.GetOrdinates().end(); ++ordIter)
		if (val <= *ordIter)
			return false;
	return true;
}
template <typename T, typename U> inline bool operator >(const ARM_T_Curve<T,U>& lhs, T val)
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 10:30:29): it's not double but T
	std::vector<double>::const_iterator ordIter;
	for(ordIter = lhs.GetOrdinates().begin(); ordIter != lhs.GetOrdinates().end(); ++ordIter)
		if (*ordIter <= val)
			return false;
	return true;
}


/// operator >=
template <typename T, typename U> inline bool operator >=(T val, const ARM_T_Curve<T,U>& rhs )
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 10:30:29): it's not double but U
	std::vector<double>::const_iterator ordIter;
	for(ordIter = rhs.GetOrdinates().begin(); ordIter != rhs.GetOrdinates().end(); ++ordIter)
		if (val < *ordIter)
			return false;
	return true;
}
template <typename T, typename U> inline bool operator >=(const ARM_T_Curve<T,U>& lhs, T val)
{
// FIXMEFRED: mig.vc8 bug (29/03/2007 11:28:41): it's not double but T
	std::vector<double>::const_iterator ordIter;
	for(ordIter = lhs.GetOrdinates.begin(); ordIter != lhs.GetOrdinates().end(); ++ordIter)
		if (*ordIter < val)
			return false;
	return true;
}

////////////////////////////////////////////////////
/// code part
////////////////////////////////////////////////////
/// function to paliate the pb of ostringstream not working in debug mode
int Sprintf( char* msg, double value );


/// function to paliate the pb of ostringstream not working in debug mode
int Sprintf( char* msg, const std::vector<double>* vec );

/// toString
template <typename T, typename U> string ARM_T_Curve<T,U>::toString( const string& indent, const string& nextIndent ) const
{
	char msg[4000];
	/*int pos = 0;
	pos  = sprintf( msg+pos, "%s ARM_T_Curve\n",		indent.c_str() ); 
	pos += sprintf( msg+pos, "%s Interpolator = %s\n",	indent.c_str(), itsInterpolator->toString().c_str()  ); 
	pos += sprintf( msg+pos, "%s Size         = %i\n",	indent.c_str(), itsAbscisses.size() ); 

	
	if( !itsAbscisses.empty() )
	{
		t_const_iterator absIter;
		u_const_iterator ordIter;
		pos += sprintf( msg+pos, "%s Abscissess\t Ordinates\t\n",  indent.c_str() );
		
		for( absIter = itsAbscisses.begin(), ordIter = itsOrdinates.begin();
			absIter != itsAbscisses.end(); ++absIter, ++ordIter )
		{
			pos += sprintf( msg+pos, "%s %f\t", indent.c_str(), *absIter );
			pos += Sprintf( msg+pos, *ordIter );
			pos += sprintf( msg+pos, "\n" );
		}
	}*/
	return string(msg);
}

template <typename T, typename U> ostream& operator<< ( ostream& os, const ARM_T_Curve<T,U>& crv )
{
	os << crv.toString();
	return os;
}

/// check routines
template <typename T, typename U> inline void ARM_T_Curve<T,U>::CheckRange(size_t i) const
{	CC_NS(ARM_Check,CheckRange)(itsAbscisses,i,"abscisses" ); }

template <typename T, typename U> inline void ARM_T_Curve<T,U>::CheckSameSize() const
{	CC_NS(ARM_Check,CheckSameArgSize)( itsAbscisses, itsOrdinates, "abscisses", "ordinates" ); };

template <typename T, typename U> inline void ARM_T_Curve<T,U>::CheckSortedAbscisses() const
{
	CheckVectorStrictlyIncreasing(itsAbscisses,"Abscisses","TCurve::CheckSortedAbsicces",__LINE__,__FILE__);
}

template <typename T, typename U>
	inline size_t ARM_T_Curve<T,U>::insert( T abscisse, U ordinate )
{
	size_t i = CC_NS(std,lower_bound)(itsAbscisses.begin(),itsAbscisses.end(), abscisse ) - itsAbscisses.begin();
	if( i < itsAbscisses.size() )
	{
		if( fabs( itsAbscisses[i]-abscisse) < K_STD_DOUBLE_TOL )
			itsOrdinates[i] = ordinate;
		else
		{
			/// shift everything by one
			itsAbscisses.resize(itsAbscisses.size()+1);
			CC_NS(ARM,ARM_GP_T_Vector)<T>::iterator tbegin = itsAbscisses.begin()+i,
				tend  = itsAbscisses.end()-1, 
				tdest = itsAbscisses.end();
			CC_NS(std,copy_backward)(tbegin,tend,tdest);

			itsOrdinates.resize(itsOrdinates.size()+1);	
			CC_NS(ARM,ARM_GP_T_Vector)<U>::iterator ubegin = itsOrdinates.begin()+i,
				uend  = itsOrdinates.end()-1,
				udest = itsOrdinates.end();
			CC_NS(std,copy_backward)(ubegin,uend,udest);

			/// insert in position i
			itsAbscisses[i] = abscisse;
			itsOrdinates[i] = ordinate;
		}
		return i;
	}
	else
	{
		push_back(abscisse, ordinate );
		return itsAbscisses.size()-1;
	}
}

template <typename T, typename U>
	inline CC_NS(std,pair)<bool,size_t> ARM_T_Curve<T,U>::contains( T abscisse ) const
{
	size_t i = CC_NS(std,lower_bound)(itsAbscisses.begin(),itsAbscisses.end() ) - itsAbscisses.begin();
	if( i < itsAbscisses.end() && fabs( itsAbscisses[i]-abscisse) < K_STD_DOUBLE_TOL )
		return CC_NS(std,pair)<bool,size_t>(true,i);
	else
		return CC_NS(std,pair)<bool,size_t>(false,0);
}

template <typename T, typename U>
inline double ARM_T_Curve<T,U>::ExactSearch(T abscisse ) const
{
    double ord = -1.0;

    size_t sz = itsAbscisses.size();

    size_t i = 0;

    while ( i < sz )
    {
        if ( fabs(itsAbscisses[i]-abscisse) <= K_STD_DOUBLE_TOL)
        {
           return(itsOrdinates[i]);
        }

        i++;
    }

    return(ord);
}



#undef CHECK_RANGE 
#undef CHECK_SORTED_ABSCISSES
#undef CHECK_SAME_SIZE
#undef UNARYOP_ONEVAL
#undef UNARYOP_CURVE
#undef K_STD_DOUBLE_TOL


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
