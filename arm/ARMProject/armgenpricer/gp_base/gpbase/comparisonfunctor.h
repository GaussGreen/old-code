/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: comparisonfunctor.h,v $
 * Revision 1.1  2004/02/02 07:52:17  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file comparisonfunctor.h
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPBASE_COMPARISONFUNCTOR_H
#define _INGPBASE_COMPARISONFUNCTOR_H

#include "port.h"
#include "functor.h"
#include <algorithm>
#include "numericconstant.h"

#include "armdef.h"

CC_BEGIN_NAMESPACE( ARM )

template<typename T>
	struct ARM_LessWithPrecision : public CC_NS( ARM_GP,BinaryFunc)<T,T,bool>
{
		ARM_LessWithPrecision( const T& precision )
			:	itsPrecision(precision)
		{}

		bool operator()(T _X, T _Y) const
		{	return _X < _Y + itsPrecision;		}
private:
	T itsPrecision;
};


/// function to return the lower bound position with a certain precision
template<typename VecType, typename ValueType >
	int lower_boundPosWithPrecision( const VecType& vec, ValueType val, ValueType precision = K_NEW_DOUBLE_TOL)
{
	return CC_NS(std,lower_bound)(vec.begin(),vec.end(),val, 
		ARM_LessWithPrecision<double>(precision) ) - vec.begin();
}

/// function to return the upper bound position with a certain precision
template<typename VecType, typename ValueType >
	int upper_boundPosWithPrecision( const VecType& vec, ValueType val, ValueType precision = K_NEW_DOUBLE_TOL)
{
	return CC_NS(std,upper_bound)(vec.begin(),vec.end(),val, 
		ARM_LessWithPrecision<double>(precision) ) - vec.begin();
}


template<typename T>
	struct ARM_LessWithDefault : public CC_NS( ARM_GP,BinaryFunc)<T,T,bool>
{
		ARM_LessWithDefault( const T& defaultvalue )
			:	itsDefaultValue(defaultvalue)
		{}

		bool operator()(T _X, T _Y) const
		{	
            return _X = _Y;
        }
private:
	T itsDefaultValue;
};


/// function to return the fist position with a certain default value
///FIX FIX To implement with STL code
template<typename VecType, typename PosType, typename DefaultType>
	int lower_boundPosWithDefault( const VecType& vec, PosType pos, DefaultType defaultValue = ARM_NumericConstants::ARM_INFINITY)
{
    while( pos> -1 && vec[pos] > defaultValue)
		--pos;
    return pos;
}
template<typename VecType, typename PosType, typename DefaultType>
	int upper_boundPosWithDefault( const VecType& vec, PosType pos, DefaultType defaultValue = ARM_NumericConstants::ARM_INFINITY)
{
    while( pos<vec.size() && vec[pos] > defaultValue)
        ++pos;
    return pos;
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/




