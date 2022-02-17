/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: vectorarithmetic.h,v $
 * Revision 1.1  2004/02/06 16:45:06  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file vectorarithmetic.h
 *
 *  \brief files to do arithmetic on vector
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_VECTORARITHMETIC_H
#define _INGPBASE_VECTORARITHMETIC_H

#include "port.h"
#include <vector>
CC_USING_NS(std,vector)

#include <algorithm>
#include <functional>

CC_BEGIN_NAMESPACE( ARM )

template <typename T> vector<T> operator+( const vector<T>& lhs, const vector<T>& rhs )
{
	CC_NS(ARM_Check,CheckSameArgSize)( lhs, rhs, "lhs", "rhs" );
	vector<T> result( lhs.size() );
	CC_NS(std,transform)(lhs.begin(),lhs.end(),rhs.begin(),result.begin(),CC_NS(std,plus)<T>() );
	return result;
}

template <typename T> vector<T> operator-( const vector<T>& lhs, const vector<T>& rhs )
{
	CC_NS(ARM_Check,CheckSameArgSize)( lhs, rhs, "lhs", "rhs" );
	vector<T> result( lhs.size() );
	CC_NS(std,transform)(lhs.begin(),lhs.end(),rhs.begin(),result.begin(),CC_NS(std,minus)<T>() );
	return result;
}

template <typename T> vector<T> operator*( const vector<T>& lhs, T val )
{
	vector<T> result( lhs.size() );
	CC_NS(std,transform)(lhs.begin(),lhs.end(),result.begin(),CC_NS(std,bind2nd)( CC_NS(std,multiplies)<T>(), val ) );
	return result;
}

template <typename T> vector<T> operator*( T val, const vector<T>& rhs )
{
	vector<T> result( rhs.size() );
	CC_NS(std,transform)(rhs.begin(),rhs.end(),result.begin(),CC_NS(std,bind2nd)( CC_NS(std,multiplies)<T>(), val ) );
	return result;
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
