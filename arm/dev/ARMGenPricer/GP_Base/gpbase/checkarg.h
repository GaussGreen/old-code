/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file checkarg.h
 *
 *  \brief header for some check
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPBASE_CHECKARG_H
#define _INGPBASE_CHECKARG_H

#include "port.h"
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM_Check )

/// Validation function that checks two vectors have the same size
template <typename T,typename U>
	void CheckSameArgSize( const T& v1, const U& v2, const char* v1Name, const char* v2Name ,
        long LINE = __LINE__, char* FILE = __FILE__)
{
	if( v1.size() != v2.size() )
	{
		char msg[255];
		sprintf( msg, "%s size = %d while %s size = %d, please advise!", v1Name, v1.size(), v2Name, v2.size() );
		throw Exception(LINE, FILE, ERR_INVALID_DATA, msg );
	}
}


/// Validation function that checks that vector have the expected size
template <typename T>
	void CheckArgSize( const T& v1, const char* v1Name, size_t expectedVectorSize )
{
	if( v1.size() != expectedVectorSize  )
	{
		char msg[255];
		sprintf( msg, "%s size = %d while expected = %d, please advise!", v1Name, v1.size(), expectedVectorSize );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
	}
}

/// Validation function that checks two vectors have the same size
template <typename T>
	void CheckRange( const T& v1, size_t i, const char* v1Name )
{
	if( i>= v1.size() )
	{
		char msg[255];
		sprintf( msg, "out of range: %s size = %d while index = %d, please advise!", v1Name, v1.size(), i );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
	}
}

template <typename T>
	void CheckNotEmpty( const T& v, const char* vName)
{
	if( v.empty())
	{
		char msg[255];
		sprintf( msg, "%s :vector is empty!",vName);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
	}
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

