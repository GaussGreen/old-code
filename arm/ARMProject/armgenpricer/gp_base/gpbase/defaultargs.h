/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file defaultargs.h
 *
 *  \brief files to have default 
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPBASE_DEFAULTARGS_H
#define _INGPBASE_DEFAULTARGS_H

/// use our macro for namespace
#include "port.h"
#include "env.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Routine: CopyArgsInPlace
///	Returns: void
///	Action : copies to a type a default args vector!
////////////////////////////////////////////////////
template <typename T, typename U>
     void CopyArgsInPlace( T& output, U args, size_t argsSize)
{
	output = T(argsSize);
	for(size_t i=0; i<argsSize; ++i )
		output[i] = args[i];
}

////////////////////////////////////////////////////
///	Routine: CopyArgsPointorInPlace
///	Returns: void
///	Action : copies to a type a default args vector!
///				warning the T* reference is a very subtil
///				programming trick to avoid more template argument
///				since the type of output is T* while the new uses
///				the type T!
////////////////////////////////////////////////////
template <typename T, typename U>
     void CopyArgsPointorInPlace( T*& output, U args, size_t argsSize)
{
	delete output;
	output = new T(argsSize);
	for(size_t i=0; i<argsSize; ++i )
		(*output)[i] = args[i];
}

////////////////////////////////////////////////////
///	Routine: CopyArgsPointorInPlace
///	Returns: void
///	Action : same as CopyArgsPointorInPlace except that we do not delete the pointor before!
////////////////////////////////////////////////////
template <typename T, typename U>
     void CopyArgsPointorInPlaceNoDelete( T*& output, U args, size_t argsSize)
{
	output = new T(argsSize);
	for(size_t i=0; i<argsSize; ++i )
		(*output)[i] = args[i];
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


