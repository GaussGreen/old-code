/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: world.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file .h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */


#ifndef _INGPINFRA_GRAMFUNCTORARGHELPER_H
#define _INGPINFRA_GRAMFUNCTORARGHELPER_H

#include "gpbase/port.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// various function to help to manipulate ARM_GramFctorArgVector
struct ARM_GramFctorArgVectorHelper
{
	/// constant to tell that the index has not been found
	static const int IndexNotFound;

	static int FindFirstVectorIndex( const ARM_GramFctorArgVector& arg, const vector<size_t>& position );

	static int FindFirstVectorIndex( const ARM_GramFctorArgVector& arg, const size_t* positionIteratorStart, const size_t* positionIteratorEnd );
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

