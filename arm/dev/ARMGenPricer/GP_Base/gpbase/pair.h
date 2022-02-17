/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: pair.h,v $
 * Revision 1.1  2003/11/14 07:51:19  ebenhamou
 * Initial revision
 *
 *
 */


/*! \file pair.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

#ifndef _INGPBASE_PAIR_H
#define _INGPBASE_PAIR_H

#include "port.h"
#include "ostringstream.h"
#include <string>

CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )

/// very similar to std::pair<T,U>
/// but cannot even derive publicly from
/// std::pair<T,U> since this is not available
/// on unix .. AAAAAAAAAAAAAAARRRRRRRRRRRRRFGGGGGGGGG!

template <typename T, typename U>
	struct ARM_Pair		

{
	T itsElem1;
	U itsElem2;
	
	ARM_Pair(T Elem1=T(), U Elem2=U()): itsElem1(Elem1), itsElem2(Elem2) {};

	string toString( const string& indent = "" ) const
	{
		CC_Ostringstream os;
		os << "(" << itsElem1 << "," << itsElem2 << ")";
		return os.str();
	}

	bool operator==( const ARM_Pair<T,U>& rhs ) const
	{
		return itsElem1 == rhs.itsElem1
			&& itsElem2 == rhs.itsElem2;
	}

	bool operator<( const  ARM_Pair<T,U>& rhs ) const
	{
		return itsElem1 < rhs.itsElem1
			|| ( itsElem1 == rhs.itsElem1 && itsElem2 < rhs.itsElem2 );
		
	}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

