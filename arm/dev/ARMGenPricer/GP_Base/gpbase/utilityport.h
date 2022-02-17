/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file utilityport.h
 *  \brief to supply for utility stl header
 * especially for max and min
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPBASE_UTILITYPORT_H
#define _INGPBASE_UTILITYPORT_H

#include "port.h"
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

/// AAAAAAAAAARRRRRRRRRRGGGGGGGGGG
/// min and max are not defined in the STL algorithm header on NT
/// could test for studio 6.0 with (_MSC_VER <= 1200 )
/// but prefer to do it on NT ... as a matter of fact
/// there is a min max template function but in xutility 
/// and defined as _cpp_max, _cpp_min

#ifdef WIN32
		template <class T>
			inline const T& CC_Min(const T& a, const T& b ){
			return b<a? b:a;
		}
// FIXMEFRED: mig.vc8 (22/05/2007 16:01:48): CC_Min -> CC_Min_c
			template <class T>
			inline T CC_Min_c(T a, T b ){
			return b<a? b:a;
		}
		template <class T>
			inline const T& CC_Max(const T& a, const T& b ){
			return b<a? a:b;
		}
// FIXMEFRED: mig.vc8 (22/05/2007 16:01:48): CC_Max -> CC_Max_c
		template <class T>
			inline T CC_Max_c(T a, T b ){
			return b<a? a:b;
		}

        template <class T>
            inline int CC_Round( T x ){
            return int(floor(x+0.5));
        }
#else
	#include <functional>
		#define CC_Min 	std::min
		#define CC_Max 	std::max

        template <class T>
            inline int CC_Round( T x ){
            return int(floor(x+0.5));
        }
#endif

template <class T>
	inline T CC_SignMult(T a, T b ){
	return b>=0.0? abs(a) : -abs(a); }


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

