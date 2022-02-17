/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: functor.h,v $
 * Revision 1.1  2004/02/02 07:52:17  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file functor.h
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPBASE_FUNCTOR_H
#define _INGPBASE_FUNCTOR_H

#include "port.h"
#include <functional>
#include <string>
CC_USING_NS( std, string )

///////////////////////////////////////////////////
/// Unary/BinaryFunc adds on top of the STL unary function
/// the operator definition in order to use derivation
/// instead of pure template design....
/// useful for run time design
///////////////////////////////////////////////////


CC_BEGIN_NAMESPACE( ARM_GP )

/// TEMPLATE ReturnFunc
template<typename T>
	struct ReturnFunc {
		virtual T operator()() const = 0;
		virtual ~ReturnFunc() {};
	};


/// TEMPLATE UnaryFunc
template<typename T,typename U=T>
	struct UnaryFunc : std::unary_function<T,U>{
		virtual U operator()(T _X) const = 0;
		virtual ~UnaryFunc() {};
	};

/// TEMPLATE BinaryFunc
template<typename T,typename U=T, typename V=U>
	struct BinaryFunc : CC_NS( std, binary_function )<T,U,V>{
		virtual V operator()(T _X, U _Y) const =0;
		virtual ~BinaryFunc() {};
	};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/




