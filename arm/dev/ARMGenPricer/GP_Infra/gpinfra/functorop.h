/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: functorop.h,v $
 * Revision 1.1  2003/10/20 07:52:17  ebenhamou
 * Initial revision
 *
 */

 /*! \file functorop.h
 *
 *  \brief simple template function for operation in place
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_FUNCTOROP_H
#define _INGPINFRA_FUNCTOROP_H

#include "gpbase/port.h"

#include <algorithm>

CC_BEGIN_NAMESPACE( ARM )
/// template for all unary function in place
template<typename UnaryOp>
	void FuncUnaryInPlace( ARM_VectorPtr& v, const UnaryOp& op ){
		CC_NS(std,transform)( (*v).begin(), (*v).end(), (*v).begin(), op );
}

/// template for all binary function in place
template <typename BinaryOp>
	void FuncBinaryInPlace( ARM_VectorPtr& v1, ARM_VectorPtr& v2,const BinaryOp& op  ){
		CC_NS(std,transform)( (*v1).begin(),(*v1).end(),(*v2).begin(),(*v1).begin(), op );
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

