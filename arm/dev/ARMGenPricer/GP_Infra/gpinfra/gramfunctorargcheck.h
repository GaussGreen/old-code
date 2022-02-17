/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorargcheck.h
 *
 *  \brief various very simple class for checking arguments size and type
 *		resolves to nothing in release mode!
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_GRAMFUNCTORARGCHECK_H
#define _INGPINFRA_GRAMFUNCTORARGCHECK_H

#include "gpbase/port.h"
#include "gramfunctorarg.h"

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_ExpNode;



struct GPAF_CheckArgSize{
	GPAF_CheckArgSize( const ARM_GramFctorArgVector& arg, int nb, const string& funcName );
};

struct GPAF_CheckArgVecSameSize{
	GPAF_CheckArgVecSameSize( const ARM_GramFctorArgVector& arg, const string& funcName );
};


struct GPAF_CheckArgType{
	GPAF_CheckArgType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type, const string& funcName );
};


struct GPAF_CheckArgReturnType{
	GPAF_CheckArgReturnType( const ARM_GramFctorArg& arg, ARM_GramFctorArg::Type, const ARM_ExpNodePtr& node, const string& funcName );
};


struct GPAF_Pb{
	static void UnexpectedType( const ARM_GramFctorArg::Type FoundType, 
		ARM_GramFctorArg::Type RequiredType, const string& funcName );
};

	
CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

