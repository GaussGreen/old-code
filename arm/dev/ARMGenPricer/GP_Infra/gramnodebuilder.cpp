/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramnodebuilder.cpp
 *
 *  \brief file to build a gramnode
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramnodebuilder.h"
#include "gpinfra/gramnode.h"


CC_BEGIN_NAMESPACE( ARM )

ARM_ExpNodeFunc* GramNodeFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate   )
{
	return	new ARM_ExpNodeFunc( Function, Args, evalDate );
}

ARM_ExpNodeFunc* GramNodeExerciseFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate )
{
	return	new ARM_ExpNodeExerciseFunc( Function, Args, evalDate );
}

ARM_ExpNodeFunc* GramNodeTriggerFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate )
{
	return	new ARM_ExpNodeTriggerFunc( Function, Args, evalDate );
}

ARM_ExpNodeFunc* GramNodeLinearPVFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate )
{
	return new ARM_ExpNodeLinearPVFunc( Function, Args, evalDate );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


