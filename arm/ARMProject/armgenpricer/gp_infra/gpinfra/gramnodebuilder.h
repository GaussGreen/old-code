/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file gramnodebuilder.h
 *
 *  \brief file for creating the corresponding gram node
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPINFRA_GRAMNODEBUILDER_H
#define _INGPINFRA_GRAMNODEBUILDER_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"


#include "gpbase/port.h"
#include "typedef.h"

/// forward declaration
class ARM_ExpNodeFunc;
class ARM_GramFunction;


CC_BEGIN_NAMESPACE( ARM )

ARM_ExpNodeFunc* GramNodeFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate   );
ARM_ExpNodeFunc* GramNodeExerciseFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate );
ARM_ExpNodeFunc* GramNodeLinearPVFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate );
ARM_ExpNodeFunc* GramNodeTriggerFuncBuilder( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate );

CC_END_NAMESPACE()

#endif 
