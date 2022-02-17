/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Revision 1.1  2003/10/08 16:44:12  ebenhamou, jmprie
 * Initial revision
 *
 */


/*! \file modelparamtype.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/modelparamtype.h"
#include "gpinfra/argconvdefault.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamType
///	Routine: GetTypeString
///	Returns: the string corresponding to the nb of the model param
///	Action : 
////////////////////////////////////////////////////

string ARM_ModelParamType::GetTypeString( ParamNb nb )
{
	return ARM_ArgConvReverse_ModelParam.GetString(nb);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

