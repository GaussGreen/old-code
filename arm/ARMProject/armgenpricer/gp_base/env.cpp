/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: env.cpp,v $
 * Revision 1.1  2003/10/16 07:52:11  ebenhamou
 * Initial revision
 *
 */


/*! \file env.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/env.h"
#include <stdlib.h>


CC_BEGIN_NAMESPACE( ARM )

string ARM_USERNAME		= "";
string ARM_COMPUTERNAME	= "";

void LOCALARM_InitEnvVariables() 
{
	ARM_USERNAME		= getenv( "USERNAME" );
	ARM_COMPUTERNAME	= getenv( "COMPUTERNAME" );
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

