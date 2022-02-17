/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file nagfunction.cpp
 *  \brief nagfunction provides some simple function for
 *		nag optimisation
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */

#include "gpcalib/nagfunction.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

void NAG_CALL NagConfunThrowExceptionIfUsed(Integer n, Integer m, Integer needc[], double x[],
	double conf[], double cjac[], Nag_Comm* comm)
{
	throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +": you should not come here!" );
}

////////////////////////////////////////////////////////////////////
/// \function NagErrorhandler
/// \brief
////////////////////////////////////////////////////////////////////
void NAG_CALL NagErrorhandler( char *strng, int code,  char *name)
{
	if ((code != NE_NOERROR)&&(code != NW_KT_CONDITIONS)&&(code != NW_NOT_CONVERGED))
	{
		string ss1("Error or warning from ");
		string ss2(name);
		string ss3(strng);
		string ss=ss1+ss2+ss3;
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +ss.c_str());
	}
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

