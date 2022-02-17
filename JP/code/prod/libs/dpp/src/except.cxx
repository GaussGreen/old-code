/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, D. Liu
 * Revision:	$Header$
 ************************************************************************/
#if defined(_WIN32)
# include <cstdio>
#endif
#include <errno.h>
#include "kstdinc.h"
#include "kerrlog.h"

extern	"C" {
#include "cerror.h"
};


//---------------------------------------------------------------
// Default Constructor 
KFailure::KFailure(int level)
{
	errLevel = level;
}


//---------------------------------------------------------------
// Constructor for {\tt KFailure}.

KFailure::KFailure(const char *fmt, ...)
{
	va_list	ap;
static	char	tmpBuf[512];


	va_start(ap, fmt);
	vsprintf(tmpBuf, fmt, ap);
	va_end(ap);

	// print to stderr 
	if (DppErrMsgStatus() == DPP_ERR_MSG_GTO)
		GtoErrMsg(tmpBuf);
	else if (DppErrMsgStatus() != DPP_ERR_MSG_OFF)
		dppErr << tmpBuf;

	errLevel = -1;
}



