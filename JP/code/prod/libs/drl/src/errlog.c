/************************************************************************
 * Module:	DRL
 * Function:	Error handling routines
 * Author:	
 * Revision:	$Header$
 ************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "drlerr.h"		/* Prototype consistency */

#if defined(WIN32) || defined(MSDOS) || defined(_WIN32)
#define vsnprintf       _vsnprintf
#endif


static	int	(*drlErrCallback)(const char *string) = NULL;
static	int	drlErrOn = 1;		/* On by default */
static	int	drlTraceLevel = 0;	/* Off by default */
static	int	(*drlTraceCallbackSet)(long) = NULL;
static	long	(*drlTraceCallbackGet)() = NULL;

/*f-------------------------------------------------------------
 * Error logging : write an error message
 *
 * <br><br>
 * Writes a formatted message to a log file.
 * The default is to write to stderr. This default can be
 * overwritten bu using DrlErrMsgAddCallback.
 * 
 */

void
DrlErrMsg(char *fmt,...)
{
	va_list         ap;
	char    	*buffer = NULL;
	int		buffTrunc = 0;	/* FALSE */
#define ERRBUFFSIZE	4096		/* Error buffer size for trunction */

	char	error_msg_append[] = "\nERROR MESSAGE TOO LONG, TRUNCATED.\n";
	long	error_msg_size;
	long	return_value;


	if (!drlErrOn) return;

	/* Allocate memory for buffer */
	error_msg_size = ERRBUFFSIZE + strlen(error_msg_append) + 1;
	buffer = (char*) malloc(error_msg_size * sizeof(char));

	if (buffer == NULL)
	{
		DrlErrMsg("DrlErrMsg: Memory allocation failed.\n");
		return;
	}

	va_start(ap, fmt);
	return_value = vsnprintf(buffer, ERRBUFFSIZE, fmt, ap);

	if (return_value > ERRBUFFSIZE ||	/* Solaris */
	    return_value < 0)			/* NT      */
		buffTrunc = 1;
	va_end(ap);
	if (buffTrunc)
		strcat(buffer, error_msg_append);


	/* Print the error message */
	if (drlErrCallback == NULL) {
		/* Default: print to stderr */
		fputs(buffer, stdout);
		fflush(stdout);
	} else {
		/* Use callback provided */
		(*drlErrCallback)(buffer);
	}

	if (buffer != NULL) free(buffer);
	return;
}


/*f-------------------------------------------------------------
 * Error logging : place a callback in the error meesage handler.
 *
 * <br><br>
 * cbfunc is a user provided error callback routine, with prototype
 *	int (*cbfunc)(const char *string),
 * where string is the error message and cbdata is optional
 * user provided data to be passed at each call.
 * This retoune can revert to the default mechanism by passing
 * NULL as arguments.
 */

int
DrlErrCallbackSet(
	int (*cbfunc)(const char*))	/* (I) callback routine */
{
	drlErrCallback = cbfunc;
	return(0);
}



/*f-------------------------------------------------------------
 * Error logging : sets error logging on/off.
 *
 * <br><br>
 */

int
DrlErrOn(int boolean)
{
	drlErrOn = boolean;
	return (0);
}


/*f-------------------------------------------------------------
 * Trace logging : sets tracing level.
 *
 * <br><br>
 */

void
DrlTraceLevelSet(long level)
{

	if (drlTraceCallbackSet) {
		(*drlTraceCallbackSet)(level);
	} else {
		drlTraceLevel = level;
	}
	return;
}


/*f-------------------------------------------------------------
 * Trace logging : get the tracing level.
 *
 * <br><br>
 */

long
DrlTraceLevelGet(void)
{
	if (drlTraceCallbackGet) {
		return (*drlTraceCallbackGet)();
	} else {
		return drlTraceLevel;
	}
}


/*f-------------------------------------------------------------
 * Trace logging : sets the tracing level.
 *
 * <br><br>
 * cbfunc is a user provided error callback routine, with prototype
 *	int (*cbfunc)(int),
 * whcich returns the trace (logging) level 
 * This retoune can revert to the default mechanism by passing
 * NULL as arguments.
 */

int	DrlTraceLevelCallbackSet(int (*cbset)(long), long (*cbget)())
{
	drlTraceCallbackSet = cbset;
	drlTraceCallbackGet = cbget;

	return (0);
}

