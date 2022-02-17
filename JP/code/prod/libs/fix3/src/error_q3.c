/************************************************************************
 * Module:	Q3TMX
 * Function:	Error handling routines
 * Author:	
 * Revision:	$Header$
 ************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "q3_tmx.h"			/* Prototype consistency */

#if defined(WIN32) || defined(MSDOS) || defined(_WIN32)
#define vsnprintf       _vsnprintf
#endif


static	int	(*q3ErrCallback)(const char *string) = NULL;


/*f-------------------------------------------------------------
 * Error logging : write an error message
 *
 * <br><br>
 * Writes a formatted message to a log file.
 * The default is to write to stderr. This default can be
 * overwritten by using Q3TMXErrMsgAddCallback.
 * 
 */

void
Q3TMXErrMsg(char *fmt,...)
{
	va_list         ap;
	char    	*buffer = NULL;
	int		buffTrunc = 0;	/* FALSE */
#define ERRBUFFSIZE	4096		/* Error buffer size for truncation */

	char	error_msg_append[] = "\nERROR MESSAGE TOO LONG, TRUNCATED.\n";
	long	error_msg_size;
	long	return_value;


	/* Allocate memory for buffer */
	error_msg_size = ERRBUFFSIZE + strlen(error_msg_append) + 1;
	buffer = (char*) malloc(error_msg_size * sizeof(char));

	if (buffer == NULL)
	{
		Q3TMXErrMsg("Q3TMXErrMsg: Memory allocation failed.\n");
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
	if (q3ErrCallback == NULL) {
		/* Default: print to stderr */
		fputs(buffer, stdout);
	} else {
		/* Use callback provided */
		(*q3ErrCallback)(buffer);
	}

	if (buffer != NULL) free(buffer);
	return;
}


/*f-------------------------------------------------------------
 * Error logging : place a callback in the error meesage handler.
 *
 * <br><br>
 * cb func is a user provided error callback routine, with prototype
 *	int (*cbfunc)(const char *string),
 * where string is the error message and cbdata is optional
 * user provided data to be passed at each call.
 * This routine can revert to the default mechanism by passing
 * NULL as arguments.
 */

int
Q3TMXErrCallbackSet(
	int (*cbfunc)(const char*))	/* (I) callback routine */
{
	q3ErrCallback = cbfunc;
	return(0);
}



