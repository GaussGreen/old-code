/*-------------------------------------------------------------------------
	Source File:	errbugg.c

	Purpose:	Redirect error handling to buffer
----------------------------------------------------------------------------*/
#include <ctype.h>                 /* toupper */
#include <math.h>		   /* fmod */
#include "bastypes.h"
#include "cgeneral.h"
#include "cerror.h"
#include "convert.h"
#include "macros.h"

#include "errbuff.h"

#define xDEBUG_ERR

static char  * gpErrorBuffer;
static char    gErrorBuffer[ DRL_MAX_BUFFER_LENGTH ];

/* Initialize error logging into internal buffer
 */
PRIVATE TBoolean  errMsgCallBackFunc
(char	* message,		/* (I) input error message  */
 void	* callBackData);	/* (I/O) call back data */


/* Initialize error logging into internal buffer
 */
GTO_EXPORT(void) drlErrorLogInit (void)
{
	GtoErrMsgOn ();

	gpErrorBuffer = gErrorBuffer;

	GtoErrMsgAddCallback ((TErrCallBackFunc*)errMsgCallBackFunc,
			      FALSE, NULL);

#ifdef DEBUG_ERR
    GtoErrMsg("Reset err pointer to the beginning of the buffer.\n");
#endif

    return;
}


/* Jim Mitsiopoulos August, 1998
 *
 * This is the error message call-back function set by
 * drlErrorLogInit().
 *
 * Note that the callback always returns TRUE so that
 * the error message will also be written to the
 * error log, if it is active.
 */

PRIVATE TBoolean  errMsgCallBackFunc
(char	* message,		/* (I) input error message  */
 void	* callBackData		/* (I/O) call back data */
)
{
    long    msglen;
    long    buflen;
    long    shift=0;

    callBackData = callBackData; /* Avoid unused variable warnings */

    if (message ISNT NULL)
    {
        msglen  = strlen (message);
        buflen  = strlen (gErrorBuffer);

        /* If message is too long, reset ptr to the
         * beginning of the buffer, and/or only print
         * the tail of the message
         */
        if (msglen + buflen > DRL_MAX_BUFFER_LENGTH)
        {
            gpErrorBuffer = gErrorBuffer;
            shift = MAX(msglen-buflen,0);
            message += shift;
            msglen -= shift;
        }

        strcpy (gpErrorBuffer, message);
        gpErrorBuffer += msglen;
    }
    return (TRUE);
}

/*f  Jim Mitsiopoulos August, 1998
 *
 * Retrieves error messages that are written to the
 * global buffer by the error message call back function.
 *
 * Returns FAILURE if outBufLen is less than the length
 * of the accumulated error message string
 */

GTO_EXPORT(int)   drlRetrieveErrorMessages
(long    outBufLen,    /* (I) length of the input buffer. */
 char  * outBuffer     /* (I/O) output buffer. */
)
{
    static int status = FAILURE;

    long   globStrlen = strlen( gErrorBuffer );

    if (outBuffer IS NULL ||
        outBufLen < globStrlen)
    {
        goto done;      /* Failed */
    }

    strcpy (outBuffer, gErrorBuffer);

    gpErrorBuffer = gErrorBuffer;

    status = SUCCESS;

  done:

    return (status);
}
