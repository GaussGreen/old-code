/******************************************************************************
 * Module:	Q3
 * Submodule:
 * File: errbuf.c
 * Function:	
 * Author:	Interest Rates DR
 * Revision:	$Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>                              /* toupper */
#include <string.h>

#include "cerror.h"                             /* GtoErrMsg */
#include "macros.h"                             /* NEW_ARRAY, FREE */
#include "dtlist.h"                             /* TDateList */
#include "busday.h"                             /* GtoHolidayLoadFromDateList */

/*#include "q3quasi.h"                            */
#include "mqbv.h"

#define Q3_MAX_BUFFER_LENGTH   8136L        /* size of internal buffer */


static char    *gpErrorBuffer; 
static char    gErrorBuffer[Q3_MAX_BUFFER_LENGTH];

static TBoolean  Q3BivarErrMsgCallBackFunc( 
     char   *message,            /* (I) input error message  */ 
     void   *callBackData);      /* (I/O) call back data */ 


/*f----------------------------------------------------------------------------
 * Q3BivarErrBuffInit
 *
 * Initialize error logging into internal buffer.
 */
void Q3BivarErrBuffInit(void)
{
    GtoErrMsgOn (); 

    gpErrorBuffer = gErrorBuffer;
    strcpy (gpErrorBuffer, ""); 

    GtoErrMsgAddCallback ((TErrCallBackFunc*) Q3BivarErrMsgCallBackFunc, 
                          FALSE, NULL); 
    return;
}

/*-----------------------------------------------------------------------------
 * Q3BivarErrMsgCallBackFunc 
 * 
 * This is the error message call-back function set by
 * Q3ErrToAlibErrInit(). 
 * 
 * Note that the callback always returns TRUE so that 
 * the error message will also be written to the 
 * error log, if it is active.  
 */ 

static TBoolean  Q3BivarErrMsgCallBackFunc(
    char   *message,        /* (I) input error message  */ 
    void   *callBackData)   /* (I/O) call back data */ 
{ 
    long    msglen; 
    long    buflen;
    long    shift=0;

    callBackData = callBackData; /* Avoid unused variable warnings */

    if (message != NULL) { 
        msglen  = strlen (message); 
        buflen  = strlen (gErrorBuffer);
 
	/* If message is too long, reset ptr to the
	 * beginning of the buffer, and/or only print
	 * the tail of the message
	 */
        if (msglen + buflen > Q3_MAX_BUFFER_LENGTH) {
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


/*f----------------------------------------------------------------------------
 * Q3BivarErrBuffRetrieve
 *
 * Retrieves error messages that are written to the 
 * global buffer by the error message call back function. 
 * 
 * Returns FAILURE if outBufLen is less than the length
 * of the accumulated error message string
 */ 

int Q3BivarErrBuffRetrieve(
    long    outBufLen,    /* (I) length of the input buffer. */ 
    char   *outBuffer)    /* (I/O) output buffer. */ 
{ 
    static int status = FAILURE;
    long   globStrlen = strlen( gErrorBuffer );

    if (outBuffer IS NULL || outBufLen < globStrlen) goto RETURN;      
    
    strcpy (outBuffer, gErrorBuffer); 
    gpErrorBuffer = gErrorBuffer;

    status = SUCCESS;

  RETURN:

    return (status); 
}


/*f----------------------------------------------------------------------------
 * Q3BivarErrBuffLength
 *
 * Retrieves the length of the internal error buffer
 */

int  Q3BivarErrBuffLength ()
{ 
    long   globStrlen = strlen(gErrorBuffer);

    return (globStrlen); 
}
 

/*f----------------------------------------------------------------------------
 * Q3QuasiHolidayLoadFromDateList
 * 
 * Loads ALIB holidays.
 */ 

int Q3BivarHolidayLoadFromDateList (
    char *holName,                  /* holiday name */
    long *datesL,                   /* counted TDate array of holiday dates */
    long satIsAlwaysHoliday,        /* 1 if Saturday is always holiday */
    long sunIsAlwaysHoliday)        /* 1 if Sunday is always holiday */
{ 
    static char   routine[] = "Q3QuasiHolidayLoadFromDateList";
    int status = FAILURE;

    TDateList *dl = NULL;
    if (datesL == NULL) goto RETURN;

    if ((dl = GtoNewDateListFromDates(
        datesL + 1,
    (int) datesL[0])) == NULL) goto RETURN;

    status = GtoHolidayLoadFromDatelist(
        holName,
        dl,
    (satIsAlwaysHoliday == 1),
    (sunIsAlwaysHoliday == 1));

  RETURN:
    
    GtoFreeDateList(dl);
    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }

    return status;
}



























