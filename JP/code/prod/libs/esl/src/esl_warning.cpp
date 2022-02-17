/****************************************************************************/
/*      Warning message handling.                                             */
/****************************************************************************/
/*      WARNING.c                                                             */
/****************************************************************************/


/*
$Header$
*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "eslhead.h"
#include "esl_error.h"

#if defined(_MSC_VER)
#define vsnprintf _vsnprintf
#endif

static IrxTErrorCallBackV warningCallBackV = NULL;

void DR_WarningSetCallBack(IrxTErrorCallBackV callback)
{
    warningCallBackV = callback;
}


/*****  DR_Warning  ***********************************************************/
/*
*       Print a warning message to stdout and warning.log
*       Note that DR_Error is deliberately called rather than sprintf so that
*       it is compatible with flows libraries which exclude error.c
*       If openNewFile = FALSE then the warning message is appended to the
*       existing warning.dat file; otherwise a new file will be created 
*/
void    DR_Warning (int openNewFile, char  *format, ...)
{

    char    FileName[MAXBUFF];
    char    WarningMsg[MAXBUFF];
    
    FILE    *stream = NULL;
    int     handled = 0;

    va_list parminfo;
    va_start(parminfo, format);

    if (warningCallBackV != NULL)
    {
       handled = (*warningCallBackV)(format, parminfo);
       (*warningCallBackV)("\n", NULL);
    }

    if (!handled)
    {
	/* Use the format string to print the variable argument list to WarningMsg buffer */
	vsnprintf(WarningMsg, MAXBUFF, format, parminfo);

	/* Open warning.log and print warning to file */
        strcpy (FileName, "warning.log");

        if (openNewFile == FALSE)
        {
            stream = fopen (FileName, "a");
        }
        else
        {
            stream = fopen (FileName, "w");
        }

        if (stream == NULL)
        {
            sprintf  (WarningMsg, "Could not open file %s!",FileName);
            DR_Error (WarningMsg);
        }
        else
        {
            fprintf (stream, "%s \n", WarningMsg);
            fclose (stream);
        }
    }
    va_end(parminfo);
}  /* DR_Warning */
