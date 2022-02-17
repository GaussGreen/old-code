/****************************************************************************/
/*      Warning message handling.                                             */
/****************************************************************************/
/*      WARNING.c                                                             */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/warning.c,v 1.3 2002/06/25 12:36:23 lcooper Exp $
*/


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "fix123head.h"





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

	va_list parminfo;

	/* Use the format string to print the variable argument list to WarningMsg buffer */
	va_start(parminfo, format);
	vsprintf(WarningMsg, format, parminfo);
	va_end(parminfo);

	/* Print to stdout stream */
	/* DR_Error(WarningMsg);  */

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

    return;

}  /* DR_Warning */
