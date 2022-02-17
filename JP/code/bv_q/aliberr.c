/******************************************************************************
 * Module:	Q3
 * Submodule:
 * File: rate.c	
 * Function:	
 * Author:	Interest Rates DR
 * Revision:	$Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>                      /* toupper */
#include <stdio.h>

#include "aliberr.h"	
#include "q3.h"	       

#include "cgeneral.h"
#include "macros.h"                     /* MAX, MIN */
#include "cerror.h"

/*-----------------------------------------------------------------------------
 * q3ToAlibErrCallBack
 */
static	int q3ToAlibErrCallback (const char *string)
{
    char *s = MALLOC(strlen(string) + 1);
    if (s == NULL) return (FAILURE);
    strcpy(s, string);

    GtoErrMsg(s);

    FREE(s);
    return(SUCCESS);
}


/*f----------------------------------------------------------------------------
 * Q3ToAlibErrInit
 *
 * Sets the callback for alib error messages
 */

void Q3ToAlibErrInit()
{
    static done = FALSE;
    if (!done) {
	Q3ErrCallbackSet(q3ToAlibErrCallback);
	done = TRUE;
    }
} /* Q3ToAlibErrInit */





