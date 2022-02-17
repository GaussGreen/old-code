/****************************************************************************/
/*      Error message handling.                                             */
/****************************************************************************/
/*      ERROR.c                                                             */
/****************************************************************************/
/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/error.c,v 1.7 2003/07/14 18:31:39 dliu Exp $
*/
#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>

#include <DRDiagnosticMsgCB.h>

/*
#undef DLL_EXPORT
#if defined(WIN32) 
#define DLL_EXPORT __declspec( dllexport ) 
#else
#define DLL_EXPORT
#endif
*/

static DRDiagnosticMsgCB ErrorCallbackFnp = NULL;
static DRDiagnosticMsgCB WarningCallbackFnp = NULL;
static DRDiagnosticMsgCB TraceCallbackFnp = NULL;

/**
 * set the error callback function as porposed by Jerry.
 */
void DR_ErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp)
{
	ErrorCallbackFnp = errorCallbackFnp;
}

/**
 * place holder for other additional functions proposed by Jerry.
 */
void DR_WarningSetCallBack(DRDiagnosticMsgCB f) 
{
	WarningCallbackFnp=f;
}
void DR_TraceSetCallBack(DRDiagnosticMsgCB f) 
{
	TraceCallbackFnp=f;
}

/*****  DR_Error  ***********************************************************/
/*
*       Print an error message.
*/
void DR_Error(char* format, ...)
{
	int handled = 0;
	va_list args;

	va_start(args, format);
	if (ErrorCallbackFnp != NULL)
        {
		handled = (*ErrorCallbackFnp)(format, args);
		(*ErrorCallbackFnp)("\n", NULL);
	}

	if (!handled)
        {
		vfprintf(stdout,format,args);
		fprintf(stdout,"\n");
		fflush(stdout);
	}
	va_end(args);
}
