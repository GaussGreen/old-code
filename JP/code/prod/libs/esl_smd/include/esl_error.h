#ifndef ESL_ERROR_DOT_H
#define ESL_ERROR_DOT_H

/** NOTE: This file should be only included through 'esl_error.c'
 */


#ifdef  __cplusplus
extern "C" {
#endif


#include "DRDiagnosticMsgCB.h"


/**
 * set the error callback function as porposed by Jerry.
 */
void DR_ErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp);

/**
 * place holder for other additional functions proposed by Jerry.
 */
void DR_WarningSetCallBack(DRDiagnosticMsgCB f);

void DR_TraceSetCallBack(DRDiagnosticMsgCB f);

/*****  DR_Error  ***********************************************************/
/**
*       Print an error message.
*/
void DR_Error(char const* format, ...);


#ifdef  __cplusplus
}
#endif


#endif



