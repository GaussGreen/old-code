#ifndef ESL_ERROR_DOT_H
#define ESL_ERROR_DOT_H

#include "irx/error.h"

#define DR_Error                irxError
#define DRDiagnosticMsgCB       IrxTErrorCallBackV

#ifdef  __cplusplus
extern "C" {
#endif

void DR_ErrorSetCallBack  (IrxTErrorCallBackV callback);
void DR_WarningSetCallBack(IrxTErrorCallBackV callback);
void DR_TraceSetCallBack  (IrxTErrorCallBackV callback);

/* esl_warning.cpp */
void    DR_Warning (int, char *, ...);

#ifdef  __cplusplus
}
#endif

#endif
