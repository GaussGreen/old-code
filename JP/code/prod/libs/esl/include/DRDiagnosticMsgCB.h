/***************************************************************************
* Header file to define types for message callback functions and declare
* message callback registration functions.
*
* Author: Jerry Cohen
***************************************************************************/ 

#ifndef _DRDiagnosticMsgCB_h
#define _DRDiagnosticMsgCB_h

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*DRDiagnosticMsgCB)(const char*, va_list);

extern void DR_ErrorSetCallBack(DRDiagnosticMsgCB);
extern void DR_WarningSetCallBack(DRDiagnosticMsgCB);
extern void DR_TraceSetCallBack(DRDiagnosticMsgCB);

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
