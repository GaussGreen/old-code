/*********************************************************************************
 * CRXERROR.H 
 * error utils
 *
 ********************************************************************************/

#ifndef __CRXERROR_H__
#define __CRXERROR_H__

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdarg.h>

#include <common/include/drmacros.h>
    
typedef int (*DRDiagnosticMsgCB)(const char*, va_list);

extern void DR_ErrorSetCallBack(DRDiagnosticMsgCB);
extern void DR_WarningSetCallBack(DRDiagnosticMsgCB);
extern void DR_TraceSetCallBack(DRDiagnosticMsgCB);

void   DR_Error(const char* format, ...);
int    DR_ErrorCallback(const char* format, va_list ap);
const  char*  DR_ErrorRetrieve();

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
    
