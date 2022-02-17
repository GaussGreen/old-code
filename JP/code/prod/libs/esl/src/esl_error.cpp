#include  "esl_error.h"

/*---------------------------------------------------------------------------
 Need these entry points to satisfy Kapital client
---------------------------------------------------------------------------*/
void DR_ErrorSetCallBack(IrxTErrorCallBackV callback)
{
    irxSetErrorCallBackV(callback);
}

void DR_TraceSetCallBack(IrxTTraceCallbackV callback)
{
    irxSetTraceCallBackV(callback);
}


