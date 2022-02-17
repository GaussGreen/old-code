#ifndef CX_ALIBFUNCLOG_H
#define CX_ALIBFUNCLOG_H

#include <alib/cgeneral.h>
#include <alib/addinlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* candidate for becoming a Gto... function - so we use TBoolean */
/* if fileName is NULL, then uses GtoErrMsg for the output       */
int CxAlibFuncLog (char* filename, TBoolean append, char *funcname, ...);

/* call this function if you are not calling from a regular alib add-in */
void CxAlibFuncLogInit (char *libname, TRegisterObjectTypesFunc regFunc);

#ifdef __cplusplus
}
#endif

#endif

