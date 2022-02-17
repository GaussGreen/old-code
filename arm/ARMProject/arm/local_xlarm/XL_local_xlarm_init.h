#ifndef XL_LOCAL_XLARM_INIT_H
#define XL_LOCAL_XLARM_INIT_H

//#include <windows.h>

#ifdef __cplusplus
extern "C"
{
#endif	/* __ cplusplus */
long XLLOCALARM_PersistentListsInit ();
long XLLOCALARM_PersistentListsClear ();
long XLLOCALARM_PersistentListsDelete ();
long XLOCALLARM_deconnexioneToolkit ();
void XLLOCALARM_InitCalendar ();
void XLLOCALARM_InitFileRead ();
void XLLOCALARM_InitEnvVariables ();
void XLLOCALARM_FreeAllObjects ();
void XLLOCALARM_EndCRMTracing ();
void XLLOCALARM_ReleaseGrandPrix ();
void XLLOCALARM_InitGrandPrix ();
#ifdef __cplusplus
}
#endif	/* __cplusplus */

#endif	/* XL_LOCAL_XLARM_INIT_H */

/* EOF %M% */