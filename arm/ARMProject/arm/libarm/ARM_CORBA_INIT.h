#ifndef ARM_CORBA_INIT_H
#define ARM_CORBA_INIT_H

#ifdef __cplusplus
extern "C" 
{
#endif	// __cplusplus

long ARM_CORBA_init ();
long ARM_CORBA_exit ();
void ARM_SetVirtualDisconnect ();
void ARM_UnsetVirtualDisconnect ();
int ARM_IsVirtualDisconnect ();
int ARM_IsCORBA_initialized ();

#ifdef __cplusplus
}
#endif	// __cplusplus

#endif	// ARM_CORBA_INIT_H

// EOF %M%
