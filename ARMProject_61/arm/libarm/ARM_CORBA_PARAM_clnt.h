/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions for ARM_CorbaRequestModule::ARM_CORBA_PARAM (for client)  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 

#ifndef ARM_CORBA_PARAM_CLNT_H 
#define ARM_CORBA_PARAM_CLNT_H


#include <CCString.h>

#include "req.h"
 

extern void PARAM_SetLongVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<long>& vecVals);
extern void PARAM_GetLongVector(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<long>& vecVals);

extern void PARAM_SetDoubleVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<double>& vecVals);
extern void PARAM_GetDoubleVector(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<double>& vecVals);

extern void PARAM_SetStringVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<CCString>& vecVals);
extern void PARAM_GetStringVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<CCString>& vecVals);

extern void PARAM_Print_message(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long idx);


#ifdef ARM_CORBA_PARAM_clnt_cpp

#include <CCmessage.h>

#include "ARM_CORBA_PARAM.h"

#endif 	// ARM_CORBA_PARAM_clnt_cpp

#endif 	// ARM_CORBA_PARAM_H


// EOF %M%
