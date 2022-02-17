/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions headers setting ARM_CORBA_REQUEST parameters (for client) */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef _ARM_CORBA_REQUEST_CLNT_H
#define _ARM_CORBA_REQUEST_CLNT_H


#include <CCString.h>

#include "ARM_CORBA_PARAM_clnt.h" 


extern void REQ_SetLongVector (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, const VECTOR<long>& vecVals, long rank);
extern void REQ_GetLongVector (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<long>& vecVals, long rank);

extern void REQ_SetDoubleVector (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, const VECTOR<double>& vecVals, long rank);
extern void REQ_GetDoubleVector (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<double>& vecVals, long rank);

extern void REQ_SetStringVector (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                                 const VECTOR<CCString>& vecVals, long rank);
extern void REQ_GetStringVector (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<CCString>& vecVals, long rank);

extern void REQ_Print_message (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);

extern long REQ_GetObjectId (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);

extern CCString REQ_GetCCString (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long rank);

extern void REQ_Delete (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);


#ifdef ARM_CORBA_REQUEST_clnt_cpp

#include <CCmessage.h>

#include "ARM_CORBA_REQUEST.h"
#include "ARM_CORBA_PARAM.h"

#endif	// ARM_CORBA_REQUEST_clnt_cpp


#endif 	// _ARM_CORBA_REQUEST_CLNT_H


// EOF %M%
