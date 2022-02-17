/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions headers setting ARM_CORBA_REQUEST parameters (for client) */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#define ARM_CORBA_REQUEST_clnt_cpp
#include "ARM_CORBA_REQUEST_clnt.h" 


void REQ_SetLongVector(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, const VECTOR<long>& vecVals, long rank)
{
	PARAM_SetLongVector(&req->paramsList[rank], vecVals);
}


void REQ_GetLongVector(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<long>& vecVals, long rank)
{
	PARAM_GetLongVector(&req->paramsList[rank], vecVals); 
}


void REQ_SetDoubleVector(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, const VECTOR<double>& vecVals, long rank)
{
	PARAM_SetDoubleVector(&req->paramsList[rank], vecVals);
}


void REQ_GetDoubleVector(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<double>& vecVals, long rank)
{
	PARAM_GetDoubleVector (&req->paramsList[rank], vecVals);
}


void REQ_SetStringVector(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, const VECTOR<CCString>& vecVals, long rank)
{
	PARAM_SetStringVector(&req->paramsList[rank], vecVals);
}


void REQ_GetStringVector(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, VECTOR<CCString>& vecVals, long rank)
{
	PARAM_GetStringVector(&req->paramsList[rank], vecVals);
}


void REQ_Print_message(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
	long i;

	MSG_printf_message (MSG_INFO, "===============> REQ");
 
	MSG_printf_message (MSG_INFO, "ReqID = %ld", req->reqId);
	MSG_printf_message (MSG_INFO, "NbPar = %ld", req->nbParams);
 
	for (i = 0; i < req->nbParams; i++)
	{
		PARAM_Print_message (&req->paramsList[i], i);
	}
 
	MSG_printf_message (MSG_INFO, "<=============== REQ");
}


long REQ_GetObjectId (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
	return REQ_GetLong (req, 0);
}


CCString REQ_GetCCString (const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long rank)
{
	CCString result;
	char* buf = PARAM_GetString(&req->paramsList[rank]);

	if(buf)
	{
		result.Set (buf);
	}

	return (result);
}

void REQ_Delete (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
	if(req)
	{
		delete req;
		req = NULL;
	}
}

// EOF %M%
