#include "ARM_result.h"


#ifndef _LOCAL_ARM
	#include "ARM_CORBA_REQUEST.h"
	#include "ARM_CORBA_REQUEST_clnt.h"








ARM_result::ARM_result (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
	set (req);
}



void ARM_result::set (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
	retCode = REQ_GetLong (req, 0);
	value = REQ_GetNthParam (req, 1)->val;
	msg = CCString ("ARM_ERR: ") + REQ_GetCCString (req, 2);
}

#endif

// EOF %M%
