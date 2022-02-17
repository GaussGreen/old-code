#include "ARM_interglob.h"
#include "ARM_barrier.h"



long ARM_CONSTBARRIER (long underlyingId,
					   long tAssetId,
					   double maturity,
					   double barrier,
					   long upDown,
					   long inOut,
					   long triggerVar,
					   double rebate,
					   double firstX,
					   ARM_result& result,
					   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_CONSTBARRIER;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CONSTBARRIER;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, underlyingId, idx);
			REQ_SetLong (reqIn, tAssetId, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(maturity), ++idx);
			REQ_SetDouble (reqIn, barrier, ++idx);
			REQ_SetLong (reqIn, upDown, ++idx);
			REQ_SetLong (reqIn, inOut, ++idx);
			REQ_SetLong (reqIn, triggerVar, ++idx);
			REQ_SetDouble (reqIn, rebate, ++idx);
			if(firstX == -1.0)
			{
				REQ_SetString (reqIn, ARM_DEFAULT_DATE, ++idx);
			}
			else
			{	
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(firstX), ++idx);
			}
						
	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  		result.set (reqOut);
        
			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}
	catch(const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}



long ARM_BARRIER (long underlyingId,
				  long tAssetId,
				  long xStyleId,
				  long refValId,
				  long upDown,
				  long inOut,
				  long triggerVar,
				  double rebate,
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 8;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_BARRIER;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_BARRIER;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, underlyingId, idx);
			REQ_SetLong (reqIn, tAssetId, ++idx);
			REQ_SetLong (reqIn, xStyleId, ++idx);
			REQ_SetLong (reqIn, refValId, ++idx);
			REQ_SetLong (reqIn, upDown, ++idx);
			REQ_SetLong (reqIn, inOut, ++idx);
			REQ_SetLong (reqIn, triggerVar, ++idx);
			REQ_SetDouble (reqIn, rebate, ++idx);
									
	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  		result.set (reqOut);
        
			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}
	catch(const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}



// EOF %M%