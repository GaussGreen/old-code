#include "ARM_interglob.h"
#include "ARM_rgnote.h"



long ARM_RNGNOTE (long swapLegId,
				  double lowerBound,
				  double upperBound,
				  long recordFreq,
				  double accruedRate,
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
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
				ARM_REQUEST_ID = RPC_SET_RNGNOTE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_RNGNOTE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetDouble (reqIn, lowerBound, ++idx);
			REQ_SetDouble (reqIn, upperBound, ++idx);
			REQ_SetLong (reqIn, recordFreq, ++idx);
			REQ_SetDouble (reqIn, accruedRate, ++idx);
						
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



long ARM_LIBORRNGNOTE (double startDate,
					   double endDate,
					   double lowerBound,
					   double upperBound,
					   long recordFreq,
					   long liborType,
					   double accruedRate,
					   double spread,
					   long resetFreq,
					   long payFreq,
					   long ccyId,
					   ARM_result& result,
					   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;
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
				ARM_REQUEST_ID = RPC_SET_LIBORRNGNOTE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORRNGNOTE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (endDate), ++idx);
			REQ_SetDouble (reqIn, lowerBound, ++idx);
			REQ_SetDouble (reqIn, upperBound, ++idx);
			REQ_SetLong (reqIn, recordFreq, idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetDouble (reqIn, accruedRate, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
									
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