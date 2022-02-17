#include "ARM_interglob.h"
#include "ARM_bdfut.h"



long ARM_BDFUT (double delivery,
			    long underIsBd,
				long underId,
				double coupon,
				double convFactor,
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
				ARM_REQUEST_ID = RPC_SETBDFUT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBDFUT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(delivery), idx);
			REQ_SetLong (reqIn, underIsBd, ++idx);
			REQ_SetLong (reqIn, underId, ++idx);
			REQ_SetDouble (reqIn, coupon, ++idx);
			REQ_SetDouble (reqIn, convFactor, ++idx);
						
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

long ARM_GetConversionFactor (long bdFutId,
							  long factId,
							  ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GETCONVERSIONFACTOR;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bdFutId, idx);
			REQ_SetLong (reqIn, factId, ++idx);
			
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

long ARM_GetCheapest (long bdFutId,
					  ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 1;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GETCHEAPEST;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bdFutId, idx);
			
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

long ARM_GILT_NOTIONNAL_BUND (const CCString &delivery,
							  long underId,
							  long notioOrGilt,
							  long market,
							  ARM_result& result,
							  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
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
				ARM_REQUEST_ID = RPC_SET_NOTIONNAL_GILT_BUND;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_NOTIONNAL_GILT_BUND;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)delivery, idx);
			REQ_SetLong (reqIn, underId, ++idx);
			REQ_SetLong (reqIn, notioOrGilt, ++idx);
			REQ_SetLong (reqIn, market, ++idx);
						
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

long ARM_GILT (const CCString &delivery,
			   long underId,
	  		   ARM_result& result,
			   long objId)
{
	return ARM_GILT_NOTIONNAL_BUND (delivery, underId, 2, 1, result, objId);
}

long ARM_NOTIONNAL (const CCString &delivery,
			        long underId,
	  				ARM_result& result,
					long objId)
{
	return ARM_GILT_NOTIONNAL_BUND (delivery, underId, 0, 0, result, objId);
}

long ARM_BUND_LIFFE (const CCString &delivery,
			         long underId,
	  				 ARM_result& result,
					 long objId)
{
	return ARM_GILT_NOTIONNAL_BUND (delivery, underId, 1, 1, result, objId);
}

long ARM_BUND_DTB (const CCString &delivery,
			       long underId,
	  			   ARM_result& result,
				   long objId)
{
	return ARM_GILT_NOTIONNAL_BUND (delivery, underId, 1, 2, result, objId);
}

long ARM_BDFUTBASKET (double delivery,
					  long underIsBd,
					  long underId,
					  double coupon,
					  double convFactor,
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
				ARM_REQUEST_ID = RPC_SETBDBASKFUT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBDBASKFUT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(delivery), idx);
			REQ_SetLong (reqIn, underIsBd, ++idx);
			REQ_SetLong (reqIn, underId, ++idx);
			REQ_SetDouble (reqIn, coupon, ++idx);
			REQ_SetDouble (reqIn, convFactor, ++idx);
						
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