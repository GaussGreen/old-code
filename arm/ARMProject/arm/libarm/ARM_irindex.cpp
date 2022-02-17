#include "ARM_interglob.h"
#include "ARM_irindex.h"




long ARM_IRINDEX(long dayCountId,
				 long frequencyId,
				 double maturity,
				 long compMethId,
				 long fwdRuleId,
				 long resetTimingId,
				 long resetGap,
				 long payTimingId,
				 long payGap,
				 long ccyId,
				 long indexType,
				 long decompFreq,
				 ARM_result& result,
				 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 12;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_IRINDEX;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_IRINDEX;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, dayCountId, idx);
			REQ_SetLong (reqIn, frequencyId, ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, compMethId, ++idx);
			REQ_SetLong (reqIn, fwdRuleId, ++idx);
			REQ_SetLong (reqIn, resetTimingId, ++idx);
			REQ_SetLong (reqIn, resetGap, ++idx);
			REQ_SetLong (reqIn, payTimingId, ++idx);
			REQ_SetLong (reqIn, payGap, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, indexType, ++idx);
			REQ_SetLong (reqIn, decompFreq, ++idx);

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



long ARM_IRINDEX_MONEY_MARKET (const CCString& mmTerm,
							   const CCString& ccy,
							   ARM_result& result,
							   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
    
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_IRINDEX_MONEY_MARKET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_IRINDEX_MONEY_MARKET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, mmTerm, idx);
			REQ_SetString (reqIn, ccy, ++idx);
			
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


long ARM_LIBOR (long liborTypeId,
				long ccyId,
				long resetFreqId,
				long payFreqId,
				ARM_result& result,
				long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SET_LIBOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, liborTypeId, idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, resetFreqId, ++idx);
			REQ_SetLong (reqIn, payFreqId, ++idx);
						
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




long ARM_CMS(long CMSType,
			 long liborTypeId,
			 long ccyId,
			 ARM_result& result,
			 long objId)
{
    long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if ( objId != -1 )
			{
			   ARM_REQUEST_ID = RPC_SET_CMS;
			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			   idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATE_CMS;
			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			   idx = 0;
			}

			REQ_SetLong (reqIn, CMSType, idx);
            REQ_SetLong (reqIn, liborTypeId, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			
						
			/*--- CORBA Server Call ---*/

 	  		CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

			result.set(reqOut);

			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}
	
    catch (const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return(ARM_KO);
}



long ARM_FixedIndex(long dayCountId,
					const CCString& ccy,
					ARM_result& result,
					long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
			
	try
	{
		ARM_CORBA_init();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if ( objId != -1 )
			{
			   ARM_REQUEST_ID = RPC_SET_FIX;
			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			   idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATE_FIX;
			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			   idx = 0;
			}

			REQ_SetLong (reqIn, dayCountId, idx);
			REQ_SetString (reqIn, ccy, ++idx);
			
						
			/*--- CORBA Server Call ---*/

 	  		CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

			result.set(reqOut);

			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}
	
    catch (const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return(ARM_KO);
}


/*---- End Of File ----*/

// EOF %M%