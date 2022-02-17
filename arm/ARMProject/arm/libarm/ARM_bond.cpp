#include "ARM_interglob.h"
#include "ARM_bond.h"



long ARM_GetBondFromSummit (const CCString& rga, ARM_result& result, 
                            long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 1;
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
				ARM_REQUEST_ID = RPC_SET_CREATED_BOND_FROM_SUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GET_BOND_FROM_SUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, rga, idx);
			
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

long ARM_bond (double issueDate, double maturityDate,
		       double firstCouponDate, double couponRate,
			   double redemptionPrice, long periodicity,
			   long dayCount, long settleGap, long couponDateFlag,
			   ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETBOND;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBOND;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(issueDate), idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(maturityDate), ++idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(firstCouponDate), ++idx);
			REQ_SetDouble (reqIn, couponRate, ++idx);
			REQ_SetDouble (reqIn, redemptionPrice, ++idx);
			REQ_SetLong (reqIn, periodicity, ++idx);
			REQ_SetLong (reqIn, dayCount, ++idx);
			REQ_SetLong (reqIn, settleGap, ++idx);
			REQ_SetLong (reqIn, couponDateFlag, ++idx);
						
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

long ARM_YTOPRICE (long bondId, double settlement, double yield, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_YIELDTOPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, yield, ++idx);
						
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

long ARM_PTOYIELD (long bondId, double settlement, double price, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_PRICETOYIELD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, price, ++idx);
						
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

long ARM_SetYield (long bondId, double yield, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_SETYIELD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetDouble (reqIn, yield, ++idx);
						
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

long ARM_CalcDuration (long bondId, double settlement, double yield, long flagCpn, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_YIELDTODURATION;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, yield, ++idx);
			REQ_SetLong (reqIn, flagCpn, ++idx);
						
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

long ARM_YTOCONVEXITY (long bondId, double settlement, double actuRate, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_YIELDTOCONVEXITY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, actuRate, ++idx);
						
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

long ARM_YTODURATION (long bondId, double settlement, double actuRate, long flagCpn, ARM_result& result)
{
	return ARM_CalcDuration (bondId, settlement, actuRate, flagCpn, result);
}

long ARM_BDFAPRICE (long bondId, double settlement,
					double actuPrice, double forwardDate, 
					double repoRate, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BONDFORWARDACTUARIELPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, actuPrice, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(forwardDate), ++idx);
			REQ_SetDouble (reqIn, repoRate, ++idx);
						
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

long ARM_BDREPORATE (long bondId, double settlement, 
					 double actuPrice, double forwardDate,
					 double forwardPrice,
					 ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BONDREPORATEFORFORWARDDATE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, bondId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
			REQ_SetDouble (reqIn, actuPrice, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(forwardDate), ++idx);
			REQ_SetDouble (reqIn, forwardPrice, ++idx);
						
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