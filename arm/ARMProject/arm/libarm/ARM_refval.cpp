#include "ARM_interglob.h"
#include "ARM_refval.h"



long ARM_CONSTREFVALUE (double value, ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_CONSTREFVALUE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CONSTREFVALUE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, value, idx);
						
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
               
long ARM_IATHREELEVREFVAL (double value,
						   double level0,
						   double amort0,
						   double level1,
						   double amort1,
						   double level2,
						   double amort2,
						   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
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
				ARM_REQUEST_ID = RPC_SET_IA3LEVREFVAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_IA3LEVREFVAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, value, idx);
			REQ_SetDouble (reqIn, level0, ++idx);
			REQ_SetDouble (reqIn, amort0, ++idx);
			REQ_SetDouble (reqIn, level1, ++idx);
			REQ_SetDouble (reqIn, amort1, ++idx);
			REQ_SetDouble (reqIn, level2, ++idx);
			REQ_SetDouble (reqIn, amort2, ++idx);
						
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

long ARM_REFVALUE (VECTOR<double>& dates,
				   VECTOR<double>& values,
				   VECTOR<double>& values2,
				   long valueType,
				   long conversion, long calcMethod,
				   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 8;
	CCString stringObjectId;
			
	try
	{
		/*--- parameters checking ---*/
		if( (dates.size() != values.size ())
			|| ((dates.size() != values2.size ()) && (values2.size () != 0))
		  )
		{
			result.setMsg ("ARM_ERR: dates and values must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_REFVALUE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_REFVALUE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> dates_str;

			for(int i = 0; i < dates.size (); i++)
			{
				dates_str.push_back (XLDATE2ARMDATE(dates[i]));
			}

			REQ_SetLong (reqIn, dates.size (), idx);
			REQ_SetStringVector (reqIn, dates_str, ++idx);
			REQ_SetDoubleVector (reqIn, values, ++idx);
			REQ_SetLong (reqIn, values2.size (), ++idx);
			REQ_SetDoubleVector (reqIn, values2, ++idx);
			REQ_SetLong (reqIn, valueType, ++idx);
			REQ_SetLong (reqIn, conversion, ++idx);
            REQ_SetLong (reqIn, calcMethod, ++idx);

						
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