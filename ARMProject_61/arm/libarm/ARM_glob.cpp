#include "ARM_interglob.h"
#include "ARM_glob.h"
#include "ARM_option.h"


char* AlreadyDefaultCountry = NULL;

double ARMDATE2XLDATE (const CCString& armdate)
{
	int y, m, d;

	sscanf ((const char*)armdate, "%2d.%2d.%4d", &d, &m, &y);
	
	return (DAT_struct_to_ssdate (y, m, d));
}



CCString XLDATE2ARMDATE (double xldate)
{
	int y, m, d;
	char buf[11];
	CCString tmp;

	long long_xldate = (long)xldate;

	DAT_ssdate_to_struct ((double)long_xldate, &y, &m, &d);

	if(d < 10)
	{
		sprintf (buf, "0%1d.", d);
	}
	else
	{
		sprintf (buf, "%2d.", d);
	}
	tmp.Set (buf);
	
	if(m < 10)
	{
		sprintf (buf, "0%1d.", m);
	}
	else
	{
		sprintf (buf, "%2d.", m);
	}
	tmp += CCString (buf);

	sprintf (buf, "%4d", y);

	tmp += CCString (buf);

	return (tmp);
}


long ARM_ARM_Price (long secId,
					long modId,
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

			ARM_REQUEST_ID = RPC_PRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
						
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



long ARM_SetPrice (long instId,
				   double price,
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

			ARM_REQUEST_ID = RPC_SETPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, instId, idx);
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



long ARM_GetFRMShortRateVols(long modId,
							 const CCString& outFile, 
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

			ARM_REQUEST_ID = RPC_FRMSHORTRATEVOLS;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong(reqIn, modId, idx);
			REQ_SetString(reqIn, outFile, ++idx);
						
	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  		result.set(reqOut);
        
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



long ARM_ExProba (long secId,
				  long modId,
				  long numEx,
				  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_EXPROBA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetLong (reqIn, numEx, ++idx);
						
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

long ARM_FreeObject (long secId,
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

			ARM_REQUEST_ID = RPC_FREEOBJECT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
						
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

long ARM_FreeAllObjects (ARM_result& result)
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

			ARM_REQUEST_ID = RPC_FREE_ALL;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, 0, idx);
						
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

long ARM_ExitArm ()
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

			ARM_REQUEST_ID = RPC_EXIT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, 0, idx);
						
	  		/*--- CORBA Server Call ---*/
			CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  		/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);
	 	}
	}
	catch(const CORBA::Exception &e)
	{
		MSG_printf_message (MSG_ERROR, CORBA_ERR_MSG);
		if(ARM_logfile)
		{
			ARM_logfile->flush ();
			*ARM_logfile << CORBA_ERR_MSG << endl;
			*ARM_logfile << &e << endl;
			ARM_logfile->flush ();
		}
		else
		{
			cout.flush ();
			cout << CORBA_ERR_MSG << endl;
			cout << &e << endl;
			cout.flush ();
		}
	}

	return ARM_OK;
}

long ARM_ARM_GetDefaultCurrency (ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 1;
	long idx = 0;
			
	try
	{
//		if (AlreadyDefaultCountry == NULL)
//		{
			ARM_CORBA_init ();

			if(CORBA_OBJECT_INTERFACE)
			{
				ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
				ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

				ARM_REQUEST_ID = RPC_GETDEFAULTCOUNTRY;

				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				
				REQ_SetLong (reqIn, 0, idx);
							
	  			/*--- CORBA Server Call ---*/

				CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  			result.set (reqOut);
        
				/*--- requests freeing ---*/
				REQ_Delete (reqOut);
				REQ_Delete (reqIn);

//				AlreadyDefaultCountry = new char[4];
//				strcpy(AlreadyDefaultCountry,result.getString());

				ARM_RESULT();
			}
/*		}
		else
		{
			result.setString(AlreadyDefaultCountry);

			result.setRetCode(ARM_OK);

			ARM_RESULT();
		}
*/
	}
	catch(const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}

long ARM_ARM_SetDefaultCurrency (const CCString& isoCCy,
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

			ARM_REQUEST_ID = RPC_SETDEFAULTCOUNTRY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)isoCCy, idx);
						
	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

	  		result.set (reqOut);

//			strcpy(AlreadyDefaultCountry,result.getString());
			
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

long ARM_NextBusinessDay (double date,
						  const CCString& cal,
						  long days,
						  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_NEXTBUSINESSDAY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetString (reqIn, (const char*)cal, ++idx);
			REQ_SetLong (reqIn, days, ++idx);
			
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

long ARM_IsBusinessDay (double date,
						const CCString& isoccyname,
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

			ARM_REQUEST_ID = RPC_ISBUSINESSDAY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetString (reqIn, (const char*)isoccyname, ++idx);
						
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

long ARM_ADJUSTTOBUSDATE (double date,
						  const CCString& currency,
						  long ruleId,
						  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_ADJUSTBUSDATE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetString (reqIn, (const char*)currency, ++idx);
			REQ_SetLong (reqIn, ruleId, ++idx);
						
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

long ARM_ARM_Accrued (long secId,
					  double fwdDate,
					  long modId,
					  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_ACCRUED;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(fwdDate), ++idx);
			REQ_SetLong (reqIn, modId, ++idx);
						
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

long ARM_Sensitivity (long secId,
					  long modId,
					  long paramId,
					  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_SENSITIVITY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetLong (reqIn, paramId, ++idx);
						
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

long ARM_CvSensitivity (long secId,
						long modId,
						long paramId,
						long viewFlagId,
						const CCString& id,
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

			ARM_REQUEST_ID = RPC_CVSENSITIVITY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetLong (reqIn, paramId, ++idx);
			REQ_SetLong (reqIn, viewFlagId, ++idx);
			REQ_SetString (reqIn, id, ++idx);
						
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

long ARM_SetMarketPrice (long id,
						 double price,
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

			ARM_REQUEST_ID = RPC_SETMARKETPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, id, idx);
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

long ARM_IMPLIEDVOL (long instId,
					 long modelId,
					 double price,
					 ARM_result& result)
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

			ARM_REQUEST_ID = RPC_IMPLIED_VOL;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, instId, idx);
			REQ_SetLong (reqIn, modelId, ++idx);
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

long ARM_ParallelShift (long secId,
						double value,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
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
				ARM_REQUEST_ID = RPC_SPARALLELSHIFT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CPARALLELSHIFT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, secId, idx);
			REQ_SetDouble (reqIn, value, ++idx);
			
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

long ARM_SetNotional (long secId,
					  long rId,
					  double percentRemainder,
					  ARM_result& result)
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

			ARM_REQUEST_ID = RPC_SETAMOUNT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, rId, ++idx);
			REQ_SetDouble (reqIn, percentRemainder, ++idx);

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

long ARM_NextBusinessDayDt (double date,
							const CCString& cal,
							long days,
							ARM_result& result)
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

			ARM_REQUEST_ID = RPC_NEXTBUSINESSDAY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (date), idx);
			REQ_SetString (reqIn, (const char*)cal, ++idx);
			REQ_SetLong (reqIn, days, ++idx);
			
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

long ARM_GetExpiry (long secId,
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

			ARM_REQUEST_ID = RPC_GETEXPIRY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, secId, idx);
			
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

long ARM_FwdPrice (long secId,
				   long modId,
				   double fwdDate,
				   ARM_result& result)
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

			ARM_REQUEST_ID = RPC_FORWARD_PRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (fwdDate), ++idx);
						
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

long ARM_BSSpot (long secId,
				 long modId,
				 double date,
				 ARM_result& result)
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

			ARM_REQUEST_ID = RPC_BSSPOT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (date), ++idx);
						
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

long ARM_SummitSwapInfo (const CCString& tradeId,
						 const CCString& curveId,
						 double valoDate,
						 const CCString& whatTo,
						 ARM_result& result)
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

			ARM_REQUEST_ID = RPC_SUMMIT_SWAP_INFO;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, (const char*)tradeId, idx);
			REQ_SetString (reqIn, (const char*)curveId, ++idx);
			if(valoDate == -1.0)
			{
				REQ_SetString (reqIn, "XX", ++idx);
			}
			else
			{
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (valoDate), ++idx);
			}
			REQ_SetString (reqIn, (const char*)whatTo, ++idx);
									
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

long ARM_FxConvert (const CCString& ccy1,
					const CCString& ccy2,
					double asOfDate,
					double amount,
					const CCString& cvname,
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

			ARM_REQUEST_ID = RPC_FXCONVERT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, ccy1, idx);
			REQ_SetString (reqIn, ccy2, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (asOfDate), ++idx); 
			REQ_SetDouble (reqIn, amount, ++idx);
			REQ_SetString (reqIn, cvname, ++idx);
						
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

long ARM_EvalSummitAsset (const CCString& tradeId,
						  const CCString& tradeType,
						  const CCString& curveId,
						  double asOfDate,
						  long assetId,
						  long isNet,
						  ARM_result &result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_EVALSUMMITASSET;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, tradeId, idx);
			REQ_SetString (reqIn, tradeType, ++idx);
			REQ_SetString (reqIn, curveId, ++idx);
			if(asOfDate == -1.0)
			{
				REQ_SetString (reqIn, "XX", ++idx);
			}
			else
			{
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (asOfDate), ++idx); 
			}
			REQ_SetLong (reqIn, assetId, ++idx);
			REQ_SetLong (reqIn, isNet, ++idx);
						
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

long ARM_KImp (long secId,
			   long modId,
			   double price,
			   long param,
			   ARM_result &result)
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

			ARM_REQUEST_ID = RPC_KIMP;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
			REQ_SetDouble (reqIn, price, ++idx);
			REQ_SetLong (reqIn, param, ++idx);
						
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

long ARM_SummitValueGreeks (const CCString& tradeId,
						    const CCString& cvId,
							const CCString& tradeType,
							double valoDate,
							const CCString& greek,
							const CCString& analytic,
							ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_SUMMIT_VALUE_GREEKS;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, tradeId, idx);
			REQ_SetString (reqIn, cvId, ++idx);
			REQ_SetString (reqIn, tradeType, ++idx);
			if(valoDate == -1.0)
			{
				REQ_SetString (reqIn, "XX", ++idx);
			}
			else
			{
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (valoDate), ++idx);
			}
			REQ_SetString (reqIn, greek, ++idx);
			REQ_SetString (reqIn, analytic, ++idx);
						
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

long ARM_SummitValueTrade (const CCString& tradeId,
						   const CCString& tradeType,
						   const CCString& curveId,
						   double date,
						   long isNetVal,
						   long version,
						   ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_SUMMIT_VALUE_TRADE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, (const char*)tradeId, idx);
            REQ_SetString (reqIn, (const char*)curveId, ++idx);
			REQ_SetString (reqIn, (const char*)tradeType, ++idx);
			
			if(date == -1.0)
			{
				REQ_SetString (reqIn, "XX", ++idx);
			}
			else
			{
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (date), ++idx);
			}
			REQ_SetLong (reqIn, isNetVal, ++idx);
			REQ_SetLong (reqIn, version, ++idx);

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

long ARM_ARM_View (long instId,
				   const CCString& sockId,
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

			ARM_REQUEST_ID = RPC_VIEW;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, instId, idx);
			REQ_SetString (reqIn, sockId, ++idx);

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


long ARM_XCccyAdjustment( long startDate,
                          long endDate, 
				          long payFreq,
                          const CCString& domCcy,
				          long forIndexTypeId,
				          const CCString& forCcy,
                          long spreadsId,
                          long zcDomId,
                          long discDomId,
                          long zcForId,
                          long discForId,
                          double FX,
						  long couponId,
                          long domDc,
						  long forDc,
                          ARM_result& result,
						  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 15;
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
				ARM_REQUEST_ID = RPC_SETXCCYADJUST;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_XCCYADJUST;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), idx);

            REQ_SetString (reqIn, XLDATE2ARMDATE(endDate), ++idx);

			REQ_SetLong   (reqIn, payFreq, ++idx);

			REQ_SetString (reqIn, domCcy, ++idx);

            REQ_SetLong   (reqIn, forIndexTypeId, ++idx);

            REQ_SetString (reqIn, forCcy, ++idx);

            REQ_SetLong   (reqIn, spreadsId, ++idx);

            REQ_SetLong   (reqIn, zcDomId, ++idx);

            REQ_SetLong   (reqIn, discDomId, ++idx);

            REQ_SetLong   (reqIn, zcForId, ++idx);

            REQ_SetLong   (reqIn, discForId, ++idx);

            REQ_SetDouble (reqIn, FX, ++idx);

            REQ_SetLong (reqIn, couponId, ++idx);

            REQ_SetLong (reqIn, domDc, ++idx);

            REQ_SetLong (reqIn, forDc, ++idx);


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


long ARM_ARM_GetPID (ARM_result& result)
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

			ARM_REQUEST_ID = RPC_DISPATCH_PID;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, 0, idx);
						
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


long ARM_ARM_SCredit (const CCString& SummitFilter,
					  const CCString& SummitCurveId,
					  const CCString& FileName,
					  ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_SCREDIT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			long oldAlgo = 0;

			REQ_SetString(reqIn, SummitFilter, idx);
			REQ_SetString(reqIn, SummitCurveId, ++idx);
			REQ_SetString(reqIn, FileName, ++idx);
			REQ_SetLong(reqIn, oldAlgo, ++idx);
						
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



long ARM_ARM_BetweenDates (long date1,
						   long date2,
						   long daycountId,
						   long isYearFrac,
						   ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BETWEENDATES;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong(reqIn, daycountId, idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(date1), ++idx);
            REQ_SetString (reqIn, XLDATE2ARMDATE(date2), ++idx);
			REQ_SetLong(reqIn, isYearFrac, ++idx);

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


long ARM_ARM_ADDYEARS  (long date,
						long nb,
						long ruleId,
						const CCString& Ccy,
						ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_ADDYEARS;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, XLDATE2ARMDATE(date), idx);
			REQ_SetLong(reqIn, nb, ++idx);
			REQ_SetLong(reqIn, ruleId, ++idx);
			REQ_SetString (reqIn, (const char*) Ccy, ++idx);

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

long ARM_ARM_ADDMONTHS (long date,
						long nb,
						long ruleId,
						const CCString& Ccy,
						ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_ADDMONTHS ;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, XLDATE2ARMDATE(date), idx);
			REQ_SetLong(reqIn, nb, ++idx);
			REQ_SetLong(reqIn, ruleId, ++idx);
			REQ_SetString (reqIn, (const char*) Ccy, ++idx);

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


long ARM_ARM_ADDPERIOD  (long date,
						 long freq,
						 const CCString& ccy,
						 long nbPeriods,
						 long adjRuleId,
						 ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
	long idx = 0;

	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_ADDPERIOD ;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetString (reqIn, XLDATE2ARMDATE(date), idx);
			REQ_SetLong(reqIn, freq, ++idx);
			REQ_SetString(reqIn, (const char*)ccy, ++idx);
			REQ_SetLong(reqIn, nbPeriods, ++idx);
			REQ_SetLong(reqIn, adjRuleId, ++idx);

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


long ARM_ARM_INTERPOL (const VECTOR<double>& vecX,
					   const VECTOR<double>& vecY,
					   double X,
					   long interpId,
					   ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
	long idx = 0;
			
	try
	{
		/*--- parameters checking ---*/
		if(vecX.size () != vecY.size ())
		{
			result.setMsg ("ARM_ERR: vector X and vector Y must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_INTERPOL ;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong(reqIn, vecX.size(), idx);
			REQ_SetDoubleVector (reqIn, vecX, ++idx);
			REQ_SetDoubleVector (reqIn, vecY, ++idx);
			REQ_SetDouble(reqIn, X, ++idx);
			REQ_SetLong(reqIn, interpId, ++idx);

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


long ARM_DiscountPriceRefvalue(long zcId,
							   long refvalId,
							   long ccyId,
							   double startDate,
							   double endDate,
							   ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_DISCOUNTPRICEREFVAL ;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong(reqIn,zcId,idx);
			REQ_SetLong(reqIn,refvalId,++idx);
			REQ_SetLong(reqIn,ccyId,++idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), ++idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(endDate), ++idx);

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

long ARM_ARM_Price_OptUnder (long secId,
							 long modId,
							 ARM_result& result)
{
	ARM_result tmpResult;

	try
	{
		long retCode = ARM_ARM_Price(secId,modId,tmpResult);

		result.setArray(tmpResult.getDouble(),0);

		retCode = ARM_GetUnderPrice(secId,tmpResult);
		result.setArray(tmpResult.getDouble(),1);
			
		ARM_RESULT();
	}

	catch(const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}


long ARM_ARM_ClonedAndSetNotional (long secId,
								   long rId,
								   double percentRemainder,
								   ARM_result& result,
								   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
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
				ARM_REQUEST_ID = RPC_SETCLONEDUPDNOTIONAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATECLONEDUPDNOTIONAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, secId, idx);
            REQ_SetLong (reqIn, rId, ++idx);
			REQ_SetDouble (reqIn, percentRemainder, ++idx);

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


long ARM_ARM_DisplayScheduleValues (long instId,
									long valuesType,
									long recId,
									long modelId,
									const CCString& sockId,
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

			ARM_REQUEST_ID = RPC_DISP_SCHED_VALUES;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, instId, idx);
			REQ_SetLong (reqIn, valuesType, ++idx);
			REQ_SetLong (reqIn, recId, ++idx);
			REQ_SetLong (reqIn, modelId, ++idx);
			REQ_SetString (reqIn, sockId, ++idx);
			
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


long ARM_ARM_DisplayScheduleDates (long instId,
								   long datesType,
								   long recId,
								   const CCString& sockId,
								   ARM_result& result)
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

			ARM_REQUEST_ID = RPC_DISP_SCHED_DATES;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong (reqIn, instId, idx);
			REQ_SetLong (reqIn, datesType, ++idx);
			REQ_SetLong (reqIn, recId, ++idx);
			REQ_SetString (reqIn, sockId, ++idx);

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

long ARM_ARM_TRIANGULARINTERPOL(const VECTOR<double>& vecX,
								const VECTOR<double>& vecY,
								const VECTOR<double>& matZ,
								double X,
								double Y,
								ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
			
	try
	{
		/*--- parameters checking ---*/
		if (matZ.size() != (vecX.size() * vecY.size()))
		{
			result.setMsg ("ARM_ERR: check your matrix dimension");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_TRIANGULAR_INTERPOL ;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetLong(reqIn, vecX.size(), idx);
			REQ_SetDoubleVector (reqIn, vecX, ++idx);
			REQ_SetLong (reqIn, vecY.size(), ++idx);
			REQ_SetDoubleVector (reqIn, vecY, ++idx);
			REQ_SetDoubleVector (reqIn, matZ, ++idx);
			REQ_SetDouble(reqIn, X, ++idx);
			REQ_SetDouble(reqIn, Y, ++idx);

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