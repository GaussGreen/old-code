#include "ARM_interglob.h"
#include "ARM_volcrv.h"




long ARM_volcurv (const VECTOR<double>& matu, const VECTOR<double>& strikes,
				  const VECTOR<double>& vols, double date,
				  long strikeType, long volType,
				  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 8;
	CCString stringObjectId;
			
	try
	{
		/*--- parameters checking ---*/
		if(((vols.size ()) != (matu.size () * strikes.size ())))
		{
			result.setMsg ("ARM_ERR: check your volatility matrix dimension");
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
				ARM_REQUEST_ID = RPC_SETVOLCURVELIN;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEVOLCURVELIN;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, matu.size (), idx);
			REQ_SetDoubleVector (reqIn, matu, ++idx);
			REQ_SetLong (reqIn, strikes.size (), ++idx);
			REQ_SetDoubleVector (reqIn, strikes, ++idx);
			REQ_SetDoubleVector (reqIn, vols, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), ++idx);
			REQ_SetLong (reqIn, strikeType, ++idx);
			REQ_SetLong (reqIn, volType, ++idx);
			
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



long ARM_GetVolFromSummit (const CCString& index, const CCString& currency, 
						   const CCString& cvName, double date, 
						   const CCString& vtype, const CCString& matuIndex,
						   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
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
				ARM_REQUEST_ID = RPC_SETVOLFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GETVOLFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)index, idx);
			REQ_SetString (reqIn, (const char*)currency, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(date), ++idx);
			REQ_SetString (reqIn, (const char*)vtype, ++idx);
			REQ_SetString (reqIn, (const char*)matuIndex, ++idx);
			
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
		CORBA_ERR()
	}

	return ARM_KO;
}




long ARM_GetVolCubeFromSummit(const CCString& index, const CCString& currency, 
						      const CCString& cvName, double date, 
						      const CCString& vtype,
							  VECTOR<CCString>& tenors,
							  const CCString&   smileOrNot,
						      ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 8;
	CCString stringObjectId;
			
	try
	{
		long nbTenors = tenors.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETVOLCUBEFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GETVOLCUBEFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString(reqIn, (const char*) index, idx);
			REQ_SetString(reqIn, (const char*) currency, ++idx);
			REQ_SetString(reqIn, (const char*) cvName, ++idx);
			REQ_SetString(reqIn, XLDATE2ARMDATE(date), ++idx);
			REQ_SetString(reqIn, (const char*) vtype, ++idx);
			REQ_SetLong (reqIn, nbTenors, ++idx);
			REQ_SetStringVector(reqIn, tenors, ++idx);
			REQ_SetString(reqIn, (const char*) smileOrNot, ++idx);
			
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
		CORBA_ERR()
	}

	return ARM_KO;
}



long ARM_ComputeVolatility(long idCurve, double matu, double strike,
						   double tenor,
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

			ARM_REQUEST_ID = RPC_COMPUTEVOLATILITY;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong(reqIn, idCurve, idx);
			REQ_SetDouble(reqIn, matu, ++idx);
			REQ_SetDouble(reqIn, strike, ++idx);
			REQ_SetDouble(reqIn, tenor, ++idx);
						
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



long ARM_volflat (double vol, double date,
				  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETVOLFLAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEVOLFLAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, vol, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), ++idx);
			
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



long ARM_ARM_GetSummitVolMatrix(const CCString& index, const CCString& currency, 
								const CCString& cvName, double date, 
								const CCString& vtype, long summitFormat,
								ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_SUMMITVOLMATRIX;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, index, idx);
			REQ_SetString (reqIn, currency, ++idx);
			REQ_SetString (reqIn, cvName, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), ++idx);
			REQ_SetString (reqIn, vtype, ++idx);
			REQ_SetLong (reqIn, summitFormat, ++idx);
						
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


long ARM_VolCube(long ATMVolId, const VECTOR<long>& volCurveIds,
						const VECTOR<double>& tenors, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETVOLCUBE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEVOLCUBE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, ATMVolId, idx);
			REQ_SetLong(reqIn, volCurveIds.size(), ++idx);
			REQ_SetLongVector(reqIn, volCurveIds, ++idx);
			REQ_SetDoubleVector (reqIn, tenors, ++idx);
			
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


long ARM_GetFXVolFromSummit(const CCString& ccy1, const CCString& ccy2, 
								 double date, const CCString& cvName, 
								 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETFXVOLFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GETFXVOLFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)ccy1, idx);
			REQ_SetString (reqIn, (const char*)ccy2, ++idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(date), ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			
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
		CORBA_ERR()
	}

	return ARM_KO;
}


long ARM_ARM_BumpVolatility(long VolId,
							double valueToBump,
							long nthLine,
							long nthCol,
							long cumulId,
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
				ARM_REQUEST_ID = RPC_SETBUMPVOLATILITY;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBUMPVOLATILITY;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, VolId, idx);
			REQ_SetDouble (reqIn, valueToBump, ++idx);
			REQ_SetLong (reqIn, nthLine, ++idx);
			REQ_SetLong (reqIn, nthCol, ++idx);
			REQ_SetLong (reqIn, cumulId, ++idx);
			
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
		CORBA_ERR()
	}

	return ARM_KO;
}


/*---- End Of File ----*/
// EOF %M%