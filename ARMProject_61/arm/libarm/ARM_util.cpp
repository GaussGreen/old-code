#include "ARM_interglob.h"
#include "ARM_util.h"






long ARM_GetPrtyByName (long objId, const CCString& prtByName, long nbArg, ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	long idx = 0;

    
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GETPRTYBYNAME;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, objId, idx);
			REQ_SetString (reqIn, (const char*)prtByName, ++idx);
			REQ_SetLong (reqIn, nbArg, ++idx);
						

	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

	  		result.set(reqOut);
        
			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}

    catch(const CORBA::Exception& e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}



long ARM_GetVolOrRatesRange(double date1,
					        double date2,
					        const CCString& ccy,
					        const CCString& index,
					        const CCString& cvName,
					        const CCString& expiry,
					        const CCString& matu,
					        long yieldOrValId,
					        long calcModId,
					        const CCString& volType,
					        const CCString& outFile,
					        ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GETFWDRATESHISTO;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date1), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date2), ++idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);
			REQ_SetString (reqIn, (const char*)index, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, (const char*)expiry, ++idx);
			REQ_SetString (reqIn, (const char*)matu, ++idx);
			REQ_SetLong (reqIn, yieldOrValId, ++idx);
			REQ_SetLong (reqIn, calcModId, ++idx);
			REQ_SetString (reqIn, (const char*)volType, ++idx);
			REQ_SetString (reqIn, (const char*)outFile, ++idx);
			
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

	return(ARM_KO);
}



long ARM_GetHistoricalVol(double date1,
					      double date2,
					      const CCString& ccy,
					      const CCString& index,
					      const CCString& cvName,
					      const CCString& expiry,
					      const CCString& matu,
					      long yieldOrValId,
					      long calcModId,
					      const CCString& volType,
					      ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_HISTORICALVOL;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date1), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date2), ++idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);
			REQ_SetString (reqIn, (const char*)index, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, (const char*)expiry, ++idx);
			REQ_SetString (reqIn, (const char*)matu, ++idx);
			REQ_SetLong (reqIn, yieldOrValId, ++idx);
			REQ_SetLong (reqIn, calcModId, ++idx);
			REQ_SetString (reqIn, (const char*)volType, ++idx);
	
			
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

	return(ARM_KO);
}



long ARM_GetAsOfVolOrRate(double date1,
					      const CCString& ccy,
					      const CCString& index,
					      const CCString& cvName,
					      const CCString& expiry,
					      const CCString& matu,
					      long yieldOrValId,
					      long calcModId,
					      const CCString& volType,
					      ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_ASOFVOLRATE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date1), idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);
			REQ_SetString (reqIn, (const char*)index, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, (const char*)expiry, ++idx);
			REQ_SetString (reqIn, (const char*)matu, ++idx);
			REQ_SetLong (reqIn, yieldOrValId, ++idx);
			REQ_SetLong (reqIn, calcModId, ++idx);
			REQ_SetString (reqIn, (const char*)volType, ++idx);

			
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

	return(ARM_KO);
}



long ARM_ComputeHistoCorrel(double date1,
					        double date2,
					        const CCString& ccy1,
					        const CCString& index1,
					        const CCString& cvName1,
					        const CCString& expiry1,
					        const CCString& matu1,
					        long yieldOrValId1,
					        long calcModId1,
					        const CCString& volType1,
                            const CCString& ccy2,
					        const CCString& index2,
					        const CCString& cvName2,
					        const CCString& expiry2,
					        const CCString& matu2,
					        long yieldOrValId2,
					        long calcModId2,
					        const CCString& volType2,
					        ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
	long idx = 0;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_HISTOCORREL;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date1), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date2), ++idx);
			REQ_SetString (reqIn, (const char*)ccy1, ++idx);
			REQ_SetString (reqIn, (const char*)index1, ++idx);
			REQ_SetString (reqIn, (const char*)cvName1, ++idx);
			REQ_SetString (reqIn, (const char*)expiry1, ++idx);
			REQ_SetString (reqIn, (const char*)matu1, ++idx);
			REQ_SetLong (reqIn, yieldOrValId1, ++idx);
			REQ_SetLong (reqIn, calcModId1, ++idx);
			REQ_SetString (reqIn, (const char*)volType1, ++idx);
	

			REQ_SetString (reqIn, (const char*)ccy2, ++idx);
			REQ_SetString (reqIn, (const char*)index2, ++idx);
			REQ_SetString (reqIn, (const char*)cvName2, ++idx);
			REQ_SetString (reqIn, (const char*)expiry2, ++idx);
			REQ_SetString (reqIn, (const char*)matu2, ++idx);
			REQ_SetLong (reqIn, yieldOrValId2, ++idx);
			REQ_SetLong (reqIn, calcModId2, ++idx);
			REQ_SetString (reqIn, (const char*)volType2, ++idx);


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

	return(ARM_KO);
}



long ARM_GetHistoFwdRates (double date1,
					       double date2,
					       const CCString& expiry,
					       const CCString& matu,
                           const CCString& ccy,
					       const CCString& outFile,
					       ARM_result& result)
{
    long rc;


    rc = ARM_GetVolOrRatesRange(date1, date2, ccy,
					            "EURIB", "MO", 
                                expiry, matu,
					            0,
					            1,
					            "IRG",
					            outFile,
					            result);

    return(rc);
}


long ARM_GetFwdRatesMatrix (double date,
		        			const CCString& ccy,
					        ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 2;
	long idx = 0;

    
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GETFWDRATESMATRIX;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);
						

	  		/*--- CORBA Server Call ---*/

			CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

	  		result.set(reqOut);
        
			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}

    catch(const CORBA::Exception& e)
	{
		CORBA_ERR();
	}

	return ARM_KO;

}
// EOF %M%