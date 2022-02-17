#include "ARM_interglob.h"
#include "ARM_capfl.h"





long ARM_CAPFLOOR (long swapLegId,
				   long capOrFloor,
				   long strikeType,
				   double strike,
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
				ARM_REQUEST_ID = RPC_SET_CAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, capOrFloor, ++idx);
			REQ_SetLong (reqIn, strikeType, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			
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



long ARM_MATCAPFLOOR (long swapLegId,
					  double annuity,
					  double initNominal,
                      long isTRI, 
                      long capOrFloor,
					  double coeff,
                      double firstTRIstrike,
					  long minStrikesId,
					  long isDigitalPayoff,
					  double increasingCoef,
                      double maxMatDate,
				      ARM_result& result,
					  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;
	CCString stringObjectId;
			
	try
	{
        char strMaxMatuDate[50];
 
        if ( maxMatDate < 0.0 )
        {
           strcpy(strMaxMatuDate, "NULL");
        }
        else
        {
           strcpy(strMaxMatuDate, XLDATE2ARMDATE(maxMatDate));
        }


		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_MATCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_MATCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, swapLegId, idx);
            REQ_SetDouble(reqIn, annuity, ++idx);
            REQ_SetDouble(reqIn, initNominal, ++idx);
            REQ_SetLong(reqIn, isTRI, ++idx); 
			REQ_SetLong(reqIn, capOrFloor, ++idx);
            REQ_SetDouble(reqIn, coeff, ++idx);
			REQ_SetDouble(reqIn, firstTRIstrike, ++idx);
            REQ_SetLong(reqIn, minStrikesId, ++idx); 
            REQ_SetLong(reqIn, isDigitalPayoff, ++idx); 
            REQ_SetDouble(reqIn, increasingCoef, ++idx);
            REQ_SetString (reqIn, strMaxMatuDate, ++idx);

			/*--- CORBA Server Call ---*/

 	  		CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

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



long ARM_GEN_CAPFLOOR (long capOrFloor, const CCString& delivery,
					   long liborType, const CCString& currency,
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
				ARM_REQUEST_ID = RPC_SET_GEN_CAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_GEN_CAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, capOrFloor, idx);
			REQ_SetString (reqIn, (const char*)delivery, ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetString (reqIn, (const char*)currency, ++idx);

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

long ARM_GEN_CAP (const CCString& delivery, long liborType, const CCString& currency,
				  ARM_result& result, long objId)
{
	return ARM_GEN_CAPFLOOR (K_CAP, delivery, liborType, currency,
							 result, objId);
}

long ARM_GEN_FLOOR (const CCString& delivery, long liborType, const CCString& currency,
					ARM_result& result, long objId)
{
	return ARM_GEN_CAPFLOOR (K_FLOOR, delivery, liborType, currency,
		                     result, objId);
}

long ARM_LIBORCF (double startDate,
				  double endDate,
				  long isItCapOrFloor,
				  long strikeType,
				  double strike,
				  long liborType,
				  double spread,
				  long resetFreq,
				  long payFreq,
				  long currencyId,
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;
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
				ARM_REQUEST_ID = RPC_SET_LIBORCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, isItCapOrFloor, ++idx);
			REQ_SetLong (reqIn, strikeType, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, currencyId, ++idx);
			
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

long ARM_FLEXCF (long swapLegId, long isItCapOrFloor, double strike,
				 long nbEx, long exerciseType,
				 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_FLEXIBLECAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_FLEXIBLECAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, isItCapOrFloor, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, nbEx, ++idx);
			REQ_SetLong (reqIn, exerciseType, ++idx);
			
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

long ARM_LIBORFLEXCF (double startDate, double endDate, long isItCapOrFloor,
					  double strike, long nbEx, long exerciseType,
					  long liborType, double spread, long resetFreq,
					  long payFreq, long currencyId,
					  ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_LIBORFLXCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORFLXCAPFLOOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, isItCapOrFloor, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, nbEx, ++idx);
			REQ_SetLong (reqIn, exerciseType, ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, currencyId, ++idx);
			
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

long ARM_EXOFLEXCF (long swapLegId, long isItCapOrFloor, long kRefValId,
					long nbEx, long exerciseType,
					ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_EXOFLEXCF;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_EXOFLEXCF;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, isItCapOrFloor, ++idx);
			REQ_SetLong (reqIn, kRefValId, ++idx);
			REQ_SetLong (reqIn, nbEx, ++idx);
			REQ_SetLong (reqIn, exerciseType, ++idx);
			
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

long ARM_CapLetPrice (long secId, long modId, long numEx,
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

			ARM_REQUEST_ID = RPC_CAPLETPRICE;

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

long ARM_STICKY (long swapLegId,
				 long capOrFloor,
				 double strike,
				 const VECTOR<double>& spreadDates,
				 const VECTOR<double>& spreadValues,
				 long kRefValId,
				 ARM_result& result,
				 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	CCString stringObjectId;
			
	try
	{
		if(spreadDates.size () != spreadValues.size ())
		{
			result.setMsg ("ARM_ERR: spread date and value array must have same size");
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
				ARM_REQUEST_ID = RPC_SET_STICKY;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_STICKY;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> spreadDates_str;
			for(int i = 0; i < spreadDates.size (); i++)
			{
				spreadDates_str.push_back (XLDATE2ARMDATE (spreadDates[i]));
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, capOrFloor, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, spreadDates.size (), ++idx);
			REQ_SetStringVector (reqIn, spreadDates_str, ++idx);
			REQ_SetDoubleVector (reqIn, spreadValues, ++idx);
			REQ_SetLong (reqIn, kRefValId, ++idx);
			
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

long ARM_RATCHET (long swapLegId,
				  long capOrFloor,
				  double strike,
				  const VECTOR<double>& spreadDates,
				  const VECTOR<double>& spreadValues,
				  const VECTOR<double>& correlDates,
				  const VECTOR<double>& correlValues,
				  const VECTOR<double>& fwdVolsDates,
				  const VECTOR<double>& fwdVolsValues,
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 12;
	CCString stringObjectId;
			
	try
	{
		if(spreadDates.size () != spreadValues.size ())
		{
			result.setMsg ("ARM_ERR: spread date and value array must have same size");
			return ARM_KO;
		}
		if(correlDates.size () != correlValues.size ())
		{
			result.setMsg ("ARM_ERR: correlation date and value array must have same size");
			return ARM_KO;
		}
		if(fwdVolsDates.size () != fwdVolsValues.size ())
		{
			result.setMsg ("ARM_ERR: forward volatility date and value array must have same size");
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
				ARM_REQUEST_ID = RPC_SET_RATCHET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_RATCHET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> spreadDates_str;
			for(int i = 0; i < spreadDates.size (); i++)
			{
				spreadDates_str.push_back (XLDATE2ARMDATE (spreadDates[i]));
			}
			VECTOR<CCString> correlDates_str;
			for(i = 0; i < correlDates.size (); i++)
			{
				correlDates_str.push_back (XLDATE2ARMDATE (correlDates[i]));
			}
			VECTOR<CCString> fwdVolsDates_str;
			for(i = 0; i < fwdVolsDates.size (); i++)
			{
				fwdVolsDates_str.push_back (XLDATE2ARMDATE (fwdVolsDates[i]));
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, capOrFloor, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, spreadDates.size (), ++idx);
			REQ_SetStringVector (reqIn, spreadDates_str, ++idx);
			REQ_SetDoubleVector (reqIn, spreadValues, ++idx);
			REQ_SetLong (reqIn, correlDates.size (), ++idx);
			REQ_SetStringVector (reqIn, correlDates_str, ++idx);
			REQ_SetDoubleVector (reqIn, correlValues, ++idx);
			REQ_SetLong (reqIn, fwdVolsDates.size (), ++idx);
			REQ_SetStringVector (reqIn, fwdVolsDates_str, ++idx);
			REQ_SetDoubleVector (reqIn, fwdVolsValues, ++idx);
			
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


long ARM_ARM_SPREADOPTION (double startDate,
						   double endDate,
						   long capOrFloorId,
						   long strike_type,
						   double strike,
						   long liborType1Id,
						   long liborType2Id,
						   double weight1,
						   double weight2,
						   long dayCountId,
						   long resetFreqId,
						   long payFreqId,
						   long resetTimingId,
						   long payTimingId,
						   long ccyId,
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

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_SPREADOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_SPREADOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, capOrFloorId, ++idx);
			REQ_SetLong (reqIn, strike_type, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, liborType1Id, ++idx);
			REQ_SetLong (reqIn, liborType2Id, ++idx);
			REQ_SetDouble (reqIn, weight1, ++idx);
			REQ_SetDouble (reqIn, weight2, ++idx);
			REQ_SetLong (reqIn, dayCountId, ++idx);
			REQ_SetLong (reqIn, resetFreqId, ++idx);
			REQ_SetLong (reqIn, payFreqId, ++idx);
			REQ_SetLong (reqIn, resetTimingId, ++idx);
			REQ_SetLong (reqIn, payTimingId, ++idx);
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


long ARM_ARM_DIGITAL (long swapLegId,
					  long isItCapOrFloorId,
					  long strikeType,
					  double strike,
					  double spread1,
					  double spread2,
					  long payoffType,
					  double payoff,
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

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_DIGITAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_DIGITAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, isItCapOrFloorId, ++idx);
			REQ_SetLong (reqIn, strikeType, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetDouble (reqIn, spread1, ++idx);
			REQ_SetDouble (reqIn, spread2, ++idx);
			REQ_SetLong (reqIn, payoffType, ++idx);
			REQ_SetDouble (reqIn, payoff, ++idx);

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