#include "ARM_interglob.h"
#include "ARM_rev.h"





long ARM_REVERSE (long structSwapLegId,
				  long classSwapLegId,
				  long ReceiveOrPay,
				  long couponId,
				  long exeId,
				  long redempId,
				  long classRedempId,
                  double dualDate,
                  double dualStrike,
                  const CCString& dualFlag,
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;
	CCString stringObjectId;
	long resetFreq = -1;
	long payFreq = -1;
			
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
				ARM_REQUEST_ID = RPC_SET_REVERSE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_REVERSE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, structSwapLegId, idx);
			REQ_SetLong (reqIn, classSwapLegId, ++idx);
			REQ_SetLong (reqIn, ReceiveOrPay, ++idx);
			REQ_SetLong (reqIn, couponId, ++idx);
			REQ_SetLong (reqIn, exeId, ++idx);
			REQ_SetLong (reqIn, redempId, ++idx);
			REQ_SetLong (reqIn, classRedempId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(dualDate), ++idx);
			REQ_SetDouble (reqIn, dualStrike, ++idx);
            REQ_SetString (reqIn, dualFlag, ++idx);


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


long ARM_REVERSE_CALENDAR (long structSwapLegId,
				  long classSwapLegId,
				  long ReceiveOrPay,
				  long couponId,
				  long exeId,
				  long redempId,
				  long classRedempId,
                  double dualDate,
                  double dualStrike,
                  const CCString& dualFlag,
				  const VECTOR<double>& dStartDates, 
				  const VECTOR<double>& dEndDates, 
				  const VECTOR<double>& dFixingDates, 
				  const VECTOR<double>& dPaymentDates, 
				  const VECTOR<double>& fStartDates, 
				  const VECTOR<double>& fEndDates, 
				  const VECTOR<double>& fFixingDates, 
				  const VECTOR<double>& fPaymentDates, 
				  ARM_result& result,
				  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 20;
	CCString stringObjectId;
	long resetFreq = -1;
	long payFreq = -1;
			
	try
	{
		if (
			(!((dStartDates.size () == dEndDates.size ()) 
             && (dEndDates.size () == dFixingDates.size ()) 
             && (dFixingDates.size () == dPaymentDates.size ())
             ))
            ||
		    (!((fStartDates.size () == fEndDates.size ()) 
             && (fEndDates.size () == fFixingDates.size ()) 
             && (fFixingDates.size () == fPaymentDates.size ())
             ))
           )	
		{
			result.setMsg ("ARM_ERR: startDate, endDate, fixingDate and paymentDates arrays must have the same size");
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
				ARM_REQUEST_ID = RPC_SET_REVERSE_CALENDAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_REVERSE_CALENDAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> dStartDates_str;
			VECTOR<CCString> dEndDates_str;
			VECTOR<CCString> dFixingDates_str;
			VECTOR<CCString> dPaymentDates_str;
			VECTOR<CCString> fStartDates_str;
			VECTOR<CCString> fEndDates_str;
			VECTOR<CCString> fFixingDates_str;
			VECTOR<CCString> fPaymentDates_str;
			
			for(int i = 0; i < dStartDates.size (); i++)
			{
				dStartDates_str.push_back (XLDATE2ARMDATE (dStartDates[i]));
				dEndDates_str.push_back (XLDATE2ARMDATE (dEndDates[i]));
				dFixingDates_str.push_back (XLDATE2ARMDATE (dFixingDates[i]));
				dPaymentDates_str.push_back (XLDATE2ARMDATE (dPaymentDates[i]));
			}
			for(i = 0; i < fStartDates.size (); i++)
			{
				fStartDates_str.push_back (XLDATE2ARMDATE (fStartDates[i]));
				fEndDates_str.push_back (XLDATE2ARMDATE (fEndDates[i]));
				fFixingDates_str.push_back (XLDATE2ARMDATE (fFixingDates[i]));
				fPaymentDates_str.push_back (XLDATE2ARMDATE (fPaymentDates[i]));
			}

			REQ_SetLong (reqIn, structSwapLegId, idx);
			REQ_SetLong (reqIn, classSwapLegId, ++idx);
			REQ_SetLong (reqIn, ReceiveOrPay, ++idx);
			REQ_SetLong (reqIn, couponId, ++idx);
			REQ_SetLong (reqIn, exeId, ++idx);
			REQ_SetLong (reqIn, redempId, ++idx);
			REQ_SetLong (reqIn, classRedempId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(dualDate), ++idx);
			REQ_SetDouble (reqIn, dualStrike, ++idx);
            REQ_SetString (reqIn, dualFlag, ++idx);
			REQ_SetStringVector (reqIn, dStartDates_str, ++idx);
			REQ_SetStringVector (reqIn, dEndDates_str, ++idx);
			REQ_SetStringVector (reqIn, dFixingDates_str, ++idx);
			REQ_SetStringVector (reqIn, dPaymentDates_str, ++idx);
			REQ_SetLong (reqIn, dStartDates.size(), ++idx);
			REQ_SetStringVector (reqIn, fStartDates_str, ++idx);
			REQ_SetStringVector (reqIn, fEndDates_str, ++idx);
			REQ_SetStringVector (reqIn, fFixingDates_str, ++idx);
			REQ_SetStringVector (reqIn, fPaymentDates_str, ++idx);
			REQ_SetLong (reqIn, fStartDates.size(), ++idx);

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


long ARM_REVERSECOUPON (double strike,
						double power,
						double callput,
						long size,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	CCString stringObjectId;
	long resetFreq = -1;
	long payFreq = -1;
			
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
				ARM_REQUEST_ID = RPC_SET_REVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_REVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, strike, idx);
			REQ_SetDouble (reqIn, power, ++idx);
			REQ_SetDouble (reqIn, callput, ++idx);
			REQ_SetLong (reqIn, size, ++idx);
			
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


/* Old version
long ARM_STRUCTREVERSECOUPON (const VECTOR<double>& strike,
							  const VECTOR<double>& power,
							  const VECTOR<double>& callput,
							  ARM_result& result,
							  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 4;
	CCString stringObjectId;
	long resetFreq = -1;
	long payFreq = -1;
			
	try
	{
		if(!((strike.size () == power.size ()) && (power.size () == callput.size ())))
		{
			result.setMsg ("ARM_ERR: strike, rate and callput array must have same size");
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
				ARM_REQUEST_ID = RPC_SET_STRUCTREVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_STRUCTREVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDoubleVector (reqIn, strike, idx);
			REQ_SetDoubleVector (reqIn, power, ++idx);
			REQ_SetDoubleVector (reqIn, callput, ++idx);
			REQ_SetLong (reqIn, strike.size (), ++idx);
			
			*--- CORBA Server Call ---*

 	  		CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

			result.set (reqOut);

			*--- requests freeing ---*
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
----*/



long ARM_STRUCTREVERSECOUPON(const VECTOR<double>& date,
							 const VECTOR<double>& strike,
							 const VECTOR<double>& power,
							 const VECTOR<double>& callput,
							 double xo,
                             const VECTOR<double>& floor,
                             const VECTOR<double>& cap,
							 ARM_result& result,
							 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 8;
	CCString stringObjectId;



	try
	{
		if (!((date.size () == strike.size ()) 
             && (strike.size () == power.size ()) 
             && (power.size () == callput.size ())
             && (callput.size () == floor.size ())
             && (floor.size () == cap.size ())
             )
           )
		{
			result.setMsg ("ARM_ERR: date, strike, rate, callput, cap and floor arrays must have the same size");
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
				ARM_REQUEST_ID = RPC_SET_STRUCTREVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_STRUCTREVERSECOUPON;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> date_str;
			for(int i = 0; i < date.size (); i++)
			{
				date_str.push_back (XLDATE2ARMDATE (date[i]));
			}

			REQ_SetStringVector (reqIn, date_str, idx);
			REQ_SetDoubleVector (reqIn, strike, ++idx);
			REQ_SetDoubleVector (reqIn, power, ++idx);
			REQ_SetDoubleVector (reqIn, callput, ++idx);
			REQ_SetDouble (reqIn, xo, ++idx);
			REQ_SetLong (reqIn, strike.size (), ++idx);
            REQ_SetDoubleVector (reqIn, floor, ++idx);
            REQ_SetDoubleVector (reqIn, cap, ++idx);
			
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



long ARM_STRUCTREVERSECOUPON2 (const VECTOR<double>& date,
							   const VECTOR<double>& strike,
							   const VECTOR<double>& power,
							   const VECTOR<double>& callput,
							   double xo,
							   ARM_result& result,
							   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	CCString stringObjectId;
	long resetFreq = -1;
	long payFreq = -1;
			
	try
	{
		if(!((date.size () == strike.size ()) && (strike.size () == power.size ()) && (power.size () == callput.size ())))
		{
			result.setMsg ("ARM_ERR: date, strike, rate and callput array must have same size");
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
				ARM_REQUEST_ID = RPC_SET_STRUCTREVERSECOUPON2;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_STRUCTREVERSECOUPON2;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> date_str;
			for(int i = 0; i < date.size (); i++)
			{
				date_str.push_back (XLDATE2ARMDATE (date[i]));
			}

			REQ_SetStringVector (reqIn, date_str, idx);
			REQ_SetDoubleVector (reqIn, strike, ++idx);
			REQ_SetDoubleVector (reqIn, power, ++idx);
			REQ_SetDoubleVector (reqIn, callput, ++idx);
			REQ_SetDouble (reqIn, xo, ++idx);
			REQ_SetLong (reqIn, strike.size (), ++idx);
			
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