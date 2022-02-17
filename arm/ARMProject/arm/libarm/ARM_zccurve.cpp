#include "ARM_interglob.h"
#include "ARM_zccurve.h"





long ARM_zclint (const VECTOR<double>& matu,
				 const VECTOR<double>& rate, 
                 long meth,
				 double aDate,
				 const CCString& ccy,
				 long interpMethId,
				 ARM_result& result,
				 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZEROCURVELIN;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZEROCURVELIN;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, matu.size (), idx);
			REQ_SetDoubleVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(aDate), ++idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);
			REQ_SetLong (reqIn, interpMethId, ++idx);

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



long ARM_GetZCFromSummit (const CCString& index, const CCString& currency, 
                          const CCString& cvName,
                          double aSdate, ARM_result& result, long objId)
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

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETCREATEDZEROCURVEFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GETZCFROMSUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)index, idx);
			REQ_SetString (reqIn, (const char*)currency, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(aSdate), ++idx);

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



long ARM_GetInitialCurveFromSummit (const CCString& index,
									const CCString& currency, 
								    const CCString& cvName,
								    double aSdate,
									long AdjOrNotId,
									const CCString& outFile,
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

			ARM_REQUEST_ID = RPC_GETINITIALCURVEFROMSUMMIT;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetString (reqIn, (const char*)index, idx);
			REQ_SetString (reqIn, (const char*)currency, ++idx);
			REQ_SetString (reqIn, (const char*)cvName, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(aSdate), ++idx);
			REQ_SetLong (reqIn, AdjOrNotId, ++idx);
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

	return ARM_KO;
}



long ARM_DiscountPrice (long idCurve, double matu, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_DISCOUNTPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, idCurve, idx);
			REQ_SetDouble (reqIn, matu, ++idx);
			
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



long ARM_DiscountYield (long idCurve, double matu, long meth, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_DISCOUNTYIELD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, idCurve, idx);
			REQ_SetDouble (reqIn, matu, ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
						
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



long ARM_CreateZCSwapFutInt (double Date, VECTOR<CCString>& matu, 
                             VECTOR<double>& rate,
							 long MMVsFut, long SwapVsFut, long Raw, 
                             long interp, const CCString& Ccy,
							 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCSWAPFUTINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCSWAPFUTINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetLong (reqIn, SwapVsFut, ++idx);
			REQ_SetLong (reqIn, Raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);

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



long ARM_CreateZCSwapInt (double Date, VECTOR<CCString>& matu, 
                          VECTOR<double>& rate,
						  long MMVsFut, long SwapVsFut, long Raw, 
                          long interp, const CCString& Ccy,
						  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCSWAPINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCSWAPINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetLong (reqIn, SwapVsFut, ++idx);
			REQ_SetLong (reqIn, Raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);

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



long ARM_C_CreateZCSwapInt (double Date, 
							CCString* C_matu, long matu_size,
                            double* C_rate, long rate_size,
						    long MMVsFut, long SwapVsFut, long Raw, 
                            long interp, const CCString& Ccy,
						    ARM_result& result, long objId)
{
	VECTOR<double> rate;
	VECTOR<CCString> matu;

	for(int i = 0; i < rate_size; i++)
	{
		rate.push_back (C_rate[i]);
	}
	for(i = 0; i < matu_size; i++)
	{
		matu.push_back (C_matu[i]);
	}

	return ARM_CreateZCSwapInt (Date, matu, rate, MMVsFut, SwapVsFut, Raw, interp, Ccy, result, objId);
}



long ARM_ForwardYield (long idCurve, double matu1, double matu2, 
                       long meth, long adjId, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_FORWARDYIELD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, idCurve, idx);
			REQ_SetDouble (reqIn, matu1, ++idx);
			REQ_SetDouble (reqIn, matu2, ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
			REQ_SetLong (reqIn, adjId, ++idx);
		
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



long ARM_ForwardPrice (long idCurve, double matu1, double matu2, 
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

			ARM_REQUEST_ID = RPC_FORWARDPRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, idCurve, idx);
			REQ_SetDouble (reqIn, matu1, ++idx);
			REQ_SetDouble (reqIn, matu2, ++idx);
						
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



long ARM_zcvsk (const VECTOR<double>& param, double date, 
                ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;

	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZEROCURVEVSK;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZEROCURVEVSK;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, param.size (), idx);
			REQ_SetDoubleVector (reqIn, param, ++idx);
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



long ARM_zcflat (const double zeroFlat, double date, ARM_result& result, 
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
				ARM_REQUEST_ID = RPC_SETZEROFLAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZEROFLAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, zeroFlat, idx);
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



long ARM_zcspl (const VECTOR<double>& param, double date, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZEROCURVESPLI;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZEROCURVESPLI;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, param.size (), idx);
			REQ_SetDoubleVector (reqIn, param, ++idx);
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





long ARM_CreateZCTAMInt (double date, VECTOR<CCString>& matu, 
					     VECTOR<double>& rate, double mean_rates, 
						 long raw, long interp, 
						 long lastBucketInt, const CCString& Ccy,
						 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCTAMINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCTAMINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetDouble (reqIn, mean_rates, ++idx);
			REQ_SetLong (reqIn, raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetLong (reqIn, lastBucketInt, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);

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



long ARM_CreateZCCashInt (double date, VECTOR<CCString>& matu,
						  VECTOR<double>& rate, VECTOR<long>& bondsId,
						  VECTOR<double>& yields, long MMVsFut,
						  const CCString& Ccy,
						  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		if(bondsId.size () != yields.size ())
		{
			result.setMsg ("ARM_ERR: bonds IDs and yields must have same size");
			return ARM_KO;
		}

		long MR_nblines = matu.size ();
		long BONDS_nblines = bondsId.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCCASHINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCCASHINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);

			REQ_SetLong (reqIn, MR_nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);

			REQ_SetLong (reqIn, BONDS_nblines, ++idx);
			REQ_SetLongVector (reqIn, bondsId, ++idx);
			REQ_SetDoubleVector (reqIn, yields, ++idx);

			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);

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



long ARM_zcsplicub (VECTOR<double>& matu, VECTOR<double>& rate,
			        long meth, double date, long lastBucket,
				    ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETSPLICUBCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATESPLICUBCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, matu.size (), idx);
			REQ_SetDoubleVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), ++idx);
			REQ_SetLong (reqIn, lastBucket, ++idx);
			
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

long ARM_zccubdiff (VECTOR<double>& matu, VECTOR<double>& rate,
					long meth, double date,
					ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
			
	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETCUBDIFFCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATECUBDIFFCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, matu.size (), idx);
			REQ_SetDoubleVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
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

long ARM_zcswapcubdiff (double date, VECTOR<CCString>& matu, VECTOR<double>& rate,
						long mmVsFut, long swapVsFut, long raw, long interp, const CCString& ccy,
						ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCSWAPCUBDIFF;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCSWAPCUBDIFF;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, mmVsFut, ++idx);
			REQ_SetLong (reqIn, swapVsFut, ++idx);
			REQ_SetLong (reqIn, raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)ccy, ++idx);

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

long ARM_zcspreaded (long zcSprId, long zcInitId, double date,
				     long MMFreq, long SwapFreq, long ccyId,
					 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;

	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETSPREADCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATESPREADCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(date), idx);
			REQ_SetLong (reqIn, zcSprId, ++idx);
			REQ_SetLong (reqIn, zcInitId, ++idx);
			REQ_SetLong (reqIn, MMFreq, ++idx);
			REQ_SetLong (reqIn, SwapFreq, ++idx);
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



long ARM_CreateZCSwapIntSmooth (double Date, VECTOR<CCString>& matu, 
								VECTOR<double>& rate,
								long MMVsFut, long SwapVsFut, long Raw, 
								long interp, const CCString& Ccy,
								double lambda, long prec,
								ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCSWAPINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCSWAPINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetLong (reqIn, SwapVsFut, ++idx);
			REQ_SetLong (reqIn, Raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);
			REQ_SetDouble (reqIn, lambda, ++idx);
			REQ_SetLong (reqIn, prec, ++idx);

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

long ARM_CreateZCSwapFutIntSmooth (double Date, VECTOR<CCString>& matu, 
								   VECTOR<double>& rate,
								   long MMVsFut, long SwapVsFut, long Raw, 
								   long interp, const CCString& Ccy,
								   double lambda, long prec,
								   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCSWAPFUTINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCSWAPFUTINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetLong (reqIn, SwapVsFut, ++idx);
			REQ_SetLong (reqIn, Raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);
			REQ_SetDouble (reqIn, lambda, ++idx);
			REQ_SetLong (reqIn, prec, ++idx);

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



long ARM_ZCINTSMOOTH (const VECTOR<double>& matu, const VECTOR<double>& rate, 
					  double aDate, long meth, double lambda, long prec,
					  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETZCINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEZCINTSMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, matu.size (), idx);
			REQ_SetDoubleVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(aDate), ++idx);
			REQ_SetLong (reqIn, meth, ++idx);
			REQ_SetDouble (reqIn, lambda, ++idx);
			REQ_SetLong (reqIn, prec, ++idx);
			
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

long ARM_ZCINTSMOOTH (long inCvId, double lambda, long prec,
					  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;

	try
	{

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_TRANS2SMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_TRANS2SMOOTH;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, inCvId, idx);
			REQ_SetDouble (reqIn, lambda, ++idx);
			REQ_SetLong (reqIn, prec, ++idx);
			
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

long ARM_CreateTOYNYZCSwapInt (double Date,
							   VECTOR<CCString>& matu,
							   VECTOR<double>& rate,
							   long MMVsFut,
							   long SwapVsFut,
							   long Raw,
							   long interp,
							   const CCString& Ccy,
							   long frq,
							   ARM_result& result,
							   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != rate.size ())
		{
			result.setMsg ("ARM_ERR: maturities and rates must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETTOYZCSWAPINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATETOYZCSWAPINT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
			REQ_SetLong (reqIn, MMVsFut, ++idx);
			REQ_SetLong (reqIn, SwapVsFut, ++idx);
			REQ_SetLong (reqIn, Raw, ++idx);
			REQ_SetLong (reqIn, interp, ++idx);
			REQ_SetString (reqIn, (const char*)Ccy, ++idx);
			REQ_SetLong (reqIn, frq, ++idx);

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


long ARM_ARM_BumpCurve (long cvId,
						VECTOR<CCString>& matu,
						VECTOR<double>& epsilon,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;

	try
	{
		/*--- parameters checking ---*/
		if(matu.size () != epsilon.size ())
		{
			result.setMsg ("ARM_ERR: maturities and bump must have same size");
			return ARM_KO;
		}

		long nblines = matu.size ();

		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_BUMPCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_BUMPCURVE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, cvId, idx);
			REQ_SetLong (reqIn, nblines, ++idx);
			REQ_SetStringVector (reqIn, matu, ++idx);
			REQ_SetDoubleVector (reqIn, epsilon, ++idx);

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