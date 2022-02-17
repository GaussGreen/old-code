#include "ARM_interglob.h"
#include "ARM_swap.h"




long ARM_LIBORSWAP (double StartDate, double EndDate,
				    long LiborType, long ReceiveOrPay,
					double FixedRate, double Spread, 
					long CcyId, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
	long resetFreq = -1;
	long payFreq = -1;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_LIBORSWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORSWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(StartDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(EndDate), ++idx);
			REQ_SetLong (reqIn, LiborType, ++idx);
			REQ_SetLong (reqIn, ReceiveOrPay, ++idx);
			REQ_SetDouble (reqIn, FixedRate, ++idx);
			REQ_SetDouble (reqIn, Spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, CcyId, ++idx);

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



long ARM_SWAP_PRICE_TO_RATE (long sId, double Date,
				             double Price, long modId,
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

			ARM_REQUEST_ID = RPC_SWAP_PRICE_TO_RATE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, sId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(Date), ++idx);
			REQ_SetDouble (reqIn, Price, ++idx);
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



long ARM_SWAPLEG (long irIndex, double startDate, double endDate,
				  long receiveOrPay, long spreadType, double spread,
				  long ccyId, long dayCount, long resetGap,
				  CCString resetCal, CCString payCal,
				  long decompPricingFlag,
                  long nxChange, long stubRuleId, double refDate,
				  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 15;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_SWAPLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_SWAPLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, irIndex, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetLong (reqIn, spreadType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, dayCount, ++idx);
			REQ_SetLong (reqIn, resetGap, ++idx);
			REQ_SetString (reqIn, (const char*)resetCal, ++idx);
			REQ_SetString (reqIn, (const char*)payCal, ++idx);
			REQ_SetLong (reqIn, decompPricingFlag, ++idx);
	        REQ_SetLong (reqIn, nxChange, ++idx);
			REQ_SetLong (reqIn, stubRuleId, ++idx);
			if (refDate != -1.0)
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(refDate), ++idx);
			else
				REQ_SetString (reqIn, "NULL", ++idx);


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



long ARM_SWAP_RATE_TO_PRICE (long sId, double date,
				             double rate, long modId,
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

			ARM_REQUEST_ID = RPC_SWAP_RATE_TO_PRICE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong(reqIn, sId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(date), ++idx);
			REQ_SetDouble(reqIn, rate, ++idx);
			REQ_SetLong(reqIn, modId, ++idx);
			
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



long ARM_FIXEDLEG (double startDate, double endDate,
				   long receiveOrPay, 
				   long fixedRateType, double fixedRate,
				   long dayCount, long freq, long decompFreq,
				   long payTiming, long intRule, long stubRule,
				   long ccyId, CCString payCalName,
                   long nxChange, double refDate,
                   ARM_result& result,
				   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 15;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_FIXEDLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_FIXEDLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString(reqIn, XLDATE2ARMDATE(startDate), idx);
			REQ_SetString(reqIn, XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong(reqIn, receiveOrPay, ++idx);
			REQ_SetLong(reqIn, fixedRateType, ++idx);
			REQ_SetDouble(reqIn, fixedRate, ++idx);
			REQ_SetLong(reqIn, dayCount, ++idx);
			REQ_SetLong(reqIn, freq, ++idx);
			REQ_SetLong(reqIn, decompFreq, ++idx);
			REQ_SetLong(reqIn, payTiming, ++idx);
			REQ_SetLong(reqIn, intRule, ++idx);
			REQ_SetLong(reqIn, stubRule, ++idx);
			REQ_SetLong(reqIn, ccyId, ++idx);
			REQ_SetString(reqIn, (const char*) payCalName, ++idx);
            REQ_SetLong (reqIn, nxChange, ++idx);
			if (refDate != -1.0)
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(refDate), ++idx);
			else
				REQ_SetString (reqIn, "NULL", ++idx);

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



long ARM_LIBORLEG (double startDate, double endDate, long liborType,
				   long receiveOrPay, long spreadType, double spread,
				   long resetFreq, long payFreq,
				   long resetTiming, long payTiming,
				   long ccyId, long intRuleId, long resetGap,
				   CCString resetCal, CCString payCal, 
                   long decompPricingFlag,
                   long nxChange, long stubRuleId, double refDate,
				   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 19;

	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_LIBORLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetLong (reqIn, spreadType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, resetTiming, ++idx);
			REQ_SetLong (reqIn, payTiming, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, intRuleId, ++idx);
			REQ_SetLong (reqIn, resetGap, ++idx);
			REQ_SetString (reqIn, (const char*)resetCal, ++idx);
			REQ_SetString (reqIn, (const char*)payCal, ++idx);
			REQ_SetLong (reqIn, decompPricingFlag, ++idx);
            REQ_SetLong (reqIn, nxChange, ++idx);
			REQ_SetLong (reqIn, stubRuleId, ++idx);
			if (refDate != -1.0)
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(refDate), ++idx);
			else
				REQ_SetString (reqIn, "NULL", ++idx);


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



long ARM_TMLEG (long tmIxType, double startDate, double endDate,
				long receiveOrPay, double spread,
				ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_TMLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_TMLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, tmIxType, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			
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



long ARM_SwapFromSummit (const CCString& swapId, long zcId,
						 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_CREATED_SWAP_FROM_SUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GET_SWAP_FROM_SUMMIT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)swapId, idx);
			REQ_SetLong (reqIn, zcId, idx);
	
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



long ARM_IMPLIEDSPREAD (long swapId, long modelId, double price,
					    long leg1Or2, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_IMPLIED_SPREAD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, swapId, idx);
			REQ_SetLong (reqIn, modelId, ++idx);
			REQ_SetDouble (reqIn, price, ++idx);
			REQ_SetLong (reqIn, leg1Or2, ++idx);
			
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



long ARM_IMPLIEDSPREADWITHMODELS (long swapId, long modelId1, long modelId2,
								  double price,
							      long leg1Or2, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_IMPLIED_SPREAD_WITH2MODELS;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, swapId, idx);
			REQ_SetLong (reqIn, modelId1, ++idx);
			REQ_SetLong (reqIn, modelId2, ++idx);
			REQ_SetDouble (reqIn, price, ++idx);
			REQ_SetLong (reqIn, leg1Or2, ++idx);
			
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



long ARM_GEN_SWAP (const CCString& delivery, long liborTypeId,
				   const CCString& ccy, ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_GEN_SWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_GEN_SWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)delivery, idx);
			REQ_SetLong (reqIn, liborTypeId, ++idx);
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



long ARM_CMTLEG(double startDate, double endDate, long cmtTypeId, long bondCouponFreq,
				long bondDayCount, long receiveOrPay, double spread,
				long yieldDecompFreq, long swapLegDayCount, long intRule,
				long resetGap, long resetFreq,
				double ntlAmount, long ccyId, 
                long resetTiming,
                ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 15;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SET_CMTLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CMTLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

		
			
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, cmtTypeId, ++idx);
			REQ_SetLong (reqIn, bondCouponFreq, ++idx);
			REQ_SetLong (reqIn, bondDayCount, ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, yieldDecompFreq, ++idx);
			REQ_SetLong (reqIn, swapLegDayCount, ++idx);
			REQ_SetLong (reqIn, intRule, ++idx);
			REQ_SetLong (reqIn, resetGap, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetDouble (reqIn, ntlAmount, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
	        REQ_SetLong (reqIn, resetTiming, ++idx);


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



long ARM_CMSLEG (double startDate, double endDate, long cmsTypeId,
				 long receiveOrPay, double spread,
				 long yieldDecompFreq, long swapLegDayCount, long resetFreq,
				 long intRule, long ccyId, long resetTiming, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 11;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_CMSLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CMSLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, cmsTypeId, ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, yieldDecompFreq, ++idx);
			REQ_SetLong (reqIn, swapLegDayCount, ++idx);
			REQ_SetLong (reqIn, intRule, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, resetTiming, ++idx);
	
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



long ARM_SWAP (long swapLeg1, long swapLeg2, long minPay,
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
				ARM_REQUEST_ID = RPC_SET_SWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_SWAP;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLeg1, idx);
			REQ_SetLong (reqIn, swapLeg2, ++idx);
			REQ_SetLong (reqIn, minPay, ++idx);
			
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



long ARM_SWAP_WITH_NOTIONNAL (long swapLeg1, long swapLeg2, long notId,
							  long minPay, ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_SWAP_WITH_NOTIONNAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_SWAP_WITH_NOTIONNAL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapLeg1, idx);
			REQ_SetLong (reqIn, swapLeg2, ++idx);
			REQ_SetLong (reqIn, notId, ++idx);
			REQ_SetLong (reqIn, minPay, ++idx);
			
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



long ARM_FIXRATES (double swapLegId,
				   VECTOR<double>& rate,
				   ARM_result& result,
				   long objId)
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

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETSET_FIXEDRATES;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATESET_FIXEDRATES;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, swapLegId, idx);
			REQ_SetLong (reqIn, rate.size (), ++idx);
			REQ_SetDoubleVector (reqIn, rate, ++idx);
	

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



long ARM_LiborAssetSwapMargin (long modelId, double startDate, double endDate,
							   double fixedRate, long fixReceiveOrPayId, 
							   long fixDayCountId, long fixFrequencyId,
							   long fixDecompFrequencyId, long fixPayTimingId, long fixIntRuleId,
							   long liborTypeId, double spread, long floatResetFreqId,
							   long floatPayFreqId, long assetGap, long vFlag,
							   double price, const CCString& discountCcy,
							   double redemptionPrice, double supplFee,
							   long solve, const CCString& id,
							   ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 22;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_LIBOR_ASSET_SWAP_MARGIN;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, modelId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetDouble (reqIn, fixedRate, ++idx);
			REQ_SetLong (reqIn, fixReceiveOrPayId, ++idx);
			REQ_SetLong (reqIn, fixDayCountId, ++idx);
			REQ_SetLong (reqIn, fixFrequencyId, ++idx);
			REQ_SetLong (reqIn, fixDecompFrequencyId, ++idx);
			REQ_SetLong (reqIn, fixPayTimingId, ++idx);
			REQ_SetLong (reqIn, fixIntRuleId, ++idx);
			REQ_SetLong (reqIn, liborTypeId, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, floatResetFreqId, ++idx);
			REQ_SetLong (reqIn, floatPayFreqId, ++idx);
			REQ_SetLong (reqIn, assetGap, ++idx);
			REQ_SetDouble (reqIn, price, ++idx);
			REQ_SetString (reqIn, discountCcy, ++idx);
			REQ_SetLong (reqIn, vFlag, ++idx);
			REQ_SetDouble (reqIn, redemptionPrice, ++idx);
			REQ_SetDouble (reqIn, supplFee, ++idx);
			REQ_SetLong (reqIn, solve, ++idx);
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



long ARM_RESTRIKABLELEG(double startDate,
						double endDate,
						long receiveOrPay,
						long payIndexId,
				        long payFreq,
						double spread, 
                        long refIndexId,
						long resetFreq,
				        long resetTiming,
						long payTiming, 
                        long stubRule,
                        double range,
						long rangeSpec,
				        long ccyId,
						long MCFreq,
						long MCInterp,
                        long decompPricingFlag,
                        ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 17;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_RESTRIKABLELEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_RESTRIKABLELEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}


			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetLong (reqIn, payIndexId, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, refIndexId, ++idx);
            REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, resetTiming, ++idx);
			REQ_SetLong (reqIn, payTiming, ++idx);
			REQ_SetLong (reqIn, stubRule, ++idx);
            REQ_SetDouble (reqIn, range, ++idx);
            REQ_SetLong (reqIn, rangeSpec, ++idx);
			REQ_SetLong (reqIn, ccyId, ++idx);
			REQ_SetLong (reqIn, MCFreq, ++idx);
			REQ_SetLong (reqIn, MCInterp, ++idx);
            REQ_SetLong (reqIn, decompPricingFlag, ++idx);

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



long ARM_CORRIDORLEG(double startDate,
					 double endDate,
					 long receiveOrPay,
					 long payIndexId,
					 long payFreq,
					 long spreadType,
					 double spread,
					 long refIndexId,
					 long resetFreq,
					 long paidRateResetTiming,
					 long refRateResetTiming,
					 long stubRule,
					 long levelDownId,
					 long downSpec,
					 long levelUpId,
					 long upSpec,
					 long ccyId,
					 long MCFreq,
					 long MCInterp,
					 long LDPricingMethod,
                     long decompPricingFlag,
					 ARM_result& result,
					 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 21;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_CORRIDORLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CORRIDORLEG;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}


			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetLong (reqIn, payIndexId, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, spreadType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, refIndexId, ++idx);
            REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, paidRateResetTiming, ++idx);
			REQ_SetLong (reqIn, refRateResetTiming, ++idx);
			REQ_SetLong (reqIn, stubRule, ++idx);

            REQ_SetLong (reqIn, levelDownId, ++idx);
            REQ_SetLong (reqIn, downSpec, ++idx);

            REQ_SetLong (reqIn, levelUpId, ++idx);
            REQ_SetLong (reqIn, upSpec, ++idx);

            REQ_SetLong (reqIn, ccyId, ++idx);
            REQ_SetLong (reqIn, MCFreq, ++idx);
            REQ_SetLong (reqIn, MCInterp, ++idx);
            REQ_SetLong (reqIn, LDPricingMethod, ++idx);

            REQ_SetLong (reqIn, decompPricingFlag, ++idx);

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


long ARM_REVERSEFLOAT(double startDate,
					  double endDate,
					  double NotionalRatio,
					  long receiveOrPay,
					  long liborType,
					  const CCString& fltccy,
					  long fltSpreadsId,
					  long fltDayCount,
					  long fxCouponsId,
					  long extraSpreadFixId,
					  const CCString& reverseCcy,
					  long payFreq,
					  long reverseIndexId,
					  long fxDayCount,
					  double multiplier,
					  double couponFloor,
					  double stubDate,
					  ARM_result& result,
					  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 17;
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_REVERSEFLOAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_REVERSEFLOAT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}



			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
            REQ_SetDouble (reqIn, NotionalRatio, ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
            REQ_SetString (reqIn, (const char*) fltccy, ++idx);
            REQ_SetLong (reqIn, fltSpreadsId, ++idx);
            REQ_SetLong (reqIn, fltDayCount, ++idx);
            REQ_SetLong (reqIn, fxCouponsId, ++idx);
            REQ_SetLong (reqIn, extraSpreadFixId, ++idx);
            REQ_SetString (reqIn, (const char*) reverseCcy, ++idx);
            REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, reverseIndexId, ++idx);
            REQ_SetLong (reqIn, fxDayCount, ++idx);
            REQ_SetDouble (reqIn, multiplier, ++idx);
            REQ_SetDouble (reqIn, couponFloor, ++idx);

            if ( stubDate > 0 ) 
            {
                REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(stubDate), ++idx);
            }
            else
            {
                REQ_SetString (reqIn, (const char*) "01/01/1970", ++idx);
            }

			
			
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


long ARM_IMPLIED_RANGE (long restrikableId,
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

			ARM_REQUEST_ID = RPC_IMPLIED_RANGE;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, restrikableId, idx);
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



long ARM_ARM_GenAmortization (long swaplegId,
							  long amortMethodId,
							  long amortFrequency,
							  double amortAmount,
							  long daycountId,
							  double legNotional,
							  double amortRate,
							  double reducedMaturity,
							  long modelId,
							  double percentOfRemainder,
							  ARM_result& result,
							  long objId)
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

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_GEN_AMORT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_GEN_AMORT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}
			
			REQ_SetLong (reqIn, swaplegId, idx);
			REQ_SetLong (reqIn, amortMethodId, ++idx);
			REQ_SetLong (reqIn, amortFrequency, ++idx);
			REQ_SetDouble (reqIn, amortAmount, ++idx);
			REQ_SetLong (reqIn, daycountId, ++idx);
			REQ_SetDouble (reqIn, legNotional, ++idx);
			REQ_SetDouble (reqIn, amortRate, ++idx);
			REQ_SetDouble (reqIn, reducedMaturity, ++idx);
			REQ_SetLong (reqIn, modelId, ++idx);
            REQ_SetDouble (reqIn, percentOfRemainder, ++idx);

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


long ARM_ARM_CUSTOMFSTCPN (long swlegId,
						   double date,
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

			ARM_REQUEST_ID = RPC_CUST_FIRST_PERIOD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, swlegId, idx);
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

/*---- End Of File ----*/

// EOF %M%
