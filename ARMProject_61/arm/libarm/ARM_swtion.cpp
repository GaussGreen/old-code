#include "ARM_interglob.h"
#include "ARM_swtion.h"



long ARM_SWAPTION (long swapId,
			       long isRecOrPay,
				   double strike,
				   double maturity,
				   long exerciseType,
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
				ARM_REQUEST_ID = RPC_SET_SWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_SWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapId, idx);
			REQ_SetLong (reqIn, isRecOrPay, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(maturity), ++idx);
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

long ARM_GEN_SWAPTION (double swapTerm,
			           double optionExpiry,
				       long liborType,
				       long isRecOrPay,
				       const CCString& ccy,
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
				ARM_REQUEST_ID = RPC_SET_GEN_SWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_GEN_SWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(swapTerm), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(optionExpiry), ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetLong (reqIn, isRecOrPay, ++idx);
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

long ARM_LIBORSWAPTION (double startDate, double endDate,
				        long receiveOrPay, double strike,
						double maturity, long liborType,
						double spread, long exerciseType,
						long resetFreq, long payFreq,
					    long CcyId, ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SET_LIBORSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_LIBORSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(maturity), ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
			REQ_SetLong (reqIn, exerciseType, ++idx);
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

long ARM_EXOSWAPTION (long swapId,
			          long isRecOrPay,
				      long xStyleId,
				      long kRefValId,
				      double swapYearTerm,
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
				ARM_REQUEST_ID = RPC_SET_EXOSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_EXOSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapId, idx);
			REQ_SetLong (reqIn, isRecOrPay, ++idx);
			REQ_SetLong (reqIn, xStyleId, ++idx);
			REQ_SetLong (reqIn, kRefValId, ++idx);
			REQ_SetDouble (reqIn, swapYearTerm, ++idx);
						
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



long ARM_EXOCFSWAPTION (long swapId,
					    long isRecOrPay,
						long isCapOrFloor,
						long xStyleId,
						long kSptionRefValId,
						long kCFloorRefValId,
						double cFloorPosition,
                        long IsBarrierCF,
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

			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SET_CFEXOSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CFEXOSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, swapId, idx);
			REQ_SetLong (reqIn, isRecOrPay, ++idx);
			REQ_SetLong (reqIn, isCapOrFloor, ++idx);
			REQ_SetLong (reqIn, xStyleId, ++idx);
			REQ_SetLong (reqIn, kSptionRefValId, ++idx);
			REQ_SetLong (reqIn, kCFloorRefValId, ++idx);
			REQ_SetDouble (reqIn, cFloorPosition, ++idx);
		    REQ_SetLong (reqIn, IsBarrierCF, ++idx);


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


                     
                     
                     
                     
                     
                     




long ARM_VARFIXSWAPTION (double startDate, double endDate,
                         long spreadsId, long exStyleId,
				        long receiveOrPay, double strike,
						double maturity, long liborType,
						double spread, 
						long resetFreq, long payFreq,
					    long ccyId, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 12;
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
				ARM_REQUEST_ID = RPC_SET_VARFIX_SWOPT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_VARFIX_SWOPT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}



			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
            REQ_SetLong (reqIn, spreadsId, ++idx);
            REQ_SetLong (reqIn, exStyleId, ++idx);
			REQ_SetLong (reqIn, receiveOrPay, ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(maturity), ++idx);
			REQ_SetLong (reqIn, liborType, ++idx);
			REQ_SetDouble (reqIn, spread, ++idx);
			REQ_SetLong (reqIn, resetFreq, ++idx);
			REQ_SetLong (reqIn, payFreq, ++idx);
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



long ARM_ARM_OPTIONALACCRUALZCBOND (double startDate,
									double endDate,
									double strike,
									long nbCurPerforAcc,
									long payFreqId,
									long ccyId,
									ARM_result& result,
									long objId)
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
				ARM_REQUEST_ID = RPC_SET_OPTIONALACCRUALZC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_OPTIONALACCRUALZC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
            REQ_SetDouble (reqIn, strike, ++idx);
            REQ_SetLong (reqIn, nbCurPerforAcc, ++idx);
            REQ_SetLong (reqIn, payFreqId, ++idx);
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


long ARM_ARM_FlexAccretSwaption (double startDate,
								 double endDate,
								 double fixedRate,
								 long nbCurPerforAcc,
								 long receiveOrPay,
								 long freqId,
								 long liborTypeId,
								 double spread,
								 long exerciseTypeId,
								 long ccyId,
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
				ARM_REQUEST_ID = RPC_SET_FLEXACCSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_FLEXACCSWAPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
            REQ_SetDouble (reqIn, fixedRate, ++idx);
            REQ_SetLong (reqIn, nbCurPerforAcc, ++idx);
            REQ_SetLong (reqIn, receiveOrPay, ++idx);
            REQ_SetLong (reqIn, freqId, ++idx);
            REQ_SetLong (reqIn, liborTypeId, ++idx);
            REQ_SetDouble (reqIn, spread, ++idx);
            REQ_SetLong (reqIn, exerciseTypeId, ++idx);
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


// EOF %M%