#include "ARM_interglob.h"
#include "ARM_option.h"



long ARM_OPTION (long underId,
			     double maturityDate,
				 double strike,
				 long optionType,
				 long exerciseType,
				 long strikeType,
				 double FstXDate,
				 ARM_result& result,
				 long objId)
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
				ARM_REQUEST_ID = RPC_SETOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, underId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(maturityDate), ++idx);
			REQ_SetDouble (reqIn, strike, ++idx);
			REQ_SetLong (reqIn, optionType, ++idx);
			REQ_SetLong (reqIn, exerciseType, ++idx);
			REQ_SetLong (reqIn, strikeType, ++idx);
			if(FstXDate == -1)
			{
				REQ_SetString (reqIn, (const char*)ARM_DEFAULT_DATE, ++idx);
			}
			else
			{
				REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(FstXDate), ++idx);
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

long ARM_EXOPTION (long underId,
				   long optionType,
				   long styleId,
				   long kRefValId,
				   ARM_result& result,
				   long objId)
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
				ARM_REQUEST_ID = RPC_SET_EXOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_EXOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, underId, idx);
			REQ_SetLong (reqIn, optionType, ++idx);
			REQ_SetLong (reqIn, styleId, ++idx);
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

long ARM_VolImp (long secId,
				 long modId,
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

			ARM_REQUEST_ID = RPC_VOLIMP;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, secId, idx);
			REQ_SetLong (reqIn, modId, ++idx);
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


long ARM_GetUnderPrice (long secId,
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

			ARM_REQUEST_ID = RPC_GETUNDPRICE;

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


long ARM_bsOption (double spot,
				   double strike,
				   double volatility,
				   double dividend,
				   double discountRate,
				   double maturity,
				   long CallPut,
				   ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BSOPTION;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetDouble (reqIn, spot/100., idx);
			REQ_SetDouble (reqIn, strike/100., ++idx);
			REQ_SetDouble (reqIn, volatility/100., ++idx);
			REQ_SetDouble (reqIn, dividend/100., ++idx);
			REQ_SetDouble (reqIn, discountRate/100., ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, CallPut, ++idx);

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

long ARM_bsDelta (double spot,
				  double strike,
				  double volatility,
				  double dividend,
				  double discountRate,
				  double maturity,
				  long CallPut,
				  ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BSDELTA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetDouble (reqIn, spot/100., idx);
			REQ_SetDouble (reqIn, strike/100., ++idx);
			REQ_SetDouble (reqIn, volatility/100., ++idx);
			REQ_SetDouble (reqIn, dividend/100., ++idx);
			REQ_SetDouble (reqIn, discountRate/100., ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, CallPut, ++idx);

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


long ARM_bsVega (double spot,
				 double strike,
				 double volatility,
				 double dividend,
				 double discountRate,
				 double maturity,
				 long CallPut,
				 ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BSVEGA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetDouble (reqIn, spot/100., idx);
			REQ_SetDouble (reqIn, strike/100., ++idx);
			REQ_SetDouble (reqIn, volatility/100., ++idx);
			REQ_SetDouble (reqIn, dividend/100., ++idx);
			REQ_SetDouble (reqIn, discountRate/100., ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, CallPut, ++idx);

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

long ARM_bsTheta (double spot,
				  double strike,
				  double volatility,
				  double dividend,
				  double discountRate,
				  double maturity,
				  long CallPut,
				  ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BSTHETA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetDouble (reqIn, spot/100., idx);
			REQ_SetDouble (reqIn, strike/100., ++idx);
			REQ_SetDouble (reqIn, volatility/100., ++idx);
			REQ_SetDouble (reqIn, dividend/100., ++idx);
			REQ_SetDouble (reqIn, discountRate/100., ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, CallPut, ++idx);

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

long ARM_bsGamma (double spot,
				  double strike,
				  double volatility,
				  double dividend,
				  double discountRate,
				  double maturity,
				  long CallPut,
				  ARM_result& result)
{

	CCString msg ("");
	double bsOpt=0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	long idx = 0;
		
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_BSGAMMA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);

			REQ_SetDouble (reqIn, spot/100., idx);
			REQ_SetDouble (reqIn, strike/100., ++idx);
			REQ_SetDouble (reqIn, volatility/100., ++idx);
			REQ_SetDouble (reqIn, dividend/100., ++idx);
			REQ_SetDouble (reqIn, discountRate/100., ++idx);
			REQ_SetDouble (reqIn, maturity, ++idx);
			REQ_SetLong (reqIn, CallPut, ++idx);

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


long ARM_ARM_OPTIONPORTFOLIO (long portfolioId,
							  long styleId,
							  long kRefValId,
							  long optionType,
							  ARM_result& result,
							  long objId)
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
				ARM_REQUEST_ID = RPC_SETPORTOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEPORTOPTION;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, portfolioId, idx);
			REQ_SetLong (reqIn, styleId, ++idx);
			REQ_SetLong (reqIn, kRefValId, ++idx);
			REQ_SetLong (reqIn, optionType, ++idx);
						
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