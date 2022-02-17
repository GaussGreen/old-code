#include "ARM_interglob.h"
#include "ARM_mod.h"



long ARM_bsmodel (double date, double spot, 
				  long dividend_type, double dividend,
				  long discrate_type, double discrate,
				  long volat_type, double volat,
				  long typstk, ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
				ARM_REQUEST_ID = RPC_SETBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(date), idx);
			REQ_SetDouble (reqIn, spot, ++idx);

			REQ_SetLong (reqIn, dividend_type, ++idx);
			REQ_SetDouble (reqIn, dividend, ++idx);

			REQ_SetLong (reqIn, discrate_type, ++idx);
			REQ_SetDouble (reqIn, discrate, ++idx);

			REQ_SetLong (reqIn, volat_type, ++idx);
			REQ_SetDouble (reqIn, volat, ++idx);

			REQ_SetLong (reqIn, typstk, ++idx);

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

long ARM_ycmod (long zeroCurveId,
				long discCurveId,
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
				ARM_REQUEST_ID = RPC_SETYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zeroCurveId, idx);
			REQ_SetLong (reqIn, discCurveId, ++idx);

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

long ARM_bsslmodel(double startDate, long zc_type, double zc, 
				   long volSpreadLock_type, double volSpreadLock,
				   long capVol_type, double capVol,
				   long indexVol_type, double indexVol,
				   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
				ARM_REQUEST_ID = RPC_SETBSSLMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBSSLMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
	
			REQ_SetLong(reqIn, zc_type, ++idx);
			REQ_SetDouble(reqIn, zc, ++idx);

			REQ_SetLong(reqIn, volSpreadLock_type, ++idx);
			REQ_SetDouble(reqIn, volSpreadLock, ++idx);

			REQ_SetLong(reqIn, capVol_type, ++idx);
			REQ_SetDouble(reqIn, capVol, ++idx);

			REQ_SetLong(reqIn, indexVol_type, ++idx);
            REQ_SetDouble(reqIn, indexVol, ++idx);


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

long ARM_GetParameter(long modId, long paramId, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_GETPARAMETER;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, modId, idx);
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


long ARM_HWTREE(long zcId, double begDate, double endDate,
                long nbSteps, double a, 
                double sigma, ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETHWTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHWTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(begDate), ++idx);
            REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong(reqIn, nbSteps, ++idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetDouble(reqIn, sigma, ++idx);

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

long ARM_HWSIGCONST(long zcId, double begDate, double endDate,
                    long nbSteps, double a, 
                    double sigma, ARM_result& result, long objId)
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

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SETHWSIGCST;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHWSIGCST;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

            REQ_SetLong(reqIn, zcId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(begDate), ++idx);
            REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong(reqIn, nbSteps, ++idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetDouble(reqIn, sigma, ++idx);


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

long ARM_GYCMODEL(long zcId,
                  double a, 
                  double sigma, 
                  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SETGYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEGYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

            REQ_SetLong(reqIn, zcId, idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetDouble(reqIn, sigma, ++idx);


			/*--- CORBA Server Call ---*/
			CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

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

long ARM_IDTREEHW(double startDate, double horizon, long nbSteps,
                  double dMeanRevSpeed, double fMeanRevSpeed,
                  double dSigma, double fSigma,
                  double prtyCorr, double prtyVol,
                  double ratesCorr, 
                  long   dZcId, long fZcId,
                  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 12;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SETIR3DTHWMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEIR3DTHWMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

		    REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(startDate), idx);
            REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
          
			REQ_SetLong(reqIn, nbSteps, ++idx);

            REQ_SetDouble(reqIn, dMeanRevSpeed, ++idx);
            REQ_SetDouble(reqIn, fMeanRevSpeed, ++idx);

            REQ_SetDouble(reqIn, dSigma, ++idx);
            REQ_SetDouble(reqIn, fSigma, ++idx);

            REQ_SetDouble(reqIn, prtyCorr, ++idx);
            REQ_SetDouble(reqIn, prtyVol, ++idx);
            REQ_SetDouble(reqIn, ratesCorr, ++idx);

			REQ_SetLong(reqIn, dZcId, ++idx);
            REQ_SetLong(reqIn, fZcId, ++idx);

            
            /*--- CORBA Server Call ---*/
			CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

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

long ARM_HWTWOFACTMOD(long zcId,
                      double a, 
                      double sigma1,
                      double b,
                      double sigma2,
                      double rho,
                      ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 6;
	CCString stringObjectId;
			
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			stringObjectId = GetLastCurCellEnvValue ();

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SETHW2FMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHW2FMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

            REQ_SetLong(reqIn, zcId, idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetDouble(reqIn, sigma1, ++idx);
            REQ_SetDouble(reqIn, b, ++idx);
            REQ_SetDouble(reqIn, sigma2, ++idx);
            REQ_SetDouble(reqIn, rho, ++idx);
			
            
            /*--- CORBA Server Call ---*/
			CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

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




long ARM_LDCANA(long zcId, 
                const VECTOR<CCString>& resetDates, 
                const VECTOR<double>& shifts,
                const VECTOR<double>& vols,
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
				ARM_REQUEST_ID = RPC_SETLOGDECANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATELOGDECANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

	        REQ_SetLong(reqIn, resetDates.size(), ++idx);
            REQ_SetStringVector(reqIn, resetDates, ++idx);
            
            REQ_SetLong(reqIn, shifts.size(), ++idx);
            REQ_SetDoubleVector(reqIn, shifts, ++idx);

            REQ_SetDoubleVector(reqIn, vols, ++idx);


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



long ARM_HWSIGVAR(long zcId, 
                  long begDate, long endDate,
                  long nbSteps, double a, long real_size,
                  const VECTOR<CCString>& dates, const VECTOR<double>& sigmas,
                  ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETHWSIGVAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHWSIGVAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(begDate), ++idx);
            REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong(reqIn, nbSteps, ++idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetLong(reqIn, real_size, ++idx);
            REQ_SetStringVector(reqIn, dates, ++idx);
            REQ_SetDoubleVector(reqIn, sigmas, ++idx);

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

long ARM_HWANALYTICSIGVAR(long zcId, 
                          double a, long real_size,
                          const VECTOR<CCString>& dates, 
                          const VECTOR<double>& sigmas,
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
				ARM_REQUEST_ID = RPC_SETHWSIGVARANALYTIC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHWSIGVARANALYTIC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

            REQ_SetLong(reqIn, zcId, idx);
			REQ_SetDouble(reqIn, a, ++idx);
            REQ_SetLong(reqIn, real_size, ++idx);
            REQ_SetStringVector(reqIn, dates, ++idx);
            REQ_SetDoubleVector(reqIn, sigmas, ++idx);

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

long ARM_HWTWOFANALYTICSIGVAR (long zcId, 
                               double a,
							   double b,
							   double sigmaRatio,
							   double rho,
							   long real_size,
							   const VECTOR<CCString>& dates, 
							   const VECTOR<double>& sigmas,
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
				ARM_REQUEST_ID = RPC_SETHW2FSIGVARANALYTIC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHW2FSIGVARANALYTIC;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

            REQ_SetLong(reqIn, zcId, idx);
			REQ_SetDouble(reqIn, a, ++idx);
			REQ_SetDouble(reqIn, b, ++idx);
			REQ_SetDouble(reqIn, rho, ++idx);
			REQ_SetDouble(reqIn, sigmaRatio, ++idx);
            REQ_SetLong(reqIn, real_size, ++idx);
            REQ_SetStringVector(reqIn, dates, ++idx);
            REQ_SetDoubleVector(reqIn, sigmas, ++idx);

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

long ARM_CRRTREE(double begDate, double endDate,
                 long nbSteps, double spot, 
                 double dividend, double discrate,
                 double volat,
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

			if ( objId != -1 )
			{
				ARM_REQUEST_ID = RPC_SETCRRTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATECRRTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(begDate), idx);
            REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong(reqIn, nbSteps, ++idx);
			REQ_SetDouble(reqIn, spot, ++idx);
            REQ_SetDouble(reqIn, dividend, ++idx);
            REQ_SetDouble(reqIn, discrate, ++idx);
            REQ_SetDouble(reqIn, volat, ++idx);


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

long ARM_GYCLSMODEL (long curveId, double a, double sigma,
					 double solv, double solvVol, double solvRCorr,
					 double lossRate,
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
				ARM_REQUEST_ID = RPC_SETGYCLSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEGYCLSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, curveId, idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetDouble (reqIn, sigma, ++idx);
			REQ_SetDouble (reqIn, solv, ++idx);
			REQ_SetDouble (reqIn, solvVol, ++idx);
			REQ_SetDouble (reqIn, solvRCorr, ++idx);
			REQ_SetDouble (reqIn, lossRate, ++idx);
			
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

long ARM_GTWOYC (double dMeanRevSpeed,
				 double dSigma,
				 long dZcId,
				 long fZcId,
				 double ratesCorr,
				 double fMeanRevSpeed,
				 double fSigma,
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
				ARM_REQUEST_ID = RPC_SETG2YCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEG2YCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, dMeanRevSpeed, idx);
			REQ_SetDouble (reqIn, dSigma, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			REQ_SetDouble (reqIn, fMeanRevSpeed, ++idx);
			REQ_SetDouble (reqIn, fSigma, ++idx);
			
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

long ARM_ARM_COVAR (long modId, double maturity1, double maturity2, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_COVARTXFWD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, modId, idx);
			REQ_SetDouble (reqIn, maturity1, ++idx);
			REQ_SetDouble (reqIn, maturity2, ++idx);
						
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

long ARM_ARM_CORREL (long modId, double maturity1, double maturity2, ARM_result& result)
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

			ARM_REQUEST_ID = RPC_CORRTXFWD;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, modId, idx);
			REQ_SetDouble (reqIn, maturity1, ++idx);
			REQ_SetDouble (reqIn, maturity2, ++idx);
						
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

long ARM_HWTWOFFNMONTECARLOSV (long zcId,
							   double horizon,
							   double a,
							   double b,
							   double sigmaRatio,
							   double rho,
							   const VECTOR<double>& sigmaDate,
							   const VECTOR<double>& sigmaVal,
							   long nbTraj,
							   ARM_result& result,
							   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 10;
	CCString stringObjectId;
			
	try
	{
		if(sigmaDate.size () != sigmaVal.size ())
		{
			result.setMsg ("ARM_ERR: dates and values array must have same size");
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
				ARM_REQUEST_ID = RPC_SETFNHW2FMONTECARLOSV;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFNHW2FMONTECARLOSV;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> sigmaDate_str;
			for(int i = 0; i < sigmaDate.size (); i++)
			{
				sigmaDate_str.push_back (XLDATE2ARMDATE (sigmaDate[i]));
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (horizon), ++idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetDouble (reqIn, b, ++idx);
			REQ_SetDouble (reqIn, sigmaRatio, ++idx);
			REQ_SetDouble (reqIn, rho, ++idx);
			REQ_SetLong (reqIn, nbTraj, ++idx);
			REQ_SetLong (reqIn, sigmaDate.size (), ++idx);
			REQ_SetStringVector (reqIn, sigmaDate_str, ++idx);
			REQ_SetDoubleVector (reqIn, sigmaVal, ++idx);

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

long ARM_HWTWOFTREE (long zcId,
					 double dateDeb,
					 double dateFin,
					 long pas,
					 double a,
					 double sig1,
					 double b,
					 double sig2,
					 double rho,
					 ARM_result& result,
					 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
				ARM_REQUEST_ID = RPC_SETHW2FTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEHW2FTREEMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(dateDeb), ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(dateFin), ++idx);
			REQ_SetLong (reqIn, pas, ++idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetDouble (reqIn, sig1, ++idx);
			REQ_SetDouble (reqIn, b, ++idx);
			REQ_SetDouble (reqIn, sig2, ++idx);
			REQ_SetDouble (reqIn, rho, ++idx);

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

long ARM_HWFNMONTECARLOSV (long zcId,
						   double horizon,
						   double a,
						   const VECTOR<double>& sigmaDate,
						   const VECTOR<double>& sigmaVal,
						   long nbTraj,
						   ARM_result& result,
						   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
	CCString stringObjectId;
			
	try
	{
		if(sigmaDate.size () != sigmaVal.size ())
		{
			result.setMsg ("ARM_ERR: dates and values array must have same size");
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
				ARM_REQUEST_ID = RPC_SETFNHWMONTECARLOSV;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFNHWMONTECARLOSV;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> sigmaDate_str;
			for(int i = 0; i < sigmaDate.size (); i++)
			{
				sigmaDate_str.push_back (XLDATE2ARMDATE (sigmaDate[i]));
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (horizon), ++idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetLong (reqIn, nbTraj, ++idx);
			REQ_SetLong (reqIn, sigmaDate.size (), ++idx);
			REQ_SetStringVector (reqIn, sigmaDate_str, ++idx);
			REQ_SetDoubleVector (reqIn, sigmaVal, ++idx);

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

long ARM_HWFNMONTECARLO (long zcId,
						 double horizon,
						 double a,
						 double sigma,
						 long nbTraj,
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
				ARM_REQUEST_ID = RPC_SETFNHWMONTECARLO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFNHWMONTECARLO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (horizon), ++idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetDouble (reqIn, sigma, ++idx);
			REQ_SetLong (reqIn, nbTraj, ++idx);
			
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

long ARM_DFBS (long dVolId,
			   long fVolId,
			   long dZcId,
			   long fZcId,
			   double fxCorr,
			   double fxVol,
			   double ratesCorr,
			   long dTypVol,
			   long fTypVol,
			   ARM_result& result,
			   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
				ARM_REQUEST_ID = RPC_SETDFBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, dVolId, idx);
			REQ_SetLong (reqIn, fVolId, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetDouble (reqIn, fxCorr, ++idx);
			REQ_SetDouble (reqIn, fxVol, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			REQ_SetLong (reqIn, dTypVol, ++idx);
			REQ_SetLong (reqIn, fTypVol, ++idx);

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

long ARM_DFFXBS (long dVolId,
			     long fVolId,
			     long dZcId,
			     long fZcId,
			     long dFxCorrId,
			     long fFxCorrId,
			     long fxVolId,
			     double ratesCorr,
				 long dBSZcId,
				 long fBSZcId,
				 double spot,
			     ARM_result& result,
			     long objId)
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
				ARM_REQUEST_ID = RPC_SETDFFXBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFFXBSMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, dVolId, idx);
			REQ_SetLong (reqIn, fVolId, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetLong (reqIn, dFxCorrId, ++idx);
			REQ_SetLong (reqIn, fFxCorrId, ++idx);
			REQ_SetLong (reqIn, fxVolId, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			REQ_SetLong (reqIn, dBSZcId, ++idx);
			REQ_SetLong (reqIn, fBSZcId, ++idx);
			REQ_SetDouble (reqIn, spot, ++idx);

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

long ARM_DFGYC (double dMeanRevSpeed,
				double fMeanRevSpeed,
				double dSigma,
				double fSigma,
				double fxCorr,
				double fxVol,
				double ratesCorr,
				long dZcId,
			    long fZcId,
			    ARM_result& result,
			    long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
				ARM_REQUEST_ID = RPC_SETDFGYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFGYCMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetDouble (reqIn, dMeanRevSpeed, idx);
			REQ_SetDouble (reqIn, fMeanRevSpeed, ++idx);
			REQ_SetDouble (reqIn, dSigma, ++idx);
			REQ_SetDouble (reqIn, fSigma, ++idx);
			REQ_SetDouble (reqIn, fxCorr, ++idx);
			REQ_SetDouble (reqIn, fxVol, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			
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

long ARM_DFFXBSPLAIN (long dVolId,
					  long fVolId,
					  long dZcId,
					  long fZcId,
					  double dFxCorr,
					  double fFxCorr,
					  double fxVol,
					  double ratesCorr,
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
				ARM_REQUEST_ID = RPC_SETDFFXBSPLAINMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFFXBSPLAINMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, dVolId, idx);
			REQ_SetLong (reqIn, fVolId, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetDouble (reqIn, dFxCorr, ++idx);
			REQ_SetDouble (reqIn, fFxCorr, ++idx);
			REQ_SetDouble (reqIn, fxVol, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			
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

long ARM_DFHWXSIGVAR (long dZcId,
					  long fZcId,
					  double dMeanRevSpeed,
					  double fMeanRevSpeed,
					  const VECTOR<double>& dDate,
					  const VECTOR<double>& fDate,
					  const VECTOR<double>& dSigma,
					  const VECTOR<double>& fSigma,
					  long dFxCorrId,
					  long fFxCorrId,
					  long fxVolId,
					  double ratesCorr,
					  ARM_result& result,
					  long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 14;
	CCString stringObjectId;
			
	try
	{
		if(dDate.size () != dSigma.size ())
		{
			result.setMsg ("ARM_ERR: domestic date and sigma array must have same size");
			return ARM_KO;
		}

		if(fDate.size () != fSigma.size ())
		{
			result.setMsg ("ARM_ERR: foreign date and sigma array must have same size");
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
				ARM_REQUEST_ID = RPC_SETDFHWXSIGVAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFHWXSIGVAR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> dDate_str;
			for(int i = 0; i < dDate.size (); i++)
			{
				dDate_str.push_back (XLDATE2ARMDATE (dDate[i]));
			}
			VECTOR<CCString> fDate_str;
			for(i = 0; i < fDate.size (); i++)
			{
				fDate_str.push_back (XLDATE2ARMDATE (fDate[i]));
			}

			REQ_SetLong (reqIn, dZcId, idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetDouble (reqIn, dMeanRevSpeed, ++idx);
			REQ_SetDouble (reqIn, fMeanRevSpeed, ++idx);
			REQ_SetLong (reqIn, dDate.size (), ++idx);
			REQ_SetLong (reqIn, fDate.size (), ++idx);
			REQ_SetStringVector (reqIn, dDate_str, ++idx);
			REQ_SetStringVector (reqIn, fDate_str, ++idx);
			REQ_SetDoubleVector (reqIn, dSigma, ++idx);
			REQ_SetDoubleVector (reqIn, fSigma, ++idx);
			REQ_SetLong (reqIn, dFxCorrId, ++idx);
			REQ_SetLong (reqIn, fFxCorrId, ++idx);
			REQ_SetLong (reqIn, fxVolId, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			
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
					  
long ARM_DFHWSIGVARTREE (double startDate,
						 double horizon,
						 long numSteps,
						 long dZcId,
						 long fZcId,
						 double dMeanRevSpeed,
						 const VECTOR<double>& dDate,
						 const VECTOR<double>& dSigma,
						 double fMeanRevSpeed,
						 const VECTOR<double>& fDate,
						 const VECTOR<double>& fSigma,
						 long dFxCorrId,
						 long fFxCorrId,
						 long fxVolId,
						 double ratesCorr,
						 double fxSpotRate,
						 ARM_result& result,
					     long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
	CCString stringObjectId;
			
	try
	{
		if(dDate.size () != dSigma.size ())
		{
			result.setMsg ("ARM_ERR: domestic date and sigma array must have same size");
			return ARM_KO;
		}

		if(fDate.size () != fSigma.size ())
		{
			result.setMsg ("ARM_ERR: foreign date and sigma array must have same size");
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
				ARM_REQUEST_ID = RPC_SETDFHWSIGVARTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEDFHWSIGVARTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			VECTOR<CCString> dDate_str;
			for(int i = 0; i < dDate.size (); i++)
			{
				dDate_str.push_back (XLDATE2ARMDATE (dDate[i]));
			}
			VECTOR<CCString> fDate_str;
			for(i = 0; i < fDate.size (); i++)
			{
				fDate_str.push_back (XLDATE2ARMDATE (fDate[i]));
			}

			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (startDate), idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (horizon), ++idx);
			REQ_SetLong (reqIn, numSteps, ++idx);
			REQ_SetLong (reqIn, dZcId, ++idx);
			REQ_SetLong (reqIn, fZcId, ++idx);
			REQ_SetDouble (reqIn, dMeanRevSpeed, ++idx);
			REQ_SetDouble (reqIn, fMeanRevSpeed, ++idx);
			REQ_SetLong (reqIn, dDate.size (), ++idx);
			REQ_SetLong (reqIn, fDate.size (), ++idx);
			REQ_SetStringVector (reqIn, dDate_str, ++idx);
			REQ_SetStringVector (reqIn, fDate_str, ++idx);
			REQ_SetDoubleVector (reqIn, dSigma, ++idx);
			REQ_SetDoubleVector (reqIn, fSigma, ++idx);
			REQ_SetLong (reqIn, dFxCorrId, ++idx);
			REQ_SetLong (reqIn, fFxCorrId, ++idx);
			REQ_SetLong (reqIn, fxVolId, ++idx);
			REQ_SetDouble (reqIn, ratesCorr, ++idx);
			REQ_SetDouble (reqIn, fxSpotRate, ++idx);
			
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

long ARM_BKIRTREE (long zcId,
				   double startDate,
				   double endDate,
				   long pas,
				   double a,
				   double sigma,
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
				ARM_REQUEST_ID = RPC_SETBKIRTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBKIRTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (startDate), ++idx);
			REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE (endDate), ++idx);
			REQ_SetLong (reqIn, pas, ++idx);
			REQ_SetDouble (reqIn, a, ++idx);
			REQ_SetDouble (reqIn, sigma, ++idx);
						
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



long ARM_LDCMC_FROM_ANA(long anaModId, double horizon,
                        long nbTraj, long mcmethod, long pricerTypeId,
                        ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 5;
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
			   ARM_REQUEST_ID = RPC_SETLOGDECMC_FROM_ANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATELOGDECMC_FROM_ANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, anaModId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
			REQ_SetLong(reqIn, nbTraj, ++idx);
			REQ_SetLong(reqIn, mcmethod, ++idx);
			REQ_SetLong(reqIn, pricerTypeId, ++idx);

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





long ARM_FRMTREE(long FrmAnaId, 
                 double horizon, 
                 long fineMonth,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETFRMTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, FrmAnaId, idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonth, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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






long ARM_FRMTREE_AUTO(long zcId, 
                 long volId, 
                 long smileId,
                 long autoMode,
                 double horizon,
                 long fineMonth,
                 double shapeDecay,
                 double shapeSlope,
                 double shapeAsymptote,
                 long nbFactor,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETFRMTREE_AUTO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMTREE_AUTO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

            REQ_SetLong(reqIn, volId, ++idx);

            REQ_SetLong(reqIn, smileId, ++idx);

            REQ_SetLong(reqIn, autoMode, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonth, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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


long ARM_FRMTREE_AUTO_G(long zcId,
						long zdId,
						long volId,
						long smileId,
						long irgvolId,
						long irgsmileId,
						long autoModeId,
						double horizon,
						long fineMonthId,
						double shapeDecay,
						double shapeSlope,
						double shapeAsymptote,
						long nbFactor,
						const VECTOR<double>& corrMatu,
						const VECTOR<double>& corrMatrix,
						const VECTOR<double>& corr2Matu,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
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
				ARM_REQUEST_ID = RPC_SETFRMTREE_AUTO_G;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMTREE_AUTO_G;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);
			REQ_SetLong(reqIn, zdId, ++idx);

            REQ_SetLong(reqIn, volId, ++idx);
            REQ_SetLong(reqIn, smileId, ++idx);

            REQ_SetLong(reqIn, irgvolId, ++idx);
            REQ_SetLong(reqIn, irgsmileId, ++idx);

            REQ_SetLong(reqIn, autoModeId, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonthId, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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


long ARM_FRMTREE_AUTO_B(long zcId,
						long zdId,
						long volId,
						long smileId,
						long autoModeId,
						double horizon,
						long fineMonthId,
						double shapeDecay,
						double shapeSlope,
						double shapeAsymptote,
						long nbFactor,
						const VECTOR<double>& corrMatu,
						const VECTOR<double>& corrMatrix,
						const VECTOR<double>& corr2Matu,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 16;
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
				ARM_REQUEST_ID = RPC_SETFRMTREE_AUTO_B;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMTREE_AUTO_B;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);
			REQ_SetLong(reqIn, zdId, ++idx);

            REQ_SetLong(reqIn, volId, ++idx);
            REQ_SetLong(reqIn, smileId, ++idx);

            REQ_SetLong(reqIn, autoModeId, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonthId, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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



long ARM_FRMANA(long zcId, 
                 const VECTOR<CCString>& resetDates, 
                 const VECTOR<double>&   spotVols, 
                 long shapeType,
                 double shapeDecay,
                 double shapeSlope,
                 double shapeAsymptote,
                 long nbFactor,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 14;
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
				ARM_REQUEST_ID = RPC_SETFRMANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

            REQ_SetLong(reqIn, resetDates.size(), ++idx);

            REQ_SetStringVector(reqIn, resetDates, ++idx);

            REQ_SetDoubleVector(reqIn, spotVols, ++idx);

            REQ_SetLong(reqIn, shapeType, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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





long ARM_FRMANA_PORT(long zcId, 
                 long portId, 
                 const VECTOR<CCString>& resetDates, 
                 double precision,
                 double min_paras,
                 double max_paras,
                 long   max_iters,
                 long shapeType,
                 double shapeDecay,
                 double shapeSlope,
                 double shapeAsymptote,
                 long nbFactor,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
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
				ARM_REQUEST_ID = RPC_SETFRMANA_PORT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMANA_PORT;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

            REQ_SetLong(reqIn, portId, ++idx);

            REQ_SetLong(reqIn, resetDates.size(), ++idx);

            REQ_SetStringVector(reqIn, resetDates, ++idx);

            REQ_SetDouble (reqIn, precision, ++idx);

            REQ_SetDouble (reqIn, min_paras, ++idx);
            
            REQ_SetDouble (reqIn, max_paras, ++idx);

            REQ_SetLong(reqIn, max_iters, ++idx);

            REQ_SetLong(reqIn, shapeType, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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


long ARM_FRMLSMC(long FrmAnaId, 
                 double horizon, 
                 long nbTraj, 
                 long fineMonth,
                 long mcMethod,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
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
				ARM_REQUEST_ID = RPC_SETFRMLSMONTECARLO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMLSMONTECARLO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, FrmAnaId, idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, fineMonth, ++idx);

            REQ_SetLong(reqIn, mcMethod, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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
                 





long ARM_BASISMCFRM(long   zcId, 
					long   baZcId,
                    long   volId, 
                    long   smileId,
                    long   productType,
                    double horizon,
                    long   nbTraj,
                    long   MCGeneratorType,
                    double shapeDecay,
                    double shapeSlope,
                    double shapeAsymptote,
                    long   nbFactor,
                    const  VECTOR<double>& indexes,
                    const  VECTOR<double>& correlatedIndex,
                    const  VECTOR<double>& corrMatrix,
					long   control,
					long   seed,
                    ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 19;
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
				ARM_REQUEST_ID = RPC_BMCFRM_SET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_BMCFRM_CREATE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

			REQ_SetLong(reqIn, baZcId, ++idx);

            REQ_SetLong(reqIn, volId, ++idx);

            REQ_SetLong(reqIn, smileId, ++idx);

            REQ_SetLong(reqIn, productType, ++idx);

			REQ_SetString(reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, MCGeneratorType, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, indexes.size(), ++idx);

            REQ_SetDoubleVector(reqIn, indexes, ++idx);

			REQ_SetLong(reqIn, correlatedIndex.size(), ++idx);

			REQ_SetDoubleVector(reqIn, correlatedIndex, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetLong(reqIn, control, ++idx);

			REQ_SetLong(reqIn, seed, ++idx);


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



long ARM_BASISMCFRM2CR(long   zcId,
					   long   baZcId,
					   long   swoptvolId,
					   long   swoptsmileId,
					   long   irgvolId,
					   long   irgsmileId,
					   long   productType,
					   double horizon,
					   long   nbTraj,
					   long   MCGeneratorType,
					   double shapeDecay,
					   double shapeSlope,
					   double shapeAsymptote,
					   long   nbFactor,
					   const  VECTOR<double>& indexes,
					   const  VECTOR<double>& correlatedIndex,
					   const  VECTOR<double>& corrMatrix,
					   long   control,
					   long   seed,
					   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 21;
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
				ARM_REQUEST_ID = RPC_BMCFRM2CR_SET;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_BMCFRM2CR_CREATE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

			REQ_SetLong(reqIn, baZcId, ++idx);

            REQ_SetLong(reqIn, swoptvolId, ++idx);

            REQ_SetLong(reqIn, swoptsmileId, ++idx);

            REQ_SetLong(reqIn, irgvolId, ++idx);

            REQ_SetLong(reqIn, irgsmileId, ++idx);

            REQ_SetLong(reqIn, productType, ++idx);

			REQ_SetString(reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, MCGeneratorType, ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, indexes.size(), ++idx);

            REQ_SetDoubleVector(reqIn, indexes, ++idx);

			REQ_SetLong(reqIn, correlatedIndex.size(), ++idx);

			REQ_SetDoubleVector(reqIn, correlatedIndex, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetLong(reqIn, control, ++idx);

			REQ_SetLong(reqIn, seed, ++idx);


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


long ARM_FRMLSMC_AUTO(long zcId, 
                 long volId, 
                 long smileId,
                 long autoModeId,
                 double horizon,
                 long fineMonthId,
                 long   nbTraj,
                 long mcMethod,
                 long noControl,
                 double shapeDecay,
                 double shapeSlope,
                 double shapeAsymptote,
                 long nbFactor,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
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
				ARM_REQUEST_ID = RPC_SETFRMLSMC_AUTO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMLSMC_AUTO;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

            REQ_SetLong(reqIn, volId, ++idx);

            REQ_SetLong(reqIn, smileId, ++idx);

            REQ_SetLong(reqIn, autoModeId, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonthId, ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, mcMethod, ++idx);

            REQ_SetLong(reqIn, noControl , ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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



long ARM_FRMLSMC_AUTO_2(long zcId, 
                 long swoptVolId, 
                 long swoptSmileId,
                 long irgVolId, 
                 long irgSmileId,
                 long autoMode,
                 double horizon,
                 long fineMonth,
                 long   nbTraj,
                 long mcMethod,
                 long noControl,
                 double shapeDecay,
                 double shapeSlope,
                 double shapeAsymptote,
                 long nbFactor,
                 const VECTOR<double>& corrMatu,
                 const VECTOR<double>& corrMatrix,
                 const VECTOR<double>& corr2Matu,
                 ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 20;
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
				ARM_REQUEST_ID = RPC_SETFRMLSMC_AUTO_2;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMLSMC_AUTO_2;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);

            REQ_SetLong(reqIn, swoptVolId, ++idx);

            REQ_SetLong(reqIn, swoptSmileId, ++idx);

            REQ_SetLong(reqIn, irgVolId, ++idx);

            REQ_SetLong(reqIn, irgSmileId, ++idx);

            REQ_SetLong(reqIn, autoMode, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonth, ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, mcMethod, ++idx);

            REQ_SetLong(reqIn, noControl , ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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


long ARM_FRMLSMC_AUTO_G(long zcId,
						long zdId,
						long swoptVolId,
						long swoptSmileId,
						long irgVolId,
						long irgSmileId,
						long intAutoMode,
						double horizon,
						long fineMonth,
						long nbTraj,
						long mcMethod,
						long noControl,
						double shapeDecay,
						double shapeSlope,
						double shapeAsymptote,
						long nbFactor,
						const VECTOR<double>& corrMatu,
						const VECTOR<double>& corrMatrix,
						const VECTOR<double>& corr2Matu,
						ARM_result& result,
						long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 21;
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
				ARM_REQUEST_ID = RPC_SETFRMLSMC_AUTO_G;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEFRMLSMC_AUTO_G;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong(reqIn, zcId, idx);
			REQ_SetLong(reqIn, zdId, ++idx);

            REQ_SetLong(reqIn, swoptVolId, ++idx);

            REQ_SetLong(reqIn, swoptSmileId, ++idx);

            REQ_SetLong(reqIn, irgVolId, ++idx);

            REQ_SetLong(reqIn, irgSmileId, ++idx);

            REQ_SetLong(reqIn, intAutoMode, ++idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

            REQ_SetLong(reqIn, fineMonth, ++idx);

            REQ_SetLong(reqIn, nbTraj, ++idx);

            REQ_SetLong(reqIn, mcMethod, ++idx);

            REQ_SetLong(reqIn, noControl , ++idx);

            REQ_SetDouble (reqIn, shapeDecay, ++idx);

            REQ_SetDouble (reqIn, shapeSlope, ++idx);

            REQ_SetDouble (reqIn, shapeAsymptote, ++idx);

            REQ_SetLong(reqIn, nbFactor, ++idx);

            REQ_SetLong(reqIn, corrMatu.size(), ++idx);

            REQ_SetLong(reqIn, corr2Matu.size(), ++idx);

            REQ_SetDoubleVector(reqIn, corrMatu, ++idx);

            REQ_SetDoubleVector(reqIn, corrMatrix, ++idx);

            REQ_SetDoubleVector(reqIn, corr2Matu, ++idx);


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



long ARM_SMILEDLDCMC(long anaModId, double horizon,
                               long nbTraj, double dVolSpot,
							   double Asymp_Row, double Asymp_Col,
							   double pUp, double pDown, long mcmethod,
                               ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
			   ARM_REQUEST_ID = RPC_SETSMILEDMCRNLDC;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATESMILEDMCRNLDC;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, anaModId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
			REQ_SetLong(reqIn, nbTraj, ++idx);
            REQ_SetDouble (reqIn, dVolSpot, ++idx);
            REQ_SetDouble (reqIn, Asymp_Row, ++idx);
            REQ_SetDouble (reqIn, Asymp_Col, ++idx);
			REQ_SetDouble (reqIn, pUp, ++idx);
			REQ_SetDouble (reqIn, pDown, ++idx);
			REQ_SetLong(reqIn, mcmethod, ++idx);
			

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


long ARM_FITTEDSMILEDLDCMC(long anaModId, double horizon,
                               long nbTraj, long volBSId,
							   long volSmileId, double strike1,
							   double strike2, double Asymp_Row, long mcmethod,
                               ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
			   ARM_REQUEST_ID = RPC_SETFITTEDSMILEDMCRNLDC;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATEFITTEDSMILEDMCRNLDC;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, anaModId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
			REQ_SetLong(reqIn, nbTraj, ++idx);
            REQ_SetLong (reqIn, volBSId, ++idx);
            REQ_SetLong (reqIn, volSmileId, ++idx);
            REQ_SetDouble (reqIn, strike1, ++idx);
			REQ_SetDouble (reqIn, strike2, ++idx);
			REQ_SetDouble (reqIn, Asymp_Row, ++idx);
			REQ_SetLong(reqIn, mcmethod, ++idx);
			

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


long ARM_SMILEDLDCANA(long anaModId, double dVolSpotUp,
						double Asymp_RowUp,
						double Asymp_ColUp,
						double pUp,
						double dVolSpotDo,
						double Asymp_RowDo,
						double Asymp_ColDo,
						double pDo,
						ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 9;
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
			   ARM_REQUEST_ID = RPC_SETSMILEDLDCANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATESMILEDLDCANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, anaModId, idx);
            REQ_SetDouble (reqIn, dVolSpotUp, ++idx);
			REQ_SetDouble (reqIn, Asymp_RowUp, ++idx);
			REQ_SetDouble (reqIn, Asymp_ColUp, ++idx);
			REQ_SetDouble (reqIn, pUp, ++idx);
            REQ_SetDouble (reqIn, dVolSpotDo, ++idx);
			REQ_SetDouble (reqIn, Asymp_RowDo, ++idx);
			REQ_SetDouble (reqIn, Asymp_ColDo, ++idx);
			REQ_SetDouble (reqIn, pDo, ++idx);
			

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


long ARM_SMILEDLDCFROMANA(long anaModId, double horizon,
							long nbTraj, long mcmethod,
							ARM_result& result, long objId)
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
			   ARM_REQUEST_ID = RPC_SETSMILEDLDCFROMANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATESMILEDLDCFROMANA;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, anaModId, idx);
			REQ_SetString(reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
			REQ_SetLong(reqIn, nbTraj, ++idx);
			REQ_SetLong(reqIn, mcmethod, ++idx);
			

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


long ARM_GLOBDFBS(long DomBSId, long DomCurrId,
				  long FrgBSId, long FrgCurrId,
				  long fxVolCrvId, long FFxCorrId,
				  long RatesCorrId,
				  ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 7;
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
			   ARM_REQUEST_ID = RPC_SETGLOBDFBSMODEL;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
			   REQ_SetLong (reqIn, objId, 0);
			
               idx = 1;
			}
			else
			{
			   ARM_REQUEST_ID = RPC_CREATEGLOBDFBSMODEL;

			   reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
               idx = 0;
			}

			REQ_SetLong(reqIn, DomBSId, idx);
			REQ_SetLong(reqIn, DomCurrId, ++idx);
			REQ_SetLong(reqIn, FrgBSId, ++idx);
			REQ_SetLong(reqIn, FrgCurrId, ++idx);
			REQ_SetLong(reqIn, fxVolCrvId, ++idx);
			REQ_SetLong(reqIn, FFxCorrId, ++idx);
			REQ_SetLong(reqIn, RatesCorrId, ++idx);
			

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



long ARM_bssmiledmodel(double date, double spot,
					   long dividend_type, double dividend,
					   long discrate_type, double discrate,
					   long volat_type, double volat, long typstk,
					   const VECTOR<double>& matu,
					   long betaType,
					   long betaObj,
					   const VECTOR<double>& beta,
					   long omegaType,
					   long omegaObj,
					   const VECTOR<double>& omega,
					   long isSABR,
					   ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 18;
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
				ARM_REQUEST_ID = RPC_SETBSSMILEDMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATEBSSMILEDMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(date), idx);
			REQ_SetDouble (reqIn, spot, ++idx);

			REQ_SetLong (reqIn, dividend_type, ++idx);
			REQ_SetDouble (reqIn, dividend, ++idx);

			REQ_SetLong (reqIn, discrate_type, ++idx);
			REQ_SetDouble (reqIn, discrate, ++idx);

			REQ_SetLong (reqIn, volat_type, ++idx);
			REQ_SetDouble (reqIn, volat, ++idx);

			REQ_SetLong (reqIn, typstk, ++idx);

			REQ_SetLong (reqIn, matu.size(), ++idx);
            REQ_SetDoubleVector(reqIn, matu, ++idx);

			REQ_SetLong (reqIn, betaType, ++idx);
            REQ_SetDoubleVector(reqIn, beta, ++idx);
			REQ_SetLong (reqIn, betaObj, ++idx);

			REQ_SetLong (reqIn, omegaType, ++idx);
            REQ_SetDoubleVector(reqIn, omega, ++idx);
			REQ_SetLong (reqIn, omegaObj, ++idx);

			REQ_SetLong (reqIn, isSABR, ++idx);

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


long ARM_BSCORRMODEL(double startDate,
					 long zcId,
					 long spreadlockId,
					 long capirgvolId,
					 long capcashvolId,
					 long indexadjvolId,
					 long spreadvolId,
					 long correlationsId,
					 long modelTypeId,
					 long volTypeId,
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
				ARM_REQUEST_ID = RPC_SET_BSCORRMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_BSCORRMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), idx);
			REQ_SetLong (reqIn, zcId, ++idx);

			REQ_SetLong (reqIn, spreadlockId, ++idx);
			REQ_SetLong (reqIn, capirgvolId, ++idx);
			REQ_SetLong (reqIn, capcashvolId, ++idx);
			REQ_SetLong (reqIn, indexadjvolId, ++idx);
			REQ_SetLong (reqIn, spreadvolId, ++idx);
			REQ_SetLong (reqIn, correlationsId, ++idx);
			REQ_SetLong (reqIn, modelTypeId, ++idx);
			REQ_SetLong (reqIn, volTypeId, ++idx);
			
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


long ARM_ARM_GetComputedSigmaSABR(long modId,
								  ARM_result& result)
{
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 1;
	long idx = 0;

    
	try
	{
		ARM_CORBA_init ();

		if (CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

			ARM_REQUEST_ID = RPC_GET_SABR_SIGMA;

			reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
			
			REQ_SetLong (reqIn, modId, idx);
            
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


long ARM_ARM_HWSIGVARFROMANA(long modId,
							 long endDate,
							 long nbSteps,
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
				ARM_REQUEST_ID = RPC_SET_HWSIGVARFROMANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_HWSIGVARFROMANA;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, modId, idx);

			REQ_SetString (reqIn, XLDATE2ARMDATE(endDate), ++idx);
			REQ_SetLong (reqIn, nbSteps, ++idx);
			
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

long ARM_ARM_CALIBRATIONHWSV(double asof,
							 long zcId,
							 long volId,
							 long secId,
							 double amin,
							 double amax,
							 double volmin,
							 double volmax,
							 const VECTOR<CCString>& dates,
							 long pfId,
							 ARM_result& result,
							 long objId)
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
				ARM_REQUEST_ID = RPC_SET_HWSIGVAR_CALIBRATOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_HWSIGVAR_CALIBRATOR;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(asof), idx);
			REQ_SetLong (reqIn, zcId, ++idx);
			REQ_SetLong (reqIn, volId, ++idx);
			REQ_SetLong (reqIn, secId, ++idx);
			REQ_SetDouble (reqIn, amin, ++idx);
			REQ_SetDouble (reqIn, amax, ++idx);
			REQ_SetDouble (reqIn, volmin, ++idx);
			REQ_SetDouble (reqIn, volmax, ++idx);

            REQ_SetLong(reqIn, dates.size(), ++idx);
            REQ_SetStringVector(reqIn, dates, ++idx);

			REQ_SetLong (reqIn, pfId, ++idx);
			
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


long ARM_ARM_CALIBRATE (long calibId,
						long calibType,
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
				ARM_REQUEST_ID = RPC_SET_HWSIGVARCALIB;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_HWSIGVARCALIB;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, calibId, idx);
			REQ_SetLong (reqIn, calibType, ++idx);
			
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


long ARM_ARM_BSSMILEDCALIBRATE(long modelId,
							   long pfId,
							   long calVolOrNotId,
							   long calRhoOrNotId,
							   long calNuOrNotId,
							   double volTenor,
							   double timeStep,
							   double minSig,
							   double maxSig,
							   double minRho,
							   double maxRho,
							   double minNu,
							   double maxNu,
							   long interpMethodId,
							   double tol,
							   long maxIter,
							   long gradCalcId,
							   double lambda,
							   long globOrBootstrapId,
							   double beginSmoothMatu,
							   ARM_result& result,
							   long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 20;
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
				ARM_REQUEST_ID = RPC_SET_CALIBRATE_BSSMILED;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CALIBRATE_BSSMILED;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, modelId, idx);
			REQ_SetLong (reqIn, pfId, ++idx);
			
			REQ_SetLong (reqIn, calVolOrNotId, ++idx);
			REQ_SetLong (reqIn, calRhoOrNotId, ++idx);
			REQ_SetLong (reqIn, calNuOrNotId, ++idx);

			REQ_SetDouble (reqIn, volTenor, ++idx);
			REQ_SetDouble (reqIn, timeStep, ++idx);

			REQ_SetDouble (reqIn, minSig, ++idx);
			REQ_SetDouble (reqIn, maxSig, ++idx);
			REQ_SetDouble (reqIn, minRho, ++idx);
			REQ_SetDouble (reqIn, maxRho, ++idx);
			REQ_SetDouble (reqIn, minNu, ++idx);
			REQ_SetDouble (reqIn, maxNu, ++idx);

			REQ_SetLong (reqIn, interpMethodId, ++idx);
			REQ_SetDouble (reqIn, tol, ++idx);
			REQ_SetLong (reqIn, maxIter, ++idx);
			REQ_SetLong (reqIn, gradCalcId, ++idx);

			REQ_SetDouble (reqIn, lambda, ++idx);
			REQ_SetLong (reqIn, globOrBootstrapId, ++idx);
			REQ_SetDouble (reqIn, beginSmoothMatu, ++idx);

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


long ARM_ARM_GETCALIBRATED_SIGRHONU(long modelId,
									long paramId,
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
				ARM_REQUEST_ID = RPC_SET_GET_SABR_SIG_RHO_NU;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_GET_SABR_SIG_RHO_NU;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, modelId, idx);
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

long ARM_ARM_FRMMARKOVTREE(double startDate,
						   double horizon,
						   long zcId,
						   long PathNumber,
						   const VECTOR<double>& Params,
						   long FRMModelId,
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
				ARM_REQUEST_ID = RPC_SET_FRMMARKOVTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_FRMMARKOVTREE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetString (reqIn, XLDATE2ARMDATE(startDate), idx);
			REQ_SetString (reqIn, XLDATE2ARMDATE(horizon), ++idx);

			REQ_SetLong (reqIn, zcId, ++idx);
			REQ_SetLong (reqIn, PathNumber, ++idx);

			REQ_SetLong (reqIn, Params.size(), ++idx);
            REQ_SetDoubleVector(reqIn, Params, ++idx);

			REQ_SetLong (reqIn, FRMModelId, ++idx);

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


long ARM_ARM_CalibrationFRMModel(long zcId,
								 long pf1Id,
								 long pf2Id,
								 long nbFactors,
								 long liborTypeId,
								 const VECTOR<double>& powers,
								 const VECTOR<double>& smoothParams,
								 const VECTOR<double>& ACalibrationSchedule,
								 const VECTOR<double>& KCalibrationSchedule,
								 long initA,
								 long initK,
								 const VECTOR<double>& meanRevA,
								 const VECTOR<double>& meanRevK,
								 const VECTOR<double>& mdec,
								 const VECTOR<double>& Bounds,
								 const VECTOR<double>& optimizerParams,
								 long calibrate,
								 ARM_result& result,
								 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 26;
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
				ARM_REQUEST_ID = RPC_SET_CALIBRATION_FRMMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_CALIBRATION_FRMMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetLong (reqIn, pf1Id, ++idx);
			REQ_SetLong (reqIn, pf2Id, ++idx);
			REQ_SetLong (reqIn, nbFactors, ++idx);
			REQ_SetLong (reqIn, liborTypeId, ++idx);

			REQ_SetLong (reqIn, powers.size(), ++idx);
            REQ_SetDoubleVector(reqIn, powers, ++idx);

			REQ_SetLong (reqIn, smoothParams.size(), ++idx);
            REQ_SetDoubleVector(reqIn, smoothParams, ++idx);

			int i;

			VECTOR<CCString> ACalibrationSchedule_str;
			for(i = 0; i < ACalibrationSchedule.size (); i++)
			{
				ACalibrationSchedule_str.push_back (XLDATE2ARMDATE (ACalibrationSchedule[i]));
			}

			VECTOR<CCString> KCalibrationSchedule_str;
			for(i = 0; i < KCalibrationSchedule.size (); i++)
			{
				KCalibrationSchedule_str.push_back (XLDATE2ARMDATE (KCalibrationSchedule[i]));
			}

			REQ_SetLong (reqIn, ACalibrationSchedule_str.size (), ++idx);
			REQ_SetStringVector (reqIn, ACalibrationSchedule_str, ++idx);

			REQ_SetLong (reqIn, KCalibrationSchedule_str.size (), ++idx);
			REQ_SetStringVector (reqIn, KCalibrationSchedule_str, ++idx);

			REQ_SetLong (reqIn, initA, ++idx);
			REQ_SetLong (reqIn, initK, ++idx);

			REQ_SetLong (reqIn, meanRevA.size(), ++idx);
            REQ_SetDoubleVector(reqIn, meanRevA, ++idx);

			REQ_SetLong (reqIn, meanRevK.size(), ++idx);
            REQ_SetDoubleVector(reqIn, meanRevK, ++idx);

			REQ_SetLong (reqIn, mdec.size(), ++idx);
            REQ_SetDoubleVector(reqIn, mdec, ++idx);

			REQ_SetLong (reqIn, Bounds.size(), ++idx);
            REQ_SetDoubleVector(reqIn, Bounds, ++idx);

			REQ_SetLong (reqIn, optimizerParams.size(), ++idx);
            REQ_SetDoubleVector(reqIn, optimizerParams, ++idx);

			REQ_SetLong (reqIn, calibrate, ++idx);

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



long ARM_ARM_BootstCalibFRMModel(long zcId,
								 long pf1Id,
								 long pf2Id,
								 long pf3Id,
								 long liborTypeId,
								 const VECTOR<double>& CalibParams,
								 const VECTOR<double>& initCurve,
								 const VECTOR<double>& mdec,
								 double meanRev,
								 long   nbfactor,
								 long   nbrows,
								 long   nbcolumns,
								 const VECTOR<double>& correlmatrix,
								 long  VolType,
								 long calibrate,  
								 ARM_result& result,
								 long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 19;
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
				ARM_REQUEST_ID = RPC_SET_BOOTSCALIB_FRMMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATE_BOOTSCALIB_FRMMODEL;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, zcId, idx);
			REQ_SetLong (reqIn, pf1Id, ++idx);
			REQ_SetLong (reqIn, pf2Id, ++idx);
			REQ_SetLong (reqIn, pf3Id, ++idx);
			REQ_SetLong (reqIn, liborTypeId, ++idx);

			REQ_SetLong (reqIn, CalibParams.size(), ++idx);
            REQ_SetDoubleVector(reqIn, CalibParams, ++idx);

			REQ_SetLong (reqIn, initCurve.size(), ++idx);
            REQ_SetDoubleVector(reqIn, initCurve, ++idx);

			REQ_SetLong (reqIn, mdec.size(), ++idx);
            REQ_SetDoubleVector(reqIn, mdec, ++idx);

			REQ_SetDouble(reqIn, meanRev, ++idx);
			REQ_SetLong (reqIn, nbfactor, ++idx);
			REQ_SetLong (reqIn, nbrows, ++idx);
			REQ_SetLong (reqIn, nbcolumns, ++idx);

			REQ_SetLong (reqIn, correlmatrix.size(), ++idx);
            REQ_SetDoubleVector(reqIn, correlmatrix, ++idx);

			REQ_SetLong (reqIn, VolType, ++idx);
			REQ_SetLong (reqIn, calibrate, ++idx);

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