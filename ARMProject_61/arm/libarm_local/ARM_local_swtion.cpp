#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_glob.h"

#include <libCCdate\CCdate.h>
#include <libCCtools++\CCstring.h>

#include <inst\swap.h>
#include "inst\swaption.h"
#include <inst\zcoptswapleg.h>
#include <inst\swaption_cf.h>
#include <inst\flexaccretswaption.h>
#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include "ARM_local_persistent.h"
#include "ARM_local_class.h"

#include <mod\qmodel.h>
#include <crv\zerocurv.h>
#include <ARM\libarm_local\ARM_local_wrapper.h>
#include "GP_Inflation\gpinflation\infswaption.h"

using ARM::ARM_InfSwaption;



/////////////////////////////////////////////////////
/// find that corresponding type depending on swaption
/////////////////////////////////////////////////////
ARM_Swaption* CorrespondingSwaption(ARM_Swap* swap,
									int isRecOrPay, 
									int exerciseType,
									long strikeType,
									double strike,
									const ARM_Date& myEndDate)
{
    if ( strikeType == 1 ) // Variable Strike
	{
		ARM_ReferenceValue* stepUpStrike = NULL;

		stepUpStrike = (ARM_ReferenceValue *)
						 LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpStrike, ARM_REFERENCE_VALUE) == 0 )
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                                 "ARM_ERR: Variable Strike is not of a good type");
		}

		if( ARM::ARM_InfSwaption::NbOfInflationLeg(swap) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "ARM_ERR: Inflation swaption with variable strike not yet implemented");
		else
			return new ARM_Swaption(swap, isRecOrPay, stepUpStrike, exerciseType, (ARM_Date&) myEndDate );

	}
	else
	{
		if( ARM::ARM_InfSwaption::NbOfInflationLeg(swap) )
			return 	new ARM_InfSwaption(swap, isRecOrPay, strike, myEndDate );
		else
			return new ARM_Swaption(swap, isRecOrPay, exerciseType, strike, (ARM_Date&) myEndDate );
	}
}



long ARMLOCAL_SWAPTION(long swapId,
					   long isRecOrPay,
					   long strikeType,
					   double strike,
					   double maturity,
					   long exerciseType,
					   ARM_result& result,
					   long objId)
{
	long swaptionId;

    ARM_Swap* swap=NULL;
    ARM_Swaption* Swaption=NULL; 
    ARM_Swaption* newSwaption=NULL;
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char myEndDate[11];
	
	CCString msg ("");

	try
    {
		Local_XLDATE2ARMDATE(maturity,myEndDate);

		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}

	    newSwaption = CorrespondingSwaption(swap, isRecOrPay, exerciseType, strikeType, strike,
                                           (ARM_Date) myEndDate);

		if (newSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swaptionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption);

			if (swaptionId == RET_KO)
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaptionId);

			return ARM_OK;
		}
		else
		{
			Swaption = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			/// a swaption can be either a swaption or an inflation swaption
			if (	LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Swaption, ARM_SWAPTION)	== 1
				||	LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Swaption, ARM_INFSWAPTION)== 1 )
			{
				if (Swaption)
				{
					delete Swaption;
					Swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
 
    catch(Exception& x)
    {
        x.DebugPrint();

		if (newSwaption)
			delete newSwaption;
		newSwaption = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_SwaptionFromExpiry(CCString& optionExpiry,
							     CCString& swapTerm,
								 long liborType,
								 long receiveOrPay,
								 double strike,
								 long spreadType,
								 double spread,
								 bool ccyIsObject,
								 const CCString& CcyName,
								 ARM_result& result,
								 long objId)
{
	long swaptionId;

	char* sCcy				  = NULL;
	ARM_Currency* ccy		  = NULL;
	ARM_Swaption* swaption    = NULL;
	ARM_Swaption* newSwaption = NULL;

	char* eur = "EUR";


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
	   return(ARM_KO);
	}


	double myStrike = strike;

	CCString msg ("");

	try
	{
		if ( CcyName == "DEFAULT" )
		{
			if (( liborType == K_PIBOR1M ) 
				||
			    ( liborType == K_PIBOR3M ) 
				||
			    (liborType == K_PIBOR6M) 
				||
			    ( liborType == K_PIBOR1Y )
			   )
			{
			   sCcy = eur;
			}
			else
			{
			   sCcy = ARM_DEFAULT_COUNTRY;
			}
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (CcyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
			   result.setMsg ("ARM_ERR: Currency is not of a good type");
				
			   return(ARM_KO);
			}

            sCcy = ccy->GetCcyName();
		}
		else
		{
		   sCcy = const_cast<char*>((const char *) CcyName);
		}
		
		if ( strike < 0 )
		{
		   myStrike =  K_MARKET_RATE;
		}

		if ( spreadType == 1 )  // Variable spread
		{
			ARM_ReferenceValue* stepUpSpread = NULL;
			stepUpSpread = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
			   result.setMsg("ARM_ERR: Variable spread is not of a good type");
				
			   return(ARM_KO);
			}

			newSwaption = new ARM_Swaption((int) receiveOrPay,
				                           (char *) optionExpiry,
										   (char *) swapTerm,
										   (ARM_INDEX_TYPE) liborType,
										   0.0, // the Spread 
										   myStrike,
										   sCcy);

			if (newSwaption)
			   newSwaption->GetFloatLeg()->SetVariableSpread(stepUpSpread);
		}
		else
		{
		   newSwaption = new ARM_Swaption((int) receiveOrPay, (char *) optionExpiry, 
			                              (char *) swapTerm, (ARM_INDEX_TYPE) liborType, 
										  spread, 
										  myStrike, 
										  sCcy);
		}
		
		if ( newSwaption == NULL )
		{
		   result.setMsg("ARM_ERR: swaption is null");
			
		   return ARM_KO;
		}
		
		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swaptionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) newSwaption);

			if ( swaptionId == RET_KO )
			{
				delete newSwaption;
				newSwaption = NULL;
				
				result.setMsg("ARM_ERR : Problem while inserting object");
				
				return ARM_KO;
			}
			result.setLong(swaptionId);
			return ARM_OK;
		}
		else
		{
			swaption = (ARM_Swaption*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swaption, ARM_SWAPTION) == 1)
			{
				delete swaption;
				swaption = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) newSwaption, objId);
			
				return ARM_OK;
			}
			else
			{
				delete newSwaption;
				newSwaption = NULL;
				
				result.setMsg("ARM_ERR : Previous object is not of a good type");
				
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newSwaption)
			delete newSwaption;
		newSwaption = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_LIBORSWAPTION (double startDate,
							 double endDate,
							 long receiveOrPay,
							 double strike,
							 double maturity,
							 long liborType,
							 long spreadType,
							 double spread,
							 long exerciseType,
							 long resetFreq,
							 long payFreq,
							 bool ccyIsObject,
							 const CCString& CcyName,
							 ARM_result& result,
							 long objId)
{
	long swaptionId;

	ARM_Currency* ccy=NULL;
	ARM_Swaption* Swaption=NULL;
	ARM_Swaption* newLiborSwaption=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char myStartDate[11];
	char myEndDate[11];
	char myMaturityDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,myStartDate);
 		Local_XLDATE2ARMDATE(endDate,myEndDate);
 		Local_XLDATE2ARMDATE(maturity,myMaturityDate);
 		
			if(CcyName == "DEFAULT")
		{
			if((liborType == K_PIBOR1M) ||
			   (liborType == K_PIBOR3M) ||
			   (liborType == K_PIBOR6M) ||
			   (liborType == K_PIBOR1Y))
			{
				ccy = ARM_FRF_CURRENCY;
			}
			else
			{
				ccy = ARM_DEFAULT_CURRENCY;
			}
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (CcyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char*) CcyName);
		}

		if( spreadType == 1) //Variable Spread
		{
			ARM_ReferenceValue* stepUpSpread = NULL;

			stepUpSpread = (ARM_ReferenceValue *)
                             LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
			//	if (newLiborSwaption)
			//		delete newLiborSwaption;
			//	newLiborSwaption = NULL;

				result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
                return(ARM_KO);
			}

			newLiborSwaption = new ARM_Swaption((ARM_Date) myStartDate,
											(ARM_Date) myEndDate,
											receiveOrPay, exerciseType,
											strike,
											(ARM_Date) myMaturityDate,
											(ARM_INDEX_TYPE) liborType, 0.0,
											0., resetFreq,
											payFreq, ccy);

			if (newLiborSwaption)
				newLiborSwaption->GetFloatLeg()->SetVariableSpread(stepUpSpread);
			
		}
		else
		{
			newLiborSwaption = new ARM_Swaption((ARM_Date) myStartDate,
											(ARM_Date) myEndDate,
											receiveOrPay, exerciseType,
											strike,
											(ARM_Date) myMaturityDate,
											(ARM_INDEX_TYPE) liborType, spread,
											0., resetFreq,
											payFreq, ccy);
		}

		if ( !(ccyIsObject)
			&&
			!(CcyName == "DEFAULT") )
		{
			delete ccy;
			ccy = NULL;
		}

		if (newLiborSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: liborSwaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swaptionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLiborSwaption);

			if (swaptionId == RET_KO)
			{
				if (newLiborSwaption)
					delete newLiborSwaption;
				newLiborSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaptionId);

			return ARM_OK;
		}
		else
		{
			Swaption = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Swaption, ARM_SWAPTION) == 1)
			{
				if (Swaption)
				{
					delete Swaption;
					Swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLiborSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (newLiborSwaption)
					delete newLiborSwaption;
				newLiborSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if ( ccyIsObject )
		{
			delete ccy;
			ccy = NULL;
		}

		if (newLiborSwaption)
			delete newLiborSwaption;
		newLiborSwaption = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_EXOSWAPTION (long swapId,
						   long isRecOrPay,
						   long xStyleId,
						   long kRefValId,
						   double swapYearTerm,
						   ARM_result& result,
						   long objId)
{
	long swaptionId;

	ARM_Swap* swap=NULL;
    ARM_Swaption* swaption=NULL;
	ARM_Swaption* newSwaption=NULL;
	ARM_ExerciseStyle *xStyle=NULL;
	ARM_ReferenceValue *kRefVal=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}

		xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(xStyleId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: xStyle is not of a good type");
			return ARM_KO;
		}

		kRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(kRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		newSwaption = new ARM_Swaption(swap, isRecOrPay, xStyle, kRefVal,
										swapYearTerm);

		if (newSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: Exotic Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swaptionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption);

			if (swaptionId == RET_KO)
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaptionId);

			return ARM_OK;
		}
		else
		{
			swaption = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swaption, ARM_SWAPTION) == 1)
			{
				if (swaption)
				{
					delete swaption;
					swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwaption)
			delete newSwaption;
		newSwaption = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_VARFIXSWAPTION (double startDate,
							  double endDate,
							  long spreadsId,
							  long exStyleId,
							  long receiveOrPay,
							  double strike,
							  double maturity,
							  long liborType,
							  double spread,
							  long resetFreq,
							  long payFreq,
							  long ccyId,
							  ARM_result& result,
							  long objId)
{
	long swtionId;

	ARM_Currency* ccy=NULL;
	ARM_Swaption* createdSwaption=NULL;
	ARM_Swaption* swaption=NULL;

	ARM_ExerciseStyle*  xStyle  = NULL;
	ARM_ReferenceValue* spreads = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];
	char sMaturity[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);
 		Local_XLDATE2ARMDATE(maturity,sMaturity);

		if ( ccyId == ARM_NULL_OBJECT )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyId == ARM_FRF_CCY_OBJECT )
		{
		   ccy = ARM_FRF_CURRENCY;
		}
		else
		{
		   ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		spreads = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spreads, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: spreads is not of a good type");
			return ARM_KO;
		}

		xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(exStyleId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: xStyle is not of a good type");
			return ARM_KO;
		}

		createdSwaption = new ARM_Swaption((ARM_Date) sStartDate,
										   (ARM_Date) sEndDate,
										   spreads, xStyle,
										   receiveOrPay,
										   strike,
										   (ARM_Date) sMaturity,
										   (ARM_INDEX_TYPE) liborType,
										   spread,
										   0., resetFreq,
										   payFreq, ccy);

		if (createdSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: VarFix Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swtionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwaption);

			if (swtionId == RET_KO)
			{
				if (createdSwaption)
					delete createdSwaption;
				createdSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swtionId);

			return ARM_OK;
		}
		else
		{
			swaption = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swaption, ARM_SWAPTION) == 1)
			{
				if (swaption)
				{
					delete swaption;
					swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (createdSwaption)
					delete createdSwaption;
				createdSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdSwaption)
			delete createdSwaption;
		createdSwaption = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (double startDate,
										 double endDate,
										 double strike,
										 long nbCurPerforAcc,
										 long payFreqId,
										 long ccyId,
										 ARM_result& result,
										 long objId)
{
	long swtionId;

	ARM_Currency* ccy=NULL;
	ARM_ZCOptSwapLeg* createdLiborSwaption=NULL;
	ARM_ZCOptSwapLeg* swaption=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);

		if ( ccyId == ARM_NULL_OBJECT )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyId == ARM_FRF_CCY_OBJECT )
		{
		   ccy = ARM_FRF_CURRENCY;
		}
		else
		{
		   ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		createdLiborSwaption = new ARM_ZCOptSwapLeg((ARM_Date) sStartDate,
                                                    (ARM_Date) sEndDate,
                                                     strike,nbCurPerforAcc,
													 K_RCV, payFreqId,
                                                     ccy);

		if (createdLiborSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swtionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdLiborSwaption);

			if (swtionId == RET_KO)
			{
				if (createdLiborSwaption)
					delete createdLiborSwaption;
				createdLiborSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swtionId);

			return ARM_OK;
		}
		else
		{
			swaption = (ARM_ZCOptSwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swaption, ARM_OPTIONALACCRUALZC) == 1)
			{
				if (swaption)
				{
					delete swaption;
					swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdLiborSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (createdLiborSwaption)
					delete createdLiborSwaption;
				createdLiborSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdLiborSwaption)
			delete createdLiborSwaption;
		createdLiborSwaption = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_EXOCFSWAPTION (long swapId,
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
	long swtcfId;

	ARM_Swap* swap = NULL;
	ARM_Swaption_CapFloor* Swaption = NULL;
	ARM_Swaption_CapFloor* newSwaption = NULL;
	ARM_ExerciseStyle *xStyle = NULL;
	ARM_ReferenceValue *kSptionRefVal = NULL;
	ARM_ReferenceValue *kCFRefVal = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}

		xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(xStyleId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: XStyle is not of a good type");
			return ARM_KO;
		}

		kSptionRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kSptionRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(kSptionRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Swaption Refvalue is not of a good type");
			return ARM_KO;
		}

		kCFRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kCFloorRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(kCFRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: CapFloor Refvalue is not of a good type");
			return ARM_KO;
		}

		newSwaption = new ARM_Swaption_CapFloor(swap, isRecOrPay, xStyle, 
												kSptionRefVal,
												kCFRefVal, isCapOrFloor, 
												cFloorPosition, IsBarrierCF);

		if (newSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swtcfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption);

			if (swtcfId == RET_KO)
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swtcfId);

			return ARM_OK;
		}
		else
		{
			Swaption = (ARM_Swaption_CapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Swaption, ARM_SWAPTION_CAPFLOOR) == 1)
			{
				if (Swaption)
				{
					delete Swaption;
					Swaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (newSwaption)
					delete newSwaption;
				newSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwaption)
			delete newSwaption;
		newSwaption = NULL;

		ARM_RESULT();
    }
}



long ARMLOCAL_FlexAccretSwaption (double startDate,
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
	long flexSwId;

	ARM_FlexAccretSwaption* newFlexSwaption = NULL;
	ARM_FlexAccretSwaption* oldFlexSwaption = NULL;
	ARM_ExerciseStyle *xStyle = NULL;
	ARM_Currency *ccy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);

		if (exerciseTypeId != ARM_NULL_OBJECT)
		{
			xStyle = (ARM_ExerciseStyle*)LOCAL_PERSISTENT_OBJECTS->GetObject(exerciseTypeId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 0)
			{
				result.setMsg ("ARM_ERR: XStyle is not of a good type");
				return ARM_KO;
			}
		}

		if ( ccyId == ARM_NULL_OBJECT )
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyId == ARM_FRF_CCY_OBJECT)
		{
			ccy = ARM_FRF_CURRENCY;
		}
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		newFlexSwaption = new ARM_FlexAccretSwaption((ARM_Date) sStartDate,
													 (ARM_Date) sEndDate,
													 fixedRate,
													 nbCurPerforAcc,
													 receiveOrPay,
													 freqId,
													 (ARM_INDEX_TYPE)liborTypeId,
													 spread,
													 xStyle,
													 ccy);

		if (newFlexSwaption == NULL)
		{
			result.setMsg ("ARM_ERR: flex Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			flexSwId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFlexSwaption);

			if (flexSwId == RET_KO)
			{
				if (newFlexSwaption)
					delete newFlexSwaption;
				newFlexSwaption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(flexSwId);

			return ARM_OK;
		}
		else
		{
			oldFlexSwaption = (ARM_FlexAccretSwaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFlexSwaption, ARM_FLEXACCRETSWAPTION) == 1)
			{
				if (oldFlexSwaption)
				{
					delete oldFlexSwaption;
					oldFlexSwaption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFlexSwaption, objId);

				return ARM_OK;
			}

			else
			{
				if (newFlexSwaption)
					delete newFlexSwaption;
				newFlexSwaption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newFlexSwaption)
			delete newFlexSwaption;
		newFlexSwaption = NULL;

		ARM_RESULT();
    }
}


double ARMLOCAL_SwaptionStickyDelta(long swaptionId,
									long modelId,
									long perturbeDiscountCurvId,
									ARM_result& result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		ARM_Swaption* origSwaption = NULL;
		if( !GetObjectFromId( &origSwaption, swaptionId, ARM_SWAPTION) )
		{
			result.setMsg ("ARM_ERR: Swaption is not of a good type");
			return ARM_KO;
		};


		ARM_Q_Model* origModel = NULL;
		if( !GetObjectFromId( &origModel, modelId, ARM_Q_MODEL) )
		{
			result.setMsg ("ARM_ERR: Modèle is not of a good type");
			return ARM_KO;
		};

		ARM_ZeroCurve* perturbeDiscountCurv = NULL;
		if( !GetObjectFromId( &perturbeDiscountCurv, perturbeDiscountCurvId, ARM_ZERO_CURVE) )
		{
			result.setMsg ("ARM_ERR: Perturbe Discount Curve is not of a good type");
			return ARM_KO;
		};

		double value = SwaptionStickyDelta(origSwaption,origModel,perturbeDiscountCurv);
                                  
		result.setDouble(value);

		return ARM_OK;

	}
		
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
