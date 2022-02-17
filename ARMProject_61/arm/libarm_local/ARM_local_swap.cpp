#pragma warning(disable : 4541)
#pragma warning(disable : 4250)

#include "firstToBeIncluded.h"
#include <libCCdate\CCdate.h>
#include <libCCtools++\CCstring.h>

#include <math.h>

#include <inst\swaption.h>
#include <inst\swap.h>
#include <inst\cmsleg.h>
#include <inst\reversefloaterswap.h>
#include <inst\corridorleg.h>
#include <inst\cmtleg.h>
#include <inst\t4mleg.h>
#include <inst\fixleg.h>
#include <inst\reversestickyleg.h>

#include <crv\zeroint.h>
#include <ccy\currency.h>
#include <mod\model.h>


#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <util\fromto.h>

#include <ARM\libarm_local\ARM_local_wrapper.h>

#include "GP_Inflation\gpinflation\infleg.h"
using ARM::ARM_AverageSwapLeg;

long ARMLOCAL_LiborAssetSwapMargin (long modelId,
									double startDate,
									double endDate,
									double fixedRate,
									long fixReceiveOrPayId,
									long fixDayCountId,
									long fixFrequencyId,
									long fixDecompFrequencyId,
									long fixPayTimingId,
									long fixIntRuleId,
									long liborTypeId,
									double spread,
									long floatResetFreqId,
									long floatPayFreqId,
									long assetGap,
									long vFlag,
									double price,
									const CCString& discountCCy,
									double redemptionPrice,
									double supplFee,
									long solve,
									char* id,
									ARM_result& result)
{	
	double theMargin (0.);
	ARM_Model* myModel;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char myStartDate[11];
	char myEndDate[11];

	ARM_Currency* aCcy = NULL;
	char* tmp = NULL;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		myModel = (ARM_Model*) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(myModel,ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		}

		ARM_Date dDateCurve = myModel->GetZeroCurve()->GetAsOfDate();

		if (dDateCurve >= (ARM_Date)myEndDate) 
		{
			result.setMsg ("ARM_ERR: the date of the curve is greater than the maturity of the asset swap");
			return ARM_KO;
		}

		tmp = (char*)discountCCy;
		if (tmp)
		{
			aCcy = new ARM_Currency(tmp);
			delete tmp;
		}
		tmp = NULL;

		theMargin = LiborAssetSwapMargin(myModel,
					(ARM_Date)myStartDate, (ARM_Date)myEndDate,
					fixedRate,
					fixReceiveOrPayId,
					fixDayCountId,
					fixFrequencyId,
					fixDecompFrequencyId,
					fixPayTimingId,
					fixIntRuleId,
					(ARM_INDEX_TYPE) liborTypeId,
					spread,
					floatResetFreqId,
					floatPayFreqId,
					assetGap,
					price,
					aCcy,
					vFlag,
					redemptionPrice,
					supplFee,
					solve,
					id);

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		result.setDouble(theMargin);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (tmp)
			delete tmp;
		tmp = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (tmp)
			delete tmp;
		tmp = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}									
}



// PL : Currency can be processed as a String or an ARM object
long ARMLOCAL_LIBORSWAP (double StartDate,
						 double EndDate,
						 long LiborType,
						 long ReceiveOrPay,
						 double FixedRate,
						 long spreadType,
						 double Spread,
						 bool CcyIsObject,
						 const CCString& CcyName,
						 long fixedDayCountId,
						 long floatingDayCountId,
						 ARM_result& result,
						 long objId)
{
	long swapId;

	ARM_Swap* swap = NULL;
	ARM_Swap* newSwap = NULL;
	ARM_Currency* ccy = NULL;

	long resetFreq = -1;
	long payFreq = -1;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char myStartDate[11];
	char myEndDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(StartDate,myStartDate);
		Local_XLDATE2ARMDATE(EndDate,myEndDate);

		if(CcyName == "DEFAULT")
		{
			if((LiborType == K_PIBOR1M) ||
			   (LiborType == K_PIBOR3M) ||
			   (LiborType == K_PIBOR6M) ||
			   (LiborType == K_PIBOR1Y))
			{
				ccy = ARM_FRF_CURRENCY;
			}
			else
			{
				ccy = ARM_DEFAULT_CURRENCY;
			}
		}
		else if (CcyIsObject)
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
		
		if(spreadType == 1) // Variable Spread
		{
			ARM_ReferenceValue* stepUpSpread = NULL;
			stepUpSpread = (ARM_ReferenceValue *)
							 LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(Spread));

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0)
			{
			//	if(newSwap)
			//		delete newSwap;
			//	newSwap = NULL;

				result.setMsg ("ARM_ERR: Variable Spread is not of good type");

				return(ARM_KO);
			}

			newSwap = new ARM_Swap((ARM_Date) myStartDate,
								(ARM_Date) myEndDate,
								(ARM_INDEX_TYPE) LiborType,
								0.0,
								FixedRate,
								ReceiveOrPay,
								resetFreq, payFreq,
								ccy, fixedDayCountId, floatingDayCountId);
			if(newSwap)
				newSwap->GetFloatLeg()->SetVariableSpread(stepUpSpread);
		}
		else //Spread fixe
		{
			newSwap = new ARM_Swap((ARM_Date) myStartDate,
								(ARM_Date) myEndDate,
								(ARM_INDEX_TYPE) LiborType,
								Spread,
								FixedRate,
								ReceiveOrPay,
								resetFreq, payFreq,
								ccy, fixedDayCountId, floatingDayCountId);
		}
        if (!(CcyIsObject)
			&&
			!(CcyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if ( newSwap == NULL )
		{
			result.setMsg ("ARM_ERR: liborSwap is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap);

			if (swapId == RET_KO)
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swapId);

			return ARM_OK;
		}
		else
		{
			swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 1)
			{
				if (swap)
				{
					delete swap;
					swap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap, objId);

				return ARM_OK;
			}
			else
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}									
}



long ARMLOCAL_FIXEDLEG (double startDate,
						double endDate,
						long receiveOrPay,
						long fixedRateType,
						double fixedRate,
						long dayCount,
						long freq,
						long decompFreq,
						long payTiming,
						long intRule,
						long stubRule,
						bool ccyIsObject,
						const CCString& ccyName,
						const CCString& payCalName,
						long nxChange,
						double refDate,
						long adjStartDateId,
						long rollDay,
						ARM_result& result,
						long objId)
{
	long swaplegId;

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

    ARM_SwapLeg* swapLeg    = NULL;
    ARM_SwapLeg* newSwapLeg = NULL;
    ARM_Currency* ccy       = NULL;
   
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
		
        return(ARM_KO);
	}

	char payCalTmp[(ARM_NB_MAX_CAL*3)+1];
    char* payCal = NULL;

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate, "NULL");
		else
		   Local_XLDATE2ARMDATE(refDate, myRefDate);

		if ( ccyName == "DEFAULT" )
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char *) ccyName);
		}

        if (strcmp((const char*)payCalName,"NULL"))
        {
		   strcpy(payCalTmp, (const char *) payCalName);
    
           payCal = payCalTmp;
        }
    
        
		if ( fixedRateType == 0 )
		{
		   newSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
										(ARM_Date) myEndDate,
										fixedRate, receiveOrPay,
										freq, dayCount, decompFreq,
										payTiming, intRule, stubRule,
										ccy, payCal, nxChange,
                                        myRefDate,adjStartDateId,
										rollDay);
		}
		else
		{
			ARM_ReferenceValue* stepUpCoupon;

			stepUpCoupon = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(fixedRate));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpCoupon, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Fixed Rate is not of a good type");
				return ARM_KO;
			}

			newSwapLeg = new ARM_FixLeg((ARM_Date) myStartDate,
									    (ARM_Date) myEndDate,
										stepUpCoupon, 
										receiveOrPay,
										freq, dayCount, decompFreq,
										payTiming, intRule, stubRule,
										ccy, payCal, nxChange,
                                        myRefDate,adjStartDateId);
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Fixed Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (swaplegId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (swapLeg)
			{
			   delete swapLeg;
			   swapLeg = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);

			return(ARM_OK);
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}									
}


long ARMLOCAL_SWAPLEG (long irIndex,
					   double startDate,
					   double endDate,
					   long receiveOrPay,
                       long spreadType,
					   double spread,
					   bool ccyIsObject,
					   const CCString& ccyName,
					   long dayCount,
					   long resetGap,
					   CCString resetCal,
					   CCString payCal,
					   long decompPricingFlag,
					   long nxChange,
					   long stubRuleId,
					   double refDate,
					   long adjStartDateId,
					   long rollDay,
					   ARM_result& result,
					   long objId)
{
	long swaplegId;

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

	ARM_IRIndex* IRI            = NULL;
	ARM_SwapLeg* createdSwapLeg = NULL;
	ARM_SwapLeg* swapLeg        = NULL;
	ARM_Currency* ccy           = NULL;



	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char* resetCalName = new char[100];
	char* payCalName   = new char[100];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate,"NULL");
		else
		   Local_XLDATE2ARMDATE(refDate,myRefDate);

		IRI = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(irIndex);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IRI, ARM_IRINDEX) == 0 )
		{
			result.setMsg("ARM_ERR: IR Index is not of a good type");

			return(ARM_KO);
		}

		if ( strcmp((const char*)resetCal, "NULL") != 0 )
		   strcpy(resetCalName, (const char*)resetCal);
		else
		{
			if (resetCalName)
				delete [] resetCalName;
			resetCalName = NULL;
		}

		if ( strcmp((const char*)payCal, "NULL") != 0 )
		   strcpy(payCalName,(const char*)payCal);
		else
		{
		   if (payCalName)
			  delete [] payCalName;
		   payCalName = NULL;
		}

		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
				result.setMsg("ARM_ERR: Currency is not of a good type");

				return(ARM_KO);
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char *) ccyName);
		}

		if ( spreadType == 0 ) // Constant Spread             
		{
           if (IRI->IsCMSIndex())
           {
              createdSwapLeg = new ARM_CMSLeg((ARM_Date) myStartDate,
											   (ARM_Date) myEndDate,
											   IRI,
											   receiveOrPay,
											   spread,
                                               ccy,
											   stubRuleId,
											   resetCalName,
											   payCalName, 
											   myRefDate,
											   adjStartDateId,
											   dayCount);
           }
           else
           {
			  createdSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
											 (ARM_Date) myEndDate,
											 IRI,
											 receiveOrPay,
											 spread,
											 stubRuleId,
											 K_COMP_PROP,  // TMP n existe pas ds macro
											 ccy, dayCount,
											 resetGap, resetCalName,
											 payCalName, decompPricingFlag,
											 nxChange,
											 myRefDate,
											 adjStartDateId,
											 rollDay);
           }
        }
        else // Variable Coupon 
        {
			if (IRI->IsCMSIndex())
			{
				ARM_ReferenceValue* stepUpSpread = NULL;

				stepUpSpread = (ARM_ReferenceValue *)
								LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

				if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
				{
					if (resetCalName)
					   delete [] resetCalName;
					resetCalName = NULL;

					if (payCalName)
					   delete [] payCalName;
					payCalName = NULL;

					result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
		            return(ARM_KO);
				}

				createdSwapLeg = new ARM_CMSLeg((ARM_Date) myStartDate,
											   (ARM_Date) myEndDate,
											   IRI,
											   receiveOrPay,
											   0.0,
											   ccy,
											   stubRuleId,
											   resetCalName,
											   payCalName, 
											   myRefDate,
											   adjStartDateId,
											   dayCount);

				createdSwapLeg->SetVariableSpread(stepUpSpread);
			}
			else if ( IRI->GetIndexType() == (ARM_INDEX_TYPE) IDXFIXED ) // Fixed Leg
			{
				ARM_ReferenceValue* stepUpCoupon = NULL;

				stepUpCoupon = (ARM_ReferenceValue *)
								LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

				if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpCoupon, ARM_REFERENCE_VALUE) == 0 )
				{
					if (resetCalName)
					   delete [] resetCalName;
					resetCalName = NULL;

					if (payCalName)
					   delete [] payCalName;
					payCalName = NULL;

					result.setMsg("ARM_ERR: Variable Coupon is not of a good type");
				
		            return(ARM_KO);
				}
            
				ARM_SwapLeg fixedTmpLeg((ARM_Date) myStartDate,
										(ARM_Date) myEndDate,
										IRI,
										receiveOrPay,
										0.0, // spread
										stubRuleId,
										K_COMP_PROP,  // TMP n existe pas ds macro
										ccy, dayCount,
										resetGap, resetCalName,
										payCalName, decompPricingFlag,
										nxChange,
										myRefDate,
										adjStartDateId,
										rollDay);

				ARM_FixLeg* fixLeg = new ARM_FixLeg();

				fixLeg->SetVarCoupons(stepUpCoupon);

				createdSwapLeg = (ARM_SwapLeg *) fixLeg;
           
				*createdSwapLeg = fixedTmpLeg;
			}
			else // Floating Leg
			{
				ARM_ReferenceValue* stepUpSpread = NULL;

				stepUpSpread = (ARM_ReferenceValue *)
								LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

				if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
				{
					if (resetCalName)
					   delete [] resetCalName;
					resetCalName = NULL;

					if (payCalName)
					   delete [] payCalName;
					payCalName = NULL;

					result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
		            return(ARM_KO);
				}

				createdSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
												 (ARM_Date) myEndDate, 
												 IRI,
												 receiveOrPay,
												 0.0,
												 stubRuleId,
												 K_COMP_PROP,  // TMP n existe pas ds macro
												 ccy, dayCount,
												 resetGap, resetCalName,
												 payCalName, decompPricingFlag,
												 nxChange,
												 myRefDate,
												 adjStartDateId,
												 rollDay);

				createdSwapLeg->SetVariableSpread(stepUpSpread);
			}
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (resetCalName)
		   delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
		   delete [] payCalName;
		payCalName = NULL;

		if ( createdSwapLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Swap Leg is NULL");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg);

			if ( swaplegId == RET_KO )
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return(ARM_OK);
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 1)
			{
				if (swapLeg)
				{
				   delete swapLeg;
				   swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg, objId);
			
				return(ARM_OK);
			}
			else
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
			
                return(ARM_KO);
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();

		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_LIBORLEG(double startDate,
					   double endDate,
					   long liborType,
					   long receiveOrPay,
                       long spreadType,
					   double spread,
					   long resetFreq,
					   long payFreq,
					   long resetTiming,
					   long payTiming,
					   bool ccyIsObject,
					   const CCString& ccyName,
					   long intRuleId,
					   long resetGap,
					   CCString resetCal,
					   CCString payCal,
					   long decompPricingFlag,
					   long nxChange,
                       long stubRuleId,
                       double refDate,
					   long adjStartDateId,
					   int couponDayCount,
                       ARM_result& result,
                       long objId)
{
	long swaplegId;
	
    ARM_SwapLeg* swapLeg = NULL;
    ARM_SwapLeg* newSwapLeg = NULL;
    ARM_Currency* ccy = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

	char* resetCalName = new char[100];
	char* payCalName = new char[100];

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate,"NULL");
		else
		   Local_XLDATE2ARMDATE(refDate, myRefDate);

		if(ccyName == "DEFAULT")
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
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);
 
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
				if (resetCalName)
					delete [] resetCalName;
				resetCalName = NULL;

				if (payCalName)
					delete [] payCalName;
				payCalName = NULL;

				result.setMsg ("ARM_ERR: Currency is not of a good type");
				
                return(ARM_KO);
			}
        }
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}

		if (strcmp((const char*)resetCal,"NULL") != 0)
		   strcpy(resetCalName,(const char*)resetCal);
		else
		{
			if (resetCalName)
				delete [] resetCalName;
			resetCalName = NULL;
		}

		if ( strcmp((const char*)payCal,"NULL") != 0 )
			 strcpy(payCalName, (const char*) payCal);
		else
		{
			if (payCalName)
				delete [] payCalName;
			payCalName = NULL;
		}

		if ( spreadType == 1 ) // Variable Spread 
        {
			ARM_ReferenceValue* stepUpSpread = NULL;

			stepUpSpread = (ARM_ReferenceValue *)
                             LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
			//	if (newSwapLeg)
			//		delete newSwapLeg;
			//	newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
                return(ARM_KO);
			}

			newSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
                                     (ARM_Date) myEndDate,
                                     (ARM_INDEX_TYPE) liborType,
                                     receiveOrPay,
                                     0.0,
                                     resetFreq,
                                     payFreq,
                                     resetTiming, payTiming,
                                     ccy, intRuleId,resetGap,
								     resetCalName,payCalName,decompPricingFlag,
									 nxChange,stubRuleId,myRefDate,adjStartDateId,
									 couponDayCount,
									 K_MOD_FOLLOWING,10000);

			if (newSwapLeg)
			   newSwapLeg->SetVariableSpread(stepUpSpread);
	
        }
        else// The default: Libor leg with a constant spread
		{
        newSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
                                     (ARM_Date) myEndDate,
                                     (ARM_INDEX_TYPE) liborType,
                                     receiveOrPay,
                                     spread,
                                     resetFreq,
                                     payFreq,
                                     resetTiming, payTiming,
                                     ccy, intRuleId,resetGap,
								     resetCalName,payCalName,decompPricingFlag,
									 nxChange,stubRuleId,myRefDate,adjStartDateId,couponDayCount,
									 K_MOD_FOLLOWING, 10000);
		}

		if (!(ccyIsObject)
			&& !(ccyName == "DEFAULT"))
		{
			delete ccy;
			ccy = NULL;
		}

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		
		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Libor Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (swaplegId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_SWAPLEG) == 1)
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);
			
				return ARM_OK;
			}
			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_SWAP_PRICE_TO_RATE (long sId,
								  double Date,
								  double Price,
								  long modId,
								  ARM_result& result)
{
    double rate;
    ARM_Swap* swap=NULL;
    ARM_Model* mod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char mydate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(Date,mydate);

		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(sId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0 &&
			LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAPTION) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}
 
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

        swap->SetModelVariable(NULL);
		swap->SetModel(mod);
        
		rate = swap->PriceToRate((ARM_Date) mydate, Price);
		
		result.setDouble(rate);
		
		return ARM_OK;
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
		
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_SWAP_RATE_TO_PRICE (long sId,
								  double date,
								  double rate,
								  long modId,
								  ARM_result& result)
{
    double price;
    ARM_Swap* swap=NULL;
    ARM_Model* mod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char mydate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,mydate);

		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(sId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}
 
		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);
		
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		swap->SetModel(mod);

        price = swap->RateToPrice((ARM_Date) mydate, rate);
		
		result.setDouble(price);
		
		return ARM_OK;
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
		
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_SWAP (long idSwapLeg1,
					long idSwapLeg2,
					long minPay,
					VECTOR<double>& fixedRates, // PL : if there are past fixings to take into account
					ARM_result& result,
					long objId)
{
	long swapId;

	ARM_Swap* swap=NULL;
	ARM_Swap* newSwap=NULL;
	ARM_SwapLeg* swapLeg1=NULL;
	ARM_SwapLeg* swapLeg2=NULL; 

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg1 = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(idSwapLeg1);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg1) == 0)
		{
			result.setMsg ("ARM_ERR: Swapleg1 is not a Swap Leg");
			return ARM_KO;
		}

		swapLeg2 = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(idSwapLeg2);
 
		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg2) == 0)
		{
			result.setMsg ("ARM_ERR: Swapleg2 is not a Swap Leg");
			return ARM_KO;
		}

		newSwap = new ARM_Swap(swapLeg1, swapLeg2, minPay);

		if (newSwap == NULL)
		{
			result.setMsg ("ARM_ERR: Swap is null");
			return ARM_KO;
		}
		
		// PL : manage in advance fixings (add every past fixings actually)
		if ( fixedRates.size() > 0 )
		{
		   ARM_Vector* fixingVect = CreateARMVectorFromVECTOR(fixedRates);
			   
		   // add a vector of fixings
		   // Rks : no effect if vector's size < nb of past fixings, since SetFixRates copies from flow 0 !
		   newSwap->SetFixRates(fixingVect);

		   delete fixingVect;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap);

			if (swapId == RET_KO)
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swapId);

			return ARM_OK;
		}
		else
		{
			swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 1)
			{
				if (swap)
				{
					delete swap;
					swap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap, objId);

				return ARM_OK;
			}
			else
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

extern long ARMLOCAL_SwapFromExpiry (CCString& expiry,
									 CCString& tenor,
									 long liborType,
									 double fixedRate,
									 bool ccyIsObject,
									 const CCString& CcyName,
									 ARM_result& result,
									 long objId = -1)
{
	
	long swapId;
	ARM_Swap* swap = NULL;
	ARM_Swap* newSwap = NULL;
	char* sCcy = NULL;
	ARM_Currency* ccy = NULL;
	char* sExpiry;
	char* sTenor;
	char* eur = "EUR";

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

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
		
		sExpiry = const_cast<char*>((const char*) expiry);
		
		sTenor = const_cast<char*>((const char*) tenor);
		
		newSwap = new ARM_Swap(sExpiry, sTenor, (ARM_INDEX_TYPE) liborType, fixedRate, sCcy); 

		if (newSwap == NULL)
		{
			result.setMsg ("ARM_ERR: Swap is null");
			return ARM_KO;
		}
		
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			swapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap);

			if (swapId == RET_KO)
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swapId);

			return ARM_OK;
		}
		else
		{
			swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 1)
			{
				if (swap)
				{
					delete swap;
					swap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwap, objId);

				return ARM_OK;
			}
			else
			{
				if (newSwap)
					delete newSwap;
				newSwap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwap)
			delete newSwap;
		newSwap = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


extern long ARMLOCAL_GETSWAPLEGFROMSWAP(
        long swapId,
        int legNumber,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_Swap* swap = NULL;
    ARM_Object* object=NULL;

	try
	{
		swap = (ARM_Swap *)(LOCAL_PERSISTENT_OBJECTS->GetObject(swapId));

		if (!swap)
		{
			result.setMsg ("ARM_ERR: Swap object is not of a good type");
			return ARM_KO;
		}

		if (legNumber == 1)
		{
			object = (ARM_SwapLeg*)(swap->Get1stLeg()->Clone());
		}
		else if (legNumber == 2)
		{
			object = (ARM_SwapLeg*)((swap->Get2ndLeg())->Clone());
		}

		/// Assign the object in ARM cache
		if(!assignObject(object, result, objId))
		{
			delete object;
			object = NULL;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		object = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_CMSLEG (double startDate,
					  double endDate,
					  long cmsTypeId,
					  long receiveOrPay,
					  long spreadType,
					  double spread,
					  long yieldDecompFreq,
					  long swapLegDayCount,
					  long resetFreq,
					  long payFreq,
					  long resetGap,		 
					  long intRule,
					  bool ccyIsObject,
					  const CCString& ccyName,
					  long resetTiming,
					  long stubRule,
					  long adjStartDate,
					  ARM_result& result,
					  long objId)
{
	long legId;

	ARM_CMSLeg* swapLeg=NULL;
	ARM_CMSLeg* newSwapLeg=NULL;
	ARM_Currency* ccy=NULL;
	ARM_INDEX_TYPE ArmCmsType;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ArmCmsType = (ARM_INDEX_TYPE) cmsTypeId;

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		if (spreadType == 0)
		{
			newSwapLeg = new ARM_CMSLeg((ARM_Date) sStartDate,
										 (ARM_Date) sEndDate,
										  ArmCmsType,
										  receiveOrPay,
										  spread,
										  yieldDecompFreq,
										  swapLegDayCount,
										  intRule,
										  resetFreq, 
										  payFreq,
										  resetGap,
										  ccy, resetTiming,stubRule,
										  NULL, NULL,
										  adjStartDate);
		}
		else
		{
			ARM_ReferenceValue* stepUpSpread;

			stepUpSpread = (ARM_ReferenceValue *)
							LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
				result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
			
		        return(ARM_KO);

			}

			newSwapLeg = new ARM_CMSLeg((ARM_Date) sStartDate,
										 (ARM_Date) sEndDate,
										  ArmCmsType,
										  receiveOrPay,
										  stepUpSpread,
										  yieldDecompFreq,
										  swapLegDayCount,
										  intRule,
										  resetFreq, 
										  payFreq,
										  resetGap,
										  ccy, resetTiming,stubRule,
										  NULL, NULL,
										  adjStartDate);
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: CMS Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (legId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(legId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_CMSLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_CMSLEG) == 1)
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);

				return ARM_OK;
			}
			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
        x.DebugPrint();

		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_REVERSEFLOAT(double startDate,
						   double endDate,
						   double NotionalRatio,
						   long receiveOrPay,
						   long liborType,
						   const CCString& sFltccy,
						   long fltSpreadsId,
						   long fltDayCount,
						   long fxCouponsId,
						   long extraCouponsId,
						   const CCString& sReverseCcy,
						   long payFreq,
						   long reverseIndexId,
						   long fxDayCount,
						   double multiplier,
						   double couponFloor,
						   double dStubDate,
						   ARM_result& result,
						   long objId)
{
	long revId;

	ARM_ReverseFloaterSwap* reverseFloaterSwap = NULL;
	ARM_ReverseFloaterSwap* createdReverseFloaterSwap = NULL;

	ARM_ReferenceValue* fltSpreads = NULL;
	ARM_ReferenceValue* fixedRates = NULL;
	ARM_ReferenceValue* extraCoupons = NULL;

	ARM_Currency* fltCcy = NULL;
	ARM_Currency* reverseCcy = NULL;

	ARM_Date* stubPointDate = NULL;
	ARM_Date  stubDate;
	ARM_Date defaultDate("01/01/1980");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* tmp = NULL;
	char sStartDate[11];
	char sEndDate[11];
	char sStubDate[11];

	CCString msg ("");

	try
	{
		tmp = (char*)sFltccy;
		if (tmp)
		{
			fltCcy = new ARM_Currency(tmp);
			free(tmp);
		}

		tmp = (char*)sReverseCcy;
		if (tmp)
		{
			reverseCcy = new ARM_Currency(tmp);
			free(tmp);
		}
		tmp = NULL;

		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

    
		if (dStubDate > 0)
			Local_XLDATE2ARMDATE(dStubDate,sStubDate);
		else
			strcpy(sStubDate,"01/01/1970");

		fltSpreads = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fltSpreadsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(fltSpreads, ARM_REFERENCE_VALUE) == 0)
		{
			if (fltCcy)
				delete fltCcy;
			fltCcy = NULL;

			if (reverseCcy)
				delete reverseCcy;
			reverseCcy = NULL;

			result.setMsg ("ARM_ERR: fltSpreads is not of a good type");
			return ARM_KO;
		}

		fltSpreads->SetCalcMethod(K_DISCRETE_REF);

		fixedRates = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fxCouponsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(fixedRates, ARM_REFERENCE_VALUE) == 0)
		{
			if (fltCcy)
				delete fltCcy;
			fltCcy = NULL;

			if (reverseCcy)
				delete reverseCcy;
			reverseCcy = NULL;

			result.setMsg ("ARM_ERR: fixedRates is not of a good type");
			return ARM_KO;
		}

		fixedRates->SetCalcMethod(K_DISCRETE_REF);

		extraCoupons = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(extraCouponsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(extraCoupons, ARM_REFERENCE_VALUE) == 0)
		{
			if (fltCcy)
				delete fltCcy;
			fltCcy = NULL;

			if (reverseCcy)
				delete reverseCcy;
			reverseCcy = NULL;

			result.setMsg ("ARM_ERR: extraCoupons is not of a good type");
			return ARM_KO;
		}

		extraCoupons->SetCalcMethod(K_DISCRETE_REF);

		stubDate = (ARM_Date) sStubDate;
        if (stubDate > defaultDate)
            stubPointDate = &stubDate;
        else
            stubPointDate = NULL;

		createdReverseFloaterSwap = new ARM_ReverseFloaterSwap((ARM_Date) sStartDate,
															   (ARM_Date) sEndDate,
															   NotionalRatio,
															   receiveOrPay,
											  (ARM_INDEX_TYPE) liborType,
															   *fltCcy,
															   fltSpreads,
															   fltDayCount,
															   fixedRates,
															   extraCoupons,
															   *reverseCcy,
															   payFreq,
											  (ARM_INDEX_TYPE) reverseIndexId,
															   fxDayCount,
															   multiplier,
															   couponFloor,
															   stubPointDate);

		if (fltCcy)
			delete fltCcy;
		fltCcy = NULL;

		if (reverseCcy)
			delete reverseCcy;
		reverseCcy = NULL;

		if (createdReverseFloaterSwap == NULL)
		{
			result.setMsg ("ARM_ERR: Reverse Floater Swap is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			revId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdReverseFloaterSwap);

			if (revId == RET_KO)
			{
				if (createdReverseFloaterSwap)
					delete createdReverseFloaterSwap;
				createdReverseFloaterSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(revId);

			return ARM_OK;
		}
		else
		{
			reverseFloaterSwap = (ARM_ReverseFloaterSwap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(reverseFloaterSwap, ARM_REVERSEFLOATERSWAP) == 1)
			{
				if (reverseFloaterSwap)
				{
					delete reverseFloaterSwap;
					reverseFloaterSwap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdReverseFloaterSwap, objId);

				return ARM_OK;
			}
			else
			{
				if (createdReverseFloaterSwap)
					delete createdReverseFloaterSwap;
				createdReverseFloaterSwap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdReverseFloaterSwap)
			delete createdReverseFloaterSwap;
		createdReverseFloaterSwap = NULL;

		if (fltCcy)
			delete fltCcy;
		fltCcy = NULL;

		if (reverseCcy)
			delete reverseCcy;
		reverseCcy = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdReverseFloaterSwap)
			delete createdReverseFloaterSwap;
		createdReverseFloaterSwap = NULL;

		if (fltCcy)
			delete fltCcy;
		fltCcy = NULL;

		if (reverseCcy)
			delete reverseCcy;
		reverseCcy = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_CORRIDORLEG(double startDate,
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
						  bool ccyIsObject,
						  const CCString& ccyName,
						  long MCFreq,
						  long MCInterp,
						  long LDPricingMethod,
                          long decompPricingFlag,
						  CCString resetCal,
						  CCString payCal,
                          long refIndexResetGap,
						  ARM_result& result,
						  long objId)
{
	long legId;

	ARM_CorridorLeg* createdCorrLeg		= NULL;
	ARM_CorridorLeg* precCorrLeg		= NULL;
	ARM_IRIndex* PayIdx					= NULL;
	ARM_IRIndex* RefIdx					= NULL;
	ARM_Currency*  ccy					= NULL;
	ARM_ReferenceValue* BarrierDown		= NULL;
	ARM_ReferenceValue* BarrierUp		= NULL;
	ARM_ReferenceValue* SpreadRefvalue	= NULL;

    ARM_ReferenceValue  BarrierDownTmp(0.0);
    ARM_ReferenceValue  BarrierUpTmp(1000.0);


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg("ARM_ERR: Pb with accessing objects");
	
        return(ARM_KO);
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		PayIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId);
		
        if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PayIdx, ARM_IRINDEX) == 0)
		{
		   result.setMsg ("ARM_ERR: Pay Index is not of a good type");
			
           return(ARM_KO);
		}

		RefIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(refIndexId);
		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RefIdx, ARM_IRINDEX) == 0 )
		{
		   result.setMsg ("ARM_ERR: Ref Index is not of a good type");
			
           return(ARM_KO);
		}

		if ( levelDownId == ARM_NULL_OBJECT )
		{
			BarrierDown = &BarrierDownTmp;
		}
		else
		{
			BarrierDown = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(levelDownId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierDown, ARM_REFERENCE_VALUE) == 0 )
			{
			   result.setMsg ("ARM_ERR: Barrier Down is not of a good type");
				
               return(ARM_KO);
			}
		}

		if ( levelUpId == ARM_NULL_OBJECT )
		{
		   BarrierUp = &BarrierUpTmp;
		}
		else
		{
		   BarrierUp = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(levelUpId);
			
           if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierUp, ARM_REFERENCE_VALUE) == 0 )
           {
			  result.setMsg("ARM_ERR: Barrier Up is not of a good type");
			
              return(ARM_KO);
           }
		}

		if ( ccyName == "DEFAULT" )
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			
            ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
			   result.setMsg("ARM_ERR: Currency is not of a good type");
				
               return(ARM_KO);
			}
		}
		else
		{
		   ccy = new ARM_Currency((const char *) ccyName);
		}

		char resetCalNameBuf[100];
		char payCalNameBuf[100];

        char* resetCalName = NULL;
        char* payCalName   = NULL;

		if ( strcmp((const char *) resetCal, "NULL") != 0 )
        {
		   strcpy(resetCalNameBuf, (const char *) resetCal);

           resetCalName = resetCalNameBuf;
        }

		if ( strcmp((const char *) payCal, "NULL") != 0 )
        {
		   strcpy(payCalNameBuf,(const char *) payCal);

           payCalName = payCalNameBuf;
        }


		if ( spreadType == 0 ) 
		{
			createdCorrLeg = new ARM_CorridorLeg((ARM_Date) sStartDate,
												 (ARM_Date) sEndDate,
												 receiveOrPay,
												 PayIdx,
												 payFreq,
												 spread,
												 RefIdx,
												 resetFreq,
												 paidRateResetTiming,
												 refRateResetTiming,
												 stubRule,
												 BarrierDown,
												 downSpec,
												 BarrierUp,
												 upSpec,
												 ccy,
												 MCFreq,
												 MCInterp,
												 LDPricingMethod,
                                                 decompPricingFlag,
												 resetCalName,
												 payCalName,
                                                 int(refIndexResetGap));
		}
		else
		{
			SpreadRefvalue = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(spread);
			
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(SpreadRefvalue, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Spread is not of a good type");
				return ARM_KO;
			}

			createdCorrLeg = new ARM_CorridorLeg((ARM_Date) sStartDate,
												 (ARM_Date) sEndDate,
												 receiveOrPay,
												 PayIdx,
												 payFreq,
												 SpreadRefvalue,
												 RefIdx,
												 resetFreq,
												 paidRateResetTiming,
												 refRateResetTiming,
												 stubRule,
												 BarrierDown,
												 downSpec,
												 BarrierUp,
												 upSpec,
												 ccy,
												 MCFreq,
												 MCInterp,
												 LDPricingMethod,
                                                 decompPricingFlag,
												 resetCalName,
												 payCalName,
                                                 int(refIndexResetGap));
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;

		   ccy = NULL;
		}

		if ( createdCorrLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Corridor Leg is null");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object * ) createdCorrLeg);

			if ( legId == RET_KO )
			{
				if (createdCorrLeg)
				   delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(legId);

			return(ARM_OK);
		}
		else
		{
			precCorrLeg = (ARM_CorridorLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(precCorrLeg, ARM_CORRIDORLEG) == 1 )
			{
				if (precCorrLeg)
				{
				   delete precCorrLeg;
					
                   precCorrLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) createdCorrLeg, objId);

				return(ARM_OK);
			}

			else
			{
				if (createdCorrLeg)
				   delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
				
                return(ARM_KO);
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdCorrLeg)
		{
			delete createdCorrLeg;
			createdCorrLeg = NULL;
		}

		ARM_RESULT();
	}

	catch(...)
	{
		if (createdCorrLeg)
		{
			delete createdCorrLeg;
		
            createdCorrLeg = NULL;
		}

		result.setMsg("ARM_ERR: unrecognized failure");
	
        return(ARM_KO);
	}
}


long ARMLOCAL_CORRIDORLEGWITHFIXINGS(const double& startDate,
									 const double& endDate,
									 const int&	payReceive,
									 const long& payIndexId,
									 const long& payFreq,
									 const double& spread,
									 const long& spreadId,
									 const long& refIndexId,
									 const long& resetFreq,
									 const long& paidRateResetTiming,
									 const long& refRateResetTiming,
									 const long& stubRule,
									 const double& levelDown,
									 const long& levelDownId,
									 const long& downSpec,
									 const double& levelUp,
									 const long& levelUpId,
									 const long& upSpec,
									 const long& ccyId,
									 const double& decompPricingFlag,
									 const double& refIndexResetGap,
									 const long& refPastFixingsId,
									 const long& payPastFixingsId,
									 ARM_result& result,
									 long objId)
{
	//Input checks

	if ( !GlobalPersistanceOk( result ) )
	{
	   return(ARM_KO);
	}

	// Used in the MACRO ARM_RESULT
	CCString msg("");
	ARM_CorridorLeg* corridorLeg		= NULL;
	
	ARM_Currency*		ccy				= NULL;
	ARM_IRIndex*		payIndex		= NULL;
	
    ARM_ReferenceValue* stepUpSpread	= NULL;
	ARM_ReferenceValue  stepUpSpreadTmp(spread);

    ARM_IRIndex*		refIndex		= NULL;

	ARM_ReferenceValue* barrierDown		= NULL;
	ARM_ReferenceValue* barrierUp		= NULL;

    ARM_ReferenceValue  barrierDownTmp;
	ARM_ReferenceValue  barrierUpTmp;

	ARM_ReferenceValue* refPastFixings	= NULL;
	ARM_ReferenceValue* payPastFixings  = NULL;
  
	try
	{
		// StartDate
		char myStartDate[20];
		Local_XLDATE2ARMDATE(startDate,	myStartDate);
	
		// EndDate
		char myEndDate[20];
		Local_XLDATE2ARMDATE(endDate, myEndDate);

		// Currency
		if ( ccyId != ARM_NULL_OBJECT )
        {
		    ccy = dynamic_cast<ARM_Currency *>(LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId));

		    if (!ccy)
		    {
			   result.setMsg("ARM_ERR: Currency is not of a good type");
			    
               return(ARM_KO);
		    }
        }
        else
        {
           result.setMsg("ARM_ERR: Currency is not of a good type");
			
           return(ARM_KO);
        }

		// PayIndex
		if ( payIndexId != ARM_NULL_OBJECT )
        {
		    payIndex = dynamic_cast<ARM_IRIndex *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId));

		    if (!payIndex)
		    {
			    result.setMsg ("ARM_ERR: Payment Index is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            result.setMsg ("ARM_ERR: Payment Index is not of a good type");
			return ARM_KO;
        }

		//RefIndex
		if(refIndexId != ARM_NULL_OBJECT)
        {
		    refIndex = dynamic_cast<ARM_IRIndex *>(LOCAL_PERSISTENT_OBJECTS->GetObject(refIndexId));

		    if (!refIndex)
		    {
			    result.setMsg ("ARM_ERR: Reference Index is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            result.setMsg ("ARM_ERR: Reference Index is not of a good type");
			return ARM_KO;
        }

		//BarrierDown
		if (levelDownId != ARM_NULL_OBJECT)
		{
			barrierDown = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(levelDownId));

		    if (!barrierDown)
		    {
			    result.setMsg ("ARM_ERR: Barrier Down is not of a good type");
			    return ARM_KO;
		    }
		}
		else
		{
            barrierDownTmp = ARM_ReferenceValue(levelDown);

			barrierDown    = &barrierDownTmp;
		}

		// BarrierUp
		if ( levelUpId != ARM_NULL_OBJECT )
		{
			barrierUp = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(levelUpId));

		    if (!barrierUp)
		    {
			    result.setMsg ("ARM_ERR: Barrier Up is not of a good type");
			    return ARM_KO;
		    }
		}
		else
		{
           barrierUpTmp = ARM_ReferenceValue(levelUp);
			
           barrierUp    = new ARM_ReferenceValue(levelUp);
		}

		// RefPastFixings
        if(refPastFixingsId != ARM_NULL_OBJECT)
        {
		    refPastFixings = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(refPastFixingsId));

		    if (!refPastFixings)
		    {
			   result.setMsg ("ARM_ERR: Reference Past Fixings is not of a good type");
			    
               return(ARM_KO);
		    }
        }
      
		// PayPastFixings
        if ( payPastFixingsId != ARM_NULL_OBJECT )
        {
		    payPastFixings = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payPastFixingsId));

		    if (!payPastFixings)
		    {
			    result.setMsg("ARM_ERR: Payment Past Fixings is not of a good type");
			    return(ARM_KO);
		    }
        }

		//StepUpSpread
		if(spreadId != ARM_NULL_OBJECT)
        {
		    stepUpSpread = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(spreadId));

		    if (!stepUpSpread)
		    {
			    result.setMsg ("ARM_ERR: Boosted Fix Rate / Spread is not of a good type");
			    return ARM_KO;
		    }
        }
		else
		{
			stepUpSpread = &stepUpSpreadTmp;
		}

		//Non used default value
		int mcFreq				= K_DEF_FREQ;
		int mcInterp			= K_LINEAR;
		int ldPricingMethod     = K_DIGITALE;

		//Calendars
		char* resetCalName		= NULL;
		char* payCalName		= NULL;

        //Create the ARM_CorridorLeg
		corridorLeg = new ARM_CorridorLeg((ARM_Date) myStartDate,
					 					  (ARM_Date) myEndDate,
										   payReceive,
										   payIndex,
										   payFreq,
										   stepUpSpread,
										   refIndex,
									 	   resetFreq,
										   paidRateResetTiming,
										   refRateResetTiming,
										   stubRule,
										   barrierDown,
										   downSpec,
										   barrierUp,
										   upSpec,
										   ccy,
										   mcFreq,
										   mcInterp,
										   ldPricingMethod,
                                           decompPricingFlag,
										   resetCalName,
										   payCalName,
										   refIndexResetGap,
										   refPastFixings,
										   payPastFixings);

		//Assign Object
		if( !assignObject( corridorLeg, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		if (corridorLeg)
		{
			delete corridorLeg;
			corridorLeg = NULL;
		}

		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_RESTRIKABLELEG(double startDate,
							 double endDate,
							 long receiveOrPay,
							 long payIndexId,
							 long spreadId,
							 double spreadValue,
							 long refIndexId,
							 long stubRule,
							 long rangeId,
							 double rangeValue,
						     bool ccyIsObject,
						     const CCString& ccyName,
							 double Alpha,
							 double Beta,
                             long decompPricingFlag,
							 long AdjStartDateFlag,
							 ARM_result& result,
							 long objId)
{
	long legId;

	ARM_CorridorLeg* createdCorrLeg = NULL;
	ARM_CorridorLeg* precCorrLeg    = NULL;
	ARM_IRIndex*     PayIdx         = NULL;
	ARM_IRIndex*     RefIdx         = NULL;
	ARM_Currency*    ccy            = NULL;
	ARM_Currency   tmpccy;
	ARM_ReferenceValue* spread		= NULL;
	ARM_ReferenceValue tmpspread;
	ARM_ReferenceValue* range		= NULL;
	ARM_ReferenceValue tmprange;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		if ( payIndexId == ARM_NULL_OBJECT )
		{
			PayIdx = NULL;
		}
		else
		{
			PayIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId);
			
            if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PayIdx, ARM_IRINDEX) == 0 )
			{
			   result.setMsg("ARM_ERR: Pay Index is not of a good type");
				
               return(ARM_KO);
			}
		}

		if ( refIndexId == ARM_NULL_OBJECT )
		{
			RefIdx = NULL;
		}
		else
		{
			RefIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(refIndexId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RefIdx, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Ref Index is not of a good type");
				return ARM_KO;
			}
		}

		if(ccyName == "DEFAULT")
		{
			ccy = ARM_DEFAULT_CURRENCY;
			tmpccy = *ccy;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}

			tmpccy = *ccy;
		}
		else
		{
			tmpccy = ARM_Currency((const char *) ccyName);
		}

		if ( spreadId == ARM_NULL_OBJECT )
		{
			tmpspread = ARM_ReferenceValue(spreadValue);
		}
		else
		{
			spread = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spread, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Spread object is not of a good type");
				return(ARM_KO);
			}

             tmpspread = *spread;
		}

		if ( rangeId == ARM_NULL_OBJECT )
		{
			tmprange = ARM_ReferenceValue(rangeValue);
		}
		else
		{
			range = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(rangeId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(range, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Range object is not of a good type");
				return(ARM_KO);
			}

			tmprange = *range;
		}

		createdCorrLeg = new ARM_CorridorLeg((ARM_Date) sStartDate, 
											 (ARM_Date) sEndDate,
											 receiveOrPay,
											 PayIdx, 
											 &tmpspread,
											 RefIdx, stubRule,
											 &tmprange, 
											 &tmpccy, 
                                             Alpha, 
                                             Beta,
                                             decompPricingFlag,
											 AdjStartDateFlag);

		if (createdCorrLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Corridor Leg is null");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCorrLeg);

			if ( legId == RET_KO )
			{
				if (createdCorrLeg)
				   delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(legId);

			return ARM_OK;
		}
		else
		{
			precCorrLeg = (ARM_CorridorLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(precCorrLeg, ARM_CORRIDORLEG) == 1)
			{
				if (precCorrLeg)
				{
					delete precCorrLeg;
					precCorrLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) createdCorrLeg, objId);

				return(ARM_OK);
			}

			else
			{
				if (createdCorrLeg)
				   delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				
                return(ARM_KO);
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdCorrLeg)
			delete createdCorrLeg;
		createdCorrLeg = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdCorrLeg)
			delete createdCorrLeg;
		createdCorrLeg = NULL;

		result.setMsg("ARM_ERR: unrecognized failure");

		return(ARM_KO);
	}
}




long ARMLOCAL_CMTLEG (double startDate,
					  double endDate,
					  long cmtTypeId,
					  long bondCouponFreq,
					  long bondDayCount,
					  long receiveOrPay,
					  double spread,
					  long yieldDecompFreq,
					  long swapLegDayCount,
					  long intRule,
					  long resetGap,
					  long resetFreq,
					  double ntlAmount,
					  bool ccyIsObject,
					  const CCString& ccyName,
					  long resetTiming,
					  long adjStartDate,
					  ARM_result& result,
					  long objId)
{
	long legId;

	ARM_CMTLeg* swapLeg = NULL;
	ARM_CMTLeg* newSwapLeg = NULL;
	ARM_Currency* ccy=NULL;

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

		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		newSwapLeg = new ARM_CMTLeg((ARM_Date) sStartDate,
										 (ARM_Date) sEndDate,
										 (ARM_INDEX_TYPE) cmtTypeId,
										 bondCouponFreq,
										 bondDayCount,
										 receiveOrPay,
										 spread,
										 yieldDecompFreq,
										 swapLegDayCount,
										 intRule,
										 resetGap,
										 resetFreq,
										 ntlAmount,
										 ccy,
										 resetTiming,
										 adjStartDate);

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: CMT leg Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (legId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(legId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_CMTLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_CMTLEG) == 1)
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);

				return ARM_OK;
			}

			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_TMLEG(long tmIxType,
					double startDate,
					double endDate,
					long receiveOrPay,
					double spread,
                    int payfreq,
                    int resetFreq,
                    int intRule,
                    int fwdRule,
					int stubRule,
                    bool ccyIsObject,
                    const CCString& ccyName,
					ARM_result& result,
					long objId)
{
	long legId;

	ARM_TMLeg* swapLeg=NULL;
	ARM_TMLeg* newSwapLeg=NULL;

    ARM_Currency* ccy = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	try
	{
        if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId(ccyName);
		
            ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);
 
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
			   result.setMsg("ARM_ERR: Currency is not of a good type");
				
               return(ARM_KO);
			}
        }
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}

		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

        double EoniaTruncPrec = 1e-5;

        newSwapLeg = new ARM_TMLeg((ARM_INDEX_TYPE) tmIxType,
                                   (ARM_Date) sStartDate,
                                   (ARM_Date) sEndDate,
                                   receiveOrPay,
                                   spread,
                                   payfreq,
                                   resetFreq,
                                   intRule,
                                   fwdRule,
								   stubRule,
                                   EoniaTruncPrec,
                                   ccy);

		if ( newSwapLeg == NULL )
		{
           if (!(ccyIsObject)
			   && 
               !(ccyName == "DEFAULT")
              )
           {
		      delete ccy;
		      ccy = NULL;
           }

		   result.setMsg ("ARM_ERR: TMleg(EONIA...) is null");
			
           return(ARM_KO);
		}

        if (!(ccyIsObject)
			&& 
            !(ccyName == "DEFAULT")
           )
		{
		   delete ccy;
		   ccy = NULL;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if ( legId == RET_KO )
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(legId);

			return(ARM_OK);
		}
		else
		{
			swapLeg = (ARM_TMLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_TMLEG) == 1 )
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);

				return ARM_OK;
			}

			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_IMPLIEDSPREAD (long swapId,
							 long modelId,
							 double price,
							 long leg1Or2,
							 ARM_result& result)
{
	double spread;
	ARM_Swap* swap=NULL;
	ARM_Model* mod=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}

		mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		swap->SetModel(mod);

		spread = swap->ComputeImpliedSpread(price, leg1Or2);

		result.setDouble(spread);
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

long ARMLOCAL_ARM_GenAmortization (long swaplegId,
								   long amortMethodId,
								   long amortFrequency,
								   double amortAmount,
								   long daycountId,
								   double legNotional,
								   double amortRate,
								   double reducedMaturity,
								   long modelId,
								   double cleanup,
								   ARM_result& result,
								   long objId)
{
	long amortId;

	ARM_SwapLeg* swapleg = NULL;
	ARM_ReferenceValue* amort = NULL;
	ARM_ReferenceValue* newAmort = NULL;
	ARM_Model* mod = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		swapleg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swaplegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapleg) == 0)
		{
			result.setMsg ("ARM_ERR: Swap Leg is not of a good type");
			return ARM_KO;
		}
	
		if (modelId != ARM_NULL_OBJECT)
		{
			mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
			{
				result.setMsg ("ARM_ERR: Model is not of a good type");
				return ARM_KO;
			}
		}

		newAmort = swapleg->GenerateAmortization(amortMethodId,
												amortFrequency,
												amortAmount,
												daycountId,
												legNotional,
												amortRate,
												reducedMaturity,
												mod,
												cleanup);

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			amortId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newAmort);

			if (amortId == RET_KO)
			{
				if (newAmort)
					delete newAmort;
				newAmort = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(amortId);

			return ARM_OK;
		}
		else
		{
			amort = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(amort, ARM_REFERENCE_VALUE) == 1)
			{
				if (amort)
				{
					delete amort;
					amort = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newAmort, objId);

				return ARM_OK;
			}

			else
			{
				if (newAmort)
					delete newAmort;
				newAmort = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newAmort)
			delete newAmort;
		newAmort = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (newAmort)
			delete newAmort;
		newAmort = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_SWAP_WITH_NOTIONNAL (long swapLeg1,
								   long swapLeg2,
								   long notId,
								   long minPay,
								   ARM_result& result,
								   long objId)
{
	long swapId;

	ARM_Swap* createdSwap = NULL;
	ARM_Swap* oldSwap = NULL;
	
	ARM_SwapLeg* swleg1 = NULL;
	ARM_SwapLeg* swleg2 = NULL;

	ARM_ReferenceValue* ref = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swleg1 = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLeg1);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swleg1) == 0)
		{
			result.setMsg ("ARM_ERR: swleg1 is not a swapleg");
			return ARM_KO;
		}

		swleg2 = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLeg2);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swleg2) == 0)
		{
			result.setMsg ("ARM_ERR: swleg2 is not a swapleg");
			return ARM_KO;
		}

		ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(notId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ref, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		createdSwap = new ARM_Swap(swleg1,swleg2,minPay);

		createdSwap->Get1stLeg()->SetAmount(ref);
		createdSwap->Get2ndLeg()->SetAmount(ref);

		if (createdSwap == NULL)
		{
			result.setMsg ("ARM_ERR: swap is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwap);

			if (swapId == RET_KO)
			{
				if (createdSwap)
					delete createdSwap;
				createdSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swapId);

			return ARM_OK;
		}
		else
		{
			oldSwap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSwap, ARM_SWAP) == 1)
			{
				if (oldSwap)
				{
					delete oldSwap;
					oldSwap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwap, objId);

				return ARM_OK;
			}

			else
			{
				if (createdSwap)
					delete createdSwap;
				createdSwap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdSwap)
			delete createdSwap;
		createdSwap = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdSwap)
			delete createdSwap;
		createdSwap = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_CUSTOMFSTCPN (long swlegId,
								double date,
								ARM_result& result)
{
	ARM_SwapLeg* swapleg = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(date,sDate);

        ARM_Date myDate(sDate);

		swapleg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swlegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapleg) == 0)
		{
			result.setMsg ("ARM_ERR: Swapleg is not of a good type");
			return ARM_KO;
		}

		swapleg->CustomizeFirstPeriod(myDate);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}




long ARMLOCAL_REVERSESTICKYLEG (double startDate,
								double endDate,
								long liborType,
								long receiveOrPay,
								long resetFreq,
								long payFreq,
								long resetTiming,
								long payTiming,
							    bool ccyIsObject,
							    const CCString& ccyName,
								long kRefValId,
								long spreadRefValId,
								double firstCoupon,
								ARM_result& result,
								long objId)
{
	long swaplegId;
	
    ARM_ReverseStickyLeg* swapLeg = NULL;
    ARM_ReverseStickyLeg* newSwapLeg = NULL;
    ARM_Currency* ccy = NULL;
   	ARM_ReferenceValue *StrikeRefVal=NULL;
   	ARM_ReferenceValue *MarginRefVal=NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char myStartDate[11];
	char myEndDate[11];

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if (ccyName == "DEFAULT")
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
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);
 
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
        }
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		StrikeRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kRefValId);
		MarginRefVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(spreadRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(StrikeRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Strikes is not of a good type");
			return ARM_KO;
		}
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(MarginRefVal, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Margins is not of a good type");
			return ARM_KO;
		}

        newSwapLeg = new ARM_ReverseStickyLeg(
                                        (ARM_Date) myStartDate,
                                        (ARM_Date) myEndDate,
                                        (ARM_INDEX_TYPE) liborType,
                                        receiveOrPay,
                                        MarginRefVal,
                                        StrikeRefVal,
                                        firstCoupon,
                                        resetFreq,
                                        payFreq,
                                        resetTiming, 
                                        payTiming,
                                        ccy);
	
        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Libor Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (swaplegId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_ReverseStickyLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_REVERSESTICKYLEG) == 1)
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);
			
				return ARM_OK;
			}
			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_IMPLIEDSPREADWITHMODELS (long swapId,
									   long modelId1,
									   long modelId2,
									   double price,
									   long leg1Or2,
									   ARM_result& result)
{
	double spread;
	ARM_Swap* swap=NULL;
	ARM_Model* mod1=NULL;
	ARM_Model* mod2=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg (" ");

	try
	{
		swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 0)
		{
			result.setMsg ("ARM_ERR: Swap is not of a good type");
			return ARM_KO;
		}

		mod1 = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId1);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod1, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		mod2 = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modelId2);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod2, ARM_MODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type");
			return ARM_KO;
		}

		spread = swap->ComputeImpliedSpreadWith2Models(price,
													   mod1,
													   mod2,
													   leg1Or2);

		result.setDouble(spread);
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_ARM_SetLastFixing(long securityId,
								double rate,
                                double beforeLastFixing,
								double asOf,
								double resetDate,
								ARM_result& result,
								long objId)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];
	char sResetDate[11];

	ARM_Security* createdSecurity = NULL;
	ARM_Security* oldSecurity = NULL;
	ARM_Security* security = NULL;

	ARM_Date dResetDate;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(asOf, sDate);

		ARM_Date AsOfDate(sDate);

		ARM_Date* eventualResetDate = NULL;

		if ( resetDate > 0 )
		{
			Local_XLDATE2ARMDATE(resetDate,sResetDate);
			dResetDate = ARM_Date(sResetDate);
	        eventualResetDate = &dResetDate;
		}

		security = (ARM_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(securityId);

		// clone the good type of security !
		if ( LocalPersistent::LOCAL_IS_SWAPLEG_OK(security) == 1 )
		{
			createdSecurity = (ARM_SwapLeg *) security->Clone();
		}
		else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CAPFLOOR) == 1 )
		{
			createdSecurity = (ARM_CapFloor *) security->Clone();
		}
		else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_SWAPTION) == 1 )
		{
			createdSecurity = (ARM_Swaption *) security->Clone();
		}
		else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_SWAP) == 1 )
		{
			createdSecurity = (ARM_Swap *) security->Clone();
		}
		else
		{
			result.setMsg ("ARM_ERR: not an ARM security\n(must be a swap, swapleg, capfloor or swaption)");
			return ARM_KO;
		}

		if ( createdSecurity == NULL )
		{
			result.setMsg ("ARM_ERR: error in cloning security");
			return ARM_KO;
		}

		createdSecurity->SetSettlement(AsOfDate);
		createdSecurity->SetLastFixing(rate, beforeLastFixing, eventualResetDate, &AsOfDate);

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			securityId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurity);

			if (securityId == RET_KO)
			{
				if (createdSecurity)
					delete createdSecurity;
				createdSecurity = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(securityId);

			return ARM_OK;
		}
		else
		{
			oldSecurity = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( (LocalPersistent::LOCAL_IS_SWAPLEG_OK(oldSecurity) == 1)
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_CAPFLOOR) == 1 )
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_SWAPTION) == 1 )
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_SWAP) == 1 ) )
			{
				if (oldSecurity)
				{
					delete oldSecurity;
					oldSecurity = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurity, objId);

				return ARM_OK;
			}
			else
			{
				if (createdSecurity)
					delete createdSecurity;
				createdSecurity = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type\n(must be a swap, swapleg, capfloor or swaption)");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdSecurity)
			delete createdSecurity;
		createdSecurity = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdSecurity)
			delete createdSecurity;
		createdSecurity = NULL;

		result.setMsg ("ARM_ERR: ARMLOCAL_ARM_SetLastFixing : unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_ARM_SetFixRates (	long securityId,
								VECTOR<double>& fixedRates,
								ARM_result& result,
								long objId = -1)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Security* createdSecurity = NULL;
	ARM_Security* oldSecurity = NULL;
	ARM_Security* security = NULL;

	CCString msg ("");

	try
	{
		// PL : manage in advance fixings (add every past fixings given in vector)
		if ( fixedRates.size() > 0 )
		{
			ARM_Vector* fixingVect = CreateARMVectorFromVECTOR(fixedRates);

			security = (ARM_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(securityId);

			// clone the good type of security !
			if ( LocalPersistent::LOCAL_IS_SWAPLEG_OK(security) == 1 )
			{
				createdSecurity = (ARM_SwapLeg *) security->Clone();
			}
			else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_CAPFLOOR) == 1 )
			{
				createdSecurity = (ARM_CapFloor *) security->Clone();
			}
			else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_SWAPTION) == 1 )
			{
				createdSecurity = (ARM_Swaption *) security->Clone();
			}
			else if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ARM_SWAP) == 1 )
			{
				createdSecurity = (ARM_Swap *) security->Clone();
			}
			else
			{
				result.setMsg ("ARM_ERR: not an ARM security\n(must be a swap, swapleg, capfloor or swaption)");
				return ARM_KO;
			}

			if ( createdSecurity == NULL )
			{
				result.setMsg ("ARM_ERR: error in cloning security");
				return ARM_KO;
			}

			createdSecurity->SetFixRates(fixingVect);

			delete fixingVect;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			securityId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurity);

			if (securityId == RET_KO)
			{
				if (createdSecurity)
					delete createdSecurity;
				createdSecurity = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(securityId);

			return ARM_OK;
		}
		else
		{
			oldSecurity = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( (LocalPersistent::LOCAL_IS_SWAPLEG_OK(oldSecurity) == 1)
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_CAPFLOOR) == 1 )
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_SWAPTION) == 1 )
				|| ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldSecurity, ARM_SWAP) == 1 ) )
			{
				if (oldSecurity)
				{
					delete oldSecurity;
					oldSecurity = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSecurity, objId);

				return ARM_OK;
			}
			else
			{
				if (createdSecurity)
					delete createdSecurity;
				createdSecurity = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdSecurity)
			delete createdSecurity;
		createdSecurity = NULL;

		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		if (createdSecurity)
			delete createdSecurity;
		createdSecurity = NULL;

		result.setMsg ("ARM_ERR: ARMLOCAL_ARM_SetFixRates: unrecognized failure");
		return ARM_KO;
	}
}



long ARMLOCAL_GENLEG(VECTOR<double> startDates,
					 VECTOR<double> endDates,
					 VECTOR<double> paymentDates,
					 VECTOR<double> resetDates,
					 VECTOR<double> intDays,
					 VECTOR<double> fwdorfixing,
					 long notionalId,
					 long indexId,
					 long receiveOrPayId,
					 long fixrateorspreadType,
					 double fixrateorspread,
					 bool ccyIsObject,
					 const CCString& ccyName,
					 long NxId, 
					 int couponDayCount,
					 ARM_result& result,
					 long objId)
{
	long swaplegId;
	
    ARM_SwapLeg* swapLeg    = NULL;
    ARM_SwapLeg* newSwapLeg = NULL;
    
	ARM_Currency* ccy       = NULL;

    
	ARM_IRIndex* index      = NULL;

    ARM_Currency theCcy;
    

    ARM_ReferenceValue* notional = NULL;

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	int size    = startDates.size();
	int fwdSize = fwdorfixing.size();

	if (( endDates.size() != size )
		|| 
        ( paymentDates.size() != size )
		|| 
        ( resetDates.size() != size )
		|| 
        ( intDays.size() != size )
	   )
	{
		result.setMsg("ARM_ERR: check the length of all your vectors");

		return(ARM_KO);
	}

    ARM_Vector vStartDates(size);
	ARM_Vector vEndDates(size);
	ARM_Vector vPaymentDates(size);
	ARM_Vector vResetDates(size);
	ARM_Vector vIntDays(size);

    ARM_Vector vFwdOrFixing(fwdSize);

	CCString msg("");


    try
    {
		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;

           theCcy = *ccy;
		}
		else if (ccyIsObject)
		{
		   long ccyId = LocalGetNumObjectId(ccyName);

		   ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);
 
		   if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
           {
			  result.setMsg("ARM_ERR: Currency is not of a good type");

			  return(ARM_KO);
           }

           theCcy = *ccy;
        }
		else
		{
           theCcy = ARM_Currency((const char *) ccyName);
		}

        ARM_IRIndex  tmpIndex;
        tmpIndex.SetCurrencyUnit(&theCcy);
        tmpIndex.UpdateCcyDependents(theCcy.GetCcyName());

		if  ( indexId != ARM_NULL_OBJECT )
		{
			index = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(indexId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(index, ARM_IRINDEX) == 0 )
			{
			   result.setMsg("ARM_ERR: IR Index is not of a good type");

			   return(ARM_KO);
			}
		}
        else
        {
           index = &tmpIndex;
        }

		if ( notionalId != ARM_NULL_OBJECT )
		{
			notional = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(notional, ARM_REFERENCE_VALUE) == 0 )
			{
			   result.setMsg ("ARM_ERR: Notional is not of a good type");
			
               return(ARM_KO);
			}
		}
		else
        {
		   notional = new ARM_ReferenceValue(100.0);
        }

		char tmpDate[11];

		for (int i = 0; i < size; i++)
		{
			Local_XLDATE2ARMDATE(startDates[i], tmpDate);
			vStartDates.Elt(i) = ((ARM_Date) tmpDate).GetJulian();
	
            Local_XLDATE2ARMDATE(endDates[i], tmpDate);
			vEndDates.Elt(i) = ((ARM_Date) tmpDate).GetJulian();
			
            Local_XLDATE2ARMDATE(paymentDates[i], tmpDate);
			vPaymentDates.Elt(i) = ((ARM_Date) tmpDate).GetJulian();
			
            Local_XLDATE2ARMDATE(resetDates[i], tmpDate);
			vResetDates.Elt(i) = ((ARM_Date) tmpDate).GetJulian();

			vIntDays.Elt(i) = intDays[i];

			if ( i < fwdSize )
			   vFwdOrFixing.Elt(i) = fwdorfixing[i];
		}

		if ( index->GetIndexType() == IDXFIXED )
		{
			if ( fixrateorspreadType == 0 ) // taux fixe constant         
			{
				newSwapLeg = new ARM_SwapLeg(&vStartDates,
											 &vEndDates,
											 &vPaymentDates,
											 &vResetDates,
											 &vIntDays,
						(( vFwdOrFixing.GetSize() == 0 ) ? NULL : &vFwdOrFixing),
											 notional,
											 index,
											 receiveOrPayId,
											 0.0,
											 fixrateorspread,
											 &theCcy,
											 NxId, 
											 1,
											 couponDayCount);
			}
			else
			{
				ARM_ReferenceValue* stepUpCoupon = NULL;

				stepUpCoupon = (ARM_ReferenceValue *)
							    LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(fixrateorspread));

				if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpCoupon, ARM_REFERENCE_VALUE) == 0)
				{
					result.setMsg("ARM_ERR: stepUpCoupon is not of a good type");
					
                    return(ARM_KO);
				}

				ARM_SwapLeg fixedTmpLeg(&vStartDates,
										&vEndDates,
										&vPaymentDates,
										&vResetDates,
										&vIntDays,
                     (( vFwdOrFixing.GetSize() == 0 ) ? NULL : &vFwdOrFixing),
										notional,
										index,
										receiveOrPayId,
										0.0,
										0.0,
										&theCcy,
										NxId, 
										1,
										couponDayCount);

				ARM_FixLeg* fixLeg = new ARM_FixLeg();

				fixLeg->SetVarCoupons(stepUpCoupon);

				newSwapLeg = (ARM_SwapLeg *) fixLeg;

				*newSwapLeg = fixedTmpLeg;
			}
		}
		else if (!(index->IsCMSIndex())) // Libor Index
		{
			newSwapLeg = new ARM_SwapLeg(&vStartDates,
										 &vEndDates,
										 &vPaymentDates,
										 &vResetDates,
										 &vIntDays,
            (( vFwdOrFixing.GetSize() == 0 ) ? NULL : &vFwdOrFixing),
										 notional,
										 index,
										 receiveOrPayId,
										 (double) fixrateorspread,
										 K_FLOAT_RATE,
										 &theCcy,
										 NxId,
										 1,
										 couponDayCount);

			if ( fixrateorspreadType == 1 ) // Variable Spread 
			{
			   ARM_ReferenceValue* stepUpSpread = NULL;

			   stepUpSpread = (ARM_ReferenceValue *)
							   LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(fixrateorspread));

			   if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0)
               {
				  result.setMsg("ARM_ERR: stepUpSpread is not of a good type");
					
                  return(ARM_KO);
               }

			   newSwapLeg->SetVariableSpread(stepUpSpread);
			}
		}
        else // Index is CMS
        {
           ARM_ReferenceValue* stepUpSpread = NULL;

           ARM_ReferenceValue theSpread((double) fixrateorspread);

 
           if ( fixrateorspreadType == 1 ) // Variable Spread 
           {
			  stepUpSpread = (ARM_ReferenceValue *)
							   LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(fixrateorspread));

			  if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0)
              {
				 result.setMsg("ARM_ERR: stepUpSpread is not of a good type");
					
                  return(ARM_KO);
              }

              theSpread = *stepUpSpread;
           }

           
           newSwapLeg = new ARM_CMSLeg(&vStartDates,
				                       &vEndDates,
				                       &vResetDates,
				                       &vPaymentDates,
				                       &vIntDays,
				                       index->GetIndexType(),
				                       notional,
				                       receiveOrPayId,
				                       &theSpread, 
				                       K_COMP_PROP,          // yieldDecompFreq = K_COMP_PROP, 
				                       index->GetDayCount(), // swapLegDayCount = K30_360,
				                       index->GetIntRule(),
				                       index->GetResetFrequency(),
				                       index->GetPayFrequency(),
				                       index->GetResetGap(),
				                       &theCcy,
				                       index->GetResetTiming());

				                       // int stubRule = K_SHORTSTART,
                                       // char* resetCal = NULL,
                                       // char* payCal = NULL);
        }

		if ( notionalId == ARM_NULL_OBJECT )
		{
		   if (notional)
			  delete notional;
		   notional = NULL;
		}

		if ( newSwapLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Libor Leg is null");
			
           return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if ( swaplegId == RET_KO )
			{
				if (newSwapLeg)
				   delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(swaplegId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

            if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 1 )
			{
				if (swapLeg)
				{
				   delete swapLeg;
				   swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSwapLeg, objId);
			
				return(ARM_OK);
			}
			else
			{
				if (newSwapLeg)
				   delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
				
                return(ARM_KO);
			}
		}

	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		if ( notionalId == ARM_NULL_OBJECT )
		{
		   if (notional)
			  delete notional;
		   notional = NULL;
		}

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
		   delete newSwapLeg;
		newSwapLeg = NULL;

		if ( notionalId == ARM_NULL_OBJECT )
		{
		   if (notional)
			  delete notional;
		   notional = NULL;
		}

		result.setMsg("ARM_ERR: unrecognized failure");
		
        return(ARM_KO);
	}
}



// Be careful this is for ARM_SwapLeg2
long ARMLOCAL_SWAPLEG(long irIndex,
					  double startDate,
					  double endDate,
					  long receiveOrPay,
                      long spreadType,
					  double spread,
					  bool ccyIsObject,
					  const CCString& ccyName,
					  long dayCount,
                      long resetTiming,
					  long resetGap,
                      long payTiming,
                      long payGap,
					  CCString resetCal,
					  CCString payCal,
					  long decompPricingFlag,
					  long nxChange,
					  long stubRuleId,
					  double refDate,
					  long adjStartDateId,
					  ARM_result& result,
					  long objId)
{
    long swaplegId;

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

	ARM_IRIndex* IRI            = NULL;
	ARM_SwapLeg* createdSwapLeg = NULL;
	ARM_SwapLeg* swapLeg        = NULL;
	ARM_Currency* ccy           = NULL;



	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char* resetCalName = new char[100];
	char* payCalName   = new char[100];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate,"NULL");
		else
		   Local_XLDATE2ARMDATE(refDate,myRefDate);

		IRI = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(irIndex);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IRI, ARM_IRINDEX) == 0 )
		{
			result.setMsg("ARM_ERR: IR Index is not of a good type");

			return(ARM_KO);
		}
        IRI->SetResetTiming(resetTiming);
        IRI->SetPayTiming(payTiming);
        IRI->SetResetGap(resetGap);
        IRI->SetPayGap(payGap);
        // update IndexStyle
        if (resetTiming == payTiming && resetTiming == K_ARREARS)
            IRI->SetIndexStyle(IN_ARREARS);
        else 
            IRI->SetIndexStyle(VANILLA);

		if ( strcmp((const char*)resetCal,"NULL") != 0 )
		   strcpy(resetCalName, (const char*)resetCal);
		else
		{
			if (resetCalName)
				delete [] resetCalName;
			resetCalName = NULL;
		}

		if ( strcmp((const char*)payCal,"NULL") != 0 )
		   strcpy(payCalName,(const char*)payCal);
		else
		{
		   if (payCalName)
			  delete [] payCalName;
		   payCalName = NULL;
		}


		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
				result.setMsg("ARM_ERR: Currency is not of a good type");

				return(ARM_KO);
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		if ( spreadType == 0 ) // Constant Spread             
		{
			createdSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
											 (ARM_Date) myEndDate,
											 IRI,
											 receiveOrPay,
											 spread,
											 stubRuleId,
											 K_COMP_PROP,  // TMP n existe pas ds macro
											 ccy, dayCount,
											 resetGap, resetCalName,
											 payCalName, decompPricingFlag,
											 nxChange,
											 myRefDate,
											 adjStartDateId);
        }
        else // Variable Coupon 
        {
			if ( IRI->GetIndexType() == (ARM_INDEX_TYPE) IDXFIXED ) // Fixed Leg
			{
				ARM_ReferenceValue* stepUpCoupon;

				stepUpCoupon = (ARM_ReferenceValue *)
								LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

				if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpCoupon, ARM_REFERENCE_VALUE) == 0 )
				{
					if (resetCalName)
					   delete [] resetCalName;
					resetCalName = NULL;

					if (payCalName)
					   delete [] payCalName;
					payCalName = NULL;

					result.setMsg ("ARM_ERR: Variable Coupon is not of a good type");
				
		            return(ARM_KO);
				}
            
				ARM_SwapLeg fixedTmpLeg((ARM_Date) myStartDate,
										(ARM_Date) myEndDate,
										IRI,
										receiveOrPay,
										0.0, // spread
										stubRuleId,
										K_COMP_PROP,  // TMP n existe pas ds macro
										ccy, dayCount,
										resetGap, resetCalName,
										payCalName, decompPricingFlag,
										nxChange,
										myRefDate,
										adjStartDateId);

				ARM_FixLeg* fixLeg = new ARM_FixLeg();

				fixLeg->SetVarCoupons(stepUpCoupon);

				createdSwapLeg = (ARM_SwapLeg *) fixLeg;
           
				*createdSwapLeg = fixedTmpLeg;
			}
			else // Floating Leg
			{
				ARM_ReferenceValue* stepUpSpread = NULL;

				stepUpSpread = (ARM_ReferenceValue *)
								LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

				if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
				{
					if (resetCalName)
					   delete [] resetCalName;
					resetCalName = NULL;

					if (payCalName)
					   delete [] payCalName;
					payCalName = NULL;

					result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
		            return(ARM_KO);
				}

				createdSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
												 (ARM_Date) myEndDate, 
												 IRI,
												 receiveOrPay,
												 0.0,
												 stubRuleId,
												 K_COMP_PROP,  // TMP n existe pas ds macro
												 ccy, dayCount,
												 resetGap, resetCalName,
												 payCalName, decompPricingFlag,
												 nxChange,
												 myRefDate,
												 adjStartDateId);

				createdSwapLeg->SetVariableSpread(stepUpSpread);
			}
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (resetCalName)
		   delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
		   delete [] payCalName;
		payCalName = NULL;

		if ( createdSwapLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Swap Leg is NULL");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg);

			if ( swaplegId == RET_KO )
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return(ARM_OK);
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_SWAPLEG) == 1)
			{
				if (swapLeg)
				{
				   delete swapLeg;
				   swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg, objId);
			
				return(ARM_OK);
			}
			else
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
			
                return(ARM_KO);
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();

		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

extern long ARMLOCAL_LIVRETALEG (double startDate,
							   double endDate,
							   long liborType,
							   long receiveOrPay,
							   long spreadType,
                               double spread,
							   long resetFreq,
							   long payFreq,
							   long resetTiming,
							   long payTiming,
							   bool ccyIsObject,
							   const CCString& ccyName,
							   long intRuleId,
							   long resetGap,
							   CCString resetCal,
							   CCString payCal,
							   long decompPricingFlag,
							   long nxChange,
							   long stubRuleId,
							   double refDate,
							   long adjStartDateId,
							   long daycountId,
							   ARM_result& result,
							   long objId)
{
	long swaplegId;
	
    ARM_SwapLeg* swapLeg = NULL;
    ARM_SwapLeg* newSwapLeg = NULL;
    ARM_Currency* ccy = NULL;


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

	char* resetCalName = new char[100];
	char* payCalName = new char[100];

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate,"NULL");
		else
		   Local_XLDATE2ARMDATE(refDate, myRefDate);

		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);
 
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
				if (resetCalName)
					delete [] resetCalName;
				resetCalName = NULL;

				if (payCalName)
					delete [] payCalName;
				payCalName = NULL;

				result.setMsg ("ARM_ERR: Currency is not of a good type");
				
                return(ARM_KO);
			}
        }
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		if (strcmp((const char*)resetCal,"NULL") != 0)
		   strcpy(resetCalName,(const char*)resetCal);
		else
		{
			if (resetCalName)
				delete [] resetCalName;
			resetCalName = NULL;
		}

		if ( strcmp((const char*)payCal,"NULL") != 0 )
			 strcpy(payCalName, (const char*) payCal);
		else
		{
			if (payCalName)
				delete [] payCalName;
			payCalName = NULL;
		}

        // The default: Libor leg with a constant spread
        newSwapLeg = new ARM_SwapLeg((ARM_Date) myStartDate,
                                     (ARM_Date) myEndDate,
                                     (ARM_INDEX_TYPE) liborType,
                                     receiveOrPay,
                                     spread,
                                     resetFreq,
                                     payFreq,
                                     resetTiming, payTiming,
                                     ccy, intRuleId,resetGap,
								     resetCalName,payCalName,decompPricingFlag,
									 nxChange,stubRuleId,myRefDate,adjStartDateId,
									 daycountId);

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

        if ( spreadType == 1 ) // Variable Spread 
        {
			ARM_ReferenceValue* stepUpSpread = NULL;

			stepUpSpread = (ARM_ReferenceValue *)
                             LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
                return(ARM_KO);
			}

			newSwapLeg->SetVariableSpread(stepUpSpread);
        }
		
		if (newSwapLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Libor Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg);

			if (swaplegId == RET_KO)
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return ARM_OK;
		}
		else
		{
			swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swapLeg, ARM_SWAPLEG) == 1)
			{
				if (swapLeg)
				{
					delete swapLeg;
					swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSwapLeg, objId);
			
				return ARM_OK;
			}
			else
			{
				if (newSwapLeg)
					delete newSwapLeg;
				newSwapLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (newSwapLeg)
			delete newSwapLeg;
		newSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_AVERAGE_LIBOR_LEG(	long			irIndex,
									double			startDate,
									double			endDate,
									long			receiveOrPay,
									long			spreadType,
									double			spread,
									bool			ccyIsObject,
									const CCString&	ccyName,
									long			dayCount,
									CCString		payCal,
									long			decompPricingFlag,
									long			nxChange,
									long			stubRuleId,
									double			refDate,
									long			adjStartDateId,
									double			couru,
									ARM_result&		result,
									long			objId)
{
	long swaplegId;

	char myStartDate[11];
	char myEndDate[11];
	char myRefDate[11];

	ARM_IRIndex* IRI					= NULL;
	ARM_AverageSwapLeg* createdSwapLeg	= NULL;
	ARM_AverageSwapLeg* swapLeg         = NULL;
	ARM_Currency* ccy					= NULL;



	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
	   result.setMsg("ARM_ERR: Pb with accessing objects");
		
       return(ARM_KO);
	}

	char* resetCalName = new char[100];
	char* payCalName   = new char[100];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,myStartDate);
		Local_XLDATE2ARMDATE(endDate,myEndDate);

		if ( refDate == -1.0 )
		   strcpy(myRefDate,"NULL");
		else
		   Local_XLDATE2ARMDATE(refDate,myRefDate);

		IRI = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(irIndex);

		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IRI, ARM_IRINDEX) == 0 )
		{
			result.setMsg("ARM_ERR: IR Index is not of a good type");

			return(ARM_KO);
		}

		if ( strcmp((const char*)payCal, "NULL") != 0 )
		   strcpy(payCalName,(const char*)payCal);
		else
		{
		   if (payCalName)
			  delete [] payCalName;
		   payCalName = NULL;
		}

		if ( ccyName == "DEFAULT" )
		{
		   ccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccyIsObject )
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0 )
			{
				result.setMsg("ARM_ERR: Currency is not of a good type");

				return(ARM_KO);
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char *) ccyName);
		}

		if ( spreadType == 0 ) // Constant Spread             
		{
          
		  createdSwapLeg = new ARM_AverageSwapLeg(	(ARM_Date) myStartDate,
													(ARM_Date) myEndDate,
													 IRI,
													 receiveOrPay,
													 spread,
													 stubRuleId,
													 K_COMP_PROP,  // TMP n existe pas ds macro
													 ccy, 
													 dayCount,
													 payCalName, 
													 decompPricingFlag,
													 nxChange,
													 myRefDate,
													 adjStartDateId,
													 couru);
        }
        else // Variable Coupon 
        {
			ARM_ReferenceValue* stepUpSpread = NULL;

			stepUpSpread = (ARM_ReferenceValue *)
							LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spread));

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stepUpSpread, ARM_REFERENCE_VALUE) == 0 )
			{
				if (resetCalName)
				   delete [] resetCalName;
				resetCalName = NULL;

				if (payCalName)
				   delete [] payCalName;
				payCalName = NULL;

				result.setMsg ("ARM_ERR: Variable Spread is not of a good type");
				
		           return(ARM_KO);
			}

			createdSwapLeg = new ARM_AverageSwapLeg(	(ARM_Date) myStartDate,
														(ARM_Date) myEndDate, 
														IRI,
														receiveOrPay,
														0.0,
														stubRuleId,
														K_COMP_PROP,  // TMP n existe pas ds macro
														ccy, 
														dayCount,
														payCalName, 
														decompPricingFlag,
														nxChange,
														myRefDate,
														adjStartDateId,
														couru);

			createdSwapLeg->SetVariableSpread(stepUpSpread);
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (payCalName)
		   delete [] payCalName;
		payCalName = NULL;

		if ( createdSwapLeg == NULL )
		{
		   result.setMsg("ARM_ERR: Swap Leg is NULL");

		   return(ARM_KO);
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swaplegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg);

			if ( swaplegId == RET_KO )
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swaplegId);

			return(ARM_OK);
		}
		else
		{
			swapLeg = (ARM_AverageSwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 1)
			{
				if (swapLeg)
				{
				   delete swapLeg;
				   swapLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwapLeg, objId);
			
				return(ARM_OK);
			}
			else
			{
				if (createdSwapLeg)
					delete createdSwapLeg;
				createdSwapLeg = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
			
                return(ARM_KO);
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();

		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		if (createdSwapLeg)
			delete createdSwapLeg;
		createdSwapLeg = NULL;

		if (resetCalName)
			delete [] resetCalName;
		resetCalName = NULL;

		if (payCalName)
			delete [] payCalName;
		payCalName = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ARMLOCAL_SWAP_FROM_SWAPTION (long swaptId,
								  long NDIndex,
								  ARM_result& result,
								  long objId)
{
	long swapId;

    ARM_Swaption* swaption	= NULL;
	ARM_Swap* swap			= NULL;
	ARM_Swap* createdSwap	= NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swaption = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(swaptId);
		
		if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swaption, ARM_SWAPTION) == 0 )
		{
			result.setMsg ("ARM_ERR: Swaption is not of a good type");
			return ARM_KO;
		}

		createdSwap = swaption->GenerateUnderlyingSwapFromBermudaNotice(NDIndex);

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			swapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwap);

			if ( swapId == RET_KO )
			{
				if (createdSwap)
					delete createdSwap;
				createdSwap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(swapId);

			return(ARM_OK);
		}
		else
		{
			swap = (ARM_Swap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(swap, ARM_SWAP) == 1 )
			{
				if (swap)
				{
				   delete swap;
				   swap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSwap, objId);
			
				return(ARM_OK);
			}
			else
			{
				if (createdSwap)
					delete createdSwap;
				createdSwap = NULL;

				result.setMsg("ARM_ERR: previous object is not of a good type");
			
                return(ARM_KO);
			}
		}		
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
		
		ARM_RESULT();
    }

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


/*---- End Of File ----*/

// EOF %M%
