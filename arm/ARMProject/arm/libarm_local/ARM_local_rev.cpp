#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <libCCdate\CCdate.h>
#include <libCCtools++\CCstring.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <ARM\libarm\ARM_result.h>

#include "ARM_local_persistent.h"
#include "ARM_local_wrapper.h"
#include "ARM_local_glob.h"

#include <inst\swapleg.h>
#include <inst\reverse.h>
#include <inst\powrev.h>
#include <gpbase\env.h>

using namespace ARM;

long ARMLOCAL_REVERSE (long structSwapLegId,
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
	long revId;

	ARM_SwapLeg* swapLegClass=NULL;
	ARM_SwapLeg* swapLegStruct=NULL;
	ARM_ReverseCoupon*  revCoupon=NULL;
	ARM_ExerciseStyle*  exeSched=NULL;
	ARM_ReferenceValue* redemp=NULL;
	ARM_ReferenceValue* classredemp=NULL ;
 
	ARM_Reverse* reverse=NULL;
	ARM_Reverse* newReverse=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	char* sDualFlag = NULL;

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(dualDate,sDate);
		sDualFlag = dualFlag.GetStr();

		swapLegClass = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(classSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegClass) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegClass is not a swapleg");
			return ARM_KO;
		}

		swapLegStruct = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(structSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegStruct) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegStruct is not a swapleg");
			return ARM_KO;
		}

		revCoupon = (ARM_ReverseCoupon*) LOCAL_PERSISTENT_OBJECTS->GetObject(couponId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(revCoupon, ARM_REVERSECOUPON) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: revCoupon is not of a good type");
			return ARM_KO;
		}

		exeSched = (ARM_ExerciseStyle*) LOCAL_PERSISTENT_OBJECTS->GetObject(exeId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(exeSched, ARM_EXERCISE_STYLE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: exeSched is not of a good type");
			return ARM_KO;
		}

		redemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(redempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(redemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: redemp is not of a good type");
			return ARM_KO;
		}

		classredemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(classRedempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(classredemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: classredemp is not of a good type");
			return ARM_KO;
		}

		newReverse = new ARM_Reverse(swapLegStruct, swapLegClass, ReceiveOrPay,
									 revCoupon,exeSched,redemp, classredemp,
									 (ARM_Date)sDate , dualStrike, sDualFlag);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		if (newReverse == NULL)
		{
			result.setMsg ("ARM_ERR: Reverse is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			revId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse);

			if (revId == RET_KO)
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(revId);

			return ARM_OK;
		}
		else
		{
			reverse = (ARM_Reverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(reverse, ARM_REVERSE) == 1)
			{
				if (reverse)
				{
					delete reverse;
					reverse = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse, objId);

				return ARM_OK;
			}

			else
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newReverse)
			delete newReverse;
		newReverse = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_STRUCTREVERSECOUPON(const VECTOR<double>& date,
								  const VECTOR<double>& strike,
								  const VECTOR<double>& power,
								  const VECTOR<double>& callput,
								  double xo,
								  const VECTOR<double>& floor,
								  const VECTOR<double>& cap,
								  ARM_result& result,
								  long objId)
{
	long revcpnId;

    ARM_ReverseCoupon* reverseCp=NULL;
    ARM_ReverseCoupon* newReverseCp=NULL;

    ARM_Vector* datevec = NULL;
    ARM_Vector* strikevec = NULL;
    ARM_Vector* powervec = NULL;
    ARM_Vector* callputflagvec = NULL;
	ARM_Vector* floorvec = NULL;
	ARM_Vector* capvec = NULL;

	long size = date.size();

	if (!((size == strike.size ()) 
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

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	VECTOR<CCString> date_str;
	char* sDate = new char[11];

	CCString msg ("");

	try
	{
		datevec = new ARM_Vector(size);
		strikevec = new ARM_Vector(size);
		powervec = new ARM_Vector(size);
		callputflagvec = new ARM_Vector(size);
		floorvec = new ARM_Vector(size);
		capvec = new ARM_Vector(size);

		for (int i = 0; i < size; i++)
		{
			Local_XLDATE2ARMDATE (date[i],sDate);
			datevec->Elt(i) = ((ARM_Date)sDate).GetJulian();
			strikevec->Elt(i) = strike[i];
			powervec->Elt(i) = power[i];
			callputflagvec->Elt(i) = callput[i];
			floorvec->Elt(i) = floor[i];
			capvec->Elt(i) = cap[i];
		}

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		newReverseCp = new ARM_ReverseCoupon(datevec, strikevec,
												 powervec, 
												 callputflagvec, xo,floorvec,capvec);

		if (datevec)
			delete datevec;
		datevec = NULL;

		if (strikevec)
			delete strikevec;
		strikevec = NULL;

		if (powervec)
			delete powervec;
		powervec = NULL;

		if (callputflagvec)
			delete callputflagvec;
		callputflagvec = NULL;

		if (floorvec)
			delete floorvec;
		floorvec = NULL;

		if (capvec)
			delete capvec;
		capvec = NULL;

		if (newReverseCp == NULL)
		{
			result.setMsg ("ARM_ERR: ReverseCoupon is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			revcpnId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverseCp);

			if (revcpnId == RET_KO)
			{
				if (newReverseCp)
					delete newReverseCp;
				newReverseCp = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(revcpnId);

			return ARM_OK;
		}
		else
		{
			reverseCp = (ARM_ReverseCoupon *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(reverseCp, ARM_REVERSECOUPON) == 1)
			{
				if (reverseCp)
				{
					delete reverseCp;
					reverseCp = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverseCp, objId);

				return ARM_OK;
			}

			else
			{
				if (newReverseCp)
					delete newReverseCp;
				newReverseCp = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newReverseCp)
			delete newReverseCp;
		newReverseCp = NULL;

		if (datevec)
			delete datevec;
		datevec = NULL;

		if (strikevec)
			delete strikevec;
		strikevec = NULL;

		if (powervec)
			delete powervec;
		powervec = NULL;

		if (callputflagvec)
			delete callputflagvec;
		callputflagvec = NULL;

		if (floorvec)
			delete floorvec;
		floorvec = NULL;

		if (capvec)
			delete capvec;
		capvec = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_REVERSE_CALENDAR (long structSwapLegId,
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
	long revId;

	int dSize = dStartDates.size();
	int fSize = fStartDates.size();

	if (
		(!((dSize == dEndDates.size ()) 
         && (dEndDates.size () == dFixingDates.size ()) 
         && (dFixingDates.size () == dPaymentDates.size ())
         ))
        ||
		(!((fSize == fEndDates.size ()) 
         && (fEndDates.size () == fFixingDates.size ()) 
         && (fFixingDates.size () == fPaymentDates.size ())
         ))
       )	
	{
		result.setMsg ("ARM_ERR: startDate, endDate, fixingDate and paymentDates arrays must have the same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	char* sDualFlag = NULL;

	int i ;
	ARM_SwapLeg*        swapLegClass=NULL;
	ARM_SwapLeg*        swapLegStruct=NULL;
	ARM_ReverseCoupon*  revCoupon=NULL;
	ARM_ExerciseStyle*  exeSched=NULL;
	ARM_ReferenceValue* redemp=NULL;
	ARM_ReferenceValue* classredemp=NULL;

	ARM_Vector * VdStartDates=NULL;
	ARM_Vector * VdEndDates=NULL;
	ARM_Vector * VdFixingDates=NULL;
	ARM_Vector * VdPaymentDates=NULL;
	ARM_Vector * VfStartDates=NULL;
	ARM_Vector * VfEndDates=NULL;
	ARM_Vector * VfFixingDates=NULL;
	ARM_Vector * VfPaymentDates=NULL;

	ARM_Reverse* reverse=NULL;
	ARM_Reverse* newReverse=NULL;

	CCString msg ("");

	try
	{
		sDualFlag = dualFlag.GetStr();
		
		swapLegClass = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(classSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegClass) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegClass is not a swapleg");
			return ARM_KO;
		}

		swapLegStruct = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(structSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegClass) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegStruct is not a swapleg");
			return ARM_KO;
		}

		revCoupon = (ARM_ReverseCoupon*) LOCAL_PERSISTENT_OBJECTS->GetObject(couponId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(revCoupon, ARM_REVERSECOUPON) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: revCoupon is not of a good type");
			return ARM_KO;
		}

		exeSched = (ARM_ExerciseStyle*) LOCAL_PERSISTENT_OBJECTS->GetObject(exeId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(exeSched, ARM_EXERCISE_STYLE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: exeSched is not of a good type");
			return ARM_KO;
		}

		redemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(redempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(redemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: redemp is not of a good type");
			return ARM_KO;
		}

		classredemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(classRedempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(classredemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: classredemp is not of a good type");
			return ARM_KO;
		}

		VdStartDates = new ARM_Vector (dSize , 0.0);
		VdEndDates = new ARM_Vector (dSize , 0.0);
		VdFixingDates = new ARM_Vector (dSize , 0.0);
		VdPaymentDates = new ARM_Vector (dSize , 0.0);
		VfStartDates = new ARM_Vector (fSize , 0.0);
		VfEndDates = new ARM_Vector (fSize , 0.0);
		VfFixingDates = new ARM_Vector (fSize , 0.0);
		VfPaymentDates = new ARM_Vector (fSize , 0.0);
		
		for (i= 0 ; i<dSize; i++)
		{
			Local_XLDATE2ARMDATE (dStartDates[i],sDate);
			VdStartDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dEndDates[i],sDate);
			VdEndDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dFixingDates[i],sDate);
			VdFixingDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dPaymentDates[i],sDate);
			VdPaymentDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
		}
		
		for (i= 0 ; i<fSize; i++)
		{
			Local_XLDATE2ARMDATE (fStartDates[i],sDate);
			VfStartDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fEndDates[i],sDate);
			VfEndDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fFixingDates[i],sDate);
			VfFixingDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fPaymentDates[i],sDate);
			VfPaymentDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
		}

		Local_XLDATE2ARMDATE (dualDate,sDate);

		newReverse = new ARM_Reverse( swapLegStruct, swapLegClass, 
									 ReceiveOrPay,revCoupon,exeSched, 
									 redemp, classredemp, 
									 (ARM_Date)sDate ,dualStrike,
									 sDualFlag,
									 VdStartDates,VdEndDates,VdFixingDates,VdPaymentDates,
									 VfStartDates,VfEndDates,VfFixingDates,VfPaymentDates);

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (VdStartDates)
			delete VdStartDates;
		VdStartDates = NULL;

		if (VdEndDates)
			delete VdEndDates;
		VdEndDates = NULL;

		if (VdFixingDates)
			delete VdFixingDates;
		VdFixingDates = NULL;

		if (VdPaymentDates)
			delete VdPaymentDates;
		VdPaymentDates = NULL;

		if (VfStartDates)
			delete VfStartDates;
		VfStartDates = NULL;

		if (VfEndDates)
			delete VfEndDates;
		VfEndDates = NULL;

		if (VfFixingDates)
			delete VfFixingDates;
		VfFixingDates = NULL;

		if (VfPaymentDates)
			delete VfPaymentDates;
		VfPaymentDates = NULL;

		if (newReverse == NULL)
		{
			result.setMsg ("ARM_ERR: Reverse Calendar is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			revId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse);

			if (revId == RET_KO)
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(revId);

			return ARM_OK;
		}
		else
		{
			reverse = (ARM_Reverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(reverse, ARM_REVERSE) == 1)
			{
				if (reverse)
				{
					delete reverse;
					reverse = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse, objId);

				return ARM_OK;
			}

			else
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newReverse)
			delete newReverse;
		newReverse = NULL;

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (VdStartDates)
			delete VdStartDates;
		VdStartDates = NULL;

		if (VdEndDates)
			delete VdEndDates;
		VdEndDates = NULL;

		if (VdFixingDates)
			delete VdFixingDates;
		VdFixingDates = NULL;

		if (VdPaymentDates)
			delete VdPaymentDates;
		VdPaymentDates = NULL;

		if (VfStartDates)
			delete VfStartDates;
		VfStartDates = NULL;

		if (VfEndDates)
			delete VfEndDates;
		VfEndDates = NULL;

		if (VfFixingDates)
			delete VfFixingDates;
		VfFixingDates = NULL;

		if (VfPaymentDates)
			delete VfPaymentDates;
		VfPaymentDates = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_REVERSENOTIONAL_CALENDAR (long structSwapLegId,
										long classSwapLegId,
										long ReceiveOrPay,
										long couponId,
										long exeId,
										long redempId,
										long classRedempId,
										long notExchFlag,
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
	long revId;

	int dSize = dStartDates.size();
	int fSize = fStartDates.size();

	if (
		(!((dSize == dEndDates.size ()) 
         && (dEndDates.size () == dFixingDates.size ()) 
         && (dFixingDates.size () == dPaymentDates.size ())
         ))
        ||
		(!((fSize == fEndDates.size ()) 
         && (fEndDates.size () == fFixingDates.size ()) 
         && (fFixingDates.size () == fPaymentDates.size ())
         ))
       )	
	{
		result.setMsg ("ARM_ERR: startDate, endDate, fixingDate and paymentDates arrays must have the same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	char* sDualFlag = NULL;

	int i ;
	ARM_SwapLeg*        swapLegClass=NULL;
	ARM_SwapLeg*        swapLegStruct=NULL;
	ARM_ReverseCoupon*  revCoupon=NULL;
	ARM_ExerciseStyle*  exeSched=NULL;
	ARM_ReferenceValue* redemp=NULL;
	ARM_ReferenceValue* classredemp=NULL;

	ARM_Vector * VdStartDates=NULL;
	ARM_Vector * VdEndDates=NULL;
	ARM_Vector * VdFixingDates=NULL;
	ARM_Vector * VdPaymentDates=NULL;
	ARM_Vector * VfStartDates=NULL;
	ARM_Vector * VfEndDates=NULL;
	ARM_Vector * VfFixingDates=NULL;
	ARM_Vector * VfPaymentDates=NULL;

	ARM_Reverse* reverse=NULL;
	ARM_Reverse* newReverse=NULL;

	CCString msg ("");

	try
	{
		sDualFlag = dualFlag.GetStr();
		
		swapLegClass = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(classSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegClass) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegClass is not a swapleg");
			return ARM_KO;
		}

		swapLegStruct = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(structSwapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLegClass) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: swapLegStruct is not a swapleg");
			return ARM_KO;
		}

		revCoupon = (ARM_ReverseCoupon*) LOCAL_PERSISTENT_OBJECTS->GetObject(couponId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(revCoupon, ARM_REVERSECOUPON) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: revCoupon is not of a good type");
			return ARM_KO;
		}

		exeSched = (ARM_ExerciseStyle*) LOCAL_PERSISTENT_OBJECTS->GetObject(exeId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(exeSched, ARM_EXERCISE_STYLE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: exeSched is not of a good type");
			return ARM_KO;
		}

		redemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(redempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(redemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: redemp is not of a good type");
			return ARM_KO;
		}

		classredemp = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(classRedempId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(classredemp, ARM_REFERENCE_VALUE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (sDualFlag)
				delete sDualFlag;
			sDualFlag = NULL;

			result.setMsg ("ARM_ERR: classredemp is not of a good type");
			return ARM_KO;
		}

		VdStartDates = new ARM_Vector (dSize , 0.0);
		VdEndDates = new ARM_Vector (dSize , 0.0);
		VdFixingDates = new ARM_Vector (dSize , 0.0);
		VdPaymentDates = new ARM_Vector (dSize , 0.0);
		VfStartDates = new ARM_Vector (fSize , 0.0);
		VfEndDates = new ARM_Vector (fSize , 0.0);
		VfFixingDates = new ARM_Vector (fSize , 0.0);
		VfPaymentDates = new ARM_Vector (fSize , 0.0);
		
		for (i= 0 ; i<dSize; i++)
		{
			Local_XLDATE2ARMDATE (dStartDates[i],sDate);
			VdStartDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dEndDates[i],sDate);
			VdEndDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dFixingDates[i],sDate);
			VdFixingDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dPaymentDates[i],sDate);
			VdPaymentDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
		}
		
		for (i= 0 ; i<fSize; i++)
		{
			Local_XLDATE2ARMDATE (fStartDates[i],sDate);
			VfStartDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fEndDates[i],sDate);
			VfEndDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fFixingDates[i],sDate);
			VfFixingDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (fPaymentDates[i],sDate);
			VfPaymentDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
		}

		Local_XLDATE2ARMDATE (dualDate,sDate);

		newReverse = new ARM_Reverse( swapLegStruct, swapLegClass, 
									 ReceiveOrPay,revCoupon,exeSched, 
									 redemp, classredemp, notExchFlag,
									 (ARM_Date)sDate ,dualStrike,
									 sDualFlag,
									 VdStartDates,VdEndDates,VdFixingDates,VdPaymentDates,
									 VfStartDates,VfEndDates,VfFixingDates,VfPaymentDates);

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (VdStartDates)
			delete VdStartDates;
		VdStartDates = NULL;

		if (VdEndDates)
			delete VdEndDates;
		VdEndDates = NULL;

		if (VdFixingDates)
			delete VdFixingDates;
		VdFixingDates = NULL;

		if (VdPaymentDates)
			delete VdPaymentDates;
		VdPaymentDates = NULL;

		if (VfStartDates)
			delete VfStartDates;
		VfStartDates = NULL;

		if (VfEndDates)
			delete VfEndDates;
		VfEndDates = NULL;

		if (VfFixingDates)
			delete VfFixingDates;
		VfFixingDates = NULL;

		if (VfPaymentDates)
			delete VfPaymentDates;
		VfPaymentDates = NULL;

		if (newReverse == NULL)
		{
			result.setMsg ("ARM_ERR: Reverse Calendar is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			revId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse);

			if (revId == RET_KO)
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(revId);

			return ARM_OK;
		}
		else
		{
			reverse = (ARM_Reverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(reverse, ARM_REVERSE) == 1)
			{
				if (reverse)
				{
					delete reverse;
					reverse = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newReverse, objId);

				return ARM_OK;
			}

			else
			{
				if (newReverse)
					delete newReverse;
				newReverse = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newReverse)
			delete newReverse;
		newReverse = NULL;

		if (sDualFlag)
			delete sDualFlag;
		sDualFlag = NULL;

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (VdStartDates)
			delete VdStartDates;
		VdStartDates = NULL;

		if (VdEndDates)
			delete VdEndDates;
		VdEndDates = NULL;

		if (VdFixingDates)
			delete VdFixingDates;
		VdFixingDates = NULL;

		if (VdPaymentDates)
			delete VdPaymentDates;
		VdPaymentDates = NULL;

		if (VfStartDates)
			delete VfStartDates;
		VfStartDates = NULL;

		if (VfEndDates)
			delete VfEndDates;
		VfEndDates = NULL;

		if (VfFixingDates)
			delete VfFixingDates;
		VfFixingDates = NULL;

		if (VfPaymentDates)
			delete VfPaymentDates;
		VfPaymentDates = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_POWERREVERSE (long initPeriodLegId,
							long realFundLegId,
							long fxUnderLegId,
							long fxNumLegId,
							const VECTOR<double>& dNoticeDates,
							const VECTOR<double>& dCancelDates,
							long FX0type,
							double FX0,
							double fxStep,
							long capValueType,
							double capValue,
							double capStepValue,
							long floorValueType,
							double floorValue,
							double floorStepValue,
							long dualOptionFlag,
							double dualOptionStrike,
							double redempNoticeDate,
							double lastLiborFixing,
							double lastFxSpotFixing,
                            ARM_result& result,
							long objId)
{
	long powRevId;

	int dSize = dNoticeDates.size();
	if (dSize != dCancelDates.size())
	{
		result.setMsg ("ARM_ERR: cancellation and notice dates must have the same size");
		return ARM_KO;
	}


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_SwapLeg*        initPeriodLeg=NULL;
	ARM_SwapLeg*        realFundLeg=NULL;
	ARM_SwapLeg*        fxUnderLeg=NULL;
	ARM_SwapLeg*        fxNumLeg=NULL;

	ARM_Vector * vNoticeDates=NULL;
	ARM_Vector * vCancelDates=NULL;

	ARM_ReferenceValue* refFX0=NULL;
	ARM_ReferenceValue* tmpRefFX0=NULL;
	ARM_ReferenceValue* refCapValue=NULL;
	ARM_ReferenceValue* tmpRefCapValue=NULL;
	ARM_ReferenceValue* refFloorValue=NULL;
	ARM_ReferenceValue* tmpRefFloorValue=NULL;

	ARM_PowerReverse* powerReverse=NULL;
	ARM_PowerReverse* newPowerReverse=NULL;

	char sDate[11];

	CCString msg ("");

	try
	{
		if (initPeriodLegId != ARM_NULL_OBJECT)
		{
			initPeriodLeg = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(initPeriodLegId);

			if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(initPeriodLeg) == 0)
			{
				result.setMsg ("ARM_ERR: initPeriodLeg is not a swapleg");
				return ARM_KO;
			}
		}

		realFundLeg = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(realFundLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(realFundLeg) == 0)
		{
			result.setMsg ("ARM_ERR: realFundLeg is not a swapleg");
			return ARM_KO;
		}

		fxUnderLeg = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(fxUnderLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(fxUnderLeg) == 0)
		{
			result.setMsg ("ARM_ERR: fxUnderLeg is not a swapleg");
			return ARM_KO;
		}

		fxNumLeg = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(fxNumLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(fxNumLeg) == 0)
		{
			result.setMsg ("ARM_ERR: fxNumLeg is not a swapleg");
			return ARM_KO;
		}

		vNoticeDates = new ARM_Vector (dSize , 0.0);
		vCancelDates = new ARM_Vector (dSize , 0.0);

		for (int i= 0 ; i<dSize; i++)
		{
			Local_XLDATE2ARMDATE (dNoticeDates[i],sDate);
			vNoticeDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
			Local_XLDATE2ARMDATE (dCancelDates[i],sDate);
			vCancelDates->Elt(i) = ((ARM_Date)sDate).GetJulian();
		}

        char noticeDate[11];

        if ( redempNoticeDate < 0 )
        {
           strcpy(noticeDate, "01/01/1981");
        }
        else
        {
            Local_XLDATE2ARMDATE(redempNoticeDate, noticeDate);
        }

		if (FX0type == 0)
		{
			ARM_Vector* vdates = fxUnderLeg->GetResetDates();

			refFX0 = GenerateStepRefValue(vdates,
										  FX0,
										  fxStep);
			tmpRefFX0 = refFX0;
		}
		else
		{
			refFX0 = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(FX0));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refFX0, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: FX0 is not of a good type");
				return ARM_KO;
			}
		}

		if (capValueType == 0)
		{
			ARM_Vector* vdates = fxUnderLeg->GetResetDates();

			refCapValue = GenerateStepRefValue(vdates,
											   capValue,
											   capStepValue);
			tmpRefCapValue = refCapValue;
		}
		else
		{
			refCapValue = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(capValue));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refCapValue, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: cap value is not of a good type");
				return ARM_KO;
			}
		}

		if (floorValueType == 0)
		{
			ARM_Vector* vdates = fxUnderLeg->GetResetDates();

			refFloorValue = GenerateStepRefValue(vdates,
												 floorValue,
												 floorStepValue);
			tmpRefFloorValue = refFloorValue;
		}
		else
		{
			refFloorValue = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(floorValue));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refFloorValue, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: floor value is not of a good type");
				return ARM_KO;
			}
		}


		newPowerReverse = new ARM_PowerReverse(initPeriodLeg,
											   realFundLeg,
											   fxUnderLeg,
											   fxNumLeg,
											   vNoticeDates,
											   vCancelDates,
											   refFX0,
											   refCapValue,
											   refFloorValue,
											   dualOptionFlag,
											   dualOptionStrike,
											   (ARM_Date) noticeDate,
											   lastLiborFixing,
											   lastFxSpotFixing);

		if (vNoticeDates)
			delete vNoticeDates;
		vNoticeDates = NULL;

		if (vCancelDates)
			delete vCancelDates;
		vCancelDates = NULL;

		if (tmpRefFX0)
			delete tmpRefFX0;
		tmpRefFX0 = NULL;

		if (tmpRefCapValue)
			delete tmpRefCapValue;
		tmpRefCapValue = NULL;

		if (tmpRefFloorValue)
			delete tmpRefFloorValue;
		tmpRefFloorValue = NULL;

		if (newPowerReverse == NULL)
		{
			result.setMsg ("ARM_ERR: Power Reverse is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			powRevId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPowerReverse);

			if (powRevId == RET_KO)
			{
				if (newPowerReverse)
					delete newPowerReverse;
				newPowerReverse = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(powRevId);

			return ARM_OK;
		}
		else
		{
			powerReverse = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(powerReverse, ARM_POWERREVERSE) == 1)
			{
				if (powerReverse)
				{
					delete powerReverse;
					powerReverse = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPowerReverse, objId);

				return ARM_OK;
			}

			else
			{
				if (newPowerReverse)
					delete newPowerReverse;
				newPowerReverse = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newPowerReverse)
			delete newPowerReverse;
		newPowerReverse = NULL;

		if (vNoticeDates)
			delete vNoticeDates;
		vNoticeDates = NULL;

		if (vCancelDates)
			delete vCancelDates;
		vCancelDates = NULL;

		if (tmpRefFX0)
			delete tmpRefFX0;
		tmpRefFX0 = NULL;

		if (tmpRefCapValue)
			delete tmpRefCapValue;
		tmpRefCapValue = NULL;

		if (tmpRefFloorValue)
			delete tmpRefFloorValue;
		tmpRefFloorValue = NULL;

		ARM_RESULT();
	}
}
/*
 * Function to give back data from a PowerReverse
 */
extern long ARMLOCAL_DatePowerReverseGetData(
	long prcsId,
	long dataType,
	VECTOR<double>& Data,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_PowerReverse* powerreverse= NULL;
		if( !GetObjectFromId( &powerreverse, prcsId, ARM_POWERREVERSE ) )
		{
			result.setMsg ("ARM_ERR: power reverse is not of a good type");
			return ARM_KO;
		};
		
		ARM_Vector* tmpData = powerreverse->GetMemberDateData( dataType );

		int i;
		/// test whether we have dates in which case, we convert this!
		if(		dataType == K_CANCEL_DATES
			||	dataType == K_NOTIF_DATES )
		{
			for( i=0; i<tmpData->size(); ++i)
			{
				Data.push_back( JulianToXLDate((*tmpData)[i]) );
			}
		}
		/// otherwise no treatment!
		else
		{
			for( i=0; i<tmpData->size(); ++i)
			{
				Data.push_back( (*tmpData)[i] );
			}
		}

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




long ARMLOCAL_GETPRCSDATA (long prcsId,
						   ARM_result& result)
{
	ARM_PowerReverse* prcs = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
    {
		prcs = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: PRCS is not of a good type");
			return ARM_KO;
		}

		prcs->GetPRCSData((CCString)"123" + (CCString)(CC_NS(ARM,ARM_USERNAME).c_str()));
    }
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }

	return ARM_OK;
}



long ARMLOCAL_DELETENEXTCALLFROMPRCS (long prcsId,
									  double asOf,
									  ARM_result& result,
									  long objId)
{
	long powRevId;

	ARM_PowerReverse* prcs = NULL;
	ARM_PowerReverse* newPrcs = NULL;
	ARM_PowerReverse* oldPrcs = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sAsOf[11];

	CCString msg ("");

    try
	{
		Local_XLDATE2ARMDATE(asOf,sAsOf);
		ARM_Date myDate(sAsOf);

		prcs = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: PRCS is not of a good type");
			return ARM_KO;
		}

		newPrcs = (ARM_PowerReverse *) prcs->Clone();

		if (newPrcs == NULL)
		{
			result.setMsg ("ARM_ERR: Power Reverse is null");
			return ARM_KO;
		}

		newPrcs->DeleteNextCall(myDate);
	
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			powRevId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPrcs);

			if (powRevId == RET_KO)
			{
				if (newPrcs)
					delete newPrcs;
				newPrcs = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(powRevId);

			return ARM_OK;
		}
		else
		{
			oldPrcs = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldPrcs, ARM_POWERREVERSE) == 1)
			{
				if (oldPrcs)
				{
					delete oldPrcs;
					oldPrcs = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPrcs, objId);

				return ARM_OK;
			}

			else
			{
				if (newPrcs)
					delete newPrcs;
				newPrcs = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
		return ARM_OK;

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		if (newPrcs)
			delete newPrcs;
		newPrcs = NULL;

		ARM_RESULT();
	}
}


//----------------------------------------------------
// PRCS access function
// if dateType = "C", get vector of Cancellation dates
// if dateType = "N", get vector of Notice dates
//----------------------------------------------------
long ARMLOCAL_GETOPTIONDATES(long prcsId,
							 const CCString& dateType,
							 ARM_result& result)
{
    CCString msg(""); // used in macro ARM_RESULT()
	ARM_PowerReverse* prcs = NULL;
	ARM_Vector* dates = NULL;
	ARM_Vector* oldVector = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

    try
	{
		prcs = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: PRCS is not of a good type");
			return ARM_KO;
		}

		if (strcmp((const char*) dateType,"C") == 0)
		{
			dates = prcs->GetItsCancelDates();
		}
		else if (strcmp((const char*) dateType,"N") == 0)
		{
			dates = prcs->GetItsNoticeDates();
		}
		else
		{
			result.setMsg ("ARM_ERR: Date type must be C (cancel) or N (notice)");
			return ARM_KO;
		}

		for (int i = 0; i < dates->size(); i++)
		{
			char dateStr[20];
			ARM_Date myDate;
			myDate.ChangeDate(dates->Elt(i));
			myDate.JulianToStrDate(dateStr);
			result.setStringInVect(dateStr);
		}
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}


//----------------------------------------------------
// PRCS access function
// if dateType = "C", get vector of Cancellation dates
// if dateType = "N", get vector of Notice dates
//----------------------------------------------------
long ARMLOCAL_GETDUALOPTIONSTRIKE(long prcsId,
								  ARM_result& result)
{
    CCString msg(""); // used in macro ARM_RESULT()
	ARM_PowerReverse* prcs = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

    try
    {
		prcs = (ARM_PowerReverse *) LOCAL_PERSISTENT_OBJECTS->GetObject(prcsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prcs, ARM_POWERREVERSE) == 0)
		{
			result.setMsg ("ARM_ERR: PRCS is not of a good type");
			return ARM_KO;
		}

		result.setDouble(prcs->GetItsDualOptionStrike());
 
		return ARM_OK;
	}
 
    catch(Exception& x)
    {
		x.DebugPrint();
 
		ARM_RESULT();
    }
}