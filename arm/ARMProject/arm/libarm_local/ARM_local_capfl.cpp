#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <ARM\libarm\ARM_result.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include "ARM_local_class.h"

#include <inst\capfloor.h>
#include <inst\flxcap.h>
#include <inst\maturitycapfl.h>
#include <inst\armsticky.h>
#include <inst\armratchet.h>
#include <inst\spreadoption.h>
#include <inst\corridorDblCondition.h>
#include <inst\armdigital.h>
#include <inst\dualcap.h>
#include <util\fromto.h>
#include <inst\globalcap.h>

#include <GP_Inflation\gpinflation\infcapfloor.h>
using ARM::ARM_InfCapFloor;


/// dummy virtual constructor
ARM_CapFloor* CorrespondingCapFloor( ARM_SwapLeg* leg, int capFloor, double strike )
{
	/// handles the null pointor
	if( leg == NULL )
		return new ARM_CapFloor( leg, capFloor, strike );

	switch( leg->GetName() )
	{
	case ARM_INFLEG:
		return new ARM_InfCapFloor( leg, capFloor, strike );
	default:
		return new ARM_CapFloor( leg, capFloor, strike );
	}
}

ARM_CLASS_NAME CorrespondingCapFloorClass( ARM_SwapLeg* leg )
{
	/// handles the null pointor
	if( leg == NULL )
		return ARM_CAPFLOOR;

	switch( leg->GetName() )
	{
	case ARM_INFLEG:
		return ARM_INFCAPFLOOR;
	default:
		return ARM_CAPFLOOR;
	}
}


//// constructor with a refValue
ARM_CapFloor* CorrespondingCapFloorwRef( ARM_SwapLeg* leg, int capFloor, ARM_ReferenceValue *strike )
{
	switch( leg->GetName() )
	{
	case ARM_INFLEG:
		return new ARM_InfCapFloor( leg, capFloor, strike );
	default:
		return new ARM_CapFloor( leg, capFloor, strike );
	}
}


/// Corresponding Cap and Floor with Id
CCString CorrespondingCapFloorClasswId( long swapLegId )
{
	ARM_SwapLeg* leg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

	if( leg == NULL )
		return LOCAL_CAPFLOOR_CLASS;
	else
	{
		switch( leg->GetName() )
		{
			case ARM_INFLEG:
				return LOCAL_INFCAPFLOOR_CLASS;
			default:
				return LOCAL_CAPFLOOR_CLASS;
		}
	}
}


long ARMLOCAL_LIBORCF (double startDate,
					   double endDate,
					   long isItCapOrFloor,
					   long strikeType,
					   double strike,
					   long liborType,
					   double spread,
					   long resetFreq,
					   long payFreq,
					   bool ccyIsObject,
					   const CCString& ccyName,
					   ARM_result& result,
					   long objId)
{
	long capflId;

    ARM_Currency* ccy=NULL;
    ARM_CapFloor* capFloor=NULL;  
    ARM_CapFloor* newLiborCapFloor=NULL;  
 
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
				result.setMsg ("ARM_ERR: ccy is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		if (strikeType == 0)
		{
			newLiborCapFloor = new ARM_CapFloor((ARM_Date) myStartDate,
												(ARM_Date) myEndDate,
												isItCapOrFloor,
												strike,
												(ARM_INDEX_TYPE) liborType,
												spread,
												resetFreq,
												payFreq,
												ccy);
		}
		else
		{
			ARM_ReferenceValue* refValue;

			refValue = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refValue, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: strike is not of a good type");
				return ARM_KO;
			}

			newLiborCapFloor = new ARM_CapFloor((ARM_Date) myStartDate,
												(ARM_Date) myEndDate,
												isItCapOrFloor,
												refValue,
												(ARM_INDEX_TYPE) liborType,
												spread,
												resetFreq,
												payFreq,
												ccy);
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (newLiborCapFloor == NULL)
		{
			result.setMsg ("ARM_ERR: Libor CapFloor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			capflId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLiborCapFloor);

			if (capflId == RET_KO)
			{
				if (newLiborCapFloor)
					delete newLiborCapFloor;
				newLiborCapFloor = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(capflId);

			return ARM_OK;
		}
		else
		{
			capFloor = (ARM_CapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(capFloor, ARM_CAPFLOOR) == 1)
			{
				if (capFloor)
				{
					delete capFloor;
					capFloor = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLiborCapFloor, objId);

				return ARM_OK;
			}

			else
			{
				if (newLiborCapFloor)
					delete newLiborCapFloor;
				newLiborCapFloor = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newLiborCapFloor)
			delete newLiborCapFloor;
		newLiborCapFloor = NULL;

		ARM_RESULT();
    } 
}



long ARMLOCAL_LIBORFLEXCF (double startDate,
						   double endDate,
						   long isItCapOrFloor,
						   double strike,
						   long nbEx,
						   long exerciseType,
						   long liborType,
						   double spread,
						   long resetFreq,
						   long payFreq,
						   bool ccyIsObject,
						   const CCString& ccyName,
						   ARM_result& result,
						   long objId)
{
	int capflId;

	ARM_Currency* ccy=NULL;
	ARM_FlexibleCapFloor* createdLiborFlexCapF=NULL;
	ARM_FlexibleCapFloor* prevLiborFlexCapF=NULL;

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
				result.setMsg ("ARM_ERR: ccy is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency ((const char*)ccyName);
		}

		createdLiborFlexCapF = new ARM_FlexibleCapFloor((ARM_Date) myStartDate,
												(ARM_Date) myEndDate,
												isItCapOrFloor, strike,
												nbEx, exerciseType,
												(ARM_INDEX_TYPE) liborType, 
												spread,
												resetFreq, payFreq,
												ccy);

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (createdLiborFlexCapF == NULL)
		{
			result.setMsg ("ARM_ERR: Libor Flexible CapFloor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			capflId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdLiborFlexCapF);

			if (capflId == RET_KO)
			{
				if (createdLiborFlexCapF)
					delete createdLiborFlexCapF;
				createdLiborFlexCapF = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(capflId);

			return ARM_OK;
		}
		else
		{
			prevLiborFlexCapF = (ARM_FlexibleCapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevLiborFlexCapF, ARM_FLEXIBLECAPFLOOR) == 1)
			{
				if (prevLiborFlexCapF)
				{
					delete prevLiborFlexCapF;
					prevLiborFlexCapF = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdLiborFlexCapF, objId);

				return ARM_OK;
			}

			else
			{
				if (createdLiborFlexCapF)
					delete createdLiborFlexCapF;
				createdLiborFlexCapF = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdLiborFlexCapF)
			delete createdLiborFlexCapF;
		createdLiborFlexCapF = NULL;

		ARM_RESULT();
    } 
}


long ARMLOCAL_CAPFLOOR (long swapLegId,
						long capOrFloor,
						long strikeType,
						double strike,
						ARM_result& result,
						long objId)
{
	int capflId;

	ARM_SwapLeg* swapLeg=NULL;
	ARM_CapFloor* capFloor=NULL;
	ARM_CapFloor* newCapFloor=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);
		ARM_CLASS_NAME CLASS_NAME = CorrespondingCapFloorClass( swapLeg );

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		if (strikeType == 0)
		{
			/// use the virtual function GetCorrespondingCap to get the appropriate cap
			/// for an inflation leg, creates the corresponding inflation cap
			newCapFloor = CorrespondingCapFloor( swapLeg, capOrFloor, strike);
		}
		else
		{
			ARM_ReferenceValue* refValue;

			refValue = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refValue, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: strike is not of a good type");
				return ARM_KO;
			}

			/// use the virtual function GetCorrespondingCap to get the appropriate cap
			/// for an inflation leg, creates the corresponding inflation cap
			newCapFloor = CorrespondingCapFloorwRef( swapLeg, capOrFloor, refValue);
		}

		if (newCapFloor == NULL)
		{
			result.setMsg ("ARM_ERR: CapFloor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			capflId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCapFloor);

			if (capflId == RET_KO)
			{
				if (newCapFloor)
					delete newCapFloor;
				newCapFloor = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(capflId);

			return ARM_OK;
		}
		else
		{
			capFloor = (ARM_CapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(capFloor, CLASS_NAME ) == 1)
			{
				if (capFloor)
				{
					delete capFloor;
					capFloor = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCapFloor, objId);

				return ARM_OK;
			}

			else
			{
				if (newCapFloor)
					delete newCapFloor;
				newCapFloor = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newCapFloor)
			delete newCapFloor;
		newCapFloor = NULL;

		ARM_RESULT();
    } 
}



long ARMLOCAL_FLEXCF (long swapLegId, long isItCapOrFloor, double strike,
						long nbEx, long exerciseType,
						ARM_result& result, long objId)
{
	long capflId;

	ARM_FlexibleCapFloor* newFxbleCapFloor=NULL;
	ARM_FlexibleCapFloor* fxbleCapFloor=NULL;
	ARM_SwapLeg* swapLeg=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

        newFxbleCapFloor = new ARM_FlexibleCapFloor(swapLeg, isItCapOrFloor,
                                   strike, nbEx, exerciseType);

		if (newFxbleCapFloor == NULL)
		{
			result.setMsg ("ARM_ERR: Flexible CapFloor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			capflId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFxbleCapFloor);

			if (capflId == RET_KO)
			{
				if (newFxbleCapFloor)
					delete newFxbleCapFloor;
				newFxbleCapFloor = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(capflId);

			return ARM_OK;
		}
		else
		{
			fxbleCapFloor = (ARM_FlexibleCapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(fxbleCapFloor, ARM_FLEXIBLECAPFLOOR) == 1)
			{
				if (fxbleCapFloor)
				{
					delete fxbleCapFloor;
					fxbleCapFloor = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newFxbleCapFloor, objId);

				return ARM_OK;
			}

			else
			{
				if (newFxbleCapFloor)
					delete newFxbleCapFloor;
				newFxbleCapFloor = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newFxbleCapFloor)
			delete newFxbleCapFloor;
		newFxbleCapFloor = NULL;

		ARM_RESULT();
    } 
}



long ARMLOCAL_MATCAPFLOOR (long swapLegId,
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
	long matcfId;

	ARM_SwapLeg*		swapLeg = NULL;
	ARM_MatCapFloor*	createdMatcf = NULL;
	ARM_MatCapFloor*	matcf = NULL;
	ARM_ReferenceValue* minStrikes = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		if (minStrikesId != ARM_NULL_OBJECT)
		{
			minStrikes = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(minStrikesId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(minStrikes,ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: min Strikes is not a Reference Value");
				return ARM_KO;
			}
		}

        if (maxMatDate < 0.0)
        {
            createdMatcf = new ARM_MatCapFloor(swapLeg,
										       annuity,
										       initNominal,
										       isTRI,
										       capOrFloor,
										       coeff,
										       firstTRIstrike,
										       minStrikes,
										       isDigitalPayoff,
										       increasingCoef);
        }
        else
        {
            char myDate[50];
            Local_XLDATE2ARMDATE(maxMatDate,myDate);
            ARM_Date MaxMatDate(myDate);
            createdMatcf = new ARM_MatCapFloor(swapLeg,MaxMatDate,
										       annuity,
										       initNominal,
										       isTRI,
										       capOrFloor,
										       coeff,
										       firstTRIstrike,
										       minStrikes,
										       isDigitalPayoff,
										       increasingCoef);
        }

		if (createdMatcf == NULL)
		{
			result.setMsg ("ARM_ERR: Maturity CapFloor is null");
			return ARM_KO;
		}
 
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			matcfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMatcf);

			if (matcfId == RET_KO)
			{
				if (createdMatcf)
					delete createdMatcf;
				createdMatcf = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(matcfId);

			return ARM_OK;
		}
		else
		{
			matcf = (ARM_MatCapFloor *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(matcf, ARM_MATCAPFLOOR) == 1)
			{
				if (matcf)
				{
					delete matcf;
					matcf = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdMatcf, objId);

				return ARM_OK;
			}

			else
			{
				if (createdMatcf)
					delete createdMatcf;
				createdMatcf = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdMatcf)
			delete createdMatcf;
		createdMatcf = NULL;

		ARM_RESULT();
    } 
}



long ARMLOCAL_STICKY (long swapLegId,
					  long capOrFloor,
					  double strike,
					  const VECTOR<double>& spreadDates,
					  const VECTOR<double>& spreadValues,
					  long kRefValId,
					  ARM_result& result,
					  long objId = -1)
{
	long stickyId;

	ARM_SwapLeg*  swapLeg = NULL;
	ARM_Sticky*   createdSticky = NULL;
	ARM_Sticky*   sticky = NULL;
	ARM_ReferenceValue* stickyStrikes = NULL;

	if(spreadDates.size () != spreadValues.size ())
	{
		result.setMsg ("ARM_ERR: spread date and value array must have same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	double dvect[ARM_NB_TERMS];
	double vvect[ARM_NB_TERMS];

	char strDate[11];

	ARM_Vector* vectSpreadDates = NULL;
	ARM_Vector* vectSpreadValues = NULL;
	ARM_ReferenceValue* spreads = NULL;
	CCString msg ("");

	try
	{

		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		for (int i = 0; i < spreadDates.size (); i++)
		{
			Local_XLDATE2ARMDATE (spreadDates[i],strDate);
			ARM_Date tmpDate(strDate);

			double dDate = tmpDate.GetJulian();

			dvect[i] = dDate;
			vvect[i] = spreadValues[i];
		}

		vectSpreadDates  = new ARM_Vector(spreadDates.size (), dvect);
		vectSpreadValues = new ARM_Vector(spreadDates.size (), vvect);
		spreads = new ARM_ReferenceValue(vectSpreadDates, vectSpreadValues);

		stickyStrikes =  (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(kRefValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stickyStrikes, ARM_REFERENCE_VALUE) == 0)
		{
			if (spreads)
				delete spreads;
			spreads = NULL;

			result.setMsg ("ARM_ERR: sticky Strikes is not of a good type");
			return ARM_KO;
		}

		createdSticky = new ARM_Sticky(swapLeg, spreads, stickyStrikes, capOrFloor,
									   strike);

		if (createdSticky == NULL)
		{
			result.setMsg ("ARM_ERR: Sticky is null");
			return ARM_KO;
		}
 
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			stickyId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSticky);

			if (stickyId == RET_KO)
			{
				if (createdSticky)
					delete createdSticky;
				createdSticky = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(stickyId);

			return ARM_OK;
		}
		else
		{
			sticky = (ARM_Sticky *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sticky, ARM_STICKY) == 1)
			{
				if (sticky)
				{
					delete sticky;
					sticky = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdSticky, objId);

				return ARM_OK;
			}

			else
			{
				if (createdSticky)
					delete createdSticky;
				createdSticky = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdSticky)
			delete createdSticky;
		createdSticky = NULL;

		if (vectSpreadDates)
			delete vectSpreadDates;
		vectSpreadDates = NULL;

		if (vectSpreadValues)
			delete vectSpreadValues;
		vectSpreadValues = NULL;

		if (spreads)
			delete spreads;
		spreads = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_SPREADOPTION (double startDate,
							double endDate,
							long capOrFloorId,
							long strike_type,
							double strike,
							long liborType1Id,
							long liborType2Id,
							double weight1,
							bool weight1IsReferenceValue,
							long weight1Id,
							double weight2,
							bool weight2IsReferenceValue,
							long weight2Id,
							long dayCountId,
							long resetFreqId,
							long payFreqId,
							long resetTimingId,
							long payTimingId,
							long ccyId,
                            long resetGap,
						    long fixing1_type,
						    long fixing1Id, 
						    VECTOR<double>& fixing1,
						    long fixing2_type,
						    long fixing2Id,
						    VECTOR<double>& fixing2,
							long intRule,
							long stubRule,
							int cptStrikeMethod,
							int computedFormula,
							VECTOR<double>& calibInfos,
							ARM_result& result,
							long objId)
{
	long sproId;
	
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	ARM_Currency *ccy = NULL;
	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;
	
	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());

	ARM_ReferenceValue* fixing2Ref = NULL;

	ARM_Vector calibInfosVect(calibInfos.size());
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	ARM_ReferenceValue* vWeight1Tmp = NULL;
	ARM_ReferenceValue* vWeight2Tmp = NULL;

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

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

        if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		ARM_ReferenceValue*	vWeight1 = NULL;
		if(weight1IsReferenceValue)
		{
			vWeight1 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(weight1Id);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vWeight1, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Weight1 is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			vWeight1Tmp = new ARM_ReferenceValue(weight1);
			vWeight1 = vWeight1Tmp;
		}

		ARM_ReferenceValue*	vWeight2 = NULL;
		if(weight2IsReferenceValue)
		{
			vWeight2 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(weight2Id);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vWeight2, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Weight2 is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			vWeight2Tmp = new ARM_ReferenceValue(weight2);
			vWeight2 = vWeight2Tmp;
		}

		if ( calibInfos.size() > 0 )
		{
			for (int i = 0; i < calibInfos.size(); i++)
			{
				calibInfosVect.Elt(i) = calibInfos[i];
			}
		}

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   vWeight1,
										   vWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   &fixing1Vect,
										   &fixing2Vect,
										   intRule,
										   stubRule,
										   cptStrikeMethod,
										   NULL,
										   NULL,
										   computedFormula,
                                           K_MOD_FOLLOWING,
										   "NULL",
										   10000, // payGap
										   &calibInfosVect);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   vWeight1,
										   vWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   fixing1Ref,
										   fixing2Ref,
										   intRule,
										   stubRule,
										   cptStrikeMethod,
										   NULL,
										   NULL,
										   computedFormula,
										   K_MOD_FOLLOWING,
										   "NULL",
										   10000, // payGap
										   &calibInfosVect);			
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");		
		}

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: SpreadOption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

		if(vWeight1Tmp)
		{
			delete	vWeight1Tmp;
		}

		if(vWeight2Tmp)
		{
			delete	vWeight2Tmp;
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if(vWeight1Tmp)
		{
			delete	vWeight1Tmp;
		}

		if(vWeight2Tmp)
		{
			delete	vWeight2Tmp;
		}

		ARM_RESULT();
    }
}

long ARMLOCAL_QUANTOSPREADOPTION (	double startDate,
									double endDate,
									long capOrFloor,
									long strike_type,
									double strikes,
									long liborIdx1Id,
									long liborIdx2Id,
									double Idx1weight,
									double Idx2weight,
									long Idx1fixingsId,
									long Idx2fixingsId,
									double idx1spread,
									double idx2spread,
									int cptStrikeMethod,
									int computedFormula,
									long CcyId,
									long dayCountId,
									long resetFreqId,
									long payFreqId,
									long resetTimingId,
									long payTimingId,
									long intRuleId,
									long stubRuleId,
									long resetGap,
									char* resetCal,
									char* payCal,
									long fwdRule,
									char* refDate,
									long notionalId,
									ARM_result& result,
									long objId)
{
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	long sproId;

	char sStartDate[11];
	char sEndDate[11];

	ARM_ReferenceValue* refStrikes = NULL;
	ARM_ReferenceValue* tmpRefStrikes = NULL;
	
	ARM_IRIndex* refIndex1 = NULL;
	ARM_IRIndex* tmpRefIndex1 = NULL;

	ARM_IRIndex* refIndex2 = NULL;
	ARM_IRIndex* tmpRefIndex2 = NULL;

	ARM_ReferenceValue* refIdx1fixings = NULL;
	ARM_ReferenceValue* refIdx2fixings = NULL;
	
	ARM_ReferenceValue* aWeight1 = NULL;
	ARM_ReferenceValue* aWeight2 = NULL;
	
	ARM_ReferenceValue* notional = NULL;

	ARM_Currency* ccy = NULL;
	ARM_Currency* ccy1 = NULL;
	ARM_Currency* ccy2 = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);
		
		if ( CcyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else if ( CcyId == ARM_FRF_CCY_OBJECT )
			ccy = ARM_FRF_CURRENCY;
		else
		{
			ccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);
			
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( liborIdx1Id < 0 )
		{
			refIndex1 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex1 = refIndex1;
		}
		else
		{
			refIndex1 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx1Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex1, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 is not of a good type");
				return ARM_KO;
			}
		}

		refIndex1->SetPayTiming(payTimingId);
		refIndex1->SetResetTiming(resetTimingId);
		refIndex1->SetResetFrequency(resetFreqId);
		refIndex1->SetPayFrequency(payFreqId);
		refIndex1->SetResetGap(resetGap);


		if ( liborIdx2Id < 0 )
		{
			refIndex2 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex2 = refIndex2;
		}
		else
		{
			refIndex2 = (ARM_IRIndex*) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex2, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index2 is not of a good type");
				return ARM_KO;
			}
		}

		refIndex2->SetResetFrequency(resetFreqId);
		refIndex2->SetPayFrequency(payFreqId);
		refIndex2->SetResetTiming(resetTimingId);
		refIndex2->SetPayTiming(payTimingId);
		refIndex2->SetResetGap(resetGap);


		if (strike_type)
		{
			refStrikes = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(strikes);
			
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrikes, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrikes = new ARM_ReferenceValue(strikes);
			tmpRefStrikes = refStrikes;
		}

		if( Idx1fixingsId >= 0 )
		{
			refIdx1fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx1fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx1fixings = NULL;
		}
		
		if( Idx2fixingsId >= 0)
		{
			refIdx2fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx2fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx2fixings = NULL;
		}

		if( notionalId >=0 )
		{
			notional = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(notional, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Notional is not of a good type");
				return ARM_KO;

			}
		}
		else
		{
			notional = NULL;
		}

		aWeight1 = new ARM_ReferenceValue(Idx1weight);
		aWeight2 = new ARM_ReferenceValue(Idx2weight);
		
		
		// vanilla spread option
		newSpro = new ARM_SpreadOption(	(ARM_Date) sStartDate, 
										(ARM_Date) sEndDate,
										capOrFloor, 
										refStrikes,
										refIndex1,  
										refIndex2,  
										aWeight1, 
										aWeight2,
										dayCountId, 
										resetFreqId, 
										payFreqId, 
										resetTimingId, 
										payTimingId, 
										ccy,
										resetGap,
										refIdx1fixings,
										refIdx2fixings,
										intRuleId,
										stubRuleId,
										cptStrikeMethod,
										resetCal,
										payCal,
										computedFormula,
										fwdRule,
										refDate,
										10000,
										NULL,
										notional);
		
	if (tmpRefIndex1)
		delete tmpRefIndex1;
	tmpRefIndex1 = NULL;

	if (tmpRefIndex2)
		delete tmpRefIndex2;
	tmpRefIndex2 = NULL;
	
	if (tmpRefStrikes)
		delete tmpRefStrikes;
	tmpRefStrikes = NULL;
	
	if (aWeight1)
		delete aWeight1;
	aWeight1 = NULL;

	if (aWeight2)
		delete aWeight2;
	aWeight2 = NULL;

	if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			
			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}
				
				sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				result.setLong(sproId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if (tmpRefIndex1)
			delete tmpRefIndex1;
		tmpRefIndex1 = NULL;

		if (tmpRefIndex2)
			delete tmpRefIndex2;
		tmpRefIndex2 = NULL;
		
		if (tmpRefStrikes)
			delete tmpRefStrikes;
		tmpRefStrikes = NULL;
		
		if (aWeight1)
			delete aWeight1;
		aWeight1 = NULL;

		if (aWeight2)
			delete aWeight2;
		aWeight2 = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
	
}


long ARMLOCAL_SPREADOPTIONWithLegs (long firstlegId,
									long secondlegId,
									long capOrFloorId,
									long strike_type,
									double strike,
									double weight1,
									bool weight1IsReferenceValue,
									long weight1Id,
									double weight2,
									bool weight2IsReferenceValue,
									long weight2Id,
									long fixing1_type,
									long fixing1Id, 
									VECTOR<double>& fixing1,
									long fixing2_type,
									long fixing2Id,
									VECTOR<double>& fixing2,
									int cptStrikeMethod,
									int computedFormula,
									ARM_result& result,
									long objId)
{
	long sproId;
	ARM_SwapLeg* firstLeg = NULL;
	ARM_SwapLeg* secondLeg = NULL;
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	ARM_Currency *ccy = NULL;
	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;
	
	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());

	ARM_ReferenceValue* fixing2Ref = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	ARM_ReferenceValue*	vWeight1Tmp = NULL;
	ARM_ReferenceValue*	vWeight2Tmp = NULL;

	try
	{
		firstLeg = (ARM_SwapLeg* )LOCAL_PERSISTENT_OBJECTS->GetObject(firstlegId);
		secondLeg = (ARM_SwapLeg* )LOCAL_PERSISTENT_OBJECTS->GetObject(secondlegId);

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

        if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		ARM_ReferenceValue*	vWeight1 = NULL;
		if(weight1IsReferenceValue)
		{
			vWeight1 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(weight1Id);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vWeight1, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Weight1 is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			vWeight1Tmp = new ARM_ReferenceValue(weight1);
			vWeight1 = vWeight1Tmp;
		}

		ARM_ReferenceValue*	vWeight2 = NULL;
		if(weight2IsReferenceValue)
		{
			vWeight2 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(weight2Id);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vWeight2, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Weight2 is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			vWeight2Tmp = new ARM_ReferenceValue(weight2);
			vWeight2 = vWeight2Tmp;
		}

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{
			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   vWeight1,
										   vWeight2,
										   &fixing1Vect,
										   &fixing2Vect,
										   cptStrikeMethod,
										   computedFormula);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   vWeight1,
										   vWeight2,
										   fixing1Ref,
										   fixing2Ref,
										   cptStrikeMethod,
										   computedFormula);			
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");		
		}

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: SpreadOption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

		if(vWeight1Tmp)
		{
			delete	vWeight1Tmp;
		}

		if(vWeight2Tmp)
		{
			delete	vWeight2Tmp;
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if(vWeight1Tmp)
		{
			delete	vWeight1Tmp;
		}

		if(vWeight2Tmp)
		{
			delete	vWeight2Tmp;
		}

		ARM_RESULT();
    }
}

long ARMLOCAL_SPREADDIGITAL(double startDate,
							double endDate,
							long capOrFloorId,
							long strike_type,
							double strike,
							long payoff_type,
							double payoff,
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
                            long resetGap,
							double spread1,
							double spread2,
						    long fixing1_type,
						    long fixing1Id, 
						    VECTOR<double>& fixing1,
						    long fixing2_type,
						    long fixing2Id,
						    VECTOR<double>& fixing2,
							long intRule,
							long stubRule,
							long slopeFlag,
							int cptStrikeMethod,
							int computedFormula,
							VECTOR<double>& calibInfos,
							ARM_result& result,
							long objId)
{
	long sproId;

	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	ARM_Currency *ccy = NULL;
	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;

	ARM_ReferenceValue *refPayoff = NULL;
	ARM_ReferenceValue *tmpRefPayoff = NULL;

	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());
	ARM_ReferenceValue* fixing2Ref = NULL;

	ARM_ReferenceValue aWeight1(weight1);
	ARM_ReferenceValue aWeight2(weight2);

	ARM_Vector calibInfosVect(calibInfos.size());

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

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

		if ( payoff_type == 1 )
		{
			refPayoff = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payoff));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refPayoff, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: PayOff is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refPayoff = new ARM_ReferenceValue(payoff);

			tmpRefPayoff = refPayoff;
		}

		if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		if ( calibInfos.size() > 0 )
		{
			for (int i = 0; i < calibInfos.size(); i++)
			{
				calibInfosVect.Elt(i) = calibInfos[i];
			}
		}

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
									   (ARM_Date)sEndDate,
									   capOrFloorId,
									   refStrike,
									   refPayoff,
									   (ARM_INDEX_TYPE)liborType1Id,
									   (ARM_INDEX_TYPE)liborType2Id,
									   &aWeight1,
									   &aWeight2,
									   dayCountId,
									   resetFreqId,
									   payFreqId,
									   resetTimingId,
									   payTimingId,
									   ccy,
                                       resetGap,
									   spread1,
									   spread2,
									   &fixing1Vect,
									   &fixing2Vect,
									   intRule,
									   stubRule,
									   NULL,
									   NULL,
									   slopeFlag,
									   cptStrikeMethod,
									   computedFormula,
									   K_MOD_FOLLOWING,
									   "NULL",
									   10000, // payGap
			calibInfosVect.GetSize() == 0 ? NULL : &calibInfosVect);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   refPayoff,
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   &aWeight1,
										   &aWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   spread1,
										   spread2,
										   fixing1Ref,
										   fixing2Ref,
										   intRule,
										   stubRule,
										   NULL,
										   NULL,
										   slopeFlag,
										   cptStrikeMethod,
										   computedFormula,
										   K_MOD_FOLLOWING,
										   "NULL",
										   10000, // payGap
				calibInfosVect.GetSize() == 0 ? NULL : &calibInfosVect);
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing type is not of a good type");		
		}

		if (tmpRefPayoff)
			delete tmpRefPayoff;
		tmpRefPayoff = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (tmpRefPayoff)
			delete tmpRefPayoff;
		tmpRefPayoff = NULL;

		ARM_RESULT();
    }
}

long ARMLOCAL_QUANTOSPREADDIGITAL ( double startDate,
									double endDate,
									long capOrFloor,
									long strike_type,
									double strikes,
									long payoff_type,
									double payoff,
									long liborIdx1Id,
									long liborIdx2Id,
									double Idx1weight,
									double Idx2weight,
									long Idx1fixingsId,
									long Idx2fixingsId,
									double spread1,
									double spread2,
									long slopeFlag,
									int cptStrikeMethod,
									int computedFormula,
									long CcyId,
									long dayCountId,
									long resetFreqId,
									long payFreqId,
									long resetTimingId,
									long payTimingId,
									long intRuleId,
									long stubRuleId,
									long resetGap,
									char* resetCal,
									char* payCal,
									long fwdRule,
									char* refDate,
									long notionalId,
									ARM_result& result,
									long objId)
{
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	long sproId;

	char sStartDate[11];
	char sEndDate[11];

	ARM_ReferenceValue* refStrikes = NULL;
	ARM_ReferenceValue* tmpRefStrikes = NULL;

	ARM_ReferenceValue* refPayoff = NULL;
	ARM_ReferenceValue* tmpRefPayoff = NULL;
	
	ARM_IRIndex* refIndex1 = NULL;
	ARM_IRIndex* tmpRefIndex1 = NULL;

	ARM_IRIndex* refIndex2 = NULL;
	ARM_IRIndex* tmpRefIndex2 = NULL;
	
	ARM_ReferenceValue* refIdx1fixings = NULL;
	ARM_ReferenceValue* refIdx2fixings = NULL;
	
	ARM_ReferenceValue* aWeight1 = NULL;
	ARM_ReferenceValue* aWeight2 = NULL;
	
	ARM_ReferenceValue* notional = NULL;

	ARM_Currency* ccy = NULL;
	ARM_Currency* ccy1 = NULL;
	ARM_Currency* ccy2 = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);
		
		if ( CcyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else if ( CcyId == ARM_FRF_CCY_OBJECT )
			ccy = ARM_FRF_CURRENCY;
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( liborIdx1Id < 0 )
		{
			refIndex1 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex1 = refIndex1;
		}
		else
		{
			refIndex1 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx1Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex1, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 is not of a good type");
				return ARM_KO;
			}
			refIndex1->SetPayTiming(payTimingId);
			refIndex1->SetResetTiming(resetTimingId);
			refIndex1->SetResetFrequency(resetFreqId);
			refIndex1->SetPayFrequency(payFreqId);
			refIndex1->SetResetGap(resetGap);
		}

		if ( liborIdx2Id < 0 )
		{
			refIndex2 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex2 = refIndex2;
		}
		else
		{
			refIndex2 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex2, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index2 is not of a good type");
				return ARM_KO;
			}
			refIndex2->SetResetFrequency(resetFreqId);
			refIndex2->SetPayFrequency(payFreqId);
			refIndex2->SetResetTiming(resetTimingId);
			refIndex2->SetPayTiming(payTimingId);
			refIndex2->SetResetGap(resetGap);
		}

		if (strike_type)
		{
			refStrikes = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(strikes);
		}
		else
		{
			refStrikes = new ARM_ReferenceValue(strikes);
			tmpRefStrikes = refStrikes;
		}
		
		if ( payoff_type)
		{
			refPayoff = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payoff));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refPayoff, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: PayOff is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refPayoff = new ARM_ReferenceValue(payoff);
			tmpRefPayoff = refPayoff;
		}

		if( Idx1fixingsId >= 0 )
		{
			refIdx1fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx1fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx1fixings = NULL;
		}
		
		if( Idx2fixingsId >= 0)
		{
			refIdx2fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx2fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx2fixings = NULL;
		}

		if( notionalId >=0 )
		{
			notional = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(notional, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Notional is not of a good type");
				return ARM_KO;

			}
		}
		else
		{
			notional = NULL;
		}

		aWeight1 = new ARM_ReferenceValue(Idx1weight);
		aWeight2 = new ARM_ReferenceValue(Idx2weight);
		
		// digital spread option
		newSpro = new ARM_SpreadOption(	(ARM_Date) sStartDate, 
										(ARM_Date) sEndDate,
										capOrFloor, 
										refStrikes,
										refPayoff,
										refIndex1,  
										refIndex2,  
										aWeight1, 
										aWeight2,
										dayCountId, 
										resetFreqId, 
										payFreqId, 
										resetTimingId, 
										payTimingId, 
										ccy,
										resetGap,
										spread1,
										spread2,
										refIdx1fixings,
										refIdx2fixings,
										intRuleId,
										stubRuleId,
										cptStrikeMethod,
										resetCal,
										payCal,
										slopeFlag,
										computedFormula,
										fwdRule,
										refDate,
										10000,
										NULL,
										notional);
		
	if (tmpRefIndex1)
		delete tmpRefIndex1;
	tmpRefIndex1 = NULL;

	if (tmpRefIndex2)
		delete tmpRefIndex2;
	tmpRefIndex2 = NULL;
	
	if (tmpRefStrikes)
		delete tmpRefStrikes;
	tmpRefStrikes = NULL;

	if (tmpRefPayoff)
		delete tmpRefPayoff;
	tmpRefPayoff = NULL;
	
	if (aWeight1)
		delete aWeight1;
	aWeight1 = NULL;

	if (aWeight2)
		delete aWeight2;
	aWeight2 = NULL;

	if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			
			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}
				
				sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				result.setLong(sproId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;
		
		if (tmpRefIndex1)
			delete tmpRefIndex1;
		tmpRefIndex1 = NULL;

		if (tmpRefIndex2)
			delete tmpRefIndex2;
		tmpRefIndex2 = NULL;
		
		if (tmpRefStrikes)
			delete tmpRefStrikes;
		tmpRefStrikes = NULL;

		if (tmpRefPayoff)
			delete tmpRefPayoff;
		tmpRefPayoff = NULL;
		
		if (aWeight1)
			delete aWeight1;
		aWeight1 = NULL;

		if (aWeight2)
			delete aWeight2;
		aWeight2 = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
	
}

long ARMLOCAL_SPREADDIGITALWithLegs(long firstlegId,
									long secondlegId,
									long capOrFloorId,
									long strike_type,
									double strike,
									long payoff_type,
									double payoff,
									double weight1,
									double weight2,
									double spread1,
									double spread2,
									long fixing1_type,
									long fixing1Id, 
									VECTOR<double>& fixing1,
									long fixing2_type,
									long fixing2Id,
									VECTOR<double>& fixing2,
									long slopeFlag,
									int cptStrikeMethod,
									int computedFormula,
									ARM_result& result,
									long objId)
{
	long sproId;

	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;

	ARM_SwapLeg* firstLeg = NULL;
	ARM_SwapLeg* secondLeg = NULL;
	
	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;

	ARM_ReferenceValue *refPayoff = NULL;
	ARM_ReferenceValue *tmpRefPayoff = NULL;

	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());
	ARM_ReferenceValue* fixing2Ref = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		firstLeg = (ARM_SwapLeg *)LOCAL_PERSISTENT_OBJECTS->GetObject(firstlegId);
		secondLeg = (ARM_SwapLeg *)LOCAL_PERSISTENT_OBJECTS->GetObject(secondlegId);

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

		if ( payoff_type == 1 )
		{
			refPayoff = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payoff));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refPayoff, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: PayOff is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refPayoff = new ARM_ReferenceValue(payoff);

			tmpRefPayoff = refPayoff;
		}

		if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		ARM_ReferenceValue vWeight1(weight1);
		ARM_ReferenceValue vWeight2(weight2);

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{
			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   refPayoff,
										   &vWeight1,
										   &vWeight2,
										   spread1,
										   spread2,
										   &fixing1Vect,
										   &fixing2Vect,
										   slopeFlag,
										   cptStrikeMethod,
										   computedFormula);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   refPayoff,
										   &vWeight1,
										   &vWeight2,
										   spread1,
										   spread2,
										   fixing1Ref,
										   fixing2Ref,
										   slopeFlag,
										   cptStrikeMethod,
										   computedFormula);
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing type is not of a good type");		
		}

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		ARM_RESULT();
    }	
}

long ARMLOCAL_SPREADCORRIDOR(double startDate,
							 double endDate,
							 long capOrFloorId,
							 long strike_type,
							 double strike,
							 long payIndexId,
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
							 long resetGap,
							 double spread1,
							 double spread2,
							 long fixing1_type,
							 long fixing1Id, 
							 VECTOR<double>& fixing1,
							 long fixing2_type,
							 long fixing2Id,
							 VECTOR<double>& fixing2,
							 long intRule,
							 long stubRule,
							 long slopeFlag,
							 int cptStrikeMethod,
							 long payMargin_type,
							 double payMargin,
							 double payWeight,
							 long fixing3_type,
							 long fixing3Id,
							 VECTOR<double>& fixing3,
							 int freezeFixing,
							 int computedFormula,
							 ARM_result& result,
							 long objId)
{
	long sproId;

	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	ARM_IRIndex* payIndex = NULL;
	ARM_Currency*      ccy = NULL;
	ARM_ReferenceValue* refStrike    = NULL;
	ARM_ReferenceValue* tmpRefStrike = NULL;

	ARM_ReferenceValue* refMargin    = NULL;
	ARM_ReferenceValue* tmpRefMargin = NULL;

	ARM_ReferenceValue* payFixedRate = NULL;
	ARM_ReferenceValue* tmpPayFixedRate = NULL;

	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing2.size());
	ARM_ReferenceValue* fixing2Ref = NULL;

	ARM_Vector fixing3Vect(fixing3.size());
	ARM_ReferenceValue* fixing3Ref = NULL;

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

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

		if ( payMargin_type == 1 )// objet
		{
			refMargin = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payMargin));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refMargin, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Margins is not of a good type");
				return ARM_KO;
			}
		}
		else//nombre
		{
			refMargin = new ARM_ReferenceValue(payMargin);

			tmpRefMargin = refMargin;
		}

		
		payFixedRate = new ARM_ReferenceValue(0.);
		tmpPayFixedRate = payFixedRate;

		if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		if ( fixing3_type == 1 )
		{
		   if (fixing3Id != ARM_NULL_OBJECT)
		   {
			   fixing3Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing3Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing3Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing3_type == 0)
		{
			if ( fixing3.size() > 1 )
			{
				for (int i = 0; i < fixing3.size(); i++)
				{
					fixing3Vect.Elt(i) = fixing3[i];
				}
			}
			else
			{
				fixing3Vect.Elt(0) = fixing3[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing3_type is not of a good type");
		}

		// PayIndex
		if(payIndexId ==  ARM_NULL_OBJECT)
		{
			payIndex = new ARM_IRIndex(LIBOR3M, payFreqId, payFreqId, ccy, dayCountId);
			payIndex->SetResetTiming(K_ADVANCE);
			payIndex->SetIntRule(K_ADJUSTED);
		}
		else
		{
			payIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(payIndex, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: PayIndex is not of a good type");
				return ARM_KO;
			}			
		}

		ARM_ReferenceValue aWeight1(weight1);
		ARM_ReferenceValue aWeight2(weight2);

		if(fixing1_type==fixing2_type && fixing1_type == 0 && fixing3_type==fixing1_type)
		{

			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   payIndex,//(ARM_INDEX_TYPE)payoffliborType1Id,
										   payFixedRate,
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   &aWeight1,
										   &aWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   spread1,
										   spread2,
										   &fixing1Vect,
										   &fixing2Vect,
										   intRule,
										   stubRule,
										   NULL,
										   NULL,
										   slopeFlag,
										   cptStrikeMethod,
										   refMargin,
										   payWeight,
										   &fixing3Vect,
										   freezeFixing,
										   computedFormula);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1 && fixing3_type==fixing1_type)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   payIndex,//(ARM_INDEX_TYPE)payoffliborType1Id, 
										   payFixedRate,
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   &aWeight1,
										   &aWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   spread1,
										   spread2,
										   fixing1Ref,
										   fixing2Ref,
										   intRule,
										   stubRule,
										   NULL,
										   NULL,
										   slopeFlag,
										   cptStrikeMethod,
										   refMargin,
										   payWeight,
										   fixing3Ref,
										   freezeFixing,
										   computedFormula);
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1, fixing2, fixing3 should be the same type");		
		}


		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (tmpPayFixedRate)
			delete tmpPayFixedRate;
		tmpPayFixedRate = NULL;

		if (tmpRefMargin)
			delete tmpRefMargin;
		tmpRefMargin = NULL;

		if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (tmpPayFixedRate)
			delete tmpPayFixedRate;
		tmpPayFixedRate = NULL;

		if (tmpRefMargin)
			delete tmpRefMargin;
		tmpRefMargin = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
}

//----------------------------------------------------------------------------------//
// Interface function for VMS corridor spread option								//
//----------------------------------------------------------------------------------//
long ARMLOCAL_SPREADCORRIDORVMS(double startDate,
							    double endDate,
							    long capOrFloorId,
							    long strike_type,
							    double strike,
							    long payIndexId,
							    long CMSIndexes1Id,
							    long CMSIndexes2Id,
							    double weight1,
							    double weight2,
							    long dayCountId,
							    long resetFreqId,
							    long payFreqId,
							    long resetTimingId,
							    long payTimingId,
							    long ccyId,
							    long resetGap,
							    double spread1,
							    double spread2,
							    long Fixing1Id, 
							    long Fixing2Id,
							    long intRule,
							    long stubRule,
							    long slopeFlag,
							    int cptStrikeMethod,
							    long payMargin_type,
							    double payMargin,
							    double payWeight,
							    long FixingPayId,
							    int freezeFixing,
							    int computedFormula,
							    ARM_result& result,
							    long objId)
{
	long sproId;

	// Parameters initialization
	//------------------------------------------------------

	ARM_SpreadOptionVMS* spro			= NULL;
	ARM_SpreadOptionVMS* newSpro		= NULL; 


	ARM_IRIndex*		payIndex		= NULL;
	ARM_Currency*		ccy				= NULL;

	ARM_ReferenceValue* refStrike		= NULL;
	ARM_ReferenceValue* tmpRefStrike	= NULL;

	ARM_ReferenceValue* refMargin		= NULL;
	ARM_ReferenceValue* tmpRefMargin	= NULL;

	ARM_ReferenceValue* payFixedRate	= NULL;
	ARM_ReferenceValue* tmpPayFixedRate	= NULL;

	ARM_ReferenceValue* CMSIndexes1		= NULL;
	ARM_ReferenceValue* CMSIndexes2		= NULL;

	ARM_ReferenceValue* Fixing1Ref		= NULL;
	ARM_ReferenceValue* Fixing2Ref		= NULL;
	ARM_ReferenceValue* FixingPayRef	= NULL;


	// First check
	//------------------------------------------------------

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
		// Construction of the start and end dates
		//------------------------------------------------------

		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);


		// Construction of the currency
		//------------------------------------------------------

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


		// Construction of the strike
		//------------------------------------------------------

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike	 = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}


		// Construction of the margin
		//------------------------------------------------------

		if ( payMargin_type == 1 )	// object
		{
			refMargin = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payMargin));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refMargin, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Margins is not of a good type");
				return ARM_KO;
			}
		}
		else						//number
		{
			refMargin = new ARM_ReferenceValue(payMargin);

			tmpRefMargin = refMargin;
		}

		
		payFixedRate    = new ARM_ReferenceValue(0.);
		tmpPayFixedRate = payFixedRate;


		// Construction of the first RefValue for CMS Indexes
		//------------------------------------------------------

	   if (CMSIndexes1Id != ARM_NULL_OBJECT)
	   {
		   CMSIndexes1 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(CMSIndexes1Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CMSIndexes1, ARM_REFERENCE_VALUE) == 0 )
		   {
			  result.setMsg ("ARM_ERR: Index types reference value is not of a good type");
				
			  return(ARM_KO);
			}
	   }


		// Construction of the second RefValue for CMS Indexes
		//------------------------------------------------------

	   if (CMSIndexes2Id != ARM_NULL_OBJECT)
	   {
		   CMSIndexes2 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(CMSIndexes2Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CMSIndexes2, ARM_REFERENCE_VALUE) == 0 )
		   {
			  result.setMsg ("ARM_ERR: Index types reference value is not of a good type");
				
			  return(ARM_KO);
			}
	   }


		// Construction of the first fixing
		//------------------------------------------------------

	   if (Fixing1Id != ARM_NULL_OBJECT)
	   {
		   Fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fixing1Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
		   {
			  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
				
			  return(ARM_KO);
			}
	   }


		// Construction of the second fixing
		//------------------------------------------------------

	   if (Fixing2Id != ARM_NULL_OBJECT)
	   {
		   Fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(Fixing2Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
		   {
			  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
				
			  return(ARM_KO);
			}
	   }


		// Construction of the third fixing
		//------------------------------------------------------

	   if (FixingPayId != ARM_NULL_OBJECT)
	   {
		   FixingPayRef = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(FixingPayId);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(FixingPayRef, ARM_REFERENCE_VALUE) == 0 )
		   {
			  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
				
			  return(ARM_KO);
			}
	   }


		// Construction of the payment index
		//------------------------------------------------------

		if(payIndexId ==  ARM_NULL_OBJECT)
		{
			payIndex = new ARM_IRIndex(LIBOR3M, payFreqId, payFreqId, ccy, dayCountId);
			payIndex->SetResetTiming(K_ADVANCE);
			payIndex->SetIntRule(K_ADJUSTED);
		}
		else
		{
			payIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(payIndex, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: PayIndex is not of a good type");
				return ARM_KO;
			}			
		}

		ARM_ReferenceValue aWeight1(weight1);
		ARM_ReferenceValue aWeight2(weight2);


		// Construction of the VMS Spread Option
		//------------------------------------------------------

		newSpro = new ARM_SpreadOptionVMS((ARM_Date)sStartDate,
										  (ARM_Date)sEndDate,
										  capOrFloorId,
										  refStrike,
										  payIndex,
										  payFixedRate,
										  CMSIndexes1,
										  CMSIndexes2,
										  &aWeight1,
										  &aWeight2,
										  dayCountId,
										  resetFreqId,
										  payFreqId,
										  resetTimingId,
										  payTimingId,
										  ccy,
										  resetGap,
										  spread1,
										  spread2,
										  Fixing1Ref,
										  Fixing2Ref,
										  intRule,
										  stubRule,
										  NULL,
										  NULL,
										  slopeFlag,
										  cptStrikeMethod,
										  refMargin,
										  payWeight,
										  FixingPayRef,
										  freezeFixing,
										  computedFormula);


		// memory release
		//------------------------------------------------------

		if (tmpRefStrike)
		{
			delete tmpRefStrike;
		}

		tmpRefStrike = NULL;

		if (tmpPayFixedRate)
		{
			delete tmpPayFixedRate;
		}

		tmpPayFixedRate = NULL;

		if (tmpRefMargin)
		{
			delete tmpRefMargin;
		}

		tmpRefMargin = NULL;


		// Second check
		//------------------------------------------------------

		if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor VMS");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
				{
					delete newSpro;
				}

				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOptionVMS *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
				{
					delete newSpro;
				}

				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			return ARM_OK;
		}
    }
 
    catch(Exception& x)
    { 
		if (newSpro)
		{
			delete newSpro;
		}

		newSpro = NULL;

		if (tmpRefStrike)
		{
			delete tmpRefStrike;
		}

		tmpRefStrike = NULL;

		if (tmpPayFixedRate)
		{
			delete tmpPayFixedRate;
		}

		tmpPayFixedRate = NULL;

		if (tmpRefMargin)
		{
			delete tmpRefMargin;
		}

		tmpRefMargin = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
}

long ARMLOCAL_QUANTOSPREADCORRIDOR (double startDate,
									double endDate,
									long capOrFloor,
									long strike_type,
									double strikes,
									long liborIdx1Id,
									long liborIdx2Id,
									long liborIdxPId,
									long pFixedRate,
									double Idx1weight,
									double Idx2weight,
									double IdxPweight,
									long Idx1fixingsId,
									long Idx2fixingsId,
									long IdxPfixingsId,
									double idx1spread,
									double idx2spread,
									long slopeFlag,
									int cptStrikeMethod,
									int computedFormula,
									long CcyId,
									long dayCountId,
									long resetFreqId,
									long payFreqId,
									long resetTimingId,
									long payTimingId,
									long intRuleId,
									long stubRuleId,
									long resetGap,
									char* resetCal,
									char* payCal,
									char* payIndexResetCal,
									long fwdRule,
									char* refDate,
									long notionalId,
									int freezeFixing,
									ARM_result& result,
									long objId)
{
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	long sproId;

	char sStartDate[11];
	char sEndDate[11];

	ARM_ReferenceValue* refStrikes = NULL;
	ARM_ReferenceValue* tmpRefStrikes = NULL;

	ARM_IRIndex* refIndex1 = NULL;
	ARM_IRIndex* tmpRefIndex1 = NULL;

	ARM_IRIndex* refIndex2 = NULL;
	ARM_IRIndex* tmpRefIndex2 = NULL;

	ARM_IRIndex* refIndexP = NULL;

	ARM_ReferenceValue* fixedRateP = NULL;

	ARM_ReferenceValue* refIdx1fixings = NULL;
	ARM_ReferenceValue* refIdx2fixings = NULL;
	ARM_ReferenceValue* refIdxPfixings = NULL;

	ARM_ReferenceValue* aWeight1 = NULL;
	ARM_ReferenceValue* aWeight2 = NULL;
	ARM_ReferenceValue* aWeightP = NULL;
	
	ARM_ReferenceValue* notional = NULL;

	ARM_Currency* ccy = NULL;
	ARM_Currency* ccy1 = NULL;
	ARM_Currency* ccy2 = NULL;
	ARM_Currency* ccyP = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);
		
		if ( CcyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else if ( CcyId == ARM_FRF_CCY_OBJECT )
			ccy = ARM_FRF_CURRENCY;
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( liborIdx1Id < 0 )
		{
			refIndex1 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex1 = refIndex1;
		}
		else
		{
			refIndex1 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx1Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex1, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 is not of a good type");
				return ARM_KO;
			}
			refIndex1->SetPayTiming(payTimingId);
			refIndex1->SetResetTiming(resetTimingId);
			refIndex1->SetResetFrequency(resetFreqId);
			refIndex1->SetPayFrequency(payFreqId);
			refIndex1->SetResetGap(resetGap);
		}

		if ( liborIdx2Id < 0 )
		{
			refIndex2 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex2 = refIndex2;
		}
		else
		{
			refIndex2 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex2, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index2 is not of a good type");
				return ARM_KO;
			}
			refIndex2->SetResetFrequency(resetFreqId);
			refIndex2->SetPayFrequency(payFreqId);
			refIndex2->SetResetTiming(resetTimingId);
			refIndex2->SetPayTiming(payTimingId);
			refIndex2->SetResetGap(resetGap);
		}

		if ( liborIdxPId < 0 )
		{
			refIndexP = NULL;
		}
		else
		{
			refIndexP = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdxPId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndexP, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: IndexP is not of a good type");
				return ARM_KO;
			}
			refIndexP->SetPayFrequency(payFreqId);
			refIndexP->SetPayTiming(payTimingId);
		}

		if (strike_type)
		{
			refStrikes = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(strikes);
		}
		else
		{
			refStrikes = new ARM_ReferenceValue(strikes);
			tmpRefStrikes = refStrikes;
		}
		
		if( pFixedRate >= 0 )
		{
			fixedRateP = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(pFixedRate);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixedRateP, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: fixed rate for paiement leg is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			fixedRateP = NULL;
		}

		if( Idx1fixingsId >= 0 )
		{
			refIdx1fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx1fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx1fixings = NULL;
		}
		
		if( Idx2fixingsId >= 0)
		{
			refIdx2fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx2fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx2fixings = NULL;
		}

		if( IdxPfixingsId >= 0 )
		{
			refIdxPfixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(IdxPfixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdxPfixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdxPfixings = NULL;
		}

		if( notionalId >=0 )
		{
			notional = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(notional, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Notional is not of a good type");
				return ARM_KO;

			}
		}
		else
		{
			notional = NULL;
		}

		aWeight1 = new ARM_ReferenceValue(Idx1weight);
		aWeight2 = new ARM_ReferenceValue(Idx2weight);
		aWeightP = new ARM_ReferenceValue(IdxPweight);

		if (refIndexP)
		{
			// Corridor spread option
			newSpro = new ARM_SpreadOption(	 (ARM_Date) sStartDate, 
										 (ARM_Date) sEndDate,
										 (int) capOrFloor, 
										 refStrikes, 
										 refIndex1, 
										 refIndex2, 
										 refIndexP,
										 fixedRateP,
										 aWeight1,
										 aWeight2,
										 aWeightP,
										 ccy,
										 refIdx1fixings,
										 refIdx2fixings,
										 refIdxPfixings,
										 idx1spread,
										 idx2spread,
										 (int) slopeFlag,
										 (int) cptStrikeMethod,
										 (int) computedFormula,
										 (int) dayCountId,
										 (int) resetFreqId,
										 (int) payFreqId,
										 (int) resetTimingId,
										 (int) payTimingId,
										 (int) intRuleId,
										 (int) stubRuleId,
										 (int) resetGap,
										 resetCal,
										 payCal,
										 payIndexResetCal,
										 (int) fwdRule,
										 refDate,
										 notional,
										 (int) freezeFixing);
		}
		else
		{
			// vanilla spread option
			newSpro = new ARM_SpreadOption(	(ARM_Date) sStartDate, 
											(ARM_Date) sEndDate,
											capOrFloor, 
											refStrikes,
											refIndex1,  
											refIndex2,  
											aWeight1, 
											aWeight2,
											dayCountId, 
											resetFreqId, 
											payFreqId, 
											resetTimingId, 
											payTimingId, 
											ccy,
											resetGap,
											refIdx1fixings,
											refIdx2fixings,
											intRuleId,
											stubRuleId,
											cptStrikeMethod,
											resetCal,
											payCal,
											computedFormula,
											fwdRule,
											refDate,
											10000,
											NULL,
											notional);
		}
		
	if (tmpRefIndex1)
		delete tmpRefIndex1;
	tmpRefIndex1 = NULL;

	if (tmpRefIndex2)
		delete tmpRefIndex2;
	tmpRefIndex2 = NULL;

	if (tmpRefStrikes)
		delete tmpRefStrikes;
	tmpRefStrikes = NULL;
	
	if (aWeight1)
		delete aWeight1;
	aWeight1 = NULL;

	if (aWeight2)
		delete aWeight2;
	aWeight2 = NULL;

	if (aWeightP)
		delete aWeightP;
	aWeightP = NULL;

	if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			
			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}
				
				sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				result.setLong(sproId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;
		
		if (tmpRefIndex1)
			delete tmpRefIndex1;
		tmpRefIndex1 = NULL;

		if (tmpRefIndex2)
			delete tmpRefIndex2;
		tmpRefIndex2 = NULL;

		if (tmpRefStrikes)
			delete tmpRefStrikes;
		tmpRefStrikes = NULL;
		
		if (aWeight1)
			delete aWeight1;
		aWeight1 = NULL;

		if (aWeight2)
			delete aWeight2;
		aWeight2 = NULL;

		if (aWeightP)
			delete aWeightP;
		aWeightP = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
	
}

long ARMLOCAL_CORRIDORDBLCONDITION(double startDate,
								   double endDate,
								   long DigitalCapFloorId,
								   long SpreadCapFloorId,
								   long digitalBarrier_Type,
								   double digitalBarrier,
								   long spreadBarrier_Type,
								   double spreadBarrier,
								   long payIndexId,
								   long payIndexMargin_Type,
								   double payIndexMargin,
								   long digitalIndex1Id,
								   long digitalIndex2Id,
								   long digitalIndex3Id,
								   long digitalIndex4Id,
								   double weight1,
								   double weight2,
								   double weight3,
								   double weight4,
								   long fixing1Id,
								   long fixing2Id,
								   long fixing3Id,
								   long fixing4Id,
								   long fixingPayId,
								   long dayCountId,
								   long resetFreqId,
								   long payFreqId,
								   long resetTimingId,
								   long payTimingId,
								   long CcyId,
								   long resetGap,
								   double spread,
								   double digitalSpread,
								   long intRuleId,
								   long stubRuleId,
								   char* resetCal,
								   char* payCal,
								   double payWeight,
								   int freezeFixing,
								   long digitalSpreadcorrelId,
								   long shiftVolCorrel_type,
								   double ShiftVolCorrelId,
								   bool isCorrelCorrection,
								   ARM_result& result,
								   long objId)
{
	long corridorId;

	ARM_Currency*	ccy = NULL;
	ARM_ReferenceValue* refDigitalBarrier = NULL;
	ARM_ReferenceValue*	refSpreadBarrier  = NULL;
	ARM_ReferenceValue* refPayFixedRate   = NULL;
	ARM_ReferenceValue* refPayIndexMargin = NULL;

	ARM_ReferenceValue* refFixing1 = NULL;
	ARM_ReferenceValue* refFixing2 = NULL;
	ARM_ReferenceValue* refFixing3 = NULL;
	ARM_ReferenceValue* refFixing4 = NULL;
	ARM_ReferenceValue* refFixingPay = NULL;
	ARM_IRIndex* refPayIndex = NULL;
	ARM_SpreadOption* oldCorridorDbl = NULL;
	ARM_SpreadOption* newCorridorDbl = NULL;
	ARM_VolLInterpol* volCurveCorrel = NULL;
	ARM_VolLInterpol* shiftCorrelVol = NULL;
	bool bDigitalBarrier = false;
	bool bSpreadBarrier	 = false;
	bool bPayFixedRate   = false;
	bool bPayIndexMargin = false;
	bool bPayIndex		 = false;
	bool bShiftVol		 = false;


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

		if ( CcyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else if ( CcyId == ARM_FRF_CCY_OBJECT )
			ccy = ARM_FRF_CURRENCY;
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		if( digitalSpreadcorrelId != ARM_NULL_OBJECT )
		{
			volCurveCorrel = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(digitalSpreadcorrelId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volCurveCorrel, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: VolCorrel is not of a good type");
				return ARM_KO;
			}
		}

		if ( digitalBarrier_Type == 1 )
		{
			refDigitalBarrier = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(digitalBarrier));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refDigitalBarrier, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Digital barrier is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refDigitalBarrier = new ARM_ReferenceValue(digitalBarrier);
			bDigitalBarrier = true;
		}

		if ( spreadBarrier_Type == 1 )
		{
			refSpreadBarrier = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(spreadBarrier));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refSpreadBarrier, ARM_REFERENCE_VALUE) == 0)
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;

				result.setMsg ("ARM_ERR: Spread barrier is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refSpreadBarrier = new ARM_ReferenceValue(spreadBarrier);
			bSpreadBarrier = true;
		}

		refPayFixedRate = new ARM_ReferenceValue(0.);
		bPayFixedRate = true;


		if ( payIndexMargin_Type == 1 )
		{
			refPayIndexMargin = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payIndexMargin));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refPayIndexMargin, ARM_REFERENCE_VALUE) == 0)
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				result.setMsg ("ARM_ERR: Payment index margin is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refPayIndexMargin = new ARM_ReferenceValue(payIndexMargin);
			bPayIndexMargin = true;
		}

		if(payIndexId ==  ARM_NULL_OBJECT)
		{
			refPayIndex = new ARM_IRIndex(LIBOR3M, payFreqId, payFreqId, ccy, dayCountId);
			refPayIndex->SetResetTiming(K_ADVANCE);
			refPayIndex->SetIntRule(K_ADJUSTED);
			bPayIndex = true;
		}
		else
		{
			refPayIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refPayIndex, ARM_IRINDEX) == 0)
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				result.setMsg ("ARM_ERR: PayIndex is not of a good type");
				return ARM_KO;
			}			
		}

		if (fixing1Id != ARM_NULL_OBJECT)
		{
		   refFixing1 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refFixing1, ARM_REFERENCE_VALUE) == 0 )
		   {
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: fixing1 reference value is not of a good type");
				return ARM_KO;
			}
		}

		if (fixing2Id != ARM_NULL_OBJECT)
		{
		   refFixing2 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refFixing2, ARM_REFERENCE_VALUE) == 0 )
		   {
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: fixing2 reference value is not of a good type");
				return ARM_KO;
			}
		}

		if (fixing3Id != ARM_NULL_OBJECT)
		{
		   refFixing3 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing3Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refFixing3, ARM_REFERENCE_VALUE) == 0 )
		   {
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: fixing3 reference value is not of a good type");
				return ARM_KO;
			}
		}

		if (fixing4Id != ARM_NULL_OBJECT)
		{
		   refFixing4 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing4Id);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refFixing4, ARM_REFERENCE_VALUE) == 0 )
		   {
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: fixing3 reference value is not of a good type");
				return ARM_KO;
			}
		}


		if (fixingPayId != ARM_NULL_OBJECT)
		{
		   refFixingPay = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixingPayId);

		   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refFixingPay, ARM_REFERENCE_VALUE) == 0 )
		   {
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: fixingPay reference value is not of a good type");
				return ARM_KO;
			}
		}

		if ( shiftVolCorrel_type == 1 )
		{
			shiftCorrelVol = (ARM_VolLInterpol *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(ShiftVolCorrelId));

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(shiftCorrelVol, ARM_VOL_CURVE) == 0)
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				result.setMsg ("ARM_ERR: ShiftVolCorrel is not of a good type");
				return ARM_KO;
			}		
		}
		else
		{
			ARM_Date tmpDate;
			shiftCorrelVol = new ARM_VolFlat( tmpDate,ShiftVolCorrelId, ccy);
			bShiftVol = true;
		}

		ARM_ReferenceValue aWeight1(weight1);
		ARM_ReferenceValue aWeight2(weight2);
		ARM_ReferenceValue aWeight3(weight3);
		ARM_ReferenceValue aWeight4(weight4);


		newCorridorDbl = new ARM_CorridorDblCondition((ARM_Date)sStartDate, (ARM_Date)sEndDate,
													 DigitalCapFloorId,
													 SpreadCapFloorId,
													 refDigitalBarrier,
													 refSpreadBarrier,
													 refPayIndex,
													 refPayFixedRate,
													 refPayIndexMargin,
													 (ARM_INDEX_TYPE)digitalIndex1Id,
													 (ARM_INDEX_TYPE)digitalIndex2Id, 
													 (ARM_INDEX_TYPE)digitalIndex3Id,
													 (ARM_INDEX_TYPE)digitalIndex4Id,
													 &aWeight1, 
													 &aWeight2,
													 &aWeight3,
													 &aWeight4,
													 refFixing1,
													 refFixing2,
													 refFixing3,
													 refFixing4,
													 refFixingPay,
													 dayCountId, 
													 resetFreqId, 
													 payFreqId, 
													 resetTimingId, 
													 payTimingId, 
													 ccy,
													 resetGap,
													 fabs(spread),
													 fabs(digitalSpread),
													 intRuleId,
													 stubRuleId,
													 resetCal,
													 payCal,		 
													 payWeight,
													 freezeFixing,
													 10000,
													 volCurveCorrel,
													 shiftCorrelVol,
													 isCorrelCorrection);

		if(bDigitalBarrier)
			delete refDigitalBarrier;
		
		if(bSpreadBarrier)
			delete refSpreadBarrier;

		if(bPayFixedRate)
			delete refPayFixedRate;

		if(bPayIndexMargin)
			delete refPayIndexMargin;

		if(bPayIndex)
			delete refPayIndex;

		if(bShiftVol)
				delete shiftCorrelVol;

		if ( newCorridorDbl == NULL )
		{
			if(bDigitalBarrier)
				delete refDigitalBarrier;
			
			if(bSpreadBarrier)
				delete refSpreadBarrier;

			if(bPayFixedRate)
				delete refPayFixedRate;

			if(bPayIndexMargin)
				delete refPayIndexMargin;

			if(bPayIndex)
				delete refPayIndex;
			
			if(bShiftVol)
				delete shiftCorrelVol;

			result.setMsg ("ARM_ERR: Unable to instanciate a corridorDblCondition");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			corridorId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newCorridorDbl);

			if (corridorId == RET_KO)
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				if(bShiftVol)
					delete shiftCorrelVol;


				delete newCorridorDbl;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(corridorId);
			return ARM_OK;
		}
		else
		{
			oldCorridorDbl = (ARM_CorridorDblCondition *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldCorridorDbl, ARM_CORRIDORDBLCONDITION) == 1)
			{
				delete oldCorridorDbl;		
				
				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCorridorDbl, objId);
				
				return ARM_OK;
			}
			else
			{
				if(bDigitalBarrier)
					delete refDigitalBarrier;
				
				if(bSpreadBarrier)
					delete refSpreadBarrier;

				if(bPayFixedRate)
					delete refPayFixedRate;

				if(bPayIndexMargin)
					delete refPayIndexMargin;

				if(bPayIndex)
					delete refPayIndex;

				if(bShiftVol)
					delete shiftCorrelVol;


				delete newCorridorDbl;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if(bDigitalBarrier)
			delete refDigitalBarrier;
		
		if(bSpreadBarrier)
			delete refSpreadBarrier;

		if(bPayFixedRate)
			delete refPayFixedRate;

		if(bPayIndexMargin)
			delete refPayIndexMargin;

		if(bPayIndex)
			delete refPayIndex;

		if(bShiftVol)
				delete shiftCorrelVol;


		delete newCorridorDbl;
		
		x.DebugPrint();

		ARM_RESULT();
    }
}

long ARMLOCAL_CHECK_CORRELS(VECTOR<double>& expiries,
							VECTOR<long>& indexes,
							long modelId,
							VECTOR<double>& resultMatrix,
							long& resultNbRows,
							long& resultNbCols,
							ARM_result& result)
{
	ARM_BSModel *model = dynamic_cast<ARM_BSModel *>(LOCAL_PERSISTENT_OBJECTS->GetObject(modelId));

	if(!model)
	{
		result.setMsg ("ARM_ERR: model is not of a good type");
		return ARM_KO;
	}

	size_t i,j,nbVal=expiries.size();
	ARM_Vector armExpiries(nbVal);
	char date[11];
	for(i=0;i<nbVal;++i)
	{
		Local_XLDATE2ARMDATE(expiries[i],date);
		armExpiries.Elt(i)=(ARM_Date(date).GetJulian()-model->GetAsOfDate().GetJulian())/K_YEAR_LEN;
	}
	nbVal=indexes.size();
	ARM_Vector armTenors(nbVal);
	for(i=0;i<nbVal;++i)
		armTenors.Elt(i)=FromIndexTypeToTerm((ARM_INDEX_TYPE)indexes[i]);

	ARM_Matrix eigenValues;
	ARM_CorridorDblCondition::CheckCorrels(armExpiries,armTenors,*model,eigenValues);

	resultNbRows = eigenValues.GetNumLines();
	resultNbCols = eigenValues.GetNumCols();
	resultMatrix.resize(resultNbRows*resultNbCols);
	for(i=0;i<resultNbRows;++i)
		for(j=0;j<resultNbCols;++j)
			resultMatrix[i*resultNbCols+j]=eigenValues.Elt(i,j);

	return ARM_OK;
}

long ARMLOCAL_CHECK_RA2_CORRELS_STATUTS(long RA2Id,
										long modelId,
										VECTOR<double>& resultMatrix,
										long& resultNbRows,
										long& resultNbCols,
										ARM_result& result)
{
	ARM_CorridorDblCondition *RA2 = dynamic_cast<ARM_CorridorDblCondition *>( LOCAL_PERSISTENT_OBJECTS->GetObject(RA2Id));
	if(!RA2)
	{
		result.setMsg ("ARM_ERR: double condition range accrual is not of a good type");
		return ARM_KO;
	}
	if(modelId != ARM_NULL_OBJECT)
	{
		/// Only check correlation consistency using RA2 schedule and tenors
		ARM_BSModel *model = dynamic_cast<ARM_BSModel *>(LOCAL_PERSISTENT_OBJECTS->GetObject(modelId));

		size_t	i,j,nbCaplets = RA2->GetNumFlows();
		double*	resetDates  = RA2->GetSpreadLeg()->GetFirstLeg()->GetResetDates()->GetElt();
		ARM_Vector armExpiries(nbCaplets);
		for(i=0;i<nbCaplets;++i)
			armExpiries[i] = (resetDates[i]-model->GetStartDate().GetJulian())/K_YEAR_LEN;
		ARM_Vector armTenors(4);
		armTenors[0] = RA2->GetUnderlyingMaturity(RA2->GetSpreadLeg()->GetFirstLeg());
		armTenors[1] = RA2->GetUnderlyingMaturity(RA2->GetSpreadLeg()->GetSecondLeg());
		armTenors[2] = RA2->GetUnderlyingMaturity(RA2->GetSpreadDigital()->GetSpreadLeg()->GetFirstLeg());
		armTenors[3] = RA2->GetUnderlyingMaturity(RA2->GetSpreadDigital()->GetSpreadLeg()->GetSecondLeg());

		ARM_Matrix eigenValues;
		ARM_CorridorDblCondition::CheckCorrels(armExpiries,armTenors,*model,eigenValues);

		resultNbRows = eigenValues.GetNumLines();
		resultNbCols = eigenValues.GetNumCols();
		resultMatrix.resize(resultNbRows*resultNbCols);
		for(i=0;i<resultNbRows;++i)
			for(j=0;j<resultNbCols;++j)
				resultMatrix[i*resultNbCols+j]=eigenValues.Elt(i,j);
	}
	else
	{
		/// Output correlation status of the input RA2
		string txt("");
		switch(RA2->GetCorrelStatus())
		{
		case ARM_CorridorDblCondition::Consistent:
			txt += "Consistent correlations";
			break;

		case ARM_CorridorDblCondition::InconsistentNotCorrected:
			txt += "Inconsistent correlations NOT CORRECTED";
			break;

		case ARM_CorridorDblCondition::InconsistentCorrected:
			txt += "Inconsistent correlations CORRECTED";
			break;
		}
		result.setString(txt.c_str());
	}

	return ARM_OK;
}

long ARMLOCAL_SPREADDIGITALFLT(	double startDate,
								double endDate,
								long capOrFloorId,
								long strike_type,
								double strike,
								long payoffliborType1Id,
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
								long resetGap,
								double spread1,
								double spread2,
							    long fixing1_type,
							    long fixing1Id, 
							    VECTOR<double>& fixing1,
							    long fixing2_type,
							    long fixing2Id,
							    VECTOR<double>& fixing2,
								long intRule,
								long stubRule,
								long slopeFlag,
								int cptStrikeMethod,
								int computedFormula,
								ARM_result& result,
								long objId)
{
	long sproId;

	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	ARM_Currency *ccy = NULL;
	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;

	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());
	ARM_ReferenceValue* fixing2Ref = NULL;

	ARM_ReferenceValue aWeight1(weight1);
	ARM_ReferenceValue aWeight2(weight2);


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

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

		if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{

			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)payoffliborType1Id, 
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   &aWeight1,
										   &aWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   spread1,
										   spread2,
										   (ARM_Vector *) &fixing1Vect,
										   (ARM_Vector *) &fixing2Vect,
										   intRule,
										   stubRule,
										   NULL,
										   NULL,
										   slopeFlag,
										   cptStrikeMethod,
                                           NULL,
                                           NULL,
										   computedFormula);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption((ARM_Date)sStartDate,
										   (ARM_Date)sEndDate,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)payoffliborType1Id, 
										   (ARM_INDEX_TYPE)liborType1Id,
										   (ARM_INDEX_TYPE)liborType2Id,
										   &aWeight1,
										   &aWeight2,
										   dayCountId,
										   resetFreqId,
										   payFreqId,
										   resetTimingId,
										   payTimingId,
										   ccy,
										   resetGap,
										   spread1,
										   spread2,
										   (ARM_ReferenceValue *) fixing1Ref,
										   (ARM_ReferenceValue *) fixing2Ref,
										   intRule,
										   stubRule,
										   NULL,
										   NULL,
										   slopeFlag,
										   cptStrikeMethod,
                                           NULL,
                                           NULL,
										   computedFormula);
		}
		else
		{
		   result.setMsg ("ARM_ERR: fixing1 type should be the same as fixing2 type ");		
		}


		if (tmpRefStrike)
		   delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		ARM_RESULT();
    }
}

extern long ARMLOCAL_SPREADDIGITALFLTWithLegs(	long firstlegId,
												long secondlegId,
												long capOrFloorId,
												long strike_type,
												double strike,
												long payoffliborType1Id,
												double weight1,
												double weight2,
												double spread1,
												double spread2,
												long fixing1_type,
												long fixing1Id, 
												VECTOR<double>& fixing1,
												long fixing2_type,
												long fixing2Id,
												VECTOR<double>& fixing2,
												long slopeFlag,
												int cptStrikeMethod,
												int computedFormula,
												ARM_result& result,
												long objId)
{
	long sproId;
	
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;

	ARM_SwapLeg* firstLeg = NULL;
	ARM_SwapLeg* secondLeg = NULL;

	ARM_ReferenceValue *refStrike = NULL;
	ARM_ReferenceValue *tmpRefStrike = NULL;

	ARM_Vector fixing1Vect(fixing1.size());
	ARM_ReferenceValue* fixing1Ref = NULL;

	ARM_Vector fixing2Vect(fixing1.size());
	ARM_ReferenceValue* fixing2Ref = NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		firstLeg = (ARM_SwapLeg *)LOCAL_PERSISTENT_OBJECTS->GetObject(firstlegId);
		secondLeg = (ARM_SwapLeg *)LOCAL_PERSISTENT_OBJECTS->GetObject(secondlegId);

		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(strike);

			tmpRefStrike = refStrike;
		}

		if ( fixing1_type == 1 )
		{
		   if (fixing1Id != ARM_NULL_OBJECT)
		   {
			   fixing1Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing1Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing1Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing1_type == 0)
		{
			if ( fixing1.size() > 1 )
			{
				for (int i = 0; i < fixing1.size(); i++)
				{
					fixing1Vect.Elt(i) = fixing1[i];
				}
			}
			else
			{
				fixing1Vect.Elt(0) = fixing1[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1_type is not of a good type");
		}

		if ( fixing2_type == 1 )
		{
		   if (fixing2Id != ARM_NULL_OBJECT)
		   {
			   fixing2Ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(fixing2Id);

			   if ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(fixing2Ref, ARM_REFERENCE_VALUE) == 0 )
			   {
				  result.setMsg ("ARM_ERR: fixing reference value is not of a good type");
					
				  return(ARM_KO);
				}
		   }
		}
		else if(fixing2_type == 0)
		{
			if ( fixing2.size() > 1 )
			{
				for (int i = 0; i < fixing2.size(); i++)
				{
					fixing2Vect.Elt(i) = fixing2[i];
				}
			}
			else
			{
				fixing2Vect.Elt(0) = fixing2[0];
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing2_type is not of a good type");
		}

		ARM_ReferenceValue vWeight1(weight1);
		ARM_ReferenceValue vWeight2(weight2);

		if(fixing1_type==fixing2_type && fixing1_type == 0)
		{

			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)payoffliborType1Id, 
										   &vWeight1,
										   &vWeight2,
										   spread1,
										   spread2,
										   &fixing1Vect,
										   &fixing2Vect,
										   slopeFlag,
										   cptStrikeMethod,
										   NULL,
										   NULL,
										   computedFormula);
		}
		else if(fixing1_type==fixing2_type && fixing1_type == 1)
		{
			newSpro = new ARM_SpreadOption(firstLeg,
										   secondLeg,
										   capOrFloorId,
										   refStrike,
										   (ARM_INDEX_TYPE)payoffliborType1Id, 
										   &vWeight1,
										   &vWeight2,
										   spread1,
										   spread2,
										   fixing1Ref,
										   fixing2Ref,
										   slopeFlag,
										   cptStrikeMethod,
										   NULL,
										   NULL,
										   computedFormula);
		}
		else
		{
			result.setMsg ("ARM_ERR: fixing1 type should be the same as fixing2 type ");		
		}


		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		if (newSpro == NULL)
		{
			result.setMsg ("ARM_ERR: Swaption is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				return ARM_OK;
			}

			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;

		if (tmpRefStrike)
			delete tmpRefStrike;
		tmpRefStrike = NULL;

		ARM_RESULT();
    }
}

long ARMLOCAL_QUANTOSPREADDIGITALFLT (	double startDate,
										double endDate,
										long capOrFloor,
										long strike_type,
										double strikes,
										long liborIdx1Id,
										long liborIdx2Id,
										long liborIdxPId,
										double Idx1weight,
										double Idx2weight,
										long Idx1fixingsId,
										long Idx2fixingsId,
										double spread1,
										double spread2,
										long slopeFlag,
										int cptStrikeMethod,
										int computedFormula,
										long CcyId,
										long dayCountId,
										long resetFreqId,
										long payFreqId,
										long resetTimingId,
										long payTimingId,
										long intRuleId,
										long stubRuleId,
										long resetGap,
										char* resetCal,
										char* payCal,
										long fwdRule,
										char* refDate,
										long notionalId,
										ARM_result& result,
										long objId)
{
	ARM_SpreadOption* spro = NULL;
	ARM_SpreadOption* newSpro = NULL;
	long sproId;

	char sStartDate[11];
	char sEndDate[11];

	ARM_ReferenceValue* refStrikes = NULL;
	ARM_ReferenceValue* tmpRefStrikes = NULL;

	ARM_IRIndex* refIndex1 = NULL;
	ARM_IRIndex* tmpRefIndex1 = NULL;

	ARM_IRIndex* refIndex2 = NULL;
	ARM_IRIndex* tmpRefIndex2 = NULL;

	ARM_IRIndex* refIndexP = NULL;
	ARM_IRIndex* tmpRefIndexP = NULL;

	ARM_ReferenceValue* refIdx1fixings = NULL;
	ARM_ReferenceValue* refIdx2fixings = NULL;
	
	ARM_ReferenceValue* aWeight1 = NULL;
	ARM_ReferenceValue* aWeight2 = NULL;
	
	ARM_ReferenceValue* notional = NULL;

	ARM_Currency* ccy = NULL;
	ARM_Currency* ccy1 = NULL;
	ARM_Currency* ccy2 = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);
		
		if ( CcyId == ARM_NULL_OBJECT )
			ccy = ARM_DEFAULT_CURRENCY;
		else if ( CcyId == ARM_FRF_CCY_OBJECT )
			ccy = ARM_FRF_CURRENCY;
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( liborIdx1Id < 0 )
		{
			refIndex1 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex1 = refIndex1;
		}
		else
		{
			refIndex1 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx1Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex1, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 is not of a good type");
				return ARM_KO;
			}
			refIndex1->SetPayTiming(payTimingId);
			refIndex1->SetResetTiming(resetTimingId);
			refIndex1->SetResetFrequency(resetFreqId);
			refIndex1->SetPayFrequency(payFreqId);
			refIndex1->SetResetGap(resetGap);
		}

		if ( liborIdx2Id < 0 )
		{
			refIndex2 = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndex2 = refIndex2;
		}
		else
		{
			refIndex2 = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdx2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex2, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index2 is not of a good type");
				return ARM_KO;
			}
			refIndex2->SetResetFrequency(resetFreqId);
			refIndex2->SetPayFrequency(payFreqId);
			refIndex2->SetResetTiming(resetTimingId);
			refIndex2->SetPayTiming(payTimingId);
			refIndex2->SetResetGap(resetGap);
		}

		if ( liborIdxPId < 0 )
		{
			refIndexP = new ARM_IRIndex (GetDefaultIndexFromCurrency("EUR"));
			tmpRefIndexP = refIndexP;
		}
		else
		{
			refIndexP = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(liborIdxPId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refIndex1, ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 is not of a good type");
				return ARM_KO;
			}
			refIndexP->SetPayTiming(payTimingId);
			refIndexP->SetResetTiming(resetTimingId);
			refIndexP->SetResetFrequency(resetFreqId);
			refIndexP->SetPayFrequency(payFreqId);
		}

		if (strike_type)
		{
			refStrikes = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(strikes);
		}
		else
		{
			refStrikes = new ARM_ReferenceValue(strikes);
			tmpRefStrikes = refStrikes;
		}

		if( Idx1fixingsId >= 0 )
		{
			refIdx1fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx1fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx1fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx1fixings = NULL;
		}
		
		if( Idx2fixingsId >= 0)
		{
			refIdx2fixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(Idx2fixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refIdx2fixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Index1 fixings is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refIdx2fixings = NULL;
		}

		if( notionalId >=0 )
		{
			notional = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(notional, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Notional is not of a good type");
				return ARM_KO;

			}
		}
		else
		{
			notional = NULL;
		}

		aWeight1 = new ARM_ReferenceValue(Idx1weight);
		aWeight2 = new ARM_ReferenceValue(Idx2weight);
		
		
		// digital spread option
		newSpro = new ARM_SpreadOption(	(ARM_Date) sStartDate, 
										(ARM_Date) sEndDate,
										capOrFloor, 
										refStrikes,
										refIndexP, 
										refIndex1,  
										refIndex2, 
										aWeight1, 
										aWeight2,
										dayCountId, 
										resetFreqId, 
										payFreqId, 
										resetTimingId, 
										payTimingId, 
										ccy,
										resetGap,
										spread1,
										spread2,
										refIdx1fixings,
										refIdx2fixings,
										intRuleId,
										stubRuleId,
										resetCal,
										payCal,
										slopeFlag,
										cptStrikeMethod,
										NULL,
										NULL,
										computedFormula,
										fwdRule,
										refDate,
										10000,
										NULL,
										notional);
		
	if (tmpRefIndex1)
		delete tmpRefIndex1;
	tmpRefIndex1 = NULL;

	if (tmpRefIndex2)
		delete tmpRefIndex2;
	tmpRefIndex2 = NULL;

	if (tmpRefIndexP)
		delete tmpRefIndexP;
	tmpRefIndexP = NULL;
	
	if (tmpRefStrikes)
		delete tmpRefStrikes;
	tmpRefStrikes = NULL;
	
	if (aWeight1)
		delete aWeight1;
	aWeight1 = NULL;

	if (aWeight2)
		delete aWeight2;
	aWeight2 = NULL;

	if ( newSpro == NULL )
		{
			result.setMsg ("ARM_ERR: Unable to instanciate a spreadoption of type corridor");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) newSpro);

			if (sproId == RET_KO)
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			
			result.setLong(sproId);

			return ARM_OK;
		}
		else
		{
			spro = (ARM_SpreadOption *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spro, ARM_SPREADOPTION) == 1)
			{
				if (spro)
				{
					delete spro;
					spro = NULL;
				}
				
				sproId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newSpro, objId);

				result.setLong(sproId);

				return ARM_OK;
		
			}
			else
			{
				if (newSpro)
					delete newSpro;
				newSpro = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    { 
		if (newSpro)
			delete newSpro;
		newSpro = NULL;
		
		if (tmpRefIndex1)
			delete tmpRefIndex1;
		tmpRefIndex1 = NULL;

		if (tmpRefIndex2)
			delete tmpRefIndex2;
		tmpRefIndex2 = NULL;

		if (tmpRefIndexP)
			delete tmpRefIndexP;
		tmpRefIndexP = NULL;
		
		if (tmpRefStrikes)
			delete tmpRefStrikes;
		tmpRefStrikes = NULL;

		if (aWeight1)
			delete aWeight1;
		aWeight1 = NULL;

		if (aWeight2)
			delete aWeight2;
		aWeight2 = NULL;

		x.DebugPrint();

		ARM_RESULT();
    }
	
}

long ARMLOCAL_CapLetPrice (long secId, long modId, long numEx,
					  ARM_result& result)
{

    ARM_Security* sec=NULL;
    ARM_Model* mod=NULL;
    double proba;
	CCString msg ("");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	try
	{
 
		    sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec, ARM_FLEXIBLECAPFLOOR) == 0)
			{

				result.setMsg ("ARM_ERR: sticky Strikes is not of a good type");
				return ARM_KO;
			}
 
			mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);
			
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
			{

				result.setMsg ("ARM_ERR: sticky Strikes is not of a good type");
				return ARM_KO;
			}

	        ((ARM_FlexibleCapFloor *) sec)->SetModel(mod);
 
			proba = ((ARM_FlexibleCapFloor *) sec)->CptCapLetPrice(numEx);
  
			result.setDouble(proba);

			return ARM_OK;
	}
	
    catch(Exception& x)
    {
        x.DebugPrint();
 
		ARM_RESULT();
    }

	return ARM_KO;
}


long ARMLOCAL_RATCHET (long swapLegId,
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
	long ratchetId;

	ARM_Ratchet* createdRatchet = NULL;
	ARM_Ratchet* oldRatchet = NULL;
	ARM_SwapLeg* swapLeg = NULL;

	int nbSpreads = spreadDates.size ();
	int nbCorrels = correlDates.size ();
	int nbFwdVols = fwdVolsDates.size ();

	int i;

	char myDate[11];
	double dDate;

	double dVectDates[ARM_NB_TERMS];
	double dVectValues[ARM_NB_TERMS];

	CCString msg (" ");

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

		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		for (i=0;i<nbSpreads;i++)
		{
			Local_XLDATE2ARMDATE(spreadDates[i],myDate);
			ARM_Date tmpDate (myDate);
			dDate = tmpDate.GetJulian();
			dVectDates[i] = dDate;
			dVectValues[i] = spreadValues[i];
		}

		ARM_Vector* vSpreadDates = new ARM_Vector(nbSpreads,dVectDates);
		ARM_Vector* vSpreadValues = new ARM_Vector(nbSpreads,dVectValues);
		ARM_ReferenceValue* vSpreads = new ARM_ReferenceValue(vSpreadDates,vSpreadValues);

		for (i=0;i<nbCorrels;i++)
		{
			Local_XLDATE2ARMDATE(correlDates[i],myDate);
			ARM_Date tmpDate (myDate);
			dDate = tmpDate.GetJulian();
			dVectDates[i] = dDate;
			dVectValues[i] = correlValues[i];
		}

		ARM_Vector* vCorrelDates = NULL;
		ARM_Vector* vCorrelValues = NULL;
		ARM_ReferenceValue* vCorrels = NULL;

		if (nbCorrels > 0)
		{
			vCorrelDates = new ARM_Vector(nbCorrels,dVectDates);
			vCorrelValues = new ARM_Vector(nbCorrels,dVectValues);
			vCorrels = new ARM_ReferenceValue(vCorrelDates,vCorrelValues);
		}

		for (i=0;i<nbFwdVols;i++)
		{
			Local_XLDATE2ARMDATE(fwdVolsDates[i],myDate);
			ARM_Date tmpDate (myDate);
			dDate = tmpDate.GetJulian();
			dVectDates[i] = dDate;
			dVectValues[i] = fwdVolsValues[i];
		}

		ARM_Vector* vFwdVolsDates = NULL;
		ARM_Vector* vFwdVolsValues = NULL;
		ARM_ReferenceValue* vFwdVols = NULL;

		if (nbFwdVols > 0)
		{
			vFwdVolsDates = new ARM_Vector(nbFwdVols,dVectDates);
			vFwdVolsValues = new ARM_Vector(nbFwdVols,dVectValues);
			vFwdVols = new ARM_ReferenceValue(vFwdVolsDates,vFwdVolsValues);
		}

		createdRatchet = new ARM_Ratchet(swapLeg,
										 vSpreads,
										 vCorrels,
										 vFwdVols,
										 capOrFloor,
										 strike);

		if (createdRatchet == NULL)
		{
			result.setMsg ("ARM_ERR: Ratchet is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			ratchetId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRatchet);

			if (ratchetId == RET_KO)
			{
				if (createdRatchet)
					delete createdRatchet;
				createdRatchet = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(ratchetId);

			return ARM_OK;
		}
		else
		{
			oldRatchet = (ARM_Ratchet *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldRatchet, ARM_RATCHET) == 1)
			{
				if (oldRatchet)
				{
					delete oldRatchet;
					oldRatchet = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRatchet, objId);

				return ARM_OK;
			}

			else
			{
				if (createdRatchet)
					delete createdRatchet;
				createdRatchet = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (createdRatchet)
			delete createdRatchet;
		createdRatchet = NULL;

		ARM_RESULT();
    }
}


long ARMLOCAL_DIGITAL (long swapLegId,
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
	long digId;

    ARM_Digital* newDigital=NULL;  
    ARM_Digital* oldDigital=NULL;  

	ARM_SwapLeg* swapLeg = NULL;

	ARM_ReferenceValue* refValueStrike = NULL;
	ARM_ReferenceValue* refValuePayoff = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swapLegId);

		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		if (strikeType == 0)
		{
			refValueStrike = new ARM_ReferenceValue(strike);
		}
		else
		{
			refValueStrike = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refValueStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: strike is not of a good type");
				return ARM_KO;
			}

		}

		if (payoffType == 0)
		{
			refValuePayoff = new ARM_ReferenceValue(payoff);
		}
		else
		{
			refValuePayoff = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(payoff));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refValuePayoff, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: payoff is not of a good type");
				return ARM_KO;
			}
		}

		newDigital = new ARM_Digital(swapLeg,isItCapOrFloorId,refValueStrike,spread1,spread2,refValuePayoff);

		if (strikeType == 0)
			delete refValueStrike;
		if (payoffType == 0)
			delete refValuePayoff;

/*
		if (refValuePayoff)
			delete refValuePayoff;
		refValuePayoff = NULL;
*/
		if (newDigital == NULL)
		{
			result.setMsg ("ARM_ERR: Digital is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			digId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newDigital);

			if (digId == RET_KO)
			{
				if (newDigital)
					delete newDigital;
				newDigital = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(digId);

			return ARM_OK;
		}
		else
		{
			oldDigital = (ARM_Digital *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldDigital, ARM_DIGITAL) == 1)
			{
				if (oldDigital)
				{
					delete oldDigital;
					oldDigital = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newDigital, objId);

				return ARM_OK;
			}

			else
			{
				if (newDigital)
					delete newDigital;
				newDigital = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newDigital)
			delete newDigital;
		newDigital = NULL;
/*
		if (refValueStrike)
			delete refValueStrike;
		refValueStrike = NULL;

		if (refValuePayoff)
			delete refValuePayoff;
		refValuePayoff = NULL;
*/
		ARM_RESULT();
    } 
}


long ARMLOCAL_DUALCAP (double startDate,
					  double endDate,
					  long isItCapOrFloorId,
					  long strikeType,
					  double strike,
					  long liborType1Id,
					  long liborType2Id,
					  long daycountId,
					  long FreqId,
					  long resetTiming1Id,
					  long resetTiming2Id,
					  double resetGap1,
					  double resetGap2,
					  long ccy1Id,
					  long ccy2Id,
					  long discountCcyId,
					  CCString resetCal1Id,
					  CCString resetCal2Id,
					  ARM_result& result,
					  long objId)
{

	long dualId;

    ARM_DualCap* newDualCap=NULL;  
    ARM_DualCap* oldDualCap=NULL;  

	ARM_ReferenceValue* refValueStrike = NULL;
	
	ARM_Currency *ccy1 = NULL;
	ARM_Currency *ccy2 = NULL;
	ARM_Currency *discountccy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];
	
	char* DresetCalName = new char[100];
	char* FresetCalName = new char[100];
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
 		Local_XLDATE2ARMDATE(endDate,sEndDate);

		if ( ccy1Id == ARM_NULL_OBJECT )
		{
			ccy1 = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccy1Id == ARM_FRF_CCY_OBJECT )
		{
			ccy1 = ARM_FRF_CURRENCY;
		}
		else
		{
			ccy1 = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccy1Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy1, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: PayIndex Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( ccy2Id == ARM_NULL_OBJECT )
		{
			ccy2 = ARM_DEFAULT_CURRENCY;
		}
		else if ( ccy2Id == ARM_FRF_CCY_OBJECT )
		{
			ccy2 = ARM_FRF_CURRENCY;
		}
		else
		{
			ccy2 = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccy2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy2, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: RangeIndex Currency is not of a good type");
				return ARM_KO;
			}
		}

		if ( discountCcyId == ARM_NULL_OBJECT )
		{
			discountccy = ARM_DEFAULT_CURRENCY;
		}
		else if ( discountCcyId == ARM_FRF_CCY_OBJECT )
		{
			discountccy = ARM_FRF_CURRENCY;
		}
		else
		{
			discountccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(discountccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: RangeIndex Currency is not of a good type");
				return ARM_KO;
			}
		}

		if (strikeType == 0)
		{
			refValueStrike = new ARM_ReferenceValue(strike);
		}
		else
		{
			refValueStrike = (ARM_ReferenceValue *)((ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(strike)))->Clone();

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refValueStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: strike is not of a good type");
				return ARM_KO;
			}

		}

		if (strcmp((const char*)resetCal1Id,"NULL") != 0)
			strcpy(DresetCalName,(const char*)resetCal1Id);
		else
		{
			if (DresetCalName)
				delete [] DresetCalName;
			DresetCalName = NULL;
		}

		if (strcmp((const char*)resetCal2Id,"NULL") != 0)
			strcpy(FresetCalName,(const char*)resetCal2Id);
		else
		{
			if (FresetCalName)
				delete [] FresetCalName;
			FresetCalName = NULL;
		}

		newDualCap = new ARM_DualCap((ARM_Date)sStartDate, (ARM_Date)sEndDate,
									isItCapOrFloorId, refValueStrike,
									(ARM_INDEX_TYPE)liborType1Id,
									(ARM_INDEX_TYPE)liborType2Id, 
									daycountId,
									FreqId,
									resetTiming1Id,
									resetTiming2Id,
									resetGap1,
									resetGap2,
									ccy1,
									ccy2,
									discountccy,
									DresetCalName,
									FresetCalName);

		
		if (refValueStrike)
			delete refValueStrike;
		refValueStrike = NULL;


		if (DresetCalName)
			delete [] DresetCalName;
		DresetCalName = NULL;

		if (FresetCalName)
			delete [] FresetCalName;
		FresetCalName = NULL;

		if (newDualCap == NULL)
		{
			result.setMsg ("ARM_ERR: Dual option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			dualId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newDualCap);

			if (dualId == RET_KO)
			{
				if (newDualCap)
					delete newDualCap;
				newDualCap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(dualId);

			return ARM_OK;
		}
		else
		{
			oldDualCap = (ARM_DualCap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldDualCap, ARM_DUALCAP) == 1)
			{
				if (oldDualCap)
				{
					delete oldDualCap;
					oldDualCap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newDualCap, objId);

				return ARM_OK;
			}

			else
			{
				if (newDualCap)
					delete newDualCap;
				newDualCap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newDualCap)
			delete newDualCap;
		newDualCap = NULL;

		if (refValueStrike)
			delete refValueStrike;
		refValueStrike = NULL;

		if (DresetCalName)
			delete [] DresetCalName;
		DresetCalName = NULL;

		if (FresetCalName)
			delete [] FresetCalName;
		FresetCalName = NULL;


		ARM_RESULT();
    }

}

long ARMLOCAL_GlobalCap(long swaplegId,
						long capOrFloor,
						double globalCapValue,
						long globalCapSpreadsId,
						long globalCapFixedRatesId,
						long globalCapBarriersId,
						double finalRatio,
						int nbIter,
						long globalCapPastFixingsId,
						ARM_result& result,
						long objId)
{
	long globalCapId;
	ARM_SwapLeg* swapLeg = NULL;
	ARM_GlobalCap* oldGlobalCap = NULL;
	ARM_GlobalCap* newGlobalCap = NULL;

	ARM_ReferenceValue* barriers = NULL;
	ARM_ReferenceValue* spreads = NULL;
	ARM_ReferenceValue* fixedRates = NULL;
	ARM_ReferenceValue* pastFixings = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		swapLeg = (ARM_SwapLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(swaplegId);
		if (LocalPersistent::LOCAL_IS_SWAPLEG_OK(swapLeg) == 0)
		{
			result.setMsg ("ARM_ERR: swapLeg is not a Swap Leg");
			return ARM_KO;
		}

		if (globalCapValue < 0.0)
		{
			result.setMsg ("ARM_ERR: global cap value must be positive");
			return ARM_KO;
		}
		if (finalRatio < 0.0)
		{
			result.setMsg ("ARM_ERR: final ratio must be positive");
			return ARM_KO;
		}
		if (nbIter < 0.0)
		{
			result.setMsg ("ARM_ERR: model iteration number must be positive");
			return ARM_KO;
		}

		spreads = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapSpreadsId);
		*spreads /= 100.0;
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(spreads, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: spreads must be a reference value");
			return ARM_KO;
		}

		fixedRates = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapFixedRatesId);
		*fixedRates /= 100;
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(fixedRates, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: fixed rates must be a reference value");
			return ARM_KO;
		}

		barriers = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapBarriersId);
		*barriers /= 100;
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(barriers, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: barriers must be a reference value");
			return ARM_KO;
		}

		if (globalCapPastFixingsId != ARM_NULL_OBJECT)
		{
			pastFixings = (ARM_ReferenceValue *)LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapPastFixingsId);
			*pastFixings /= 100;
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pastFixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: past fixings must be a reference value");
				return ARM_KO;
			}
		}

		newGlobalCap = new ARM_GlobalCap(swapLeg, 
										 capOrFloor,
										 globalCapValue,
										 spreads,
										 fixedRates,
										 barriers,
										 finalRatio,
										 nbIter,
										 pastFixings);

		if (newGlobalCap == NULL)
		{
			result.setMsg ("ARM_ERR: newGlobalCap is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			globalCapId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newGlobalCap);

			if (globalCapId == RET_KO)
			{
				if (newGlobalCap)
					delete newGlobalCap;
				newGlobalCap = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(globalCapId);

			return ARM_OK;
		}
		else
		{
			oldGlobalCap = (ARM_GlobalCap *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldGlobalCap, ARM_GLOBALCAP ) == 1)
			{
				if (oldGlobalCap)
				{
					delete oldGlobalCap;
					oldGlobalCap = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newGlobalCap, objId);

				return ARM_OK;
			}

			else
			{
				if (newGlobalCap)
					delete newGlobalCap;
				newGlobalCap = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (newGlobalCap)
			delete newGlobalCap;
		newGlobalCap = NULL;

		ARM_RESULT();
    } 
}
