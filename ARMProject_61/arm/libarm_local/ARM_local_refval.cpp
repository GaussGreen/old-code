#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_glob.h"
#include <libCCdate\CCdate.h>
#include <libCCtools++\CCstring.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <ARM\libarm\ARM_result.h>

#include "ARM_local_persistent.h"

#include <util\refvalue.h>
#include <util\ia3levrefval.h>


long ARMLOCAL_REFVALUE (VECTOR<double>& dates,
						VECTOR<double>& values,
						VECTOR<double>& values2,
						long valueType,
						long conversion,
						long calcMethod,
						ARM_result& result,
						long objId)
{
	long refvalId;

	ARM_ReferenceValue *refVal=NULL, *newRefVal=NULL;

	ARM_Vector* vDates = NULL;
	ARM_Vector* vValues = NULL;
	ARM_Vector* vValues2 = NULL;

	int i;
	int size = dates.size ();
	
	if( (size != values.size ()) || ((size != values2.size ()) && (values2.size () != 0)) )
	{
		result.setMsg ("ARM_ERR: dates and values must have same size");
		return ARM_KO;
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sDate[11];

	CCString msg ("");

	try
	{
		vDates = new ARM_Vector(size);
		vValues = new ARM_Vector(size);

		if (values2.size() != 0)
			vValues2 = new ARM_Vector(size);

		for(i = 0; i < size; i++)
		{
			Local_XLDATE2ARMDATE(dates[i],sDate);

			vDates->Elt(i) = ((ARM_Date) sDate).GetJulian();
			vValues->Elt(i) =  values[i];
			
			if (vValues2)
				vValues2->Elt(i) =  values2[i];
		}

		if (vValues2 == NULL)
			newRefVal = new ARM_ReferenceValue(vDates, vValues, valueType,
												conversion);
		else
			newRefVal = new ARM_ReferenceValue(vDates, vValues, vValues2, valueType,
												conversion);

		newRefVal->SetCalcMethod(calcMethod);

		// Pas de destruction car vecteurs non clonés à la construction											
/*
		if (vDates)
			delete vDates;
		vDates = NULL;

		if (vValues)
			delete vValues;
		vValues = NULL;
*/
		if (newRefVal == NULL)
		{
			result.setMsg ("ARM_ERR: refValue is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			refvalId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal);

			if (refvalId == RET_KO)
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(refvalId);

			return ARM_OK;
		}
		else
		{
			refVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refVal, ARM_REFERENCE_VALUE) == 1)
			{
				if (refVal)
				{
					delete refVal;
					refVal = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal, objId);

				return ARM_OK;
			}

			else
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newRefVal)
			delete newRefVal;
		newRefVal = NULL;

		if (vDates)
			delete vDates;
		vDates = NULL;

		if (vValues)
			delete vValues;
		vValues = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_CONSTREFVALUE (double value,
							 ARM_result& result,
							 long objId)
{
	long refvalId;

	ARM_ReferenceValue* refVal = NULL; 
	ARM_ReferenceValue* newRefVal = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		newRefVal = new ARM_ReferenceValue(value);
	
		if (newRefVal == NULL)
		{
			result.setMsg ("ARM_ERR: refValue is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			refvalId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal);

			if (refvalId == RET_KO)
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(refvalId);

			return ARM_OK;
		}
		else
		{
			refVal = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refVal, ARM_REFERENCE_VALUE) == 1)
			{
				if (refVal)
				{
					delete refVal;
					refVal = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newRefVal, objId);

				return ARM_OK;
			}

			else
			{
				if (newRefVal)
					delete newRefVal;
				newRefVal = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newRefVal)
			delete newRefVal;
		newRefVal = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_IATHREELEVREFVAL (double value,
								double level0,
								double amort0,
								double level1,
								double amort1,
								double level2,
								double amort2,
								ARM_result& result,
								long objId = -1)
{
	long refvalId;

	ARM_IA3LevRefVal* createdRefValue = NULL;
	ARM_IA3LevRefVal* refval = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		createdRefValue = new ARM_IA3LevRefVal(value, level0, amort0,
												level1, amort1,
												level2, amort2);
	
		if (createdRefValue == NULL)
		{
			result.setMsg ("ARM_ERR: refValue is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			refvalId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRefValue);

			if (refvalId == RET_KO)
			{
				if (createdRefValue)
					delete createdRefValue;
				createdRefValue = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(refvalId);

			return ARM_OK;
		}
		else
		{
			refval = (ARM_IA3LevRefVal *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval, ARM_IAREFVAL) == 1)
			{
				if (refval)
				{
					delete refval;
					refval = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRefValue, objId);

				return ARM_OK;
			}

			else
			{
				if (createdRefValue)
					delete createdRefValue;
				createdRefValue = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdRefValue)
			delete createdRefValue;
		createdRefValue = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_CptRefValue (long refvalId,
						   double date,
						   ARM_result& result)
{
	double dResult;
	ARM_ReferenceValue* refval=NULL;

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

		refval = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Reference Value is not of a good type");
			return ARM_KO;
		}

		dResult = refval->CptReferenceValue(myDate);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


long ARMLOCAL_SumRefValue (long refval1Id,
						   long refval2Id,
						   double coef,
						   ARM_result& result,
						   long objId)
{
	long refvalId;

	ARM_ReferenceValue* createdRefValue = NULL;
	ARM_ReferenceValue* oldRefValue = NULL;

	ARM_ReferenceValue* refval1 = NULL;
	ARM_ReferenceValue* refval2 = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		refval1 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(refval1Id);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval1,ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refvalue1 is not of a good type");
			return ARM_KO;
		}

		if (!(refval2Id == ARM_NULL_OBJECT))
		{
			refval2 = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(refval2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refval2,ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: refvalue2 is not of a good type");
				return ARM_KO;
			}
		}

		if (refval2 == NULL)
		{
				createdRefValue = (ARM_ReferenceValue*) refval1->Clone();
				*createdRefValue *= coef;
		}
		else
		{
			if (refval1->GetCalcMethod() == K_CONSTANT_REF)
			{
				createdRefValue = (ARM_ReferenceValue*) refval2->Clone();
				*createdRefValue *= coef;
				*createdRefValue += refval1->GetDiscreteValues()->Elt(0);
			}
			else if (refval2->GetCalcMethod() == K_CONSTANT_REF)
			{
				createdRefValue = (ARM_ReferenceValue*) refval1->Clone();
				*createdRefValue += coef * refval2->GetDiscreteValues()->Elt(0);
			}
			else if (refval1->GetCalcMethod() == refval2->GetCalcMethod())
			{
				ARM_Vector* vectDates1 = refval1->GetDiscreteDates();
				ARM_Vector* vectDates2 = refval2->GetDiscreteDates();

				ARM_Vector* newDateslonger;
				ARM_Vector* newDatesshorter;

				if (vectDates1->GetSize() >= vectDates2->GetSize())
				{
					newDateslonger = vectDates1;
					newDatesshorter = vectDates2;
				}
				else
				{
					newDateslonger = vectDates2;
					newDatesshorter = vectDates1;
				}

				int j = 0;

				for (int i = 0; i < newDatesshorter->GetSize(); i++)
				{
					while ( (j < newDateslonger->GetSize()) && (newDatesshorter->Elt(i) > newDateslonger->Elt(j)) )
						j++;

					if ( (j == newDateslonger->GetSize()) || (newDateslonger->Elt(j) != newDatesshorter->Elt(i))) 
					{
						result.setMsg ("ARM_ERR: check the dates of refval1 and refval2 ");
						return ARM_KO;
					}
				}

				if (newDateslonger == vectDates1)
				{
					createdRefValue = (ARM_ReferenceValue*) refval1->Clone();
				}
				else
				{
					createdRefValue = (ARM_ReferenceValue*) refval2->Clone();
				}

				for (i = 0; i < createdRefValue->GetDiscreteDates()->GetSize(); i++)
				{
					if (newDateslonger == vectDates1)
					{
						createdRefValue->GetDiscreteValues()->Elt(i) += coef * refval2->CptReferenceValue(createdRefValue->GetDiscreteDates()->Elt(i));
					}
					else
					{
						createdRefValue->GetDiscreteValues()->Elt(i) *= coef;
						createdRefValue->GetDiscreteValues()->Elt(i) += refval1->CptReferenceValue(createdRefValue->GetDiscreteDates()->Elt(i));
					}
				}
			}
			else
			{
				result.setMsg ("ARM_ERR: refval1 and refval2 have not the same calculation method");
				return ARM_KO;
			}
		}

		if (createdRefValue == NULL)
		{
			result.setMsg ("ARM_ERR: refValue is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			refvalId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRefValue);

			if (refvalId == RET_KO)
			{
				if (createdRefValue)
					delete createdRefValue;
				createdRefValue = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(refvalId);

			return ARM_OK;
		}
		else
		{
			oldRefValue = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldRefValue, ARM_REFERENCE_VALUE) == 1)
			{
				if (oldRefValue)
				{
					delete oldRefValue;
					oldRefValue = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdRefValue, objId);

				return ARM_OK;
			}

			else
			{
				if (createdRefValue)
					delete createdRefValue;
				createdRefValue = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdRefValue)
			delete createdRefValue;
		createdRefValue = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_DisplayRefValue(long refvalId,
							  bool isDate,
							  ARM_result& result)
{
	ARM_ReferenceValue* refval=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		refval = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(refval, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

        char buf[30];
		double dDate;

		if (refval->GetDiscreteDates() == NULL)
		{
			result.setLong(1);
			// datetoday
			ARM_Date date;
			
			if (isDate)
			{
				date.JulianToStrDate(buf);
				dDate = Local_ARMDATE2XLDATE(buf);
			}
			else
			{
				dDate = 0;
			}

			result.setArray(dDate,0);
			result.setArray(refval->GetDiscreteValues()->Elt(0),1);
		}
		else
		{
			result.setLong(refval->GetDiscreteValues()->GetSize());
			for (int i=0;i<refval->GetDiscreteValues()->GetSize();i++)
			{
				if (isDate)
				{
					ARM_Date date(refval->GetDiscreteDates()->Elt(i));
					date.JulianToStrDate(buf);
					dDate = Local_ARMDATE2XLDATE(buf);
				}
				else
				{
					dDate = refval->GetDiscreteDates()->Elt(i);
				}
				
				result.setArray(dDate,i);
				result.setArray(refval->GetDiscreteValues()->Elt(i),refval->GetDiscreteValues()->GetSize()+i);
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
