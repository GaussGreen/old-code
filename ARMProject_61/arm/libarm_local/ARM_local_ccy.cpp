#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_glob.h"
#include "ARM_local_persistent.h"

#include <ccy\currency.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>


long ARMLOCAL_ISOCCY (const CCString& cname,
						ARM_result& result, long objId)
{
	long curId;

	ARM_Currency* createdCcy = NULL;
	ARM_Currency* prevCcy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg("");

	try
	{
        createdCcy = new ARM_Currency((const char*) cname);

		if (createdCcy == NULL)
		{
			result.setMsg ("ARM_ERR: currency is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			curId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCcy);

			if (curId == RET_KO)
			{
				if (createdCcy)
					delete createdCcy;
				createdCcy = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curId);

			return ARM_OK;
		}
		else
		{
			prevCcy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevCcy, ARM_CURRENCY) == 1)
			{
				if (prevCcy)
				{
					delete prevCcy;
					prevCcy = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCcy, objId);

				return ARM_OK;
			}
			else
			{
				if (createdCcy)
					delete createdCcy;
				createdCcy = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
 
	catch(Exception& x)
	{
        x.DebugPrint();

		if (createdCcy)
			delete createdCcy;
		createdCcy = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_CCY (const CCString& name, long idCurve,
				     double crossValue, long daycount,
					 ARM_result& result, long objId)
{
	long ccyId;

	ARM_Currency* newCcy=NULL;
	ARM_Currency* ccy=NULL;
	ARM_ZeroCurve* zc=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* tmp = NULL;
	CCString msg("");

	try
	{
		tmp = (char *) name;

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		newCcy = new ARM_Currency(tmp, zc);

		if (tmp)
			free(tmp);
		tmp = NULL;

		newCcy->SetCrossValue(crossValue);

		if ( daycount < 0 )
		{
		   ARM_Currency ccy((const char *) name);
			
           newCcy->SetMMDayCount(ccy.GetMMDayCount());
		}
		else
		{
			newCcy->SetMMDayCount(daycount);
		}

		if (newCcy == NULL)
		{
			result.setMsg ("ARM_ERR: currency is null");
			return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			ccyId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCcy);

			if ( ccyId == RET_KO )
			{
				if (newCcy)
				   delete newCcy;
				newCcy = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(ccyId);

			return ARM_OK;
		}
		else
		{
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 1)
			{
				if (ccy)
				{
					delete ccy;
					ccy = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCcy, objId);

				return ARM_OK;
			}
			else
			{
				if (newCcy)
					delete newCcy;
				newCcy = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();
 
		if (newCcy)
			delete newCcy;
		newCcy = NULL;
		
		if (tmp)
			free(tmp);
		tmp = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_GetSpotDays(const CCString&  ccy,
						  ARM_result& result)
{
	double dResult;
	ARM_Currency* myCcy=NULL;

	char * sCCY = (char*) ccy;
	CCString msg ("");

	try
	{
		myCcy = new ARM_Currency(sCCY);

		if (myCcy)
		{
			dResult = myCcy->GetSpotDays();
			result.setDouble(dResult);

			if (sCCY)
				delete sCCY;
			sCCY = NULL;

			if (myCcy)
				delete myCcy;
			myCcy = NULL;

			return ARM_OK;
		}

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (myCcy)
			delete myCcy;
		myCcy = NULL;

		return ARM_KO;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (myCcy)
			delete myCcy;
		myCcy = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_GetInfoFromCcy(long ccyId,
							 const CCString& type,
							 ARM_result& result)
{
	ARM_Currency* ccy = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ccy = (ARM_Currency*) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
		{
			result.setMsg ("ARM_ERR: Currency is not of a good type");
			return ARM_KO;
		}

		if (strcmp((const char*) type,"RESETCAL") == 0)
		{
			char resCal[5];
			ccy->CalcResetCal(resCal);
			result.setString(resCal);
		}
		else if (strcmp((const char*) type,"PAYCAL") == 0)
		{
			char payCal[5];
			ccy->CalcFloatPayCal(payCal);
			result.setString(payCal);
		}
		else
		{
			result.setMsg("Wrong data type. Valid is RESETCAL or PAYCAL");
			return ARM_KO;
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
