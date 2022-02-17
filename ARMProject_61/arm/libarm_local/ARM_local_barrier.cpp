#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <inst\barrier.h>


long ARMLOCAL_CONSTBARRIER (long underlyingId,
					   long tAssetId,
					   double maturity,
					   double dbarrier1,
					   double dbarrier2,
					   long upDownDouble,
					   long inOut,
					   long triggerVar,
					   double rebate,
					   double firstX,
					   int isInArrear,
					   ARM_result& result,
					   long objId)
{
	long barId;

	ARM_Security* sec=NULL;
	ARM_Security* asset=NULL;

	ARM_Barrier* barrier=NULL;
	ARM_Barrier* newbarrier=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sMaturDate[11];
	char sFirstX[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(maturity,sMaturDate);
		Local_XLDATE2ARMDATE(firstX,sFirstX);

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underlyingId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		asset = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(tAssetId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(asset, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

			newbarrier = new ARM_Barrier(sec,
										 asset,
										 (ARM_Date)sMaturDate,
										 dbarrier1,
										 upDownDouble,
										 inOut,
										 triggerVar,
										 rebate,
										 (ARM_Date)sFirstX,
										 isInArrear,
										 dbarrier2);

        if (newbarrier == NULL)
		{
			result.setMsg ("ARM_ERR: barrier is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
            CREATE_GLOBAL_OBJECT();

            barId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbarrier);

			if (newbarrier == NULL)
			{
				result.setMsg ("ARM_ERR: Option is null");
				return ARM_KO;
			}

			if (barId == RET_KO)
			{
				if (newbarrier)
					delete newbarrier;
				newbarrier = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(barId);

			return ARM_OK;
		}
		else
		{
			barrier = (ARM_Barrier *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
 
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(barrier, ARM_BARRIER) == 1)
			{
				if (barrier)
					delete barrier;
                barrier = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbarrier, objId);

				return ARM_OK;
			}
   			else
			{
				if (newbarrier)
					delete newbarrier;
				newbarrier = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				

			result.setLong(objId);

			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

long ARMLOCAL_BARRIER (long underlyingId,
					   long tAssetId,
					   long xStyleId,
					   long refVal1Id,
					   long refVal2Id,
					   long upDownDouble,
					   long inOut,
					   long triggerVar,
					   double rebate,
					   int isInArrear,
					   ARM_result& result,
					   long objId)
{
	long barId;

	ARM_Security* sec=NULL;
	ARM_Security* asset=NULL;

	ARM_Barrier* barrier=NULL;
	ARM_Barrier* newbarrier=NULL;

    ARM_ExerciseStyle* barrierStyle = NULL; 
    ARM_ReferenceValue* Refbarrier1 = NULL;
	ARM_ReferenceValue* Refbarrier2 = NULL; 


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underlyingId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		asset = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(tAssetId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(asset, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		Refbarrier1 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refVal1Id);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Refbarrier1, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Barrier Refvalue 1 is not of a good type");
			return ARM_KO;
		}

		if (upDownDouble == K_DOUBLE)
		{
			Refbarrier2 = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refVal2Id);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Refbarrier2, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Barrier Refvalue 2 is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			Refbarrier2 = NULL;
		}

		barrierStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(xStyleId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(barrierStyle, ARM_EXERCISE_STYLE) == 0)
		{
			result.setMsg ("ARM_ERR: Exercise Style is not of a good type");
			return ARM_KO;
    	}

		
        newbarrier = new ARM_Barrier(sec,  asset, barrierStyle, Refbarrier1,
                                  upDownDouble, inOut, triggerVar, rebate, isInArrear, Refbarrier2);

       if (newbarrier == NULL)
		{
			result.setMsg ("ARM_ERR: barrier is null");
			return ARM_KO;
		}


		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			barId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbarrier);

			if (barId == RET_KO)
			{
				if (newbarrier)
					delete newbarrier;
				newbarrier = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(barId);

			return ARM_OK;
		}
		else
		{
		    barrier = (ARM_Barrier *) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
	
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(barrier, ARM_BARRIER) == 1)
			{
				if (barrier)
					delete barrier;
				barrier = NULL;

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbarrier, objId);

				return ARM_OK;
			}
			else
			{
				if (newbarrier)
					delete newbarrier;
				newbarrier = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newbarrier)
			delete newbarrier;
		newbarrier = NULL;

		ARM_RESULT();
	}

}



// EOF %M%