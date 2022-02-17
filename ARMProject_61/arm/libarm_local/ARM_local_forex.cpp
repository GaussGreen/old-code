#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_glob.h"
#include "ARM_local_persistent.h"


#include <inst\forex.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

long ARMLOCAL_FOREX (long LeftCcyId,
					 long RightCcyId,
                     double SpotValue,
					 ARM_result& result,
					 long objId)
{
	long forexId;

	ARM_Currency* c1 = NULL;
	ARM_Currency* c2 = NULL;

	ARM_Forex* createdFx = NULL;
	ARM_Forex* oldFx = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		c1 = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(LeftCcyId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(c1, ARM_CURRENCY) == 0)
		{
			result.setMsg ("ARM_ERR: Left Currency is not of a good type");
			return ARM_KO;
		}

		c2 = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(RightCcyId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(c1, ARM_CURRENCY) == 0)
		{
			result.setMsg ("ARM_ERR: Right Currency is not of a good type");
			return ARM_KO;
		}

		createdFx = new ARM_Forex(c1, c2, SpotValue);

		if (createdFx == NULL)
		{
			result.setMsg ("ARM_ERR: Forex is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			forexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFx);

			if (forexId == RET_KO)
			{
				if (createdFx)
					delete createdFx;
				createdFx = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(forexId);

			return ARM_OK;
		}
		else
		{
			oldFx = (ARM_Forex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFx, ARM_FOREX) == 1)
			{
				if (oldFx)
				{
					delete oldFx;
					oldFx = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFx, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFx)
					delete createdFx;
				createdFx = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();
 
		if (createdFx)
			delete createdFx;
		createdFx = NULL;
		
		ARM_RESULT();
	}
}

long ARMLOCAL_FOREX (const CCString& LeftCcyName,
					 const CCString& RightCcyName,
                     double SpotValue,
					 ARM_result& result,
					 long objId)
{
	long forexId;

	ARM_Currency* c1 = NULL;
	ARM_Currency* c2 = NULL;

	ARM_Forex* createdFx = NULL;
	ARM_Forex* oldFx = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		ARM_Currency c1((const char *) LeftCcyName);
		ARM_Currency c2((const char *) RightCcyName);

		createdFx = new ARM_Forex(&c1, &c2, SpotValue);
		
		if ( createdFx == NULL )
		{
		   result.setMsg ("ARM_ERR: Forex is null");
		
		   return ARM_KO;
		}

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();
			
			forexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFx);

			if ( forexId == RET_KO )
			{
				if (createdFx)
					delete createdFx;
				createdFx = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(forexId);

			return ARM_OK;
		}
		else
		{
			oldFx = (ARM_Forex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldFx, ARM_FOREX) == 1)
			{
				if (oldFx)
				{
					delete oldFx;
					oldFx = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdFx, objId);

				return ARM_OK;
			}
			else
			{
				if (createdFx)
					delete createdFx;
				createdFx = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

    catch(Exception& x)
    {
		x.DebugPrint();
 
		if (createdFx)
			delete createdFx;
		createdFx = NULL;

		ARM_RESULT();
	}
}

