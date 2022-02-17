#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <ARM\libarm\ARM_result.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <inst\iasecurity.h>


long ARMLOCAL_IASEC (long underlyingId,
					 long iaCtrlId,
					 long refValId,
					 long iaCtrlType,
					 ARM_result& result,
					 long objId)
{
	long idxId;

	ARM_IdxAmortSec* createdIaSec = NULL;
	ARM_IdxAmortSec* iaSec = NULL;

	ARM_Security *Und = NULL;
	ARM_Security *CtrlAsset = NULL;
	ARM_IARefVal *RefValue = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg("");

	try
	{
		Und = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(underlyingId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Und, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Underlying is not of a good type");
			return ARM_KO;
		}

		CtrlAsset = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(iaCtrlId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(CtrlAsset, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Ctrl Asset is not of a good type");
			return ARM_KO;
		}

		RefValue = (ARM_IARefVal *) LOCAL_PERSISTENT_OBJECTS->GetObject(refValId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RefValue, ARM_IAREFVAL) == 0)
		{
			result.setMsg ("ARM_ERR: IA RefValue is not of a good type");
			return ARM_KO;
		}

		createdIaSec = new ARM_IdxAmortSec(Und, CtrlAsset, RefValue, 
										   iaCtrlType);

		if (createdIaSec == NULL)
		{
			result.setMsg ("ARM_ERR: IA Security is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			idxId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIaSec);

			if (idxId == RET_KO)
			{
				if (createdIaSec)
					delete createdIaSec;
				createdIaSec = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(idxId);

			return ARM_OK;
		}
		else
		{
			iaSec = (ARM_IdxAmortSec *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(iaSec, ARM_IDXAMORTSEC) == 1)
			{
				if (iaSec)
				{
					delete iaSec;
					iaSec = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIaSec, objId);

				return ARM_OK;
			}

			else
			{
				if (createdIaSec)
					delete createdIaSec;
				createdIaSec = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIaSec)
			delete createdIaSec;
		createdIaSec = NULL;

		ARM_RESULT();
	}
}


		
