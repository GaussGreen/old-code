#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <libCCtools++\CCstring.h>

#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <inst\irfuture.h>


long ARMLOCAL_IRFUT (double delivery, long idUnderlying,
					   ARM_result& result, long objId)
{
	long irfutId;

	ARM_IRIndex* IrUnder=NULL;
    ARM_IRFuture* IrFut=NULL;
    ARM_IRFuture* createdIrFut=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char * sDelivery = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(delivery,sDelivery);

		IrUnder = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(idUnderlying);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IrUnder, ARM_IRINDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Underlying is not of a good type");
			return ARM_KO;
		}

		if (objId == -1)
		{
			createdIrFut = new ARM_IRFuture((ARM_Date) sDelivery, IrUnder);

			if (sDelivery)
				delete [] sDelivery;
			sDelivery = NULL;

			if (createdIrFut == NULL)
			{
				result.setMsg ("ARM_ERR: IrFut is null");
				return ARM_KO;
			}
			
			CREATE_GLOBAL_OBJECT();
			
			irfutId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrFut);

			if (irfutId == RET_KO)
			{
				if (createdIrFut)
					delete createdIrFut;
				createdIrFut = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irfutId);

			return ARM_OK;
		}
		else
		{
			IrFut = (ARM_IRFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IrFut, ARM_IRFUTURE) == 1)
			{
				IrFut->Set((ARM_Date) sDelivery, IrUnder);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrFut)
			delete createdIrFut;
		createdIrFut = NULL;

		if (sDelivery)
			delete [] sDelivery;
		sDelivery = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_THREE_MONTH_FUT (const CCString& delivery, long market, const CCString& ccy,
							     ARM_result& result, long objId)
{
	long irfutId;

    ARM_IRFuture* IrFut=NULL;
    ARM_IRFuture* createdIrFut=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sCCY = NULL;
	char * sDelivery = NULL;

	CCString msg ("");

	try
	{
		char* sCCY = ccy.GetStr();		
		char * sDelivery = delivery.GetStr();

        createdIrFut = new ARM_IRFuture(sDelivery, (ARM_MARKET) market, sCCY);

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (sDelivery)
			delete sDelivery;
		sDelivery = NULL;

		if (createdIrFut == NULL)
		{
			result.setMsg ("ARM_ERR: IrFut is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			irfutId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrFut);

			if (irfutId == RET_KO)
			{
				if (createdIrFut)
					delete createdIrFut;
				createdIrFut = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irfutId);

			return ARM_OK;
		}
		else
		{
			IrFut = (ARM_IRFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IrFut, ARM_IRFUTURE) == 1)
			{
				if (IrFut)
				{
					delete IrFut;
					IrFut = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrFut, objId);

				return ARM_OK;
			}
			else
			{
				if (createdIrFut)
					delete createdIrFut;
				createdIrFut = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrFut)
			delete createdIrFut;
		createdIrFut = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (sDelivery)
			delete sDelivery;
		sDelivery = NULL;

		ARM_RESULT();
	}
}
