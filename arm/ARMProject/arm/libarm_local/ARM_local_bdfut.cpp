#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_bdfut.h"

#include <inst\bond.h>
#include <inst\portfolio.h>
#include <inst\bdfuture.h>

long ARMLOCAL_BDFUT (double delivery,
					 long underIsBd,
					 long underId,
					 double coupon,
					 double convFactor,
					 ARM_result& result,
					 long objId)
{
	CCString stringObjectId;
	ARM_Bond* BdUnder=NULL;
	ARM_Portfolio* pf=NULL;
	ARM_BondFuture* BdFut=NULL;
	double futurePrice;
	int id;
	ARM_BondFuture* createdBdFut=NULL;
	char * date_delivery = new char[11];
	ARM_Container* Cont = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(delivery,date_delivery);

		if ( underIsBd == 1 )
		{
			BdUnder = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdUnder, ARM_BOND) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
      	}
		else
		{
			pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			futurePrice = convFactor; // in this case convFactor is futurePrice

			Cont = pf->PfToCont();

		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			if ( underIsBd == 1 )
			{
				createdBdFut = new ARM_BondFuture((ARM_Date) date_delivery,
												  BdUnder,
												  coupon,
												  convFactor);
			}
			else
			{
				createdBdFut = new ARM_BondFuture(futurePrice,
												  (ARM_Date) date_delivery,
												  Cont,
												  coupon);

			}
			
			id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdBdFut);
			result.setLong(id);
		}
		else
		{
			BdFut = (ARM_BondFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdFut, ARM_BONDFUTURE) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			if ( underIsBd == 1 )
				BdFut->Set((ARM_Date) date_delivery, BdUnder, coupon, convFactor);
			else
				BdFut->Set(futurePrice, (ARM_Date) date_delivery, Cont, coupon);
		}

		if (date_delivery)
			free(date_delivery);
		date_delivery = NULL;

		return ARM_OK;
	}
	catch(Exception& x)
	{
        x.DebugPrint();
		ARM_RESULT();

		if (date_delivery)
			delete date_delivery;
		date_delivery = NULL;

	}
	return ARM_KO;
}

long ARMLOCAL_GetConversionFactor (long bdFutId,
							  long factId,
							  ARM_result& result)
{
    ARM_BondFuture* bdFut = NULL;
    double convFact = 0.0; 
    ARM_Vector* factors = NULL;
	CCString msg ("");
	
	try
	{
		bdFut = (ARM_BondFuture*) LOCAL_PERSISTENT_OBJECTS->GetObject(bdFutId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bdFut, ARM_BONDFUTURE) != 1)
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}

        factors = bdFut->GetConversionFactors();

		if (!factors) 
		{ 
		  convFact = bdFut->GetCTDConversionFactor();
		  result.setDouble(convFact);
		  return ARM_OK;
		}
 
        if (( factId >= 0 ) 
            &&
            ( factors != NULL ) 
            && 
            ( factId < factors->GetSize() )
           )
        {
           convFact = (*factors)[factId];
        }
        else
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Invalid range for conversion factor or NULL factors");
        }

		result.setDouble(convFact);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_KO;
}

long ARMLOCAL_GetCheapest (long bdFutId,
						   ARM_result& result)
{

    ARM_BondFuture* bdFut=NULL;
    ARM_Bond* cheapest=NULL;
    int cheapestId;
	CCString msg ("");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
 
    bdFut = (ARM_BondFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(bdFutId);
 
	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bdFut, ARM_BONDFUTURE) != 1)
	{
		result.setMsg ("ARM_ERR: previous object is not of a good type");
		return ARM_KO;
	}

	try
	{
        cheapest = bdFut->GetCheapest();
 
        cheapestId = LOCAL_PERSISTENT_OBJECTS->GetObjectIdFromPointer(cheapest);

		result.setDouble(cheapestId); //dam

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
        x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_KO;
}

long ARMLOCAL_GILT_NOTIONNAL_BUND (const CCString &delivery,
								   long underId,
								   long notioOrGilt,
								   long market,
								   ARM_result& result,
								   long objId)
{
    int id;
    ARM_Bond* BdUnder=NULL;
    ARM_BondFuture* BdFut=NULL;
    ARM_BondFuture* newBdFut=NULL;
	CCString msg ("");
	char * date_delivery = new char[11];

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	try
	{
		date_delivery= (char*) delivery;

		BdUnder = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdUnder, ARM_BOND) != 1)
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}

		if ( notioOrGilt == BDF_NOTIONNAL )
		{
			newBdFut = new ARM_BondFuture(date_delivery, BdUnder, MATIF,
										 BDF_NOTIONNAL);
		}
		else if ( notioOrGilt == BDF_GILT )
		{
			newBdFut = new ARM_BondFuture(date_delivery, BdUnder, LIFFE,
										 BDF_GILT);
		}
		else
		{
			newBdFut = new ARM_BondFuture(date_delivery, BdUnder, 
									 (ARM_MARKET) market, BDF_BUND);
		}

		if (date_delivery)
			free(date_delivery);
		date_delivery = NULL;

		if (newBdFut == NULL)
		{
			result.setMsg ("ARM_ERR: bdfut is null");
			return ARM_KO;
		}

		if (objId != -1)
		{
			BdFut = (ARM_BondFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdFut, ARM_BONDFUTURE) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
 
			(void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(newBdFut, objId);

			if (BdFut)
				delete BdFut;
			BdFut = NULL;
			
		}
		else
		{
			CREATE_GLOBAL_OBJECT();

			id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(newBdFut);
			
			result.setLong(id);			
		}

		return ARM_OK;

	}
	catch(Exception& x)
	{
		if (date_delivery)
			free(date_delivery);
		date_delivery = NULL;

        x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_KO;
}

long ARMLOCAL_GILT (const CCString &delivery,
					long underId,
					ARM_result& result,
					long objId)
{
	return ARMLOCAL_GILT_NOTIONNAL_BUND (delivery, underId, 2, 1, result, objId);
}

long ARMLOCAL_NOTIONNAL (const CCString &delivery,
						 long underId,
						 ARM_result& result,
						 long objId)
{
	return ARMLOCAL_GILT_NOTIONNAL_BUND (delivery, underId, 0, 0, result, objId);
}

long ARMLOCAL_BUND_LIFFE (const CCString &delivery,
						  long underId,
						  ARM_result& result,
						  long objId)
{
	return ARMLOCAL_GILT_NOTIONNAL_BUND (delivery, underId, 1, 1, result, objId);
}

long ARMLOCAL_BUND_DTB (const CCString &delivery,
						long underId,
						ARM_result& result,
						long objId)
{
	return ARMLOCAL_GILT_NOTIONNAL_BUND (delivery, underId, 1, 2, result, objId);
}

long ARMLOCAL_BDFUTBASKET (double delivery,
						   long underIsBd,
						   long underId,
						   double coupon,
						   double convFactor,
						   ARM_result& result,
						   long objId)
{
	CCString msg ("");
	ARM_Bond* BdUnder=NULL;
	ARM_Portfolio* pf=NULL;
    ARM_BondFuture* BdFut=NULL;
    double futurePrice;
	int id;
    ARM_BondFuture* createdBdFut=NULL;
	char * date_delivery = new char[11];
	ARM_Container* Cont = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	try
	{
		Local_XLDATE2ARMDATE(delivery,date_delivery);

		if ( underIsBd == 1 )
		{
			BdUnder = (ARM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdUnder, ARM_BOND) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(underId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}		

			futurePrice = convFactor; // in this case convFactor is futurePrice

			ARM_Container* Cont = pf->PfToCont();
		}
	
		if (objId != -1)
		{
			BdFut = (ARM_BondFuture *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BdFut, ARM_BONDFUTURE) != 1)
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}

			if ( underIsBd == 1 )
			{
				BdFut->Set((ARM_Date) date_delivery, BdUnder, coupon, convFactor);
			}
			else
			{
		        BdFut->Set(futurePrice, (ARM_Date) date_delivery, Cont, coupon);
			}
		}
		else
		{
			if ( underIsBd == 1 )
			{
				createdBdFut = new ARM_BondFuture((ARM_Date) date_delivery,
												  BdUnder,
												  coupon,
												  convFactor);
			}
			else
			{
				createdBdFut = new ARM_BondFuture(futurePrice,
												  (ARM_Date) date_delivery,
												  Cont,
												  coupon);

			}

		    CREATE_GLOBAL_OBJECT();

			id = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdBdFut);
 
			result.setLong(id);
		}

		if (date_delivery)
			free(date_delivery);
		date_delivery = NULL;

		return ARM_OK;
	}
	catch(Exception& x)
	{
		if (date_delivery)
			free(date_delivery);
		date_delivery = NULL;

        x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}



// EOF %M%