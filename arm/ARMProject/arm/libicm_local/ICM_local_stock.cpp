

#pragma warning(disable :4786 4541 4800 4250)

#include "firstToBeIncluded.h"

#include "ICM_local_stock.h"

#include "CCdate.h"
#include "CCstring.h"

#include <math.h>

#include <ICMKernel\glob\icm_enums.h>
#include <ICMKernel\inst\icm_stock.h>
#include <ICMKernel\mod\defaulttree.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

extern long ICMLOCAL_STOCK   (double	FirmVolIn,
					  double	DividendYieldIn,
					  long TreeModelId,
					  double SpotIn,
					  ARM_result& result,
					  long objId)
{
	long stockId;

	ICM_Stock* stock = NULL;
	ICM_Stock* newstock = NULL;

	ICM_DefaultTree* tree=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		tree = (ICM_DefaultTree*) LOCAL_PERSISTENT_OBJECTS->GetObject(TreeModelId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(tree, ICM_DEFAULTTREEMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: Model is not of a good type, Default Tree expected !");
			return ARM_KO;
		}

		newstock = new ICM_Stock(FirmVolIn,
								DividendYieldIn);

		newstock->SetSpot(SpotIn,tree);

		if (newstock == NULL)
		{
			result.setMsg ("ARM_ERR: Stock is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			stockId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newstock);

			if (stockId == RET_KO)
			{
				if (newstock)
					delete newstock;
				newstock = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(stockId);

			return ARM_OK;
		}
		else
		{
			stock = (ICM_Stock*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stock, ICM_STOCK) == 1)
			{
				if (stock)
				{
					delete stock;
					stock = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newstock, objId);

				return ARM_OK;
			}
			else
			{
				if (newstock)
					delete newstock;
				newstock = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newstock)
			delete newstock;
		newstock = NULL;

		ARM_RESULT();
	}
};

extern long ICMLOCAL_STOCKCALLOPTION  (long	StockId,
					  double		StrikePriceIn,
					  double		FirstDateIn,
					  double		LastDateIn,
					  ARM_result& result,
					  long objId)
{
	long CallOptionId;

	ICM_StockCallOption* option = NULL;
	ICM_StockCallOption* newoption = NULL;

	ICM_Stock* stock=NULL;

	char*pFirstDate=new char[11];
	char*pLastDate=new char[11];

	Local_XLDATE2ARMDATE(FirstDateIn,pFirstDate);
	Local_XLDATE2ARMDATE(LastDateIn,pLastDate);

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		stock = (ICM_Stock*) LOCAL_PERSISTENT_OBJECTS->GetObject(StockId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(stock, ICM_STOCK) == 0)
		{
			result.setMsg ("ARM_ERR: Underlying stock is not of a good type, Defaultable stock expected !");
			return ARM_KO;
		}

		newoption = new ICM_StockCallOption(stock,
											StrikePriceIn,
											(ARM_Date)pFirstDate,
											(ARM_Date)pLastDate);

		if (newoption == NULL)
		{
			result.setMsg ("ARM_ERR: Stock Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			CallOptionId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newoption);

			if (CallOptionId == RET_KO)
			{
				if (newoption)
					delete newoption;
				newoption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(CallOptionId);

			return ARM_OK;
		}
		else
		{
			option = (ICM_StockCallOption*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(option, ICM_STOCKCALLOPTION) == 1)
			{
				if (option)
				{
					delete option;
					option = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newoption, objId);

				return ARM_OK;
			}
			else
			{
				if (newoption)
					delete newoption;
				newoption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newoption)
			delete newoption;
		newoption = NULL;

		ARM_RESULT();
	}
}
// EOF %M%
