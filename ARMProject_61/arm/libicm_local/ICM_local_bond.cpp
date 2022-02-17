
#pragma warning(disable :4786 4541 4800 4250)

#include "firstToBeIncluded.h"
#include "CCdate.h"
#include "CCstring.h"

#include <math.h>

#include <ICMKernel\glob\icm_enums.h>
#include <ARMKernel\inst\swap.h>
#include <ICMKernel\inst\icm_bond.h>
#include <ICMKernel\inst\icm_frn.h>
#include <ARMKernel\ccy\currency.h>
#include <ARMKernel\mod\model.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>



long ICMLOCAL_BOND   (double	CouponRateIn,
					  double	Int_Accrual_DateIn,
					  double	MaturityIn,
					  int		frequencyIn,
					  double	First_Period_Reference_DateIn ,
					  double	NotionalAmountIn ,
					  	qPAYMENT_PREMIUM_LEG  AccOnDef,
					  CCString	Currency,
					  CCString  PayCalendar,
					  int		DayCount,
					  int		AccruedDayCount,
					  int		SettlementGap,
					  double	RedemptionValue,	
					  ARM_result& result,
					  long objId)
{
	long bondId;

	ICM_Bond* bond = NULL;
	ICM_Bond* newbond = NULL;
	// ARM_Currency* ccy = NULL;
	
	long resetFreq = -1;
	long payFreq = -1;
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pcurrentdate=new char[11];
	char* pMaturityIn=new char[11];
	//char* pRefDate=new char[11];

	CCString msg ("");

	//char* payCal = NULL;

	try
	{
		// Creation des ccy
		//char* tmp = (char*)Currency;
		//if (tmp)
		//{
		//	ccy = new ARM_Currency(tmp);
		//}
		//tmp = NULL;

		//if (strcmp((const char*)PayCalendar,"NULL"))
		//{
		//	payCal = new char[(ARM_NB_MAX_CAL*3)+1];
		//	strcpy(payCal,PayCalendar);
		//}


		Local_XLDATE2ARMDATE(MaturityIn,pMaturityIn);
		Local_XLDATE2ARMDATE(Int_Accrual_DateIn,pcurrentdate);

		ARM_Date refDate; 
		if (First_Period_Reference_DateIn == -1.0)
		{	
			//strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		Local_XLDATE2ARMDATE(First_Period_Reference_DateIn,refDate);

		
		newbond = new ICM_Bond( (CouponRateIn/100.),
								(ARM_Date) pcurrentdate,
								(ARM_Date) pMaturityIn,
								First_Period_Reference_DateIn==-1?0:&refDate, 
								 frequencyIn,
								 // pRefDate,
								 NotionalAmountIn,
								 (qPAYMENT_PREMIUM_LEG)AccOnDef,
								 CCSTringToSTLString(Currency),
								 CCSTringToSTLString(PayCalendar),
								 DayCount,
								 AccruedDayCount,
								 SettlementGap,
								 stubrule,
								 RedemptionValue);

//		if (ccy)
//			delete ccy;
//		ccy = NULL;
		
		if (pMaturityIn)
			delete [] pMaturityIn;
		pMaturityIn = NULL;

		if (pcurrentdate)
			delete [] pcurrentdate;
		pcurrentdate = NULL;

		// if (pRefDate)
		// 	delete [] pRefDate;
		// pRefDate = NULL;

		if (newbond == NULL)
		{
			result.setMsg ("ARM_ERR: liborSwap is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			bondId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbond);

			if (bondId == RET_KO)
			{
				if (newbond)
					delete newbond;
				newbond = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(bondId);

			return ARM_OK;
		}
		else
		{
			bond = (ICM_Bond *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(bond, ICM_BOND) == 1)
			{
				if (bond)
				{
					delete bond;
					bond = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newbond, objId);

				return ARM_OK;
			}
			else
			{
				if (newbond)
					delete newbond;
				newbond = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		if (newbond)
			delete newbond;
		newbond = NULL;

		if (pcurrentdate)
			delete [] pcurrentdate;
		pcurrentdate = NULL;

		if (pMaturityIn)
			delete [] pMaturityIn;
		pMaturityIn = NULL;

//		if (pRefDate)
//			delete [] pRefDate;
//		pRefDate = NULL;


		ARM_RESULT();
	}
}


long ICMLOCAL_YTOPRICE (long bondId, double settlement, double yield, ARM_result& result)
{
	double price;
	ARM_Security* security = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sSettlement = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		security = (ARM_Security*) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ICM_BOND) != 0)
		{
			ICM_Bond* bond = (ICM_Bond*) security;
			price = bond->YieldToPrice((ARM_Date) sSettlement, yield);

			if (sSettlement)
				delete [] sSettlement;
			sSettlement = NULL;
		}
		else
		{

			if (sSettlement)
				delete [] sSettlement;
			sSettlement = NULL;

			result.setMsg ("ARM_ERR: Asset is not of a good type");
			return ARM_KO;
		}

		result.setDouble(price);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sSettlement)
			delete [] sSettlement;
		sSettlement = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_PTOYIELD (long bondId, double settlement, double price, ARM_result& result)
{
	double yield;
	ARM_Security* security=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sSettlement = new char [11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(settlement,sSettlement);

		security = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(bondId);
 
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(security, ICM_BOND) != 0)
		{
			ICM_Bond* bond = (ICM_Bond*)security;

			yield = bond->PriceToYield((ARM_Date) sSettlement, price);

			if (sSettlement)
				delete [] sSettlement;
			sSettlement = NULL;
		}
		else
		{
			if (sSettlement)
				delete [] sSettlement;
			sSettlement = NULL;

			result.setMsg ("ARM_ERR: asset is not of a good type");
			return ARM_KO;
		}

		result.setDouble(yield);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sSettlement)
			delete [] sSettlement;
		sSettlement = NULL;

		ARM_RESULT();
	}
}




/*---- End Of File ----*/

// EOF %M%