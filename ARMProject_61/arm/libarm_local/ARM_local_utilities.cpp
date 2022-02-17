#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"


#include <ARM\libarm_local\ARM_local_utilities.h>

#include <util/fromto.h>
#include <GP_Base/gpbase/gpmatrix.h>
#include <GP_Infra/gpinfra/gensecurity.h>
#include <GP_Infra/gpinfra/pricingadviser.h>
#include <GP_Infra/gpinfra/pricingmodel.h>
#include <GP_Calib/gpcalib/calibmethod.h>

#include <GP_Calculators/gpcalculators/csoquantocalculator.h>

using namespace ARM;



long ARMLOCAL_CQSOCREATE(double dStartDate, 
						 double dEndDate,
						 int iFrequency,
						 int iDayCounter,
						 int iIsAdjusted,
						 int iResetType,
						 int iResetGap,
						 CCString sResetCalendar,
						 CCString sPaymtCalendar,						
						 CCString domCurrency,
						 CCString forCurrency,
						 CCString index1,
						 CCString index2,						 
						 const ARM_Curve& nominal, 
						 const ARM_Curve& strikes, 
						 const ARM_Curve& leveragesShort, 
						 const ARM_Curve& leveragesLong, 
						 const ARM_Curve& cpnMin, 
						 const ARM_Curve& cpnMax, 
						 const ARM_Curve& margins,
						 const ARM_Curve& fees,
						 ARM_result& result,
						 long& objId)
{
	char sStartDate[11];
	char sEndDate[11];
	Local_XLDATE2ARMDATE(dStartDate,sStartDate);
 	Local_XLDATE2ARMDATE(dEndDate,sEndDate);
	ARM_Date startDate=ARM_Date (sStartDate);
	ARM_Date endDate=ARM_Date (sEndDate);

	CallableQuantoSpreadOptionCreator cqsoCreator(startDate, endDate); 

	//Indices information
	ARM_Currency cpnCcy(forCurrency.c_str()); 
	ARM_Currency fundCcy(domCurrency.c_str()); 
	long liborType1Id = ARM_ConvIrType (index1);
	long liborType2Id = ARM_ConvIrType (index2);
	ARM_INDEX_TYPE shortIndex = ARM_INDEX_TYPE(liborType2Id); 
	ARM_INDEX_TYPE longIndex = ARM_INDEX_TYPE(liborType1Id); 
	ARM_INDEX_TYPE fundIndex = EURIBOR1Y;

	cqsoCreator.setUnderlyings(cpnCcy, shortIndex,longIndex, fundCcy,fundIndex);

	// Schedule
	ScheduleInformation schedule;
	schedule.frequency   = iFrequency; 
	schedule.isAdjusted  = iIsAdjusted; 
	schedule.dayCounter  = iDayCounter; 
	schedule.resetType	 = iResetType; 
	schedule.resetGap	 = iResetGap;
	schedule.resetCalendar   = sResetCalendar;
	schedule.paymentCalendar = sPaymtCalendar;

	cqsoCreator.setScheduleInformation(schedule);

	// Payoff 
	cqsoCreator.setPayoffInformation(nominal,strikes,leveragesShort, 
								leveragesLong, cpnMin,cpnMax,margins);

	int exerciseFreq=1, noticeGap=2; 
	cqsoCreator.setCallabilityInformation(exerciseFreq, noticeGap, fees);

	ARM_GenSecurityPtr security = cqsoCreator.create();  

	ARM_GenSecurity* newCqso = (ARM_GenSecurity*) security->Clone();

	// Memory utilisation: why do we need to complicate everywhere it's possible?
	long cqsoId; 
	if (objId == -1)
	{
		CREATE_GLOBAL_OBJECT();
		cqsoId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCqso);
		if (cqsoId == RET_KO)
		{
			if (newCqso)
				delete newCqso;
			newCqso=0;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return ARM_KO;
		}
		result.setLong(cqsoId);
		objId=cqsoId; 
	}
	else 
	{
		ARM_GenSecurity * cqso = (ARM_GenSecurity *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cqso, ARM_GENSECURITY) == 1)
		{
			if (cqso)
			{
				delete cqso;
				cqso = 0;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newCqso, objId);
			return ARM_OK;
		}
		else
		{
			if (newCqso) {
				delete newCqso;
				newCqso = 0;
			}
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}

	}

	return ARM_OK;
}
