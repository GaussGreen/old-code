#include "firstToBeIncluded.h"
#include "CCdate.h"
#include "CCstring.h"

#include <ICMKernel\glob\icm_enums.h>
#include <ICMKernel\inst\icm_cds.h>
#include <ICMKernel\inst\icm_credit_index.h>
#include <ICMKernel\inst\icm_corridorleg.h>
#include <ICMKernel\pricer\icm_pricer.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <ARM\libicm_local\icm_local_leg.h>
#include <ICMKernel\util\icm_schedule_info.h>
#include <ARMKernel\crv\volflat.h>
#include <ICMKernel\crv\icm_defaultcurve.h>


long ICMLOCAL_FIXEDLEG (double		startDateIn, 
					    double		endDateIn,
						double		fixedRateIn,
						qPAYMENT_PREMIUM_LEG			AccruedOnDefaultIn,
						int			AccruedDayCountIn,
						double		LastIndexFixingIn,
						int			rcvOrPayIn,
						int			freqIn,
						int			dayCountIn,
						int			decompFreqIn,
						int			payTimingIn,
						int			intRuleIn,
						int			stubRuleIn,
						long		discountCcyIn,
						CCString	payCalNameIn,
						int			nxChangeIn,
						double		refDateIn,
						ARM_result& result,
						long		objId)
{
	long LegId;

	ICM_Leg* leg = NULL;
	ICM_Leg* newLeg = NULL;
	
	ARM_Currency* discountCcy = NULL;
	
	long resetFreq = -1;
	long payFreq = -1;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* startDate=new char[11];
	char* endDate=new char[11];

	CCString msg ("");

	try
	{
		// Creation des ccy
		if ( discountCcyIn == ARM_NULL_OBJECT )
		{
			discountCcy = ARM_DEFAULT_CURRENCY;
		}
		else
		{
			discountCcy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyIn);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(discountCcy, ARM_CURRENCY) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (endDate)
					delete [] endDate;
				endDate = NULL;


				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}


		ARM_Date refDate ; 
		if (refDateIn!=-1) Local_XLDATE2ARMDATE(refDateIn,refDate); 


		Local_XLDATE2ARMDATE(startDateIn,startDate);
		Local_XLDATE2ARMDATE(endDateIn,endDate);

		newLeg = new ICM_Leg((ARM_Date)startDate, 
							 (ARM_Date)endDate, 
							 refDateIn==-1 ? 0 : &refDate,
							 0,
							 (fixedRateIn/100.),
							 AccruedOnDefaultIn,
							 AccruedDayCountIn,
							 LastIndexFixingIn,
							 rcvOrPayIn , 
							 freqIn , 
							 dayCountIn , 
							 decompFreqIn ,
							 payTimingIn ,
							 intRuleIn ,
							 stubRuleIn ,
							 discountCcy->GetCcyName(),
							 CCSTringToSTLString(payCalNameIn),
							 nxChangeIn,
								EXCLUDE_MATURITY,//const bool& includematurity /*= EXCLUDE_MATURITY*/,
							K_ADJUSTED,// const int& adjStartDate /*= K_ADJUSTED*/,
							qRunning_Leg, // const qCredit_Leg_Type& LegType /*= qRunning_Leg*/,
							CREDIT_DEFAULT_VALUE,// const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
							ISSUER_UNDEFINE // const string& name /*= ISSUER_UNDEFINE*/) 
							 ) ; 
							 // refDate);


		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;


		if (newLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			LegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg);

			if (LegId == RET_KO)
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(LegId);

			return ARM_OK;
		}
		else
		{
			leg = (ICM_Leg*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(leg, ICM_LEG) == 1)
			{
				if (leg)
				{
					delete leg;
					leg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg, objId);

				return ARM_OK;
			}
			else
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (discountCcy)
			delete discountCcy;
		discountCcy = NULL;


		if (newLeg)
			delete newLeg;
		newLeg = NULL;

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_FIXEDLEGGEN (double		startDateIn, 
							double		endDateIn,
							double		fixedRateIn,
							int			AccruedDayCountIn,
							int			freqIn,
							int			dayCountIn,
							int			payTimingIn,
							int			intRuleIn,
							int			stubRuleIn,
							long		discountCcyIn,
							CCString	payCalNameIn,
							double		refDateIn,
							ARM_result& result,
							long		objId)
{
	long LegId;

	ICM_Leg* leg = NULL;
	ICM_Leg* newLeg = NULL;
	
	ARM_Currency* discountCcy = NULL;
	
	long resetFreq = -1;
	long payFreq = -1;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* startDate=new char[11];
	char* endDate=new char[11];
	// char* refDate=new char[11];

	CCString msg ("");

	try
	{
		// Creation des ccy
		if ( discountCcyIn == ARM_NULL_OBJECT )
		{
			discountCcy = ARM_DEFAULT_CURRENCY;
		}
		else
		{
			discountCcy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(discountCcyIn);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(discountCcy, ARM_CURRENCY) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (endDate)
					delete [] endDate;
				endDate = NULL;


				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}

		ARM_Date refDate ;
		if (refDateIn!=-1) Local_XLDATE2ARMDATE(refDateIn,refDate);


		Local_XLDATE2ARMDATE(startDateIn,startDate);
		Local_XLDATE2ARMDATE(endDateIn,endDate);


		newLeg = new ICM_Leg((ARM_Date)startDate, 
							 (ARM_Date)endDate, 
							 refDateIn==-1 ? 0 : &refDate,
							 0,
							 (fixedRateIn/100.),
							 AccruedDayCountIn,
							 freqIn , 
							 dayCountIn , 
							 payTimingIn ,
							 intRuleIn ,
							 stubRuleIn ,
							 discountCcy->GetCcyName(),
							 CCSTringToSTLString(payCalNameIn),
							EXCLUDE_MATURITY, // const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
							K_ADJUSTED, // const int& adjStartDate /* = K_ADJUSTED*/ ,
							qRunning_Leg, // const qCredit_Leg_Type& LegType/* LegType = qRunning_Leg*/ ,
							CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
							ISSUER_UNDEFINE // const string& name /* = ISSUER_UNDEFINE*/ ) 
							 ) ;


		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;


		if (newLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			LegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg);

			if (LegId == RET_KO)
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(LegId);

			return ARM_OK;
		}
		else
		{
			leg = (ICM_Leg*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(leg, ICM_LEG) == 1)
			{
				if (leg)
				{
					delete leg;
					leg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg, objId);

				return ARM_OK;
			}
			else
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (discountCcy)
			delete discountCcy;
		discountCcy = NULL;


		if (newLeg)
			delete newLeg;
		newLeg = NULL;

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;



		ARM_RESULT();
	}
}


long ICMLOCAL_SetVariableSpread (long secId,
								long rId,
								ARM_result& result)
{
	ICM_Leg* sec = NULL;
	ARM_ReferenceValue* ref = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(sec, ICM_LEG) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		ref = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(rId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ref, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		sec->SetVariableSpread(ref);

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


long ICMLOCAL_GenLeg (const double&	startDateIn, 
					  const double&	endDateIn,
					  const double&	fixedRate,
					  const double&	fixedNotional,
					  const long&	VarNotId,
					  const long&	VarRateId,
					  const long&	ExchangeNotId,
 					  const int  &	frequency,
					  const int &	daycount,
					  const int	&	payTiming,
					  const int &	intRule,
					  const int	&	stubRule,
					  const long &	CcyId,
					  const CCString& payCal,
					  const double&	RefDate,
					  const bool&	includematurity,
					  const int&	adjuststartDate,
					  const int&	legtype,
					  const int&	IndexId,
					  const int&	creditlag,
					  const double& binary,
					  const CCString& name,
					  const int&	nxchange,
					  qPAYMENT_PREMIUM_LEG	accruedOnDef,
					  ARM_result& result,
					  long	objId)
{
	long LegId;

	ICM_Leg* leg = NULL;
	ICM_Leg* newLeg = NULL;
	
	ARM_Currency* discountCcy = NULL;
	// char* payCalName = NULL;

	ARM_ReferenceValue* Notional = NULL;
	ARM_ReferenceValue* Rates = NULL;
	ARM_ReferenceValue* Exchange = NULL;
	ARM_IRIndex* pIndex = NULL;
	
	long Frequency = -1;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* startDate=new char[11];
	char* endDate=new char[11];

	CCString msg ("");

	try
	{
		// Creation des ccy
		if ( CcyId == -1 )
		{
			discountCcy = ARM_DEFAULT_CURRENCY;
		}
		else
		{
			discountCcy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(CcyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(discountCcy, ARM_CURRENCY) == 0)
			{
				if (startDate)
					delete [] startDate;
				startDate = NULL;

				if (endDate)
					delete [] endDate;
				endDate = NULL;

				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}


		ARM_Date refDate ;
		if (RefDate!=-1) 	Local_XLDATE2ARMDATE(RefDate,refDate);

		Local_XLDATE2ARMDATE(startDateIn,startDate);
		Local_XLDATE2ARMDATE(endDateIn,endDate);

		if (VarNotId>=0)
		{
		Notional = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(VarNotId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Notional, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Notional is not a reference value");
			return ARM_KO;
		}
		}

		if (VarRateId>=0)
		{
		Rates = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(VarRateId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Rates, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Notional is not a reference value");
			return ARM_KO;
		}
		}

		if (ExchangeNotId>=0)
		{
		Exchange = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(ExchangeNotId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Exchange, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Exchange is not a reference value");
			return ARM_KO;
		}
		}

		pIndex = dynamic_cast<ARM_IRIndex*> (LOCAL_PERSISTENT_OBJECTS->GetObject(IndexId));


		newLeg = new ICM_Leg((ARM_Date)startDate, 
							 (ARM_Date)endDate, 
							 RefDate == -1 ? 0 : &refDate,
							 NULL,
							 (fixedRate/10000.),
							 fixedNotional,	
							 Notional,
							 Rates,
							 Exchange,
							 frequency , 
							 daycount , 
							 payTiming,
							 intRule ,
							 stubRule ,
							 discountCcy->GetCcyName(),
							 CCSTringToSTLString(payCal),
							 (qCredit_Leg_Type) legtype,
							 includematurity,
							 adjuststartDate,
							 pIndex,
							 binary,
							 CCSTringToSTLString(name),
							 nxchange,
							 accruedOnDef);


		newLeg->SetCreditLag(creditlag);

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;

		// if (refDate)
		//	delete [] refDate;
		//refDate = NULL;

		if (newLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			LegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg);

			if (LegId == RET_KO)
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(LegId);

			return ARM_OK;
		}
		else
		{
			leg = (ICM_Leg*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(leg, ICM_LEG) == 1)
			{
				if (leg)
				{
					delete leg;
					leg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newLeg, objId);

				return ARM_OK;
			}
			else
			{
				if (newLeg)
					delete newLeg;
				newLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (discountCcy)
			delete discountCcy;
		discountCcy = NULL;

//		if (payCalName)
//			delete payCalName;
//		payCalName = NULL;

		if (newLeg)
			delete newLeg;
		newLeg = NULL;

		if (startDate)
			delete [] startDate;
		startDate = NULL;

		if (endDate)
			delete [] endDate;
		endDate = NULL;

		// if (refDate)
		//	delete [] refDate;
		// refDate = NULL;


		ARM_RESULT();
	}
}


long ICMLOCAL_SetLeg (long secId,
					  long LegId,
					  long option,	
					  ARM_result& result)
{
	ICM_Cds* sec = NULL;
	ICM_Leg* Leg = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: object is not of a good type");
			return ARM_KO;
		}

		Leg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(LegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Leg, ICM_LEG) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		switch (option)
		{
		case qCMPFEELEGPV : 
			{
			sec->SetFeeLeg(Leg);
			break;
			}
		case qCMPDEFLEGPV : 
			{
			sec->SetDefLeg(Leg);
			break;
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

long ICMLOCAL_Index (const string& IndexName, 
					 const vector<string>&labels,
					 int Basis,
					 int ResetFreq,
					 int PayFreq,
					 ARM_Vector yearterm,
					 ARM_Vector Spread,
					 const string& currency,
					 long Method,
					 int DefaultCurveId,
					 int fwdRule,
					 int resetTiming,
					 int resetGap,
					 int payTiming,
					 int payGap,
					 int intRule,
					 qCDS_ADJ AdjCalType,
					 int cm_resetWeekDay,
					 int cm_resetOccur,
					 long objId)
{
	//string IndexName;
	//vector<string>* labels = NULL;

	ICM_Credit_Index* index = NULL;
	long IndexId = 0;
	string ccy("");
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_ERR: Pb with accessing objects"); 
	
	//char* pRefDate=new char[11];
	ICM_DefaultCurve* pDefaultCurve = NULL;
	ARM_VolFlat* pVolatility = NULL;
	

	// Creation des ccy
	if (currency == "NONE") ccy = ARM_DEFAULT_COUNTRY/*ARM_DEFAULT_CURRENCY*/;
	else ccy = currency;
	
	ARM_Date refDate; 
	if (DefaultCurveId != -1)
	{
		pDefaultCurve = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DefaultCurveId);	
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pDefaultCurve, ICM_DEFAULTCURVE) == 0)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_ERR: The Default Curve is not of a good type"); 
	}

	
	index = new ICM_Credit_Index(Basis, 
										 ResetFreq, 
										 PayFreq, 
										yearterm,
										 ccy,
										 IndexName,
										 labels,
										 (qINDEX_CMPT_METHOD)Method,
										 Spread,
										 pDefaultCurve,
										 fwdRule,
										 intRule,
										 resetTiming,
										 resetGap,
										 payTiming,
										 payGap,
										 AdjCalType,
										 cm_resetWeekDay,
										 cm_resetOccur);
		
	long QId = LocalPersistent::get().adopt(index,objId );
	return QId; 
}



long ICMLOCAL_SetCoupons(long secId,
						 long couponId,
						 long styleId,
						 long PartRatesId,
						 ARM_result& result)
{
	ICM_Cds* sec = NULL;
	ARM_ReferenceValue* refcpn = NULL;
	ARM_ReferenceValue* refstyle = NULL;
	ARM_ReferenceValue* partrates = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		refcpn = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(couponId);

		if (refcpn) {
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refcpn, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}}

		refstyle = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(styleId);

		if (refstyle) {
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refstyle, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}}

		partrates = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(PartRatesId);

		if (partrates){
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(partrates, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}}

		sec->SetCoupons(refcpn,refstyle,partrates);

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


long ICMLOCAL_SetRiskyProfile(long secId,
							  long typeId,
							  long dateId,
							  ARM_result& result)
{
	ICM_Cds* sec = NULL;
	ARM_ReferenceValue* refcpn = NULL;
	ARM_ReferenceValue* refstyle = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		sec = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		refcpn = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(dateId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refcpn, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		refstyle = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(typeId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refstyle, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: refValue is not of a good type");
			return ARM_KO;
		}

		sec->GetFeeLeg()->SetRefRiskyDate(refcpn,refstyle);

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

long ICMLOCAL_CORRIDORLEG_SCHE(const string& name,
								 double Notional,
								 long recieveOrPay,
								 long RefValueSpreads,
								 long floatingIdx,
								 double leverageFloatIdx,
								 long creditIdx,
								 long refvalueKUP,
								 long refvalueKDW,
								 long scheduleId,
								 qPAYMENT_PREMIUM_LEG accondef,
								 const string& disc_ccy,
								 ARM_result& result,
								 long objId)
{
	
	ICM_CorridorLeg* createdCorrLeg=NULL;
	ARM_IRIndex* PayIdx=NULL;
	ICM_Credit_Index* RefIdx=NULL;
	ICM_Schedule_Info* pSche=NULL;
	
	ARM_ReferenceValue *BarrierDown=NULL;
	ARM_ReferenceValue *BarrierUp=NULL;
	ARM_ReferenceValue *SpreadRefvalue=NULL;

		pSche = (ICM_Schedule_Info*) LOCAL_PERSISTENT_OBJECTS->GetObject(scheduleId);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pSche, ICM_SCHEDULE_INFO) == 0)
		{
			result.setMsg ("ARM_ERR: Schedule_info is not of a good type");
			return ARM_KO;
		}
		
		RefIdx = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(creditIdx);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RefIdx, ICM_CREDIT_INDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Ref Index is not of a good type");
			return ARM_KO;
		}

		BarrierDown = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalueKDW);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierDown, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Barrier Down is not of a good type");
			return ARM_KO;
		}

		PayIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(floatingIdx);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PayIdx, ARM_IRINDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Pay Index is not of a good type");
			return ARM_KO;
		}

		BarrierUp = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalueKUP);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierUp, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Barrier Up is not of a good type");
			return ARM_KO;
		}
	
			SpreadRefvalue = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(RefValueSpreads);
			
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(SpreadRefvalue, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Spread is not of a good type");
			return ARM_KO;
		}

		createdCorrLeg = new ICM_CorridorLeg(pSche->itsEffectiveDate, 
											pSche->itsMaturityDate, 
											pSche->GetpReferenceDate(),
											pSche->GetpFirstCpnEffDate(),
											SpreadRefvalue,
											Notional,
											recieveOrPay,
											PayIdx,
											leverageFloatIdx,
											RefIdx,
											BarrierUp,
											BarrierDown,
											accondef /*= qACCRUED_SETTLED */,
											pSche->itsAccDayCount /*= KACTUAL_365*/,
											pSche->itsPayFrequency /*= K_ANNUAL*/, 
											pSche->itsResetFreq /*= K_ANNUAL*/, 
											pSche->itsDayCount/*= K30_360*/, 
											pSche->itsPayTiming /*= K_ARREARS*/,
											pSche->itsIntRule /*= K_ADJUSTED*/,
								            pSche->itsStubrule /*= K_SHORTSTART*/,
											disc_ccy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
											pSche->itsPayCalName/*= NULL*/);
		// from SwpaLeg
		createdCorrLeg->SetSingleName(name);
		long QId = LocalPersistent::get().adopt(createdCorrLeg,objId );

		return QId; 

}

long ICMLOCAL_CORRIDORLEG(		const CCString& name,
								 double startDate,
								 double endDate,
								 double refdate,
								 double fstcpneffdate,
								 long RefValueSpreads,
								 //double Notional,
								 long floatingIdx,
								 double leverageFloatIdx,
								 long creditIdx,
								 long refvalueKUP,
								 long refvalueKDW,
								 qPAYMENT_PREMIUM_LEG accondef,
								 long accdaycount,
								 long payfreq,
								 long resetfreq,
								 long daycount,
								 long paytiming,
								 long intrule,
								 long stubrule,
								 const CCString& disc_ccy,
								 const CCString& paycalname,
								 ARM_result& result,
								 long objId)
{
	long legId;

	ICM_CorridorLeg* createdCorrLeg=NULL;
	ICM_CorridorLeg* precCorrLeg=NULL;
	ARM_IRIndex* PayIdx=NULL;
	ICM_Credit_Index* RefIdx=NULL;
	
	ARM_ReferenceValue *BarrierDown=NULL;
	ARM_ReferenceValue *BarrierUp=NULL;
	ARM_ReferenceValue *SpreadRefvalue=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char sStartDate[11];
	char sEndDate[11];

	CCString msg ("");

	double Notional = 100;

	try
	{
		Local_XLDATE2ARMDATE(startDate,sStartDate);
		Local_XLDATE2ARMDATE(endDate,sEndDate);

		ARM_Date D_refDate ; 
		if (refdate!=-1) Local_XLDATE2ARMDATE(refdate,D_refDate); 

		ARM_Date D_fstcpneffdate; 
		if (fstcpneffdate!=-1) Local_XLDATE2ARMDATE(fstcpneffdate,D_fstcpneffdate); 

		PayIdx = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(floatingIdx);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(PayIdx, ARM_IRINDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Pay Index is not of a good type");
			return ARM_KO;
		}

		RefIdx = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(creditIdx);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RefIdx, ICM_CREDIT_INDEX) == 0)
		{
			result.setMsg ("ARM_ERR: Ref Index is not of a good type");
			return ARM_KO;
		}

		BarrierDown = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalueKDW);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierDown, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Barrier Down is not of a good type");
			return ARM_KO;
		}


		BarrierUp = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(refvalueKUP);
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(BarrierUp, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Barrier Up is not of a good type");
			return ARM_KO;
		}

		SpreadRefvalue = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(RefValueSpreads);
			
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(SpreadRefvalue, ARM_REFERENCE_VALUE) == 0)
		{
			result.setMsg ("ARM_ERR: Spread is not of a good type");
			return ARM_KO;
		}

		// qPAYMENT_PREMIUM_LEG AccOnDefault = (qPAYMENT_PREMIUM_LEG) accondef;

		createdCorrLeg = new ICM_CorridorLeg((ARM_Date)sStartDate, 
											(ARM_Date)sEndDate, 
											refdate==-1 ? 0 : &D_refDate,
											fstcpneffdate==-1 ? 0 : &D_fstcpneffdate,
											SpreadRefvalue,
											//Notional,
											PayIdx,
											leverageFloatIdx,
											RefIdx,
											BarrierUp,
											BarrierDown,
											accondef /*= qACCRUED_SETTLED */,
											accdaycount /*= KACTUAL_365*/,
											payfreq /*= K_ANNUAL*/, 
											resetfreq /*= K_ANNUAL*/, 
											daycount/*= K30_360*/, 
											paytiming /*= K_ARREARS*/,
											intrule /*= K_ADJUSTED*/,
								            stubrule /*= K_SHORTSTART*/,
											CCSTringToSTLString(disc_ccy) /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
											CCSTringToSTLString(paycalname) /*= NULL*/);

		createdCorrLeg->SetSingleName(CCSTringToSTLString(name));
		if (createdCorrLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Corridor Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			legId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCorrLeg);

			if (legId == RET_KO)
			{
				if (createdCorrLeg)
					delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(legId);

			return ARM_OK;
		}
		else
		{
			precCorrLeg = (ICM_CorridorLeg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(precCorrLeg, ICM_CORRIDORLEG) == 1)
			{
				if (precCorrLeg)
				{
					delete precCorrLeg;
					precCorrLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdCorrLeg, objId);

				return ARM_OK;
			}

			else
			{
				if (createdCorrLeg)
					delete createdCorrLeg;
				createdCorrLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdCorrLeg)
			delete createdCorrLeg;
		createdCorrLeg = NULL;

		ARM_RESULT();
	}

	catch(...)
	{
		if (createdCorrLeg)
			delete createdCorrLeg;
		createdCorrLeg = NULL;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}


long ICMLOCAL_IRLEGTOCREDITLEG(int SwapLegId,
							   int LegType,
							   int creditindexId,
							   int PricerId,
							   ARM_result& result,
							   long	objId)
{
	long LegId;

	ICM_Leg* pLeg = NULL;
	ICM_Leg* pPrevLeg = NULL;
	ARM_SwapLeg* pSwapLeg = NULL;
	ICM_Credit_Index* pIndex = NULL;
	ICM_Pricer* pPricer = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		pSwapLeg = (ARM_SwapLeg*) LOCAL_PERSISTENT_OBJECTS->GetObject(SwapLegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pSwapLeg, ARM_SWAPLEG) == 0) 
		{
			if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pSwapLeg, ARM_INFLEG) == 0) &&  
				(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pSwapLeg, ARM_CMSLEG) == 0))
			{
			result.setMsg ("ARM_ERR: SwapLeg is not of a good type");
			return ARM_KO;
			}
		}

		pIndex = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(creditindexId);

		if (pIndex)
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0) 
		{
			result.setMsg ("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		}

		pPricer = (ICM_Pricer*) LOCAL_PERSISTENT_OBJECTS->GetObject(PricerId);
		if (pPricer)
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pPricer, ICM_PRICER) == 0) 
		{
			result.setMsg ("ARM_ERR: Pricer is not of a good type");
			return ARM_KO;
		}

		pLeg = new ICM_Leg(*pSwapLeg,(qCredit_Leg_Type)LegType,pIndex);
		if (pPricer) pLeg->SetRatesPricer(pPricer);

		if (pLeg == NULL)
		{
			result.setMsg ("ARM_ERR: Leg is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			LegId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pLeg);

			if (LegId == RET_KO)
			{
				if (pLeg)
					delete pLeg;
				pLeg = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(LegId);

			return ARM_OK;
		}
		else
		{
			pPrevLeg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pPrevLeg, ICM_LEG) == 1)
			{
				if (pPrevLeg)
				{
					delete pPrevLeg;
					pPrevLeg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pLeg, objId);

				return ARM_OK;
			}
			else
			{
				if (pLeg)
					delete pLeg;
				pLeg = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (pLeg)
			delete pLeg;
		pLeg = NULL;

		ARM_RESULT();
	}
}


/*---- End Of File ----*/

// EOF %M%
