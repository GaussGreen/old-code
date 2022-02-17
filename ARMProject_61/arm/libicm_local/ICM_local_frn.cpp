
#pragma warning(disable :4786 4541 4800 4250)

#include "firstToBeIncluded.h"
#include "CCdate.h"
#include "CCstring.h"

#include <math.h>

#include <ICMKernel\glob\icm_enums.h>
#include <ARMKernel\inst\swap.h>
#include <ICMKernel\inst\icm_frn.h>
#include <ICMKernel\inst\icm_cln.h>
#include <ARMKernel\ccy\currency.h>
#include <ARMKernel\mod\model.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>



long ICMLOCAL_FRN   (double	Spread,
					 double	Int_Accrual_DateIn,
					 double	MaturityIn,
					 long IrIndex,
					 double InitialRate,
					 double LastIndexFixing,
					 qPAYMENT_PREMIUM_LEG AccOnDef,
					 double	NotionalAmountIn ,
					 int		DayCount,
					 int	AccruedDayCount,
					 int		SettlementGap,	
					 CCString	DiscCurrency,
					 CCString   ResetCalendar,
					 CCString   PayCalendar,
					 double	First_Period_Reference_DateIn,
					 ARM_result& result,
					 long objId)
{
	long frnId;

	ICM_Frn* frn = NULL;
	ICM_Frn* newfrn = NULL;
	
	long resetFreq = -1;
	long payFreq = -1;
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pStartDatee=new char[11];
	char* pMaturityIn=new char[11];

	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(Int_Accrual_DateIn,pStartDatee);
		Local_XLDATE2ARMDATE(MaturityIn,pMaturityIn);

		ARM_Date refDate; 
		if (First_Period_Reference_DateIn == -1.0)
		{	
			// strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		Local_XLDATE2ARMDATE(First_Period_Reference_DateIn,refDate);

		ARM_IRIndex* IRI = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(IrIndex);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IRI, ARM_IRINDEX) == 0)
		{

			if (pStartDatee)
				delete [] pStartDatee;
			pStartDatee = NULL;

			if (pMaturityIn)
				delete [] pMaturityIn;
			pMaturityIn = NULL;

			result.setMsg ("ARM_ERR: IR Index is not of a good type");
			return ARM_KO;
		}

		
		newfrn = new ICM_Frn((Spread/10000.),
							(ARM_Date)pStartDatee,
							(ARM_Date)pMaturityIn, 
							First_Period_Reference_DateIn == -1.0 ? 0 : &refDate,
							IRI,
							(InitialRate/100.) ,
							LastIndexFixing,
							 AccOnDef ,
							NotionalAmountIn , 
							DayCount , 
							AccruedDayCount, 
							SettlementGap , 
							stubrule,
							CCSTringToSTLString(DiscCurrency) ,
							CCSTringToSTLString(ResetCalendar) , 
							CCSTringToSTLString(PayCalendar),
							K_NX_END); 
							// pRefDate);

		


			if (pStartDatee)
				delete [] pStartDatee;
			pStartDatee = NULL;

			if (pMaturityIn)
				delete [] pMaturityIn;
			pMaturityIn = NULL;



		if (newfrn == NULL)
		{
			result.setMsg ("ARM_ERR: FRN is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			frnId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newfrn);

			if (frnId == RET_KO)
			{
				if (newfrn)
					delete newfrn;
				newfrn = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(frnId);

			return ARM_OK;
		}
		else
		{
			frn = (ICM_Frn *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(frn, ICM_FRN) == 1)
			{
				if (frn)
				{
					delete frn;
					frn = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newfrn, objId);

				return ARM_OK;
			}
			else
			{
				if (newfrn)
					delete newfrn;
				newfrn = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();


			if (pStartDatee)
				delete [] pStartDatee;
			pStartDatee = NULL;

			if (pMaturityIn)
				delete [] pMaturityIn;
			pMaturityIn = NULL;



		ARM_RESULT();
	}
}


long ICMLOCAL_CLN(const double& startdate,
				  const double& enddate,
				  const double& refdate,
				  const double& fstcpnrefdate,
				  const long&	IrIndex,
				  const double& spread,
				  const double&	Notional,
				  qPAYMENT_PREMIUM_LEG	AccOnDef,
				  const int&	DayCount,
				  const int&	decompfreq,
				  const int&	stubrule,
				  const int&	resetgap,
				  CCString	DiscCurrency,
				  CCString   ResetCalendar,
				  CCString   PayCalendar,
				  const int&	nxchange,
				  const bool&   includematurity,
				  const int&	adjstartdate,
				  const int&	LegType,
				  const double& binary,
				  ARM_result& result,
				  long objId)
{
	long clnId;

	ICM_Cln* cln = NULL;
	ICM_Cln* newcln = NULL;
	// ARM_Currency* ccy = NULL;
	
	int p_stubrule = stubrule;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{	result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;	}
	
	char p_startdate[11];
	char p_enddate[11];
	char p_refdate[11];
	char p_fstcpnrefdate[11];

	CCString msg ("");

	//char* payCal = NULL;
	//char* ResetCal = NULL;

	try
	{
		// Creation des ccy
		//char* tmp = (char*)DiscCurrency;
		//if (tmp)
		//{	ccy = new ARM_Currency(tmp); }
		//tmp = NULL;

		//if (strcmp((const char*)ResetCalendar,"NULL"))
		//{
		//	ResetCal = new char[(ARM_NB_MAX_CAL*3)+1];
		//	strcpy(ResetCal,ResetCalendar);
		//}

		//if (strcmp((const char*)PayCalendar,"NULL"))
		//{
		//	payCal = new char[(ARM_NB_MAX_CAL*3)+1];
		//	strcpy(payCal,PayCalendar);
		//}

		Local_XLDATE2ARMDATE(startdate,p_startdate);
		Local_XLDATE2ARMDATE(enddate,p_enddate);

		ARM_Date D_refdate;
		if (refdate == -1.0)
		{	p_stubrule = 1;} // stub rule is shortstart 
		else 
		{	Local_XLDATE2ARMDATE(refdate,p_refdate);
			D_refdate = (ARM_Date)p_refdate;}

		ARM_Date D_fstcpnrefdate;
		if (fstcpnrefdate == -1.0)
		{	p_stubrule = 1;} // stub rule is shortstart 
		else 
		{	Local_XLDATE2ARMDATE(fstcpnrefdate,p_fstcpnrefdate);
			D_fstcpnrefdate = (ARM_Date)p_fstcpnrefdate;}

		ARM_IRIndex* IRIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(IrIndex);

		if (IRIndex)
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(IRIndex, ARM_IRINDEX) == 0)
		{
			//if (ResetCal)
			//	delete [] ResetCal;
			//ResetCal = NULL;

			//if (payCal)
			//	delete [] payCal;
			//payCal = NULL;

			result.setMsg ("ARM_ERR: IR Index is not of a good type");
			return ARM_KO;
		}

		newcln = new ICM_Cln((ARM_Date)p_startdate, 
							(ARM_Date)p_enddate,
							refdate==-1 ? 0 : &D_refdate,
							fstcpnrefdate==-1 ? 0 : &D_fstcpnrefdate,
							IRIndex,
							(spread/1.e4),
							Notional,
							(PAYMENT_PREMIUM_LEG)AccOnDef,
							DayCount,
							K_PAY, 
							DayCount, 
							decompfreq,
							p_stubrule,
							resetgap,
							CCSTringToSTLString(ResetCalendar), 
							CCSTringToSTLString(DiscCurrency),
							CCSTringToSTLString(PayCalendar),
							nxchange,
							includematurity,
							adjstartdate,
							(qCredit_Leg_Type)LegType,
							binary,
							ISSUER_UNDEFINE);

		

			//if (ResetCal)
			//	delete [] ResetCal;
			//ResetCal = NULL;

			//if (payCal)
			//	delete [] payCal;
			//payCal = NULL;


		if (newcln == NULL)
		{
			result.setMsg ("ARM_ERR: CLN is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			clnId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcln);

			if (clnId == RET_KO)
			{
				if (newcln)
					delete newcln;
				newcln = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(clnId);

			return ARM_OK;
		}
		else
		{
			cln = (ICM_Cln *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cln, ICM_CLN) == 1)
			{
				if (cln)
				{
					delete cln;
					cln = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcln, objId);

				return ARM_OK;
			}
			else
			{
				if (newcln)
					delete newcln;
				newcln = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		//if (ResetCal)
		//	delete [] ResetCal;
		//ResetCal = NULL;

		//if (payCal)
		//	delete [] payCal;
		//payCal = NULL;

		ARM_RESULT();
	}
}


/*---- End Of File ----*/

// EOF %M%