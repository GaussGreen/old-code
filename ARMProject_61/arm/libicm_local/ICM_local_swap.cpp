#include "firstToBeIncluded.h"

#pragma warning(disable :4786 4541 4800 4250)

#include "CCdate.h"
#include "CCstring.h"

#include <math.h>
#include <ICMKernel\glob\icm_enums.h>
#include <ARMKernel\inst\swap.h>
#include <ICMKernel\inst\icm_cds.h>
#include <ICMKernel\inst\icm_cmcds.h>
#include <ICMKernel\inst\icm_ftd.h>
#include <ICMKernel\inst\icm_nthtd.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_cppi.h>
#include <ICMKernel\inst\icm_cdo2.h>
#include <ICMKernel\inst\icm_credit_index.h>
#include <ICMKernel\inst\icm_option.h>
#include <ICMKernel\inst\icm_spreadoption.h>
#include <ICMKernel\inst\icm_customized_cdo.h>
#include "ICMKernel\inst\icm_corridorleg.h"
#include "ICMKernel\inst\icm_gen.h"
#include "ICMKernel\inst\icm_lss_gap_option.h"
#include "ICMKernel\inst\icm_cpdo.h"
#include "ICMKernel\inst\icm_option_tranche.h"
#include <crv\volflat.h>


/// remove warnings on va_start and va_end
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
/// and do not change the order as this one has to be precisely here
/// if you do not know what you are doing, please ask someone who does!
#include <ARM\libarm_local\undef_va_vars.h>
#include <ARMKernel\ccy\currency.h>
#include <ARMKernel\mod\model.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\libicm_local\icm_local_swap.h>
#include "ICMKernel\inst\icm_collateral.h"


long ICMLOCAL_CDS   (double Spread,
					 double	EffectiveDate,
					 double MaturityDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 int	Frequency,
					 int	DayCount,
					 double FixedPayerAmount,
					 double FloatingPayerAmount,
					 int	Stubrule,
					 CCString Currency,
					 int	intRule,
					 int	CreditLag,
					 bool	IncludeMaturity,
					 double ProtectionStartDate,
					 double ProtectionEndDate,
					 CCString name,
					 double binary,
					 int	intStartAdj,
					 qPAYMENT_PREMIUM_LEG accruedOnDef,
					 long l_NotionalEch_Type, 
					 long l_NotionalEchange,
					 ARM_result& result,
					 long	objId)
{
	long cdsId;

	ICM_Cds* cds = NULL;
	ICM_Cds* newcds = NULL;
	
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	if (EffectiveDate>=MaturityDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	char* pProtectionStartDate=new char[11];
	char* pProtectionEndDate=new char[11];

	ARM_Date dProtectionStartDate;
	ARM_Date dProtectionEndDate;

	CCString msg ("");

	Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
	Local_XLDATE2ARMDATE(MaturityDate,pEndDate);
	if (ProtectionStartDate<0) 
		dProtectionStartDate = (ARM_Date)pEffectiveDate;
	else
	{
		Local_XLDATE2ARMDATE(ProtectionStartDate,pProtectionStartDate);
		dProtectionStartDate = (ARM_Date)pProtectionStartDate;
	}
	if (ProtectionEndDate<0) 
		dProtectionEndDate = (ARM_Date)pEndDate;
	else
	{
		Local_XLDATE2ARMDATE(ProtectionEndDate,pProtectionEndDate);
		dProtectionEndDate = (ARM_Date)pProtectionEndDate;
	}
	ARM_Date refDate ;
	if (ReferenceDate == -1.0)
	{	
		stubrule = 1; // stub rule is shortstart
	}
	else
	{
		Local_XLDATE2ARMDATE(ReferenceDate, refDate);
		stubrule = Stubrule;
	}		
	ARM_Date FstCpnEffDate;
	if (FirstCpnEffDate != -1.0)
		Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);

	int FrequencyDefLeg = -1;
	ARM_ReferenceValue * pRefValueNotionalExchange = NULL;
	pRefValueNotionalExchange = dynamic_cast<ARM_ReferenceValue *> (LOCAL_PERSISTENT_OBJECTS->GetObject(l_NotionalEchange));

	newcds = new ICM_Cds((ARM_Date) pEffectiveDate,
							   (ARM_Date) pEndDate,
							   ReferenceDate==-1 ? 0 : &refDate,
							   FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   dProtectionStartDate,
							   dProtectionEndDate,
							   (Spread/10000.),
							   FixedPayerAmount,FloatingPayerAmount,// 0,0,	// notionals
							   Frequency,
							   DayCount,
							   accruedOnDef,
							   CCSTringToSTLString(Currency), 
							   stubrule,
							   CreditLag,
							   FrequencyDefLeg,
							   intRule,
							   IncludeMaturity,
							   intStartAdj, //K_ADJUSTED, 
							   CCSTringToSTLString(Currency),
							   qRunning_Leg,
							   qStandart_Recovery_Leg,
							   CCSTringToSTLString(name),
							   binary,
							   l_NotionalEch_Type,
							   pRefValueNotionalExchange
							   );

 
	if (pEffectiveDate)
		delete [] pEffectiveDate;
	pEffectiveDate = NULL;

	if (pEndDate)
		delete [] pEndDate;
	pEndDate = NULL;

	if (pProtectionStartDate)
		delete [] pProtectionStartDate;
	pProtectionStartDate = NULL;

	if (pProtectionEndDate)
		delete [] pProtectionEndDate;
	pProtectionEndDate = NULL;

	long QId = LocalPersistent::get().adopt(newcds,objId );

	return QId; 

}
// ------------------------------------------------------------
// CMCDS
// ------------------------------------------------------------
long ICMLOCAL_CMCDS (double PartRate,
					 double	EffectiveDate,
					 double MaturityDate,
					 double FirstCpnEffDate,
					 int	IndexId,
					 double ReferenceDate,
					 int	Frequency,
					 int	DayCount,
					 double FixedPayerAmount,
					 double FloatingPayerAmount,
					 int	Stubrule,
					 CCString Currency,
					 int	intRule,
					 int	CreditLag,
					 bool	IncludeMaturity,
					 double ProtectionStartDate,
					 double ProtectionEndDate,
					 ARM_result& result,
					 long	objId)
{
	long cdsId;

	ICM_Cds* cds = NULL;
	ICM_Cds* newcds = NULL;
	// ARM_Currency* ccy = NULL;
	ICM_Credit_Index* pIndex = NULL;
	qSecurity_TYPE typecds = qCM_CDS;
	int adjStartDate = 1;
	// char* payCalName = NULL;
	
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}


	if (EffectiveDate>=MaturityDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		


	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];
	char* pProtectionStartDate=new char[11];
	char* pProtectionEndDate=new char[11];

	ARM_Date dProtectionStartDate;
	ARM_Date dProtectionEndDate;

	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(MaturityDate,pEndDate);

		if (ProtectionStartDate<0) 
			dProtectionStartDate = (ARM_Date)pEffectiveDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionStartDate,pProtectionStartDate);
			dProtectionStartDate = (ARM_Date)pProtectionStartDate;
		}

		if (ProtectionEndDate<0) 
			dProtectionEndDate = (ARM_Date)pEndDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionEndDate,pProtectionEndDate);
			dProtectionEndDate = (ARM_Date)pProtectionEndDate;
		}

		ARM_Date refDate; 
		if (ReferenceDate == -1.0)
		{	
			stubrule = 1; // stub rule is shortstart
		}
		else
		{
			Local_XLDATE2ARMDATE(ReferenceDate,refDate);
			stubrule = Stubrule;
		}

		ARM_Date FstCpnEffDate;
		if (FirstCpnEffDate != -1.0)
			Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);


		int FrequencyDefLeg = -1;

		pIndex = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0) 
		{
			result.setMsg ("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		}


		// On utilise le constructeur de CDS prenant en parametre un Index.
		newcds = new ICM_Cds((ARM_Date) pEffectiveDate,
							   (ARM_Date) pEndDate,
							   ReferenceDate==-1 ? 0 : &refDate,
							   FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   dProtectionStartDate,
							   dProtectionEndDate,
							   PartRate,
							   FixedPayerAmount,FloatingPayerAmount,0,0,
							   pIndex,
							   Frequency,
							   DayCount,
							   qACCRUED_SETTLED,
							   CCSTringToSTLString(Currency), 
							   stubrule,
							   CreditLag,
							   FrequencyDefLeg,
							   intRule,
							   IncludeMaturity,
							   adjStartDate,
							   std::string(), // payCalName ,
							   typecds,
								ISSUER_UNDEFINE, // const string& name  ,
								CREDIT_DEFAULT_VALUE // const double& Binary  );
							   );


		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pProtectionStartDate)
			delete [] pProtectionStartDate;
		pProtectionStartDate = NULL;

		if (pProtectionEndDate)
			delete [] pProtectionEndDate;
		pProtectionEndDate = NULL;

		if (newcds == NULL)
		{
			result.setMsg ("ARM_ERR: CMCDS is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds);

			if (cdsId == RET_KO)
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdsId);

			return ARM_OK;
		}
		else
		{
			cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 1)
			{
				if (cds)
				{
					delete cds;
					cds = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds, objId);

				return ARM_OK;
			}
			else
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newcds)
			delete newcds;
		newcds = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		ARM_RESULT();
	}
}



long ICMLOCAL_CAPFLOORCMCDS (double PartRate,
							 double	EffectiveDate,
							 double MaturityDate,
							 double CapLevel,
							 double FloorLevel,
							 int	IndexId,
							 double ReferenceDate,
							 int	Frequency,
							 int	DayCount,
							 double FixedPayerAmount,
							 double FloatingPayerAmount,
							 int	Stubrule,
							 CCString Currency,
							 int	intRule,
							 int	CreditLag,
							 bool	IncludeMaturity,
							 double ProtectionStartDate,
							 double ProtectionEndDate,
							 ARM_result& result,
							 long objId )
{
	long cmcdsId;

	ICM_Cmcds* cmcds = NULL;
	ICM_Cmcds* newcmcds = NULL;
	// ARM_Currency* ccy = NULL;
	ICM_Credit_Index* pIndex = NULL;
	qSecurity_TYPE typecds = qCM_CDS;
	int adjStartDate = 1;
	// char* payCalName = NULL;
	
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}


	if (EffectiveDate>=MaturityDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		


	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];
	char* pProtectionStartDate=new char[11];
	char* pProtectionEndDate=new char[11];

	ARM_Date dProtectionStartDate;
	ARM_Date dProtectionEndDate;

	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(MaturityDate,pEndDate);

		if (ProtectionStartDate<0) 
			dProtectionStartDate = (ARM_Date)pEffectiveDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionStartDate,pProtectionStartDate);
			dProtectionStartDate = (ARM_Date)pProtectionStartDate;
		}

		if (ProtectionEndDate<0) 
			dProtectionEndDate = (ARM_Date)pEndDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionEndDate,pProtectionEndDate);
			dProtectionEndDate = (ARM_Date)pProtectionEndDate;
		}

		ARM_Date refDate ;
		if (ReferenceDate == -1.0)
		{	
			// strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		{
			Local_XLDATE2ARMDATE(ReferenceDate,refDate);
			stubrule = Stubrule;
		}

		int FrequencyDefLeg = -1;

		pIndex = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0) 
		{
			result.setMsg ("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		}


		// On utilise le constructeur de CMCDS

		newcmcds = new ICM_Cmcds((ARM_Date) pEffectiveDate,
								(ARM_Date) pEndDate,
								ReferenceDate==-1 ? 0 : &refDate,
								0,
								dProtectionStartDate,
								dProtectionEndDate,
								CapLevel,
								FloorLevel,
								PartRate,
								pIndex,
								// pRefDate,
								Frequency,
								DayCount,
								FixedPayerAmount, 
								qACCRUED_SETTLED,
								CCSTringToSTLString(Currency), 
								FloatingPayerAmount,
								stubrule,
								CreditLag,
								FrequencyDefLeg,
								intRule,
								IncludeMaturity,
								adjStartDate,
								std::string(), // payCalName,
								typecds
								);


		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pProtectionStartDate)
			delete [] pProtectionStartDate;
		pProtectionStartDate = NULL;

		if (pProtectionEndDate)
			delete [] pProtectionEndDate;
		pProtectionEndDate = NULL;

		if (newcmcds == NULL)
		{
			result.setMsg ("ARM_ERR: CAPFLOORCMCDS is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cmcdsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcmcds);

			if (cmcdsId == RET_KO)
			{
				if (newcmcds)
					delete newcmcds;
				newcmcds = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cmcdsId);

			return ARM_OK;
		}
		else
		{
			cmcds = (ICM_Cmcds *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cmcds, ICM_CAPFLOORCMCDS) == 1)
			{
				if (cmcds)
				{
					delete cmcds;
					cmcds = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcmcds, objId);

				return ARM_OK;
			}
			else
			{
				if (newcmcds)
					delete newcmcds;
				newcmcds = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();


		if (newcmcds)
			delete newcmcds;
		newcmcds = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		ARM_RESULT();
	}
}

// ------------------------------------------------------------

long ICMLOCAL_CDSGEN(int FeelegId,
					 int DeflegId,
					 double RcvFee,
					 double TradedNotional,
					 ARM_result& result,
					 long	objId)
{
	long cdsId;

	ICM_Cds* cds = NULL;
	ICM_Cds* newcds = NULL;
	ICM_Leg* pFeeleg = NULL;
	ICM_Leg* pDefleg = NULL;
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		pFeeleg = dynamic_cast<ICM_Leg *>(LOCAL_PERSISTENT_OBJECTS->GetObject(FeelegId));

		if (	( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pFeeleg, ICM_LEG) == 0 )
			&&  (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pFeeleg, ICM_CORRIDORLEG) == 0) 
			) 
		{
			result.setMsg ("ARM_ERR: FeeLeg is not of a good type");
			return ARM_KO;
		}


		pDefleg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(DeflegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pDefleg, ICM_LEG) == 0) 
		{
			result.setMsg ("ARM_ERR: Defleg is not of a good type");
			return ARM_KO;
		}

		newcds = new ICM_Cds(pFeeleg,pDefleg);
		newcds->SetPorS(RcvFee);
		double InitNot=0.;
		if ((pDefleg->GetCreditLegType() == qNone_Leg) || (pFeeleg->GetCreditInfosRef().GetNotionals().Elt(0)))
		{
			InitNot=fabs(pFeeleg->GetCreditInfosRef().GetNotionals().Elt(0));
			newcds->SetTradedCoef(TradedNotional/InitNot);
		}
		else
		if ((pFeeleg->GetCreditLegType() == qNone_Leg) || (pDefleg->GetCreditInfosRef().GetNotionals().Elt(0)))
		{	
			InitNot=fabs(pDefleg->GetCreditInfosRef().GetNotionals().Elt(0));
			newcds->SetTradedCoef(TradedNotional/InitNot);
		}
		else
			newcds->SetTradedCoef(1.);

		if (newcds == NULL)
		{
			result.setMsg ("ARM_ERR: CDS is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds);

			if (cdsId == RET_KO)
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdsId);

			return ARM_OK;
		}
		else
		{
			cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 1)
			{
				if (cds)
				{
					delete cds;
					cds = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds, objId);

				return ARM_OK;
			}
			else
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newcds)
			delete newcds;
		newcds = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_FTD   (double Spread,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 VECTOR<CCString> labels,
					 int	Frequency,
					 int	DayCount,
					 double IssuerNotional,
					 qPAYMENT_PREMIUM_LEG	AccruedOnDefault,
					 CCString Currency,
					 int CreditLag,
					 int stub,
					 int intRule,
					 int startAdj,
					 ARM_result& result,
					 long objId)
{
	long ftdId;

	ICM_Ftd* ftd = NULL;
	ICM_Ftd* newftd = NULL;

	int il=0;
	// int nbissuers = 0;
	//char** plabels = NULL;
	
	int stubrule = stub;// stub rule is shortend

	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];

	char* pPayCreditLag=new char[11];

	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate ;
		if (ReferenceDate != -1) Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		// nbissuers = labels.size();
		vector<string> labels_(labels.size()); 
		for(il=0;il<labels_.size() ;il++) 
			labels_[il]=labels[il]; 


		ARM_Date FstCpnEffDate;
		if (FirstCpnEffDate != -1.0)
			Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);

		newftd = new ICM_Ftd((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
							   FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   (Spread/10000.),
							   	intRule, // K_ADJUSTED,  intRule
								startAdj, //	K_ADJUSTED, adjStartDate	
							   // nbissuers,
							   labels_,	
							   Frequency,
							   DayCount,
							   IssuerNotional, 
							   AccruedOnDefault,
							   CCSTringToSTLString(Currency), 
							   IssuerNotional,
							   stubrule,
							DEFAULT_CREDIT_LAG, // const double& CreditLag  ,
							DEFAULT_FRQ_DEFLEG, // const int& FreqDefLeg  ,
							CREDIT_DEFAULT_VALUE, // const double& Binary  ,
							std::string(), // const std::string& payCalName  ,
							INCLUDE_MATURITY // const bool& IncludeMaturity  
							   );

		// if (plabels)
		// {
		// 	for (il=0; il<nbissuers; il++)
		// 		delete plabels[il];
		// 	delete[] plabels;
		// }


		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;


		if (newftd == NULL)
		{
			result.setMsg ("ARM_ERR: FTD is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			ftdId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newftd);

			if (ftdId == RET_KO)
			{
				if (newftd)
					delete newftd;
				newftd = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(ftdId);

			return ARM_OK;
		}
		else
		{
			ftd = (ICM_Ftd *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ftd, ICM_FTD) == 1)
			{
				if (ftd)
				{
					delete ftd;
					ftd = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newftd, objId);

				return ARM_OK;
			}
			else
			{
				if (newftd)
					delete newftd;
				newftd = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		// if (plabels)
		// {
		// 	for (il=0; il<nbissuers; il++)
		// 		delete plabels[il];
		// 	delete[] plabels;
		// }


		if (newftd)
			delete newftd;
		newftd = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_NTHTD   (double Spread,
					 double		EffectiveDate,
					 double		EndDate,
					 double		ReferenceDate,
					 double		FirstCpnEffDate,
					 int		FirstNumDefault,
					 int		LastNumDefault,
					 VECTOR<CCString> labels,
					 int		Frequency,
					 int		DayCount,
					 double		IssuerNotional,
					 qPAYMENT_PREMIUM_LEG		AccruedOnDefault,
					 CCString	Currency,
					 int		CreditLag,
					 int		stub,
					 int		freqdefleg,
					 double		Binary,
					 CCString	PayCal,
					 double		RcvFee,
					 double		TradedNotional,
					 bool		IncludeMaturity,
					 int intRule,
					 int startAdj,
					 ARM_result& result,
					 long objId)
{
	long nthId;

	ICM_Nthtd* nth = NULL;
	ICM_Nthtd* newnth = NULL;
	// ARM_Currency* ccy = NULL;

	// char* payCal = NULL;

	int il=0;
	// int nbissuers = 0;
	// char** plabels = NULL;
	int stubrule = stub;// stub rule is shortend

	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];

	char* pPayCreditLag=new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate ;
		if (ReferenceDate!=-1) 	Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		// nbissuers = labels.size();

		// if (nbissuers)
		// {
		// 	plabels = new char*[nbissuers];
// 
// 				for (il=0; il<nbissuers; il++)
// 			{
// 			plabels[il] = new char[60];
// 			strcpy(plabels[il],(char*)labels[il]);
// 			}

		//}
		
		// double* pIssuersNotionals = new double[nbissuers];
		// for (int il=0; il<nbissuers; il++) pIssuersNotionals[il]=IssuerNotional;

		vector<string> labels_(labels.size()); 
		for(il=0;il<labels_.size() ;il++) 
			labels_[il]=labels[il]; 
		vector<double> notionals(labels.size(),IssuerNotional); 


		ARM_Date FstCpnEffDate;
		if (FirstCpnEffDate != -1.0)
			Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);

		
		newnth = new ICM_Nthtd((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
							  FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   (Spread/10000.),
							   intRule, //	//K_ADJUSTED, intRule
							   startAdj , //, K_ADJUSTED,	// startAdj
							   FirstNumDefault,
							   LastNumDefault,
							   // nbissuers,
							   labels_,
							   notionals,
							   Frequency,
							   DayCount,
							   RcvFee*IssuerNotional, 
							   AccruedOnDefault,
							   CCSTringToSTLString(Currency), 
							   RcvFee*IssuerNotional,
							   stubrule,
							   CreditLag,
							   freqdefleg,
							   Binary,
							   CCSTringToSTLString(PayCal),
							   IncludeMaturity
							   );



		newnth->SetTradedCoef(TradedNotional/IssuerNotional);

		// if (plabels)
		// {
		// 	for (il=0; il<nbissuers; il++)
		// 		delete plabels[il];
		// 	delete[] plabels;
		// }

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		if (newnth == NULL)
		{
			result.setMsg ("ARM_ERR: NTHTD is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			nthId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newnth);

			if (nthId == RET_KO)
			{
				if (newnth)
					delete newnth;
				newnth = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(nthId);

			return ARM_OK;
		}
		else
		{
			nth = (ICM_Nthtd *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(nth, ICM_NTD) == 1)
			{
				if (nth)
				{
					delete nth;
					nth = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newnth, objId);

				return ARM_OK;
			}
			else
			{
				if (newnth)
					delete newnth;
				newnth = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		// if (plabels)
		// {
		// 	for (il=0; il<nbissuers; il++)
		// 		delete plabels[il];
		// 	delete[] plabels;
		// }


		if (newnth)
			delete newnth;
		newnth = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		ARM_RESULT();
	}
}



long ICMLOCAL_MEZZANINE   (double Spread,
					 double		EffectiveDate,
					 double		EndDate,
					 double		ReferenceDate,
					 double		FirstCpnEffDate,
					 double		MezzAmount,
					 double		SubAmount,
					 VECTOR<CCString> labels,
					 VECTOR<double> notionals,	
					 int		FreqFeeLeg,
					 int		FreqDefLeg,
					 int		DayCount,
					 double		FixedPayerAmount,
					 double		FloatingPayerAmount,
					 qPAYMENT_PREMIUM_LEG		AccruedOnDefault,
					 CCString	Currency,
					 int		CreditLag,
					 int		stub,
					 double		Binary,
					 CCString	PayCal,
					 double		RcvFee,
					 double		TradedNotional,
					 long		TypeFeeLeg,
					 long		TypeDefLeg,
					 bool		IncludeMaturity,
					 int		intRule,
					 int		adjStartDate,	
					 ARM_result& result,
					 long objId)
{
	long mezId;

	ICM_Mez* mez = NULL;
	ICM_Mez* newmez = NULL;
	//ARM_Currency* ccy = NULL;
	//int nbissuers =0;	
	int stubrule = stub;// stub rule is shortend

	//char* payCal = NULL;

	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];
	int il=0;
	// char** plabels = NULL;
	std::vector<std::string> plabels; 
	// double* pnotionals = NULL;
	
	char* pPayCreditLag=new char[11];

	CCString msg ("");

	try
	{
		// Creation des ccy

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate ;
		if (ReferenceDate!=-1) Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		// nbissuers = notionals.size();
		plabels.resize(labels.size()); // plabels = new char*[nbissuers];
		// pnotionals = new double[nbissuers];

		for (int i=0; i<plabels.size(); i++)
		{
			plabels[i]=labels[i];
			// pnotionals[i] = notionals[i];
		}

		ARM_Date FstCpnEffDate;
		if (FirstCpnEffDate != -1.0)
			Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);

		newmez = new ICM_Mez((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
   							  FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   (Spread/10000.),
							   intRule,  //	K_ADJUSTED  intRule
							   adjStartDate, //K_ADJUSTED  adjStartDate	
							   MezzAmount,
							   SubAmount,
							   plabels,
							   notionals,
							   // nbissuers,
							   // pRefDate,
							   FreqFeeLeg,
							   DayCount,
							   RcvFee*SubAmount, 
							   AccruedOnDefault,
							   CCSTringToSTLString(Currency), 
							   RcvFee*SubAmount,
							   stubrule,
							   CreditLag,
							   FreqDefLeg,
							   Binary,
							   CCSTringToSTLString(PayCal),
							   (qCredit_Leg_Type)TypeFeeLeg,
							   (qCredit_Leg_Type)TypeDefLeg,
							   IncludeMaturity);


		newmez->SetTradedCoef(TradedNotional/SubAmount);

 
		// if (pnotionals)
		// 	delete[] pnotionals;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;
		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		if (newmez == NULL)
		{
			result.setMsg ("ARM_ERR: MEZ is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			mezId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newmez);

			if (mezId == RET_KO)
			{
				if (newmez)
					delete newmez;
				newmez = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mezId);

			return ARM_OK;
		}
		else
		{
			mez = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mez, ICM_MEZ) == 1)
			{
				if (mez)
				{
					delete mez;
					mez = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newmez, objId);

				return ARM_OK;
			}
			else
			{
				if (newmez)
					delete newmez;
				newmez = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

 
		// if (pnotionals)
		// 	delete[] pnotionals;


		if (newmez)
			delete newmez;
		newmez = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;


		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_CMTranche (double	EffectiveDate,
					 double		EndDate,
					 double		ReferenceDate,
					 double		MezzAmount,
					 double		SubAmount,
					 VECTOR<CCString> labels,
					 VECTOR<double> notionals,	
					 int		IndexId,
					 double		PaticipationRate, 	
					 int		FreqFeeLeg,
					 int		FreqDefLeg,
					 int		DayCount,
					 double		FixedPayerAmount,
					 double		FloatingPayerAmount,
					 qPAYMENT_PREMIUM_LEG		AccruedOnDefault,
					 CCString	Currency,
					 int		CreditLag,
					 int		stub,
					 double		Binary,
					 CCString	PayCal,
					 double		RcvFee,
					 double		TradedNotional,
					 double		FwdFixedDate,
					 bool		IncludeMaturity,
					 ARM_result& result,
					 long objId)
{
	long mezId;

	ICM_Mez* mez = NULL;
	ICM_Mez* newmez = NULL;
	// ARM_Currency* ccy = NULL;
	int nbissuers =0;	
	int stubrule = 3;// stub rule is shortend

	// char* payCal = NULL;

	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	//char* pRefDate=new char[11];
	// char* pFwdFixedDate=new char[11];
	int il=0;
	// char** plabels = NULL;
	std::vector<std::string> plabels; 
	// double* pnotionals = NULL;
	
	char* pPayCreditLag=new char[11];

	CCString msg ("");
	ICM_Credit_Index* pIndex = NULL;

	try
	{

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate ;
		if (ReferenceDate == -1.0)
		{	
			// strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		ARM_Date fwdFixedDate ;
		if (FwdFixedDate == -1.0) 
		{	
		// 	strcpy(pFwdFixedDate,"NULL");
		}
		else
		Local_XLDATE2ARMDATE(FwdFixedDate,fwdFixedDate);

		nbissuers = notionals.size();
		plabels.resize(nbissuers); //  = new char*[nbissuers];

		// pnotionals = new double[nbissuers];

		for (int i=0; i<nbissuers; i++)
		{
			plabels[i] =labels[i];
			// pnotionals[i] = notionals[i];
		}

		pIndex = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0) 
		{
			result.setMsg ("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		}


		newmez = new ICM_Mez((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
							  0,
							  	K_ADJUSTED,	// intRule
								K_ADJUSTED,	// adjStartDate
							   MezzAmount,
							   SubAmount,
							   plabels,
							   notionals,
							   // nbissuers,
							   pIndex,
							   PaticipationRate,
							   FreqFeeLeg,
							   DayCount,
							   RcvFee*SubAmount, 
							   AccruedOnDefault,
							   CCSTringToSTLString(Currency), 
							   RcvFee*SubAmount,
							   stubrule,
							   CreditLag,
							   FreqDefLeg,
							   Binary,
							   CCSTringToSTLString(PayCal),
							   FwdFixedDate==-1?0:&fwdFixedDate,
							   IncludeMaturity);


		newmez->SetTradedCoef(TradedNotional/SubAmount);

 		//if (pnotionals)
			// delete[] pnotionals;


		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;


		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		if (newmez == NULL)
		{
			result.setMsg ("ARM_ERR: MEZ is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			mezId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newmez);

			if (mezId == RET_KO)
			{
				if (newmez)
					delete newmez;
				newmez = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mezId);

			return ARM_OK;
		}
		else
		{
			mez = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mez, ICM_MEZ) == 1)
			{
				if (mez)
				{
					delete mez;
					mez = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newmez, objId);

				return ARM_OK;
			}
			else
			{
				if (newmez)
					delete newmez;
				newmez = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

   
		if (newmez)
			delete newmez;
		newmez = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		// if (pRefDate)
		// 	delete [] pRefDate;
		// pRefDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_CDO2   (double	Spread,
					 double		EffectiveDate,
					 double		EndDate,
					 double		ReferenceDate,
					 double		FirstCpnEffDate,
					 double		MezzAmount,
					 double		SubAmount,
					 int		FreqFeeLeg,
					 int		DayCount,
					 qPAYMENT_PREMIUM_LEG		AccruedOnDefault,
					 CCString	Currency,
					 long		PtfId, 	
					 int		CreditLag,
					 int		stub,
					 int		FreqDefLeg,
					 double		Binary,
					 CCString	PayCal,
					 double		RcvFee,
					 double		TradedNotional,
					 bool		CrossSubordination,
					 bool		IncludeMaturity,
					 ARM_result& result,
					 long objId)
{
	long cdo2Id;
	ICM_Cdo2* cdo2 = NULL;
	ICM_Cdo2* newcdo2 = NULL;
	// ARM_Currency* ccy = NULL;
	int stubrule = 3;// stub rule is shortend

	// char* payCal = NULL;

	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];
	int il=0;
	
	char* pPayCreditLag=new char[11];

	CCString msg ("");
	ICM_Portfolio* ptf = NULL;

	try
	{
		// Creation des ccy
		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		ARM_Date refDate; 
		if (ReferenceDate == -1.0)
		{	
			// strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		Local_XLDATE2ARMDATE(ReferenceDate,refDate);

		ptf = (ICM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(PtfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ptf, ICM_PF) == 0) 
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}

		ARM_Date FstCpnEffDate;
		if (FirstCpnEffDate != -1.0)
			Local_XLDATE2ARMDATE(FirstCpnEffDate,FstCpnEffDate);

		newcdo2 = new ICM_Cdo2((ARM_Date) pEffectiveDate,
							  (ARM_Date) pEndDate,
							  ReferenceDate==-1 ? 0 : &refDate,
							  FirstCpnEffDate==-1 ? 0 : &FstCpnEffDate,
							   (Spread/10000.),
							   	K_ADJUSTED,	// intRule
								K_ADJUSTED,	// adjStartDate
							   MezzAmount,
							   RcvFee*SubAmount,
							   //pRefDate,
							   FreqFeeLeg,
							   DayCount,
							   AccruedOnDefault,
							   CCSTringToSTLString(Currency), 
							   stubrule,
							   ptf,
							   CreditLag,
							   FreqDefLeg,
							   Binary,
							   CCSTringToSTLString(PayCal),
							   CrossSubordination,
							   IncludeMaturity);

		newcdo2->SetTradedCoef(TradedNotional/SubAmount);

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;

		// if (pRefDate)
		// 	delete [] pRefDate;
		// pRefDate = NULL;

		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		if (newcdo2 == NULL)
		{
			result.setMsg ("ARM_ERR: CDO2 is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdo2Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcdo2);

			if (cdo2Id == RET_KO)
			{
				if (newcdo2)
					delete newcdo2;
				newcdo2 = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdo2Id);

			return ARM_OK;
		}
		else
		{
			cdo2 = (ICM_Cdo2 *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cdo2, ICM_CDO2) == 1)
			{
				if (cdo2)
				{
					delete cdo2;
					cdo2 = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcdo2, objId);

				return ARM_OK;
			}
			else
			{
				if (newcdo2)
					delete newcdo2;
				newcdo2 = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();


		if (newcdo2)
			delete newcdo2;
		newcdo2 = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;


		if (pPayCreditLag)
			delete [] pPayCreditLag;
		pPayCreditLag = NULL;

		ARM_RESULT();
	}
}

long ICMLOCAL_CDSIndex(double	Spread,
					   int		IndexId,
					   double	EffectiveDate,
					   double	EndDate,
					   double	ReferenceDate,
					   int		Frequency,
					   int		DayCount,
					   double	FixedPayerAmount,
					   double	FloatingPayerAmount,
					   int		Stubrule,
					   CCString Currency,
					   int		intRule,
					   int		CreditLag,
					   bool		IncludeMaturity,
					   double	ProtectionStartDate,
					   double	ProtectionEndDate,
					   ARM_result& result,
					   long objId )
{							  
	long cdsId;

	ICM_Cds* cds = NULL;
	ICM_Cds* newcds = NULL;
	// ARM_Currency* ccy = NULL;
	
	int stubrule = 3;// stub rule is shortend

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}


	if (EffectiveDate>=EndDate)
	{
		result.setMsg ("Effective Date>=End Date");				
		return ARM_KO;
	}		

	ICM_Credit_Index* pIndex = NULL;

	char* pEffectiveDate=new char[11];
	char* pEndDate=new char[11];
	// char* pRefDate=new char[11];
	char* pProtectionStartDate=new char[11];
	char* pProtectionEndDate=new char[11];

	ARM_Date dProtectionStartDate;
	ARM_Date dProtectionEndDate;

	CCString msg ("");

	try
	{

		Local_XLDATE2ARMDATE(EffectiveDate,pEffectiveDate);
		Local_XLDATE2ARMDATE(EndDate,pEndDate);

		if (ProtectionStartDate<0) 
			dProtectionStartDate = (ARM_Date)pEffectiveDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionStartDate,pProtectionStartDate);
			dProtectionStartDate = (ARM_Date)pProtectionStartDate;
		}

		if (ProtectionEndDate<0) 
			dProtectionEndDate = (ARM_Date)pEndDate;
		else
		{
			Local_XLDATE2ARMDATE(ProtectionEndDate,pProtectionEndDate);
			dProtectionEndDate = (ARM_Date)pProtectionEndDate;
		}

		ARM_Date refDate ;
		if (ReferenceDate == -1.0)
		{	
			//strcpy(pRefDate,"NULL");
			stubrule = 1; // stub rule is shortstart
		}
		else
		{
			Local_XLDATE2ARMDATE(ReferenceDate,refDate);
			stubrule = Stubrule;
		}

		int FrequencyDefLeg = -1;

		pIndex = (ICM_Credit_Index *) LOCAL_PERSISTENT_OBJECTS->GetObject(IndexId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0) 
		{
			result.setMsg ("ARM_ERR: Index is not of a good type");
			return ARM_KO;
		}


		newcds = new ICM_Cds((ARM_Date) pEffectiveDate,
							   (ARM_Date) pEndDate,
							   ReferenceDate==-1 ? 0 : &refDate,
							   0,
							   dProtectionStartDate,
							   dProtectionEndDate,
							   (Spread/10000.),
								FixedPayerAmount,FloatingPayerAmount,0,0,
							   pIndex,
							   Frequency,
							   DayCount,
							   qCONTINUE_TO_MATURITY,
							   CCSTringToSTLString(Currency), 
							   stubrule,
							   CreditLag,
							   FrequencyDefLeg,
							   intRule,
							   IncludeMaturity,
								 K_ADJUSTED, // const int& adjStartDate /  ,
								 std::string(), // const std::string& payCalName   ,
								 qCDS_INDEX, // const qSecurity_TYPE& cdstype   ,
								 ISSUER_UNDEFINE, // const string& name   ,
								 CREDIT_DEFAULT_VALUE // const double& Binary );
							   );


		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;


		if (pProtectionStartDate)
			delete [] pProtectionStartDate;
		pProtectionStartDate = NULL;

		if (pProtectionEndDate)
			delete [] pProtectionEndDate;
		pProtectionEndDate = NULL;

		if (newcds == NULL)
		{
			result.setMsg ("ARM_ERR: CDS is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds);

			if (cdsId == RET_KO)
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdsId);

			return ARM_OK;
		}
		else
		{
			cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 1)
			{
				if (cds)
				{
					delete cds;
					cds = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcds, objId);

				return ARM_OK;
			}
			else
			{
				if (newcds)
					delete newcds;
				newcds = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();


		if (newcds)
			delete newcds;
		newcds = NULL;

		if (pEffectiveDate)
			delete [] pEffectiveDate;
		pEffectiveDate = NULL;

		if (pEndDate)
			delete [] pEndDate;
		pEndDate = NULL;


		ARM_RESULT();
	}
}




extern long ICMLOCAL_Option(string UnderlyingMaturity,
							ARM_Date UnderlyingMaturityDate,
							double Expiry,
							string ccy,
							qCDS_ADJ cds_adj,
							bool endAdj,
							double Strike,
							int OptionType,
							qDEF_MAT KO,
							double Notional,
							ARM_result& result,
							long objId)
{	
	ICM_Option* Option = NULL ;	
	char penddate[11];
	Local_XLDATE2ARMDATE(Expiry,penddate);
	ARM_Date expiry(penddate);

	Option = new ICM_Option(UnderlyingMaturity,UnderlyingMaturityDate,expiry,ccy,cds_adj,endAdj,
							Strike,OptionType,KO, Notional) ;
	long QId = LocalPersistent::get().adopt(Option,objId );
	return QId; 
}


extern long ICMLOCAL_SpreadOption(long idCDS,
							double Strike,
							bool isCall,
							int  koStyle,
							int  accelerationStyle,
							const std::vector<double>& exerciseDates,
							int exerciseStyle,
							int exerciseFrequency,
							int  matuStyle,
							ARM_result& result,
							long objId)
{
	long optId;
	


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	try
	{
		ICM_Cds* sec =dynamic_cast<ICM_Cds*>(LOCAL_PERSISTENT_OBJECTS->GetObject(idCDS)); 
		if (!sec) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find object id "<<idCDS); 

		ARM_Vector localExerciseDates(unconst(exerciseDates)); 
		for(unsigned int i=0;i<localExerciseDates.size();i++) 
			localExerciseDates[i]=XLDateToJulian(localExerciseDates[i]); 

		std::auto_ptr<ICM_SpreadOption> item(
			ICM_SpreadOption::build(Strike,isCall,(qKoStyle)koStyle,(qAccelerationStyle)accelerationStyle,localExerciseDates,exerciseStyle,exerciseFrequency,(qUnderlying_Maturity_Style)matuStyle,*sec) 
			); 
		
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();			
			optId=LOCAL_PERSISTENT_OBJECTS->SetPersistent(item.get()); 
			if ( optId== RET_KO ) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_ERR: Pb with inserting object"); 
			item.release(); 
			result.setLong(optId);
			return ARM_OK;
		}
		else
		{
			ICM_Security* tmp=dynamic_cast<ICM_Security*>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId)); 
			if (tmp) delete tmp; 
			LOCAL_PERSISTENT_OBJECTS->SetPersistent(item.release(), objId);
			result.setLong(objId); 
		}
	}

	catch(Exception& x)
    {
		CCString msg ("");
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_OK; 
}


// ---------------- CPPI 
extern long ICMLOCAL_CPPI  (double	StartDate,
							double MaturityDate,
							long idSecurity,
							CCString Currency,
							vector<double> Min,
							vector<double> Max,
							vector<double> ValueMin,
							vector<double> ValueMax,
							double Notional,
							double ProtectedAmount,
							double ManagementCost,
							double AdditionalLeverage,
							double DesactivateCushion,
							CCString CorrelName,
							ARM_result& result,
							long objId )
{
	long cppiId;
		
	int i = 0;
	
	ICM_Cppi* cppi = NULL;
	ICM_Cppi* newcppi = NULL;
	ARM_Currency* ccy = NULL;
	ICM_RangeFactor* riskyFactor = NULL;
	ICM_Security* sec = NULL;
	
	
	if (StartDate >= MaturityDate)
	{
		result.setMsg ("Start Date >= Maturity Date");				
		return ARM_KO;
	}		

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pStartDate = new char[11];
	char* pMaturityDate = new char[11];

	CCString msg ("");
	
	try
	{
		// Creation des ccy
		char* tmp = (char*)Currency;
		if (tmp)
		{
			ccy = new ARM_Currency(tmp);
		}
		tmp = NULL; 

		// Recuperation des dates
		Local_XLDATE2ARMDATE(StartDate,pStartDate);
		Local_XLDATE2ARMDATE(MaturityDate,pMaturityDate);

		// Creation du RangeFactor
		riskyFactor = new ICM_RangeFactor(Min,Max,ValueMin,ValueMax);

		// Recuperation de la security
		sec = (ICM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(idSecurity);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0) 
		{
			result.setMsg ("ARM_ERR: security is not of a good type");
			return ARM_KO;
		}

		newcppi = new ICM_Cppi(pStartDate,
							   pMaturityDate,
							   sec,
							   Notional,
							   ProtectedAmount,
							   AdditionalLeverage,
							   riskyFactor,
							   DesactivateCushion,
							   ManagementCost,
							   ccy,
							   CCSTringToSTLString(CorrelName)
							   );


		if (ccy)
			delete ccy;
		ccy = NULL;

		if (pStartDate)
			delete [] pStartDate;
		pStartDate = NULL;

		if (pMaturityDate)
			delete [] pMaturityDate;
		pMaturityDate = NULL;

		if (riskyFactor)
			delete riskyFactor;
		riskyFactor = NULL;

		if (newcppi == NULL)
		{
			result.setMsg ("ARM_ERR: CPPI is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cppiId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcppi);

			if (cppiId == RET_KO)
			{
				if (newcppi)
					delete newcppi;
				newcppi = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cppiId);

			return ARM_OK;
		}
		else
		{
			cppi = (ICM_Cppi *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cppi, ICM_CPPI) == 1)
			{
				if (cppi)
				{
					delete cppi;
					cppi = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcppi, objId);

				return ARM_OK;
			}
			else
			{
				if (newcppi)
					delete newcppi;
				newcppi = NULL;

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

		if (newcppi)
			delete newcppi;
		newcppi = NULL;

		if (pStartDate)
			delete [] pStartDate;
		pStartDate = NULL;

		if (pMaturityDate)
			delete [] pMaturityDate;
		pMaturityDate = NULL;

		if (riskyFactor)
			delete riskyFactor;
		riskyFactor = NULL;

		if (sec)
			delete sec;
		sec = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_NTDGEN(const int& CdsId,
					 const int& FirstNumDefault,
					 const int& LastNumDefault,
					 const int& CollateralId,
					 const double& Binary,
					 const double& RcvFee,
					 ARM_result& result,
					 long	objId)
{
	long NtdId;

	ICM_Cds* cds = NULL;
	ICM_Nthtd* ntd = NULL;
	ICM_Nthtd* prev_ntd = NULL;
	ICM_Collateral* collateral = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(CdsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 0) 
		{
			result.setMsg ("ARM_ERR: CDS is not of a good type");
			return ARM_KO;
		}

		collateral = (ICM_Collateral *) LOCAL_PERSISTENT_OBJECTS->GetObject(CollateralId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(collateral, ICM_COLLATERAL) == 0) 
		{
			result.setMsg ("ARM_ERR: Collateral is not of a good type");
			return ARM_KO;
		}

		ntd = new ICM_Nthtd(cds,FirstNumDefault,LastNumDefault,*collateral,Binary);
		ntd->SetPorS(RcvFee);

		if (ntd == NULL)
		{
			result.setMsg ("ARM_ERR: Ntd is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			NtdId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)ntd);

			if (NtdId == RET_KO)
			{
				if (ntd)
					delete ntd;
				ntd = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(NtdId);

			return ARM_OK;
		}
		else
		{
			prev_ntd = (ICM_Nthtd *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prev_ntd, ICM_NTD) == 1)
			{
				if (prev_ntd)
				{
					delete prev_ntd;
					prev_ntd = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)ntd, objId);

				return ARM_OK;
			}
			else
			{
				if (ntd)
					delete ntd;
				ntd = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (ntd)
			delete ntd;
		ntd = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_CDOGEN(const int& CdsId,
					 const double& SubAmountId,
					 const int strike_type,
					 const int& CollateralId,
					 const double& Binary,
					 const double& RcvFee,
					 ARM_result& result,
					 long	objId)
{
	long cdoId;

	ICM_Cds* cds = NULL;
	ICM_Mez* cdo = NULL;
	ICM_Mez* prev_cdo = NULL;
	ICM_Collateral* collateral = NULL;

	ARM_ReferenceValue *refStrike = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		
		if ( strike_type == 1 )
		{
			refStrike = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(DBL_INT(SubAmountId));
			refStrike = (ARM_ReferenceValue *) refStrike->Clone();

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(refStrike, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			refStrike = new ARM_ReferenceValue(SubAmountId);
		}

		cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(CdsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 0) 
		{
			result.setMsg ("ARM_ERR: CDS is not of a good type");
			return ARM_KO;
		}

		collateral = (ICM_Collateral *) LOCAL_PERSISTENT_OBJECTS->GetObject(CollateralId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(collateral, ICM_COLLATERAL) == 0) 
		{
			result.setMsg ("ARM_ERR: Collateral is not of a good type");
			return ARM_KO;
		}

		cdo = new ICM_Mez(cds,refStrike,*collateral,Binary);
		cdo->SetPorS(RcvFee);

		if (refStrike)
			delete refStrike;
		refStrike = NULL;

		if (cdo == NULL)
		{
			result.setMsg ("ARM_ERR: Cdo is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdoId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cdo);

			if (cdoId == RET_KO)
			{
				if (cdo)
					delete cdo;
				cdo = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdoId);

			return ARM_OK;
		}
		else
		{
			prev_cdo = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prev_cdo, ICM_MEZ) == 1)
			{
				if (prev_cdo)
				{
					delete prev_cdo;
					prev_cdo = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cdo, objId);

				return ARM_OK;
			}
			else
			{
				if (cdo)
					delete cdo;
				cdo = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (cdo)
			delete cdo;
		cdo = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_CDO_SQUARE_GEN(const int& CdsId,
							const double& SubAmount,
							const int& PortfolioId,
							const double& Binary,
							const double& RcvFee,
							ARM_result& result,
							long	objId)
{
	long cdoId;

	ICM_Cds* cds = NULL;
	ICM_Mez* cdo = NULL;
	ICM_Mez* prev_cdo = NULL;
	ICM_Portfolio* pf = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		cds = (ICM_Cds *) LOCAL_PERSISTENT_OBJECTS->GetObject(CdsId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cds, ICM_CDS) == 0) 
		{
			result.setMsg ("ARM_ERR: CDS is not of a good type");
			return ARM_KO;
		}

		pf = (ICM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(PortfolioId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ICM_PF) == 0) 
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}

		cdo = new ICM_Cdo2(cds,SubAmount,pf,Binary);
		cdo->SetPorS(RcvFee);

		if (cdo == NULL)
		{
			result.setMsg ("ARM_ERR: Cdo is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cdoId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cdo);

			if (cdoId == RET_KO)
			{
				if (cdo)
					delete cdo;
				cdo = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cdoId);

			return ARM_OK;
		}
		else
		{
			prev_cdo = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prev_cdo, ICM_CDO2) == 1)
			{
				if (prev_cdo)
				{
					delete prev_cdo;
					prev_cdo = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)cdo, objId);

				return ARM_OK;
			}
			else
			{
				if (cdo)
					delete cdo;
				cdo = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (cdo)
			delete cdo;
		cdo = NULL;

		ARM_RESULT();
	}
}

//	-------------------------------------------------------------------------------------
long 
ICMLOCAL_getLastFixingDate(long instId,double xlAsOfDate,ARM_result&result)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	ICM_Cds * item=dynamic_cast<ICM_Cds*>(LOCAL_PERSISTENT_OBJECTS->GetObject(instId)); 
	if (!item) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO; 
	}
	ARM_Date fixingDate; 
	char toto[11] ;
	Local_XLDATE2ARMDATE(xlAsOfDate,toto) ; 
	if ( ! item->GetFeeLeg()->GetCreditInfos()->getLastFixingDate(ARM_Date(toto),fixingDate) )
	{
		result.setMsg("NOFIXING"); 
		return ARM_OK; 
	}
	char dateStr[10]; 
	fixingDate.JulianToStrDate(dateStr); 
	result.setString(dateStr); 
	return ARM_OK; 
}

/* for corridor Leg
long 
ICMLOCAL_getFirstCreditFixingDate(long instId,double AsOfDate,ARM_result&result)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	ICM_CorridorLeg * item=dynamic_cast<ICM_CorridorLeg*>(LOCAL_PERSISTENT_OBJECTS->GetObject(instId)); 
	if (!item) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects security");
		return ARM_KO; 
	}
	ARM_Date fixingDate; 
	char toto[11] ;
	Local_XLDATE2ARMDATE(AsOfDate,toto) ; 
	if ( item->getFirstCreditFixingDate(ARM_Date(toto),fixingDate) )
	{
		result.setMsg("NOFIXING"); 
		return ARM_OK; 
	}
	char dateStr[10]; 
	fixingDate.JulianToStrDate(dateStr); 
	result.setString(dateStr); 
	return ARM_OK; 
}
*/
//	-------------------------------------------------------------------------------------
long 
ICMLOCAL_setPastFixing(long instId,double xlResetDate,double fixingValue,ARM_result&result)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	ICM_Cds * item=dynamic_cast<ICM_Cds*>(LOCAL_PERSISTENT_OBJECTS->GetObject(instId)); 
	if (!item) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO; 
	}
	char toto[11]; 
	Local_XLDATE2ARMDATE(xlResetDate,toto); 
	// item->GetFeeLeg()->GetCreditInfos()->setPastFixing(ARM_Date(toto),fixingValue); 
	item->GetFeeLeg()->setPastFixing(ARM_Date(toto),fixingValue); 
	result.setLong(instId); 
	return ARM_OK; 
}

long 
ICMLOCAL_SetPricerForRatesComputation(long legId,long pricerId,ARM_result& result)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ICM_Leg * leg=dynamic_cast<ICM_Leg*>(LOCAL_PERSISTENT_OBJECTS->GetObject(legId)); 
	if (!leg) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO; 
	}

	ICM_Pricer* pricer=dynamic_cast<ICM_Pricer*>(LOCAL_PERSISTENT_OBJECTS->GetObject(pricerId)); 
	if (!pricer) 
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO; 
	}

	leg->SetRatesPricer(pricer);
	
	result.setMsg("DONE"); 
	return ARM_OK; 
}

extern long ICMLOCAL_GetBounds(long CdoId,
							   double& low,
							   double& up,
							   ARM_result& result)
{
	ICM_Ftd* ftd=dynamic_cast<ICM_Ftd*>(LOCAL_PERSISTENT_OBJECTS->GetObject(CdoId)); 
	ftd->Bounds(low,up,ftd->GetFeeLeg()->GetStartDate());

	result.setMsg("DONE"); 
	return ARM_OK; 
}


extern long ICMLOCAL_Customized_CDO(
							const VECTOR<CCString>&	Labels,
							const VECTOR<double>&	notionals,
							const CCString			Currency,
							const long&				CP_DefaultLegId,
							const long&				CP_PremiumLegId,
							const long&				CP_ProductParametersId,
							ARM_result&			result,
							long				objId)
{
	double dResult=0.;
	CCString msg ("");

	long	Cutomized_CDO_Id;

	ICM_Customized_CDO*	Cutomized_CDO	=	NULL;
	ICM_Customized_CDO*	new_Customized_CDO	=	NULL;

	ICM_GenCF* TheCashFlows = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Default = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Premium = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_ProductParameters = NULL;

	ICM_Parameters* Default_Parameters = NULL;
	ICM_Parameters* Premium_Parameters = NULL;
	ICM_Parameters* Data_Parameters = NULL;

//	char* pEffectiveDate=new char[11];

	char** pLabels = NULL;
	double* pnotionals = NULL;
	
	try
	{
		// ---------------------------------------------------------------------------------
		// Effective Date
		// ---------------------------------------------------------------------------------
//		Local_XLDATE2ARMDATE(EffectiveDate, pEffectiveDate);
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Labels
		// ---------------------------------------------------------------------------------
		int	nbissuers = Labels.size();

		if (nbissuers)
		{
			pLabels = new char*[nbissuers];
			pnotionals = new double[nbissuers];

			for (int il=0; il<nbissuers; il++)
			{
				pLabels[il] = new char[60];
				strcpy(pLabels[il], (char*)Labels[il]);
				pnotionals[il] = notionals[il];
			}
		}

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
//		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_DefaultLegId); 

//		if (TheCashFlows)	Matrix_Default = TheCashFlows->GetMatrix();
		Default_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(CP_DefaultLegId)); 
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
//		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_PremiumLegId); 

//		if (TheCashFlows)	Matrix_Premium = TheCashFlows->GetMatrix();
		Premium_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(CP_PremiumLegId)); 
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Pricing Parameters, series of Cash Flows
		// ---------------------------------------------------------------------------------
//		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_ProductParametersId); 

//		if (TheCashFlows)	Matrix_ProductParameters = TheCashFlows->GetMatrix();
		Data_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(CP_ProductParametersId)); 
		// ---------------------------------------------------------------------------------

		new_Customized_CDO	= new ICM_Customized_CDO(
//										Matrix_Premium,
//										Matrix_Default,
//										Matrix_ProductParameters,
										Default_Parameters,
										Premium_Parameters,
										Data_Parameters,
										nbissuers,
										pLabels,
										pnotionals,
										CCSTringToSTLString(Currency)
											);


		if (new_Customized_CDO == NULL)
		{
			result.setMsg ("ARM_ERR: Cutomized CDO is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			Cutomized_CDO_Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Customized_CDO);

			if (Cutomized_CDO_Id == RET_KO)
			{
				if (new_Customized_CDO)
					delete new_Customized_CDO;
				new_Customized_CDO = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object Customized CDO");				
				return ARM_KO;
			}

			result.setLong(Cutomized_CDO_Id);

			return ARM_OK;
		}
		else
		{
			Cutomized_CDO = (ICM_Customized_CDO *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Cutomized_CDO, ICM_CUSTOMIZED_CDO) == 1)
			{
				if (Cutomized_CDO)
				{
					delete Cutomized_CDO;
					Cutomized_CDO = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Customized_CDO, objId);

				return ARM_OK;
			}
			else
			{
				if (new_Customized_CDO)
					delete new_Customized_CDO;
				new_Customized_CDO = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type - Customized CDO");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		ARM_RESULT();
	}
}

/* 
extern long ICMLOCAL_Option_Index_Gen(
							const long&		ScheduleParametersId,
							const long&		DataParametersId,
							ARM_result&		result,
							long			objId)
{
	double dResult=0.;
	CCString msg ("");

	long	Option_Index_Gen_Id;

	ICM_Option_Index_Gen*	Option_Index_Gen	=	NULL;
	ICM_Option_Index_Gen*	new_Option_Index_Gen	=	NULL;

	ICM_Parameters* Schedule_Parameters = NULL;
	ICM_Parameters* Data_Parameters = NULL;
	
	try
	{
		// ---------------------------------------------------------------------------------
		Schedule_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(ScheduleParametersId)); 
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		Data_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(DataParametersId)); 
		// ---------------------------------------------------------------------------------

		new_Option_Index_Gen	= new ICM_Option_Index_Gen(
										Schedule_Parameters,
										Data_Parameters
											);


		if (new_Option_Index_Gen == NULL)
		{
			result.setMsg ("ARM_ERR: Index Gen Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			Option_Index_Gen_Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Option_Index_Gen);

			if (Option_Index_Gen_Id == RET_KO)
			{
				if (new_Option_Index_Gen)
					delete new_Option_Index_Gen;
				new_Option_Index_Gen = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object Index Gen Option");				
				return ARM_KO;
			}

			result.setLong(Option_Index_Gen_Id);

			return ARM_OK;
		}
		else
		{
			Option_Index_Gen = (ICM_Option_Index_Gen *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Option_Index_Gen, ICM_OPTION_INDEX_GEN) == 1)
			{
				if (Option_Index_Gen)
				{
					delete Option_Index_Gen;
					Option_Index_Gen = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Option_Index_Gen, objId);

				return ARM_OK;
			}
			else
			{
				if (new_Option_Index_Gen)
					delete new_Option_Index_Gen;
				new_Option_Index_Gen = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type - Index Gen Option");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ICMLOCAL_Corridor_Index_Gen(
							const long&		ScheduleParametersId,
							const long&		DataParametersId,
							const long&		SubScheduleParametersId,
							ARM_result&		result,
							long			objId)
{
	double dResult=0.;
	CCString msg ("");

	long	Corridor_Index_Gen_Id;

	ICM_Corridor_Index_Gen*	Corridor_Index_Gen	=	NULL;
	ICM_Corridor_Index_Gen*	new_Corridor_Index_Gen	=	NULL;

	ICM_Parameters* Schedule_Parameters = NULL;
	ICM_Parameters* Data_Parameters = NULL;
	ICM_Parameters* Sub_Schedule_Parameters = NULL;
	
	try
	{
		// ---------------------------------------------------------------------------------
		Schedule_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(ScheduleParametersId)); 
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		Data_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(DataParametersId)); 
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		Sub_Schedule_Parameters	=	dynamic_cast<ICM_Parameters*> (LOCAL_PERSISTENT_OBJECTS->GetObject(SubScheduleParametersId)); 
		// ---------------------------------------------------------------------------------

		new_Corridor_Index_Gen	= new ICM_Corridor_Index_Gen(
										Schedule_Parameters,
										Data_Parameters,
										Sub_Schedule_Parameters
											);

		if (new_Corridor_Index_Gen == NULL)
		{
			result.setMsg ("ARM_ERR: Index Gen Corridor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			Corridor_Index_Gen_Id = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Corridor_Index_Gen);

			if (Corridor_Index_Gen_Id == RET_KO)
			{
				if (new_Corridor_Index_Gen)
					delete new_Corridor_Index_Gen;
				new_Corridor_Index_Gen = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object Index Gen Corridor");				
				return ARM_KO;
			}

			result.setLong(Corridor_Index_Gen_Id);

			return ARM_OK;
		}
		else
		{
			Corridor_Index_Gen = (ICM_Corridor_Index_Gen *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Corridor_Index_Gen, ICM_CORRIDOR_INDEX_GEN) == 1)
			{
				if (Corridor_Index_Gen)
				{
					delete Corridor_Index_Gen;
					Corridor_Index_Gen = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)new_Corridor_Index_Gen, objId);

				return ARM_OK;
			}
			else
			{
				if (new_Corridor_Index_Gen)
					delete new_Corridor_Index_Gen;
				new_Corridor_Index_Gen = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type - Index Gen Corridor");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		ARM_RESULT();
	}
}
*/ 
/*
extern long ICMLOCAL_Get_Instrument_Data(
							long				InstrumentId,
							VECTOR<double*>&	OutputMatrix,
							VECTOR<CCString>&	OutputLabels,
							int&				OutputNbRows,
							int&				OutputNbCols,
							ARM_result&			result)
{
	ARM_Security*	AnInstrument	=	NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	int	j;
	vector<string> AllLabels;

	try
	{
		// ---------------------------------------------------------------------------------
		AnInstrument = dynamic_cast<ARM_Security*> (  LOCAL_PERSISTENT_OBJECTS->GetObject(InstrumentId) );

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(AnInstrument, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Instrument is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		ARM_CLASS_NAME	Instrument_Name;
		Instrument_Name	=	AnInstrument->GetName();

		switch (Instrument_Name)
		{
		case ICM_CORRIDORLEG:

			((ICM_CorridorLeg*) AnInstrument)->GetSubSchedulesDataFromLabel(OutputMatrix, AllLabels, OutputNbRows, OutputNbCols);
			
			break;

			// ---------------------------------------------------------------------------------

		default:
			return	ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Labels
		OutputLabels.clear();
		for (j=0; j<OutputNbCols; j++)
			OutputLabels.push_back((AllLabels[j]).c_str());

		// ---------------------------------------------------------------------------------
		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;

}
*/


extern long ICMLOCAL_LssGapOption(long MezzId,
							long nblines_spread_triggers,
							long nbcols_spread_triggers,
							vector<double> spread_triggers,
							long nblines_default_triggers,
							long nbcols_default_triggers,
							vector<double> default_triggers,
							long nblines_mtm_triggers,
							long nbcols_mtm_triggers,
							vector<double> mtm_triggers,
							double	mtm_single_cond,
							ARM_result& result,
							long objId)
{
	long optId;
	
	ICM_Mez* mez = NULL;
	ICM_Lss_Gap_Option* OldOption = NULL ;
	ICM_Lss_Gap_Option* NewOption = NULL ;

	ICM_QMatrix<double> MAT_spread_triggers(nblines_spread_triggers,nbcols_spread_triggers);
	ICM_QMatrix<double> MAT_default_triggers(nblines_default_triggers,nbcols_default_triggers);
	ICM_QMatrix<double> MAT_mtm_triggers(nblines_mtm_triggers,nbcols_mtm_triggers);
	
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

	try
	{
		mez = (ICM_Mez*) LOCAL_PERSISTENT_OBJECTS->GetObject(MezzId);

		int i,j;

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mez, ICM_MEZ) == 0)
		{
			result.setMsg ("ARM_ERR: Mezz is not of a good type");
			return ARM_KO;
		}

		for ( i=0;i<nblines_spread_triggers;i++)
		for ( j=0;j<nbcols_spread_triggers;j++)
		{
			MAT_spread_triggers(i,j)=spread_triggers[i*nbcols_spread_triggers+j];
		}

		for ( i=0;i<nblines_default_triggers;i++)
		for ( j=0;j<nbcols_default_triggers;j++)
		{
			MAT_default_triggers(i,j)=default_triggers[i*nbcols_default_triggers+j];
		}

		for ( i=0;i<nblines_mtm_triggers;i++)
		for ( j=0;j<nbcols_mtm_triggers;j++)
		{
			MAT_mtm_triggers(i,j)=mtm_triggers[i*nbcols_mtm_triggers+j];
		}


		NewOption = new ICM_Lss_Gap_Option(mez,										 
										 &MAT_spread_triggers,
										 &MAT_default_triggers,
										 &MAT_mtm_triggers,
										 mtm_single_cond) ;
	
		if (NewOption == NULL)
		{
			result.setMsg ("ARM_ERR: Gap Option is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			optId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)NewOption);

			if (optId == RET_KO)
			{
				if (NewOption)
					delete NewOption;
				NewOption = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(optId);

			return ARM_OK;
		}
		else
		{
			OldOption = (ICM_Lss_Gap_Option *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(OldOption, ICM_LSS_GAP_OPTION) == 1)
			{
				if (OldOption)
				{
					delete OldOption;
					OldOption = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)NewOption, objId);

				return ARM_OK;
			}
			else
			{
				if (NewOption)
					delete NewOption;
				NewOption = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (NewOption)
			delete NewOption;
		NewOption = NULL;

		ARM_RESULT();
	}
}

long ICMLOCAL_CPDO (const long&			RiskyLegId,
					 const long&		RollLegId,
					 const long&		NoRiskyLegId,
					 const double&		CPDOMaturity,
					 const double&		InitialValo,
					 const double&		Target,
					 const CCString&	CouponType,
					 const double&		UFFees,
					 const double&		RunningFees,
					 const double&		VExpo,
					 const double&		V0Expo,
					 const double&		Alpha,
					 const double&		Beta,
					 const double&		Desactivation,
					 const double&		NbAssets,	
					 ARM_result& result,
					 long objId)
{
	long cpdoId;

	ICM_cpdo* cpdo = NULL;
	ICM_cpdo* newcpdo = NULL;
	ICM_Leg* RiskyLeg = NULL;
	ICM_Leg* RollLeg = NULL;
	ICM_Leg* NoRiskyLeg = NULL;
	//ARM_Currency* ccy = NULL;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	char* pMaturityDate=new char[11];

	// char** plabels = NULL;
	// std::vector<std::string> plabels; 
	// double* pnotionals = NULL;
	
	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE(CPDOMaturity, pMaturityDate);
		
		RiskyLeg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(RiskyLegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RiskyLeg, ICM_LEG) == 0) 
		{
			result.setMsg ("ARM_ERR: Leg is not of a good type");
			return ARM_KO;
		}

		RollLeg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(RollLegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(RollLeg, ICM_LEG) == 0) 
		{
			result.setMsg ("ARM_ERR: Leg is not of a good type");
			return ARM_KO;
		}

		NoRiskyLeg = (ICM_Leg *) LOCAL_PERSISTENT_OBJECTS->GetObject(NoRiskyLegId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(NoRiskyLeg, ICM_LEG) == 0) 
		{
			result.setMsg ("ARM_ERR: Leg is not of a good type");
			return ARM_KO;
		}

		newcpdo = new ICM_cpdo(RiskyLeg,
							RollLeg,
							NoRiskyLeg,
   							InitialValo,
							Target,
							(ARM_Date) pMaturityDate,  
							CCSTringToSTLString(CouponType), 
							UFFees,
							RunningFees,
							VExpo,
							V0Expo,
							Alpha,
							Beta,
							Desactivation,
							NbAssets);

		if (pMaturityDate)
			delete [] pMaturityDate;
		pMaturityDate = NULL;

		if (newcpdo == NULL)
		{
			result.setMsg ("ARM_ERR: CPDO is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			cpdoId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcpdo);

			if (cpdoId == RET_KO)
			{
				if (newcpdo)
					delete newcpdo;
				newcpdo = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(cpdoId);

			return ARM_OK;
		}
		else
		{
			cpdo = (ICM_cpdo *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(cpdo, ICM_CPDO) == 1)
			{
				if (cpdo)
				{
					delete cpdo;
					cpdo = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newcpdo, objId);

				return ARM_OK;
			}
			else
			{
				if (newcpdo)
					delete newcpdo;
				newcpdo = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
    {
		x.DebugPrint();

		if (newcpdo)
			delete newcpdo;
		newcpdo = NULL;

		if (pMaturityDate)
			delete [] pMaturityDate;
		pMaturityDate = NULL;


		ARM_RESULT();
	}
}
extern long ICMLOCAL_RESTRIKABLE_CDO(const ARM_Date* TriggerStartDate,
								const ARM_Date* Expiry,
								double Strike,
								double InitSpread,
								int OptionType,
								long idSecurity,
								double Rehauss,
								double Frequency,
								int DiffCDO,
								int IsCMSpread,
								double CMSpreadMatu,
								long objId)
{

	ICM_Mez*	Underlying;

	LocalPersistent::get().convert(idSecurity,Underlying);

	ICM_Option_Tranche* Restrikable_CDO ;

	int ExecType = K_BERMUDAN;

	Restrikable_CDO = new ICM_Option_Tranche(*TriggerStartDate,*Expiry,Strike,InitSpread,	OptionType,Underlying,Rehauss,Frequency,ExecType,DiffCDO,IsCMSpread,CMSpreadMatu);

	if( !Restrikable_CDO)

		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_RESTRIKABLE_CDO: Can't create Restrikable Tranche."); 

	return LocalPersistent::get().adopt(Restrikable_CDO,objId);
}


/*---- End Of File ----*/

// EOF %M%
