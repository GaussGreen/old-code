#include "ARMKernel\glob\firsttoinc.h"
#include "ARMKernel\util\merge.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\inst\icm_cds.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ARMKernel\util\interpol.h"

ICM_Cds::ICM_Cds(const ICM_Cds& cds) : ICM_Security(cds)
{ 
	Init();
	BitwiseCopy(&cds);  
	
 	int nbflows = cds.GetFeeLeg().GetNumFlows(); 

}	
// virtual 
ICM_Cds::~ICM_Cds()
{
	if (itsFeeLeg)
		delete itsFeeLeg;
	itsFeeLeg = NULL;

	if (itsDefLeg)
		delete itsDefLeg;
	itsDefLeg = NULL;

	if (itsSchedule)
		delete itsSchedule;
	itsSchedule = NULL;

	if (itsRiskSchedule)
		delete itsRiskSchedule;
	itsRiskSchedule = NULL;

	//if (itsYearTermSearch)
	// 	delete itsYearTermSearch;
	// itsYearTermSearch = NULL;
}	
void ICM_Cds::ResetSchedules()
{
	if (itsSchedule) delete itsSchedule;
	itsSchedule=NULL;
	if (itsRiskSchedule) delete itsRiskSchedule;
	itsRiskSchedule=NULL;
	// if (itsYearTermSearch) delete itsYearTermSearch;
	// itsYearTermSearch=NULL;
}

qCredit_Leg_Type ICM_Cds::DeduceFeeLegType(const qSecurity_TYPE& sectype)
{
	switch (sectype)
	{
	case qRUNNING:
		return qRunning_Leg;
	case qFUNDED:
		return qFunded_Leg;
	case qINFINE:
		return qInfine_Leg;
	case qZEROCOUPON:
		return qZeroCoupon_Leg;
	case qCM_TRANCHE:
		return qConstant_Maturity_Leg;
	case qCDS_INDEX:
		return qConstant_Maturity_Leg;
	case qCM_CDS:
		return qCMCDS_Leg;
	case qBASKETasMEZ:
		return qRunning_Leg;
	default:
	throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Unknown LegType ");
	}

	return qNone_Leg;
}

qCredit_Leg_Type ICM_Cds::DeduceDefLegType(const qSecurity_TYPE& sectype)
{
	switch (sectype)
	{
	case qRUNNING:
	case qFUNDED:
	case qINFINE:
	case qZEROCOUPON:
	case qCM_TRANCHE:
	case qCDS_INDEX:
	case qCM_CDS:
	case qBASKETasMEZ:
		return qStandart_Recovery_Leg;
	default:
	throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Unknown LegType ");
	}

	return qNone_Leg;
}

void ICM_Cds::DisplayScheduleDates(int datesType,char* id, FILE* fOut)
{
	itsFeeLeg->DisplayScheduleDates(datesType,K_YES,id,fOut);
}

void ICM_Cds::DisplayScheduleValues(int valuesType, char* id, FILE* fOut)
{
	itsFeeLeg->DisplayScheduleValues(valuesType,id,fOut);
}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Cds::BitwiseCopy(const ARM_Object* srccds)
{
    ICM_Cds* cds = (ICM_Cds *) srccds;

	if (cds->itsFeeLeg)
	{
		if (itsFeeLeg)
			delete itsFeeLeg;
		itsFeeLeg = NULL;

		itsFeeLeg = dynamic_cast<ICM_Leg*>(cds->itsFeeLeg->Clone());
	}

	if (cds->itsDefLeg)
	{
		if (itsDefLeg)
			delete itsDefLeg;
		itsDefLeg = NULL;

		itsDefLeg = dynamic_cast<ICM_Leg*>(cds->itsDefLeg->Clone());
	}


	if (itsSchedule)
		delete itsSchedule;
	itsSchedule = NULL;

	// if (itsYearTermSearch)
	// 	delete itsYearTermSearch;
	// itsYearTermSearch = NULL;

	if (itsRiskSchedule)
		delete itsRiskSchedule;
	itsRiskSchedule = NULL;

	itsTradedCoef=cds->itsTradedCoef ; 
}

// -------------
//	Copy Method 
// -------------
void ICM_Cds::Copy(const ARM_Object* srccds)
{
	ICM_Security::Copy(srccds);
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Cds::Clone(void)
{
     ICM_Cds* theClone = new ICM_Cds();

     theClone->Copy(this);
 
     return(theClone);
}

ARM_Object* ICM_Cds::Clone() const 
{
     ICM_Cds* theClone = new ICM_Cds();
     theClone->Copy(this);
     return(theClone);
}

ICM_Cds::ICM_Cds(ICM_Leg* Feeleg,
				 ICM_Leg* Defleg,
				 const int& creditlag,
				 const string& name)
{
	Init();
	
	Set(Feeleg,Defleg,creditlag,name);
}

void ICM_Cds::Set(ICM_Leg* Feeleg,
			 ICM_Leg* Defleg,
			 const int& creditlag,
			 const string& name)
{
	Init();
	
	if (itsFeeLeg)
		delete itsFeeLeg;
	itsFeeLeg = dynamic_cast<ICM_Leg*>(Feeleg->Clone());
	
	if (itsDefLeg)
		delete itsDefLeg;
	itsDefLeg = dynamic_cast<ICM_Leg*>(Defleg->Clone());

	itsDefLeg->SetCreditLag(creditlag);
	SetSingleName(name);

	ICM_Security::Set(Feeleg->GetCreditInfosRef().GetStartDateNA(),
					  Feeleg->GetCreditInfosRef().GetEndDateNA(),
					  Feeleg->GetNumFlows(),
					  Feeleg->GetNotionalAmount(),
					  Feeleg->GetDayCount(),
					  Feeleg->GetCurrencyUnit()->GetCcyName());
					  
}

// ---------------------
//	Init of members data
// ---------------------

void ICM_Cds::Init()
{
	SetName(ICM_CDS);
	SetModel(NULL);

	itsFeeLeg = NULL;
	itsDefLeg = NULL;


	itsSchedule		=	NULL;
	itsRiskSchedule	=	NULL;
	// itsYearTermSearch	=	NULL;

	itsTradedCoef=1 ;
}

// ----------------------------------------------
//	Constructor of CDS (Reference Date is char*)
// ----------------------------------------------

ICM_Cds::ICM_Cds(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				const ARM_ReferenceValue& premiumNotionals,
				const ARM_ReferenceValue& defaultNotionals,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& Ccy, 
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const int& intRule,
				const bool& includematurity,
				const int& adjStartDate,
				const std::string& payCalName,
				const qCredit_Leg_Type& TypeFeeLeg,
				const qCredit_Leg_Type& TypeDefLeg,
				const string& name,
				const double& Binary,
				const long l_NotionalEch_Type ,
				const ARM_ReferenceValue* l_NotionalEchange )
{

	Init();

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		ProtectionStartDate,
		ProtectionEndDate,
		FixedRate,
		premiumNotionals,
		defaultNotionals,
		FrequencyFeeLeg,
		DayCountBasis,
		AccruedOnDefault,
		Ccy, 
		stubrule,
		CreditLag,
		FrequencyDefLeg,
		intRule,
		includematurity,
		adjStartDate,
		payCalName,
		TypeFeeLeg,
		TypeDefLeg,
		name,
		Binary,
		l_NotionalEch_Type,
		l_NotionalEchange
		);
}
// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------
void ICM_Cds::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				const ARM_ReferenceValue& premiumNotionals,
				const ARM_ReferenceValue& defaultNotionals,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy, 
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const int& intRule,
				const bool& includematurity,
				const int& adjStartDate,
				const std::string& payCalName,
				const qCredit_Leg_Type& TypeFeeLeg,
				const qCredit_Leg_Type& TypeDefLeg,
				const string& name,
				const double& Binary,
				const long l_NotionalEch_Type ,
				const ARM_ReferenceValue* l_NotionalEchange )
{

	int DefFreq = 0;

	if (FrequencyDefLeg == K_DEF_FREQ ) 
		DefFreq = FrequencyFeeLeg;
	else
		DefFreq = FrequencyDefLeg;

	// Premium Leg ---------------------------------------
	if (itsFeeLeg)
		delete itsFeeLeg;
	itsFeeLeg = new ICM_Leg(EffectiveDate,				// startDate, 
							  ScheduleTerminateDate,	// endDate, 
							  FirstPeriodReferenceDate, // For break period
							  FstCpnEffDate,
							  FixedRate,				// fixedRate,
							  AccruedOnDefault,			// AccruedOnDefault = qACCRUED_SETTLED,
							  KACTUAL_365,				// AccruedDayCount
							  0.,						// LastIndexFixing
							  K_RCV,					//rcvOrPay
							  FrequencyFeeLeg,			// freq = K_ANNUAL, 
							  DayCountBasis,			// dayCount= K30_360, 
							  K_COMP_PROP,				// decompFreq
							  K_ARREARS,				// payTiming
							  intRule,					// intRule
							  stubrule,					// stubRule
							  ccy,						// discountCcy = ARM_DEFAULT_CURRENCY,
						      payCalName,				// char* payCalName
							  K_NX_NONE,				// nxChange							  
							  includematurity,			//include maturity for last coupon computation 
							  adjStartDate,				//Adjusted startDate
							  TypeFeeLeg,				//Leg Type
							  Binary,					//Binary Recovery
							  name);					//Single Issuer Name

	itsFeeLeg->SetAmount(&unconst(premiumNotionals),100); 


	
	// Default Leg ---------------------------------------
	if (itsDefLeg)
		delete itsDefLeg;
	itsDefLeg = new ICM_Leg(ProtectionStartDate,				// startDate, 
							  ProtectionEndDate,	// endDate, 
							  FirstPeriodReferenceDate,
							  FstCpnEffDate, 
							  1.,						// fixedRate,
							  qACCRUED_SETTLED,			// AccruedOnDefault = qACCRUED_SETTLED,
							  KACTUAL_365,				// AccruedDayCount
							  0.,						// LastIndexFixing
							  K_PAY,					//rcvOrPay
							  DefFreq,					// freq = K_ANNUAL, 
							  DayCountBasis,			// dayCount= K30_360, 
							  K_COMP_PROP,				// decompFreq
							  K_ARREARS,				// payTiming
							  getDefLegInterestRule(intRule),	// intRule
							  stubrule,					// stubRule
							  ccy,						// discountCcy = ARM_DEFAULT_CURRENCY,
						      payCalName,				// char* payCalName
							  l_NotionalEch_Type,		//K_NX_NONE	// nxChange
							  false,					//include maturity for last coupon computation 
							  adjStartDate,				//Adjusted startDate
							  TypeDefLeg,				//Leg Type
							  Binary,					//Binary Recovery
							  name);					//Single Issuer Name

	itsDefLeg->SetAmount(&unconst(defaultNotionals),-100); 
	itsDefLeg->SetCreditLag(CreditLag);
	ARM_ReferenceValue* refValueCloned = NULL;
	if(l_NotionalEchange){
		refValueCloned = (ARM_ReferenceValue*)(((ARM_ReferenceValue*)l_NotionalEchange)->Clone());	
	}else {// creation d'une refValue
		/*ARM_Vector* pDates = new ARM_Vector(1);
		pDates[0] = ProtectionStartDate.GetJulian();
		ARM_Vector* pValues = new ARM_Vector(1);
		pValues[0] = unconst(defaultNotionals);*/
		refValueCloned = (ARM_ReferenceValue*)(((ARM_ReferenceValue)defaultNotionals).Clone());
		//refValueCloned = new ARM_ReferenceValue(pDates,pValues);
	}
	itsDefLeg->SetNotExchange(refValueCloned);
	//	JLA: à quoi ça sert ? 
	ICM_Security::Set(EffectiveDate,
					ScheduleTerminateDate,
					itsFeeLeg->GetNumFlows(),
					*itsFeeLeg->GetAmount(),// FixedPayerAmount,
					DayCountBasis,
					ccy); 
}


// -------------------------------------------------------------------------
// Constructor For CDS Index Product
// -------------------------------------------------------------------------

ICM_Cds::ICM_Cds(const ARM_Date& EffectiveDate,
				const ARM_Date& MaturityDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				double premiumNotional,
				double defaultNotional,
				const ARM_ReferenceValue* premiumNotionals,
				const ARM_ReferenceValue* defaultNotionals,
				ICM_Credit_Index* Index,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy, 
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const int& intRule,
				const bool& includematurity,
				const int& adjStartDate,
				const std::string& payCalName,
				const qSecurity_TYPE& cdstype,
				const string& name,
				const double& Binary)
{

	Init();

	Set(EffectiveDate,
		MaturityDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		ProtectionStartDate,
		ProtectionEndDate,
		FixedRate,
		premiumNotional,
		defaultNotional,
		premiumNotionals,
		defaultNotionals,
		Index,
		FrequencyFeeLeg,
		DayCountBasis,
		AccruedOnDefault,
		ccy, 
		stubrule,
		CreditLag,
		FrequencyDefLeg,
		intRule,
		includematurity,
		adjStartDate,
		payCalName,
		cdstype,
		name,
		Binary);

}


// ----------------------------
//	Set Method for CDS Index Product
// ----------------------------

void ICM_Cds::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& MaturityDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				double premiumNotional,
				double defaultNotional,
				const ARM_ReferenceValue* premiumNotionals,
				const ARM_ReferenceValue* defaultNotionals,
				ICM_Credit_Index* Index,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string & ccy, 
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const int& intRule,
				const bool& includematurity,
				const int& adjStartDate,
				const std::string& payCalName,
				const qSecurity_TYPE& cdstype,
				const string& name,
				const double& Binary)
{
	

	double spread = 0.0;
	double RedemptionValueFix = 100.0;
	double RedemptionValueFloat = -100.0;

	int DefFreq = 0;

	if (FrequencyDefLeg == -1) 
		DefFreq = FrequencyFeeLeg;
	else
		DefFreq = FrequencyDefLeg;

	// Premium Leg ---------------------------------------
	if (itsFeeLeg)
		delete itsFeeLeg;
	itsFeeLeg = new ICM_Leg(EffectiveDate,				// startDate, 
							  MaturityDate,				// endDate, 
							  FirstPeriodReferenceDate,
							  FstCpnEffDate,
							  FixedRate,				// fixedRate,
							  AccruedOnDefault,			// AccruedOnDefault = qACCRUED_SETTLED,
							  KACTUAL_365,				// AccruedDayCount
							  0.,						// LastIndexFixing
							  K_RCV,					//rcvOrPay
							  FrequencyFeeLeg,			// freq = K_ANNUAL, 
							  DayCountBasis,			// dayCount= K30_360, 
							  K_COMP_PROP,				// decompFreq
							  K_ARREARS,				// payTiming
							  intRule,					// intRule
							  stubrule,					// stubRule
							  ccy,						// discountCcy = ARM_DEFAULT_CURRENCY,
						      payCalName,				// char* payCalName
							  K_NX_NONE,				// nxChange
							  includematurity,			//include maturity for last coupon computation 
							  adjStartDate,				//Adjusted startDate
							  DeduceFeeLegType(cdstype),//Leg Type
							  Binary,					//Binary Recovery
							  name);					//Single Issuer Name


	if (premiumNotionals) 
		itsFeeLeg->SetAmount(&unconst(*premiumNotionals)); 
	else 
	{
		ARM_ReferenceValue refvalfix(premiumNotional, 1 /* price */, 0 /* K_CONSTANT */);
		itsFeeLeg->SetAmount(&refvalfix,RedemptionValueFix);
	}

	itsFeeLeg->SetIRIndex(Index);
	itsFeeLeg->CptCashFlowDates();
	itsFeeLeg->SetCreditIndex((ICM_Credit_Index*)Index->Clone());
	itsFeeLeg->SetIndexStyle(Index->GetIndexStyle());
	itsFeeLeg->SetRefSpread(FixedRate);

	// Default Leg ---------------------------------------
	if (itsDefLeg)
		delete itsDefLeg;
	itsDefLeg = new ICM_Leg(ProtectionStartDate,		// startDate, 
							  ProtectionEndDate,		// endDate, 
							  FirstPeriodReferenceDate,
							  FstCpnEffDate,
							  1.,						// fixedRate,
							  qACCRUED_SETTLED,			// AccruedOnDefault = qACCRUED_SETTLED,
							  KACTUAL_365,				// AccruedDayCount
							  0.,						// LastIndexFixing
							  K_PAY,					//rcvOrPay
							  DefFreq,					// freq = K_ANNUAL, 
							  DayCountBasis,			// dayCount= K30_360, 
							  K_COMP_PROP,				// decompFreq
							  K_ARREARS,				// payTiming
							  getDefLegInterestRule(intRule),					// intRule
							  stubrule,					// stubRule
							  ccy,						// discountCcy = ARM_DEFAULT_CURRENCY,
						      payCalName,				// char* payCalName
							  K_NX_NONE,				// nxChange
							  false,					//include maturity for last coupon computation 
							  adjStartDate,				//Adjusted startDate
							  DeduceDefLegType(cdstype),//Leg Type
							  Binary,					//Binary Recovery
							  name);					//Single Issuer Name

	if (defaultNotionals)
		itsDefLeg->SetAmount(&unconst(*defaultNotionals)); 
	else 
	{
		ARM_ReferenceValue refvalfloat(defaultNotional, 1 /* price */, 0 /* K_CONSTANT */);
		itsDefLeg->SetAmount(&refvalfloat,RedemptionValueFloat);
	}
	itsDefLeg->SetCreditLag(CreditLag);

 
	ICM_Security::Set(EffectiveDate,
					  MaturityDate,
					  itsFeeLeg->GetNumFlows(),
					  *itsFeeLeg->GetAmount(),
					  // 0,
					  DayCountBasis,
					  // DayCountBasis,
					  ccy); // Ccy->GetCcyName());
					  
	
	SetSecurityType(cdstype);

}



void ICM_Cds::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char  fOutName[200];

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
    else
    {
       fOut = ficOut;
    } 

	// ICM_Security::View(id, fOut);
	fprintf(fOut, "\t\t\t ----------------- Basis Informations ----------------- \n");
	fprintf(fOut,"Purchase Or Sale :%f\n",GetPorS());
	fprintf(fOut,"Traded Coeficient :%f\n",GetTradedCoef());
	fprintf(fOut,"Single Issuer Name :%s\n",GetSingleName().c_str());

    fprintf(fOut, "\t\t\t ----------------- Premium Leg ----------------- \n");
	GetFeeLeg()->View(id, fOut);
    fprintf(fOut, "\t\t\t ----------------- Default Leg ----------------- \n");
	GetDefLeg()->View(id, fOut);

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}

// -------------------------------------------------------------------------------
//Compute the schedule following start&end periods in year fractions starting at 0
// -------------------------------------------------------------------------------
void	ICM_Cds::GenerateSchedule(const ARM_Date& TheDate, const int& stepnbdays)
{
	if (itsSchedule != NULL) return ;

	ICM_Security* securityfee = GetFeeLeg()->GetCreditInfos();
	ICM_Security* securitydef = GetDefLeg()->GetCreditInfos();

	ARM_Vector* Schedule1 = NULL;
	ARM_Vector*	Schedule2 = NULL;
	ARM_Vector*	Schedule3 = NULL;
	ARM_Vector*	Schedule4 = NULL;

	ICM_Security* security = NULL;
		
	security = GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(TheDate);
	MergeDates(&Schedule1,security->GetYFAccStartDates(),security->GetYFAccEndDates());

	security = GetDefLeg()->GetCreditInfos();
	security->ComputeYF(TheDate);
	MergeDates(&Schedule2,security->GetYFAccStartDates(),security->GetYFAccEndDates());

	MergeDates(&Schedule3,Schedule1,Schedule2);	

	ARM_Vector* VTheDate = new ARM_Vector(1,0.);
	MergeDates(&Schedule4,VTheDate,Schedule3);

	for (int k=0; k<Schedule4->GetSize(); k++)
		if (Schedule4->Elt(k)>0.) break;
	
	itsSchedule = new ARM_Vector(Schedule4->GetSize() - k +1,0.);
	
	for (int l=0; l<itsSchedule->GetSize(); l++)
		itsSchedule->Elt(l)=Schedule4->Elt(l + k -1);

	if (Schedule1) delete Schedule1;
	if (Schedule2) delete Schedule2;
	if (Schedule3) delete Schedule3;
	if (Schedule4) delete Schedule4;
	if (VTheDate) delete VTheDate;

	//For Searches
	// if (itsYearTermSearch != NULL)
	// 	delete itsYearTermSearch;
	
	// ARM_Vector* vector = new ARM_Vector(itsSchedule->GetSize(),0.);

	// for (k=0;k<vector->GetSize();k++)
	// 	vector->Elt(k)=(double)k;

	// if (itsYearTermSearch)
	// 	delete itsYearTermSearch;
	// itsYearTermSearch = new ARM_ReferenceValue((ARM_Vector*)itsSchedule->Clone(),vector,K_STEPUP_LEFT,0);
	// itsYearTermSearch->SetCalcMethod(K_STEPUP_LEFT);

	if (stepnbdays>0)
	{
		//ICMMSG(WARN,"Using stepnbdays "<<stepnbdays); 

		int nbflowsfee = securityfee->GetAccStartDates().size();
		int nbflowsdef = securitydef->GetAccStartDates().size();
		ARM_Date StartFee = securityfee->GetAccStartDates().Elt(0);
		ARM_Date EndFee = securityfee->GetAccEndDates().Elt(nbflowsfee-1);
		ARM_Date StartDef = securitydef->GetAccStartDates().Elt(0);
		ARM_Date EndDef = securitydef->GetAccEndDates().Elt(nbflowsdef-1);

		double start = MAX(TheDate.GetJulian(),MIN(StartFee.GetJulian(),StartDef.GetJulian()));
		double end = MAX(EndFee.GetJulian(),EndDef.GetJulian());

		//conversion en julian days
		for (l=0; l<itsSchedule->GetSize(); l++)
		{itsSchedule->Elt(l)=TheDate.GetJulian()+(int)(itsSchedule->Elt(l)*365.+0.5);}

		ARM_Vector p(itsSchedule);
		if (itsSchedule) delete itsSchedule;
		itsSchedule = GenerateIntSchInludingOtherSchedule(start,end,&p,stepnbdays);

		//conversion en yearfractions
		for (l=0; l<itsSchedule->GetSize(); l++)
		{itsSchedule->Elt(l)=(itsSchedule->Elt(l)-TheDate.GetJulian())/365.;}
	}

}

// -------------------------------------------------------------------------------------
//Compute the schedule following Risky start&end periods in year fractions starting at 0
// -------------------------------------------------------------------------------------
void ICM_Cds::GenerateRiskSchedule(const ARM_Date& TheDate, const int& stepnbdays)
{
	if (itsRiskSchedule != NULL) return ;

	ICM_Security* securityfee = GetFeeLeg()->GetCreditInfos();
	ICM_Security* securitydef = GetDefLeg()->GetCreditInfos();
	
	int l=0;
	ARM_Vector* Schedule1 = NULL;
	ARM_Vector*	Schedule2 = NULL;
	ARM_Vector*	Schedule3 = NULL;
	ARM_Vector*	Schedule4 = NULL;

	securityfee->ComputeYF(TheDate);
	MergeDates(&Schedule1,securityfee->GetYFStartRiskDates(),securityfee->GetYFEndRiskDates());

	securitydef->ComputeYF(TheDate);
	MergeDates(&Schedule2,securitydef->GetYFStartRiskDates(),securitydef->GetYFEndRiskDates());
	MergeDates(&Schedule3,Schedule1,Schedule2);	

	ARM_Vector* VTheDate = new ARM_Vector(1,0.);
	MergeDates(&Schedule4,VTheDate,Schedule3);

	for (int k=0; k<Schedule4->GetSize(); k++) {if (Schedule4->Elt(k)>0.) break;}
	itsRiskSchedule = new ARM_Vector(Schedule4->GetSize() - k +1,0.);
	
	for (l=0; l<itsRiskSchedule->GetSize(); l++)
	{itsRiskSchedule->Elt(l)=Schedule4->Elt(l + k -1);}

	if (Schedule1) delete Schedule1;
	if (Schedule2) delete Schedule2;
	if (Schedule3) delete Schedule3;
	if (Schedule4) delete Schedule4;
	if (VTheDate) delete VTheDate;

	if (stepnbdays>0)
	{
		//ICMMSG(WARN,"Using stepnbdays "<<stepnbdays); 

		int nbflowsfee = securityfee->GetAccStartDates().size();
		int nbflowsdef = securitydef->GetAccStartDates().size();
		ARM_Date StartFee = securityfee->GetStartRiskDates().Elt(0);
		ARM_Date EndFee = securityfee->GetEndRiskDates().Elt(nbflowsfee-1);
		ARM_Date StartDef = securitydef->GetStartRiskDates().Elt(0);
		ARM_Date EndDef = securitydef->GetEndRiskDates().Elt(nbflowsdef-1);

		double start = MAX(TheDate.GetJulian(),MIN(StartFee.GetJulian(),StartDef.GetJulian()));
		double end = MAX(EndFee.GetJulian(),EndDef.GetJulian());

		//conversion en julian days
		for (l=0; l<itsRiskSchedule->GetSize(); l++)
		{itsRiskSchedule->Elt(l)=TheDate.GetJulian()+(int)(itsRiskSchedule->Elt(l)*365.+0.5);}

		ARM_Vector p(itsRiskSchedule);
		if (itsRiskSchedule) delete itsRiskSchedule;
		itsRiskSchedule = GenerateIntSchInludingOtherSchedule(start,end,&p,stepnbdays);

		//conversion en yearfractions
		for (l=0; l<itsRiskSchedule->GetSize(); l++)
		{itsRiskSchedule->Elt(l)=(itsRiskSchedule->Elt(l)-TheDate.GetJulian())/365.;}
	}
}



int	ICM_Cds::GetYearTermSearchInSchedule(const double& yearterm) const 
{
	// returns -1 on failure.
	return locateIndex(unconst(itsSchedule),yearterm); 
	// 
	// 	if (itsYearTermSearch == NULL)
	// 		return -1;
	// 	return itsYearTermSearch->CptReferenceValue(yearterm);
}
