#include "ARMKernel\glob\firsttoinc.h"
#include "ARMKernel\util\merge.h"
#include "ICMKernel\inst\icm_ftd.h"
#include "ICMKernel\inst\icm_collateral.h"

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Ftd::BitwiseCopy(const ARM_Object* srccds)
{
    ICM_Ftd* cds = (ICM_Ftd *) srccds;

	int i = 0;
	if (itsCollateral)
		delete itsCollateral;
	if (cds->itsCollateral) itsCollateral = (ICM_Collateral*) cds->itsCollateral->Clone();

}
//	-------------------------------------------------------------------
ICM_Ftd::~ICM_Ftd() 
{
if (itsCollateral) delete itsCollateral;itsCollateral=0; 
}
// -------------
//	Copy Method 
// -------------
void ICM_Ftd::Copy(const ARM_Object* srccds)
{
     ICM_Cds::Copy(srccds);
 
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Ftd::Clone(void)
{
     ICM_Ftd* theClone = new ICM_Ftd();

     theClone->Copy(this);
 
     return(theClone);
}

// ---------------------
//	Init of members data
// ---------------------

void ICM_Ftd::Init()
{
	SetName(ICM_FTD);

	// itsProportions = NULL;
	itsCollateral = NULL; 

}

// ----------------------------------------------
//	Constructor of FTD (Reference Date is char*)
// ----------------------------------------------

ICM_Ftd::ICM_Ftd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& FreqFeeLeg,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy/* ARM_Currency *Ccy */ , 
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string& payCalName/* char* payCalName */ ,
				const qCredit_Leg_Type& TypeFeeLeg,
				const qCredit_Leg_Type& TypeDefLeg,
				const bool& IncludeMaturity) 
{

	Init();

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		FixedRate,
		intRule,
		adjStartDate,
		IssuersLabels,
		IssuersNotionals,
		FreqFeeLeg,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FreqDefLeg,
		Binary,
		payCalName,
		TypeFeeLeg,
		TypeDefLeg,
		IncludeMaturity) ;
}


// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

void ICM_Ftd::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& FreqFeeLeg,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy,
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string& payCalName, // char* payCalName,
				const qCredit_Leg_Type& TypeFeeLeg,
				const qCredit_Leg_Type& TypeDefLeg,
				const bool& IncludeMaturity)
{
	
	ICM_Cds::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				EffectiveDate,
				ScheduleTerminateDate,
				FixedRate,
				ARM_ReferenceValue(FixedPayerAmount), 
				ARM_ReferenceValue(FloatingPayerAmount), // 0,0,	// notionals 
				FreqFeeLeg,
				DayCountBasis,
				AccruedOnDefault,
				ccy, 
				stubrule,
				CreditLag,
				FreqDefLeg,
				intRule, // K_ADJUSTED,/*intRule*/
				IncludeMaturity,
				adjStartDate, // 1, /*adjStartDate*/
				payCalName,
				TypeFeeLeg,
				TypeDefLeg,
				ISSUER_UNDEFINE,
				Binary);

	if (itsCollateral)
		delete itsCollateral;
	itsCollateral = new ICM_Collateral(IssuersLabels,ARM_Vector(IssuersNotionals));
}


// ----------------------------------------------
//	Constructor of FTD (Reference Date is char*)
// ----------------------------------------------

ICM_Ftd::ICM_Ftd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels ,
				const int& FreqFeeLeg ,
				const int& DayCountBasis ,
				const double& FixedPayerAmount , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault ,
				const std::string& ccy , 
				const double& FloatingPayerAmount ,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string& PayCalName,
				const bool& IncludeMaturity)
{

	Init();

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		FixedRate,
		intRule,
		adjStartDate,
		// NbIssuers,
		IssuersLabels,
		// FirstPeriodReferenceDate,
		FreqFeeLeg,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FreqDefLeg,
		Binary,
		PayCalName,
		IncludeMaturity);
}


// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

void ICM_Ftd::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>&IssuersLabels, 
				const int& FreqFeeLeg ,
				const int& DayCountBasis ,
				const double& FixedPayerAmount , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault ,
				const std::string& ccy, 
				const double& FloatingPayerAmount ,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string& payCalName,
				const bool& IncludeMaturity)
{
	
	ICM_Cds::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				EffectiveDate,
				ScheduleTerminateDate,
				FixedRate,
				ARM_ReferenceValue(FixedPayerAmount), 
				ARM_ReferenceValue(FloatingPayerAmount),// 0,0,	// notionals 
				FreqFeeLeg,
				DayCountBasis,
				AccruedOnDefault,
				ccy, 
				stubrule,
				CreditLag,
				FreqDefLeg,
				intRule,// K_ADJUSTED,/*intRule*/
				IncludeMaturity, 
				adjStartDate, // 1, /*adjStartDate*/
				payCalName,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				ISSUER_UNDEFINE,
				Binary);


	if (itsCollateral)delete itsCollateral;
	
	// std::vector<std::string> issuers(NbIssuers); 
	ARM_Vector notios(IssuersLabels.size()); 
	for(int i=0;i<IssuersLabels.size();i++) 
	{
		notios[i]=fabs(FixedPayerAmount);
	}
	itsCollateral = new ICM_Collateral(IssuersLabels,notios);

	SetBinary(Binary);
}

 bool ICM_Ftd::SearchBoundsForStepUp(const double& yf,const ARM_Date& Asof,vector<double>& Odates)
{
	bool output = false;

	return (output);

}
void ICM_Ftd::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Basket Viewer ----------------- \n");
    fprintf(fOut, "\tBinary Value :%f\n",GetBinary());

	ICM_Cds::View(id, fOut);

	// if (itsProportions) itsProportions->View(id, fOut);
	if (itsCollateral) itsCollateral->View(id, fOut);

	fprintf(fOut, "\n");

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}



// virtual 
void ICM_Ftd::ExcludeIssuer(const std::string&label) 
{ 
	itsCollateral->ExcludeIssuer(label);
}


void ICM_Ftd::SetCollateral(const ICM_Collateral& value) 
{
	// if (value== NULL) return;
	if (itsCollateral) delete itsCollateral;itsCollateral=0; 
	itsCollateral = dynamic_cast<ICM_Collateral*>(unconst(value).Clone());
}
// 
void ICM_Ftd::SetIssuersInfos(const std::vector<std::string>&IssuersLabels, const ARM_Vector&IssuersNotionals)
{
	itsCollateral->SetIssuersInfos(IssuersLabels,IssuersNotionals);
}	
