#include "ARMKernel\glob\firsttoinc.h"
#include "ARMKernel\util\merge.h"
#include "ICMKernel\inst\icm_nthtd.h"
#include "ICMKernel\inst\icm_collateral.h"

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Nthtd::BitwiseCopy(const ARM_Object* src)
{
    ICM_Nthtd* ntd = (ICM_Nthtd *) src;

	itsFirstNumDefault = ntd->itsFirstNumDefault;
	itsLastNumDefault = ntd->itsLastNumDefault;

}

// -------------
//	Copy Method 
// -------------
void ICM_Nthtd::Copy(const ARM_Object* srccds)
{
     ICM_Ftd::Copy(srccds);
 
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Nthtd::Clone(void)
{
     ICM_Nthtd* theClone = new ICM_Nthtd();

     theClone->Copy(this);
 
     return(theClone);
}

ICM_Nthtd::ICM_Nthtd(ICM_Ftd* ftd)
{
	Init();

	ICM_Ftd::Copy(ftd);

	itsFirstNumDefault=0;
	itsLastNumDefault=1;

}

// ---------------------
//	Init of members data
// ---------------------

void ICM_Nthtd::Init()
{
	SetName(ICM_NTD);

	itsFirstNumDefault = 0;
	itsLastNumDefault = 0;
}

// ----------------------------------------------
//	Constructor of NTD (Reference Date is char*)
// ----------------------------------------------

ICM_Nthtd::ICM_Nthtd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date*FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				// double*	IssuersNotionals,
				const int& Frequency,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& Ccy, 
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string&  Paycal,
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
		FirstNumDefault,
		LastNumDefault,
		// NbIssuers,
		IssuersLabels,
		IssuersNotionals,
		Frequency,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		Ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FreqDefLeg,
		Binary,
		Paycal,
		IncludeMaturity); 

}


/** 
ICM_Nthtd::ICM_Nthtd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date*FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				// char** IssuersLabels_,
				// double*	IssuersNotionals,
				const std::vector<std::string>&IssuersLabels,
				const std::vector<double>&IssuersNotionals,
				const int& Frequency,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& Ccy, 
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string&  Paycal,
				const bool& IncludeMaturity) 
{

	Init();

	// std::vector<std::string> IssuersLabels(NbIssuers); 
	// int i;
	// for(i=0;i<NbIssuers;i++) 
	// 	IssuersLabels[i]=IssuersLabels_[i]; 

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		FixedRate,
		intRule,
		adjStartDate,
		FirstNumDefault,
		LastNumDefault,
		// NbIssuers,
		IssuersLabels,
		IssuersNotionals,
		Frequency,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		Ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FreqDefLeg,
		Binary,
		Paycal,
		IncludeMaturity); 

}

**/ 
// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

void ICM_Nthtd::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				// double*	IssuersNotionals,
				const int& Frequency,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& Ccy, 
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string&  Paycal,
				const bool& IncludeMaturity) 
{
	
	ICM_Ftd::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				FixedRate,
				intRule,
				adjStartDate,
				//NbIssuers,
				IssuersLabels,
				IssuersNotionals,
				Frequency,
				DayCountBasis,
				FixedPayerAmount, 
				AccruedOnDefault,
				Ccy, 
				FloatingPayerAmount,
				stubrule,
				CreditLag,
				FreqDefLeg,
				Binary,
				Paycal,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				IncludeMaturity); 

	itsFirstNumDefault = FirstNumDefault;
	itsLastNumDefault = LastNumDefault;
}


// ----------------------------------------------
//	Constructor of NTD (Reference Date is char*)
// ----------------------------------------------

ICM_Nthtd::ICM_Nthtd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				// char**	IssuersLabels,
				const std::vector<std::string>& IssuersLabels_,
				const int& Frequency ,
				const int& DayCountBasis ,
				const double& FixedPayerAmount , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault ,
				const std::string& Ccy , 
				const double& FloatingPayerAmount ,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string&  Paycal,
				const bool& IncludeMaturity)
{

	Init();

//	std::vector<std::string> IssuersLabels(NbIssuers); 
//	int i ;
//	for(i=0;i<NbIssuers;i++) IssuersLabels[i]=IssuersLabels_[i]; 

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		FixedRate,
		intRule,
		adjStartDate,
		FirstNumDefault,
		LastNumDefault,
		// NbIssuers,
		IssuersLabels_,
		Frequency,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		Ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FreqDefLeg,
		Binary,
		Paycal,
		IncludeMaturity); 
}

// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

void ICM_Nthtd::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date*FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				// char**	IssuersLabels_,
				const std::vector<std::string>& IssuersLabels_,
				const int& Frequency ,
				const int& DayCountBasis ,
				const double& FixedPayerAmount , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault ,
				const std::string& Ccy, 
				const double& FloatingPayerAmount ,
				const int& stubrule,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string&  Paycal,
				const bool& IncludeMaturity)
{

	// std::vector<std::string> IssuersLabels(NbIssuers); 
	// int i ;
	// for(i=0;i<NbIssuers;i++) IssuersLabels[i]=IssuersLabels_[i]; 

	ICM_Ftd::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				FixedRate,
				intRule,
				adjStartDate,
				//NbIssuers,
				IssuersLabels_,
				Frequency,
				DayCountBasis,
				FixedPayerAmount, 
				AccruedOnDefault,
				Ccy, 
				FloatingPayerAmount,
				stubrule,
				CreditLag,
				FreqDefLeg,
				Binary,
				Paycal,
				IncludeMaturity); 

	itsFirstNumDefault = FirstNumDefault;
	itsLastNumDefault = LastNumDefault;
}

// ----------------------------------------------------------------
// Rebuild NTD after default
// ----------------------------------------------------------------
void ICM_Nthtd::RebuildAfterDefault(ICM_ModelMultiCurves* mod)
{
	double nbdefault = 0;
	int i=0;

	vector<string> issuersindefault;
	GetCollateral()->GetIssuersInDefault(issuersindefault);
	nbdefault = issuersindefault.size();

	if (nbdefault==0) return; //no default

	itsFirstNumDefault -= nbdefault;
	itsFirstNumDefault = MAX(itsFirstNumDefault,0);
	itsLastNumDefault -= nbdefault;
	itsLastNumDefault = MAX(itsLastNumDefault,0);
	
	for (i=0;i<nbdefault;i++)
	{	GetCollateral()->TransferIssuerInDefaultCollateral((char*)(const char*)issuersindefault[i].c_str()); }

	if ((itsFirstNumDefault == itsLastNumDefault) &&
		(itsFirstNumDefault == 0))
	{
	GetFeeLeg()->SetCreditLegType(qNone_Leg);
	GetFeeLeg()->SetCreditLegStyle(qStyle_None_Leg);
	GetDefLeg()->SetCreditLegType(qNone_Leg);
	GetDefLeg()->SetCreditLegStyle(qStyle_None_Leg);
	}
}