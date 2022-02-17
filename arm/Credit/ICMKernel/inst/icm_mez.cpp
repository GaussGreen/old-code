
#include "ARMKernel/glob/firsttoinc.h"
#include "ARMKernel/util/merge.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\modelmulticurves.h"

// virtual
ICM_Mez::~ICM_Mez()
{
	if (itsSubAmount) delete itsSubAmount ; 
	itsSubAmount=NULL; 
}
// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Mez::BitwiseCopy(const ARM_Object* srccds)
{
	int i = 0;
    ICM_Mez* mezz = (ICM_Mez *) srccds;

	if (itsSubAmount) delete itsSubAmount;
	itsSubAmount =	(ARM_ReferenceValue*) mezz->itsSubAmount->Clone();

}

// -------------
//	Copy Method 
// -------------
void ICM_Mez::Copy(const ARM_Object* srccds)
{
     ICM_Ftd::Copy(srccds);
 
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Mez::Clone(void)
{
     ICM_Mez* theClone = new ICM_Mez();

     theClone->Copy(this);
 
     return(theClone);
}

// ---------------------
//	Init of members data
// ---------------------

void ICM_Mez::Init()
{
	SetName(ICM_MEZ);
	if (itsSubAmount) delete itsSubAmount; 
	itsSubAmount =NULL;
}

// Ce constructeur ne differe du precedent que par l'utilisation de vecteur de double et de string
ICM_Mez::ICM_Mez(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const vector<string>& IssuersLabels,	//char**	IssuersLabels,
				const vector<double>& IssuersNotionals,	//double*	IssuersNotionals,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy, 
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const double& Binary,
				const std::string& payCalName,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
				const bool& IncludeMaturity) 
{
	itsSubAmount=NULL; 
	Init();

		Set(EffectiveDate,
			ScheduleTerminateDate,
			FirstPeriodReferenceDate,
			FstCpnEffDate,
			FixedRate,
			intRule,
			adjStartDate,
			SubAmount,
			MezzAmount,
			IssuersLabels,
			IssuersNotionals,
			FrequencyFeeLeg,
			DayCountBasis,
			FixedPayerAmount, 
			AccruedOnDefault,
			ccy, 
			FloatingPayerAmount,
			stubrule,
			CreditLag,
			FrequencyDefLeg,
			Binary,
			payCalName,
			TypeFeeLeg,
			TypeDefLeg,
			IncludeMaturity);

}

// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

/** 

// Ce set ne differe du precedant que par l'utilisation de vector<double> et de vector<string>
void ICM_Mez::Set(const ARM_Date& EffectiveDate,
				  const ARM_Date& ScheduleTerminateDate,
				  const ARM_Date* pcFirstPeriodReferenceDate,
  				  const ARM_Date* FstCpnEffDate,
				  const double& FixedRate,
				  int intRule,
				  int adjStartDate,
				  const double& SubAmount,
				  const double& MezzAmount,
				  const vector<string>& IssuersLabels,			//char**	IssuersLabels,
				  const vector<double>& IssuersNotionals,		//double*	IssuersNotionals,
				  const int& FrequencyFeeLeg,
				  const int& DayCountBasis,
				  const double& FixedPayerAmount, 
				  const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				  const std::string& ccy, 
				  const double& FloatingPayerAmount,
				  const int& stubrule,
				  const double& CreditLag,
				  const int& FrequencyDefLeg,
				  const double& Binary,
				  const std::string& payCalName,
				  const bool& IncludeMaturity)
{
	
	ICM_Ftd::Set(EffectiveDate,
				ScheduleTerminateDate,
				pcFirstPeriodReferenceDate,
				FstCpnEffDate,
				FixedRate,
				intRule,
				adjStartDate,
				IssuersLabels,
				IssuersNotionals,
				FrequencyFeeLeg,
				DayCountBasis,
				FixedPayerAmount, 
				AccruedOnDefault,
				ccy, 
				FloatingPayerAmount,
				stubrule,
				CreditLag,
				FrequencyDefLeg,
				Binary,
				payCalName,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				IncludeMaturity);

	SetSubAmount(SubAmount);

}

  */ 
void ICM_Mez::Set(const ARM_Date& EffectiveDate,
				  const ARM_Date& ScheduleTerminateDate,
				  const ARM_Date* pcFirstPeriodReferenceDate,
  				  const ARM_Date* FstCpnEffDate,
				  const double& FixedRate,
				  int intRule,
				  int adjStartDate,
				  const double& SubAmount,
				  const double& MezzAmount,
				  const vector<string>& IssuersLabels,			//char**	IssuersLabels,
				  const std::vector<double>& IssuersNotionals,
				  const int& FrequencyFeeLeg,
				  const int& DayCountBasis,
				  const double& FixedPayerAmount, 
				  const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				  const std::string& ccy, 
				  const double& FloatingPayerAmount,
				  const int& stubrule,
				  const double& CreditLag,
				  const int& FrequencyDefLeg,
				  const double& Binary,
				  const std::string& payCalName,
				const qCredit_Leg_Type& TypeFeeLeg ,
				const qCredit_Leg_Type& TypeDefLeg ,
				  const bool& IncludeMaturity)
{
	
	ICM_Ftd::Set(EffectiveDate,
				ScheduleTerminateDate,
				pcFirstPeriodReferenceDate,
				FstCpnEffDate,
				FixedRate,
				intRule,
				adjStartDate,
				// NbIssuers,
				IssuersLabels,
				IssuersNotionals,
				FrequencyFeeLeg,
				DayCountBasis,
				FixedPayerAmount, 
				AccruedOnDefault,
				ccy, 
				FloatingPayerAmount,
				stubrule,
				CreditLag,
				FrequencyDefLeg,
				Binary,
				payCalName,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				IncludeMaturity);

	SetSubAmount(SubAmount);


}

// ----------------------------------------------
//	Constructor of Constant Maturity Tranche
// ----------------------------------------------

ICM_Mez::ICM_Mez(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				// const int& NbIssuers,
				ICM_Credit_Index* Index,
				const double& ParticipationRate,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy/** ARM_Currency *Ccy **/ ,
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const double& Binary,
				const std::string& payCalName,
				const ARM_Date* FwdFixedDate,
				const bool& IncludeMaturity)
{

	itsSubAmount=NULL; 
	Init();

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		intRule,
		adjStartDate,
		SubAmount,
		MezzAmount,
		IssuersLabels,
		IssuersNotionals,
		// NbIssuers,
		Index,
		ParticipationRate,
		FrequencyFeeLeg,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FrequencyDefLeg,
		Binary,
		payCalName,
		FwdFixedDate,
		IncludeMaturity);

}


// -----------------------------------------------
//	Set Method of Constant Maturity Tranche
// -----------------------------------------------

void ICM_Mez::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				ICM_Credit_Index* Index,
				const double& ParticipationRate,
				const int& FrequencyFeeLeg,
				const int& DayCountBasis,
				const double& FixedPayerAmount, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string& ccy,/* ARM_Currency *Ccy,*/  
				const double& FloatingPayerAmount,
				const int& stubrule,
				const double& CreditLag,
				const int& FrequencyDefLeg,
				const double& Binary,
				const std::string&payCalName /* char* payCalName*/ ,
				const ARM_Date* FwdFixedDate/* char* FwdFixedDate*/ ,
				const bool& IncludeMaturity)
{

	
	ICM_Ftd::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				0.00,
				intRule,
				adjStartDate,
				IssuersLabels,
				IssuersNotionals,
				FrequencyFeeLeg,
				DayCountBasis,
				FixedPayerAmount, 
				AccruedOnDefault,
				ccy, 
				FloatingPayerAmount,
				stubrule,
				CreditLag,
				FrequencyDefLeg,
				Binary,
				payCalName,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				IncludeMaturity);

	SetSubAmount(SubAmount);

	GetFeeLeg()->SetCreditLegType(DeduceFeeLegType(qCM_TRANCHE));
	GetDefLeg()->SetCreditLegType(DeduceDefLegType(qCM_TRANCHE));
	GetFeeLeg()->SetCreditIndex((ICM_Credit_Index*)Index->Clone());
	GetFeeLeg()->SetIndexStyle(Index->GetIndexStyle());
	//	The participation Rate is set as a ReferenceValue over the FeeLeg... 
	// GetFeeLeg()->SetRefSpread(ParticipationRate);
	GetFeeLeg()->SetRefPartRate(&ARM_ReferenceValue(ParticipationRate)); 
	// by default: the coupon will only pay FWD spread 
	GetFeeLeg()->SetFwdCalcTypes(&ARM_ReferenceValue(qCPN_FWD)); 
	GetFeeLeg()->CptCashFlowDates();

	GetFeeLeg()->SetFwdFixedDate(FwdFixedDate);
	SetSecurityType(qCM_TRANCHE);

}

void ICM_Mez::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Synthetic Tranche ----------------- \n");
    fprintf(fOut, "\tSubordination Amount :%f\n",GetSubAmount(GetFeeLeg()->GetStartDate()));
    fprintf(fOut, "\tPercent down :%f\n",GetPercentLow(GetFeeLeg()->GetStartDate()));
	fprintf(fOut, "\tPercent up :%f\n",GetPercentLow(GetFeeLeg()->GetStartDate())+GetPercentHight(GetFeeLeg()->GetStartDate()));

	ICM_Ftd::View(id, fOut);

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


void ICM_Mez::RebuildAfterDefault(ICM_ModelMultiCurves* mod)
{
	double Sens = 1.;
	double subamount2 = GetSubAmount(GetFeeLeg()->GetStartDate());
	double InitNotional = GetFeeLeg()->GetCreditInfos()->GetNotionals().Elt(0); // GetInitialNotional();
	double notional2 = InitNotional;
	double nbdefault = 0;
	double SumRecXIssuerNot=0.;
	int i=0;

	//recuperation du sens de l'operation
	Sens=(notional2>0.) ? 1.:-1.;

	vector<string> issuersindefault;
	//recuperation des issuers en defaut
	GetCollateral()->GetIssuersInDefault(issuersindefault);
	nbdefault = issuersindefault.size();

	if (nbdefault==0) return; //no default

	double SizeAttackNot = 0.;
	
	//on retranche du montant de subordination (1-R)*NotIssuer
	for (i=0;i<nbdefault;i++)
	{ 
		int issuerRank = GetCollateral()->getIssuerPosition(issuersindefault[i]); 
		double IssuerNot = GetCollateral()->GetIssuersNotional(issuerRank);

		//cas binaire
		if (GetBinaryFlg()) {
			subamount2 -= IssuerNot*(1.-GetBinary()); 
			SumRecXIssuerNot += IssuerNot*GetBinary();
		}//cas leverage recovery
		else if (GetCollateral()->GetRecovCoefFlg()){
			subamount2 -= IssuerNot*(1.-MIN(mod->GetRecoveryRate(issuersindefault[i])*GetCollateral()->GetRecovCoef(),1));
			SumRecXIssuerNot += IssuerNot*MIN(mod->GetRecoveryRate(issuersindefault[i])*GetCollateral()->GetRecovCoef(),1);
		}//cas général
		else {
			subamount2 -= IssuerNot*(1.-mod->GetRecoveryRate(issuersindefault[i]));
			SumRecXIssuerNot += IssuerNot*mod->GetRecoveryRate(issuersindefault[i]);
		}
	}

	//on trancfert les issuers en defaut dans le collateral de defauts
	for (i=0;i<nbdefault;i++)
	{ GetCollateral()->TransferIssuerInDefaultCollateral(issuersindefault[i]);}

	if (subamount2<0.) SizeAttackNot = -subamount2;
	subamount2 = MAX(0.,subamount2);
	notional2 -= Sens*SizeAttackNot;//si la tranche est attaquée, on impacte le notional de la tranche

	//dans le cas d'une super senior avec un point d'attachement à 100%, on impacte le notionel de la tranche 
	if (fabs(subamount2+notional2)/GetCollateral()->SumNotionals(GetStartDateNA())>1.) {notional2 -= Sens*SumRecXIssuerNot;}

	if (Sens>0.) notional2 = MAX(0.,notional2);
	else notional2 = MIN(0.,notional2);

	
	//case everything is dead
	if ((CHECK_NULL(subamount2))&& (CHECK_NULL(notional2)))
	{
	GetFeeLeg()->SetCreditLegType(qNone_Leg);
	GetFeeLeg()->SetCreditLegStyle(qStyle_None_Leg);
	GetDefLeg()->SetCreditLegType(qNone_Leg);
	GetDefLeg()->SetCreditLegStyle(qStyle_None_Leg);
	}

	SetSubAmount(subamount2);
	GetFeeLeg()->GetCreditInfos()->SetFixedNotional(notional2);
	// 14514 GetFeeLeg()->GetCreditInfos()->SetInitialNotional(notional2);
	GetDefLeg()->GetCreditInfos()->SetFixedNotional(notional2);
	// 14514 GetDefLeg()->GetCreditInfos()->SetInitialNotional(notional2);
	SetFixedNotional(notional2);
	// 14514 SetInitialNotional(notional2);
 
}

bool ICM_Mez::SearchBoundsForStepUp(const double& yf,const ARM_Date& Asof,vector<double>& Odates)
{
	Odates.clear();
	double valueMAX = 5.e6; 
	ARM_Date date = yf*365. + Asof.GetJulian();
	vector<double> dates;
	vector<double> dates1;
	vector<double> dates2;
	ARM_ReferenceValue* ref = itsSubAmount;
	bool output = true;
	ICM_Leg* leg = GetDefLeg();
	ARM_ReferenceValue* amount = leg->GetAmount();
	int i=0;

	//cas ou on dispose d'une seule date
	if ((ref->GetDiscreteDates()==NULL) && (amount->GetDiscreteDates()==NULL))
	{ Odates.push_back(date.GetJulian()); 
		return (false);
	}
	else if (ref->GetDiscreteDates()==NULL) {}
	else if (amount->GetDiscreteDates()==NULL) {}

	if ((ref->GetDiscreteDates()->GetSize()== amount->GetDiscreteDates()->GetSize()) && 
	   (ref->GetDiscreteDates()->GetSize()== 1))
	{	if CHECK_EQUAL(ref->GetDiscreteDates()->Elt(0),amount->GetDiscreteDates()->Elt(0))
		{Odates.push_back(date.GetJulian()); 
		return (false);}
	}

	dates1.push_back(0.);
	for (i=0; i<ref->GetDiscreteDates()->GetSize();i++)
	{	dates1.push_back(ref->GetDiscreteDates()->Elt(i)); }
	dates1.push_back(valueMAX);

	dates2.push_back(0.);
	for (i=0; i<amount->GetDiscreteDates()->GetSize();i++)
	{	dates2.push_back(amount->GetDiscreteDates()->Elt(i)); }
	dates2.push_back(valueMAX);
	
	MergeVector(dates1,dates2,dates);

	double inf_= FlatVectorInterpol(dates,dates,date.GetJulian(),true);
	double sup_= FlatVectorInterpol(dates,dates,date.GetJulian(),false);
/*
	if CHECK_EQUAL(inf_,sup_)
	{	inf_= FlatVectorInterpol(dates,dates,date.GetJulian()-1.,true);
		sup_= FlatVectorInterpol(dates,dates,date.GetJulian()-1.,false);}
*/

	//cas ou la date inf correspond à la 1ere date de la référence value
	if CHECK_EQUAL(inf_,MIN(ref->GetDiscreteDates()->Elt(0),amount->GetDiscreteDates()->Elt(0)))
	{	Odates.push_back(date.GetJulian()); 
		return (false);	}

	if CHECK_EQUAL(inf_,0.)
	{	Odates.push_back(date.GetJulian()); 
		return (false);
	}
	else //if CHECK_EQUAL(sup_,valueMAX)
	{	Odates= dates;
		Odates.resize(Odates.size()-1);
		for (i=0; i<Odates.size()-1;i++)
		{Odates[i]=Odates[i+1];};
		Odates.resize(Odates.size()-1);
		Odates.push_back(date.GetJulian()); 
	}
/*	else
	{
		for (i=0; i<dates.size();i++)
		{	if (dates[i]<=inf_) 
			{Odates.push_back(dates[i]);}}
	}
*/	

	return (output);

}

// virtual 
double ICM_Mez::GetPercentLow(const ARM_Date  &date) 
{	double totnot = GetCollateral()->SumNotionals(date);
	return (GetSubAmount(date)/totnot);	
}


// virtual 
double ICM_Mez::GetPercentHight(const ARM_Date &date) 
{	double totnot = GetCollateral()->SumNotionals(date);
	double sizetrch = fabs(GetMezzAmount(date));
	return (sizetrch/totnot);
}

void ICM_Mez::SetMezzAmount(const double& MezzAmount)
{  

	double RedemptionValueFix = 100.0;
	double RedemptionValueFloat = -100.0;

	SetFixedNotional(MezzAmount);
	// JLA14514 SetInitialNotional(MezzAmount) ;

	ARM_ReferenceValue* refvalfix = new ARM_ReferenceValue(MezzAmount , 1 /* price */, 0 /* K_CONSTANT */);

 	GetFeeLeg()->SetAmount(refvalfix,RedemptionValueFix);
	GetFeeLeg()->ResetCreditInfos(); 

	GetDefLeg()->SetAmount(refvalfix,RedemptionValueFloat);
	GetDefLeg()->ResetCreditInfos(); 

	if (refvalfix)
		delete refvalfix;
	refvalfix = NULL;
}
double ICM_Mez::GetMezzAmount(const ARM_Date &Asof) const
{  
	double MezzAmount = 0.;

	//if (Asof){
	int i = GetFeeLeg().GetCreditInfos()->PeriodIndex(Asof) ;	// index of the period including AsOf (-1 if no period) 
	MezzAmount = ABS(GetFeeLeg().GetCreditInfos()->GetNotionals().Elt(i));
	return (MezzAmount); 
}