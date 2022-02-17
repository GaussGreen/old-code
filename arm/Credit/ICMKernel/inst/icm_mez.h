
#if !defined(_ICM_MEZ_H_)
#define _ICM_MEZ_H_

#include "ICMKernel\inst\icm_ftd.h"


/*********************************************************************************/
/*! \class  ICM_Mez icm_mez.h "icm_mez.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  <B> Mezzanine Default Swap </B> */
/***********************************************************************************/

class ICM_Credit_Index;

class ICM_Mez : public ICM_Ftd  
{

private: 
	ARM_ReferenceValue*		itsSubAmount;		 // Subordination Amount

public:

	ICM_Mez() {itsSubAmount=NULL; Init();  }	

	ICM_Mez(ICM_Cds* cds,
			double SubAmount,
			const ICM_Collateral& Collateral,
			const double& Binary = CREDIT_DEFAULT_VALUE) : ICM_Ftd(cds,Collateral,Binary)
	{ 
		itsSubAmount=NULL; 
		Init();
 		SetSubAmount(SubAmount);
	}	

	ICM_Mez(ICM_Cds* cds,
			ARM_ReferenceValue* SubAmount,
			const ICM_Collateral& Collateral,
			const double& Binary = CREDIT_DEFAULT_VALUE) : ICM_Ftd(cds,Collateral,Binary)
	{ 
		itsSubAmount=NULL; 
		Init();
 		SetSubAmount(SubAmount);
	}	



	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& FrequencyFeeLeg /* = K_QUARTERLY*/,
				const int& DayCountBasis /* = KACTUAL_360*/,
				const double& FixedPayerAmount /* = 100.*/, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/, 
				const double& FloatingPayerAmount /* = 0.0*/,
				const int& stubrule /* = K_SHORTSTART*/,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/,
				const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				const std::string&  payCalName /* = NULL*/,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);


// Meme constructeur avec un const vecteur<double>& à la place du double* passant les issuers notionals

	ICM_Mez(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const vector<string>& IssuersLabels,		 //char**	IssuersLabels,
				const vector<double>& IssuersNotionals,		 //double*	IssuersNotionals,
				const int& FrequencyFeeLeg /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				const double& FloatingPayerAmount /* = 0.0*/ ,
				const int& stubrule /* = K_SHORTSTART*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
				const std::string&  payCalName /* = NULL*/ ,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ); 

// JLA void Set(const ARM_Date& EffectiveDate,
// JLA 			const ARM_Date& ScheduleTerminateDate,
// JLA 			const ARM_Date* FirstPeriodReferenceDate,
// JLA 				const ARM_Date* FstCpnEffDate,
// JLA 				const double& FixedRate,
// JLA 				int intRule,
// JLA 				int adjStartDate,
// JLA 				const double& SubAmount,
// JLA 				const double& MezzAmount,
// JLA 				const vector<string>& IssuersLabels,// //char**	IssuersLabels,
// JLA 				const vector<double>& IssuersNotionals, // //double*	IssuersNotionals,
// JLA 				const int& FrequencyFeeLeg /* = K_QUARTERLY*/ ,
// JLA 				const int& DayCountBasis /* = KACTUAL_360*/ ,
// JLA 				const double& FixedPayerAmount /* = 100.*/ , 
// JLA 				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
// JLA 				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
// JLA 				const double& FloatingPayerAmount /* = 0.0*/ ,
// JLA 				const int& stubrule /* = K_SHORTSTART*/ ,
// JLA 				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
// JLA 				const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
// JLA 				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
// JLA 				const std::string& payCalName /* = NULL*/ ,
// JLA 				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ );

	ICM_Mez(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>&	IssuersNotionals,
				ICM_Credit_Index* Index,
				const double& ParticipationRate,
				const int& FrequencyFeeLeg /*= K_QUARTERLY*/ ,
				const int& DayCountBasis /* KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& ccy /* = ARM_DEFAULT_CURRENCY */ , 
				const double& FloatingPayerAmount /*= 0.0 */ ,
				const int& stubrule /* = K_SHORTSTART*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /*  CREDIT_DEFAULT_VALUE*/ ,
				const std::string& /* char* payCalName = NULL*/ ,
				const ARM_Date* FwdFixedDate /* = "NULL"*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ );

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date*FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>&	IssuersNotionals,
				// const int& NbIssuers,
				ICM_Credit_Index* Index,
				const double& ParticipationRate,
				const int& FrequencyFeeLeg /* = K_QUARTERLY */ ,
				const int& DayCountBasis /*= KACTUAL_360 */ ,
				const double& FixedPayerAmount /*= 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& ccy /*  ARM_DEFAULT_CURRENCY */ , 
				const double& FloatingPayerAmount /*= 0.0*/ ,
				const int& stubrule/* = K_SHORTSTART */ ,
				const double& CreditLag /*  DEFAULT_CREDIT_LAG*/ ,
				const int& FrequencyDefLeg /* DEFAULT_FRQ_DEFLEG */ ,
				const double& Binary /*  CREDIT_DEFAULT_VALUE */ ,
				const std::string& payCalName /*  NULL */ ,
				const ARM_Date* FwdFixedDate /* = "NULL" */ ,
				const bool& IncludeMaturity /*= INCLUDE_MATURITY*/ );

	virtual ~ICM_Mez() ; 

	void Init();

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	inline double GetSubAmount(const ARM_Date &Asof ) const
	{  
		return itsSubAmount->CptReferenceValue(unconst(Asof));
	}

	void SetSubAmount(const double& SubAmount)
	{  
		if (itsSubAmount)
			delete itsSubAmount;
		itsSubAmount = NULL;
		itsSubAmount = new ARM_ReferenceValue(SubAmount); 
	}

	void SetSubAmount(ARM_ReferenceValue* refval)
	{  
		if (itsSubAmount)
			delete itsSubAmount;
		itsSubAmount = NULL;
		itsSubAmount = (ARM_ReferenceValue*) refval->Clone(); 
	}

	double GetMezzAmount(const ARM_Date &Asof) const ;
 
	void SetMezzAmount(const double& MezzAmount);
	 

	void View(char* id, FILE* ficOut);

	virtual double GetPercentLow(const ARM_Date&  Asof) ;
	virtual double GetPercentHight(const ARM_Date& Asof) ;

	virtual void RebuildAfterDefault(ICM_ModelMultiCurves* mod);
	virtual void Bounds(double& low,double& up,const ARM_Date & AsOf) 
	{
		low=GetPercentLow(AsOf);
		up=GetPercentHight(AsOf);
	}

	inline ARM_ReferenceValue*	 GetSubAmount() { return itsSubAmount;}
	bool SearchBoundsForStepUp(const double& yf,const ARM_Date& Asof, vector<double>& dates);
private:
	ICM_Mez(const ICM_Mez&) ;//NA
	ICM_Mez& operator=(const ICM_Mez&); //NA 
};

#endif 
