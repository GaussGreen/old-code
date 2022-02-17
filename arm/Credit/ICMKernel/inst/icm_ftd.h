
#if !defined(_ICM_FTD_H_)
#define _ICM_FTD_H_

#include "ICMKernel\inst\icm_cds.h"

class ICM_Collateral ;
class ICM_ModelMultiCurves;

/*********************************************************************************/
/*! \class  ICM_ProportionsInfo icm_ftd.h "icm_ftd.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  <B> Proportions informations for summit </B> */
/***********************************************************************************/

/*********************************************************************************/
/*! \class  ICM_Ftd icm_ftd.h "icm_ftd.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_ftd.h
 *	\brief  <B> First to Default Definition </B> */
/***********************************************************************************/


class ICM_Ftd : public ICM_Cds  
{
private: 

	// ICM_ProportionsInfo* itsProportions;
	ICM_Collateral*		 itsCollateral;			// 0-1, aggreg.

public:

	ICM_Ftd() { Init();}	

	ICM_Ftd(ICM_Cds* cds,
			const ICM_Collateral& Collateral,
			const double& Binary = CREDIT_DEFAULT_VALUE) : ICM_Cds(*cds)
	{ 
		Init();
		Set(Binary,Collateral);
	}	

	void Set(const double& Binary,
			const ICM_Collateral& Collateral) 
	{ 
		SetBinary(Binary);
		SetCollateral(Collateral);
	}	

/* ***************************************************************************************************************** */

	ICM_Ftd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels,
				const int& FreqFeeLeg /* = K_QUARTERLY*/,
				const int& DayCountBasis /* = KACTUAL_360*/,
				const double& FixedPayerAmount /* = 100.*/, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
				const std::string& ccy, 
				const double& FloatingPayerAmount /* = 0.0*/,
				const int& stubrule /* = K_SHORTSTART*/,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				const std::string& payCalName /* = NULL*/,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/
				);

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>&IssuersLabels , 
				const int& FreqFeeLeg /* = K_QUARTERLY*/,
				const int& DayCountBasis /* = KACTUAL_360*/,
				const double& FixedPayerAmount /* = 100.*/, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
				const std::string& ccy, 
				const double& FloatingPayerAmount /* = 0.0*/,
				const int& stubrule /* = K_SHORTSTART*/,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				const std::string& payCalName /* = NULL*/,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);




/* ***************************************************************************************************************** */

	ICM_Ftd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& FreqFeeLeg /* = K_QUARTERLY */ ,
				const int& DayCountBasis /* = KACTUAL_360 */ ,
				const double& FixedPayerAmount /* = 100. */ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED */ ,
				const std::string& ccy /* =  */ , 
				const double& FloatingPayerAmount /* = 0.0*/ ,
				const int& stubrule /* = K_SHORTSTART */ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG */ ,
				const int& FreqDefLeg /*  DEFAULT_FRQ_DEFLEG */ ,
				const double& Binary /*  CREDIT_DEFAULT_VALUE*/ ,
				const std::string& /* char* payCalName = NULL */ ,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg */ ,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY */ ); 

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& FreqFeeLeg /* K_QUARTERLY */ ,
				const int& DayCountBasis /*  KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100. */ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED */ ,
				const std::string& /*  */ , 
				const double& FloatingPayerAmount /* = 0.0 */ ,
				const int& stubrule /* = K_SHORTSTART */ ,
				const double& CreditLag /*  DEFAULT_CREDIT_LAG */ ,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG */ ,
				const double& Binary /*  CREDIT_DEFAULT_VALUE */ ,
				const std::string& /* char* payCalName = NULL */ ,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg */ ,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/ ,
				const bool& IncludeMaturity /*  INCLUDE_MATURITY*/ ) ;


	
	virtual ~ICM_Ftd() ; 

	void Init();

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	virtual void ExcludeIssuer(const std::string&label) ; 

	void View(char* id, FILE* ficOut);


	virtual double GetPercentLow(const ARM_Date&AsOf) {return 0.;}
	virtual double GetPercentHight(const ARM_Date&AsOf) {return 0.;}


	ICM_Collateral* GetCollateral() { return itsCollateral;}
	const ICM_Collateral* GetCollateral() const { return itsCollateral;}

	void SetCollateral(const ICM_Collateral& value) ;

	void SetIssuersInfos(const std::vector<std::string>&IssuersLabels, const ARM_Vector&IssuersNotionals) ;
	

	virtual bool SearchBoundsForStepUp(const double& yf,const ARM_Date& Asof,vector<double>& Odates);
	virtual void RebuildAfterDefault(ICM_ModelMultiCurves* mod) {}
	virtual void Bounds(double& low,double& up, const ARM_Date &AsOf) 
	{
		low=0.;
		up=1.;
	}
private:
	ICM_Ftd(const ICM_Ftd&) ;//NA
	ICM_Ftd& operator=(const ICM_Ftd&); //NA 
};

#endif 
