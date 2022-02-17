
#if !defined(_ICM_NTHTD_H_)
#define _ICM_NTHTD_H_

#include "ICMKernel\inst\icm_ftd.h"


/*********************************************************************************/
/*! \class  ICM_Nthtd icm_nthtd.h "icm_nthtd.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  <B> Nth to Default Definition </B> */
/***********************************************************************************/

class ICM_Nthtd : public ICM_Ftd  
{
private: 

		int itsFirstNumDefault;
		int itsLastNumDefault;

public:

	ICM_Nthtd()	{ Init(); }	

	ICM_Nthtd(ICM_Cds* cds,
			const int&  FirstNumDefault,
			const int&  LastNumDefault,
			const ICM_Collateral& Collateral,
			const double& Binary = CREDIT_DEFAULT_VALUE) : ICM_Ftd(cds,Collateral,Binary)
	{ 
		Init();
		SetFirstNumDefault(FirstNumDefault);
		SetLastNumDefault(LastNumDefault);
	}

/* ***************************************************************************************************************** */
/* ***************************************************************************************************************** */

	ICM_Nthtd(const ARM_Date& EffectiveDate,
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
				const std::vector<std::string>& IssuersLabels,
				const int& Frequency /* = K_QUARTERLY */,
				const int& DayCountBasis /* = KACTUAL_360*/,
				const double& FixedPayerAmount /* = 100.*/, 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/, 
				const double& FloatingPayerAmount /* = 0.0*/,
				const int& stubrule /* = K_SHORTSTART*/,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				const std::string&  Paycal /* = NULL*/,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/) ;

	void Set(const ARM_Date& EffectiveDate,
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
				const std::vector<std::string>& IssuersLabels,
				const int& Frequency /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& Ccy/* = ARM_DEFAULT_CURRENCY*/ , 
				const double& FloatingPayerAmount /* = 0.0*/ ,
				const int& stubrule /* = K_SHORTSTART*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
				const std::string&  /* = NULL*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ) ;



/* ***************************************************************************************************************** */
	
/* ***************************************************************************************************************** */

	ICM_Nthtd(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				// char**  IssuersLabels,
				// double*	IssuersNotionals,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& Frequency /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				const double& FloatingPayerAmount/*  = 0.0*/ ,
				const int& stubrule /* = K_SHORTSTART*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary/*  = CREDIT_DEFAULT_VALUE*/ ,
				const std::string&  /* = NULL*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ) ;
	
	//JLA ICM_Nthtd(const ARM_Date& EffectiveDate,
	//JLA 			const ARM_Date& ScheduleTerminateDate,
	//JLA 			const ARM_Date* FirstPeriodReferenceDate,
	//JLA 			const ARM_Date* FstCpnEffDate,
	//JLA 			const double& FixedRate,
	//JLA 			int intRule,
	//JLA 			int adjStartDate,
	//JLA 			const int& FirstNumDefault,
	//JLA 			const int& LastNumDefault,
	//JLA 			// const int& NbIssuers,
	//JLA 			// const std::vector<std::string> & IssuersLabels,
	//JLA 			// double*	IssuersNotionals,
	//JLA 			const std::vector<std::string>& IssuersLabels,
	//JLA 			const std::vector<std::string>& IssuersNotionals,
	//JLA 			const int& Frequency /* = K_QUARTERLY*/ ,
	//JLA 			const int& DayCountBasis /* = KACTUAL_360*/ ,
	//JLA 			const double& FixedPayerAmount /* = 100.*/ , 
	//JLA 			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
	//JLA 			const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
	//JLA 			const double& FloatingPayerAmount/*  = 0.0*/ ,
	//JLA 			const int& stubrule /* = K_SHORTSTART*/ ,
	//JLA 			const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
	//JLA 			const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
	//JLA 			const double& Binary/*  = CREDIT_DEFAULT_VALUE*/ ,
	//JLA 			const std::string&  /* = NULL*/ ,
	//JLA 			const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ) ;

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const int& FirstNumDefault,
				const int& LastNumDefault,
				// const int& NbIssuers,
				// const std::vector<std::string>& IssuersLabels,
				// double*	IssuersNotionals,
				const std::vector<std::string>& IssuersLabels,
				const std::vector<double>& IssuersNotionals,
				const int& Frequency /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				const double& FloatingPayerAmount /* = 0.0*/ ,
				const int& stubrule /* = K_SHORTSTART*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
				const std::string&  Paycal /* = NULL*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ); 

	
	ICM_Nthtd(ICM_Ftd* ftd);

	virtual ~ICM_Nthtd() {}

	void Init();

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	int GetFirstNumDefault(void) { return itsFirstNumDefault; }
	int GetLastNumDefault(void) { return itsLastNumDefault; }

	void View(char* id, FILE* ficOut)
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

	    fprintf(fOut, "\t\t\t ----------------- Nth to default ----------------- \n");
		fprintf(fOut, "\tFirst Num Default :%f\n", (double)itsFirstNumDefault);
		fprintf(fOut, "\tLast Num Default :%f\n", (double)itsLastNumDefault);

		ICM_Ftd::View(id, fOut);

	    if ( ficOut == NULL )
	    {
		fclose(fOut);
		}
	}

	void SetFirstNumDefault(const int& num) {itsFirstNumDefault = num;}
	void SetLastNumDefault(const int& num) {itsLastNumDefault = num;}

	virtual void RebuildAfterDefault(ICM_ModelMultiCurves* mod);
	virtual void Bounds(double& low,double& up,ARM_Date AsOf) {low=(double)itsFirstNumDefault;up=(double)itsLastNumDefault;}
private:
	ICM_Nthtd(const ICM_Nthtd&) ;//NA
	ICM_Nthtd& operator=(const ICM_Nthtd&); //NA 
};

#endif 
