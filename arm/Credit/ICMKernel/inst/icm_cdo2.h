
#if !defined(_ICM_CDO2_H_)
#define _ICM_CDO2_H_

#include "ICMKernel\inst\icm_pf.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel/util/icm_qmatrix.h"

/*********************************************************************************/
/*! \class  ICM_Cdo2 icm_cdo2.h "icm_cdo2.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   April 2004
 *	\file   icm_cdo2.h
 *	\brief  <B> Cdo Square </B> */
/***********************************************************************************/



class ICM_Cdo2 : public ICM_Mez  
{
private: 

	ICM_Portfolio* itsPortfolio;			//Portfolio of Underlyings Cdos	
	bool		   itsIsCrossSubordinate;   //For Cdo² with Cross Subordination

public:

	ICM_Cdo2()
	{Init();}	

	ICM_Cdo2(ICM_Cds* cds,
			double SubAmount,
			ICM_Portfolio* pf,
			const double& Binary = CREDIT_DEFAULT_VALUE) ; 

	ICM_Cdo2(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate ,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const int& Frequency /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string & Ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				const int& stubrule /* = K_SHORTSTART*/ ,
				ICM_Portfolio* ptf /* = NULL*/ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FreqDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
				const std::string &  payCalName /* = NULL*/ ,
				const bool& CrossSub /* = false*/ ,
				const bool& IncludeMaturity /* = INCLUDE_MATURITY */ );

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate ,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const int& Frequency /*= K_QUARTERLY*/ ,
				const int& DayCountBasis /*= KACTUAL_360*/ ,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED*/ ,
				const std::string & Ccy /*= ARM_DEFAULT_CURRENCY*/ , 
				const int& stubrule /*= K_SHORTSTART*/ ,
				ICM_Portfolio* ptf /*= NULL*/ ,
				const double& CreditLag /*= DEFAULT_CREDIT_LAG*/ ,
				const int& FreqDefLeg /*= DEFAULT_FRQ_DEFLEG*/ ,
				const double& Binary /*= CREDIT_DEFAULT_VALUE*/ ,
				const std::string &  payCalName /*= NULL*/ ,
				const bool& CrossSub /*= false*/ ,
				const bool& IncludeMaturity /*= INCLUDE_MATURITY*/ );


// ------------------------------------------------------------------------
	virtual ~ICM_Cdo2() 
	{
		if (itsPortfolio)
			delete itsPortfolio;
		itsPortfolio = NULL;
	}

	void Init();

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	ICM_Portfolio* GetPortfolio(void) {return itsPortfolio;}
	void SetPortfolio(ICM_Portfolio* Portfolio)	
	{
		itsPortfolio = Portfolio;
		CptIssuersData();
	}

	void View(char* id, FILE* ficOut);

	void CptIssuersData(void) ;

	ICM_QMatrix<int>* GenAppMatrix(std::vector<std::string>&_UnionIssuers);
	ICM_QMatrix<int>* AppMatrix(const std::vector<std::string>& vUnionIssuers);
	bool IsCrossSubordinate() { return itsIsCrossSubordinate; }
private:
	ICM_Cdo2(const ICM_Cdo2&ref) ; // NA 
	ICM_Cdo2& operator=(const ICM_Cdo2&ref); // NA 
};

#endif 
