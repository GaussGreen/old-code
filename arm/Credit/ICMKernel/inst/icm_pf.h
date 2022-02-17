
#if !defined(_ICM_PORTFOLIO_H_)
#define _ICM_PORTFOLIO_H_

#include "ICMKernel\inst\icm_ftd.h"
#include "ICMKernel\util\icm_utils.h"


/*********************************************************************************/
/*! \class  ICM_Portfolio icm_pf.h "icm_pf.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_pf.h
 *	\brief  Description of a Portfolio */
/***********************************************************************************/

class ICM_Parameters; 

class ICM_Portfolio : public ARM_Security  
{
private: 

	int							itsNbSec;
	ICM_Parameters*				itsParams;
	ARM_Security**				itsSecurities;  

private: 
	
	void Init();

public:


	ICM_Portfolio()
	{
		Init();
	}	

/* ***************************************************************************************************************** */
/*!	\fn ICM_Portfolio(ARM_Security** securities,
					  int nbofsecurities)

	\brief Constructor of an <B> ISDA CDS </B>
	\param EffectiveDate Effective Date as in ISDA definition
	\param ScheduleTerminationDate Maturity Date Unadjusted
	\param FixedRate Premiun in unit of the leg
	\param FirstPeriodReferenceDate of the fixed leg if <B> 1 </B> computed from TerminationDate. 
	\param Ccy Currency of the leg: The default currency is a function of the site 
	\param Frequency Fixed Leg Frequency default value is K_QUATERLY
	\param DayCountBasis Fixed Leg Day Count Basis default value is KACTUAL_360
	\param FixedPayerAmount Fixed payer amount default value is 100
	\param AccruedOnDefault Fixed Payer accrued in case of default
	\param Ccy Currency of the Fixed and Floating leg
	\param FloatingPayerAmount: Amount of the floating leg: If <B> 0 </B> it should be equal to the fixed leg amount
	\note This constructor will make a FixedLeg and Floating Leg 
*/	
/* ***************************************************************************************************************** */

	ICM_Portfolio(ARM_Security** securities,
				  int nbofsecurities,
				  ICM_Parameters* params = NULL) : ARM_Security()
	{
		Init();

		Set(securities, nbofsecurities,params);
	}

	void Set(ARM_Security** securities,
				  int nbofsecurities,
				  ICM_Parameters* params = NULL);


	virtual ~ICM_Portfolio() ;


	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	void SetParams(ICM_Parameters* params) ;

	ICM_Parameters* GetParams() { return itsParams;}

	ARM_Object* Clone(void);

	//void View(char* id, FILE* ficOut);

	ARM_Security** GetSecurities()
	{
		return itsSecurities;
	}

	ARM_Security* GetSecurity(int n)
	{
		if (n>=0 && n<itsNbSec)
			return itsSecurities[n];
		else
			return NULL;
	}

	int GetNbSec() {return itsNbSec;};
	void View(char* id, FILE* ficOut) ;
	//void GetIssuersDatas(char** &names,ARM_Vector* &notional, int& nbissuers) ;
	void GetIssuersDatas(std::vector<std::string>& names,ARM_Vector& notional); 

private:
	ICM_Portfolio(const ICM_Portfolio&ref); // NA 
	ICM_Portfolio& operator=(const ICM_Portfolio&ref) ; //NA 
};


#endif 
