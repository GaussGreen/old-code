
#ifndef _ICM_PRICER_ADVISOR_H
#define _ICM_PRICER_ADVISOR_H


/*********************************************************************************/
/*! \class  ICM_Pricer_Advisor icm_pricer_adviser.h "icm_pricer_adviser.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004
 *	\brief  Pricer Adviser */
/***********************************************************************************/

#include <ICMKernel/pricer/icm_pricer.h>
#include <set>

class ICM_PRICER_TYPE
{
	public :
		ARM_CLASS_NAME m_security_type;
		ARM_CLASS_NAME m_pricer_type;
	public :
		ICM_PRICER_TYPE(const ARM_CLASS_NAME& security_type,
						const ARM_CLASS_NAME& pricer_type = ARM_OBJECT):m_security_type(security_type),m_pricer_type(pricer_type){}
		ICM_PRICER_TYPE(){m_pricer_type=ARM_OBJECT;}
		bool operator < (const ICM_PRICER_TYPE & rhs) const {return ((long)m_security_type < (long)rhs.m_security_type);}
		void Set(const ARM_CLASS_NAME& security_type,
				const ARM_CLASS_NAME& pricer_type=ARM_OBJECT){m_security_type = security_type; m_pricer_type =pricer_type;}
} ;

class ICM_Pricer_Advisor : public ARM_Object
{

private:
	std::set<ICM_PRICER_TYPE> itsMapAdvisor;

	
public :
	ICM_Pricer_Advisor(void) { Init();}

	void Init()  ;


	inline ARM_CLASS_NAME GetPricerType(ARM_Security* sec) ;

	virtual ICM_Pricer* GeneratePricer(ARM_Security* sec, 
										ARM_Object * mod, //JLA: not really a model !
										ARM_CLASS_NAME InputPricerType /** = ARM_OBJECT **/ ,
										int nbpaths /** = CREDIT_DEFAULT_VALUE **/ ,
										const ICM_Parameters* parameters ,
										// JLA: useless ICM_MktDataMng* MktDataManager ,
										const ARM_Date& asof ) ; 
   

	virtual void CheckForForwardPricer(ICM_Pricer*& initpricer) ;


};
#endif
