
//////////////////////////////////////////////////////////////////////
// icm_leg_basket.h: interface for the ICM_Leg_Basket class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ICM_LEG_BASKET_H_)
#define _ICM_LEG_BASKET_H_

#include "ICMKernel\inst\icm_pf.h"


/*********************************************************************************/
/*! \class  ICM_Leg_Basket icm_leg_basket.h "icm_leg_basket.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   Febuary 2005
 *	\brief  Description of Basket of legs */
/***********************************************************************************/

class ICM_Leg_Basket : public ICM_Portfolio
{


private: 
	

public:

	void Init() {}

	ICM_Leg_Basket() {}


	ICM_Leg_Basket(ARM_Security** securities,
				  int nbofsecurities,
				  ICM_Parameters* params = NULL) : ICM_Portfolio(securities,
															 nbofsecurities,
															 params)
	{
		Init();

		Set(securities, nbofsecurities,params);
	}

	void Set(ARM_Security** securities,
				  int nbofsecurities,
				  ICM_Parameters* matrix = NULL)
	{
	}


	virtual ~ICM_Leg_Basket()
	{
	}	

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg){}

	ARM_Object* Clone(void){return NULL;}

};


#endif 
