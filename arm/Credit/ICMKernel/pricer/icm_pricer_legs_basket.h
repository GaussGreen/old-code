
#ifndef _ICM_PRICER_LEGS_BASKET_H
#define _ICM_PRICER_LEGS_BASKET_H

#include "ICMKernel\pricer\icm_pricer_adviser.h"
// #include "ICMKernel\inst\icm_leg_basket.h"
// #include "ICMKernel\mod\icm_meta_model.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Legs_Basket icm_pricer_legs_basket.h "icm_pricer_legs_basket.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   Febuary 2005
 *	\brief  Pricer for Generic Basket of Legs */
/***********************************************************************************/

class ICM_DefaultCurveModel; 

class ICM_Pricer_Legs_Basket : public ICM_Pricer
{
private:

	vector<ARM_Security*>			itsSecurities;
	vector<ICM_DefaultCurveModel*>	itsModels;
	vector<ARM_CLASS_NAME>			itsPricerType;
	vector<ICM_Pricer*>				itsPricers;

public :
	ICM_Pricer_Legs_Basket(void) { Init();}

	ICM_Pricer_Legs_Basket(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
	{
		Init();

		Set(sec, mod,parameters,asof);
	}

	~ICM_Pricer_Legs_Basket()
	{
		itsSecurities.clear();
		itsModels.clear();
		itsPricerType.clear();

		for (int i=0; i<itsPricers.size(); i++) 
		{
			if (itsPricers[i])delete itsPricers[i];
			itsPricers[i] = NULL;
		}
		
		itsPricers.clear();
	}

	void Init()
	{
		SetName(ICM_DIFFUSIONPRICER);

		itsSecurities.clear();
		itsModels.clear();
		itsPricerType.clear();
		itsPricers.clear();
	}

	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof);


	virtual double ComputePrice(qCMPMETH);

	virtual void Reset(void)
	{
		ICM_Pricer::Reset();
	}

	void View(char* id, FILE* ficOut){}

protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma = 0) ; 
private:
	double Accrued(void) ; 
	double FeeLegPV(void) ;
	double DefLegPV(void) ;
	double Compute_Fwd_Spread(const ARM_Date &,const ARM_Date &, double& dur); 
	double ComputeDuration(void) ;
	double ComputeSpread(const double &) ;
	double ComputeImpliedVol(const double &) ;
};

#endif

