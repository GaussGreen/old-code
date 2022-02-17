
#ifndef _ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD_H
#define _ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD_H

#include "ICMKernel/pricer/icm_pricer_homogeneous_smile.h"
#include <vector>

/*********************************************************************************/
/*! \class  ICM_Pricer_Distrib_Smile_Collat_Fwd 
 *  \author C. Nicklaus
 *	\version 1.0
 *	\date   Sept 2005
 *	\file   icm_pricer_homogeneous_smile_collat_fwd.h
 *	\brief  Random Factors */
/***********************************************************************************/

class ICM_Pricer_Distrib_Smile_Collat_Fwd : public ICM_Pricer_Distrib_Smile
{

public:

	ICM_Pricer_Distrib_Smile_Collat_Fwd(void) { Init();}

	/** ICM_Pricer_Distrib_Smile_Collat_Fwd(ARM_Security *sec, 
								   ARM_Model *mod, 
								   const ICM_Parameters& parameters ,
								   const ARM_Date&asof)
	{
		Init();
		Set(sec, mod, parameters,asof);
	}
	**/ 
	void Init() 
	{
		SetName(ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD);
	}

	~ICM_Pricer_Distrib_Smile_Collat_Fwd(void) {}

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters&parameters ,const ARM_Date&asof);

	double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
	double ExpectedLossTrancheFullHomogene(const double& yearterm,vector<double>& losses);
	void View(char* id, FILE* ficOut);

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Distrib_Smile_Collat_Fwd* theClone = new ICM_Pricer_Distrib_Smile_Collat_Fwd; 
	theClone->Set(sec,mod,GetParameters(),GetAsOfDate());
    return(theClone);
	}
};

#endif

