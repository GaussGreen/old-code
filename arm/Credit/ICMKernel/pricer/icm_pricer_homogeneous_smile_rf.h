
#ifndef _ICM_PRICER_HOMOGENEOUS_SMILE_RANDOM_FACTORS_H
#define _ICM_PRICER_HOMOGENEOUS_SMILE_RANDOM_FACTORS_H

#include "ICMKernel/pricer/icm_pricer_homogeneous_smile.h"
#include <vector>

/*********************************************************************************/
/*! \class  ICM_Pricer_Distrib_RandFactors icm_pricer_homogeneous_smile_rf.h "icm_pricer_homogeneous_smile_rf.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_pricer_homogeneous_smile_rf.h
 *	\brief  Random Factors */
/***********************************************************************************/

class ICM_Pricer_Distrib_RandFactors : public ICM_Pricer_Distrib_Smile
{

private:
	std::vector<double> itsRFParameters; //Parameters for distribution
	qCAL_INDEX_CORR_TYPE itsModelType; //Model Type

public:
	
	std::vector<double>& GetRFParameters()	{return itsRFParameters;}
	void	SetRFParameters(const std::vector<double>& data)	{itsRFParameters = data;}

public:

	ICM_Pricer_Distrib_RandFactors(void) { Init();}

	/** ICM_Pricer_Distrib_RandFactors(ARM_Security *sec, 
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
		SetName(ICM_PRICER_HOMOGENEOUS_SMILE_RF);
		itsRFParameters.clear();
	}

	~ICM_Pricer_Distrib_RandFactors(void) {}

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters ,const ARM_Date&asof);

	double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
	double ExpectedLossTrancheFullHomogene(const double& yearterm);
	void View(char* id, FILE* ficOut);

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Distrib_RandFactors* theClone = new ICM_Pricer_Distrib_RandFactors; 
	theClone ->Set(sec,mod,GetParameters(),GetAsOfDate());
    return(theClone);
	}
};

#endif

