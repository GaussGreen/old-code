
#ifndef _ICM_Pricer_IndexOption_H
#define _ICM_Pricer_IndexOption_H

#include "ICMKernel\pricer\icm_pricer_cdsoption.h"


/*********************************************************************************/
/*! \class  ICM_Pricer_IndexOption icm_pricer_indexoption.h "icm_pricer_indexoption.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   October 2004
 *	\file   ICM_Pricer_IndexOption.h
 *	\brief  Pricer for Index Options */
/***********************************************************************************/

class ICM_Pricer_IndexOption : public ICM_Pricer_CdsOption
{
private:


public :
	ICM_Pricer_IndexOption(void) { Init();}

	~ICM_Pricer_IndexOption(){}

	void Init()
	{
		SetName(ICM_PRICER_INDEXOPTION);
		this->unsetFlgs();
	}

	void ICM_Pricer_IndexOptionvoid(ARM_Security *option, ARM_Model *mod,const ICM_Parameters&params,
								const ARM_Date&asof);
	void Set(ARM_Security *option, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);

	virtual void Reset(void) ;


	double ComputePriceBS();
	double ComputePriceBS_LN_SPREADS();
	double ComputePriceBS_OU();
	double Payoff_OU(double x,void* params);

	//virtual double ComputeImpliedVol(double Price);
	virtual void Compute_Fwd_Value(void); //const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur);


	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon,  double epsilonGamma = 0); 
	double CptDefaultPL(const std::string&  label, double epsvalue);
};

#endif

