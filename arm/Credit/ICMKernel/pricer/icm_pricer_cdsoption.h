
#ifndef _ICM_Pricer_CdsOption_H
#define _ICM_Pricer_CdsOption_H

#include "ICMKernel/pricer/icm_pricer_option.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_CdsOption icm_pricer_cdsoption.h "icm_pricer_cdsoption.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   October 2004
 *	\file   ICM_Pricer_CdsOption.h
 *	\brief  Pricer for Cds Options */
/***********************************************************************************/

class ICM_Pricer_CdsOption : public ICM_Pricer_Option
{
	
public :
	ICM_Pricer_CdsOption(void) { Init();}
	~ICM_Pricer_CdsOption(){}

	void Init()	{
		SetName(ICM_PRICER_CDSOPTION);
	}
	ICM_Pricer_CdsOption( ARM_Security *option,  ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);
	void Set( ARM_Security *option,  ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);
	virtual void Reset(void)
	{
		ICM_Pricer_Option::Reset();
	}

	virtual double ComputePriceBS();
	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 

	public:
	virtual void Compute_Fwd_Values(void);//const ARM_Date &Mty, const ARM_Date & CDS_ExpiryDate,double& dur);
	double CptDefaultPL(const std::string&  label, double epsvalue);
	
};

#endif

