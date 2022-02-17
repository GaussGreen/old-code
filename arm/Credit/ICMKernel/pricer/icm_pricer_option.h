
#ifndef _ICM_Pricer_Option_H
#define _ICM_Pricer_Option_H

#include "ICMKernel\inst\icm_option.h"
#include "ICMKernel\pricer\icm_pricer_index.h"
#include "ARMKernel\crv\volflat.h"
#include "ICMKernel\mod\icm_vasicek_diffusion.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Option icm_pricer_option.h "icm_pricer_option.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   September 2004
 *	\file   ICM_Pricer_Option.h
 *	\brief  Pricer for Spread Options */
/***********************************************************************************/

class ICM_Pricer_Option : public ICM_Pricer
{

public :
	ICM_Pricer_Option(void) { Init();}

	~ICM_Pricer_Option(){}

	void Init();

	ICM_Pricer_Option( ARM_Security *option,  ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);
	void Set( ARM_Security *option,  ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);
	virtual void Reset(void);



	virtual double ComputeBSGreeks (const int& Type);
	virtual void Compute_Fwd_Values(void) {throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <ComputePriceBS> method");}; //const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur);
	//double CptDefaultPL(const std::string&  label, double epsvalue);	
	void View(char* id, FILE* ficOut);

private:

	virtual double Accrued() { ICMTHROW(ERR_INVALID_ARGUMENT,"Accrued : Not implemented"); }
	virtual double FeeLegPV () { ICMTHROW(ERR_INVALID_ARGUMENT,"FeeLegPV : Not implemented"); }
	virtual double DefLegPV () { ICMTHROW(ERR_INVALID_ARGUMENT,"DefLegPV : Not implemented"); }
protected:
	virtual double ComputeDuration(void) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); }
	virtual double ComputeSpread(const double& MtM = 0.) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); }

	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon ,  double epsilonGamma =0 ) =0; 

	virtual double DoPrice(qCMPMETH measure)
	{
		if(!getFlg(measure))
			ComputePriceBS() ;
		return getValue(measure) ;
	}

	virtual double ComputePriceBS()
	{
		throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <ComputePriceBS> method");
	}

	virtual double ComputeImpliedVol(const double& Price);

	
	double GetBSDelta()
	{
		ComputePrice(qCMP_OPT_AN_DELTA);
		return this->getValue(qCMP_OPT_AN_DELTA);
	}

	double GetBSGamma()
	{
		ComputePrice(qCMP_OPT_AN_GAMMA);
		return this->getValue(qCMP_OPT_AN_GAMMA);
	}
	 
	double GetBSVega()
	{
		ComputePrice(qCMP_OPT_AN_VEGA);
		return this->getValue(qCMP_OPT_AN_VEGA);
	}

	double GetBSTheta()
	{
		double sensi = ComputeSensitivity(ICM_THETA_Type,"NONE","NONE",0.01);
		return sensi;
	}

	double GetBSRho()
	{
		
		double sensi = ComputeSensitivity(ICMIRCURVE_TYPE,"NONE","NONE",0.01);
		return sensi;
	}

	void	SetDiffusionType(qDIFFUSION_TYPE& type) {its_diffusion_type=type;}
	qDIFFUSION_TYPE	GetDiffusionType() {return its_diffusion_type;}

	void	SetDiffusion(ICM_Diffusion* value) {its_diffusion=value;}
	ICM_Diffusion*	GetDiffusion() {return its_diffusion;}

private :

	qDIFFUSION_TYPE	its_diffusion_type;
	ICM_Diffusion*	its_diffusion;
};

#endif

