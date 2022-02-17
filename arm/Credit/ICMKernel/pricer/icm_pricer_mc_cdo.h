
#ifndef _ICM_PRICER_MC_CDO_H
#define _ICM_PRICER_MC_CDO_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/inst/icm_mez.h"

/*********************************************************************************/
/*! \class  ICM_PRICER_MC_CDO icm_pricer_mc_cdo.h "icm_pricer_mc_cdo.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   February 2007 
 *	\file   icm_pricer_mc_cdo.h
 *	\brief  Monte Carlo Cdo Pricer */
/***********************************************************************************/

typedef double (__stdcall *FUNCPTR)();

static FUNCPTR globalFunc;
void __stdcall Register(long address);

class _TESTGEN{};

class ICM_Pricer_MC_CDO : public ICM_Pricer
{
private:

	long itsnbsimuls;

public :
	ICM_Pricer_MC_CDO(void) 
	{ Init(); }

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters, const ARM_Date&AsOf,long nbsimuls=-1)
	{	
		ICM_Pricer::Set(sec,mod,parameters,&AsOf);

		if (nbsimuls>0) 
		{itsnbsimuls = nbsimuls;}
	}

	~ICM_Pricer_MC_CDO()
	{
	}

	void Init()
	{
		SetName(ICM_PRICER_MC_CDO);
		itsnbsimuls = 100;
	}

	virtual double ComputePrice(qCMPMETH mode);
	double PayOff();

	virtual double ComputeDuration() {return CREDIT_DEFAULT_VALUE;}	
	virtual double ComputeSpread(const double& MtM = 0.) {return CREDIT_DEFAULT_VALUE;}	
private:
	virtual double Accrued() {return CREDIT_DEFAULT_VALUE;}	
	virtual double FeeLegPV () {return CREDIT_DEFAULT_VALUE;}	
	virtual double DefLegPV () {return CREDIT_DEFAULT_VALUE;}	
public:
	virtual	double ComputeImpliedVol(const double& Price) {return CREDIT_DEFAULT_VALUE;}	
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
									const std::string& plot,
									const std::string& label,
									double epsilon, double epsvalueGamma = 0) {return CREDIT_DEFAULT_VALUE;}	

};


#endif

