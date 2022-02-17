
#ifndef _VASICEK_DIFFUSION_MODEL_H
#define _VASICEK_DIFFUSION_MODEL_H

#include "ARMKernel\mod\model.h"
#include "ARMKernel\crv\volcurv.h"
#include "ICMKernel\crv\icm_defaultcurve.h"

/*********************************************************************************/
/*! \class  ICM_VasicekDiffusion icm_defcurvemodel.h "icm_defcurvemodel.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   May 2004
 *	\brief  Vasicek diffusion */
/***********************************************************************************/

class ICM_Diffusion{};
class ICM_VasicekDiffusion : public ICM_Diffusion
{        
private :
	ARM_Vector	its_yf_schedule;
	double	its_MeanReversion;
	double	its_LongRate;
	ARM_Vector	its_Volatility;
	double its_R0;

public:

	//YF_SCHEDULE_VASICEK	schedule in year fractions for volatility
	//LONG_RATE_VASICEK		long rate 
	//MEANREV_VASICEK		mean reversion
	//VOLATILITY_VASICEK	volatility vector
	//R0_VASICEK	volatility vector
	ICM_VasicekDiffusion(const ICM_Parameters& params)
	{ 
		Init(); 

		ARM_Vector* YF_SCHEDULE_VASICEK = params.GetColVect("YF_SCHEDULE_VASICEK");
		if (YF_SCHEDULE_VASICEK) {its_yf_schedule = *YF_SCHEDULE_VASICEK;}

		ARM_Vector* VOLATILITY_VASICEK = params.GetColVect("VOLATILITY_VASICEK");
		if (VOLATILITY_VASICEK) {its_Volatility = *VOLATILITY_VASICEK;}

		ARM_Vector* LONG_RATE_VASICEK = params.GetColVect("LONG_RATE_VASICEK");
		if (LONG_RATE_VASICEK) {its_LongRate = LONG_RATE_VASICEK->Elt(0);}

		ARM_Vector* MEANREV_VASICEK = params.GetColVect("MEANREV_VASICEK");
		if (MEANREV_VASICEK) {its_MeanReversion = MEANREV_VASICEK->Elt(0);}

		ARM_Vector* R0_VASICEK = params.GetColVect("R0_VASICEK");
		if (R0_VASICEK) {its_R0 = R0_VASICEK->Elt(0);}
	}

	ICM_VasicekDiffusion(ARM_Vector& yf_schedule,
						double& MeanReversion,
						double& LongRate,
						ARM_Vector& Volatility,
						double R0)
	{ 
		Init(); 

		its_yf_schedule=yf_schedule;
		its_MeanReversion=MeanReversion;
		its_LongRate=LongRate;
		its_Volatility=Volatility;	
		its_R0 = R0;
	}

	void Init()
	{
		its_MeanReversion=0.;
		its_LongRate=0.;
		its_R0=0.;
	}

	ICM_VasicekDiffusion(const ICM_VasicekDiffusion& in)
	{
		its_yf_schedule=in.its_yf_schedule;
		its_MeanReversion=in.its_MeanReversion;
		its_LongRate=in.its_LongRate;
		its_Volatility=in.its_Volatility;
		its_R0=in.its_R0;
	}

	double FwdZeroCoupon(double t, double T, double Rt);

};

#endif
