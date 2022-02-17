
#ifndef _ICM_PRICER_IR_INFCRV_H
#define _ICM_PRICER_IR_INFCRV_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "gpinflation/infleg.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Security icm_pricer_ir_infcurve.h "icm_pricer_ir_infcurve.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   May 2005 
 *	\file   icm_pricer_ir_infcurve.h
 *	\brief  Pricer for all security using ICM_Security */
/***********************************************************************************/


class ICM_Pricer_IR_InfCrv : public ICM_Pricer
{
private:

public :
	ICM_Pricer_IR_InfCrv(void) ;
	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& params,const ARM_Date& asof) ;
	virtual ~ICM_Pricer_IR_InfCrv(){}

	void Init() ;
	

	
	virtual void Reset(void) {}
	virtual void View(char* id, FILE* ficOut) {}

	virtual double ComputePrice(qCMPMETH measure) ;


	virtual void GenerateRates(ARM_Vector& rates) ;
protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon,  double epsilonGamma = 0)  ;
 

private:
	double Accrued(void)  
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double FeeLegPV(void)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double DefLegPV(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double Compute_Fwd_Spread(const  ARM_Date &,const  ARM_Date &, double& dur) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double ComputeDuration(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double ComputeSpread(const double &) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double ComputeImpliedVol(const double &) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}

};


#endif

