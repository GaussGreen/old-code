
#ifndef _ICM_PRICER_IR_YCMOD_H
#define _ICM_PRICER_IR_YCMOD_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "ARMKernel/pricer/ipricer.h"
#include "gpcalculators\pricerfactory.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_IR_YcMod icm_pricer_ir_ycmod.h "icm_pricer_ir_ycmod.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   May 2005 
 *	\file   icm_pricer_ir_ycmod.h
 *	\brief  Pricer for all security using ICM_Security */
/***********************************************************************************/
// using ARM::ARM_PricerFactory;

class ICM_Pricer_IR_YcMod : public ICM_Pricer
{
private:

public :
	ICM_Pricer_IR_YcMod(void) ;

 
	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);
	;

	~ICM_Pricer_IR_YcMod();

	void Init();
	

	
	virtual void Reset(void) ;
	void View(char* id, FILE* ficOut) ;

	virtual void GenerateRates(ARM_Vector& rates); 

	virtual double ComputeSpread(const double& MtM = 0.); 
	


	virtual double ComputeDuration(void);
  
private:
	virtual double FeeLegPV () ;
    virtual double DefLegPV () ;
	virtual double Accrued() ;
public:
	 

	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
										const std::string& plot, 
										const std::string& label , 
										double  epsvalue ,  double epsilonGamma = 0);


private:
	double Compute_Fwd_Spread(const  ARM_Date &,const  ARM_Date &, double& dur)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
	double ComputeImpliedVol(const double &) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
};


#endif


