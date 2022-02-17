
#ifndef _ICM_Pricer_Analytic_CDO2_H
#define _ICM_Pricer_Analytic_CDO2_H

#include "ICMKernel/pricer/icm_pricer_homogeneous.h"
// #include "ICMKernel/inst/icm_mez.h"


/*********************************************************************************/
/*! \class  ICM_Pricer_Analytic_Cdo2 icm_pricer_analytic_cdo2.h "icm_pricer_analytic_cdo2.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  Analytic pricing of CDO² */
/***********************************************************************************/

class ICM_Pricer_Analytic_Cdo2 : public ICM_Pricer_Distrib
{
private:

	ICM_ModelMultiCurves* itsIntermediateModel0; //Initial Intermediate Model	
	ICM_ModelMultiCurves* itsIntermediateModel;	 //Intermediate Model for sensitivities

	ARM_Security*		  itsIntermediateSecurity0; //Initial Intermediate Security			
	ARM_Security*		  itsIntermediateSecurity; //Intermediate Security			

public :
	ICM_Pricer_Analytic_Cdo2(void) { Init();}

	/** ICM_Pricer_Analytic_Cdo2(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters,const ARM_Date&asof)
	{
		Init();

		Set(sec, mod, parameters,asof);
	}
	**/ 
	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof);

	virtual void BeforePrice(const std::string& label,qSENSITIVITY_TYPE type /** = ICMSPREAD_TYPE**/ ) ; 

	void Init() ; 


	virtual void Reset(void); 

	~ICM_Pricer_Analytic_Cdo2(void) ;
	 

	ICM_DefaultCurve* CptDefaultCurve(const ARM_Vector& VectorYF,const ARM_Vector& spreads,ICM_ModelMultiCurves* model,const std::string&label);
	void CptTranchesModel();

	ICM_ModelMultiCurves* GetIntermediateModel() {return itsIntermediateModel;}	
	void SetIntermediateModel(ICM_ModelMultiCurves* model) 
	{
		itsIntermediateModel = model;
	}	

	ARM_Security* GetIntermediateSecurity() { return itsIntermediateSecurity;}			
	void SetIntermediateSecurity(ARM_Security* security) 
	{ 
		itsIntermediateSecurity=security;
	}	

	ICM_ModelMultiCurves* GetIntermediateModel0() {return itsIntermediateModel0;}	
	void SetIntermediateModel0(ICM_ModelMultiCurves* model) 
	{
		itsIntermediateModel0 = model;
	}	

	ARM_Security* GetIntermediateSecurity0() { return itsIntermediateSecurity0;}			
	void SetIntermediateSecurity0(ARM_Security* security) 
	{ 
		itsIntermediateSecurity0=security;
	}	

	void View(char* id, FILE* ficOut);

	virtual double ComputePrice(qCMPMETH measure/**const ARM_Date& ExecutionDate,
								const int& mode=0, 
								const int& settleGap=0**/ );

	void ResetTranchesModel(const std::string& label);

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Analytic_Cdo2* theClone = new ICM_Pricer_Analytic_Cdo2 ;
	theClone->Set(sec,mod,GetParameters(),GetAsOfDate());
    return(theClone);
	}

};



#endif

