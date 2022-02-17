
#pragma warning( disable : 4150 ) 
#ifndef _ICM_Pricer_Analytic_CDO2_STRIKE_H
#define _ICM_Pricer_Analytic_CDO2_STRIKE_H

#include "ICMKernel/pricer/icm_pricer_homogeneous_smile.h"
#include "ICMKernel/glob/icm_smile_correlation.h"
#include "ICMKernel\glob\icm_correlation_sector.h"


/*********************************************************************************/
/*! \class  ICM_Pricer_Analytic_Cdo2 icm_pricer_analytic_cdo2.h "icm_pricer_analytic_cdo2.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_pricer_analytic_cdo2.h
 *	\brief  Analytic pricing of CDO² */
/***********************************************************************************/

class ICM_CorrMatrix; 
class ICM_Pricer_Analytic_Cdo2_Smile : public ICM_Pricer_Distrib_Smile
{
private:

	ICM_ModelMultiCurves* itsIntermediateModel0; //Initial Intermediate Model	
	ICM_ModelMultiCurves* itsIntermediateModel;	 //Intermediate Model for sensitivities

	ARM_Security*		  itsIntermediateSecurity; //Intermediate Security			

	vector<ICM_Smile_Correlation*> itsImplicitUnderlyingCorrelations; //Implicit Underlying correlations
	ICM_Correlation*	  itsForcedCorrelationForCdo2; //Forced Correlation

	int					  itsIntegrationStepUnderlying; //Integration step for underlying
	int					  itsIntegrationMethodUnderlying; //Integration method for underlying
	int					  itsNumSimulationsForCorrelation; //Num Simulation step for underlying
	int					  itsForcedRescalingForUnderlyings; //Forced rescaling for underlyings

	double				  itsAverageCorrTpsDefault_Up;
	double				  itsAverageCorrTpsDefault_Down;	
	bool				  itsAverageCorrTpsDefaultFlg;
	double				  itsAverageCorrTpsDefault_positive_bump_Up;
	double				  itsAverageCorrTpsDefault_positive_bump_Down;
	double				  itsAverageCorrTpsDefault_negative_bump_Up;
	double				  itsAverageCorrTpsDefault_negative_bump_Down;

	// Sectorial Correlation
	qCORRELATION_FIT_TYPE	its_Correlation_Fit_Type;
	//ICM_Correlation_Sector*		its_Sectorial_Correlation;

	ICM_CorrMatrix*			itsCorrelTpsDefaut ;
private:
	void setCorrelTpsDefaut(const ICM_CorrMatrix*); 
	
public :
	ICM_Pricer_Analytic_Cdo2_Smile(void) { Init();}

	/**ICM_Pricer_Analytic_Cdo2_Smile(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters ,const ARM_Date&asof)
	{
		Init();

		Set(sec, mod, parameters,asof);
	}
	**/ 
	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters ,const ARM_Date&asof);

	virtual void BeforePrice(const std::string& label  ,qSENSITIVITY_TYPE type /** = ICMSPREAD_TYPE **/ ) ;

	void Init() ;
	 

	virtual void Reset(void) ; 

	~ICM_Pricer_Analytic_Cdo2_Smile(void) ;

	void SetForcedCorrelationForCdo2(ICM_Correlation*  correl); // adoption : can be null. 
	void SetImplicitUnderlyingCorrelations(const std::vector<ICM_Smile_Correlation*>&ref) ; // adoption

	ICM_Correlation* GetForcedCorrelationForCdo2() {return itsForcedCorrelationForCdo2;}
	
	bool GetAverageCorrTpsDefaultFlg() const { return itsAverageCorrTpsDefaultFlg;}
	double GetAverageCorrTpsDefault_Up() const {return itsAverageCorrTpsDefault_Up;}
	double GetAverageCorrTpsDefault_Down() const {return itsAverageCorrTpsDefault_Down;}
	
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

	void View(char* id, FILE* ficOut);

	virtual double ComputePrice(qCMPMETH measure);
	virtual double DoPrice(qCMPMETH measure);
	virtual void ComputeGetAverageCorrTpsDefault();

	void ResetTranchesModel(const std::string& label ,qSENSITIVITY_TYPE type /** = ICMSPREAD_TYPE **/ );
	void CptProportionsIndexName(int& sizeprop,ARM_Vector**& Proportions_IndexName,std::vector<std::string>& IndexName);
	void CptForcedCorrelationForCdo2(const ICM_Parameters& parameters);
	void ComputeImplicitUnderlyingCorrelations(vector<ICM_Smile_Correlation*>& UnderlyingCorrelations,const bool& Reset=true);
	double ComputeImplicitSpread(const string& name, const string& tenor);

	virtual	void	GetDataFromLabel(string TheLabel, double& TheValue) 
	{
		Price(qCMPPRICE);
		if (TheLabel=="AVGCORRDEF") TheValue = itsAverageCorrTpsDefault_Down;
		if (TheLabel=="AVGCORRDEF_UP") TheValue = itsAverageCorrTpsDefault_Up;
		if (TheLabel=="AVGCORRDEF_DOWN") TheValue = itsAverageCorrTpsDefault_Down;
		else if (TheLabel=="AVGCORRDEF_POS_BUMP_UP") TheValue = itsAverageCorrTpsDefault_positive_bump_Up;
		else if (TheLabel=="AVGCORRDEF_POS_BUMP_DOWN") TheValue = itsAverageCorrTpsDefault_positive_bump_Down;
		else if (TheLabel=="AVGCORRDEF_NEG_BUMP_UP") TheValue = itsAverageCorrTpsDefault_negative_bump_Up;
		else if (TheLabel=="AVGCORRDEF_NEG_BUMP_DOWN") TheValue = itsAverageCorrTpsDefault_negative_bump_Down;
	}

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Analytic_Cdo2_Smile* theClone = new ICM_Pricer_Analytic_Cdo2_Smile; 
	theClone->Set(sec,mod,GetParameters(),GetAsOfDate());
    return(theClone);
	}

	virtual double ComputeSpread(const double& MtM = 0.);
	virtual double ComputeDuration(void);

	void	Set_Correlation_Fit_Type(const qCORRELATION_FIT_TYPE& value){its_Correlation_Fit_Type = value;}
	void	Get_Correlation_Fit_Type(qCORRELATION_FIT_TYPE& value){value = its_Correlation_Fit_Type;}

	virtual void MarketData(ARM_Security* sec,vector<string>& DMKT) ;

protected:
	virtual void ResetLU(void)  
	{
		ResetLossUnit() ;
	}

};



#endif

