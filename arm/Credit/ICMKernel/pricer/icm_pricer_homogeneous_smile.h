
#ifndef _ICM_Pricer_HOMOGENEOUS_SMILE_H
#define _ICM_Pricer_HOMOGENEOUS_SMILE_H

#include "ICMKernel/pricer/icm_pricer_homogeneous.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Distrib icm_pricer_homogeneous.h "icm_pricer_homogeneous.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  ce pricer est basé sur l'utilisation d'un cache de distribution */
/***********************************************************************************/

class ICM_Distribution ;

class ICM_Pricer_Distrib_Smile : public ICM_Pricer_Distrib
{

private:
		
	int		itsCopulaType;					  //Copula Type	
	int		itsIntegrationMethod;			  //Hermite, Legendre	
	int		itsFreedomDegree;				  //FreedomDegree	
	int		itsIntegrationStep1;			  //integration step	
	qRescalType	itsRescalType;				  //Rescaling Type	
	int     itsApproxComputation;			  //fast distribution 

	ICM_Distribution*	itsDistribution;

public:
	
	int		GetCopulaType() const { return itsCopulaType; }
	void	GetCopulaType(int& data)	{data = itsCopulaType;}
	void	SetCopulaType(const int& data)	{itsCopulaType = data;}

	void	GetFreedomDegree(int& data)	{data = itsFreedomDegree;}
	void	SetFreedomDegree(const int& data)	{itsFreedomDegree = data;}

	void	GetIntegrationMethod(int& data)	{data = itsIntegrationMethod;}
	void	SetIntegrationMethod(const int& data)	{itsIntegrationMethod = data;}

	double	GetApproxComputation()	{return itsApproxComputation;}
	void	SetApproxComputation(const double& data)	{itsApproxComputation = data;}

	int		GetIntegrationStep1()	const { return itsIntegrationStep1;}
	void	GetIntegrationStep1(int& data)	{data = itsIntegrationStep1;}
	void	SetIntegrationStep1(const int& data)	{itsIntegrationStep1 = data;}

	ICM_Distribution* GetDistribution(void)	{return itsDistribution;}
	void	SetDistribution(ICM_Distribution* distrib)	;

public:

	ICM_Pricer_Distrib_Smile(void) { Init();}

	/** ICM_Pricer_Distrib_Smile(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters ,const ARM_Date&asof)
	{
		Init();

		Set(sec, mod, parameters,asof);
	}
	**/ 
	void Init() ;
	

	virtual ~ICM_Pricer_Distrib_Smile(void) ;

	void Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters ,const ARM_Date&asof);

	double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
	// double ExpectedLossTrancheFullHomogeneous(const double& yearterm,vector<double>& losses);
// 17783 	double ExpectedLossTrancheForFastSpreadHedge(const double& yearterm);

	void	ExtractParameters();
	qIntegratorChoice GetGaussLegendreMethodFromIntegrationStep(const int& TheIntegrationStep);
	inline qRescalType GetRescalType() {return itsRescalType;}
	inline void SetRescalType(const qRescalType& value) {itsRescalType=value;}
	
	void View(char* id, FILE* ficOut);

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Distrib_Smile* theClone = new ICM_Pricer_Distrib_Smile; 
	theClone->Set(sec,mod,GetParameters(),GetAsOfDate());
    return(theClone);
	}


};



#endif

