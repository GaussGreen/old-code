/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		icm_Homogene_FastLossDistrib_RF.H
	PROJECT:	MOD
	
	DESCRIPTION:	Loss Computation with random factors (mixture copula)

  -----------------------------------------------------------------

 	CREATION:	June, 2005

	LAST MODIF:	June, 2005
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by CN

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef _ICM_Gauss1FLossDistrib_RF_H__
#define _ICM_Gauss1FLossDistrib_RF_H__

#include "ICMKernel\mod\icm_distribution.h"
#include "ICMKernel\util\icm_integrator.h"
#include "ICMKernel\util\icm_CorrelFunction.h"


class ICM_Gauss1FLossDistrib_RF  : public ICM_Distribution  
{

//Members

private:

	// Gauss Legendre Integration
	int					itsCopulaType;
	int					itsIntegrationStep;
	int					itsFreedomDegree;
	qIntegratorChoice	itsIntegrationMethod;

	//Correlation Function
	qCAL_INDEX_CORR_TYPE	itsCorrelationType;
	std::vector<double>		itsCorrelFunctionParam;
	ICM_CorrelFunction*		itsCorrelFunction;

	// Combinaison : Used in Full Homogene pricing
	ICM_QMatrix<double> itsCombinations;

	// Integrator
	ICM_Integrator	TheIntegrator;

	// Pre-computed parameters : including Beta coef (Issuers * integrationStep)
	ICM_QMatrix<double>* its_coeff_a;				 
	ICM_QMatrix<double>* its_coeff_b;				


public:

	//-----  Constructors/destructors

	void Init();

	ICM_Gauss1FLossDistrib_RF()
	{
		Init();
	}	

	~ICM_Gauss1FLossDistrib_RF()
	{
		if (itsCorrelFunction) delete itsCorrelFunction;
		if (its_coeff_a) delete its_coeff_a;
		if (its_coeff_b) delete its_coeff_b;
	}	


	ICM_Gauss1FLossDistrib_RF(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector&  beta,
							const ARM_Vector&  LossRates,
							const std::vector<int>& SortedIndice,
							qCAL_INDEX_CORR_TYPE CorrelType,
							const std::vector<double>& CorrelParam,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);
	
	//Full Homogene constructor
	ICM_Gauss1FLossDistrib_RF(const int& nbnames,
							qCAL_INDEX_CORR_TYPE CorrelType,
							const std::vector<double>& CorrelParam,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);


	
	void Set(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  beta,
			 const ARM_Vector&  LossRates,
			 const std::vector<int>& SortedIndice,
			 qCAL_INDEX_CORR_TYPE CorrelType,
			 const std::vector<double>& CorrelParam,
			 const int& discretizationstep = 20,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);

	//Full Homogene Set
	void Set(const int& nbnames,
			 qCAL_INDEX_CORR_TYPE CorrelType,
			 const std::vector<double>& CorrelParam,
			 const int& discretizationstep = 20,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);
	
	void BitwiseCopy(const ARM_Object* src);

	void Copy(const ARM_Object* src);

	ARM_Object* Clone(void);


	double compute_expectedlosstranche(const double& tranche_up, 
									   const double& tranche_down, 
									   const double& lossunit,
									   //ICM_QMatrix<double>* ShiftMatrix = NULL,
									   const int& Tenor = -1,
									   const int& Issuer = -1);
	
	double compute_expectedlosstrancheFullHomog(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef);

	double compute_expectedlosstrancheFullHomogLegendre(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef);

	double compute_expectedlosstrancheFullHomogHermite(const double& LossUnit,
												const double& tranche_down,
												const double& tranche_up,
												const double& Pdef);

public:

	//----- Assessors
	inline void SetCopulaType(const int& value){itsCopulaType = value;}
	inline void GetCopulaType(int& value){value = itsCopulaType;}

	inline void SetIntegrationMethod(const qIntegratorChoice& value){itsIntegrationMethod = value;}
	inline void GetIntegrationMethod(qIntegratorChoice& value){value = itsIntegrationMethod ;}
	inline qIntegratorChoice GetIntegrationMethod(){return itsIntegrationMethod ;}

	inline void SetIntegrationStep(const int& value){itsIntegrationStep = value;}
	inline void GetIntegrationStep(int& value){value = itsIntegrationStep ;}

	inline void SetFreedomDegree(const int& value){itsFreedomDegree = value;}
	inline void GetFreedomDegree(int& value){value = itsFreedomDegree;}

	//Homogenous Pdef
	inline bool check_homogeneous_Pdef(const std::vector<double>& pdef)
	{ 
		bool res =true;
		for(int k=0;k<pdef.size();k++)
		{
			if(fabs(pdef[k]-pdef[0]) > DB_TOL) res = false;
		}
		return res;
	}

	//----- Utilities
	void compute_barrier();
	void precompute_coeffs(); // barrier computation

	inline void		SetProbCond_Elt(const int& x, const int& y,const int& z, const double& value) {its_ProbCond->SetElt(x, y, z,value);}
	inline double	GetProbCond_Elt(const int& x,const int& y, const int& z) {return its_ProbCond->Elt(x, y, z);}
	
	// 17783 inline void		SetProbCond_Perturb_Elt(const int& x,const int& y,const int& z,const double& value) {its_ProbCond_Perturb->SetElt(x, y, z,value);}
	// 17783 inline double	GetProbCond_Perturb_Elt(const int& x,const int& y,const int& z) {return its_ProbCond_Perturb->Elt(x, y, z);}
	
	inline double	GetCoeff_a(const int& issuer, const int& intstep)	{return (*its_coeff_a)(issuer, intstep);}
	inline double	GetCoeff_b(const int& issuer, const int& intstep)	{return (*its_coeff_b)(issuer, intstep);}
	
	double compute_cond_distribZeroLoss_Gaussian(double x);
	double compute_cond_distrib_Gaussian(double x);
	
	// No Smile
	//virtual void compute_distrib(const int& lup);
	// Error au link
	//virtual void compute_distrib(const int& lup, const int& ldown);

	


	//Long Short CDO
	/*********************************************************************/
	double compute_expectedlosstranche_LongShort(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  // const double& minloss
//												  ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  );
	//Distribution
	virtual void compute_distrib_LongShort(const int& lup, const int& lmin);
	
	void compute_distrib_LongShortLegendre(const int& lup, const int& lmin);
	void compute_distrib_LongShortHermite(const int& lup, const int& lmin);
	
	//Utilities
	void InitProbCondPortSize0(const int& lmin, const int& ind_x);
	void ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname);
	void ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname);	
	
	/*********************************************************************/

	
	// INTEGRATOR
	void	SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, const int& TheStep = 0, const double& lbound = -6., const double& ubound = 6.);

	// View
	void View(char* id, FILE* ficOut);

};
#endif 
