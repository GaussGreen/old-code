/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		icm_Homogene_FastLossDistrib_Gaussian_MF.H
	PROJECT:	MOD
	
	DESCRIPTION:	Loss Computation with Smile

  -----------------------------------------------------------------

 	CREATION:	May, 2006

	LAST MODIF:	May, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Christophe Nicklaus and Valentin Rodriguez

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef __ICM_Gaussian_LossDistrib_MF_H__
#define __ICM_Gaussian_LossDistrib_MF_H__

#include <string>
#include <vector>
#include <deque>
#include "ICMKernel\mod\icm_distribution.h"
#include "ICMKernel\util\icm_integrator.h"


class ICM_Gaussian_LossDistrib_MF  : public ICM_Distribution  
{

	//-----  Constructors/destructors
private:

	// Gauss Legendre Integration
	int		itsCopulaType;
	int		itsFreedomDegree;
	int		itsNbStrikes;
	int		itsNbSectors;
	vector<int> itsIndicSectors;

	ICM_QMatrix<double>*	itsBetasMatrix;

	//Type de corrélation
	qTWO_FACTORS_CORRELATION_TYPE	its_CorrelationType;

	// Integrator
	ICM_Integrator	TheIntegrator;
	qIntegratorChoice	itsIntegrationMethod;
	int		itsIntegrationStep;

	// Prise en compte du Smile : structure dupliquées pour le strike Down
	double	its_taildistrib_Down;			//Queue de distribution
	vector<double> its_beta_Down;			//Beta
	vector<double> its_lambda;			//Beta
	vector<double> its_lambda_Down;			//Beta
	vector<double> its_lossdistrib;	//For assessor, m_dens[k] = P(N(T)=k) 
	vector<double> its_lossdistrib_Down;	//For assessor, m_dens[k] = P(N(T)=k) 
	ICM_QCubix<double>* its_ProbCond;	//Stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QCubix<double>* its_ProbCond_Down;	//Stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QCubix<double>* its_SectorLoss;		//Stockage distribution de pertes conditionnelles au facteur systémique pour tous les secteur
	ICM_QCubix<double>* its_SectorLoss_Down;	//Stockage distribution de pertes conditionnelles au facteur systémique pour tous les secteur
	
	// Valeurs des coef stockés pour le calcul des prob conditionnelles de défaut
	vector<double>	its_coeff_a;
	vector<double>	its_coeff_b;
	vector<double>	its_coeff_c;
	vector<double>	its_coeff_a_down;
	vector<double>	its_coeff_b_down;
	vector<double>	its_coeff_c_down;

	vector<double>	itsStrikes;
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/	itsIsShortName;
	
public:

	void Init();

	ICM_Gaussian_LossDistrib_MF()
	{
		Init();
	}	

	~ICM_Gaussian_LossDistrib_MF()
	{
		if (itsBetasMatrix)
			delete itsBetasMatrix;
		itsBetasMatrix = NULL;

		if (its_ProbCond)
			delete its_ProbCond;
		its_ProbCond = NULL;

		if (its_ProbCond_Down)
			delete its_ProbCond_Down;
		its_ProbCond_Down = NULL;

		if (its_SectorLoss)
			delete its_SectorLoss;
		its_SectorLoss = NULL;

		if (its_SectorLoss_Down)
			delete its_SectorLoss_Down;
		its_SectorLoss_Down = NULL;

		itsIndicSectors.clear();

		its_beta_Down.clear();		 
		its_lambda.clear();
		its_lambda_Down.clear();
		its_lossdistrib.clear(); 
		its_lossdistrib_Down.clear(); 

		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_c.clear();
		its_coeff_a_down.clear();
		its_coeff_b_down.clear();
		its_coeff_c_down.clear();

		itsStrikes.clear();
		itsIsShortName.clear();

	}	


	ICM_Gaussian_LossDistrib_MF(const int& nbnames,
										   const ARM_Vector& pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector&  LossRates,
										   const int& NbSec,
										   const std::vector<int>& IndSec,
										   const int& CopulaType ,						//= qNO_COPULA
										   const qIntegratorChoice&	IntegrationMethod ,	//= qGAUSS_HERMITE
										   const int& IntegrationStep);					//= 20


	//LongShort
	/*********************************************************************/
	ICM_Gaussian_LossDistrib_MF(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  beta,
							const ARM_Vector&  LossRates,
							const int& NbSec,
							const std::vector<int>& IndSec,
							const std::vector<int>& SortedIndice,
							const int& CopulaType ,								//=qNO_COPULA
							const qIntegratorChoice&	IntegrationMethod,		//= qGAUSS_HERMITE
							const int& IntegrationStep);						//= 20
	/*********************************************************************/

	void Set(const int& nbnames,
			 const ARM_Vector& ,
			 const ARM_Vector& ,
			 const ARM_Vector&  LossRates,
			 const int& NbSec,
			 const std::vector<int>& IndSec,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);


	//LongShort
	/*********************************************************************/
	void Set(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  beta,
			 const ARM_Vector&  LossRates,
			 const int& NbSec,
			 const std::vector<int>& IndSec,
			 const std::vector<int>& SortedIndice,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);
	/*********************************************************************/

	void BitwiseCopy(const ARM_Object* src);

	void Copy(const ARM_Object* src);

	ARM_Object* Clone(void);


	double compute_expectedlosstranche(const double& tranche_up, 
									   const double& tranche_down, 
									   const double& lossunit,
									   //ICM_QMatrix<double>* ShiftMatrix = NULL,
									   const int& Tenor = -1,
									   const int& Issuer = -1);
	
	double compute_expectedlosstranche_LongShort(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  // const double& minloss
												  //ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  );
//----- Assessors
public:

	inline void SetCopulaType(const int& value){itsCopulaType = value;}
	inline void GetCopulaType(int& value){value = itsCopulaType;}

	inline void SetIntegrationMethod(const qIntegratorChoice& value){itsIntegrationMethod = value;}
	inline void GetIntegrationMethod(qIntegratorChoice& value){value = itsIntegrationMethod ;}
	inline qIntegratorChoice GetIntegrationMethod(){return itsIntegrationMethod ;}

	inline void SetIntegrationStep(const int& value){itsIntegrationStep = value;}
	inline void GetIntegrationStep(int& value){value = itsIntegrationStep ;}

	inline void SetFreedomDegree(const int& value){itsFreedomDegree = value;}
	inline void GetFreedomDegree(int& value){value = itsFreedomDegree;}

	inline void SetNbStrikes(const int& value){itsNbStrikes = value;}
	inline void GetNbStrikes(int& value){value = itsNbStrikes ;}

	void SetSmileParameters(const int& NbStrikes,
					   const vector<double>& Strikes,
					   ICM_QMatrix<double>* BetaMatrix);
	
	void UpdateSmileCorrelation(const double& tranche_down, 
								const double& tranche_up,
								const bool& HedgeFlag = false);
	



//----- Utilities

	//Estimation des probas conditionnelles
	double compute_cond_distribZeroLoss_Gaussian(double x,void* params);
	double compute_cond_distrib_Gaussian(double x, void* params);

	void InitProbCondPortSize0(const int& lmin, const int& ind_x);
/*	void InitAllProbCondPortSize0(const int& lmin, const int& ind_x, const int& lmax, const int& nbnames);*/
	void ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk);
	void ComputeProbCondDown(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk_down);
	void ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk);
	void ComputeProbCondDownShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk_down);
	double Compute_pk(const int& ind_x, const int& ind_sortedname);
	double Compute_pk_down(const int& ind_x, const int& ind_sortedname);

	// Estimation de la distribution de perte conditionnelle
	virtual void compute_distrib(const int& lup);
	virtual void compute_distrib(const int& lup, const int& ldown);

	virtual void compute_distrib_LongShort(const int& lup, const int& lmin);
	virtual void compute_distrib_LongShort(const int& lup, const int& ldown, const int& lmin);

	void integrate_SectorLoss(const int& lmax, const int& nb_fact_cond,const int& nb_name_secteur,const int& ind_x, const int& ind_Sect);
	void integrate_SectorLossDown(const int& lmax, const int& nb_fact_cond,const int& nb_name_secteur,const int& ind_x, const int& ind_Sect);

	void integrate_LossDistrib(const int& lmax, const int& nb_fact_cond);
	void integrate_LossDistribDown(const int& lmax, const int& nb_fact_cond);

	void initialize_SectorLoss(const int& lmax, const int& nb_fact_cond);
	void initialize_SectorLoss0(const int& lmin, const int& lmax, const int& nb_fact_cond);
	void initialize_SectorLossDown(const int& lmax, const int& nb_fact_cond);
	void initialize_SectorLoss0Down(const int& lmin, const int& lmax, const int& nb_fact_cond);


	//Fonctions intermédiaires 
	void compute_barrier();
	void precompute_coeffs(); // barrier computation
	void UpdateCorrelation(const double& beta_down, 
												  const double& beta_up,
												  const double& lambda_down, 
												  const double& lambda_up,
												  const bool& HedgeFlag = false);

	void UpdateCorrelation(std::vector<double>& V_beta_down,
												  std::vector<double>& V_beta_up,
												  std::vector<double>& V_lambda_down,
												  std::vector<double>& V_lambda_up, 
												  bool HedgeFlag = false);

	//Utilities
	void View(char* id, FILE* ficOut){};
	void ViewF(char* id, FILE* ficOut,const int& ind_x);
	void ViewS(char* id, FILE* ficOut,const int& ind_x);

	// INTEGRATOR
	void	SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, const int& TheStep = 0);

	//Assesseusr

	inline void	SetProbCond_Elt(const int& x, const int& y,const int& z, const double& value) {its_ProbCond->SetElt(x, y, z,value);}
	inline void	SetProbCond_Down_Elt(const int& x, const int& y, const int& z, const double& value) {its_ProbCond_Down->SetElt(x, y, z,value);}
	inline double	GetProbCond_Elt(const int& x,const int& y, const int& z) {return its_ProbCond->Elt(x, y, z);}
	inline double	GetProbCond_Down_Elt(const int& x,const int& y,const int& z) {return its_ProbCond_Down->Elt(x, y, z);}
	inline double	GetCoeff_a(const int& Num)	{return its_coeff_a[Num];}
	inline double	GetCoeff_b(const int& Num)	{return its_coeff_b[Num];}
	inline double	GetCoeff_c(const int& Num)	{return its_coeff_c[Num];}
	inline double	GetCoeff_a_down(const int& Num)	{return its_coeff_a_down[Num];}
	inline double	GetCoeff_b_down(const int& Num)	{return its_coeff_b_down[Num];}
	inline double	GetCoeff_c_down(const int& Num)	{return its_coeff_c_down[Num];}
	inline double	GetNbNameSector(const int& Num)	
	{
		int res=0;
		for(int i=0; i<itsIndicSectors.size(); i++)
			if (itsIndicSectors[i]==Num) res++;
		return res;
	}
	inline void SetIsShortName(const ARM_Vector&lossrates)
	{ 
		itsIsShortName.resize(lossrates.size()); 
		for(int k=0;k<lossrates.size();k++) 
		{
			if(lossrates[k] < 0) 
				itsIsShortName[k]=true;
			else
				itsIsShortName[k]=false;
		}
	}
};
#endif 


