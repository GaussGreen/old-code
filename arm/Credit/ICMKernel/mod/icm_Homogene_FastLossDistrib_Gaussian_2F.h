/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		icm_Homogene_FastLossDistrib_Gaussian_2F.H
	PROJECT:	MOD
	
	DESCRIPTION:	Loss Computation with Smile

  -----------------------------------------------------------------

 	CREATION:	February, 2006

	LAST MODIF:	February, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef __ICM_Gaussian_LossDistrib_2F_H__
#define __ICM_Gaussian_LossDistrib_2F_H__

#include <string>
#include <vector>

#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\mod\icm_homogene_fastlossdistrib.h"

#include "ICMKernel\util\icm_integrator_dim_2.h"



class ICM_Gaussian_LossDistrib_2F  : public ICM_Gauss1FLossDistrib  
{
//-----  Constructors/destructors
protected:

	// Sectorial Numerical Integration 
	qCopula_TYPE		its_CopulaType;		// practical value: GAUSSIAN

	int		its_IntegrationStep_Common_Factor;
	int		its_IntegrationStep_Sector_Factor;

	qIntegratorChoice	itsIntegrationMethod_Common_Factor;
	qIntegratorChoice	itsIntegrationMethod_Sector_Factor;

	qTWO_FACTORS_CORRELATION_TYPE	its_CorrelationType;

	// --------------------------------------
	// some pre-computed parameters, to be optimized, when dealing with a single correlation value,
	// instead of a vector
	vector<double>	its_coeff_a;
	vector<double>	its_coeff_b;
	vector<double>	its_coeff_c;
	vector<double>	its_coeff_den;

	// later on for perturbation...

	// --------------------------------------
	// Integrator
	ICM_Integrator_Dim_2	TheIntegrator;

	// --------------------------------------
	// SECTORIAL APPROACH
	vector<int>		its_sector_membership;

	double			its_Single_intra_sector_correlation;
	double			its_Single_inter_sector_correlation;

	vector<double>	its_Betas;
	vector<double>	its_Lambdas;

protected:

	// INTERNAL
	int				its_Nb_Sectors;

	double			its_SQRT_OneMinus_intra_sector_correlation;
	double			its_SQRT_inter_sector_correlation;
	double			its_SQRT_inter_minus_intra_sector_correlation;

	DoubleVector	its_Barriers;
	
	// INTEGRATION
	
	// Current state with factors (common and sector)
	// ids...
	int			its_factor_state;
	int			its_sector_state;
	
	int			its_sector_id;		// SECTOR NUMBER ID

	// values
	double		its_factor_value;
	double		its_sector_value;


	int				its_index_loss;		// current state for loss

	DoubleVector	its_Current_DefProb;	// for index_z

	// output of 'LossProbabilityDistribution' method
	DoubleVector	its_lossdistrib;		// for assessor, m_dens[k] = P(N(T)=k)
	double			its_taildistrib;

	bool		its_first_passing_flag;

public:

	void Init();

	ICM_Gaussian_LossDistrib_2F()
	{
		Init();
	}	

	~ICM_Gaussian_LossDistrib_2F()
	{
		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_c.clear();
		its_coeff_den.clear();

		its_sector_membership.clear();

		its_Betas.clear();
		its_Lambdas.clear();

		its_Barriers.clear();
		its_Current_DefProb.clear();

		its_lossdistrib.clear();
	}	


	ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& betas_sector,
							const ARM_Vector& lambdas_sector,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qGAUSSIAN,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& betas_sector,
							double	intra_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qGAUSSIAN,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	ICM_Gaussian_LossDistrib_2F(const int& nbnames,
							const ARM_Vector& pdefault,
							double	intra_sector_correl,
							double	inter_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qGAUSSIAN,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	void Set(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& betas_sector,
							const ARM_Vector& lambdas_sector,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	void Set(const int& nbnames,
							const ARM_Vector&  pdefault,
							const ARM_Vector&  betas_sector,
							double	intra_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	void Set(const int& nbnames,
							const ARM_Vector& pdefault,
							double	intra_sector_correl,
							double	inter_sector_correl,
							int*	sector_membership,
							int		nb_sectors,
							const ARM_Vector&  LossRates,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20,
							const qIntegratorChoice&	IntegrationMethod_2 = qGAUSS_HERMITE,
							const int& IntegrationStep_2 = 20
							);

	void BitwiseCopy(const ARM_Object* src);

	void Copy(const ARM_Object* src);

	ARM_Object* Clone(void);


	double compute_expectedlosstranche(double& tranche_down, 
									   double& tranche_up, 
									   const double& lossunit,
									   // ICM_QMatrix<double>* ShiftMatrix = NULL,
									   const int& Tenor = -1,
									   const int& Issuer = -1);
//----- Assessors
public:

	inline void SetCopulaType(const qCopula_TYPE& value){its_CopulaType = value;}
	inline void GetCopulaType(qCopula_TYPE& value){value = its_CopulaType;}

	inline void SetCorrelationType(const qTWO_FACTORS_CORRELATION_TYPE& value){its_CorrelationType = value;}
	inline void GetCorrelationType(qTWO_FACTORS_CORRELATION_TYPE& value){value = its_CorrelationType;}
	
	inline void SetIntegrationMethod_Common_Factor(const qIntegratorChoice& value){itsIntegrationMethod_Common_Factor = value;}
	inline void GetIntegrationMethod_Common_Factor(qIntegratorChoice& value){value = itsIntegrationMethod_Common_Factor;}
	inline qIntegratorChoice GetIntegrationMethod_Common_Factor(){return itsIntegrationMethod_Common_Factor;}

	inline void SetIntegrationMethod_Sector_Factor(const qIntegratorChoice& value){itsIntegrationMethod_Sector_Factor = value;}
	inline void GetIntegrationMethod_Sector_Factor(qIntegratorChoice& value){value = itsIntegrationMethod_Sector_Factor;}
	inline qIntegratorChoice GetIntegrationMethod_Sector_Factor(){return itsIntegrationMethod_Sector_Factor;}

	inline void SetIntegrationStep_Common_Factor(const int& value){its_IntegrationStep_Common_Factor = value;}
	inline void GetIntegrationStep_Sector_Factor(int& value){value = its_IntegrationStep_Sector_Factor;}

	inline void SetNbSectors(const int& value){its_Nb_Sectors = value;}
	inline void GetNbSectors(int& value){value = its_Nb_Sectors;}



	//	INTEGRATION PROCESS
	//	CURRENT STATE 
	inline void SetFactorState(const int& value) {its_factor_state = value;}
	inline void GetFactorState(int& value) {value = its_factor_state;}

	inline void SetSectorState(const int& value) {its_sector_state = value;}
	inline void GetSectorState(int& value) {value = its_sector_state;}

	inline void SetFactorValue(const double& value) {its_factor_value = value;}
	inline void GetFactorValue(double& value) {value = its_factor_value;}

	inline void SetSectorValue(const double& value) {its_sector_value = value;}
	inline void GetSectorValue(double& value) {value = its_sector_value;}

	inline void SetSectorId(const int& value) {its_sector_id = value;}
	inline void GetSectorId(int& value) {value = its_sector_id;}

//----- Utilities

public:
	
	void	ComputeConditionalDefaultProbabilities(double& value_factor, double& value_sector);
	double	ComputeExpectedLossTranche(const double& LossMin, const double& LossMax);

	void	Compute_Conditional_Probabilities_Coeffs();
	void	Update_Conditional_Probabilities_Coeffs();

private:

	void	LossProbabilityDistribution(const int& lup);

	double	Compute_Conditional_Distribution_Loss_Zero();
	double	Compute_Conditional_Distribution_Losses();

/*
	double compute_cond_distribZeroLoss_Gaussian(double x,void* params);
	double compute_cond_distribZeroLoss_Student(const double& x);
	double compute_cond_distrib_Gaussian(double x, void* params);
	double compute_cond_distrib_Student(const double& x);

	// No Smile
	virtual void compute_distrib(const int& lup);
	virtual void compute_distrib_Gauss_Legendre(const int& lup);

	// Smile
	virtual void compute_distrib(const int& lup, const int& ldown);

	void precompute_coeffs(); // barrier computation
	void precompute_coeffs_perturb(); // barrier computation


	/*********************************************************************/
	//Utilities
/*	void InitProbCondPortSize0(const int& lmin, const int& ind_x);
	void ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname);
	void ComputeProbCondDown(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname);
	void ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname);
	void ComputeProbCondDownShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname);

	/*********************************************************************/

	
	// -----------------------------------------------------------------------------------
	//	UTILITIES
	// -----------------------------------------------------------------------------------
	// Generates all Betas and Lambdas values

	void	Compute_Barriers();

	void	SetIntegratorType(qIntegratorChoice&	TheIntegratorType, const int&	TheStep);

	// -----------------------------------------------------------------------------------
	// HEDGES
/*	void compute_expectedlosstranche_fast_spread_hedge(const double& tranche_up, 
											   const double& tranche_down, 
											   const double& lossunit,
											   ICM_QMatrix<double>* ShiftMatrix);

 	void compute_expectedlosstranche_perturb(const double& tranche_up, 
											 const double& tranche_down, 
											 const double& lossunit,
											ICM_QMatrix<double>* ShiftMatrix);


	double compute_cond_distribZeroLoss_shift_k_Gaussian(double x);
	double compute_cond_distrib_shift_k_Gaussian(double x);


	// Smile
	virtual void compute_distrib_perturb(const int& lup, const int& ldown);


	// No Smile
	virtual void compute_distrib_perturb(const int& lup);
	virtual void compute_distrib_perturb_Gauss_Legendre(const int& lup);
	virtual void compute_distrib_perturb_Hermite(const int& lup);
*/	
	// -----------------------------------------------------------------------------------
	// ACCESS TO ITEMS

	inline	double	GetCoeff_a(const int& Num)	{return its_coeff_a[Num];}
	inline	double	GetCoeff_b(const int& Num)	{return its_coeff_b[Num];}
	inline	double	GetCoeff_c(const int& Num)	{return its_coeff_c[Num];}
	inline	double	GetCoeff_den(const int& Num)	{return its_coeff_den[Num];}
	inline	int		GetSectorMemberShip(const int& Num) {return its_sector_membership[Num];}

protected:

	inline	bool	DoesThisNameBelongToSector_s(int& Num, int& s) {return (its_sector_membership[Num] == s);}

	// -----------------------------------------------------------------------------------

public:

	void	Set_First_Passing() {its_first_passing_flag = true;}
	void	Set_Not_First_Passing() {its_first_passing_flag = false;}
	bool	Is_First_Passing() {return its_first_passing_flag;}
	

public:

	void View(char* id, FILE* ficOut);

};
#endif 
