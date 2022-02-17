/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		icm_Homogene_FastLossDistrib_Gaussian.H
	PROJECT:	MOD
	
	DESCRIPTION:	Loss Computation with Smile

  -----------------------------------------------------------------

 	CREATION:	September, 2004

	LAST MODIF:	October, 2004
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef _ICM_Gauss1FLossDistrib_LJ_H__
#define _ICM_Gauss1FLossDistrib_LJ_H__

#include <string>
#include <vector>
#include "ICMKernel\mod\icm_distribution.h"
#include "ICMKernel\util\icm_integrator.h"


// 17783 extern "C" void 
// 17783 compute_cond_distribZeroLoss_Gaussian_Smile_TSR(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down);

// 17783 extern "C" void 
// 17783 compute_cond_distrib_Gaussian_Smile_TSR(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down);

// 17783 extern "C" void 
// 17783 compute_cond_distrib_Gaussian_Integrator_TSR(void* Param, double x, double& res);

// 17783 extern "C" void 
// 17783 compute_cond_distribZeroLoss_Gaussian_Integrator_TSR(void* Param, double x, double& res);

// 17783 extern "C" void 
// 17783 compute_cond_distrib_Gaussian_Integrator_TSR(void* Param, double x, double& res);

class ICM_Gauss1FLossDistrib_LJ  : public ICM_Distribution  
{
//-----  Constructors/destructors
private:

	// Gauss Legendre Integration
	int		itsCopulaType;
	int		itsIntegrationStep;
	int		itsFreedomDegree;
	int		itsNbStrikes;

	double	its_taildistrib_Down;

	qIntegratorChoice	itsIntegrationMethod;

	ICM_QMatrix<double>*	itsBetasMatrix;

	// SMILE: parameters for Down Barrier
	ICM_QCubix<double>* its_ProbCond_Down; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QCubix<double>* its_ProbCond_Perturb_Down; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QMatrix<double>* its_lossdistrib_perturb_Down; // for assessor, m_dens[k] = P(N(T)=k)

	vector<double> its_beta_Down;		//beta strike dw (à its_t1) for maturity
	vector<double> its_beta_Down_t2;	//beta strike dw (à its_t2) for maturity	 

	vector<double> its_lossdistrib_Down; // for assessor, m_dens[k] = P(N(T)=k)
	vector<double> its_taildistrib_perturb_Down;

	// --------------------------------------
	// some pre-computed parameters, to be optimized, when dealing with a single correlation value,
	// instead of a vector
	vector<double>	its_coeff_a;
	vector<double>	its_coeff_b;
	vector<double>	its_coeff_a_down;
	vector<double>	its_coeff_b_down;
	vector<double>	its_coeff_a_perturb;
	vector<double>	its_coeff_a_down_perturb;
	vector<double>	itsStrikes;

	// --------------------------------------
	// Integrator
	ICM_Integrator	TheIntegrator;

	//Forward Collateral --------------------------------------------------------------
	double	its_collat_fwd_taildistrib_Down;

	// SMILE: parameters for Down Barrier
	ICM_QCubix<double>* its_collat_fwd_ProbCond_Down; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QCubix<double>* its_collat_fwd_ProbCond_Perturb_Down; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	ICM_QMatrix<double>* its_collat_fwd_lossdistrib_perturb_Down; // for assessor, m_dens[k] = P(N(T)=k)

	vector<double> its_collat_fwd_beta_Down;		//beta strike dw (à its_collat_t1) for start date
	vector<double>  its_collat_fwd_beta_Down_t2;	//beta strike dw (à its_collat_t2) for start date

	vector<double> its_collat_fwd_lossdistrib_Down; // for assessor, m_dens[k] = P(N(T)=k)
	vector<double> its_collat_fwd_taildistrib_perturb_Down;

	vector<double>	its_collat_fwd_coeff_a;
	vector<double>	its_collat_fwd_coeff_b;
	vector<double>	its_collat_fwd_coeff_a_down;
	vector<double>	its_collat_fwd_coeff_b_down;
	vector<double>	its_collat_fwd_coeff_a_perturb;
	vector<double>	its_collat_fwd_coeff_a_down_perturb;
	vector<double>	its_collat_fwdStrikes;

	// --------------------------------------------------------------------------------
	// Term structure review
	// --------------------------------------------------------------------------------

	vector<double>	its_coeff_b_perturb;
	vector<double>	its_coeff_b_down_perturb;

	vector<double>	its_coeff_a_t2;
	vector<double>	its_coeff_b_t2;
	vector<double>	its_coeff_a_down_t2;
	vector<double>	its_coeff_b_down_t2;

	vector<double>	its_coeff_a_perturb_t2;
	vector<double>	its_coeff_a_down_perturb_t2;
	vector<double>	its_coeff_b_perturb_t2;
	vector<double>	its_coeff_b_down_perturb_t2;

	// --------------------------------------------------------------------------------
	// step up subordination
	// --------------------------------------------------------------------------------

	ICM_QCubix<double>* its_ProbCond_stepup;	//stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l


public:

	void Init();

	ICM_Gauss1FLossDistrib_LJ()
	{
		Init();
	}	

	~ICM_Gauss1FLossDistrib_LJ()
	{
		if (itsBetasMatrix)
			delete itsBetasMatrix;
		itsBetasMatrix = NULL;

		if (its_ProbCond_Down)
			delete its_ProbCond_Down;
		its_ProbCond_Down = NULL;

		if (its_ProbCond_Perturb_Down)
			delete its_ProbCond_Perturb_Down;
		its_ProbCond_Perturb_Down = NULL;

		if (its_lossdistrib_perturb_Down)
			delete its_lossdistrib_perturb_Down;
		its_lossdistrib_perturb_Down = NULL;

		its_beta_Down.clear();		 
		its_lossdistrib_Down.clear(); 
		its_taildistrib_perturb_Down.clear();
		its_coeff_a.clear();
		its_coeff_b.clear();
		its_coeff_a_down.clear();
		its_coeff_b_down.clear();
		its_coeff_a_perturb.clear();
		its_coeff_a_down_perturb.clear();
		itsStrikes.clear();

		//Forward Collateral --------------------------------------------------------------
		if (its_collat_fwd_ProbCond_Down)
			delete its_collat_fwd_ProbCond_Down;
		its_collat_fwd_ProbCond_Down = NULL;

		if (its_collat_fwd_ProbCond_Perturb_Down)
			delete its_collat_fwd_ProbCond_Perturb_Down;
		its_collat_fwd_ProbCond_Perturb_Down = NULL;

		if (its_collat_fwd_lossdistrib_perturb_Down)
			delete its_collat_fwd_lossdistrib_perturb_Down;
		its_collat_fwd_lossdistrib_perturb_Down = NULL;

		its_collat_fwd_beta_Down.clear();
		its_collat_fwd_lossdistrib_Down.clear();
		its_collat_fwd_taildistrib_perturb_Down.clear();

		its_collat_fwd_coeff_a.clear();
		its_collat_fwd_coeff_b.clear();
		its_collat_fwd_coeff_a_down.clear();
		its_collat_fwd_coeff_b_down.clear();
		its_collat_fwd_coeff_a_perturb.clear();
		its_collat_fwd_coeff_a_down_perturb.clear();
		its_collat_fwdStrikes.clear();

		its_beta_Down_t2.clear();		 
		its_collat_fwd_beta_Down_t2.clear();
		its_coeff_b_perturb.clear();
		its_coeff_b_down_perturb.clear();
		its_coeff_a_t2.clear();
		its_coeff_b_t2.clear();
		its_coeff_a_down_t2.clear();
		its_coeff_b_down_t2.clear();
		its_coeff_a_perturb_t2.clear();
		its_coeff_a_down_perturb_t2.clear();
		its_coeff_b_perturb_t2.clear();
		its_coeff_b_down_perturb_t2.clear();	

		if (its_ProbCond_stepup)
			delete its_ProbCond_stepup;
		its_ProbCond_stepup = NULL;

	}	

	// ---------------------------------------------------------------------
	// Generic constructor
	// ---------------------------------------------------------------------
	ICM_Gauss1FLossDistrib_LJ(CtxtDistrib* c);
	void Set(CtxtDistrib* c);


	ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& beta,
							const ARM_Vector&LossRates,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);


	// ---------------------------------------------------------------------
	//LongShort
	// ---------------------------------------------------------------------
	ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& beta,
							const ARM_Vector&LossRates,
							const std::vector<int>& SortedIndice,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);

	// ---------------------------------------------------------------------
	//Collat Fwd
	// ---------------------------------------------------------------------
	ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& pdef_start,
							const ARM_Vector& beta,
							const ARM_Vector&   LossRates,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);
	
	// ---------------------------------------------------------------------
	void Set(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  beta,
			 const ARM_Vector&LossRates,
			 const int& discretizationstep = 20,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);

	void SetFwdCollat(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  beta,
			 const ARM_Vector& LossRates,
			 const int& discretizationstep = 20,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);

	// ---------------------------------------------------------------------
	//	LongShort
	// ---------------------------------------------------------------------
	void Set(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  beta,
			 const ARM_Vector&  LossRates,
			 const std::vector<int>& SortedIndice,
			 const int& discretizationstep = 20,
			 const int& CopulaType = qNO_COPULA,
			 const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
			 const int& IntegrationStep = 20);

	// ---------------------------------------------------------------------
	//Collat Fwd
	// ---------------------------------------------------------------------
	void Set(const int& nbnames,
			 const ARM_Vector&  pdefault,
			 const ARM_Vector&  pdef_start,
			 const ARM_Vector&  beta,
			 const ARM_Vector& LossRates,
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
									   // ICM_QMatrix<double>* ShiftMatrix,
									   // const int& Tenor,
									   // const int& Issuer,
									   vector<double>& losses,
									   vector<double>& ProbasLossesUp,
									   vector<double>& ProbasLossesDown);

	//-------------------------------------------------------------------
	// Generic EL computation
	//-------------------------------------------------------------------
	double ComputeEL(CtxtDistrib* ctxt);

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
	
	void UpdateCorrelation(const double& beta_down, 
							const double& beta_up, 
							const bool& HedgeFlag = false);

	void UpdateCorrelation(const ARM_Vector& V_beta_down, 
						   const ARM_Vector& V_beta_up, 
						   bool HedgeFlag = false);

	void Collat_fwd_UpdateCorrelation(std::vector<double>& V_beta_down, 
									std::vector<double>& V_beta_up, 
									bool HedgeFlag = false);


//----- Utilities

	double compute_cond_distribZeroLoss_Gaussian(double x,void* params);
	double compute_cond_distribZeroLoss_Student(const double& x);
	double compute_cond_distrib_Gaussian(double x, void* params);
	double compute_cond_distrib_Student(const double& x);

	void compute_cond_distribZeroLoss_Student_Smile(const double& x, 
													double& ProbCond_Down,
													double& ProbCond_Up) {};

	void compute_cond_distrib_Student_Smile(const double& x, 
											const double& ProbCond_Down, 
											const double& ProbCond_Up) {};

	// No Smile
	virtual void compute_distrib(const int& lup);
	virtual void compute_distrib_Gauss_Legendre(const int& lup);
	//virtual void compute_distrib_Hermite(const int& lup);

	// Smile
	virtual void compute_distrib(const int& lup, const int& ldown);

	void compute_barrier();

	void precompute_coeffs();				// barrier computation
	void precompute_coeffs_perturb();		// barrier computation perturb


	// ---------------------------------------------------------------------
	//Long Short CDO
	// ---------------------------------------------------------------------
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
	virtual void compute_distrib_LongShort(const int& lup, const int& ldown, const int& lmin);

	//Utilities
	void InitProbCondPortSize0(const int& lmin, const int& ind_x);
	void ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk);
	void ComputeProbCondDown(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk_down);
	void ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk);
	void ComputeProbCondDownShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk_down);
	double Compute_pk(const int& ind_x, const int& ind_sortedname);
	double Compute_pk_down(const int& ind_x, const int& ind_sortedname);

	void View(char* id, FILE* ficOut);
	
	// ---------------------------------------------------------------------
	//Collat Fwd Start
	// ---------------------------------------------------------------------
	double compute_expectedlosstranche_CollatFwd(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  //ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  );
	//Distribution
	virtual void compute_distrib_CollatFwd(const int& lup);
	virtual void compute_distrib_CollatFwd(const int& lup, const int& ldown);
	
	// ---------------------------------------------------------------------
	// HEDGES
	// -----------------------------------------------------------------------------------
// 17783 	void compute_expectedlosstranche_fast_spread_hedge(const double& tranche_up, 
// 17783 											   const double& tranche_down, 
// 17783 											   const double& lossunit,
// 17783 											   ICM_QMatrix<double>* ShiftMatrix);

// 17783  	void compute_expectedlosstranche_perturb(const double& tranche_up, 
// 17783 											 const double& tranche_down, 
// 17783 											 const double& lossunit,
// 17783 											ICM_QMatrix<double>* ShiftMatrix);


	// 17783 double compute_cond_distribZeroLoss_shift_k_Gaussian(double x);
	// 17783 double compute_cond_distrib_shift_k_Gaussian(double x);


	// Smile
	// 17783 virtual void compute_distrib_perturb(const int& lup, const int& ldown);


	// No Smile
	// 17783 virtual void compute_distrib_perturb(const int& lup);
	// 17783 virtual void compute_distrib_perturb_Gauss_Legendre(const int& lup);
	// 17783 virtual void compute_distrib_perturb_Hermite(const int& lup);
	
	// -----------------------------------------------------------------------------------


	// INTEGRATOR
	void	SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, const int& TheStep = 0);


//-------------------

	// no tests

	inline void	SetProbCond_Elt(const int& x, const int& y,const int& z, const double& value) {its_ProbCond->SetElt(x, y, z,value);}
	inline void	SetProbCond_Down_Elt(const int& x, const int& y, const int& z, const double& value) {its_ProbCond_Down->SetElt(x, y, z,value);}
	inline double	GetProbCond_Elt(const int& x,const int& y, const int& z) {return its_ProbCond->Elt(x, y, z);}
	inline double	GetProbCond_Down_Elt(const int& x,const int& y,const int& z) {return its_ProbCond_Down->Elt(x, y, z);}
	// 17783 inline void	SetProbCond_Perturb_Elt(const int& x,const int& y,const int& z,const double& value) {its_ProbCond_Perturb->SetElt(x, y, z,value);}
	inline void	SetProbCond_Perturb_Down_Elt(const int& x,const int& y, const int& z,const double& value) {its_ProbCond_Perturb_Down->SetElt(x, y, z,value);}
	// 17783 inline double	GetProbCond_Perturb_Elt(const int& x,const int& y,const int& z) {return its_ProbCond_Perturb->Elt(x, y, z);}
	inline double	GetProbCond_Perturb_Down_Elt(const int& x, const int& y,const int& z) {return its_ProbCond_Perturb_Down->Elt(x, y, z);}
	inline double	GetCoeff_a(const int& Num)	{return its_coeff_a[Num];}
	inline double	GetCoeff_b(const int& Num)	{return its_coeff_b[Num];}
	inline double	GetCoeff_a_down(const int& Num)	{return its_coeff_a_down[Num];}
	inline double	GetCoeff_b_down(const int& Num)	{return its_coeff_b_down[Num];}
	inline double	GetCoeff_Perturb_a(const int& Num)	{return its_coeff_a_perturb[Num];}
	inline double	GetCoeff_Perturb_a_down(const int& Num)	{return its_coeff_a_down_perturb[Num];}

	//Forward Collateral
	inline double	Get_collat_fwdCoeff_a(const int& Num)	{return its_collat_fwd_coeff_a[Num];}
	inline double	Get_collat_fwdCoeff_b(const int& Num)	{return its_collat_fwd_coeff_b[Num];}
	inline double	Get_collat_fwdCoeff_a_down(const int& Num)	{return its_collat_fwd_coeff_a_down[Num];}
	inline double	Get_collat_fwdCoeff_b_down(const int& Num)	{return its_collat_fwd_coeff_b_down[Num];}
	inline double	Get_collat_fwdCoeff_Perturb_a(const int& Num)	{return its_collat_fwd_coeff_a_perturb[Num];}
	inline double	Get_collat_fwdCoeff_Perturb_a_down(const int& Num)	{return its_collat_fwd_coeff_a_down_perturb[Num];}

	//---------------------------------------------------------------------------------
	//cas term structure review
	//---------------------------------------------------------------------------------

	ICM_Gauss1FLossDistrib_LJ(const double& t1_corr,
							const double& t2_corr,
							const int& nbnames,
							const ARM_Vector& pdefault_maturity_up_t1,
							const ARM_Vector&  pdefault_maturity_up_t2,
							const ARM_Vector&  pdefault_maturity_dw_t1,
							const ARM_Vector&  pdefault_maturity_dw_t2,
							const ARM_Vector&  beta_t1,
							const ARM_Vector&  beta_t2,
							const ARM_Vector & LossRates,
							const int& discretizationstep = 20,
							const int& CopulaType = qNO_COPULA,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE,
							const int& IntegrationStep = 20);

	void Set(const double& t1_corr,
							const double& t2_corr,
							const int& nbnames,
							const ARM_Vector&  pdefault_maturity_up_t1,
							const ARM_Vector&  pdefault_maturity_up_t2,
							const ARM_Vector&  pdefault_maturity_dw_t1,
							const ARM_Vector&  pdefault_maturity_dw_t2,
							const ARM_Vector&  beta_t1,
							const ARM_Vector&  beta_t2,
							const ARM_Vector & LossRates,
							const int& discretizationstep = 20 ,
							const int& CopulaType = qNO_COPULA ,
							const qIntegratorChoice&	IntegrationMethod = qGAUSS_HERMITE ,
							const int& IntegrationStep = 20);

	void UpdateCorrelation_TSR(const ARM_Vector & V_beta_down_t1, 
								const ARM_Vector & V_beta_down_t2, 
								const ARM_Vector & V_beta_up_t1, 
								const ARM_Vector & V_beta_up_t2);

	void precompute_coeffs_TSR(); // barrier computation

	//---------------------------------------------------------------------------------
	//cas term structure review : collat forward
	//---------------------------------------------------------------------------------

	ICM_Gauss1FLossDistrib_LJ(const double& t1_corr,
													 const double& t2_corr,
													 const double& start_t1_corr,
													 const double& start_t2_corr,
													 const int& nbnames,
													 const ARM_Vector&  pdef,
													 const ARM_Vector&  pdef_t2,
													 const ARM_Vector&  pdef_start,
													 const ARM_Vector&  pdef_start_t2,
													 const ARM_Vector&  beta_t1,
													 const ARM_Vector&  beta_t2,
													 const ARM_Vector&  start_beta_t1,
													 const ARM_Vector&  start_beta_t2,
													 const ARM_Vector & LossRates,
													 const int& discretizationstep  = 20,
													 const int& CopulaType = qNO_COPULA,
													 const qIntegratorChoice&	IntegrationMethod  = qGAUSS_HERMITE,
													 const int& IntegrationStep= 20);


	void Set(const double& t1_corr,
													 const double& t2_corr,
													 const double& start_t1_corr,
													 const double& start_t2_corr,
													 const int& nbnames,
													 const ARM_Vector&  pdef,
													 const ARM_Vector&  pdef_t2,
													 const ARM_Vector&  pdef_start,
													 const ARM_Vector&  pdef_start_t2,
													 const ARM_Vector&  beta_t1,
													 const ARM_Vector&  beta_t2,
													 const ARM_Vector&  start_beta_t1,
													 const ARM_Vector&  start_beta_t2,
													 const ARM_Vector & LossRates,
													 const int& discretizationstep  = 20,
													 const int& CopulaType = qNO_COPULA,
													 const qIntegratorChoice&	IntegrationMethod  = qGAUSS_HERMITE,
													 const int& IntegrationStep= 20);


	void SetFwdCollat_TSR(const int& nbnames,
													 const ARM_Vector&  pdef,
													 const ARM_Vector&  pdef_t2,
													 const ARM_Vector&  beta_t1,
													 const ARM_Vector&  beta_t2,
													 const ARM_Vector&  start_beta_t1,
													 const ARM_Vector&  start_beta_t2,
													 const ARM_Vector & LossRates,
													 const int& discretizationstep  = 20,
													 const int& CopulaType = qNO_COPULA,
													 const qIntegratorChoice&	IntegrationMethod  = qGAUSS_HERMITE,
													 const int& IntegrationStep= 20);

	void Collat_fwd_UpdateCorrelation_TSR(const ARM_Vector & V_beta_down,
													const ARM_Vector  & V_beta_down_t2,
													const ARM_Vector & V_beta_up,
													const ARM_Vector & V_beta_up_t2,
													const ARM_Vector & V_collat_beta_down,
													const ARM_Vector & V_collat_beta_down_t2,
													const ARM_Vector & V_collat_beta_up,
													const ARM_Vector & V_collat_beta_up_t2);

	inline double	GetCoeff_a_t2(const int& Num)	{return its_coeff_a_t2[Num];}
	inline double	GetCoeff_b_t2(const int& Num)	{return its_coeff_b_t2[Num];}
	inline double	GetCoeff_a_down_t2(const int& Num)	{return its_coeff_a_down_t2[Num];}
	inline double	GetCoeff_b_down_t2(const int& Num)	{return its_coeff_b_down_t2[Num];}
	inline double	GetCoeff_Perturb_a_t2(const int& Num)	{return its_coeff_a_perturb_t2[Num];}
	inline double	GetCoeff_Perturb_a_down_t2(const int& Num)	{return its_coeff_a_down_perturb_t2[Num];}

	inline double	GetCoeff_Perturb_b_t2(const int& Num)	{return its_coeff_b_perturb_t2[Num];}
	inline double	GetCoeff_Perturb_b_down_t2(const int& Num)	{return its_coeff_b_down_perturb_t2[Num];}
	inline double	GetCoeff_Perturb_b(const int& Num)	{return its_coeff_b_perturb[Num];}
	inline double	GetCoeff_Perturb_b_down(const int& Num)	{return its_coeff_b_down_perturb[Num];}

	void precompute_coeffs_perturb_TSR();	// barrier computation perturb TSR

	// ---------------------------------------------------------------------
	// Step up computation
	// ---------------------------------------------------------------------
	double ComputeStepUp(CtxtDistrib* ctxt);
	double ZeroLoss_TSR_stepup(double x,void* params);

	// ---------------------------------------------------------------------
	// Notionels variables
	// ---------------------------------------------------------------------
	void Set_VN(CtxtDistrib* c);
	void compute_distrib_VN(const int& lup);
	double compute_cond_distribZeroLoss_Gaussian_VN(double x, void* params);
	double compute_cond_distrib_Gaussian_VN(double x,void* params);
	double compute_cond_distribZeroLoss_Student_VN(const double& x);
	double compute_cond_distrib_Student_VN(const double& x);
	void compute_distrib_VN(const int& lup,const int& ldown);
	void compute_distrib_Gauss_Legendre_VN(const int& lup);
	double compute_expectedlosstranche_VN(CtxtDistrib* c);
	void SetIntegratorType_VN(const qIntegratorChoice&	TheIntegratorType, 
														const int&	TheStep);

	// -------------------------------------------------------------------------
	// Approximation
	// -------------------------------------------------------------------------
	void Set_APPROX(CtxtDistrib* c);
	double compute_cond_distrib_Gaussian_APPROX(double x,void* params);
	double compute_distrib_APPROX(const int& lup,const int& ldown);
	double compute_expectedlosstranche_APPROX(CtxtDistrib* c);
	void SetIntegratorType_APPROX(const qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep);

	// ---------------------------------------------------------------------
	// Notionels variables
	// ---------------------------------------------------------------------
	void Set_STD(CtxtDistrib* c);
	void compute_distrib_STD(const int& lup);
	double compute_cond_distribZeroLoss_Gaussian_STD(double x, void* params);
	double compute_cond_distrib_Gaussian_STD(double x,void* params);
	double compute_cond_distribZeroLoss_Student_STD(const double& x);
	double compute_cond_distrib_Student_STD(const double& x);
	void compute_distrib_STD(const int& lup,const int& ldown);
	void compute_distrib_Gauss_Legendre_STD(const int& lup);
	double compute_expectedlosstranche_STD(CtxtDistrib* c);
	void SetIntegratorType_STD(const qIntegratorChoice&	TheIntegratorType, const int&	TheStep);

};	


#endif 
