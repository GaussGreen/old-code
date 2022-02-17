#if !defined(AFX_ICM_Gauss1FLossDistrib_H__)
#define AFX_ICM_Gauss1FLossDistrib_H__


#include <string>
#include <vector>
#include "ICMKernel\mod\icm_distribution.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\util\icm_HermiteIntegration.h"


/*********************************************************************************/
/*! \class  ICM_Gauss1FLossDistrib icm_Homogene_FastLossDistrib.h "icm_Homogene_FastLossDistrib.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004
 *	\brief  Cette classe définit une distribution Gaussienne à 1 Facteur */
/***********************************************************************************/

class ICM_Gauss1FLossDistrib  : public ICM_Distribution  
{

//-----  Constructors/destructors
private:

ICM_QMatrix<double>* its_ProbCond_indep;

double its_proportion_full_correl;
double its_proportion_independant;

//Full Homogeneous
int its_FH_i;
double its_FH_Barrier;
double its_FH_Beta;

public:

static  ICM_QMatrix<double> itsCombinations; //Pascal Combinations for full homogeneous case

	ICM_Gauss1FLossDistrib(const int& nbnames)
	{
		Init();
		Set(nbnames);
	}

	void Set(const int& nbnames)
	{
		ICM_Distribution::Set(nbnames);

		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,20,1,0.);
	}

	void Init();

	ICM_Gauss1FLossDistrib() {Init();}	

	ICM_Gauss1FLossDistrib(const int& nbnames,
					     const ARM_Vector& pdefault,
						 const ARM_Vector& beta,
						 const ARM_Vector& LossRates);

	void Set(const int& nbnames,
		     const ARM_Vector& pdef,
			 const ARM_Vector& beta,
			 const ARM_Vector& LossRates);

	void BitwiseCopy(const ARM_Object* src);

	void Copy(const ARM_Object* src);


	ARM_Object* Clone(void);

	virtual ~ICM_Gauss1FLossDistrib();

//----- Assessors
public:

	inline void setPercentIndep(const double& value){ its_proportion_independant = value;}
	inline void setPercentFullCorrel(const double& value){ its_proportion_full_correl = value;}


	virtual double compute_cond_distribZeroLoss(const double& x);
	// 17783 virtual double compute_cond_distribZeroLoss_shift_k(const double& x);

	double compute_cond_distrib(const double& x);
	// 17783 double compute_cond_distrib_shift_k(const double& x);

	virtual void compute_distrib(const int& lup);
	// 17783 void compute_distrib_perturb(const int& lup);
	void compute_distrib_independant_case(); // dirty seulement dans le cas full homogène

	virtual double compute_expectedlosstranche(const double& tranche_up, 
											   const double& tranche_down, 
											   const double& lossunit,
											   //ICM_QMatrix<double>* ShiftMatrix,
											   // const int& Tenor,
											   // const int& Issuer,
											   vector<double>& losses);

// 17783 	virtual void compute_expectedlosstranche_fast_spread_hedge(const double& tranche_up, 
// 17783 											   const double& tranche_down, 
// 17783 											   const double& lossunit,
// 17783 											   ICM_QMatrix<double>* ShiftMatrix);

// 17783  	virtual void compute_expectedlosstranche_perturb(const double& tranche_up, 
// 17783 											  const double& tranche_down, 
// 17783 											  const double& lossunit,
// 17783 											  ICM_QMatrix<double>* ShiftMatrix);

	double ZeroLossFullHomog(const double& x);

	double ComputeEL_FullHomog(const double& LossUnit,
							   const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   const int& IntStep,
							   vector<double>& losses);


};

//----------------------------------------------------------------
// Homogeneous case only
//----------------------------------------------------------------
class ICM_LightDistrib : public ICM_Distribution  
{

//-----  Constructors/destructors
private:

bool   its_FullProba;

//Full Homogeneous
int its_nbnames;
int its_FH_i;
double its_FH_Barrier;
double its_FH_Beta;
double its_sqrt2;
double its_normPi;

//Full Homogeneous Collat Fwd
double its_FH_collat_Barrier;
double its_FH_collat_Beta;

//term structure
double its_FH_Barrier_T;
double its_FH_Beta_T;
double its_T1;
double its_T2;

//Full Homogeneous Collat Fwd
double its_FH_collat_Barrier_T;
double its_FH_collat_Beta_T;
double its_collat_T1;
double its_collat_T2;

//new barriers
double its_FH_Barrier_down;
double its_FH_Barrier_down_T;
double its_FH_collat_Barrier_down;
double its_FH_collat_Barrier_down_T;
bool   its_IsUp;

//step up
ICM_QCubix<double>* its_ProbCond_stepup;	//stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
int		its_ind_x_stepup;						//common factor

public:

static  ICM_QMatrix<double> itsCombinations; //Pascal Combinations for full homogeneous case

	ICM_LightDistrib() {Init();}
	
	ICM_LightDistrib(const int& nbnames)
	{
		Init();
		its_nbnames= nbnames;
	}

	void Init()
	{
		its_nbnames=0;
		its_FH_i=-1;
		its_FH_Barrier=0.;
		its_FH_Beta=0.;
		its_sqrt2 = sqrt(2.);
		its_normPi = 1./sqrt(PI);
		its_FH_Barrier_T=0.;
		its_FH_Beta_T=0.;
		its_T1=0.;
		its_T2=0.;
		its_FH_collat_Barrier=0.;
		its_FH_collat_Beta=0.;
		its_FH_collat_Barrier_T=0.;
		its_FH_collat_Beta_T=0.;
		its_collat_T1=0.;
		its_collat_T2=0.;
		//new barriers
		its_FH_Barrier_down=0.;
		its_FH_Barrier_down_T=0.;
		its_FH_collat_Barrier_down=0.;
		its_FH_collat_Barrier_down_T=0.;
		its_IsUp=true;
		its_ProbCond_stepup=NULL;
		its_ind_x_stepup=-1;
		its_FullProba=true;
	}	

public:

	double ZeroLossFullHomog(const double& x);
	double ZeroLossFullHomog_TSR(const double& x);

	//-------------------------------------------------------------------
	// Generic EL computation
	//-------------------------------------------------------------------
	double ComputeEL_FullHomog(CtxtDistrib* ctxt);

	//-------------------------------------------------------------------
	// Full Homogeneous case with integration
	//-------------------------------------------------------------------
	 double ComputeEL_FullHomog(const double& LossUnit,
							   const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   const int& IntStep,
							   vector<double>& losses);

	//-------------------------------------------------------------------
	// Full Homogeneous case with integration TSR
	//-------------------------------------------------------------------
	double ComputeEL_FullHomog_TSR(const double& LossUnit,
										const double& tranche_down,
										const double& tranche_up,
										const double& beta_down_T1,
										const double& beta_down_T2,
										const double& beta_up_T1,
										const double& beta_up_T2,
										const double& Pdef_T1,
										const double& Pdef_T2,
										const double& Pdef_down_T1,
										const double& Pdef_down_T2,
										const double& T1,
										const double& T2,
										const int& IntStep,
										vector<double>& losses);

	//-------------------------------------------------------------------
	// LHP
	//-------------------------------------------------------------------

	double ComputeEL_LHP(const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   double recovery,
							   vector<double>& losses);

	//-------------------------------------------------------------------
	// Full Homogeneous case with integration
	//-------------------------------------------------------------------
	double ComputeEL_FullHomog_collat(const double& LossUnit,
							   const double& tranche_down,
							   const double& tranche_up,
							   const double& beta_down,
							   const double& beta_up,
							   const double& Pdef,
							   const double& collat_beta_down,
							   const double& collat_beta_up,
							   const double& collat_Pdef,
							   const int& IntStep,
							   vector<double>& losses);

	//-------------------------------------------------------------------
	// Full Homogeneous case with integration TSR
	//-------------------------------------------------------------------
	double ComputeEL_FullHomog_collat_TSR(const double& LossUnit,
										const double& tranche_down,
										const double& tranche_up,
										const double& beta_down_T1,
										const double& beta_down_T2,
										const double& beta_up_T1,
										const double& beta_up_T2,
										const double& Pdef_T1,
										const double& Pdef_T2,
										const double& Pdef_down_T1,
										const double& Pdef_down_T2,
										const double& T1,
										const double& T2,
										const double& collat_beta_down_T1,
										const double& collat_beta_down_T2,
										const double& collat_beta_up_T1,
										const double& collat_beta_up_T2,
										const double& collat_Pdef_T1,
										const double& collat_Pdef_T2,
										const double& collat_Pdef_down_T1,
										const double& collat_Pdef_down_T2,
										const double& collat_T1,
										const double& collat_T2,
										const int& IntStep,
										vector<double>& losses);

	// ---------------------------------------------------------------------
	// Step up computation
	// ---------------------------------------------------------------------
	double ComputeStepUp(CtxtDistrib* ctxt);
	double ZeroLossFullHomog_TSR_stepup(const double& x);
};

#endif 
