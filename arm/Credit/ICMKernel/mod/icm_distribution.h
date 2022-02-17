
#if !defined(AFX_DISTRIBUTION_H__)
#define AFX_DISTRIBUTION_H__

#include <string>
#include <vector>
#include "ARMKernel\mod\model.h"
#include "ICMKernel\util\icm_qmatrix.h"
#include "ICMKernel\util\icm_pgcd.h"
#include "ICMKernel\util\icm_utils.h"

/*********************************************************************************/
/*! \class  ICM_Distribution icm_distribution.h "icm_distribution.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004
 *	\brief  Cette classe définit la classe Root des distributions */
/***********************************************************************************/

#ifndef MIN
#define   MIN(X,Y) ((X)<(Y)?(X):(Y))
#endif


#ifndef PI
#define PI		3.14159265359
#endif

#ifndef RHO_DIST
//#define   RHO_DIST(X,Y) sqrt(X/Y)
#define   RHO_DIST(X,Y) 1
#endif

class ICM_Pricer;
// void CptStrikesLosses(double k1,double k2,ICM_Pricer* p,int nblosses, double& Outk1,double& Outk2);

// -------------------------------------------------------------------------
// context pour le calcul des distributions
// -------------------------------------------------------------------------
class CtxtDistrib;
class CtxtDistrib
{
public:
	bool					m_IsUsed;
	bool					m_IsHomogeneous;
	bool					m_IsBegin;

	double					m_YearTerm;
	
	int m_nbnames;
	double					m_LossUnit;
	ARM_Vector				m_LossRates;
	int						m_discretizationstep;
	int						m_CopulaType;
	qIntegratorChoice		m_IntegrationMethod;
	int						m_IntegrationStep;
	qDISTRIB_TYPE			m_DistribType;
	double					m_TrancheDown;
	double					m_TrancheUp;
	double					m_TrancheDown_stepup;
	double					m_TrancheUp_stepup;

	bool					m_IsLongShort;

	vector<double> 			m_OutLosses;
	vector<double> 			m_ProbasLossesUp;
	vector<double> 			m_ProbasLossesDown;
	double					m_TailProbaLosseUp;
	double					m_TailProbaLosseDown;

	//cubes de stockage des probabilités conditionelles
	ICM_QCubix<double>*		m_ProbaCond_Up;	 //stockage des probas conditionnelles (nbnames,nbsteps,lup)
	ICM_QCubix<double>*		m_ProbaCond_Down; //stockage des probas conditionnelles (nbnames,nbsteps,ldown)

	//standard distrib parameters
	ARM_Vector m_pdef_maturity;
	ARM_Vector m_barriers_maturity;
	
	ARM_Vector m_beta_up_maturity;
	ARM_Vector m_beta_dw_maturity;

	std::vector<double> m_a_up_maturity;
	std::vector<double> m_b_up_maturity;

	std::vector<double> m_a_dw_maturity;
	std::vector<double> m_b_dw_maturity;

	//collat
	ARM_Vector m_collat_pdef_start;
	ARM_Vector m_collat_barriers_start;

	ARM_Vector m_collat_beta_up_start;
	ARM_Vector m_collat_beta_dw_start;

	std::vector<double> m_collat_a_up_start;
	std::vector<double> m_collat_b_up_start;

	std::vector<double> m_collat_a_dw_start;
	std::vector<double> m_collat_b_dw_start;

	//term structure review
	ARM_Vector m_ts_pdef_up_maturity_t1;
	ARM_Vector m_ts_pdef_up_maturity_t2;

	ARM_Vector m_ts_pdef_dw_maturity_t1;
	ARM_Vector m_ts_pdef_dw_maturity_t2;
	
	ARM_Vector m_ts_beta_up_maturity_t1;
	ARM_Vector m_ts_beta_up_maturity_t2;

	ARM_Vector  m_ts_barriers_up_maturity_t1;
	ARM_Vector  m_ts_barriers_up_maturity_t2;

	ARM_Vector m_ts_beta_dw_maturity_t1;
	ARM_Vector m_ts_beta_dw_maturity_t2;

	ARM_Vector   m_ts_barriers_dw_maturity_t1;
	ARM_Vector  m_ts_barriers_dw_maturity_t2;

	std::vector<double> m_ts_a_up_maturity_t1;
	std::vector<double> m_ts_b_up_maturity_t1;

	std::vector<double> m_ts_a_up_maturity_t2;
	std::vector<double> m_ts_b_up_maturity_t2;

	std::vector<double> m_ts_a_dw_maturity_t1;
	std::vector<double> m_ts_b_dw_maturity_t1;

	std::vector<double> m_ts_a_dw_maturity_t2;
	std::vector<double> m_ts_b_dw_maturity_t2;

	double	m_ts_t1_corr;
	double	m_ts_t2_corr;

	//collateral forward + ts

	ARM_Vector m_ts_collat_pdef_up_start_t1;
	ARM_Vector m_ts_collat_pdef_up_start_t2;

	ARM_Vector m_ts_collat_pdef_dw_start_t1;
	ARM_Vector m_ts_collat_pdef_dw_start_t2;

	ARM_Vector m_ts_collat_beta_up_start_t1;
	ARM_Vector m_ts_collat_beta_up_start_t2;

	ARM_Vector m_ts_collat_barriers_up_start_t1;
	ARM_Vector m_ts_collat_barriers_up_start_t2;

	ARM_Vector m_ts_collat_beta_dw_start_t1;
	ARM_Vector m_ts_collat_beta_dw_start_t2;

	ARM_Vector m_ts_collat_barriers_dw_start_t1;
	ARM_Vector m_ts_collat_barriers_dw_start_t2;

	std::vector<double> m_ts_collat_a_up_start_t1;
	std::vector<double> m_ts_collat_b_up_start_t1;

	std::vector<double> m_ts_collat_a_dw_start_t1;
	std::vector<double> m_ts_collat_b_dw_start_t1;

	std::vector<double> m_ts_collat_a_up_start_t2;
	std::vector<double> m_ts_collat_b_up_start_t2;

	std::vector<double> m_ts_collat_a_dw_start_t2;
	std::vector<double> m_ts_collat_b_dw_start_t2;

	double	m_ts_collat_start_t1_corr;
	double	m_ts_collat_start_t2_corr;

	//step up
	double	m_losses_reset_beg;
	double	m_losses_reset_end;

	//term structure review step up
	ARM_Vector m_ts_stepup_pdef_up_maturity_t1;
	ARM_Vector m_ts_stepup_pdef_up_maturity_t2;

	ARM_Vector m_ts_stepup_pdef_dw_maturity_t1;
	ARM_Vector m_ts_stepup_pdef_dw_maturity_t2;
	
	ARM_Vector m_ts_stepup_beta_up_maturity_t1;
	ARM_Vector m_ts_stepup_beta_up_maturity_t2;

	ARM_Vector m_ts_stepup_barriers_up_maturity_t1;
	ARM_Vector m_ts_stepup_barriers_up_maturity_t2;

	ARM_Vector m_ts_stepup_beta_dw_maturity_t1;
	ARM_Vector m_ts_stepup_beta_dw_maturity_t2;

	ARM_Vector m_ts_stepup_barriers_dw_maturity_t1;
	ARM_Vector m_ts_stepup_barriers_dw_maturity_t2;

	std::vector<double> m_ts_stepup_a_up_maturity_t1;
	std::vector<double> m_ts_stepup_b_up_maturity_t1;

	std::vector<double> m_ts_stepup_a_up_maturity_t2;
	std::vector<double> m_ts_stepup_b_up_maturity_t2;

	std::vector<double> m_ts_stepup_a_dw_maturity_t1;
	std::vector<double> m_ts_stepup_b_dw_maturity_t1;

	std::vector<double> m_ts_stepup_a_dw_maturity_t2;
	std::vector<double> m_ts_stepup_b_dw_maturity_t2;

	double	m_ts_stepup_t1_corr;
	double	m_ts_stepup_t2_corr;

	//Notionel variable
	vector<CtxtDistrib>	m_vn_contexts;

	//fast computation
	std::vector<double> m_Recoveries;
	std::vector<double> m_Notionals;
	double				m_TotalNotionals;


	CtxtDistrib() {Init();}

	~CtxtDistrib() ; 


	void Init() ;


	CtxtDistrib& operator= (const CtxtDistrib& ref) ;

	CtxtDistrib(const CtxtDistrib& i) ; 


	void ComputeBarriers(int& size,const ARM_Vector&,ARM_Vector& output);
	void ComputeALLBarriers();
	void ComputeCoefs(const ARM_Vector& barriers,
					  const ARM_Vector& betas,
					  std::vector<double>& coefsA,
					  std::vector<double>& coefsB);
	void ComputeALLCoefs();
	void ComputeAll(); 
	void IsHomogeneous()
	{ for(int k=0;k<m_LossRates.size();k++) if(m_LossRates[k]!=m_LossRates[0]) m_IsHomogeneous=false;}

	void ResetAll();
	
	friend bool operator > (const CtxtDistrib &, const CtxtDistrib &);
};


// -------------------------------------------------------------------------
// classe distrib
// -------------------------------------------------------------------------
class ICM_Distribution  : public ARM_Object  
{
//-----  Attributes
protected: //For compatibility only !!!!

	CtxtDistrib* its_ctxt;		// global context, association.

	bool   its_IsUp;			//barrier Up

	bool its_ishomogeneous;		//are the Loss Rates homogeneous ?
	bool itsTSR;				//Term structure mode

	unsigned int its_ind_name;		//Index for issuer name
	unsigned int its_lup;		// indice max de loss: dépend de la tranche up

	double  its_min_pdef;		// Min probabilty for default probability
	double	its_taildistrib;	// proba de loss superieur à la tranche up
	double	its_lossunit;		// loss_unit=(1-recov)*Nominal_One_Name (cadre homogène)

	// 17783 ICM_QCubix<double>* its_ProbCond_Perturb; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	// 17783 ICM_QMatrix<double>* its_lossdistrib_perturb;
	// 17783 ICM_QMatrix<double>* itsShifts;

	std::vector<double>	its_unique_beta;	//beta strike up (à its_t1) for maturity
	std::vector<double>	its_unique_beta_t2; //beta strike up (à its_t2) for maturity

	std::vector<double>	its_pdef_at_maturity;		// default probabilty at maturity (matu-)
	std::vector<double>	its_pdef_at_maturity_t2;	// default probabilty at maturity (matu+)
	std::vector<double>	its_pdef_perturb;			// default probabilty at fwd start (matu-)
	std::vector<double>	its_pdef_perturb_t2;		// default probabilty at fwd start (matu+)

	std::vector<double>	its_barrier;			// barriere associée à its_pdef_at_maturity
	std::vector<double>	its_barrier_t2;			// barriere associée à its_pdef_at_maturity_t2
	std::vector<double>	its_barrier_perturb;	// barriere associée à its_pdef_perturb
	std::vector<double>	its_barrier_perturb_t2; // barriere associée à its_pdef_perturb_t2

	std::vector<double> its_lossdistrib; // for assessor, m_dens[k] = P(N(T)=k)
	// 17783 std::vector<double> its_taildistrib_perturb; // for assessor, m_dens[k] = P(N(T)=k)
	
	std::vector<double> its_dbl_lossrates;	// for assessor, m_dens[k] = P(N(T)=k)
	std::vector<int>	its_int_lossrates;	// for assessor, m_dens[k] = P(N(T)=k)

	ICM_QCubix<double>* its_ProbCond; //stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l
	int its_ind_x;			//common factor
	unsigned int its_ind_Loss;		//Index for Losses
	unsigned int its_nbnames;		//Number of issuers

	// ----------------------------------------------------------------------
	//	LongShort
	// ----------------------------------------------------------------------
	
	unsigned int its_nbnames_short;			//Number of short issuers
	double its_LSminloss;					// sum of LS positives.
	std::vector<int> its_sorted_indices;	//Vector of sorted indices

	// ----------------------------------------------------------------------
	//collateral forward
	// ----------------------------------------------------------------------
	std::vector<double>	its_collat_fwd_unique_beta;		// beta strike up (à its_collat_t1) for start fwd
	std::vector<double>	its_collat_fwd_unique_beta_t2;	// beta strike up (à its_collat_t2) for start fwd

	// ----------------------------------------------------------------------
	//Term Structure
	// ----------------------------------------------------------------------

	double 	its_t1;			//date correl borne inf for maturity
	double	its_t2;			//date correl borne sup for maturity
	double	its_collat_t1;	//date correl borne inf for start fwd
	double	its_collat_t2;	//date correl borne sup for start fwd

	// ----------------------------------------------------------------------
	// StepUp
	// ----------------------------------------------------------------------
	ICM_QCubix<double>* its_ProbCond_stepup;	//stockage des probas conditionnelles pour un basket de dim k, de facteur commun x et de loss l

	int its_ind_x_stepup;						//common factor

//-----  Constructors/destructors
public:

	//Forward collateral
	static bool itsFwdCollatTScase;

	// For normalisation
	double its_normPi;
	double its_2Pi;
	double its_sqrt2;

	ICM_Distribution(const int& nbnames)
	{
		Init();
		Set(nbnames);
	}

	void Set(const int& nbnames) ;


	ICM_Distribution()
	{
		Init();
	}

	virtual ~ICM_Distribution() ;

	//----- Utilities

	void Init() ;

	virtual double compute_expectedlosstranche(double tranche_up, 
											   double tranche_down, 
											   double lossunit,
											   //ICM_QMatrix<double>* ShiftMatrix = NULL,
											   int TenorShift = -1,
											   int IssuerShift = -1)
	{return 0.;}

	// 17783 void SetShifts(ICM_QMatrix<double>* value) 
	// 17783 { 
	// 17783 	if (itsShifts)
	// 17783 		delete itsShifts;
	// 17783 	itsShifts = value; 
	// 17783 }

	// 17783 ICM_QMatrix<double>* GetShifts(void) { return itsShifts;}

	void View(char* id, FILE* ficOut);

	inline bool IsHomogeneous(void) { return its_ishomogeneous;}	
	inline void SetNbNames(const unsigned int& value) { its_nbnames = value;}	
	inline unsigned int GetNbNames(void) { return its_nbnames;}	
	inline void SetLup(const unsigned int& value) { its_lup = value;}	
	inline unsigned int GetLup(void) { return its_lup;}	
	inline void SetProbCond(ICM_QCubix<double>* value) 
	{
		if (its_ProbCond)
			delete its_ProbCond;
		its_ProbCond = (ICM_QCubix<double>*)value->Clone();
	}
	inline ICM_QCubix<double>* GetProbCond(void) { return its_ProbCond;}
	// 17783 inline void SetProbCond_Perturb(ICM_QCubix<double>* value) 
	// 17783 {
	// 17783 	if (its_ProbCond_Perturb)
	// 17783 		delete its_ProbCond_Perturb;
	// 17783 	its_ProbCond_Perturb = (ICM_QCubix<double>*)value->Clone();
	// 17783 }
	// 17783 inline ICM_QCubix<double>* GetProbCond_Perturb(void) { return its_ProbCond_Perturb;}
	// 17783 inline void SetLossDistrib_Perturb(ICM_QMatrix<double>* value) 
	// 17783 {
	// 17783 	if (its_lossdistrib_perturb)
	// 17783 		delete its_lossdistrib_perturb;
	// 17783 	its_lossdistrib_perturb = (ICM_QMatrix<double>*) value->Clone();
	// 17783 }
	// 17783 inline ICM_QMatrix<double>* GetLossDistrib_Perturb(void) { return its_lossdistrib_perturb;}
	inline const std::vector<double>& GetUniqueBeta(void) const {return its_unique_beta;} 
	inline const std::vector<double>& GetPdefAtMaturity(void) const {return its_pdef_at_maturity;} 
	// 17783 inline const std::vector<double>& GetPdefPerturb(void) const {return its_pdef_perturb;} 
	inline const std::vector<double>& GetBarrier(void) const {return its_barrier;}
	// 17783 inline const std::vector<double>& GetBarrierPerturb(void) const {return its_barrier_perturb;} 
	inline const std::vector<int>& GetIntLossRates(void) const {return its_int_lossrates;} 
	inline const std::vector<double>& GetDblLossRates(void) const {return its_dbl_lossrates;} 
	inline void SetIntLossRates(const ARM_Vector& value) 
	{
		its_int_lossrates.resize(value.size());
		its_dbl_lossrates.resize(value.size());
		// for(int k=0;k<value.size();k++) its_int_lossrates[k] = floor(round(value[k]));
		for(int k=0;k<value.size();k++) 
		{
			its_dbl_lossrates[k] = value[k];
			its_int_lossrates[k] = floor(value[k]);
		}
	}

	inline void check_homogeneous(const ARM_Vector& lossrates)
	{ for(int k=0;k<lossrates.size();k++) if(lossrates[k]!=lossrates[0]) its_ishomogeneous=false;}

	inline void SetPdefPerturb(const ARM_Vector& value) 
	{
		its_pdef_perturb.resize(value.size());
		for(int k=0;k<value.size();k++) its_pdef_perturb[k] = value[k];
	}

	inline void SetPdefAtMaturity(const ARM_Vector& value) 
	{
		its_pdef_at_maturity.resize(value.size());
		for(int k=0;k<value.size();k++) its_pdef_at_maturity[k] = value[k];
	}

	inline void SetUniqueBeta(const ARM_Vector& value) 
	{
		its_unique_beta.resize(value.size());
		for(int k=0;k<value.size();k++) its_unique_beta[k] = value[k];
	}

	void BitwiseCopy(const ARM_Object* src);
	void Copy(const ARM_Object* src);
	ARM_Object* Clone(void);

	inline void compute_min_pdef(const vector<double>& vpdef)
	{
	its_min_pdef=vpdef[0];
	for(int k=1;k<vpdef.size();k++) its_min_pdef=MIN(its_min_pdef,vpdef[k]);
	}

	virtual void compute_barrier();
	void compute_barrier_perturb();

	//Long Short
	/*********************************************************************/
	inline void ComputeNbNameShort(const ARM_Vector& lossrates)
	{ 
		its_nbnames_short=0;
		double sum=0; 
		for(int k=0;k<lossrates.size();k++) 
		{
			if(lossrates[k] < 0) { its_nbnames_short++; sum += lossrates[k]; }
		}
		its_LSminloss = sum; 
	}
	inline void SetSortedIndices(const std::vector<int>& value)
	{
		//Size test
		if (value.size() != its_nbnames)
		{
			its_sorted_indices.resize(0);
			return;
		}
		else
		{
			its_sorted_indices.resize(value.size());
			for(int k=0;k<value.size() ;k++) 
				its_sorted_indices[k] = value[k];
			return;
		}
	}
	/*********************************************************************/

	inline void SetFwdCollatTScase(bool status) {itsFwdCollatTScase=status;}
	inline bool GetFwdCollatTScase() {return itsFwdCollatTScase;}

	inline const std::vector<double>& GetUniqueBeta_t2(void) {return its_unique_beta_t2;} 

	inline const std::vector<double>& GetPdefAtMaturity_t2(void) {return its_pdef_at_maturity_t2;} 

	inline const std::vector<double>& GetBarrier_t2(void) {return its_barrier_t2;}

	inline bool IsTSR() {return itsTSR;}
	inline void SetTSR(bool value) {itsTSR=value;}

	inline void SetPdefAtMaturity_t2(const ARM_Vector& value) 
	{
		its_pdef_at_maturity_t2.resize(value.size());
		for(int k=0;k<value.size();k++) its_pdef_at_maturity_t2[k] = value[k];
	}

	inline void SetUniqueBeta_t2(const ARM_Vector& value) 
	{
		its_unique_beta_t2.resize(value.size());
		for(int k=0;k<value.size();k++) its_unique_beta_t2[k] = value[k];
	}

	inline double GetT1() {return its_t1;}
	inline double GetT2() {return its_t2;}

	inline void SetT1(const double& t1) {its_t1=t1;}
	inline void SetT2(const double& t2) {its_t2=t2;}

	inline double GetCollatT1() {return its_collat_t1;}
	inline double GetCollatT2() {return its_collat_t2;}

	inline void SetCollatT1(const double& t1) {its_collat_t1=t1;}
	inline void SetCollatT2(const double& t2) {its_collat_t2=t2;}

	inline void SetPdefPerturb_t2(const ARM_Vector&value) 
	{
		its_pdef_perturb_t2.resize(value.size());
		for(int k=0;k<value.size();k++) its_pdef_perturb_t2[k] = value[k];
	}

	inline void SetCollatUniqueBeta_t2(const ARM_Vector& value ) 
	{
		its_collat_fwd_unique_beta_t2.resize(value.size());
		for(int k=0;k<value.size();k++) its_collat_fwd_unique_beta_t2[k] = value[k];
	}

	inline void SetCollatUniqueBeta(const ARM_Vector& value) 
	{
		its_collat_fwd_unique_beta.resize(value.size());
		for(int k=0;k<value.size();k++) its_collat_fwd_unique_beta[k] = value[k];
	}

	inline void SetProbCond_stepup(ICM_QCubix<double>* value) 
	{
		if (its_ProbCond_stepup)
			delete its_ProbCond_stepup;
		its_ProbCond_stepup = (ICM_QCubix<double>*)value->Clone();
	}
	inline ICM_QCubix<double>* GetProbCond_stepup(void) { return its_ProbCond_stepup;}

	void SetCtxt(CtxtDistrib* ctxt) {ctxt->m_IsUsed=true;its_ctxt = ctxt;}
	CtxtDistrib* GetCtxt() {return its_ctxt;}

	// 17783 inline const std::vector<double>& GetPdefPerturb_t2(void) {return its_pdef_perturb_t2;} 

};
#endif 
