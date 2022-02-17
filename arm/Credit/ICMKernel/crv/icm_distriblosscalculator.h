
#ifndef _ICM_DISTRIBLOSS_CALCULATOR_H_ 
#define _ICM_DISTRIBLOSS_CALCULATOR_H_ 

#include "ICMKernel/inst/icm_security.h"
#include "ICMKernel/pricer/icm_pricer.h"
class ICM_ModelMultiCurves ;
class ICM_DefaultCurve; 
class ICM_Correlation; 
class ICM_Credit_Index;
//	1st stage... 
class _barrier_ctxt
{
	public:
	double m_barrier;
	double m_correl;
	double m_pdef;
	_barrier_ctxt(){m_pdef=m_barrier=m_correl=0.;}
};

/** static void __stdcall objfun_barrier(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);
				   **/ 

//	Signature should not change.. 

double cpt_elt_default(const double& yf,ICM_Pricer& , ARM_Model& ,ARM_Security&,vector<double>& ) ;
double cpt_elt_pricer_distrib(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_fullhomogeneous(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_fullhomogeneous_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
// 17783 double cpt_elt_pricer_distrib_fast(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_); 
// 17783 double cpt_elt_pricer_distrib_smile_fast(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_); 
double cpt_elt_pricer_distrib_MF(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_) ;

//Fwd start collateral
double cpt_elt_pricer_distrib_smile_collat_fwd(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_collat_fwd_TSR(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd_TSR(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;

//stepup subordination
double cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;
double cpt_elt_pricer_distrib_smile_stepup_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses) ;

//Variable notional
double cpt_elt_pricer_distrib_smile_VN(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses);
double cpt_elt_pricer_distrib_smile_APPROX(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses);

// 2F for Sectorial CDO Square
double cpt_elt_pricer_distrib_sector(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_) ;

double cpt_quick_el_full_homog(double pdefault,double recovery,int nbnames,double strikedw,
							   double strikeup,double correldw,double correlup,vector<double>& losses,
							   int intstep=40,bool LHP=false);

double FastCDOPricing(ARM_Vector* YFstartdates,ARM_Vector* YFenddates,double rate,
					  double strikedw,double strikeup,double correldw,double correlup,
					  double notfeeleg,double notdefleg,vector<double>& pdef,
					  vector<double>& disc,double recovery,int nbnames,int intstep,
					  double& feepv,double& defpv,bool LHP = false);

double FastCDOPricing(ARM_Date& AsOf,ARM_Date& startdate,ARM_Date& enddate,
					  int frequency,double rate, double strikedw,double strikeup,
					  double correldw, double correlup, double notfeeleg,double notdefleg,
					  ICM_DefaultCurve* pdef, ARM_ZeroCurve* disc, double recovery,
					  int nbnames, int intstep, double& feepv, double& defpv, bool LHP);

double FastCDOPricing(ARM_Vector& schedule, const int& begin,const int& end,  double rate,
					  const vector<double>& el,double notfeeleg,double notdefleg, 
					  ARM_ZeroCurve* zc, double& feepv,double& defpv);

double CptPtf_ELoss(ICM_Pricer* pricer,int payfreq = K_QUARTERLY);
double CptPtf_ELoss_Index(ICM_ModelMultiCurves* mod,ICM_Credit_Index* index ,ARM_Date& Maturity,int payfreq = K_QUARTERLY);
double CptPtf_ELoss_Index(ICM_ModelMultiCurves* mod,const std::string& indexlabel,ARM_Date& Maturity,int payfreq = K_QUARTERLY);

//term structure review
double cpt_quick_el_full_homog_TSR(double pdefault_T1,double pdefault_T2,double recovery,int nbnames,double strikedw,
							   double strikeup,double correldw_T1,double correldw_T2,double correlup_T1,double correlup_T2,vector<double>& losses,
							   int intstep=40,bool LHP=false);

double FastCDOPricing_TSR(ARM_Vector* YFstartdates,ARM_Vector* YFenddates,double rate,
					  double strikedw,double strikeup,ARM_VolCurve* correl,
					  double notfeeleg,double notdefleg,vector<double>& pdef_T1,vector<double>& pdef_T2,
					  vector<double>& disc,double recovery,int nbnames,int intstep,
					  double& feepv,double& defpv,bool LHP = false);

double FastCDOPricing_TSR(ARM_Date& AsOf,ARM_Date& startdate,ARM_Date& enddate,
					  int frequency,double rate, double strikedw,double strikeup,
					  double notfeeleg,double notdefleg,ICM_DefaultCurve* pdef, 
					  ARM_ZeroCurve* disc, double recovery,int nbnames,ARM_VolCurve* correl,
					   int intstep, double& feepv, double& defpv, bool LHP);

double ComputeBarrierTSR(ICM_Correlation* correl,
						 const ICM_DefaultCurve* defcurve,
						 const double& yf_t1,
						 const double& yf_t2,
						 const double& yf,
						 const qCorrel_By_Strike& striketype);

#endif // _ICM_DISTRIBLOSS_CALCULATOR_H_