#include "ARMKernel/glob/firsttoinc.h" 



#include "icm_distriblosscalculator.h" 
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"

#include <ICMKernel/glob/icm_maths.h>
#include "ICMKernel\util\icm_RootFinderND.h"


// --------------------------------------------------------------------------------------------
// loss distribution for collateral with base correlation computation : Term Structure Review
// --------------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	int i=0;
	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Mez& ftd=dynamic_cast<ICM_Mez&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	//case step up subordination
	if (isStepUp)
	{	return cpt_elt_pricer_distrib_smile_TSR(yf,pricer_,model_,sec_,losses); }

	// -----------------------------------------------------------
	// JUST TO TAKE INTO ACCOUNT A SECTORIAL APPROACH
	ICM_Correlation* tmp_correlation = model.GetCorrelation();

	if (tmp_correlation->GetName() == ICM_CORRELATION_SECTOR)
		return	cpt_elt_pricer_distrib_sector(yf, pricer_, model_, sec_);
		//return	cpt_elt_pricer_distrib_MF(yf, pricer_, model_, sec_);
	// -----------------------------------------------------------

	//case full homogeneous 
	if (collat->IsFullHomog()) 
		return cpt_elt_pricer_distrib_smile_fullhomogeneous_TSR(yf,pricer_,model_,sec_,losses) ; 
		// return pricer.ExpectedLossTrancheFullHomogeneous(yf);
	if (yf<=0) return 0.;

	ARM_Date Maturity = ftd.GetEndDateNA(); 
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_LJ* Distrib = NULL;

	if (!collat->SumNotionals(ftd.GetStartDateNA())) return 0.;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	if (tranche_down>=tranche_up) return 0.; 

	int nbnames = collat->GetNbIssuers();

	// double* beta = NULL;
	// ARM_Vector* VBeta = NULL;
	// beta = new double[nbnames];	
	ARM_Vector beta(nbnames); 

	double ELT = 0.;

	// int integrationStep1; 
	// pricer.GetIntegrationStep1(integrationStep1) ;
	qIntegratorChoice	TheIntegratorType;
	if (! pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		// ntegrationStep1 = 40;
		pricer.SetIntegrationStep1(40) ;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());
	
	double	current_Pdefault_t1=0.,current_Pdefault_t2=0.;
	//double Pdefault_t1=0.,Pdefault_t2=0.;
	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;

	DeduceYearTermsForTSR(tmp_correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);

	/** double* Pdefault_t1 = new double[nbnames];
	double* Pdefault_t2 = new double[nbnames];
	double* Pdefault_down_t1 = new double[nbnames];
	double* Pdefault_down_t2 = new double[nbnames]; **/ 
	ARM_Vector Pdefault_t1 (nbnames); 
	ARM_Vector Pdefault_t2 (nbnames); 
	ARM_Vector Pdefault_down_t1 (nbnames); 
	ARM_Vector Pdefault_down_t2 (nbnames); 

  	// Smile Parameters
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	ICM_Correlation* correlation = model.GetCorrelation();
	for (i=0; i<nbnames; i++) 
	{	
		Pdefault_t1[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt1_pdef,qStrike_UP);
		Pdefault_t2[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt2_pdef,qStrike_UP);
		Pdefault_down_t1[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt1_pdef,qStrike_LOW);
		Pdefault_down_t2[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt2_pdef,qStrike_LOW);
	}

	if (pricer.GetDistribution() == 0)
	{

			Distrib = new ICM_Gauss1FLossDistrib_LJ(yt1_corr,yt2_corr,nbnames,
				Pdefault_t1,Pdefault_t2,Pdefault_down_t1,Pdefault_down_t2,beta,beta,
				// pricer.getLossRate(),
				lossUnits.getLossRates(date),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1());

			pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );

		Distrib->Set(yt1_corr,yt2_corr,nbnames,Pdefault_t1,Pdefault_t2,Pdefault_down_t1,Pdefault_down_t2,beta,beta,
				lossUnits.getLossRates(date),pricer.GetIntegrationStep1(),pricer.GetCopulaType(),
				TheIntegratorType,pricer.GetIntegrationStep1());
	}

	// if (Pdefault_t1) delete[] Pdefault_t1;
	// if (Pdefault_t2) delete[] Pdefault_t2;
	// if (Pdefault_down_t1) delete[] Pdefault_down_t1;
	// if (Pdefault_down_t2) delete[] Pdefault_down_t2;

  	// Smile Parameters
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector corr_low_t1(nbnames);
	ARM_Vector  corr_Hight_t1(nbnames);
	ARM_Vector  corr_low_t2(nbnames);
	ARM_Vector  corr_Hight_t2(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_NONE);
	}

	correlation->SetForcedStrikeType(qStrike_NONE);
	Distrib->UpdateCorrelation_TSR(corr_low_t1,corr_low_t2,corr_Hight_t1,corr_Hight_t2);

	// New...
	Distrib->SetCopulaType(qGAUSSIAN);
	Distrib->SetIntegratorType(TheIntegratorType,pricer.GetIntegrationStep1());
	Distrib->SetIntegrationStep(pricer.GetIntegrationStep1());


	try
	{
			vector<double> ProbaslossesUp;
			vector<double> ProbaslossesDown;

			ELT =Distrib->compute_expectedlosstranche(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date),
// 17783 															pricer.GetTenorShift(),
// 17783 															pricer.GetIssuerShift(),
													losses,ProbaslossesUp,ProbaslossesDown); 

	}
	catch (std::exception&e)
	{	ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());	}
	catch(...)
	{	ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib_smile :: Unable to compute Loss Distribution for maturity="<<yf);	}
	
	if (ELT<0.) ELT = 0.;
	// if (beta) delete[] beta;
	
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}


// --------------------------------------------------------------------------------------------
// Loss distribution with base correlation computation & Fwd Start Collat : Term Structure Review
// --------------------------------------------------------------------------------------------
double cpt_elt_pricer_distrib_smile_collat_fwd_TSR(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_, vector<double>& losses)
{
	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile_Collat_Fwd&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
ICM_Collateral*collat =ftd.GetCollateral(); 
const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 
	
	if ( (yf<=0)||((yf-yf_fwd_start)<=0) ) return 0.;

	if (collat->IsFullHomog())
	{ return cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd(yf,yf_fwd_start,pricer_,model_,sec_,losses);}
	
	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_LJ* Distrib = NULL;

	if (!collat->SumNotionals(ftd.GetStartDateNA())) return 0.;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	if (tranche_down>=tranche_up) return 0.; 

	int nbnames = collat->GetNbIssuers();

	// double* beta = NULL;
	// ARM_Vector* VBeta = NULL;
	// beta = new double[nbnames];	
	ARM_Vector beta(nbnames); 

	double ELT = 0.;

	qIntegratorChoice	TheIntegratorType;
	if (! pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		pricer.SetIntegrationStep1(40) ;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());
	
	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;

	double start_beta_down_t1=0.,start_beta_up_t1=0.,start_beta_down_t2=0.,start_beta_up_t2=0.;
	double start_yt1_corr=0.,start_yt2_corr=0.;
	double start_yt1_pdef=0.,start_yt2_pdef=0.;

	DeduceYearTermsForTSR(correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);
	DeduceYearTermsForTSR(correlation,yf_fwd_start,start_yt1_corr,start_yt2_corr,start_yt1_pdef,start_yt2_pdef);

	/** double* Pdefault_t1 = new double[nbnames];
	double* PdefaultStart_t1 = new double[nbnames];
	double* Pdefault_t2 = new double[nbnames];
	double* PdefaultStart_t2 = new double[nbnames]; **/ 
	ARM_Vector Pdefault_t1(nbnames);
	ARM_Vector PdefaultStart_t1(nbnames);
	ARM_Vector Pdefault_t2(nbnames);
	ARM_Vector PdefaultStart_t2(nbnames);

	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<nbnames; i++) 
	{
		//Estimation des Probas de défaut sur la période forward
		Pdefault_t1[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yt1_pdef);
		PdefaultStart_t1[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(start_yt1_pdef);

		Pdefault_t2[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yt2_pdef);
		PdefaultStart_t2[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(start_yt2_pdef);

		if (Pdefault_t1[i]>1. || PdefaultStart_t1[i]>1. || Pdefault_t2[i]>1. || PdefaultStart_t2[i]>1.)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
		if (Pdefault_t1[i]<0.0 || PdefaultStart_t1[i]<0.0 || Pdefault_t2[i]<0.0 || PdefaultStart_t2[i]<0.0)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
	}

	if (pricer.GetDistribution() == 0)
	{
		Distrib = new ICM_Gauss1FLossDistrib_LJ(yt1_corr,yt2_corr,start_yt1_corr,start_yt2_corr,
				nbnames,Pdefault_t1,Pdefault_t2,PdefaultStart_t1,PdefaultStart_t2,beta,beta,beta,beta,
				lossUnits.getLossRates(date),pricer.GetIntegrationStep1(),pricer.GetCopulaType(),
				TheIntegratorType,pricer.GetIntegrationStep1());

		pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );
		Distrib->Set(yt1_corr,yt2_corr,start_yt1_corr,start_yt2_corr,
				nbnames,Pdefault_t1,Pdefault_t2,PdefaultStart_t1,PdefaultStart_t2,beta,beta,beta,beta,
				lossUnits.getLossRates(date),pricer.GetIntegrationStep1(),pricer.GetCopulaType(),
				TheIntegratorType,pricer.GetIntegrationStep1());
	}

	// if (Pdefault_t1) delete[] Pdefault_t1;
	// if (Pdefault_t2) delete[] Pdefault_t2;
	// if (PdefaultStart_t1) delete[] PdefaultStart_t1;
	// if (PdefaultStart_t2) delete[] PdefaultStart_t2;

	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector corr_low_t1(nbnames);
	ARM_Vector corr_Hight_t1(nbnames);
	ARM_Vector corr_low_t2(nbnames);
	ARM_Vector corr_Hight_t2(nbnames);

	ARM_Vector corr_collat_low_t1(nbnames);
	ARM_Vector corr_collat_Hight_t1(nbnames);
	ARM_Vector corr_collat_low_t2(nbnames);
	ARM_Vector corr_collat_Hight_t2(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_collat_low_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),start_yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_collat_Hight_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),start_yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_collat_low_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),start_yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_collat_Hight_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),start_yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_NONE);
	}
	correlation->SetForcedStrikeType(qStrike_NONE);

	Distrib->Collat_fwd_UpdateCorrelation_TSR(corr_low_t1,corr_low_t2,corr_Hight_t1,corr_Hight_t2,
										corr_collat_low_t1,corr_collat_low_t2,corr_collat_Hight_t1,corr_collat_Hight_t2);
	//Distrib->precompute_coeffs_perturb(); //Update aussi les coef Perturb pour modifier les coef Start

	// New...
	Distrib->SetCopulaType(qGAUSSIAN);
	Distrib->SetIntegratorType(TheIntegratorType,pricer.GetIntegrationStep1());
	Distrib->SetIntegrationStep(pricer.GetIntegrationStep1());

	try
	{
		//Do not deal with Long Short
		if (ftd.GetCollateral()->IsLongShort())
		{
			ELT =0.;
		}
		else
		{
			ELT =Distrib->compute_expectedlosstranche_CollatFwd(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date)
													); 

		}
	}
	catch (std::exception&e)
	{ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distribCollatFwd:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());	}
	catch(...)
	{ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib_smileCollatFwd:: Unable to compute Loss Distribution for maturity="<<yf);}
	
	if (ELT<0.) ELT = 0.;
// 	if (beta) delete[] beta;

	ELT /= (tranche_up-tranche_down);

	return (ELT);
}

static void __stdcall objfun_barrier(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	double result = 0.;
    _barrier_ctxt* p = (_barrier_ctxt*)comm->p ;

	result = NAG_cumul_normal(-p->m_barrier);
	result -= NAG_cumul_normal(x[0]);
	result += NAG_bivariate_normal_dist(p->m_barrier,x[0],p->m_correl);
	result -= (1. - p->m_pdef);

	*objf = result*result;
}

double ComputeBarrierTSR(ICM_Correlation* correl,
						 const ICM_DefaultCurve* defcurve,
						 const double& yf_t1,
						 const double& yf_t2,
						 const double& yf,
						 const qCorrel_By_Strike& striketype)
{
	using namespace OptimTools;
	double output = 0.;

	vector<double> X;X.resize(1);
	vector<double> bound_inf_X;bound_inf_X.resize(1);bound_inf_X[0]=-10.;
	vector<double> bound_sup_X;bound_sup_X.resize(1);bound_sup_X[0]=10.;

	double barrier=defcurve->DefaultProba(yf);

	if (yf_t2==0.)
		{return barrier;}

	ARM_Vector BaseCorrelSchedule; 
	correl->GetCorrelationTerms(BaseCorrelSchedule);

	double icorr0 = BaseCorrelSchedule.Elt(0);

	for (int i=0;i<BaseCorrelSchedule.GetSize();i++)
	{
		double yt_bc=BaseCorrelSchedule.Elt(i);
		double yt_bc_prev=0.;
		if (i)
			{yt_bc_prev=BaseCorrelSchedule.Elt(i-1);}

		correl->SetForcedStrikeType(striketype);
		double beta_t1 = correl->GetBeta("NOWAY",yt_bc_prev,CREDIT_DEFAULT_VALUE,yf);
		double beta_t2 = correl->GetBeta("NOWAY",yt_bc,CREDIT_DEFAULT_VALUE,yf);

		double rho = RHO_DIST(yt_bc_prev,yt_bc);
		//rho=1.;

		_barrier_ctxt c;
		c.m_barrier=barrier;
		c.m_correl=rho*beta_t1*beta_t2+rho*sqrt((1.-beta_t1*beta_t1)*(1.-beta_t2*beta_t2));
		//bound_inf_X[0]=barrier;

		if (yt_bc>yf_t2)
		{break;}
		else if (yt_bc<=yf_t1)
		{
			if (yt_bc_prev &&  (!CHECK_EQUAL(yf,icorr0)) ){
			c.m_pdef=defcurve->DefaultProba(yt_bc);
			//X[0]=(NAG_deviates_normal_dist(c.m_pdef)+bound_sup_X[0])/2.;
			X[0]=MAX(NAG_deviates_normal_dist(c.m_pdef),barrier);

			output = OptimTools::NagMultiOptimisator((void*) &c,
												objfun_barrier,
												X,
												bound_inf_X,
												bound_sup_X,1.E-6,200);

			#ifdef _DEBUG
			if (X[0]<barrier)
			{	ICMLOG("yf_corr_t1: "<<yf_t1<<"yf_corr_t2: "<<yf_t2<<" yt: "<<yf);
				ICMLOG("Searching barrier : "<<"barrier++="<<X[0]<<" < barrier="<<barrier);}
			#endif

			barrier = X[0];}
			else if (!CHECK_EQUAL(yf,icorr0))
			{barrier=NAG_deviates_normal_dist(defcurve->DefaultProba(yt_bc));}
			else
			{	// yf = maturity min correl
				return defcurve->DefaultProba(yf);
			}
		}
		else if (!CHECK_EQUAL(yf,icorr0))
		{
			c.m_pdef=defcurve->DefaultProba(yf);
			//X[0]=(NAG_deviates_normal_dist(c.m_pdef)+bound_sup_X[0])/2.;
			X[0]=MAX(NAG_deviates_normal_dist(c.m_pdef),barrier);

			output = OptimTools::NagMultiOptimisator((void*) &c,
												objfun_barrier,
												X,
												bound_inf_X,
												bound_sup_X,1.E-6,200);

			#ifdef _DEBUG
			if (X[0]<barrier)
			{	ICMLOG("yf_corr_t1: "<<yf_t1<<"yf_corr_t2: "<<yf_t2<<" yt: "<<yf);
				ICMLOG("Searching barrier : "<<"barrier++="<<X[0]<<" < barrier="<<barrier);}
			#endif

			barrier = X[0];
		}
		else
		{   // yf = maturity min correl
			return defcurve->DefaultProba(yf);}
	}

	return (NAG_cumul_normal(barrier));
}


// --------------------------------------------------------------------------------------------
// Step Up subordination
// loss distribution for collateral with base correlation computation : Term Structure Review
// --------------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_stepup_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	int i=0;
	if (yf<=0) return 0.;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Mez& ftd=dynamic_cast<ICM_Mez&>(sec_); 
ICM_Collateral*collat =ftd.GetCollateral(); 
const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	// -----------------------------------------------------------
	// JUST TO TAKE INTO ACCOUNT A SECTORIAL APPROACH
	ICM_Correlation* tmp_correlation = model.GetCorrelation();

	if (tmp_correlation->GetName() == ICM_CORRELATION_SECTOR)
		return	cpt_elt_pricer_distrib_sector(yf, pricer_, model_, sec_);

	//case full homogeneous 
	if (collat->IsFullHomog()) 
		return cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(yf,pricer_,model_,sec_,losses) ; 

	ARM_Date Maturity = ftd.GetEndDateNA(); 
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_LJ* Distrib = NULL;

	if (!collat->SumNotionals(ftd.GetStartDateNA())) return 0.;

	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	vector<double> tranche_down;tranche_down.resize(stepup_dates.size());
	vector<double> tranche_up;tranche_up.resize(stepup_dates.size());

	for (i=0;i<stepup_dates.size();i++)
	{
		tranche_down[i]= ABS(round(pricer.GetTranche_Down(stepup_dates[i])));
		tranche_up[i]= ABS(round(pricer.GetTranche_Up(stepup_dates[i])));
		if (tranche_down[i]>=tranche_up[i]) return 0.; 
	}

	int nbnames = collat->GetNbIssuers();

	// double* beta = NULL;
	// ARM_Vector* VBeta = NULL;
	// beta = new double[nbnames];	

	double ELT = 0.;

	qIntegratorChoice	TheIntegratorType;
	if (! pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		pricer.SetIntegrationStep1(40) ;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());
	
	double	current_Pdefault_t1=0.,current_Pdefault_t2=0.;
	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;
	double yt_stepup = 0.;
	double yt1_corr_stepup=0.,yt2_corr_stepup=0.;
	double yt1_pdef_stepup=0.,yt2_pdef_stepup=0.;


	DeduceYearTermsForTSR(tmp_correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);
	if (stepup_dates.size()>1)
	{	
		yt_stepup = (stepup_dates[stepup_dates.size()-2]-model.GetStartDate().GetJulian())/365.;
		DeduceYearTermsForTSR(tmp_correlation,yt_stepup,yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,yt2_pdef_stepup);	
	}

	if (CHECK_EQUAL(yf,yt_stepup))
	{isStepUp=false;}

	/** double* Pdefault_t1 = new double[nbnames];
	double* Pdefault_t2 = new double[nbnames];
	double* Pdefault_down_t1 = new double[nbnames];
	double* Pdefault_down_t2 = new double[nbnames];

	double* Pdefault_t1_stepup = new double[nbnames];
	double* Pdefault_t2_stepup = new double[nbnames];
	double* Pdefault_down_t1_stepup = new double[nbnames];
	double* Pdefault_down_t2_stepup = new double[nbnames];
	**/ 
	ARM_Vector Pdefault_t1(nbnames); 
	ARM_Vector Pdefault_t2(nbnames); 
	ARM_Vector Pdefault_down_t1(nbnames); 
	ARM_Vector Pdefault_down_t2(nbnames); 

	ARM_Vector Pdefault_t1_stepup(nbnames); 
	ARM_Vector Pdefault_t2_stepup(nbnames); 
	ARM_Vector Pdefault_down_t1_stepup(nbnames); 
	ARM_Vector Pdefault_down_t2_stepup(nbnames); 

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	const std::vector<int> &collatRank = lossUnits.getCollatRank(Maturity); 
	for (i=0; i<nbnames; i++) 
	{	
		Pdefault_t1[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt1_pdef,qStrike_UP);
		Pdefault_t2[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt2_pdef,qStrike_UP);
		Pdefault_down_t1[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt1_pdef,qStrike_LOW);
		Pdefault_down_t2[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr,yt2_corr,yt2_pdef,qStrike_LOW);
		Pdefault_t1_stepup[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,qStrike_UP);
		Pdefault_t2_stepup[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr_stepup,yt2_corr_stepup,yt2_pdef_stepup,qStrike_UP);
		Pdefault_down_t1_stepup[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,qStrike_LOW);
		Pdefault_down_t2_stepup[i]=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i])),yt1_corr_stepup,yt2_corr_stepup,yt2_pdef_stepup,qStrike_LOW);
	}

	CtxtDistrib c;
	c.m_nbnames = collat->GetNbIssuers();
	c.m_DistribType=qDISTRIB_STD_TSR;
	c.m_LossUnit=lossUnits.getLossUnit(Maturity); ;
	c.m_IntegrationStep=pricer.GetIntegrationStep1();
	c.m_TrancheDown=tranche_down[tranche_down.size()-3];
	c.m_TrancheUp=tranche_up[tranche_up.size()-3];
	c.m_TrancheDown_stepup=tranche_down[tranche_down.size()-1];
	c.m_TrancheUp_stepup=tranche_up[tranche_up.size()-1];

	c.m_ts_pdef_up_maturity_t1=Pdefault_t1;
	c.m_ts_pdef_up_maturity_t2=Pdefault_t2;
	c.m_ts_pdef_dw_maturity_t1=Pdefault_down_t1;
	c.m_ts_pdef_dw_maturity_t2=Pdefault_down_t2;

	c.m_ts_t1_corr=yt1_corr;
	c.m_ts_t2_corr=yt2_corr;
	c.m_LossRates=lossUnits.getLossRates(Maturity);
	c.m_IntegrationMethod=TheIntegratorType;
	c.m_CopulaType=pricer.GetCopulaType();

	c.m_ts_beta_dw_maturity_t1=std::vector<double>(c.m_nbnames,0.);
	c.m_ts_beta_dw_maturity_t2=std::vector<double>(c.m_nbnames,0.);
	c.m_ts_beta_up_maturity_t1=std::vector<double>(c.m_nbnames,0.);
	c.m_ts_beta_up_maturity_t2=std::vector<double>(c.m_nbnames,0.);

	c.m_ts_stepup_pdef_up_maturity_t1=Pdefault_t1_stepup;
	c.m_ts_stepup_pdef_up_maturity_t2=Pdefault_t2_stepup;
	c.m_ts_stepup_pdef_dw_maturity_t1=Pdefault_down_t1_stepup;
	c.m_ts_stepup_pdef_dw_maturity_t2=Pdefault_down_t2_stepup;

	if (pricer.GetDistribution() == 0)
	{
		Distrib = new ICM_Gauss1FLossDistrib_LJ(&c);
		pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );
		Distrib->Set(&c);
	}

	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector  corr_low_t1(nbnames);
	ARM_Vector  corr_Hight_t1(nbnames);
	ARM_Vector  corr_low_t2(nbnames);
	ARM_Vector  corr_Hight_t2(nbnames);

	std::vector<double> corr_low_t1_stepup;corr_low_t1_stepup.resize(nbnames);
	std::vector<double> corr_Hight_t1_stepup;corr_Hight_t1_stepup.resize(nbnames);
	std::vector<double> corr_low_t2_stepup;corr_low_t2_stepup.resize(nbnames);
	std::vector<double> corr_Hight_t2_stepup;corr_Hight_t2_stepup.resize(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t1[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t2[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_NONE);

		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t1_stepup[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t1_stepup[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt1_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low_t2_stepup[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight_t2_stepup[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),yt2_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
		correlation->SetForcedStrikeType(qStrike_NONE);
	}

	correlation->SetForcedStrikeType(qStrike_NONE);
	Distrib->UpdateCorrelation_TSR(corr_low_t1,corr_low_t2,corr_Hight_t1,corr_Hight_t2);

	// New...
	Distrib->SetCopulaType(qGAUSSIAN);
	Distrib->SetIntegratorType(TheIntegratorType,pricer.GetIntegrationStep1());
	Distrib->SetIntegrationStep(pricer.GetIntegrationStep1());

	c.m_ts_beta_dw_maturity_t1=corr_low_t1;
	c.m_ts_beta_dw_maturity_t2=corr_low_t2;
	c.m_ts_beta_up_maturity_t1=corr_Hight_t1;
	c.m_ts_beta_up_maturity_t2=corr_Hight_t2;

	c.m_ts_stepup_beta_dw_maturity_t1=corr_low_t1_stepup;
	c.m_ts_stepup_beta_dw_maturity_t2=corr_low_t2_stepup;
	c.m_ts_stepup_beta_up_maturity_t1=corr_Hight_t1_stepup;
	c.m_ts_stepup_beta_up_maturity_t2=corr_Hight_t2_stepup;

	if (!isStepUp) //cas non stepup dans l'intervalle de temps [0;T0[
	{
		try
		{	ELT =Distrib->ComputeEL(&c);
			losses = c.m_OutLosses;	}
		catch (std::exception&e)
		{	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case : "<<e.what());	}
		catch (...)
		{	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");	}

		ELT /= (c.m_TrancheUp-c.m_TrancheDown);
	}
	else 
	{// cas stepup
		ELT = Distrib->ComputeStepUp(&c);
	}

	if (ELT<0.) ELT = 0.;
	//if (beta) delete[] beta;
/*
	if (Pdefault_t1) delete[] Pdefault_t1;
	if (Pdefault_t2) delete[] Pdefault_t2;
	if (Pdefault_down_t1) delete[] Pdefault_down_t1;
	if (Pdefault_down_t2) delete[] Pdefault_down_t2;
	if (Pdefault_t1_stepup) delete[] Pdefault_t1_stepup;
	if (Pdefault_t2_stepup) delete[] Pdefault_t2_stepup;
	if (Pdefault_down_t1_stepup) delete[] Pdefault_down_t1_stepup;
	if (Pdefault_down_t2_stepup) delete[] Pdefault_down_t2_stepup;
*/
	return (ELT);
}
