#include "ARMKernel/glob/firsttoinc.h" 
#include "icm_distriblosscalculator.h" 


#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "ICMKernel\inst\icm_ftd.h"

// ------------------------------------------------------------------
// loss distribution for collateral with base correlation computation
// ------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_APPROX(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	int k=0;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd & ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	ARM_Date Maturity = ftd.GetEndDateNA(); 
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_LJ* Distrib = NULL;

	if (!ftd.GetCollateral()->SumNotionals(ftd.GetStartDateNA())) return 0.;

	ARM_Date AsOf = (ARM_Date)(model.GetStartDate().GetJulian());
	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down)); 
	tranche_up = ABS(round(tranche_up));

	if (tranche_down>=tranche_up) return 0.; 

	int nbnames = ftd.GetCollateral()->GetNbIssuers();

	// double* beta = NULL;
	// ARM_Vector* VBeta = NULL;
	// beta = new double[nbnames];	
	ARM_Vector beta (nbnames); 

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
	
	double	current_Pdefault;

	// double* Pdefault = new double[nbnames];
	ARM_Vector Pdefault(nbnames); 
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	std::vector<double> Recoveries; Recoveries.resize(nbnames);
	std::vector<double> Notionals;Notionals.resize(nbnames);

	for (int i=0; i<nbnames; i++) 
	{	
		current_Pdefault	=	model.GetDefaultCurve(ftd.GetCollateral()->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);
		
		if (CHECK_EQUAL(current_Pdefault, 0.0))
		{Pdefault[i] = 0.0;	}
		else if (CHECK_EQUAL(current_Pdefault, 1.0))
		{Pdefault[i] = 1.0;}
		else if (current_Pdefault > 1.)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<ftd.GetCollateral()->GetIssuersLabels(collatRank[i]));}
		else if (current_Pdefault < 0.0)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<ftd.GetCollateral()->GetIssuersLabels(collatRank[i]));}
		else Pdefault[i] = current_Pdefault;

		Recoveries[i]= model.GetDefaultCurve(ftd.GetCollateral()->GetIssuersLabels(collatRank[i]))->GetRecovery();
		Notionals[i] = ftd.GetCollateral()->GetIssuersNotional(collatRank[i]);
	}

	CtxtDistrib c_;
	c_.m_IsUsed=true;
	c_.m_DistribType=qDISTRIB_APPROX;
	c_.m_nbnames=nbnames;
	c_.m_LossUnit=lossUnits.getLossUnit(AsOf);
	c_.m_TrancheDown=tranche_down;
	c_.m_LossRates=lossUnits.getLossRates(AsOf);
	c_.m_CopulaType=pricer.GetCopulaType();
	c_.m_IntegrationMethod=TheIntegratorType;
	c_.m_TrancheUp=tranche_up;
	c_.m_pdef_maturity=Pdefault;
	c_.m_IntegrationStep=pricer.GetIntegrationStep1();
	c_.m_Recoveries=Recoveries;
	c_.m_Notionals=Notionals;
	ARM_Vector Nots(Notionals.size(),Notionals.begin());
	c_.m_TotalNotionals = Nots.sum();

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	std::vector<double> corr_low;corr_low.resize(nbnames);
	std::vector<double> corr_Hight;corr_Hight.resize(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low[i] = correlation->GetBeta(ftd.GetCollateral()->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight[i] = correlation->GetBeta(ftd.GetCollateral()->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	}

	c_.m_beta_dw_maturity=corr_low;
	c_.m_beta_up_maturity=corr_Hight;

	CtxtDistrib c=c_;
	c_.m_pdef_maturity.Resize(0) ; 

	correlation->SetForcedStrikeType(qStrike_NONE);

	if (pricer.GetDistribution() == 0)
	{
		Distrib = new ICM_Gauss1FLossDistrib_LJ(&c);
		pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );
		Distrib->Set_APPROX(&c);
	}

	try
	{
		ELT =Distrib->ComputeEL(&c);
	}
	catch (std::exception&e)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());
	}
	catch(...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib_smile :: Unable to compute Loss Distribution for maturity="<<yf);
	}
	

	if (ELT<0.) ELT = 0.;
	
	ELT /= (tranche_up-tranche_down);
//	ELT = 1. - ELT;

	return (ELT);
}