#include "ARMKernel/glob/firsttoinc.h" 
#include "icm_distriblosscalculator.h" 
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib.h"

// -------------------------------------------------------------------------------------
// loss distribution for a Full Homogeneous collateral with base correlation computation 
// --------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_fullhomogeneous(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	//calcul de l'expected loss dans le cas Term Structure Review
	if ((pricer.GetTermStructurePricing()==qTermStructureR) || (isStepUp))
	{return cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(yf,pricer_,model_,sec_,losses);}
	else if (pricer.GetTermStructurePricing()==qTermStructureR)
	{return cpt_elt_pricer_distrib_smile_fullhomogeneous_TSR(yf,pricer_,model_,sec_,losses);}

	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_LightDistrib* Distrib = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	int nbnames = collat->GetNbIssuers();
	double ELT = 0.,Pdefault=0.,beta_down=0.,beta_up=0.;

	Pdefault= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yf);

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down = correlation->GetBeta(collat->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up = correlation->GetBeta(collat->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_NONE);

	CtxtDistrib c;
	c.m_DistribType=qDISTRIB_STD;
	c.m_nbnames=nbnames;
	c.m_LossRates=lossUnits.getLossRates(date); // pricer.getLossRate();
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_LossUnit=lossUnits.getLossUnit(date);
	c.m_TrancheDown=tranche_down;
	c.m_TrancheUp=tranche_up;
	c.m_beta_dw_maturity=std::vector<double>(nbnames,beta_down);
	c.m_beta_up_maturity=std::vector<double>(nbnames,beta_up);
	// double* VPdefault = new double[1];VPdefault[0]=Pdefault;
	c.m_pdef_maturity.Resize(1,Pdefault); 
	c.m_IntegrationStep=pricer.GetIntegrationStep1();

	if (pricer.GetDistribution()== 0)
	{
		Distrib = new ICM_LightDistrib(nbnames);
		pricer.SetDistribution( Distrib);
	}
	else
		Distrib = dynamic_cast<ICM_LightDistrib*>(pricer.GetDistribution());

	try
	{
		ELT =Distrib->ComputeEL_FullHomog(&c);
		losses=c.m_OutLosses;
	}
	catch (std::exception&e)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case : "<<e.what());
	}
	catch (...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");
	}

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}

// -------------------------------------------------------------------------------------
// collat forward : loss distribution for a Full Homogeneous collateral with base correlation computation 
// --------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd(const double& yf, const double& yf_fwd_start,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	//calcul de l'expected loss dans le cas Term Structure Review
	if (pricer.GetTermStructurePricing()==qTermStructureR)
	{return cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd_TSR(yf,yf_fwd_start,pricer_,model_,sec_,losses);}

	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_LightDistrib* Distrib = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	int nbnames = collat->GetNbIssuers();
	double ELT = 0.,Pdefault=0.,beta_down=0.,beta_up=0.,Pdefault_start=0.;

	Pdefault= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yf);
	Pdefault_start= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yf_fwd_start);

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down = correlation->GetBeta(collat->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up = correlation->GetBeta(collat->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_NONE);

	CtxtDistrib c;
	c.m_DistribType=qDISTRIB_COLLFWD;
	c.m_nbnames=nbnames;
	c.m_LossUnit=lossUnits.getLossUnit(date);
	c.m_LossRates=lossUnits.getLossRates(date);
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_TrancheDown=tranche_down;
	c.m_TrancheUp=tranche_up;
	c.m_IntegrationStep=pricer.GetIntegrationStep1();

	c.m_beta_dw_maturity=std::vector<double>(nbnames,beta_down);
	c.m_beta_up_maturity=std::vector<double>(nbnames,beta_up);
	// double* VPdefault = new double[1];VPdefault[0]=Pdefault;
	c.m_pdef_maturity.Resize(1,Pdefault) ; // VPdefault;
	// double* VPdefault_start = new double[1];VPdefault_start[0]=Pdefault_start;
	c.m_collat_pdef_start.Resize(1,Pdefault_start);
	c.m_collat_beta_dw_start=c.m_beta_dw_maturity;
	c.m_collat_beta_up_start=c.m_beta_up_maturity;

	if (pricer.GetDistribution()== 0)
	{
		Distrib = new ICM_LightDistrib(nbnames);
		pricer.SetDistribution( Distrib);
	}
	else
		Distrib = dynamic_cast<ICM_LightDistrib*>(pricer.GetDistribution());

	try
	{
		ELT =Distrib->ComputeEL_FullHomog(&c);
		losses = c.m_OutLosses;

	}
	catch (std::exception&e)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case : "<<e.what());
	}
	catch (...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");
	}

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}

// -------------------------------------------------------------------------------------
// collat forward TSR: loss distribution for a Full Homogeneous collateral with base correlation computation 
// --------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd_TSR(const double& yf, const double& yf_fwd_start,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_LightDistrib* Distrib = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	int nbnames = collat->GetNbIssuers();
	double ELT = 0.;

	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;

	double start_beta_down_t1=0.,start_beta_up_t1=0.,start_beta_down_t2=0.,start_beta_up_t2=0.;
	double start_yt1_corr=0.,start_yt2_corr=0.;
	double start_yt1_pdef=0.,start_yt2_pdef=0.;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();

	DeduceYearTermsForTSR(correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);
	DeduceYearTermsForTSR(correlation,yf_fwd_start,start_yt1_corr,start_yt2_corr,start_yt1_pdef,start_yt2_pdef);

	/*
	double Pdefault_t1= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yt1_pdef);
	double Pdefault_t2= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yt2_pdef);
	double start_Pdefault_t1= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(start_yt1_pdef);
	double start_Pdefault_t2= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(start_yt2_pdef);
	*/
	double Pdefault_t1=0.,Pdefault_t2=0.,start_Pdefault_t1=0.,start_Pdefault_t2=0.;
	double Pdefault_down_t1=0.,Pdefault_down_t2=0.,start_Pdefault_down_t1=0.,start_Pdefault_down_t2=0.;

	Pdefault_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_UP);
	Pdefault_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_UP);
	Pdefault_down_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_LOW);
	Pdefault_down_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_LOW);

	start_Pdefault_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),start_yt1_corr,start_yt2_corr,start_yt1_pdef,qStrike_UP);
	start_Pdefault_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),start_yt1_corr,start_yt2_corr,start_yt2_pdef,qStrike_UP);
	start_Pdefault_down_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),start_yt1_corr,start_yt2_corr,start_yt1_pdef,qStrike_LOW);
	start_Pdefault_down_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),start_yt1_corr,start_yt2_corr,start_yt2_pdef,qStrike_LOW);

	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_LOW);

	start_beta_down_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),start_yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	start_beta_up_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),start_yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_LOW);
	start_beta_down_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),start_yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	start_beta_up_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),start_yt2_corr,CREDIT_DEFAULT_VALUE,yf);

	correlation->SetForcedStrikeType(qStrike_NONE);

	if (pricer.GetDistribution()== 0)
	{
		Distrib = new ICM_LightDistrib(nbnames);
		pricer.SetDistribution( Distrib);
	}
	else
		Distrib = dynamic_cast<ICM_LightDistrib*>(pricer.GetDistribution());

	CtxtDistrib c;
	c.m_DistribType=qDISTRIB_COLLFWD_TSR;
	c.m_nbnames=nbnames;
	c.m_LossUnit=lossUnits.getLossUnit(date);
c.m_LossRates=lossUnits.getLossRates(date);	
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_TrancheDown=tranche_down;
	c.m_TrancheUp=tranche_up;
	c.m_IntegrationStep=pricer.GetIntegrationStep1();

	c.m_ts_beta_dw_maturity_t1=std::vector<double>(nbnames,beta_down_t1);
	c.m_ts_beta_dw_maturity_t2=std::vector<double>(nbnames,beta_down_t2);
	c.m_ts_beta_up_maturity_t1=std::vector<double>(nbnames,beta_up_t1);
	c.m_ts_beta_up_maturity_t2=std::vector<double>(nbnames,beta_up_t2);
	// double* VPdefault_t1 = new double[1];VPdefault_t1[0]=Pdefault_t1;
	c.m_ts_pdef_up_maturity_t1.Resize(1,Pdefault_t1) ;
	// double* VPdefault_t2 = new double[1];VPdefault_t2[0]=Pdefault_t2;
	// c.m_ts_pdef_up_maturity_t2=VPdefault_t2;
	c.m_ts_pdef_up_maturity_t2.Resize(1,Pdefault_t2) ;
	// double* VPdefault_down_t1 = new double[1];VPdefault_down_t1[0]=Pdefault_down_t1;
	// c.m_ts_pdef_dw_maturity_t1=VPdefault_down_t1;
	c.m_ts_pdef_dw_maturity_t1.Resize(1,Pdefault_down_t1) ;
	// double* VPdefault_down_t2 = new double[1];VPdefault_down_t2[0]=Pdefault_down_t2;
	// c.m_ts_pdef_dw_maturity_t2=VPdefault_down_t2;
	c.m_ts_pdef_dw_maturity_t2.Resize(1,Pdefault_down_t2);
	c.m_ts_t1_corr=yt1_corr;
	c.m_ts_t2_corr=yt2_corr;

	c.m_ts_collat_beta_dw_start_t1=std::vector<double>(nbnames,start_beta_down_t1);
	c.m_ts_collat_beta_dw_start_t2=std::vector<double>(nbnames,start_beta_down_t2);
	c.m_ts_collat_beta_up_start_t1=std::vector<double>(nbnames,start_beta_up_t1);
	c.m_ts_collat_beta_up_start_t2=std::vector<double>(nbnames,start_beta_up_t2);
	//double* Vstart_Pdefault_t1 = new double[1];Vstart_Pdefault_t1[0]=start_Pdefault_t1;
	// c.m_ts_collat_pdef_up_start_t1=Vstart_Pdefault_t1;
	c.m_ts_collat_pdef_up_start_t1.Resize(1,start_Pdefault_t1); 
	//double* Vstart_Pdefault_t2 = new double[1];Vstart_Pdefault_t2[0]=start_Pdefault_t2;
	// c.m_ts_collat_pdef_up_start_t2=Vstart_Pdefault_t2;
	c.m_ts_collat_pdef_up_start_t2.Resize(1,start_Pdefault_t2);
	//double* Vstart_Pdefault_down_t1 = new double[1];Vstart_Pdefault_down_t1[0]=start_Pdefault_down_t1;
	// c.m_ts_collat_pdef_dw_start_t1=Vstart_Pdefault_down_t1;
	c.m_ts_collat_pdef_dw_start_t1.Resize(1,start_Pdefault_down_t1);
	//double* Vstart_Pdefault_down_t2 = new double[1];Vstart_Pdefault_down_t2[0]=start_Pdefault_down_t2;
	// c.m_ts_collat_pdef_dw_start_t2=Vstart_Pdefault_down_t2;
	c.m_ts_collat_pdef_dw_start_t2.Resize(1,start_Pdefault_down_t2);
	c.m_ts_collat_start_t1_corr=start_yt1_corr;
	c.m_ts_collat_start_t2_corr=start_yt2_corr;

	try
	{
	/*	ELT =Distrib->ComputeEL_FullHomog_collat_TSR(
				pricer.GetLossUnit(),tranche_down,tranche_up,
				beta_down_t1,beta_down_t2,beta_up_t1,beta_up_t2,
				Pdefault_t1,Pdefault_t2,Pdefault_down_t1,Pdefault_down_t2,
				yt1_corr,yt2_corr,start_beta_down_t1,start_beta_down_t2,
				start_beta_up_t1,start_beta_up_t2,
				start_Pdefault_t1,start_Pdefault_t2,start_Pdefault_down_t1,start_Pdefault_down_t2,
				start_yt1_corr,start_yt2_corr,pricer.GetIntegrationStep1(),losses);*/

		ELT =Distrib->ComputeEL_FullHomog(&c);
		losses=c.m_OutLosses;

	}
	catch (std::exception&e)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case : "<<e.what());
	}
	catch (...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");
	}

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}


// -------------------------------------------------------------------------------------
// loss distribution for a Full Homogeneous collateral with base correlation computation TSR 
// --------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_fullhomogeneous_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Mez& ftd=dynamic_cast<ICM_Mez&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	// cas stepup
	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	//calcul de l'expected loss dans le cas Term Structure Review
	if (isStepUp)
	{return cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(yf,pricer_,model_,sec_,losses);}


	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_LightDistrib* Distrib = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	int nbnames = collat->GetNbIssuers();
	double ELT = 0.,Pdefault_t1=0.,Pdefault_t2=0.;
	double Pdefault_down_t1=0.,Pdefault_down_t2=0.;
	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	DeduceYearTermsForTSR(correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);

	/*
	Pdefault_t1= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yt1_pdef);
	Pdefault_t2= model.GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yt2_pdef);*/

	Pdefault_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_UP);
	Pdefault_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_UP);
	Pdefault_down_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_LOW);
	Pdefault_down_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_LOW);

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_NONE);

	CtxtDistrib c;
	c.m_DistribType=qDISTRIB_STD_TSR;
	c.m_nbnames=nbnames;
	c.m_LossUnit=lossUnits.getLossUnit(date);
	c.m_LossRates=lossUnits.getLossRates(date);
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_TrancheDown=tranche_down;
	c.m_TrancheUp=tranche_up;
	c.m_IntegrationStep=pricer.GetIntegrationStep1();
	c.m_ts_beta_dw_maturity_t1=std::vector<double>(1,beta_down_t1);
	c.m_ts_beta_dw_maturity_t2=std::vector<double>(1,beta_down_t2);
	c.m_ts_beta_up_maturity_t1=std::vector<double>(1,beta_up_t1);
	c.m_ts_beta_up_maturity_t2=std::vector<double>(1,beta_up_t2);
	//double* VPdefault_t1 = new double[1];VPdefault_t1[0]=Pdefault_t1;
	// c.m_ts_pdef_up_maturity_t1=VPdefault_t1;
	c.m_ts_pdef_up_maturity_t1.Resize(1,Pdefault_t1);
	//double* VPdefault_t2 = new double[1];VPdefault_t2[0]=Pdefault_t2;
	// c.m_ts_pdef_up_maturity_t2=VPdefault_t2;
	c.m_ts_pdef_up_maturity_t2.Resize(1,Pdefault_t2); 
	//double* VPdefault_down_t1 = new double[1];VPdefault_down_t1[0]=Pdefault_down_t1;
	// c.m_ts_pdef_dw_maturity_t1=VPdefault_down_t1;
	c.m_ts_pdef_dw_maturity_t1.Resize(1,Pdefault_down_t1); 
	//double* VPdefault_down_t2 = new double[1];VPdefault_down_t2[0]=Pdefault_down_t2;
	// c.m_ts_pdef_dw_maturity_t2=VPdefault_down_t2;
	c.m_ts_pdef_dw_maturity_t2.Resize(1,Pdefault_down_t2); 
	c.m_ts_t1_corr=yt1_corr;
	c.m_ts_t2_corr=yt2_corr;

	if (pricer.GetDistribution()== 0)
	{
		Distrib = new ICM_LightDistrib(nbnames);
		pricer.SetDistribution( Distrib);
	}
	else
		Distrib = dynamic_cast<ICM_LightDistrib*>(pricer.GetDistribution());

	try
	{
		ELT =Distrib->ComputeEL_FullHomog(&c);
		losses = c.m_OutLosses;

	}
	catch (std::exception&e)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case : "<<e.what());
	}
	catch (...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");
	}

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}


// -------------------------------------------------------------------------------------
// loss distribution for a Full Homogeneous collateral with base correlation computation TSR 
// --------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_smile_fullhomogeneous_stepup_TSR(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;
	int i=0;

	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Mez& ftd=dynamic_cast<ICM_Mez&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_LightDistrib* Distrib = NULL;

	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	vector<double> tranche_down;tranche_down.resize(stepup_dates.size());
	vector<double> tranche_up;tranche_up.resize(stepup_dates.size());

	for (i=0;i<stepup_dates.size();i++)
	{
		tranche_down[i]= ABS(round(pricer.GetTranche_Down(stepup_dates[i])));
		tranche_up[i]= ABS(round(pricer.GetTranche_Up(stepup_dates[i])));
	}
	
	int nbnames = collat->GetNbIssuers();
	double ELT = 0.,Pdefault_t1=0.,Pdefault_t2=0.;
	double Pdefault_down_t1=0.,Pdefault_down_t2=0.;
	double beta_down_t1=0.,beta_up_t1=0.,beta_down_t2=0.,beta_up_t2=0.;
	double yt1_corr=0.,yt2_corr=0.;
	double yt1_pdef=0.,yt2_pdef=0.;

	double yt_stepup = 0.;
	double yt1_corr_stepup=0.,yt2_corr_stepup=0.;
	double yt1_pdef_stepup=0.,yt2_pdef_stepup=0.;
	double Pdefault_t1_stepup = 0.,Pdefault_t2_stepup=0.;
	double Pdefault_down_t1_stepup = 0.,Pdefault_down_t2_stepup = 0.;
	double beta_down_t1_stepup=0.,beta_up_t1_stepup=0.,beta_down_t2_stepup=0.,beta_up_t2_stepup=0.;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	DeduceYearTermsForTSR(correlation,yf,yt1_corr,yt2_corr,yt1_pdef,yt2_pdef);
	if (stepup_dates.size()>1)
	{	yt_stepup = (stepup_dates[stepup_dates.size()-2]-model.GetStartDate().GetJulian())/365.;
		DeduceYearTermsForTSR(correlation,yt_stepup,yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,yt2_pdef_stepup);	
	}

	if (CHECK_EQUAL(yf,yt_stepup))
	{isStepUp=false;}

	Pdefault_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_UP);
	Pdefault_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_UP);
	Pdefault_down_t1=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt1_pdef,qStrike_LOW);
	Pdefault_down_t2=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr,yt2_corr,yt2_pdef,qStrike_LOW);

	Pdefault_t1_stepup=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,qStrike_UP);
	Pdefault_t2_stepup=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr_stepup,yt2_corr_stepup,yt2_pdef_stepup,qStrike_UP);
	Pdefault_down_t1_stepup=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr_stepup,yt2_corr_stepup,yt1_pdef_stepup,qStrike_LOW);
	Pdefault_down_t2_stepup=ComputeBarrierTSR(correlation,model.GetDefaultCurve(collat->GetIssuersLabels(0)),yt1_corr_stepup,yt2_corr_stepup,yt2_pdef_stepup,qStrike_LOW);

	//correlation
	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t1 = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t2 = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr,CREDIT_DEFAULT_VALUE,yf);
	correlation->SetForcedStrikeType(qStrike_NONE);

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t1_stepup = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t1_stepup = correlation->GetBeta(collat->GetIssuersLabels(0),yt1_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down_t2_stepup = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up_t2_stepup = correlation->GetBeta(collat->GetIssuersLabels(0),yt2_corr_stepup,CREDIT_DEFAULT_VALUE,yt_stepup);
	correlation->SetForcedStrikeType(qStrike_NONE);

	CtxtDistrib c;
	c.m_nbnames = collat->GetNbIssuers();
	c.m_DistribType=qDISTRIB_STD_TSR;
	c.m_LossUnit=lossUnits.getLossUnit(Maturity);
	c.m_LossRates=lossUnits.getLossRates(Maturity);
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_TrancheDown=tranche_down[tranche_down.size()-3];
	c.m_TrancheUp=tranche_up[tranche_up.size()-3];
	c.m_TrancheDown_stepup=tranche_down[tranche_down.size()-1];
	c.m_TrancheUp_stepup=tranche_up[tranche_up.size()-1];
	c.m_IntegrationStep=pricer.GetIntegrationStep1();
	c.m_ts_t1_corr=yt1_corr;
	c.m_ts_t2_corr=yt2_corr;
//	c.m_pricer=&pricer;

	//term structure
	c.m_ts_beta_dw_maturity_t1=std::vector<double>(1,beta_down_t1);
	c.m_ts_beta_dw_maturity_t2=std::vector<double>(1,beta_down_t2);
	c.m_ts_beta_up_maturity_t1=std::vector<double>(1,beta_up_t1);
	c.m_ts_beta_up_maturity_t2=std::vector<double>(1,beta_up_t2);
	// double* VPdefault_t1 = new double[1];VPdefault_t1[0]=Pdefault_t1;
	// c.m_ts_pdef_up_maturity_t1=VPdefault_t1;
	c.m_ts_pdef_up_maturity_t1.Resize(1,Pdefault_t1); 
	// double* VPdefault_t2 = new double[1];VPdefault_t2[0]=Pdefault_t2;
	// c.m_ts_pdef_up_maturity_t2=VPdefault_t2;
	c.m_ts_pdef_up_maturity_t2.Resize(1,Pdefault_t2); 
	//double* VPdefault_down_t1 = new double[1];VPdefault_down_t1[0]=Pdefault_down_t1;
	// c.m_ts_pdef_dw_maturity_t1=VPdefault_down_t1;
	c.m_ts_pdef_dw_maturity_t1.Resize(1,Pdefault_down_t1); 
	// double* VPdefault_down_t2 = new double[1];VPdefault_down_t2[0]=Pdefault_down_t2;
	// c.m_ts_pdef_dw_maturity_t2=VPdefault_down_t2;
	c.m_ts_pdef_dw_maturity_t2.Resize(1,Pdefault_down_t2); 

	//step up subordination
	c.m_ts_stepup_beta_dw_maturity_t1=std::vector<double>(1,beta_down_t1_stepup);
	c.m_ts_stepup_beta_dw_maturity_t2=std::vector<double>(1,beta_down_t2_stepup);
	c.m_ts_stepup_beta_up_maturity_t1=std::vector<double>(1,beta_up_t1_stepup);
	c.m_ts_stepup_beta_up_maturity_t2=std::vector<double>(1,beta_up_t2_stepup);
	// double* VPdefault_t1_stepup = new double[1];VPdefault_t1_stepup[0]=Pdefault_t1_stepup;
	// c.m_ts_stepup_pdef_up_maturity_t1=VPdefault_t1_stepup;
	c.m_ts_stepup_pdef_up_maturity_t1.Resize(1,Pdefault_t1_stepup);
	// double* VPdefault_t2_stepup = new double[1];VPdefault_t2_stepup[0]=Pdefault_t2_stepup;
	// c.m_ts_stepup_pdef_up_maturity_t2=VPdefault_t2_stepup;
	c.m_ts_stepup_pdef_up_maturity_t2.Resize(1,Pdefault_t2_stepup); 
	// double* VPdefault_down_t1_stepup = new double[1];VPdefault_down_t1_stepup[0]=Pdefault_down_t1_stepup;
	// c.m_ts_stepup_pdef_dw_maturity_t1=VPdefault_down_t1_stepup;
	c.m_ts_stepup_pdef_dw_maturity_t1.Resize(1,Pdefault_down_t1_stepup); 
	// double* VPdefault_down_t2_stepup = new double[1];VPdefault_down_t2_stepup[0]=Pdefault_down_t2_stepup;
	// c.m_ts_stepup_pdef_dw_maturity_t2=VPdefault_down_t2_stepup;
	c.m_ts_stepup_pdef_dw_maturity_t2.Resize(1,Pdefault_down_t2_stepup); 

	if (pricer.GetDistribution()== 0)
	{
		Distrib = new ICM_LightDistrib(nbnames);
		pricer.SetDistribution( Distrib);
	}
	else
		Distrib = dynamic_cast<ICM_LightDistrib*>(pricer.GetDistribution());

	if (!isStepUp) //cas non stepup dans l'intervalle de temps [0;T0[
	{
		try
		{	ELT =Distrib->ComputeEL_FullHomog(&c);
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

	return (ELT);
}