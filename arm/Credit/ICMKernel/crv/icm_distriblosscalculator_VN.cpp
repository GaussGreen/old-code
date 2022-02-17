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
cpt_elt_pricer_distrib_smile_VN(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	if (yf<=0) return 0.;

	int k=0,j=0;

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

	ARM_Vector beta (nbnames); 

	double ELT = 0.;

	qIntegratorChoice	TheIntegratorType;
	if (! pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		pricer.SetIntegrationStep1(40) ;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());
	
	double	current_Pdefault;

	ARM_Vector Pdefault(nbnames); 
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	std::vector<int> ReverseCollatRank (collatRank.size()); 
	for (j=0;j<collatRank.size(); j++)  {ReverseCollatRank[collatRank[j]]=j;} 

	for (int i=0; i<nbnames; i++) 
	{	
		current_Pdefault	=	model.GetDefaultCurve(ftd.GetCollateral()->GetIssuersLabels(i))->DefaultProba(yf);
		
		if (CHECK_EQUAL(current_Pdefault, 0.0))
		{Pdefault[i] = 0.0;	}
		else if (CHECK_EQUAL(current_Pdefault, 1.0))
		{Pdefault[i] = 1.0;}
		else if (current_Pdefault > 1.)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<ftd.GetCollateral()->GetIssuersLabels(i));}
		else if (current_Pdefault < 0.0)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<ftd.GetCollateral()->GetIssuersLabels(i));}
		else Pdefault[i] = current_Pdefault;
	}

	CtxtDistrib c_;
	c_.m_IsUsed=true;
	c_.m_DistribType=qDISTRIB_VN;
	c_.m_nbnames=nbnames;
	c_.m_LossUnit=lossUnits.getLossUnit(AsOf);
	c_.m_TrancheDown=tranche_down;
	c_.m_LossRates.Resize(nbnames);
	
	const ARM_Vector& LR1=lossUnits.getLossRates(AsOf);
	for (j=0;j<nbnames;j++)
	{c_.m_LossRates[j] = LR1.Elt(ReverseCollatRank[j]);}
	
	c_.m_CopulaType=pricer.GetCopulaType();
	c_.m_IntegrationMethod=TheIntegratorType;
	c_.m_TrancheUp=tranche_up;
	c_.m_pdef_maturity=Pdefault;
	c_.m_IntegrationStep=pricer.GetIntegrationStep1();

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	std::vector<double> corr_low;corr_low.resize(nbnames);
	std::vector<double> corr_Hight;corr_Hight.resize(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low[i] = correlation->GetBeta(ftd.GetCollateral()->GetIssuersLabels(ReverseCollatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight[i] = correlation->GetBeta(ftd.GetCollateral()->GetIssuersLabels(ReverseCollatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	}

	c_.m_beta_dw_maturity=corr_low;
	c_.m_beta_up_maturity=corr_Hight;

	CtxtDistrib c=c_;
	c_.m_pdef_maturity.Resize(0) ; 

	// ----------------------------------------------------------------------------------------------
	// BEGIN : Determination du vecteur de dates pour la Loss Unit
	// ----------------------------------------------------------------------------------------------
	int szmax=lossUnits.size() ;
	int min_idx=0;

	vector<double> LUDates;
	ARM_Date FirstDate = (ARM_Date) model.GetStartDate();

	//if (FirstDate>=lossUnits.getDate(idxl)) 
	if (FirstDate>=lossUnits.getDate(0)) 
	{LUDates.push_back(FirstDate.GetJulian());}

	//detection des dates < Asof
	for (k=0;k<szmax;k++)
	{	
		ARM_Date date_ = (ARM_Date) lossUnits.getDate(k);
		double year = (date_.GetJulian()-model.GetStartDate().GetJulian())/K_YEAR_LEN; 

		if ((date_>model.GetStartDate()) && !CHECK_EQUAL(date_.GetJulian(),model.GetStartDate().GetJulian()) &&	(year<yf) && (!lossUnits.IsLossUnitNull(date_)))
			{LUDates.push_back(date_.GetJulian());}
	}

	if (LUDates.size()) 
	{
		ARM_Date LDate =(ARM_Date)LUDates[LUDates.size()-1];
		ARM_Date LastDate = (ARM_Date) model.GetStartDate().GetJulian() + K_YEAR_LEN*yf;

		if ((LastDate>=lossUnits.getDate(0)) && (LastDate!=LDate))
			{LUDates.push_back(LastDate.GetJulian());}
	}

	// ----------------------------------------------------------------------------------------------
	// END : Determination du vecteur de dates pour la Loss Unit
	// ----------------------------------------------------------------------------------------------

	c.m_vn_contexts.resize(LUDates.size());

	if (!LUDates.size()) 
		{return 0;}

	for (k=0;k<LUDates.size();k++)
	{
		CtxtDistrib ctemp=c_;
		
		ARM_Date date_ = (ARM_Date) LUDates[k];
		ctemp.m_YearTerm= (date_.GetJulian()-model.GetStartDate().GetJulian())/K_YEAR_LEN; 
		date_ -= 1.; //on retranche d'un jour

		const std::vector<int> &collatRank_ = lossUnits.getCollatRank(date_); 
		ReverseCollatRank.resize(collatRank_.size()); 
		for (j=0;j<collatRank_.size(); j++)  {ReverseCollatRank[collatRank_[j]]=j;} 


		ARM_Vector Pdefault (nbnames);
		for (int i=0; i<nbnames; i++) 
		{	
		current_Pdefault =	model.GetDefaultCurve(ftd.GetCollateral()->GetIssuersLabels(i))->DefaultProba(ctemp.m_YearTerm);
		
		if (CHECK_EQUAL(current_Pdefault, 0.0))
		{Pdefault[i] = 0.0;	}
		else if (CHECK_EQUAL(current_Pdefault, 1.0))
		{Pdefault[i] = 1.0;}
		else if (current_Pdefault > 1.)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<ftd.GetCollateral()->GetIssuersLabels(i));}
		else if (current_Pdefault < 0.0)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<ftd.GetCollateral()->GetIssuersLabels(i));}
		else Pdefault[i] = current_Pdefault;
		}

		ctemp.m_pdef_maturity = Pdefault;
		ctemp.m_LossUnit=lossUnits.getLossUnit(date_);
		const ARM_Vector& LR=lossUnits.getLossRates(date_);

		for (j=0;j<nbnames;j++)
		{ctemp.m_LossRates[j] = LR.Elt(ReverseCollatRank[j]);}

		if (k==0) {ctemp.m_IsBegin=true;}

		c.m_vn_contexts[k]=ctemp;

/*		if (k)
		{
			if (!(c.m_vn_contexts[k-1]>c.m_vn_contexts[k]))
				{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: invalid ordination for context yt:"<<ctemp.m_YearTerm);}
		}
*/	}

	correlation->SetForcedStrikeType(qStrike_NONE);

	if (pricer.GetDistribution() == 0)
	{
		Distrib = new ICM_Gauss1FLossDistrib_LJ(&c);
		pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );
		Distrib->Set_VN(&c);
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

	return (ELT);
}