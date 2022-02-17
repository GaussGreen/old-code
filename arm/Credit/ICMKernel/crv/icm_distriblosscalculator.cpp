#include "ARMKernel/glob/firsttoinc.h" 
#include "icm_distriblosscalculator.h" 
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\glob\icm_correlation_sector.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_2F.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_MF.h"
#include "ICMKernel\inst\icm_mez.h"
// #include "ICMKernel\pricer\icm_pricer_cds.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\util\icm_RootFinderND.h"
#include "ICMKernel\util\icm_brentsolver.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
//#include "ICMKernel\mod\icm_Homogene_FastLossDistrib.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"

//	---------------------------------------------------------------------------------
//	helper to compute the "SortedInd". 
std::vector<int> 
CptSortedInd(const ICM_Collateral& collat,const ICM_LossUnits&lossUnits,const ARM_Date&date,const ICM_Correlation&corr) 
{
	const ARM_Vector& LSLossRates = lossUnits.getLSLossRates(date);
	const std::vector<int> & collatRank = lossUnits.getCollatRank(date);
	int nbNames = LSLossRates.size(); 
	std::vector<int> sorted_ind_temp(nbNames); 	
	int nbShort=0 ; 
	//	counting how many short names.
	//	this does not depend on collat order 
	for(int i=0;i<nbNames;i++) if (LSLossRates.Elt(i)<0) nbShort++; 
	int iShort=0,iLong=0; 
	//	aggregating positions 
	//	sorted_ind referenced the ordered collat
	for(i=0;i<nbNames;i++) 
	{
		if (LSLossRates.Elt(i)<0) 
			sorted_ind_temp[iShort++]= i; // collatRank[i] ; 
		else sorted_ind_temp[nbShort+iLong++]=i; // collatRank[i] ; 
	}
	if (corr.GetName() != ICM_CORRELATION_SECTOR ) return sorted_ind_temp ;
	//
	const ICM_Correlation_Sector&	correlation	=	dynamic_cast<const ICM_Correlation_Sector&>(corr);
	std::vector<int> Ref_Sector_Membership ; 
	int NbSectors; 
	correlation.Get_Sector_Membership(Ref_Sector_Membership);
	correlation.GetNbSectors(NbSectors);
	std::vector<int> sorted_ind(nbNames); 
	int FinalRank=0; 
	for (int iSector = 0; iSector < NbSectors; iSector++)
	{
		for (int iTotal = 0; iTotal < nbNames; iTotal ++)
		{
			int InitialRank;
			int Rank;;
			InitialRank = sorted_ind_temp[iTotal];
			const std::string& Label = collat.GetIssuersLabels(InitialRank); 
			Rank =  correlation.Get_IssuersRankFromLabel(Label);
			if (Ref_Sector_Membership[Rank] == iSector)
			{	
				sorted_ind [FinalRank]=InitialRank;
				FinalRank++;
			}
		}
	}
	return sorted_ind; 
}


// ----------------------------------------------
// loss distribution for single name product
// ----------------------------------------------
double 
cpt_elt_default(const double& yf, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	
	ICM_ModelMultiCurves& model=dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	const ICM_DefaultCurve* defaultcurve = model.GetDefaultCurves(0); 
	if (!defaultcurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_default: Can't get Default Curve 0"); 
	double ELT = defaultcurve->DefaultProba(yf);
	return (ELT);
	
}

// ---------------------------------------------------------------------------------------------
// loss distribution for collateral with beta correlation (use only with ntd or ftd securities)
// ---------------------------------------------------------------------------------------------
double 
cpt_elt_pricer_distrib(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	ICM_Ftd & ftd= dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Pricer_Distrib& pricer = dynamic_cast<ICM_Pricer_Distrib&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Collateral* collat=ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits=  pricer.getLossUnits(); 

	double ELT = 0.;
	if (yf<=0) return 0.;

// 17783 	if (pricer.IsActivateShift()) 
// 17783 		return cpt_elt_pricer_distrib_fast(yf,pricer_,model_,sec_); 

	int nbnames = collat->GetNbIssuers();
		
	ARM_Vector* Betas = pricer.GetBetas();
	// double* beta = new double[nbnames];
	// memcpy(beta,Betas->GetElt(),sizeof(double)*nbnames);
	// double* Pdefault = new double[nbnames];
	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	ARM_Vector Pdefault (nbnames); 
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<nbnames; i++)
	{
		Pdefault[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);
		if (Pdefault[i]>1.)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
		else if (Pdefault[i]<0.0)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
	}
	
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	//Prise en compte des portefeuilles LongShort dans un modèle avec Beta différents
	//Très sale : a modifier
	if (ftd.GetCollateral()->IsLongShort())
	{
		//GaussHermite et 40 obligatoire si Long Short
		qIntegratorChoice	TheIntegratorType;
		TheIntegratorType =qGAUSS_HERMITE;
		
		ICM_Gauss1FLossDistrib_LJ Distrib(nbnames,Pdefault,*Betas,
					// pricer.getLongShort_LossRate(), 
					lossUnits.getLSLossRates(date),
					// *pricer.GetSortedInd(),
					CptSortedInd(*collat,lossUnits,date,*model.GetCorrelation()), 
					40,
					qGAUSSIAN,
					TheIntegratorType,
					40);
		
		// if (Pdefault) delete[] Pdefault;

		ARM_Vector corr_low(nbnames);
		ARM_Vector corr_Hight(nbnames);

		//Reconstruction du vecteur des betas
		for (i=0; i<nbnames; i++) 
		{ 
			corr_low[i] = Betas->Elt(i);
			corr_Hight[i] = Betas->Elt(i);
		}

		Distrib.UpdateCorrelation(corr_low,corr_Hight);

		// New...
		Distrib.SetCopulaType(qGAUSSIAN);
		Distrib.SetIntegratorType(TheIntegratorType,40);
		Distrib.SetIntegrationStep(40);
		
		try
		{
			//Insertion du vecteur d'indice triés
			ELT =Distrib.compute_expectedlosstranche_LongShort(tranche_up,
														tranche_down,
														lossUnits.getLossUnit(date)
														); 
			
		}
		catch (std::exception&e)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());
		}
		catch(...)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib_smile :: Unable to compute Loss Distribution for maturity="<<yf);
		}


	}
	else
	{
		ICM_Gauss1FLossDistrib Distrib(nbnames,Pdefault,*Betas,lossUnits.getLossRates(date));
		// if (Pdefault) delete[] Pdefault;

		Distrib.setPercentIndep(model.GetIndependantPart());
		Distrib.setPercentFullCorrel(model.GetFullCorrelPart());

		try 
		{
			ELT =	Distrib.compute_expectedlosstranche(tranche_up,
															 tranche_down,
															 lossUnits.getLossUnit(date),
// 17783 																	 pricer.GetTenorShift(),
// 17783 																	 pricer.GetIssuerShift(),
															 losses);
		}
		catch (std::exception&e)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());
		}
		catch(...)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib:: Unable to compute Loss Distribution for maturity="<<yf);
		}
	}

	if (ELT<0) ELT =0.;
	// if (beta) delete[] beta;

	ELT /= (tranche_up-tranche_down);

	return (ELT);
}

// ------------------------------------------------------------------
// loss distribution for collateral with base correlation computation
// ------------------------------------------------------------------

//qDISTRIB_TYPE TypeDistrib = qDISTRIB_STD;

double 
cpt_elt_pricer_distrib_smile(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd & ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral *collat=ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits=  pricer.getLossUnits(); 

	vector<double> stepup_dates;
	bool isStepUp= ftd.SearchBoundsForStepUp(yf,(ARM_Date)model.GetStartDate(),stepup_dates);

	//Fast Approximation
	if (pricer.GetApproxComputation())
	{return cpt_elt_pricer_distrib_smile_APPROX(yf,pricer_,model_,sec_,losses);}

	//cas loss unit variable following time
	if (lossUnits.size()>1)
	{return cpt_elt_pricer_distrib_smile_VN(yf,pricer_,model_,sec_,losses);}

	//calcul de l'expected loss dans le cas Term Structure Review
	if ((pricer.GetTermStructurePricing()==qTermStructureR) && (isStepUp))
	{return cpt_elt_pricer_distrib_smile_stepup_TSR(yf,pricer_,model_,sec_,losses);}

	if (pricer.GetTermStructurePricing()==qTermStructureR)
	{return cpt_elt_pricer_distrib_smile_TSR(yf,pricer_,model_,sec_,losses);}

	// -----------------------------------------------------------
	// Sector Correlation for CDO2 if needed
	ICM_Correlation* tmp_correlation = model.GetCorrelation();

	if (tmp_correlation->GetName() == ICM_CORRELATION_SECTOR)
		//return	cpt_elt_pricer_distrib_sector(yf, pricer_, model_, sec_);
		return	cpt_elt_pricer_distrib_MF(yf, pricer_, model_, sec_);
	// -----------------------------------------------------------

	//case full homogeneous 
	if (collat->IsFullHomog()) 
		return cpt_elt_pricer_distrib_smile_fullhomogeneous(yf,pricer_,model_,sec_,losses) ; 
		// return pricer.ExpectedLossTrancheFullHomogeneous(yf);
	if (yf<=0) return 0.;

	//for sensivities
// 17783 	if (pricer.IsActivateShift())  
// 17783 		return cpt_elt_pricer_distrib_smile_fast(yf,pricer_,model_,sec_) ; 

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
	
	double	current_Pdefault;

	// double* Pdefault = new double[nbnames];
	ARM_Vector Pdefault(nbnames); 
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<nbnames; i++) 
	{	

		current_Pdefault	=	model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);
		
		if (CHECK_EQUAL(current_Pdefault, 0.0))
		{Pdefault[i] = 0.0;	}
		else if (CHECK_EQUAL(current_Pdefault, 1.0))
		{Pdefault[i] = 1.0;}
		else if (current_Pdefault > 1.)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<collat->GetIssuersLabels(collatRank[i]));}
		else if (current_Pdefault < 0.0)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<collat->GetIssuersLabels(collatRank[i]));}
		else Pdefault[i] = current_Pdefault;
	}

	CtxtDistrib c;
	c.m_IsUsed=true;
	c.m_DistribType=qDISTRIB_STD;
	//c.m_DistribType=TypeDistrib;
	c.m_YearTerm = yf;
	c.m_nbnames=nbnames;
	c.m_LossUnit=lossUnits.getLossUnit(date);
	c.m_TrancheDown=tranche_down;
	c.m_LossRates=lossUnits.getLossRates(date);
	c.m_CopulaType=pricer.GetCopulaType();
	c.m_IntegrationMethod=TheIntegratorType;
	c.m_TrancheUp=tranche_up;
	c.m_pdef_maturity=Pdefault;
	c.m_IntegrationStep=pricer.GetIntegrationStep1();
	c.m_beta_dw_maturity=beta ;
	c.m_beta_up_maturity=beta; 

	if (pricer.GetDistribution() == 0)
	{
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Pour l'instant GaussHermite obligatoire si Long Short
			TheIntegratorType =qGAUSS_HERMITE;
			// integrationStep1= 40;
			pricer.SetIntegrationStep1(40); 
			Distrib = new ICM_Gauss1FLossDistrib_LJ(nbnames,Pdefault,beta,
				lossUnits.getLSLossRates(date), // pricer.getLongShort_LossRate(), 
				CptSortedInd(*collat,lossUnits,date,*tmp_correlation),
				// *pricer.GetSortedInd(),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1()
				);
		}
		else
		/*	Distrib = new ICM_Gauss1FLossDistrib_LJ(nbnames,Pdefault,beta,
				pricer.getLossRate(),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1()); */
			
			Distrib = new ICM_Gauss1FLossDistrib_LJ(&c);

	pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gauss1FLossDistrib_LJ*> ( pricer.GetDistribution() );
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Pour l'instant GaussHermite obligatoire si Long Short
			TheIntegratorType =qGAUSS_HERMITE;
			// integrationStep1= 40;
			pricer.SetIntegrationStep1(40) ;
			Distrib->Set(nbnames,Pdefault,beta,
				lossUnits.getLSLossRates(date), // pricer.getLongShort_LossRate(), 
				// *pricer.GetSortedInd() , 
				CptSortedInd(*collat,lossUnits,date,*tmp_correlation),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1());
		}
		else
			/*
			Distrib->Set(nbnames,Pdefault,beta,
				pricer.getLossRate(),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1()); */
			{Distrib->Set(&c);}
	}

//	if (Pdefault) delete[] Pdefault;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector corr_low(nbnames);
	ARM_Vector corr_Hight(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	}

	c.m_beta_dw_maturity=corr_low;
	c.m_beta_up_maturity=corr_Hight;

	correlation->SetForcedStrikeType(qStrike_NONE);
	Distrib->UpdateCorrelation(corr_low,corr_Hight);

	// New...
	Distrib->SetCopulaType(qGAUSSIAN);
	Distrib->SetIntegratorType(TheIntegratorType,pricer.GetIntegrationStep1());
	Distrib->SetIntegrationStep(pricer.GetIntegrationStep1());


	try
	{
		// do not deal with two vectors of betas 
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Insertion du vecteur d'indice triés
			ELT =Distrib->compute_expectedlosstranche_LongShort(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date)
													); 
		}
		else
		{
			
			ELT =Distrib->ComputeEL(&c);

		}
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
//	if (beta) delete[] beta;
	
	ELT /= (tranche_up-tranche_down);

	return (ELT);
}

// ------------------------------------------------------------------
// loss distribution for collateral with Multi factor 
// ------------------------------------------------------------------
double 
cpt_elt_pricer_distrib_MF(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_)
{
	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits=  pricer.getLossUnits(); 

	ICM_Correlation* tmp_correlation = model.GetCorrelation();
	ICM_Correlation_Sector*	correlation	=	(ICM_Correlation_Sector*) tmp_correlation;

	qTWO_FACTORS_CORRELATION_TYPE	CorrelationType;
	correlation->GetCorrelationType(CorrelationType);

	int	Nb_Sectors;
	correlation->GetNbSectors(Nb_Sectors);

	vector<int>	Ref_Sector_Membership;
	correlation->Get_Sector_Membership(Ref_Sector_Membership);

	vector<double>	Ref_Betas;
	correlation->Get_Betas(Ref_Betas);

	vector<double>	Ref_Lambdas;
	correlation->Get_Lambdas(Ref_Lambdas);

	vector<double>	Ref_Betas_Down;
	correlation->Get_Betas_Down(Ref_Betas_Down);

	vector<double>	Ref_Lambdas_Down;
	correlation->Get_Lambdas_Down(Ref_Lambdas_Down);
	

	if (yf<=0) return 0.;

	//for sensivities
	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;
	ICM_Gaussian_LossDistrib_MF* Distrib = NULL;

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
	vector<double>	VBetas(nbnames,0.);
	vector<double>	VBetas_Down(nbnames,0.);

	int	inam;
	int isec;
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (inam=0; inam<nbnames; inam++)
	{
		int RankCorrel;
		RankCorrel = correlation->Get_IssuersRankFromLabel(ftd.GetCollateral()->GetIssuersLabels(collatRank[inam]));
		isec=Ref_Sector_Membership[RankCorrel];
		VBetas[collatRank[inam]]	=	Ref_Betas[isec];
		VBetas_Down[collatRank[inam]]	=	Ref_Betas_Down[isec];
	}

	double ELT = 0.;

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
	ARM_Vector Pdefault (nbnames); 
	for (int i=0; i<nbnames; i++) 
	{	
		current_Pdefault	=	model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);
		
		if (CHECK_EQUAL(current_Pdefault, 0.0))
		{
			Pdefault[i] = 0.0;
		}
		else if (CHECK_EQUAL(current_Pdefault, 1.0))
		{
			Pdefault[i] = 1.0;
		}
		else if (current_Pdefault > 1.)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
		else if (current_Pdefault < 0.0)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
		else
			Pdefault[i] = current_Pdefault;
	}

	if (pricer.GetDistribution() == 0)
	{
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Pour l'instant GaussHermite obligatoire si Long Short
			TheIntegratorType =qGAUSS_HERMITE;
			// integrationStep1= 40;
			//pricer.SetIntegrationStep1(40); 

			Distrib=new ICM_Gaussian_LossDistrib_MF(nbnames,
									Pdefault,
									beta,
									// pricer.getLongShort_LossRate(),
									lossUnits.getLSLossRates(date),
									Nb_Sectors,
									Ref_Sector_Membership,
									// *pricer.GetSortedInd(),
									CptSortedInd(*collat,lossUnits,date,*correlation),
									qNO_COPULA,
									qGAUSS_HERMITE,
									pricer.GetIntegrationStep1()); 

		}
		else
			Distrib=new ICM_Gaussian_LossDistrib_MF(nbnames,
												Pdefault,
												beta,
												// pricer.getLossRate(),
												lossUnits.getLossRates(date),
												Nb_Sectors,
												Ref_Sector_Membership,
												// *pricer.GetSortedInd(),
												CptSortedInd(*collat,lossUnits,date,*correlation),
												qNO_COPULA,
												qGAUSS_HERMITE,
												pricer.GetIntegrationStep1()); 

		pricer.SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = dynamic_cast<ICM_Gaussian_LossDistrib_MF*> ( pricer.GetDistribution() );
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Pour l'instant GaussHermite obligatoire si Long Short
			TheIntegratorType =qGAUSS_HERMITE;
			Distrib->Set(nbnames,Pdefault,beta,
				// pricer.getLongShort_LossRate(),
				lossUnits.getLSLossRates(date),
				Nb_Sectors,
				Ref_Sector_Membership,
				CptSortedInd(*collat,lossUnits,date,*correlation),
				// *pricer.GetSortedInd(),
				qNO_COPULA,	qGAUSS_HERMITE,	pricer.GetIntegrationStep1());
		}
		else
			Distrib->Set(nbnames,Pdefault,beta,
				// pricer.getLossRate(),
				lossUnits.getLSLossRates(date),
				Nb_Sectors,
				Ref_Sector_Membership,
				// *pricer.GetSortedInd(),
				CptSortedInd(*collat,lossUnits,date,*correlation),
				qNO_COPULA,	qGAUSS_HERMITE,	pricer.GetIntegrationStep1());
	}

	//if (Pdefault) delete[] Pdefault;

	Distrib->UpdateCorrelation(VBetas_Down,VBetas,Ref_Lambdas_Down,Ref_Lambdas);

	try
	{
		// do not deal with two vectors of betas 
		if (ftd.GetCollateral()->IsLongShort())
		{
			//Insertion du vecteur d'indice triés
			ELT =Distrib->compute_expectedlosstranche_LongShort(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date)
													); 
		}
		else
		{
			ELT =Distrib->compute_expectedlosstranche(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date)
// 17783 															pricer.GetTenorShift(),
// 17783 															pricer.GetIssuerShift()
													); 

		}
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
	// if (beta) delete[] beta;

	ELT /= (tranche_up-tranche_down);

	return (ELT);
}



// ---------------------------------------------------------------
// Loss distribution with base correlation computation & Fwd Start Collat
// ---------------------------------------------------------------
double cpt_elt_pricer_distrib_smile_collat_fwd(const double& yf, const double& yf_fwd_start, ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_,vector<double>& losses)
{
	ICM_Pricer_Distrib_Smile& pricer=dynamic_cast<ICM_Pricer_Distrib_Smile_Collat_Fwd&>(pricer_); 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(sec_); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits=  pricer.getLossUnits(); 

	if (pricer.GetTermStructurePricing()==qTermStructureR) 
	{return cpt_elt_pricer_distrib_smile_collat_fwd_TSR(yf,yf_fwd_start,pricer_,model_,sec_,losses);}

	if (collat->IsFullHomog())
	{return cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd(yf,yf_fwd_start,pricer_,model_,sec_,losses);}
	
	if ( (yf<=0) || ((yf-yf_fwd_start)<=0) ) return 0.;

	
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

	qIntegratorChoice	TheIntegratorType;
	if (! pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		pricer.SetIntegrationStep1(40) ;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());
	
	// double* Pdefault = new double[nbnames];
	// double* PdefaultStart = new double[nbnames];
	ARM_Vector Pdefault (nbnames); 
	ARM_Vector PdefaultStart  (nbnames); 

	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<nbnames; i++) 
	{
		//Estimation des Probas de défaut sur la période forward
		Pdefault[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);
		PdefaultStart[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf_fwd_start);
		if (Pdefault[i]>1. || PdefaultStart[i]>1.)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability >1. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
		if (Pdefault[i]<0.0 || PdefaultStart[i]<0.0)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"Error: Default Probability <0. for "<<collat->GetIssuersLabels(collatRank[i]));
		}
	}

	if (pricer.GetDistribution() == 0)
	{
		Distrib = new ICM_Gauss1FLossDistrib_LJ(nbnames,Pdefault,PdefaultStart,beta,
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
		Distrib->Set(nbnames,Pdefault,PdefaultStart,beta,
				lossUnits.getLossRates(date),// pricer.getLossRate(),
				pricer.GetIntegrationStep1(),
				pricer.GetCopulaType(),
				TheIntegratorType,
				pricer.GetIntegrationStep1());
	}

	// if (Pdefault) delete[] Pdefault;
	// if (PdefaultStart) delete[] PdefaultStart;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector corr_low(nbnames);
	ARM_Vector corr_Hight(nbnames);

	std::vector<double> corr_low_start;corr_low_start.resize(nbnames);
	std::vector<double> corr_Hight_start;corr_Hight_start.resize(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		if (pricer.GetTermStructurePricing()) 
		{
			corr_low_start[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf_fwd_start);
		}

		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		if (pricer.GetTermStructurePricing()) 
		{
			corr_Hight_start[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf_fwd_start);
		}
	}

	correlation->SetForcedStrikeType(qStrike_NONE);
	Distrib->UpdateCorrelation(corr_low,corr_Hight);
	if (pricer.GetTermStructurePricing()) 
	{
		Distrib->Collat_fwd_UpdateCorrelation(corr_low_start,corr_Hight_start);
	}
	Distrib->precompute_coeffs_perturb(); //Update aussi les coef Perturb pour modifier les coef Start

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
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distribCollatFwd:: Unable to compute Loss Distribution for maturity="<<yf<<" : "<<e.what());
	}
	catch(...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"cpt_elt_pricer_distrib_smileCollatFwd:: Unable to compute Loss Distribution for maturity="<<yf);
	}
	

	if (ELT<0.) ELT = 0.;
	//if (beta) delete[] beta;

	ELT /= (tranche_up-tranche_down);

	return (ELT);
}


// ---------------------------------------------------------------
// SPECIAL CDO SQUARE
// ---------------------------------------------------------------

// 2F for Sectorial CDO Square

double cpt_elt_pricer_distrib_sector(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_)
{
	return 0.;
}


