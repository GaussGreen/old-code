#include "firsttoinc.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile_rf.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_RF.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\inst\icm_ftd.h"
#include "ICMKernel\mod\modelmulticurves.h"


//	temp here
std::vector<int> 
CptSortedInd(const ICM_Collateral& collat,const ICM_LossUnits&lossUnits,const ARM_Date&date,const ICM_Correlation&corr) ; 


void ICM_Pricer_Distrib_RandFactors::Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters,
										 const ARM_Date&asof)
{

	//la matrice parameters du pricer doit contenir les champs suivant :
	//RF_PARAMS
	//RF_MODEL

	ARM_Vector* INTEGRATION_METHOD = NULL;
	ICM_Security* security = (ICM_Security*) sec;

	ICM_Pricer_Distrib_Smile::Set(sec, mod,parameters,asof);
	ARM_Vector* RFParams = parameters.GetColVect("RF_PARAMS");
	itsRFParameters.resize(RFParams->GetSize());
	for (int i=0; i<RFParams->GetSize();i++) {itsRFParameters[i] = RFParams->Elt(i);}

	ICM_QMatrix<string>* RFModel = parameters.GetColVectStr("STR_RF_MODEL");
	
	if ((*RFModel)(0,0)=="PWL") itsModelType = qCAL_PWL_CORREL;
	if ((*RFModel)(0,0)=="PWC") itsModelType = qCAL_PWC_CORREL;
	
}


// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************
double ICM_Pricer_Distrib_RandFactors::ExpectedLossTranche(const double& yearterm,vector<double>& losses)
{
	if (yearterm<=0) return 0.;

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();
	const ICM_LossUnits& lossUnits= getLossUnits(); 
	//Full Homogene Case
	if (ftd->GetCollateral()->IsFullHomog()) return ExpectedLossTrancheFullHomogene(yearterm);

	//for sensivities
// 17783 	if (IsActivateShift()) return ExpectedLossTrancheForFastSpreadHedge(yearterm);

	ARM_Date Maturity = ftd->GetEndDateNA();
	double YF_Maturity = (Maturity-model->GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_RF* Distrib = NULL;

	if (!ftd->GetCollateral()->SumNotionals(ftd->GetStartDateNA())) return 0.;

	ARM_Date date = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yearterm);
	double tranche_down = GetTranche_Down(date);
	double tranche_up = GetTranche_Up(date);

	tranche_down = fabs(round(tranche_down));
	tranche_up = fabs(round(tranche_up));

	if (tranche_down>=tranche_up) return 0.; 

	int nbnames = ftd->GetCollateral()->GetNbIssuers();
	// double* beta = new double[nbnames];	
	ARM_Vector* VBeta = GetBetas(); //For correlation Matrix compatibility
	// if (VBeta) memcpy(beta,VBeta->GetElt(),sizeof(double)*nbnames);

	double ELT = 0.;

	qIntegratorChoice	TheIntegratorType;
	int IntegrationStep;GetIntegrationStep1(IntegrationStep);
	int CopulaType;GetCopulaType(CopulaType);

	if (!IntegrationStep) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		SetIntegrationStep1(40);
	}
	else	
		TheIntegratorType	=	GetGaussLegendreMethodFromIntegrationStep(IntegrationStep);
	
	ARM_Vector Pdefault (nbnames); // = new double[nbnames];
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<nbnames; i++) 
	{	Pdefault[i] = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(collatRank[i]))->DefaultProba(yearterm);
		if (Pdefault[i]>1.)
			ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability >1. for "<<ftd->GetCollateral()->GetIssuersLabels(collatRank[i])); 
		if (Pdefault[i]<0.0)
			ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability <0. for "<<ftd->GetCollateral()->GetIssuersLabels(collatRank[i])); 
	}

	if (GetDistribution() == NULL)
	{
		//Pour l'instant GaussHermite obligatoire si Long Short
		TheIntegratorType =qGAUSS_HERMITE;
		Distrib = new ICM_Gauss1FLossDistrib_RF(nbnames,Pdefault,*VBeta,
			getLossUnits().getLSLossRates(date),// getLongShort_LossRate(), 
			CptSortedInd(*ftd->GetCollateral(),getLossUnits(),date,*model->GetCorrelation()), // *GetSortedInd(),
			itsModelType,itsRFParameters,IntegrationStep,
												CopulaType,TheIntegratorType,IntegrationStep);


	SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = (ICM_Gauss1FLossDistrib_RF*)GetDistribution();
		//Pour l'instant GaussHermite obligatoire si Long Short
		TheIntegratorType =qGAUSS_HERMITE;
		Distrib->Set(nbnames,Pdefault,*VBeta,
				getLossUnits().getLSLossRates(date),// getLongShort_LossRate(), 
				CptSortedInd(*ftd->GetCollateral(),getLossUnits(),date,*model->GetCorrelation()) ,//*GetSortedInd(),
					itsModelType,itsRFParameters,IntegrationStep,
					CopulaType,TheIntegratorType,IntegrationStep);
		
	}

	// if (Pdefault) delete[] Pdefault;
	// if (!its_lossunit_flg) 	computelossunit();	


	//Insertion du vecteur d'indice triés
	ELT =Distrib->compute_expectedlosstranche_LongShort(tranche_up,
													tranche_down,
													getLossUnits().getLossUnit(date)// GetLossUnit(),
													// GetMaxShortLoss()
//													GetShiftMatrix(),
													// GetTenorShift(),
													// GetIssuerShift()
													); 

	if (ELT<0.) ELT = 0.;
	//if (beta) delete[] beta;

	ELT /= (tranche_up-tranche_down);

	return (ELT);
}


// *************************************************************
// Computing of Expected Loss Tranche Cas Full Homogene
// *************************************************************
double ICM_Pricer_Distrib_RandFactors::ExpectedLossTrancheFullHomogene(const double& yearterm)
{
	if (yearterm<=0) return 0.;

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	//for sensivities
// 17783 	if (IsActivateShift()) return ExpectedLossTrancheForFastSpreadHedge(yearterm);

	ARM_Date Maturity = ftd->GetEndDateNA();
	double YF_Maturity = (Maturity-model->GetStartDate())/365.;
	ICM_Gauss1FLossDistrib_RF* Distrib = NULL;

	if (!ftd->GetCollateral()->SumNotionals(ftd->GetStartDateNA())) return 0.;

	ARM_Date date = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yearterm);
	double tranche_down = GetTranche_Down(date);
	double tranche_up = GetTranche_Up(date);

	tranche_down = fabs(round(tranche_down));
	tranche_up = fabs(round(tranche_up));

	if (tranche_down>=tranche_up) return 0.; 

	int nbnames = ftd->GetCollateral()->GetNbIssuers();
	double ELT = 0.;

	//Integrator
	qIntegratorChoice	TheIntegratorType;
	int IntegrationStep;GetIntegrationStep1(IntegrationStep);
	int CopulaType;GetCopulaType(CopulaType);

	if (!IntegrationStep) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		SetIntegrationStep1(40);
	}
	else	
		TheIntegratorType	=	GetGaussLegendreMethodFromIntegrationStep(IntegrationStep);
	
	//Default Probability
	double Pdef = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(0))->DefaultProba(yearterm);
	if (Pdef>1.)
		ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability >1. for "<<ftd->GetCollateral()->GetIssuersLabels(0)); 
	if (Pdef<0.0)
		ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability <0. for "<<ftd->GetCollateral()->GetIssuersLabels(0)); 

	//Distribution
	if (GetDistribution() == NULL)
	{
		Distrib = new ICM_Gauss1FLossDistrib_RF(nbnames,
												itsModelType,itsRFParameters,IntegrationStep,
												CopulaType,TheIntegratorType,IntegrationStep);
		SetDistribution((ICM_Distribution*)Distrib);
	}
	else
	{
		Distrib = (ICM_Gauss1FLossDistrib_RF*)GetDistribution();
		Distrib->Set(nbnames,
					itsModelType,itsRFParameters,IntegrationStep,
					CopulaType,TheIntegratorType,IntegrationStep);
	}

	//	if (!its_lossunit_flg) 	computelossunit();	
	// SetLossUnit(CheckLossUnit_(its_lossunit,tranche_up,tranche_down));

	ELT = Distrib->compute_expectedlosstrancheFullHomog(getLossUnits().getLossUnit(date),
												tranche_down,
												tranche_up,
												Pdef);

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);
	return (ELT);
}



void ICM_Pricer_Distrib_RandFactors::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 

	int size =0;

	fprintf(fOut, "\t\t\t ---------------------------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t ------------ Homogeneous Strike Security Pricer Range Factor --------- \n");
	fprintf(fOut, "\t\t\t ---------------------------------------------------------------------- \n");

	switch (itsModelType)
	{	
	case qCAL_PWC_CORREL:
		fprintf(fOut, "\tModel: PWC_CORREL\n");
		break;
	case qCAL_PWL_CORREL:
		fprintf(fOut, "\tModel: PWL_CORREL\n");
		break;
	default:
		fprintf(fOut, "\tIntegrationMethod: UNKNOWN !\n");
	}
	
	ICM_Pricer_Distrib_Smile::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
