#include "firsttoinc.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile.h"
#include "ICMKernel\util\ICM_pgcd.h"
#include "ICMKernel/crv/icm_distriblosscalculator.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel/mod/icm_distribution.h"

ICM_Pricer_Distrib_Smile::~ICM_Pricer_Distrib_Smile(void) 
	{
		if (itsDistribution)
			delete itsDistribution;
		itsDistribution = NULL;
	}
void ICM_Pricer_Distrib_Smile::Init() 
{
		ICM_Pricer_Distrib::Init();
		SetName(ICM_PRICER_HOMOGENEOUS_SMILE);

		itsCopulaType = qGAUSSIAN;
		itsIntegrationMethod = qGAUSS_HERMITE;
		itsIntegrationStep1 = INTEGRATION_STEP;
		itsFreedomDegree = 0;
		itsDistribution = NULL;
		itsRescalType = qRescal_Std_Maturity;
		itsApproxComputation=0.;
}

void ICM_Pricer_Distrib_Smile::Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters,
								   const ARM_Date&asof)
{

	//la matrice parameters du pricer doit contenir les champs suivant :
	//COPULA = 0                          'Integration Method
    //INTEGRATION_METHOD= 0
    //FREEDOM_DEGREE = 0
    //INTEGRATION_STEP_1 = 0

	ARM_Vector* INTEGRATION_METHOD = NULL;
	ICM_Security* security = (ICM_Security*) sec;

	ICM_Pricer_Distrib::Set(sec, mod,parameters,asof);

	computelossunit();


	if (!parameters.empty()) 
	{
	// if (parameters != NULL)
	// {
		// SetParameters((ICM_Parameters*)parameters->Clone());	
		SetParameters(parameters); 
		// Get Parameters
		ExtractParameters();

		INTEGRATION_METHOD = parameters.GetColVect("INTEGRATION_METHOD");
	}
	// }

	// how to say whether I have Hermite or Legendre with n: odd or even!
	if (INTEGRATION_METHOD == NULL)
		itsIntegrationMethod = GetGaussLegendreMethodFromIntegrationStep(itsIntegrationStep1);

}


qIntegratorChoice ICM_Pricer_Distrib_Smile::GetGaussLegendreMethodFromIntegrationStep(const int& TheIntegrationStep)
{

	switch(itsIntegrationMethod)
	{
	case 2:
		return qTRAPEZE;
		break;
	case 0:
	case 1:
	default:
		if (TheIntegrationStep % 2 == 0) // modulo
			// even number --> General Hermite
			return qGAUSS_HERMITE;
		else
			// odd number --> General Legendre
			return qGAUSS_LEGENDRE;
		break;
	}
}

void ICM_Pricer_Distrib_Smile::ExtractParameters()
{
	const ICM_Parameters& FlowsMatrix = GetParameters();

	// if (FlowsMatrix == NULL) return;

	// ------------------------------------------------------
	// COPULA TYPE: GAUSSIAN, STUDENT
	ARM_Vector* COPULA = FlowsMatrix.GetColVect("COPULA");

	if (COPULA) itsCopulaType = (int) COPULA->Elt(0);
	// ------------------------------------------------------

	// Afin de distinguer les differentes méthodes numériques, il faut que INTEGRATION_METHOD soit renseigné
	// ------------------------------------------------------
	// INTEGRATION METHOD: GAUSS-LEGENDRE, HERMITE, TRAPEZE
	ARM_Vector* INTEGRATION_METHOD = FlowsMatrix.GetColVect("INTEGRATION_METHOD");

	if (INTEGRATION_METHOD) itsIntegrationMethod = (int) INTEGRATION_METHOD->Elt(0);

	// ------------------------------------------------------
	// ------------------------------------------------------
	// INTEGRATION STEPS: for Numerical Integrator
	ARM_Vector* INTEGRATION_STEP_1 = FlowsMatrix.GetColVect("INTEGRATION_STEP_1");

	if (INTEGRATION_STEP_1) itsIntegrationStep1 = (int) INTEGRATION_STEP_1->Elt(0);
	// ------------------------------------------------------

	// On verifie que le nombre de pas et la méthode employée sont compatibles
	// et on change la méthode s'il le faut.
	if ( (itsIntegrationStep1 % 2 == 0) && (itsIntegrationMethod == 0))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep1 % 2 == 1) && (itsIntegrationMethod == 1))
		itsIntegrationMethod = qGAUSS_LEGENDRE;

	ARM_Vector* VResc = FlowsMatrix.GetColVect("TERMS_RESCALING");
	if (VResc) 
	{
		qTERM_STRUCTURE type;
		bool ret;
		ICM_EnumsCnv::cnv(VResc->Elt(0),type , ret);
		SetTermStructurePricing(type);
	}	

	// ------------------------------------------------------
	// FREEDOM_DEGREE: for Numerical Integrator
	ARM_Vector* RESCAL_TYPE = FlowsMatrix.GetColVect("RESCAL_TYPE");

	if (RESCAL_TYPE)
	{	
		int value = (int)RESCAL_TYPE->Elt(0);
		bool res;
		ICM_EnumsCnv::cnv(value,itsRescalType, res);
	}else {	itsRescalType = qRescal_Std_Maturity;}
	// ------------------------------------------------------

	ARM_Vector* DISTRIB_APPROX = FlowsMatrix.GetColVect("DISTRIB_APPROX");

	if (DISTRIB_APPROX)
	{	
		double value = DISTRIB_APPROX->Elt(0);
		SetApproxComputation(value);
	}

}


// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************
double ICM_Pricer_Distrib_Smile::ExpectedLossTranche(const double& yearterm,vector<double>& losses)
{
	return cpt_elt_pricer_distrib_smile(yearterm,*this,*GetModel(),*GetSecurity(),losses); 
};


// *************************************************************
// Computing of Expected Loss Tranche in case Full Homogeneous
// *************************************************************
/**
double ICM_Pricer_Distrib_Smile::ExpectedLossTrancheFullHomogeneous(const double& yearterm,vector<double>& losses)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	ARM_Date Maturity = ftd->GetEndDateNA();
	double YF_Maturity = (Maturity-model->GetStartDate())/365.;
	ICM_Gauss1FLossDistrib* Distrib = NULL;

	ARM_Date date = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yearterm);
	double tranche_down = GetTranche_Down(date);
	double tranche_up = GetTranche_Up(date);

	tranche_down = fabs(round(tranche_down));
	tranche_up = fabs(round(tranche_up));

	int nbnames = ftd->GetNbIssuers();
	double ELT = 0.,Pdefault=0.,beta_down=0.,beta_up=0.;

	Pdefault= model->GetDefaultCurve(ftd->GetIssuersLabels(0))->DefaultProba(yearterm);

  	// Smile Parameters
	ICM_Correlation* correlation = model->GetCorrelation();
	if (!GetFaster()) correlation->ComputeStrikesEq(this);

	correlation->SetForcedStrikeType(qStrike_LOW);
	beta_down = correlation->GetBeta(ftd->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,CREDIT_DEFAULT_VALUE);
	correlation->SetForcedStrikeType(qStrike_UP);
	beta_up = correlation->GetBeta(ftd->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,CREDIT_DEFAULT_VALUE);
	correlation->SetForcedStrikeType(qStrike_NONE);

	if (itsDistribution == NULL)
	{
		Distrib = new ICM_Gauss1FLossDistrib(nbnames);
		SetDistribution((ICM_Distribution*)Distrib);
	}
	else
		Distrib = (ICM_Gauss1FLossDistrib*)GetDistribution();

	if (!its_lossunit_flg) 	computelossunit();	

	try{
	ELT =Distrib->ComputeEL_FullHomog(its_lossunit,tranche_down,tranche_up,beta_down,beta_up,Pdefault,itsIntegrationStep1,losses);
	}
	catch (...)
		{ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Distrib_Smile :: unable to compute Expected Loss in Full Homogeneous case");}

	if (ELT<0.) ELT = 0.;
	ELT /= (tranche_up-tranche_down);

	return (ELT);
};
**/ 

// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************

// 17783 double ICM_Pricer_Distrib_Smile::ExpectedLossTrancheForFastSpreadHedge(const double& yearterm)
// 17783 {
// 17783 	return cpt_elt_pricer_distrib_smile_fast(yearterm,*this,*GetModel(),*GetSecurity()); 
// 17783 }


void ICM_Pricer_Distrib_Smile::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Strike Security Pricer ----------------- \n");

	fprintf(fOut, "\tCopulaType: %i\n",itsCopulaType);
	fprintf(fOut, "\tFreedomDegree: %i\n",itsFreedomDegree);
	fprintf(fOut, "\tIntegrationStep1: %i\n",itsIntegrationStep1);
	string resS(""); ICM_EnumsCnv::toString(itsRescalType, resS);
	int resI=0; ICM_EnumsCnv::toInt(itsRescalType, resI);
	fprintf(fOut, "\tRescalType: %s : corresponds to %d \n\n", resS.c_str(), resI);
	switch (itsIntegrationMethod)
	{	
	case 0:
		fprintf(fOut, "\tIntegrationMethod: GAUSS LEGENDRE\n");
		break;
	case 1:
		fprintf(fOut, "\tIntegrationMethod: GAUSS HERMITE\n");
		break;
	case 2:
		fprintf(fOut, "\tIntegrationMethod: TRAPEZE\n");
		break;
	default:
		fprintf(fOut, "\tIntegrationMethod: UNKNOWN !\n");
	}
	
	ICM_Pricer_Distrib::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}

void	ICM_Pricer_Distrib_Smile::SetDistribution(ICM_Distribution* distrib)	
	{
		if (itsDistribution)
			delete itsDistribution;
		itsDistribution = distrib;
	}
