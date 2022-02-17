#include "firsttoinc.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile_collat_fwd.h"
#include "ICMKernel/crv/icm_distriblosscalculator.h"
#include "ICMKernel/inst/icm_ftd.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include <string>

void ICM_Pricer_Distrib_Smile_Collat_Fwd::Set(ARM_Security *sec, ARM_Model *mod, const ICM_Parameters& parameters,
											  const ARM_Date&asof)
{
	ARM_Vector* INTEGRATION_METHOD = NULL;
	ICM_Security* security = (ICM_Security*) sec;

	ICM_Pricer_Distrib_Smile::Set(sec, mod,parameters,asof);
}


// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************

double ICM_Pricer_Distrib_Smile_Collat_Fwd::ExpectedLossTranche(const double& yearterm,vector<double>& losses)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	//Date de fwd start en fonction de l'effective date de la tranche
	ARM_Date Fwd_Start = ftd->GetStartDateNA();
	double yf_fwd_start = (Fwd_Start - model->GetStartDate())/365.;
	if (yf_fwd_start <=0.) yf_fwd_start = 0.;

	//Estimation de la loss distribution
	vector<double> vectlosses;
	return cpt_elt_pricer_distrib_smile_collat_fwd(yearterm, yf_fwd_start, *this,*GetModel(),*GetSecurity(),vectlosses); 
}


// *************************************************************
// Computing of Expected Loss Tranche Cas Full Homogene
// *************************************************************

double ICM_Pricer_Distrib_Smile_Collat_Fwd::ExpectedLossTrancheFullHomogene(const double& yearterm,vector<double>& losses)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	//Date de fwd start en fonction de l'effective date de la tranche
	ARM_Date Fwd_Start  = ftd->GetStartDateNA();
	double yf_fwd_start = (Fwd_Start - model->GetStartDate())/365.;
	if (yf_fwd_start <=0.) yf_fwd_start = 0.;

	//Estimation de la loss distribution
	vector<double> vectlosses;
	return cpt_elt_pricer_distrib_smile_fullhomogeneous_collat_fwd(yearterm, yf_fwd_start, *this,*GetModel(),*GetSecurity(),vectlosses);
}



void ICM_Pricer_Distrib_Smile_Collat_Fwd::View(char* id, FILE* ficOut)
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
    fprintf(fOut, "\t\t\t ------------ Homogeneous Strike Security Pricer Collat Fwd Start ----- \n");
	fprintf(fOut, "\t\t\t ---------------------------------------------------------------------- \n");

	ICM_Pricer_Distrib_Smile::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
