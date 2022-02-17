#include "ICMKernel\inst\icm_cpdo.h"


ICM_cpdo::ICM_cpdo(ICM_Leg* RiskyLeg, 
			ICM_Leg* RollLeg,
			ICM_Leg* NoRiskyLeg,
			const double& InitialValo,
			const double& Target,
			const ARM_Date& Maturity,
			const std::string& CpnType,
			const double& UFFees,
			const double& RunningFees,
			const double& VExpo,
			const double& V0Expo,
			const double& Alpha,
			const double& Beta,
			const double& Desactivation,
			const int& NbAssets)
{
	Init();
	Set(RiskyLeg, 
		RollLeg,
		NoRiskyLeg,
		InitialValo,
		Target,
		Maturity,
		CpnType,
		UFFees,
		RunningFees,
		VExpo,
		V0Expo,
		Alpha,
		Beta,
		Desactivation, 
		NbAssets);
}

void ICM_cpdo::Init()
{
	SetName(ICM_CPDO);

	itsRiskyLeg = NULL;
	itsRollLeg = NULL;
	itsNoRiskyLeg = NULL;
	
	// Caracteristique du produit
	itsInitialValo = 100.0;
	itsTarget = 2.0;
	
	itsUFFees = 0.01;
	itsRunningFees = 0.0025;

	itsVExpo= 0.25;
	itsV0Expo = 0.2;

	itsAlpha = 5.;
	itsBeta = 0.;

	itsDesactivation = 0.1;

	itsNbAssets = 125; 

}

void ICM_cpdo::Set(ICM_Leg* RiskyLeg, 
			ICM_Leg* RollLeg,
			ICM_Leg* NoRiskyLeg,
			const double& InitialValo,
			const double& Target,
			const ARM_Date& Maturity,
			const std::string& CpnType,
			const double& UFFees,
			const double& RunningFees,
			const double& VExpo,
			const double& V0Expo,
			const double& Alpha,
			const double& Beta,
			const double& Desactivation,
			const int& NbAssets)
{

	if (itsRiskyLeg)
		delete itsRiskyLeg;
	itsRiskyLeg = dynamic_cast<ICM_Leg*>(RiskyLeg->Clone());

	if (itsRollLeg)
		delete itsRollLeg;
	itsRollLeg = dynamic_cast<ICM_Leg*>(RollLeg->Clone());

	if (itsNoRiskyLeg)
		delete itsNoRiskyLeg;
	itsNoRiskyLeg = dynamic_cast<ICM_Leg*>(NoRiskyLeg->Clone());

	itsInitialValo  = InitialValo;
	itsTarget = Target;
	itsCPDOMaturity = Maturity; 
	itsCpnType = CpnType;
	itsUFFees = UFFees;
	itsRunningFees = RunningFees;
	itsVExpo = VExpo;
	itsV0Expo = V0Expo;
	itsAlpha = Alpha;
	itsBeta = Beta;
	itsDesactivation = Desactivation;
	itsNbAssets = NbAssets; 

}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_cpdo::BitwiseCopy (const ARM_Object* srcCPDO)
{
	ICM_cpdo* cpdo = (ICM_cpdo*) srcCPDO;

	itsInitialValo	= cpdo->itsInitialValo;
	itsTarget		= cpdo->itsTarget;
	itsCPDOMaturity	= cpdo->itsCPDOMaturity; 
	itsCpnType		= cpdo->itsCpnType;

	itsUFFees		= cpdo->itsUFFees;
	itsRunningFees	= cpdo->itsRunningFees;
	
	itsVExpo	= cpdo->itsVExpo;
	itsV0Expo = cpdo->itsV0Expo;
	
	itsAlpha = cpdo->itsAlpha;
	itsBeta = cpdo->itsBeta;
	
	itsDesactivation = cpdo->itsDesactivation;
	
	itsNbAssets = cpdo->itsNbAssets; 

	if (cpdo->itsRiskyLeg)
	{
		if (itsRiskyLeg)
			delete itsRiskyLeg;
		itsRiskyLeg = NULL;

		itsRiskyLeg = dynamic_cast<ICM_Leg*>(cpdo->itsRiskyLeg->Clone());
	}

	if (cpdo->itsRollLeg)
	{
		if (itsRollLeg)
			delete itsRollLeg;
		itsRollLeg = NULL;

		itsRollLeg = dynamic_cast<ICM_Leg*>(cpdo->itsRollLeg->Clone());
	}

	if (cpdo->itsNoRiskyLeg)
	{
		if (itsNoRiskyLeg)
			delete itsNoRiskyLeg;
		itsNoRiskyLeg = NULL;

		itsNoRiskyLeg = dynamic_cast<ICM_Leg*>(cpdo->itsNoRiskyLeg->Clone());
	}


}

// -------------
//	Copy Method 
// -------------
void ICM_cpdo::Copy(const ARM_Object* srcCPDO)
{
	ICM_Security::Copy(srcCPDO);
	BitwiseCopy(srcCPDO);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_cpdo::Clone(void)
{
	ICM_cpdo* theClone = new ICM_cpdo();

	theClone->Copy(this);

	return(theClone);
}

// --------------
//	View Method
// --------------
void ICM_cpdo::View(char* id, FILE* ficOut)
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

	char d1[20];

	itsCPDOMaturity.JulianToStrDate(d1);


	fprintf(fOut, "\t\t\t ------------------ CPDO Informations ----------------- \n\n\n");

	fprintf(fOut, " Initial valo : %f \n",itsInitialValo);
	fprintf(fOut, " Target : %f \n",itsTarget);
	fprintf(fOut, " Maturity :  %14s \n",d1);
	fprintf(fOut, " Coupon Type : %s \n",itsCpnType.c_str());
	fprintf(fOut, " Up Front fees: %f \n",itsUFFees);
	fprintf(fOut, " Running fees : %f \n",itsRunningFees);
	fprintf(fOut, " Exposition on V : %f \n",itsVExpo);
	fprintf(fOut, " Exposition on V0 : %f \n",itsV0Expo);
	fprintf(fOut, " Alpha : %f \n",itsAlpha);
	fprintf(fOut, " Beta : %f \n",itsBeta);
	fprintf(fOut, " Desactivation : %f \n",itsDesactivation);
	fprintf(fOut, " Nb of Assets : %d \n",itsNbAssets);

	
	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ----------------- Risky Leg Informations ----------------- \n\n");
	itsRiskyLeg->View(id,fOut);

	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ----------------- Roll Leg Informations ----------------- \n\n");
	itsRollLeg->View(id,fOut);

	fprintf(fOut, "\n\n");
	fprintf(fOut, "\t\t\t ----------------- No Risky Leg Informations ----------------- \n\n");
	itsNoRiskyLeg->View(id,fOut);

	 if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
