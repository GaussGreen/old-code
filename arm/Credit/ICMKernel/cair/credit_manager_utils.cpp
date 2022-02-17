/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_UTILS.CPP
	PROJECT:	CAIR
	
	DESCRIPTION:	Misc tools mainly for conversion

   -----------------------------------------------------------------
   
	CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */
#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\util\icm_integrator_dim_2.h"
#include "ICMKernel\glob\icm_corrmatrix.h"
#include "ICMKernel\inst\icm_collateral.h"
// ------------------------------------------------
// Fill		-->		Reset
// Update	-->		Switch according to Models		
// ------------------------------------------------

void	CreditManager :: Price_CreditProduct(ICM_Mez* Current_Product, CreditManager_Model The_Pricing_And_Copula_Choice, DoubleVector& TheParameters, DoubleVector& TheCorrelationParameters, DoubleVector& TheProductParameters)
{
	// ---------------------------
	// PART 0: CREDIT DATA
	// ---------------------------
	Fill_CreditDataParameters();
	Update_CreditDataParameters();

	// ---------------------------
	// PART 1: MODEL DATA
	// ---------------------------
	Fill_CreditDataDescription();
	Fill_CreditModelParameters();
	Fill_CreditCorrelationParameters();

	Update_CreditDataDescription();
	Update_CreditModelParameters(The_Pricing_And_Copula_Choice, TheParameters);
	Update_CreditCorrelationParameters(TheCorrelationParameters);

	// ---------------------------
	// PART 2: PRODUCT DATA - DEFAULT LEG
	// ---------------------------
	Fill_CreditProductDefault();
	Update_CreditProductDefault(Current_Product);

	// ---------------------------
	// PART 3: PRODUCT DATA - PREMIUM LEG
	// ---------------------------
	Fill_CreditProductPemium();
	Update_CreditProductPremium(Current_Product);

	// ---------------------------
	// PART 4: PRODUCT DATA - PARAMETERS
	// ---------------------------
	Fill_CreditProductPricingParameters();
	Update_CreditProductPricingParameters(Current_Product, TheProductParameters);

	// ---------------------------
	// PART 5: INTERNAL DATA - PARAMETERS
	// ---------------------------
	Fill_CreditManagerInternalData();

	// ---------------------------
	// PART 6: PRICING!!!
	// ---------------------------
	// a lot of things are redondant
	Price();

}



// ------------------------------------------------------
//	CREDIT DATA PARAMETERS
// ------------------------------------------------------

void	CreditManager :: Fill_CreditDataParameters()
{
	// ------------------------------------------------------
	// ROLL_DATE:	BOOLEAN

	its_RollDateFlag = true;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// HEDGES_TYPE:	INT --> ENUM
	
	its_Bump_Choice = CHB_NO;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD:	DOUBLE

	its_BumpSpread = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD_TYPE:	int
	
	its_BumpSpread_Type = BT_ADD;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_RECOVERY:	DOUBLE

	its_BumpRecovery = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_CORRELATION:	DOUBLE

	its_BumpCorrelation = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQUARE_NB_UNDERLYINGS:	INT

	its_CDO_Square_Nb_Underlyings = 0;
	// ------------------------------------------------------
}

void	CreditManager :: Update_CreditDataParameters()
{
	// a priori, nothing to do
	return;
}

//--------------------------------------------
//	DATA DESCRIPTION
//--------------------------------------------

void	CreditManager::Fill_CreditDataDescription()
{
	return;
}


void	CreditManager::Update_CreditDataDescription()
{
	return;
}


//--------------------------------------------
//	CREDIT MODEL PARAMETERS
//--------------------------------------------

void	CreditManager::Fill_CreditModelParameters()
{
	// ------------------------------------------------------
	// DEFAULT VALUES, THEN SWITCH
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MODEL_TYPE:	INT --> ENUM
	
	its_CreditModelType = CMT_ANALYTIC_LHP;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// KEEP_CALIBRATION:	BOOLEAN

	its_KeepCalibrationFlag = true;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CORRELATION_TYPE:	INT --> ENUM

	its_CorrelationType = CT_FLAT;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COPULA_TYPE:	INT --> ENUM
	
	its_CopulaType = CCT_GAUSSIAN;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_SIMUL:	DOUBLE

	its_NbSimul = 1000;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// its_FreedomDegree:	INT

	its_FreedomDegree = 4;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// its_NIntegration_1F:	INT

	its_NIntegration_1F = 51;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// its_N_FFT:	INT

	its_N_FFT = 512;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION:	INT --> ENUM

	its_MC_Variance_Reduction	= CMCVR_NONE;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_THETA:	DOUBLE

	its_IS_Theta = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_MU:	DOUBLE

	its_IS_Mu = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_LOSS_LEVEL:	DOUBLE

	its_IS_Loss_Level = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_LOSS_MATURITY:	DOUBLE

	its_IS_Loss_Maturity = 5.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_CHOICE:	INT --> ENUM
	
	its_Recursive_1F_Loss_Unit_Choice	= COFRLUC_PGCD;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_MIN:	DOUBLE

	its_Recursive_1F_Loss_Unit_Min = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_NB_STEP:	INT

	its_Recursive_1F_NbLossStep = 0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TARGET_DEFAULT_PROB:	DOUBLE

	its_Target_Default_Prob = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION_THETA_CHOICE:	INT --> ENUM

	its_Theta_Choice	= CVRI_IMPOSED;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION_MU_CHOICE:	INT --> ENUM

	its_Mu_Choice	= CVRF_TARGET_DEFPROB;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_MATURITY:	DOUBLE

	its_LHP_Maturity = 5.0 * 365.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_SPREAD:	DOUBLE

	its_LHP_Spread = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_SPREAD:	DOUBLE

	its_LHP_Recovery = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LP_DEGREE:	INT

	its_LP_Degree = 0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_ALPHA:	DOUBLE

	its_NIG_Alpha	= 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_BETA:	DOUBLE

	its_NIG_Beta	= 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_MU:	DOUBLE

	its_NIG_Mu = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_DELTA:	DOUBLE

	its_NIG_Delta = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_RHO:	DOUBLE

	its_NIG_Rho = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_1:	DOUBLE

	its_SC_Rho1	= 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_2:	DOUBLE

	its_SC_Rho2	= 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_PROB:	DOUBLE

	its_SC_Prob	= 0.0;
	// ------------------------------------------------------
}

//--------------------------------------------
//	CREDIT MODEL PARAMETERS
//--------------------------------------------

void	CreditManager::Fill_CreditCorrelationParameters()
{
	its_CorrelationValue	=	0.0;

	its_Beta.clear();
	its_Base_Correlation_Strikes.clear();
	its_Base_Correlation_Values.clear();

	if (its_CorrelationMatrix)
		delete its_CorrelationMatrix;
	its_CorrelationMatrix	=	NULL;
	
	if (its_CholeskyMatrix)
		delete its_CholeskyMatrix;
	its_CholeskyMatrix	=	NULL;

	if (itsCorrelation)
		delete itsCorrelation;
	itsCorrelation	=	NULL;

	its_FL_Alpha.clear();
	its_FL_Beta1.clear();
	its_FL_Beta2.clear();

	its_Used_Beta.clear();
	its_Used_SQRT_OneMinusBetaSquare.clear();
}


void	CreditManager::Update_CreditCorrelationParameters(DoubleVector& TheCorrelationParameters)
{
	double	TheParam;

	CreditManager_Correlation_Type	TheCorrelationType;

	TheParam	=	TheCorrelationParameters[0];

	TheCorrelationType	=	(CreditManager_Correlation_Type) ((int) TheParam);
	switch (TheCorrelationType)
	{
	case CMCT_COMPOUND_CORRELATION:

		its_CorrelationValue	=	TheCorrelationParameters[1];
		break;

	case CMCT_BASE_CORRELATION:
		// to do
			//TheCorrelationParameters[1]
			//TheCorrelationParameters[2]
			//TheCorrelationParameters[3];
		break;
		
	case CMCT_NIG_CORRELATION:

		// Retrieve Model Data
		its_NIG_Beta	=	TheCorrelationParameters[2];

		its_NIG_Alpha	=	sqrt(exp(TheCorrelationParameters[1]) + its_NIG_Beta*its_NIG_Beta);

		its_NIG_Rho		=	cos(2.0 * TheCorrelationParameters[3] - 1.0);

		break;

	}
}


void	CreditManager::Update_CreditModelParameters(CreditManager_Model The_Pricing_And_Copula_Choice, DoubleVector& TheParameters)
{
	switch (The_Pricing_And_Copula_Choice)
	{
	case CMM_LHP_GAUSSIAN:

		its_CreditModelType	=	CMT_ANALYTIC_LHP;
//		its_CopulaType		=	CCT_GAUSSIAN;			// redondant with Default Settings

		its_LHP_Maturity	=	TheParameters[0];
		its_LHP_Spread		=	TheParameters[1];
		its_LHP_Recovery	=	TheParameters[2];
		
		break;
		
	case CMM_LHP_JPM_GAUSSIAN:
		its_CreditModelType	=	CMT_ANALYTIC_LHP_JPM;
		its_CopulaType		=	CCT_GAUSSIAN;			// redondant with Default Settings

		its_LHP_Maturity	=	TheParameters[0];
		its_LHP_Spread		=	TheParameters[1];
		its_LHP_Recovery	=	TheParameters[2];
		break;
		
	case CMM_LHP_JPM_NIG:
		its_CreditModelType	=	CMT_ANALYTIC_LHP_JPM;
		its_CopulaType		=	CCT_NIG;

		its_LHP_Maturity	=	TheParameters[0];
		its_LHP_Spread		=	TheParameters[1];
		its_LHP_Recovery	=	TheParameters[2];
		break;
		
	case CMM_HT_JPM_GAUSSIAN:
			break;

	case CMM_RSCB:
		break;
		
	case CMM_ANDERSEN_GAUSSIAN:
		its_CreditModelType	=	CMT_ANALYTIC_RECURSIVE_1F;
		break;
	}
}


void	CreditManager::Fill_CreditProductDefault()
{
	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	RELATIVE DATE

	its_DL_CreditWindowLow = 0.0;	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	RELATIVE DATE

	its_DL_CreditWindowUp = 5.0 * 365.0;	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MIN: DOUBLE

	its_DL_LossMin = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MAX: DOUBLE

	its_DL_LossMax = 0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_TYPE:	RELATIVE DATE

	its_DL_PaymentType = CEP_ATDEFAULTDATE;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_DATE:	RELATIVE DATE

	its_DL_PaymentDate = 0.0;	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT LAG: INT

	its_DL_PaymentLag = 0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MIN : INT

	its_DL_NbDefMin = 0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MAX : INT

	its_DL_NbDefMax = 0;
	// ------------------------------------------------------
}


void	CreditManager::Update_CreditProductDefault(ICM_Mez* Current_Product)
{
	// no NULL product
	ICM_Leg* TheDefLeg;
	ARM_Date Asof = GetAsOfDate();

	TheDefLeg	=	Current_Product->GetDefLeg();

//	its_DL_CreditWindowDown	=	0.0;
	its_DL_CreditWindowUp	=	TheDefLeg->GetEndDate() - Asof;

	its_DL_LossMin	=	Current_Product->GetPercentLow(Current_Product->GetFeeLeg()->GetStartDate());
	// HIGH!!!
	its_DL_LossMax	=	its_DL_LossMin + Current_Product->GetPercentHight(Current_Product->GetFeeLeg()->GetStartDate());
}



void	CreditManager::Fill_CreditProductPemium()
{
	int	size;
	size	=	0;

	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	VECTOR OF RELATIVE DATE

	its_PL_NbFlows	=	0;	
	
	its_PL_CreditWindowLows.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	VECTOR OF RELATIVE DATE

	its_PL_CreditWindowUps.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// START DATES:	VECTOR OF RELATIVE DATE

	its_PL_StartDates.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// END DATES:	VECTOR OF RELATIVE DATE

	its_PL_EndDates.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT DATES:	VECTOR OF RELATIVE DATE

	its_PL_PaymentDates.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS MINS:	VECTOR OF DOUBLE

	its_PL_LossMins.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS MAXS:	VECTOR OF DOUBLE

	its_PL_LossMaxs.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB MINS:	VECTOR OF INT

	its_PL_NbDefMins.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB MAXS:	VECTOR OF INT

	its_PL_NbDefMaxs.resize(size);
	// ------------------------------------------------------
	
	// ------------------------------------------------------
	// RATIOS:	VECTOR OF DOUBLE

	its_PL_Ratios.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NOTIOS:	VECTOR OF DOUBLE

	its_PL_Notios.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREADS:	VECTOR OF DOUBLE

	its_PL_Spreads.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREAD CAPS:	VECTOR OF DOUBLE

	its_PL_CreditSpreadCaps.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// REDEMPTION:	VECTOR OF DOUBLE
	
	its_PL_Redemptions.resize(size);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT FLAGS:	VECTOR OF DOUBLE

	its_PL_CreditFlags.resize(size);
	// ------------------------------------------------------
}


void	CreditManager::Update_CreditProductPremium(ICM_Mez* Current_Product)
{
	// no NULL product
	ICM_Leg* TheFeeLeg;
	ARM_Date Asof = GetAsOfDate();

	TheFeeLeg	=	Current_Product->GetFeeLeg();

	const ARM_Vector*		TheStartDates;
	const ARM_Vector*		TheEndDates;
	const ARM_Vector*		TheRatios;
	const ARM_Vector*		TheNotionals;
	const ARM_Vector*		TheSpreads;

	ICM_Security* TheSecurity;
	
	TheSecurity		=	TheFeeLeg->GetCreditInfos();

	TheStartDates	=	&TheSecurity->GetAccStartDates();
	TheEndDates		=	&TheSecurity->GetAccEndDates();
	TheNotionals	=	&TheSecurity->GetNotionals();
	TheSpreads		=	&TheSecurity->GetCouponRates();

	// To be Improved?
	// if (! TheSecurity->IsYFAlreadyComputed())
	TheSecurity->ComputeYF(Asof);
	
	TheRatios	=	&TheSecurity->GetYFInterestDays();

	int	i, size;

	size	=	TheStartDates->GetSize();
	
	// size
	its_PL_NbFlows	=	size;

	its_PL_CreditWindowLows.resize(size);
	its_PL_CreditWindowUps.resize(size);
	its_PL_StartDates.resize(size);
	its_PL_EndDates.resize(size);
	its_PL_PaymentDates.resize(size);
	its_PL_LossMins.resize(size);
	its_PL_LossMaxs.resize(size);

	its_PL_NbDefMins.resize(size);
	its_PL_NbDefMaxs.resize(size);
	its_PL_Ratios.resize(size);
	its_PL_Notios.resize(size);
	its_PL_Spreads.resize(size);
	its_PL_CreditSpreadCaps.resize(size);
	its_PL_Redemptions.resize(size);
	its_PL_CreditFlags.resize(size);


	double	PctLow;
	double	PctUp;

	PctLow	=	Current_Product->GetPercentLow(Current_Product->GetFeeLeg()->GetStartDate());
	PctUp	=	PctLow + Current_Product->GetPercentHight(Current_Product->GetFeeLeg()->GetStartDate());

	double	AsOf_Julian;
	AsOf_Julian	=	Asof.GetJulian();

	for (i=0; i<size; i++)
	{
		its_PL_CreditWindowLows[i]	=	0.0;
		its_PL_CreditWindowUps[i]	=	(RelativeDate) ((*TheEndDates)[i] - AsOf_Julian);		//	relative to AsOf
		its_PL_StartDates[i]		=	(RelativeDate) ((*TheStartDates)[i] - AsOf_Julian);	//	relative to AsOf
		its_PL_EndDates[i]			=	(RelativeDate) ((*TheEndDates)[i] - AsOf_Julian);		//	relative to AsOf
		its_PL_PaymentDates[i]		=	(RelativeDate) ((*TheEndDates)[i] - AsOf_Julian);		//	relative to AsOf

		its_PL_LossMins[i]			=	PctLow;
		its_PL_LossMaxs[i]			=	PctUp;

		its_PL_NbDefMins[i]			=	0;	
		its_PL_NbDefMaxs[i]			=	0;

		its_PL_Ratios[i]			=	(*TheRatios)[i];		

		its_PL_Notios[i]			=	(*TheNotionals)[i];		

		its_PL_Spreads[i]			=	(*TheSpreads)[i];	// 100.0		

		its_PL_CreditSpreadCaps[i]	=	0.0;
		its_PL_Redemptions[i]		=	0.0;
		
		its_PL_CreditFlags[i]		=	CPLT_OUTSTANDING;

	}

}


void	CreditManager::Fill_CreditProductPricingParameters()
{
	
	// ------------------------------------------------------
	// PRICING_LEGS:	INT --> ENUM

	its_PricingLegsType	=	CBN_DEFAULTMINUSPREMIUM;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_OBSERVATION:	INT --> ENUM

	its_CreditObservationType	=	CO_LOSSES;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_FORMAT:	INT --> ENUM

	its_ATMDataFlag	=	CPLADT_PURESPREAD;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_PREMIUM:	DOUBLE

	its_UpFront_Premium	=	0.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_RUNNING_SPREAD:	DOUBLE

	its_UpFront_RunningSpread	=	0.0;	// 500.0 bps.
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED:	INT --> ENUM

	its_CreditPremiumLegAccrued	=	CAP_NON_PRORATA;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED_PAYMENT:	INT --> ENUM

	its_PL_PaymentType	=	CEP_ATDEFAULTDATE;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TIME_STEP_1F_PRORATA_NB_DAYS:	DOUBLE

	its_Time_Step_Prorata_1F	=	30.0;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COMPUTE_HEDGES:	INT: 0 or 1

	its_HedgesRunning	=	false;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRODUCT_TYPE:	INT --> ENUM

	its_CDO_Type	=	CPT_STD_CDO;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_TYPE:	INT --> ENUM		// not yet in interface

	its_NPV_Type	=	CNPV_STANDARD;
	// ------------------------------------------------------

}


void	CreditManager::Update_CreditProductPricingParameters(ICM_Mez* Current_Product, DoubleVector& TheProductParameters)
{
	int		UpFrontFlag;
	double	TheValue;

	TheValue	=	TheProductParameters[0];
	UpFrontFlag	=	(int) TheValue;

	if (UpFrontFlag)
	{
		its_NPV_Type	=	CNPV_WITH_UP_FRONT_PREMIUM;
		its_UpFront_Premium	=	TheProductParameters[1];

		its_ATMDataFlag	=	CPLADT_UF_UPFRONT;
		its_UpFront_RunningSpread	=	its_PL_Spreads[0];
	}
	else
		its_ATMDataFlag	=	CPLADT_PURESPREAD;

/*	if ()
	{
		its_ATMDataFlag		=	CPLADT_UF_UPFRONT;
		its_UpFront_Premium	=	0.0;

		its_UpFront_RunningSpread	=	500.0;
	}
*/
}


void	CreditManager::Fill_CreditManagerInternalData()
{
	its_IsActivateShift	=	false;
}


void	CreditManager::Fill_DataFrom_ModelMultiCurves(ARM_Security *option, ARM_Model *mod)
{
	if ((mod == NULL) || (option == NULL)) return;

	ICM_Ftd*	The_Security;
	ICM_ModelMultiCurves*	The_Model;

	The_Security	=	dynamic_cast<ICM_Ftd*> (option);
	if (!The_Security)
		ICMTHROW(ERR_INVALID_DATA," Collateral product ");
	The_Model		=	((ICM_ModelMultiCurves*) mod);
	
	// GET ALL DATA
	const std::vector<std::string>& The_Labels	=	The_Security->GetCollateral()->GetIssuersLabels();
	// Nb credits
	double NbCredits	=	The_Security->GetCollateral()->GetNbIssuers();
	if (NbCredits == 0)
		ICMTHROW(ERR_INVALID_DATA,"No Credit! in Fill_DataFrom_ModelMultiCurves!");
	
	const ICM_DefaultCurve*	The_Def_Curve	=	The_Model->GetDefaultCurve(0);

	if (The_Def_Curve == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL Default Curve! in Fill_DataFrom_ModelMultiCurves!");
	
	// ROLL DATE
	its_RollDateFlag	=	true;

	// CDO SQUARE
	its_CDO_Square_Nb_Underlyings	=	0;

	// HEDGES PARAMETERS TO BE SET ELSEWHERE
	its_BumpSpread		=	10.0;
	its_BumpSpread_Type	=	BT_ADD;
	its_BumpRecovery	=	0.1;
	its_BumpCorrelation	=	0.1;

	// NOT REALLY USED
	its_Categories.clear();
	its_Currencies.clear();
	its_Accrueds.clear();
	its_Recoveries.clear();
	its_Notionals.clear();
	its_Input_Losses.clear();
	its_DefaultDates.clear();
	its_AmortizationDates.clear();

	its_Categories.resize(NbCredits);
	its_Currencies.resize(NbCredits);
	its_Accrueds.resize(NbCredits);
	its_Recoveries.resize(NbCredits);
	its_Notionals.resize(NbCredits);
	its_Input_Losses.resize(NbCredits);
	its_DefaultDates.resize(NbCredits);
	its_AmortizationDates.resize(NbCredits);

	int	i;
	// char* Current_Label;
	
	double	TheRecovery;
	double	TheNotional;
	double	TheValue;

	ARM_Vector theLosses(NbCredits);
	// theLosses	=	(double *)malloc(Get_NbCredits()*sizeof(double));

	// double* notionals	=	The_Security->GetIssuersNotionals();

	if (its_ArrayDefaultCrv)
		for (int i=0; i<its_CreditsLabels.size(); i++){	
			if(its_ArrayDefaultCrv[i]) delete 
				its_ArrayDefaultCrv[i];
		}
		delete [] its_ArrayDefaultCrv;
	its_ArrayDefaultCrv	=	NULL;

	its_ArrayDefaultCrv = new ICM_DefaultCurve*[NbCredits];
	std::vector<const ICM_DefaultCurve*> items(NbCredits); 
	for (i=0; i<NbCredits; i++)
	{
		// Label
		const std::string& Current_Label	=The_Labels[i];

		TheRecovery	=	The_Model->GetRecoveryRate(Current_Label);
		TheNotional	=	The_Security->GetCollateral()->GetIssuersNotional(i); // notionals[i];
		
		its_Notionals[i]	=	TheNotional;
		its_Recoveries[i]	=	TheRecovery;

		TheValue			=	(1.0 - TheRecovery);
		its_Input_Losses[i]		=	TheValue;
		TheValue			*=	TheNotional;
		theLosses[i]		=	TheValue;
		
		The_Def_Curve		=	The_Model->GetDefaultCurve(Current_Label);

		items[i] = its_ArrayDefaultCrv[i] = (ICM_DefaultCurve*) The_Def_Curve->Clone();

	}

	Set_CreditsLabels(The_Labels);
/*
	// ugly
	char	The_Terms[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];

	data	=	NULL;
	its_Credit_Market_Data->Get_CreditDataMaturitiesAsChar(data);

	Set_CreditDataMaturitiesAsChar(The_Terms);

	// Spreads
	its_CreditDataSpreads	=	its_Credit_Market_Data->Get_CreditDataSpreads();

	GenerateAllDefaultCurves();
*/
	// ------------------------------------------------------------

	// Create MultiCurveModel
	if (its_ModelMultiCurves != NULL)
		delete its_ModelMultiCurves;

	its_ZeroCurve	=	The_Model->GetZeroCurve();
	itsCorrelation	=	(ICM_Correlation*) (The_Model->GetCorrelation()->Clone());

	
	its_ModelMultiCurves = new ICM_ModelMultiCurves(items,
											// Get_NbCredits(),
											// its_ArrayDefaultCrv,
											its_ZeroCurve, 
											theLosses,
											itsCorrelation);

	// if (theLosses != NULL)
	//	free(theLosses);

}



/** 
17783 
extern "C" void 
dummy_function_to_integrate(void* Param, double x, double y, double& Result)
{
	qIntegratorChoice	*TheIntegrationMethod;
	
	TheIntegrationMethod	=	(qIntegratorChoice*)((*(AddressVector*)Param)[0]);

	if (*TheIntegrationMethod	==	qGAUSS_HERMITE)
	{
		x	*=	SQRT_TWO;
		y	*=	SQRT_TWO;
	}

	Result = x * y;		//	1.0;
}
*/ 

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
/** 
17783 
void	CreditManager::Test_Integration_Dim()
{
	ICM_Integrator_Dim_2	TheIntegrator;

	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;

	double	x_min	=	-6.0;
	double	x_max	=	6.0;
	int		x_nbsteps	=	21;

	double	y_min	=	-6.0;
	double	y_max	=	6.0;
	int		y_nbsteps	=	21;

	qIntegratorChoice	TheIntegratorChoice;

	TheIntegratorChoice	=	qGAUSS_LEGENDRE;

	TheIntegrator.SetIntegrationType(TheIntegratorChoice);
	TheIntegrator.SetIntegrationStep(x_nbsteps);
	TheIntegrator.SetLowBound(x_min);
	TheIntegrator.SetUpBound(x_max);

	TheIntegrator.SetIntegrationType_2(TheIntegratorChoice);
	TheIntegrator.SetIntegrationStep_2(y_nbsteps);
	TheIntegrator.SetLowBound_2(y_min);
	TheIntegrator.SetUpBound_2(y_max);

	double	Value;

	TheParameterVector.Reset();
	TheParameterVector.Append(&TheIntegratorChoice);

	TheIntegrator.Integrate(x_min, x_max, y_min, y_max, dummy_function_to_integrate, &TheParameterVector, Value);

	ICMLOG("Result for Integrator - GAUSS LEGENDRE : " << Value);

	TheIntegratorChoice	=	qGAUSS_HERMITE;

	TheParameterVector.Reset();
	TheParameterVector.Append(&TheIntegratorChoice);

	x_nbsteps	=	20;
	TheIntegrator.SetIntegrationType(TheIntegratorChoice);
	TheIntegrator.SetIntegrationStep(x_nbsteps);
	TheIntegrator.SetLowBound(x_min);
	TheIntegrator.SetUpBound(x_max);

	y_nbsteps	=	20;
	TheIntegrator.SetIntegrationType_2(TheIntegratorChoice);
	TheIntegrator.SetIntegrationStep_2(y_nbsteps);
	TheIntegrator.SetLowBound_2(y_min);
	TheIntegrator.SetUpBound_2(y_max);

	TheIntegrator.Integrate(x_min, x_max, y_min, y_max, dummy_function_to_integrate, &TheParameterVector, Value);

	ICMLOG("Result for Integrator - GAUSS HERMITE : " << Value);}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


17783 */ 