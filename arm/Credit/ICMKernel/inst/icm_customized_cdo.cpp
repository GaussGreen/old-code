/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CUSTOMIZED_CDO.H
	PROJECT:	INST
	
	DESCRIPTION:	this class provides a Cash Flow CDO Description


   -----------------------------------------------------------------
   
	ICM CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICMKernel\inst\icm_customized_cdo.h"
#include "ICMKernel\inst\icm_collateral.h"

void ICM_Customized_CDO::Init()
{
	SetName(ICM_CUSTOMIZED_CDO);
	SetCcy("EUR");

	its_Matrix_PL_Data		=	NULL;
	its_Matrix_DL_Data		=	NULL;
	its_Matrix_Parameters_Data	=	NULL;

	// DEFAULT LEG DATA
	its_DL_CreditWindowLow	=	0.0;
	its_DL_CreditWindowUp	=	5.0 * 365.0;
	its_DL_LossMin = 0.0;
	its_DL_LossMax = 0.0;
	its_DL_PaymentType = CEP_ATDEFAULTDATE;
	its_DL_PaymentDate = 0.0;
	its_DL_PaymentLag = 0;
	its_DL_NbDefMin = 0;
	its_DL_NbDefMax = 0;

	// ------------------------------------------------------
	// PREMIUM LEG DATA
	its_PL_CreditWindowLows.clear();
	its_PL_CreditWindowUps.clear();
	its_PL_StartDates.clear();
	its_PL_EndDates.clear();
	its_PL_PaymentDates.clear();

	its_PL_LossMins.clear();
	its_PL_LossMaxs.clear();
	its_PL_NbDefMins.clear();
	its_PL_NbDefMaxs.clear();
	its_PL_Ratios.clear();
	its_PL_Notios.clear();
	its_PL_Spreads.clear();

	its_PL_CreditFlags.clear();
	its_PL_CreditSpreadCaps.clear();
	its_PL_Redemptions.clear();

	// PRODUCT PRICING DATA
	its_PricingLegsType		=	CBN_DEFAULTMINUSPREMIUM;
	its_CreditObservationType	=	CO_LOSSES;
	its_ATMDataFlag			=	CPLADT_PURESPREAD;
	its_UpFront_Premium		=	0.0;
	its_UpFront_RunningSpread	=	0.0;	// 500.0 bps.
	its_CreditPremiumLegAccrued	=	CAP_NON_PRORATA;
	its_PL_PaymentType			=	CEP_ATDEFAULTDATE;
	its_Time_Step_Prorata	=	30.0;

	its_CDO_Type				=	CPT_STD_CDO;
	its_NPV_Type				=	CNPV_STANDARD;
	// ------------------------------------------------------

	its_ValDate			=	"01/01/2001";
	its_ValDateAsDouble	=	0.0;

	its_PL_NbFlows	=	0;
}


ICM_Customized_CDO::ICM_Customized_CDO(
										ICM_Parameters*	Default_Parameters,
										ICM_Parameters* Premium_Parameters,
										ICM_Parameters* Data_Parameters,
										const int& NbIssuers,
									   	char**	IssuersLabels,
										double*	IssuersNotionals,
										const std::string& ccy, 
										ARM_Date* EffectiveDate)
{
	Init();

	// ------------------------------------------------------------
	// Parameters:
	ICM_Matrix<ARM_Vector>* DL_matrix = NULL;
	ICM_Matrix<ARM_Vector>* PL_matrix = NULL;
	ICM_Matrix<ARM_Vector>* Parameters_matrix = NULL;

	if (Default_Parameters != NULL)
		DL_matrix	=	Default_Parameters->GetDblParams();

	if (Premium_Parameters != NULL)
		PL_matrix	=	Premium_Parameters->GetDblParams();

	if (Data_Parameters != NULL)
		Parameters_matrix	=	Data_Parameters->GetDblParams();
	
	Set(PL_matrix, DL_matrix, Parameters_matrix, NbIssuers, IssuersLabels, IssuersNotionals, ccy, EffectiveDate);
	
}

ICM_Customized_CDO::ICM_Customized_CDO(
										ICM_Matrix<ARM_Vector>* PL_matrix,
										ICM_Matrix<ARM_Vector>* DL_matrix,
										ICM_Matrix<ARM_Vector>* Parameters_matrix,
										const int& NbIssuers,
									   	char**	IssuersLabels,
										double*	IssuersNotionals,
										const std::string& ccy, 
										ARM_Date* EffectiveDate)
{
	Init();

	Set(PL_matrix, DL_matrix, Parameters_matrix, NbIssuers, IssuersLabels, IssuersNotionals, ccy, EffectiveDate);
	
}

void ICM_Customized_CDO::Set(	
								ICM_Matrix<ARM_Vector>* PL_matrix,
								ICM_Matrix<ARM_Vector>* DL_matrix,
								ICM_Matrix<ARM_Vector>* Parameters_matrix,
								const int& NbIssuers,
								char**	IssuersLabels,
								double*	IssuersNotionals,
								const std::string& ccy, 
								ARM_Date* EffectiveDate)
{
	SetName(ICM_CUSTOMIZED_CDO);

	if (EffectiveDate)
	{
		its_ValDate	=	(*EffectiveDate);
		its_ValDateAsDouble	=	EffectiveDate->GetJulian();
	}

	// It has to be done first
	if (its_Matrix_Parameters_Data)
		delete its_Matrix_Parameters_Data;

	if (Parameters_matrix)
		its_Matrix_Parameters_Data = (ICM_Matrix<ARM_Vector>*) Parameters_matrix->Clone();

	if (its_Matrix_PL_Data)
		delete its_Matrix_PL_Data;

	if (PL_matrix)
		its_Matrix_PL_Data = (ICM_Matrix<ARM_Vector>*) PL_matrix->Clone();

	if (its_Matrix_DL_Data)
		delete its_Matrix_DL_Data;

	if (DL_matrix)
		its_Matrix_DL_Data = (ICM_Matrix<ARM_Vector>*) DL_matrix->Clone();

	//Ccy
	SetCcy(ccy);
	
	// COLLATERAL
 

 	
	std::vector<std::string> issuers(NbIssuers); 
	ARM_Vector notios(NbIssuers); 
	for(int i=0;i<NbIssuers;i++) 
	{
		if (IssuersNotionals) notios[i]=IssuersNotionals[i]; 
		else notios[i]=10000.0;
		issuers[i]=IssuersLabels[i]; 
	}
	
 	ICM_Collateral TmpCollateral (issuers, notios);
	SetCollateral(TmpCollateral);
}


void ICM_Customized_CDO::BitwiseCopy(const ARM_Object* src)
{
    ICM_Customized_CDO* gen = (ICM_Customized_CDO*) src;

	if (gen->its_Matrix_PL_Data)
	{
		if (its_Matrix_PL_Data)
			delete its_Matrix_PL_Data;
		its_Matrix_PL_Data = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_PL_Data->Clone();
	}

	if (gen->its_Matrix_DL_Data)
	{
		if (its_Matrix_DL_Data)
			delete its_Matrix_DL_Data;
		its_Matrix_DL_Data = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_DL_Data->Clone();
	}

	if (gen->its_Matrix_Parameters_Data)
	{
		if (its_Matrix_Parameters_Data)
			delete its_Matrix_Parameters_Data;
		its_Matrix_Parameters_Data = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_Parameters_Data->Clone();
	}

	// DEFAULT LEG DATA
	its_DL_CreditWindowLow	=	gen->its_DL_CreditWindowLow;
	its_DL_CreditWindowUp	=	gen->its_DL_CreditWindowUp;
	its_DL_LossMin			=	gen->its_DL_LossMin;
	its_DL_LossMax			=	gen->its_DL_LossMax;
	its_DL_NbDefMin			=	gen->its_DL_NbDefMin;
	its_DL_NbDefMax			=	gen->its_DL_NbDefMax;
	its_DL_PaymentType		=	gen->its_DL_PaymentType;
	its_DL_PaymentDate		=	gen->its_DL_PaymentDate;
	its_DL_PaymentLag		=	gen->its_DL_PaymentLag;

	// PREMIUM LEG DATA
	its_PL_CreditWindowLows	=	gen->its_PL_CreditWindowLows;
	its_PL_CreditWindowUps	=	gen->its_PL_CreditWindowUps;
	its_PL_StartDates		=	gen->its_PL_StartDates;
	its_PL_EndDates			=	gen->its_PL_EndDates;
	its_PL_PaymentDates		=	gen->its_PL_PaymentDates;

	its_PL_LossMins			=	gen->its_PL_LossMins;
	its_PL_LossMaxs			=	gen->its_PL_LossMaxs;
	its_PL_NbDefMins		=	gen->its_PL_NbDefMins;
	its_PL_NbDefMaxs		=	gen->its_PL_NbDefMaxs;
	its_PL_Ratios			=	gen->its_PL_Ratios;
	its_PL_Notios			=	gen->its_PL_Notios;
	its_PL_Spreads			=	gen->its_PL_Spreads;

	its_PL_CreditFlags		=	gen->its_PL_CreditFlags;
	its_PL_CreditSpreadCaps	=	gen->its_PL_CreditSpreadCaps;
	its_PL_Redemptions		=	gen->its_PL_Redemptions;

	// PRODUCT PRICING DATA
	its_PricingLegsType		=	gen->its_PricingLegsType;
	its_DefNPVFlag			=	gen->its_DefNPVFlag;
	its_PremNPVFlag			=	gen->its_PremNPVFlag;

	its_CreditObservationType	=	gen->its_CreditObservationType;
	its_CreditPremiumLegAccrued	=	gen->its_CreditPremiumLegAccrued;
	
	its_ATMDataFlag			=	gen->its_ATMDataFlag;
	its_NPV_Type			=	gen->its_NPV_Type;

	its_ValDate				=	gen->its_ValDate;
	its_ValDateAsDouble		=	gen->its_ValDateAsDouble;

}


void ICM_Customized_CDO::Copy(const ARM_Object* src)
{
     ICM_Ftd::Copy(src);
 
     BitwiseCopy(src);
}


ARM_Object* ICM_Customized_CDO::Clone(void)
{
     ICM_Customized_CDO* theClone = new ICM_Customized_CDO();

     theClone->Copy(this);
 
     return(theClone);
}

// *************************************************************
// View Matrix 
// *************************************************************
void ICM_Customized_CDO::View(char* id, FILE* ficOut)
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

 	// FTD Viewer
	//fprintf(fOut, "\t\t\t ----------------- FTD Viewer ----------------- \n");

	//ICM_Ftd::View(id, fOut);

	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Premium Leg Matrix ----------------- \n");

	if (Get_Matrix_PL_Data())
		Get_Matrix_PL_Data()->View(id, fOut);

	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Default Leg Matrix ----------------- \n");

	if (Get_Matrix_DL_Data())
		Get_Matrix_DL_Data()->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- Pricing Data Matrix ----------------- \n");

	if (Get_Matrix_Parameters_Data())
		Get_Matrix_Parameters_Data()->View(id, fOut);

	if ( ficOut == NULL )
	{
	fclose(fOut);
	}
}


// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------
//	GET DEFAULT LEG INFORMATION
// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------

void	ICM_Customized_CDO::SetCreditProductDefaultMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	int	TheValue;
	if (parameters == NULL) return;

	ARM_Vector* DATA	=	NULL;
	
	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_LOW");
	if (!DATA)
		its_DL_CreditWindowLow	=	0.0;		
	else
		its_DL_CreditWindowLow = (RelativeDate) (DATA->Elt(0) - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_UP");
	if (!DATA)
		its_DL_CreditWindowUp	=	5.0 * 365.0;	//	relative to AsOf
	else
		its_DL_CreditWindowUp = (RelativeDate) (DATA->Elt(0) - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MIN: DOUBLE

	DATA	= parameters->GetColVect("LOSS_MIN");
	if (!DATA)
		its_DL_LossMin = 0.0;
	else
		its_DL_LossMin = (double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MAX: DOUBLE

	DATA	= parameters->GetColVect("LOSS_MAX");
	if (!DATA)
		its_DL_LossMax = 0.0;
	else
		its_DL_LossMax = (double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_TYPE:	RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_TYPE");
	if (!DATA)
		its_DL_PaymentType = CEP_ATDEFAULTDATE;
	else
	{
		TheValue	=	(int)	DATA->Elt(0);
		if ((TheValue < CEP_ATDEFAULTDATE) || (TheValue > CEP_ATFIXEDDATE))
			its_DL_PaymentType = CEP_ATDEFAULTDATE;
		else
			its_DL_PaymentType = (CreditEventPayment) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_DATE:	RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_DATE");
	if (!DATA)
		its_DL_PaymentDate = 0.0;	//	relative to AsOf
	else
		its_DL_PaymentDate = (RelativeDate) (DATA->Elt(0) - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT LAG: INT

	DATA	= parameters->GetColVect("PAYMENT_LAG");
	if (!DATA)
		its_DL_PaymentLag = 0;
	else
		its_DL_PaymentLag = (int) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MIN : INT

	DATA	= parameters->GetColVect("NB_MIN");
	if (!DATA)
		its_DL_NbDefMin = 0;
	else
		its_DL_NbDefMin = (int) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MAX : INT

	DATA	= parameters->GetColVect("NB_MAX");
	if (!DATA)
		its_DL_NbDefMax = 0;
	else
		its_DL_NbDefMax = (int) DATA->Elt(0);
	// ------------------------------------------------------
}



void	ICM_Customized_CDO::SetCreditProductPremiumMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;

	ARM_Vector* DATA	=	NULL;
	int	i, size, tmp_size;
	int	TheValue;

	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_LOW");
	if (!DATA)
		size	=	0;
	else
		size	=	DATA->GetSize();

	its_PL_NbFlows	=	size;
	its_PL_CreditWindowLows.resize(size);

	for (i=0;i<size;i++)
		its_PL_CreditWindowLows[i]	=	(RelativeDate) ((*DATA)[i] - its_ValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_UP");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_CreditWindowUps.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_CreditWindowUps[i]	=	(RelativeDate) ((*DATA)[i] - its_ValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// START DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("START_DATE");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_StartDates.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_StartDates[i]	=	(RelativeDate) ((*DATA)[i] - its_ValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// END DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("END_DATE");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_EndDates.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_EndDates[i]	=	(RelativeDate) ((*DATA)[i] - its_ValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_DATE");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_PaymentDates.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_PaymentDates[i]	=	(RelativeDate) ((*DATA)[i] - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS MINS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("LOSS_MIN");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_LossMins.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_LossMins[i]	=	(*DATA)[i];		

	// ------------------------------------------------------
	// LOSS MAXS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("LOSS_MAX");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_LossMaxs.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_LossMaxs[i]	=	(*DATA)[i];		

	// ------------------------------------------------------
	// NB MINS:	VECTOR OF INT

	DATA	= parameters->GetColVect("NB_MINS");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_NbDefMins.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_NbDefMins[i]	=	(int) (*DATA)[i];		

	// ------------------------------------------------------
	// NB MAXS:	VECTOR OF INT

	DATA	= parameters->GetColVect("NB_MAXS");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_NbDefMaxs.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_NbDefMaxs[i]	=	(int) (*DATA)[i];		

	// ------------------------------------------------------
	// RATIOS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("RATIOS");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_Ratios.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_Ratios[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// NOTIOS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("NOTIOS");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_Notios.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_Notios[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREADS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("SPREADS");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_Spreads.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_Spreads[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREAD CAPS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("SPREAD_CAP");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_CreditSpreadCaps.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_CreditSpreadCaps[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// REDEMPTION:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("REDEMPTION");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_Redemptions.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
		its_PL_Redemptions[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT FLAGS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("CREDIT_FLAG");
	if (!DATA)
		tmp_size	=	0;
	else
		tmp_size	=	DATA->GetSize();

	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_PL_CreditFlags.resize(tmp_size);

	for (i=0;i<tmp_size;i++)
	{
		TheValue	=	(int) (*DATA)[i];
		if ((TheValue < CPLT_GUARANTEED) || (TheValue > CPLT_CUMULREDEMPTION)) 
			TheValue	=	CPLT_OUTSTANDING;

		its_PL_CreditFlags[i]	=	(CreditPremiumLegType) TheValue;		
	}

	// ------------------------------------------------------
}


void	ICM_Customized_CDO::SetCreditProductPricingParametersMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	int	TheValue;

	ARM_Vector* DATA	=	NULL;
	
	// ------------------------------------------------------
	// DOUBLE TYPE: ValDate from Excel
	
	DATA	= parameters->GetColVect("EXCEL_VALDATE");
	if (!DATA)
		its_ValDateAsDouble	=	0.0;

	its_ValDateAsDouble = (double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICING_LEGS:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICING_LEGS");
	if (!DATA)
		its_PricingLegsType	=	CBN_DEFAULTMINUSPREMIUM;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CBN_DEFAULTMINUSPREMIUM) || (TheValue > CBN_DEFAULTLEGONLY))
			its_PricingLegsType	=	CBN_DEFAULTMINUSPREMIUM;
		else	
			its_PricingLegsType	=	(CreditBasketNPV) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_OBSERVATION:	INT --> ENUM

	DATA	= parameters->GetColVect("CREDIT_OBSERVATION");
	if (!DATA)
		its_CreditObservationType	=	CO_LOSSES;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CO_NBDEFAULTS) || (TheValue > CO_LOSSES))
			its_CreditObservationType	=	CO_LOSSES;
		else	
			its_CreditObservationType	=	(CreditObservation) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_FORMAT:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICE_FORMAT");
	if (!DATA)
		its_ATMDataFlag	=	CPLADT_PURESPREAD;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CPLADT_PURESPREAD) || (TheValue > CPLADT_UF_SPREAD))
			its_ATMDataFlag	=	CPLADT_PURESPREAD;
		else	
			its_ATMDataFlag	=	(CreditPremiumLegATMData) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_PREMIUM:	DOUBLE

	DATA	= parameters->GetColVect("UP_FRONT_PREMIUM");
	if (!DATA)
		its_UpFront_Premium	=	0.0;
	else	
		its_UpFront_Premium	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_RUNNING_SPREAD:	DOUBLE

	DATA	= parameters->GetColVect("UP_FRONT_RUNNING_SPREAD");
	if (!DATA)
		its_UpFront_RunningSpread	=	0.0;
	else
		its_UpFront_RunningSpread	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED:	INT --> ENUM

	DATA	= parameters->GetColVect("PREMIUM_ACCRUED");
	if (!DATA)
		its_CreditPremiumLegAccrued	=	CAP_NON_PRORATA;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CAP_PRORATA) || (TheValue > CAP_NON_PRORATA))
			its_CreditPremiumLegAccrued	=	CAP_NON_PRORATA;
		else
			its_CreditPremiumLegAccrued	=	(CreditAccruedPayment) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED_PAYMENT:	INT --> ENUM

	DATA	= parameters->GetColVect("PREMIUM_ACCRUED_PAYMENT");
	if (!DATA)
		its_PL_PaymentType	=	CEP_ATDEFAULTDATE;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CEP_ATDEFAULTDATE) || (TheValue > CEP_ATFIXEDDATE))
			its_PL_PaymentType	=	CEP_ATDEFAULTDATE;
		else		
			its_PL_PaymentType	=	(CreditEventPayment) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TIME_STEP_1F_PRORATA_NB_DAYS:	DOUBLE

	DATA	= parameters->GetColVect("TIME_STEP_1F_PRORATA_NB_DAYS");
	if (!DATA)
		its_Time_Step_Prorata	=	30.0;
	else
		its_Time_Step_Prorata	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COMPUTE_HEDGES:	INT: 0 or 1
/*
	DATA	= parameters->GetColVect("COMPUTE_HEDGES");
	if (!DATA)
		its_HedgesRunning	=	false;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if (TheValue)
			its_HedgesRunning	=	true;
		else
			its_HedgesRunning	=	false;
	}
*/	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRODUCT_TYPE:	INT --> ENUM

	DATA	= parameters->GetColVect("PRODUCT_TYPE");
	if (!DATA)
		its_CDO_Type	=	CPT_STD_CDO;
	else
	{
		TheValue = (int) DATA->Elt(0);
		if ((TheValue < CPT_STD_CDO) || (TheValue > CPT_CDO_SQUARE))
			its_CDO_Type	=	CPT_STD_CDO;
		
		its_CDO_Type	=	(CreditProductType) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_TYPE:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICE_TYPE");
	if (!DATA)
		TheValue =	CNPV_STANDARD;
	else
	{
		TheValue = (int) DATA->Elt(0);

		if ((TheValue < CNPV_STANDARD) || (TheValue > CNPV_WITH_RUNNING_SPREAD_AND_UP_FRONT_PREMIUM))
			its_NPV_Type	=	CNPV_STANDARD;
		else
			its_NPV_Type	=	(CreditNPV) TheValue;
	}
	// ------------------------------------------------------

}