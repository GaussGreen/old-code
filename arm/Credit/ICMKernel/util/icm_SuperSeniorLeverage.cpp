	

#include "ICMKernel/util/icm_superseniorleverage.h"

#include "ICMKernel/pricer/icm_pricer_homogeneous_smile.h"

#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include "ICMKernel/crv/icm_defaultcurve.h"
#include "ICMKernel/glob/icm_maths.h"
// #include <nag.h>

#define	NB_ITER_MAX	50

// ------------------------------------------------------
ICM_SuperSeniorLeverage::~ICM_SuperSeniorLeverage() 
{ 
	if (its_Pricer)
		delete its_Pricer;
	its_Pricer = NULL;

	if (its_Matrix_Multiples)
		delete its_Matrix_Multiples;
	its_Matrix_Multiples = NULL;

	if (its_MultiCurvesModel)
		delete its_MultiCurvesModel;
	its_MultiCurvesModel = NULL;

	if (its_Min_Matrix_Multiples)
		delete its_Min_Matrix_Multiples;
	its_Min_Matrix_Multiples = NULL;

	if (its_Max_Matrix_Multiples)
		delete its_Max_Matrix_Multiples;
	its_Max_Matrix_Multiples = NULL;
}
void ICM_SuperSeniorLeverage::BitwiseCopy(const ARM_Object * src)
{
	ICM_SuperSeniorLeverage* SSL = (ICM_SuperSeniorLeverage *) src;

	its_leverage		=	SSL->Get_Leverage();
	its_trigger_pct		=	SSL->Get_Trigger_Pct();
	its_trigger_tol		=	SSL->Get_Trigger_Tolerance();
	its_target_spread	=	SSL->Get_Target_Spread();
	its_loss_level		=	SSL->Get_Loss_Level();
	its_loss_step		=	SSL->Get_Loss_Step();
	its_maturity_step	=	SSL->Get_Maturity_Step();
	its_min_down_multiple	=	SSL->Get_Min_Down_Multiple();
	its_max_down_multiple	=	SSL->Get_Max_Down_Multiple();
	its_min_high_multiple	=	SSL->Get_Min_High_Multiple();
	its_max_high_multiple	=	SSL->Get_Max_High_Multiple();
	its_Pricing_Type	=	SSL->Get_PricingType();
	its_Pricing_Method	=	SSL->Get_PricingMethod();
	its_cutoff_loss		=	SSL->Get_CutOff_Loss();
	its_cutoff_maturity	=	SSL->Get_CutOff_Maturity();
	
//JLA: no clone method on pricers... 		its_Pricer				=	(ICM_Pricer*) (SSL->GetPricer())->Clone();
//		its_Matrix_Multiples	=	(ICM_Pricer*) (SSL->GetMatrix_Multiples())->Clone();

//		Vect_of_MaturitiesInYF	=	SSL->GetVectOfMaturities()
}
void	ICM_SuperSeniorLeverage::SetPricer(ICM_Pricer* data)
{
	if (its_Pricer)
		delete its_Pricer;
//		its_Pricer = (ICM_Pricer*) data->Clone();
	its_Pricer =	data;
}
void ICM_SuperSeniorLeverage::Set_Matrix_Flags(ICM_QMatrix<double>* value) 
{ 
	if (its_Matrix_Flags)
		delete its_Matrix_Flags;
	its_Matrix_Flags = value; 
}

void ICM_SuperSeniorLeverage::Set_Multiples_Input(ICM_QMatrix<double>* value) 
{ 
	if (its_Matrix_MultiplesInput)
		delete its_Matrix_MultiplesInput;
	its_Matrix_MultiplesInput = value; 
}
ICM_QMatrix<double>* ICM_SuperSeniorLeverage::GetMatrix_Outputs()
	{
		switch (its_Computation_Type)
		{
		case SSL_CC_MULTIPLES:
			return its_Matrix_Multiples;

		case SSL_CC_TRIGGERS:
			return its_Matrix_Triggers;

		default:
			return	NULL;
		}
	}
// ------------------------------------------------------
// virtual 
void ICM_SuperSeniorLeverage::SetTriggeredCorrelation(ICM_Correlation* correl) 
{
	if (its_Triggered_Correlation)
		delete its_Triggered_Correlation;
	its_Triggered_Correlation = (ICM_Correlation*) correl->Clone();
}

void	ICM_SuperSeniorLeverage::SetDataParameters(ICM_Matrix<ARM_Vector>* parameters) 
{
	int	TheValue;

	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;

	// ------------------------------------------------------
	// LEVERAGE:	DOUBLE

	THE_VECT	= parameters->GetColVect("LEVERAGE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad LEVERAGE");

	its_leverage = (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TRIGGER_PCT:		DOUBLE

	THE_VECT	= parameters->GetColVect("TRIGGER_PCT");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad TRIGGER_PCT");

	its_trigger_pct = (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TRIGGER_TOL:		DOUBLE

	THE_VECT	= parameters->GetColVect("TRIGGER_TOL");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad TRIGGER_TOL");

	its_trigger_tol = (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TARGET_SPREAD:		DOUBLE

	THE_VECT	= parameters->GetColVect("TARGET_SPREAD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad TARGET_SPREAD");

	its_target_spread	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_LEVEL:		DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS_LEVEL");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad LOSS_LEVEL");

	its_loss_level	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_STEP:		DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS_STEP");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad LOSS_STEP");

	its_loss_step	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MATURITY_STEP:		DOUBLE

	THE_VECT	= parameters->GetColVect("MATURITY_STEP");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad MATURITY_STEP");

	its_maturity_step	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MIN_DOWN_MULTIPLE:		DOUBLE

	THE_VECT	= parameters->GetColVect("MIN_DOWN_MULTIPLE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad MIN_DOWN_MULTIPLE");

	its_min_down_multiple	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MAX_DOWN_MULTIPLE:		DOUBLE

	THE_VECT	= parameters->GetColVect("MAX_DOWN_MULTIPLE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad MAX_DOWN_MULTIPLE");

	its_max_down_multiple	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MIN_HIGH_MULTIPLE:		DOUBLE

	THE_VECT	= parameters->GetColVect("MIN_HIGH_MULTIPLE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad MIN_HIGH_MULTIPLE");

	its_min_high_multiple	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MAX_HIGH_MULTIPLE:		DOUBLE

	THE_VECT	= parameters->GetColVect("MAX_HIGH_MULTIPLE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad MAX_HIGH_MULTIPLE");

	its_max_high_multiple	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICING_TYPE:		ENUM

	THE_VECT	= parameters->GetColVect("PRICING_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad PRICING_TYPE");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < SSL_PV) || (TheValue > SSL_SPREAD)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PRICING_TYPE");
	
	its_Pricing_Type = (SSL_PricingType) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CUTOFF_LOSS:		INT

	THE_VECT	= parameters->GetColVect("CUTOFF_LOSS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad CUTOFF_LOSS");

	its_cutoff_loss = (int)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CUTOFF_MATURITY:		INT

	THE_VECT	= parameters->GetColVect("CUTOFF_MATURITY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad CUTOFF_MATURITY");

	its_cutoff_maturity = (int)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICING METHOD: DICHOTOMY or BRENT

	THE_VECT	= parameters->GetColVect("PRICING_METHOD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad PRICING_METHOD");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < SSL_PM_DICHOTOMY) || (TheValue > SSL_PM_BRENT)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PRICING_TYPE");
	
	its_Pricing_Method = (SSL_PricingMethod) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COMPUTATION_TYPE : MULTIPLES or LEVERAGES

	THE_VECT	= parameters->GetColVect("COMPUTATION_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad PRICING_METHOD");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < SSL_CC_MULTIPLES) || (TheValue > SSL_CC_TRIGGERS)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad COMPUTATION_TYPE");
	
	its_Computation_Type = (SSL_ComputationChoice) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_TRIGGER:		DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS_TRIGGER");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad LOSS_TRIGGER");

	its_loss_trigger	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CORRELATION_TRIGGERED:		DOUBLE

	THE_VECT	= parameters->GetColVect("CORRELATION_TRIGGERED");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad CORRELATION_TRIGGERED");

	its_correlation_triggered	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// INITIAL_MTM:		DOUBLE

	THE_VECT	= parameters->GetColVect("INITIAL_MTM");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad INITIAL_MTM");

	its_initial_mtm	= (double)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_FLAGS_ROWS:		DOUBLE

	THE_VECT	= parameters->GetColVect("NB_FLAGS_ROWS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad NB_FLAGS_ROWS");

	its_Nb_Flags_Rows	= (int)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_FLAGS_COLS:		DOUBLE

	THE_VECT	= parameters->GetColVect("NB_FLAGS_COLS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad NB_FLAGS_COLS");

	its_Nb_Flags_Cols	= (int)	THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// GET_DATA : MULTIPLES or LEVERAGES

	THE_VECT	= parameters->GetColVect("GET_DATA");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  bad GET_DATA");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < SSL_GD_NONE) || (TheValue > SSL_GD_TRIGGERS)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad GET_DATA");
	
	its_GetData_Type = (SSL_GetData) TheValue;
	// ------------------------------------------------------

}


double	ICM_SuperSeniorLeverage::ComputeLeverage()
{
	int	i;
	int	iloss;

	int	NbSpreadsMaturities;

	const ARM_Vector	*Current_Spreads	=	NULL;

	// -------------------------------------------
//FILE *stream = fopen("c:\\test\\ObjectViewer.txt", "w+");
	// -------------------------------------------

	// SOMEHOW AD-HOC
	if (its_GetData_Type != SSL_GD_NONE)
		return 1.0;

	if (its_Pricer == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL PRICER for Super Senior Leverage Pricing!");

	// Curves Spreads Shifts
	const ICM_DefaultCurve*		Current_DefaultCurve = NULL;

	its_MultiCurvesModel	=	(ICM_ModelMultiCurves*) its_Pricer->GetModel();
	if (its_MultiCurvesModel == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL MULTI CURVE MODEL for Super Senior Leverage Pricing!");

	its_NbCredits	=	its_MultiCurvesModel->GetNbDefCurves();
	if (its_NbCredits == 0)
		ICMTHROW(ERR_INVALID_DATA,"No Credit in MULTI CURVE MODEL for Super Senior Leverage Pricing!");

	// Store all initial Spreads
	Current_DefaultCurve	=	its_MultiCurvesModel->GetDefaultCurve(0);
	if (Current_DefaultCurve == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL DEFAULT CURVE for Super Senior Leverage Pricing!");
	NbSpreadsMaturities		=	Current_DefaultCurve->GetSize();
	if (NbSpreadsMaturities == 0)
		ICMTHROW(ERR_INVALID_DATA,"No Spreads Maturities in DEFAULT CURVE for Super Senior Leverage Pricing!");

	// -------------------------------------------------------------------------
	// Keep them in a Matrix

	its_Initial_Spreads	=	new ICM_QMatrix<double>(its_NbCredits, NbSpreadsMaturities);

	for (i=0; i<its_NbCredits; i++)
	{
		Current_DefaultCurve	=	its_MultiCurvesModel->GetDefaultCurve(i);		 

		Current_Spreads			=	Current_DefaultCurve->GetRates();
		its_Initial_Spreads->SetRow(i, Current_Spreads);
	}

	// ---------------
	// VIEW
//	its_Initial_Spreads->View("",stream);
	// ---------------

	// -------------------------------------------------------------------------
	// Products building
	its_ValDate	=	its_Pricer->GetAsOfDate();

	ICM_Mez*	MyProduct;
	MyProduct		=	(ICM_Mez*)	its_Pricer->GetSecurity();
	if (MyProduct == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL PRODUCT for Super Senior Leverage Pricing!");

	its_Notional	=	MyProduct->GetCollateral()->SumNotionals(MyProduct->GetStartDateNA());

	// The Target is the spread (in bps)
	MyProduct->SetCreditSpread(its_target_spread / 100);

	// Tranche Size
	its_Tranche_Size	=	MyProduct->GetMezzAmount(MyProduct->GetFeeLeg()->GetStartDate());

	// ------------------------------------------------------------------------------------------------
	// PRODUCT MODIFICATION
	// ------------------------------------------------------------------------------------------------
	double	current_loss;
	double	SubAmount;
	double	MezzAmount;

	its_Init_Detach	=	MyProduct->GetPercentLow(MyProduct->GetFeeLeg()->GetStartDate());
	its_Init_Attach	=	MyProduct->GetPercentHight(MyProduct->GetFeeLeg()->GetStartDate());		// sorry, just High!

	current_loss	=	0.0;
	SubAmount	=	FMAX(its_Notional * (its_Init_Detach - current_loss), 0.0);
	MezzAmount	=	its_Notional * (FMAX(its_Init_Attach + its_Init_Detach - current_loss, 0.0) - FMAX(its_Init_Detach - current_loss, 0.0));

	MyProduct->SetSubAmount(SubAmount);
	MyProduct->SetMezzAmount(MezzAmount);

	// ------------------------------------------------------------------------------------------------
	// Price
	its_Pricer->ResetRootPricer();
	((ICM_Pricer_Distrib_Smile*) its_Pricer)->clearDistribLoss();
//	current_price	=	((ICM_Pricer_Distrib_Smile*) its_Pricer)->ComputePrice(its_ValDate);

	double	multiple;
	double	min_multiple;
	double	max_multiple;

	double	f_multiple;

	int	nbiter;

	PreciseChrono chrono;

	// ----------------------------------------------------------------------------
//	if ((its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+")) == NULL) return -9.9;
//	fprintf(its_fOut, " ----------------- Super Senior Leverage ----------------- \n");
//	fclose(its_fOut);
	// ----------------------------------------------------------------------------

	// -------------------------------------------------------------------------------------------------
	// MATURITIES
	// -------------------------------------------------------------------------------------------------

	ICM_Leg*	MyDefLeg	=	NULL;

	switch (its_Pricing_Type)
	{
	case SSL_PV:
		its_PricingMode	=	qCMPPRICE;

		break;
	
	case SSL_DEF_LEG:
		its_PricingMode	=	qCMPDEFLEGPV;

		break;
	
	case SSL_DEF_LEG_ZC:
		its_PricingMode	=	qCMPDEFLEGPV;

		MyDefLeg	=	MyProduct->GetDefLeg();
		if (MyDefLeg == NULL)
			ICMTHROW(ERR_INVALID_DATA,"NULL DEF LEG for Super Senior Leverage Pricing!");

		MyDefLeg->SetPaymentFreq(K_ZEROCOUPON);
		MyDefLeg->CptCashFlowDatesCredit();

		break;
	
	case SSL_SPREAD:
		its_PricingMode	=	qCMPSPREAD;

		break;
	}

	ICM_Mez*	MyCloneProduct;	

	ARM_Date	The_Maturity;
	ARM_Date	Current_Maturity;
	char*		Matu= new char[10];

	The_Maturity	=	MyProduct->GetEndDateNA();

	Current_Maturity	=	The_Maturity;
	
	string	MaturityLabel;
	char	buffer[5];

	itoa(-its_maturity_step, buffer, 10);
	MaturityLabel.append(buffer);
	MaturityLabel.append("M");

	strcpy(Matu, MaturityLabel.c_str());

	ARM_Currency	*MyCurrency;
	MyCurrency	=	MyProduct->GetCurrencyUnit();
	if (MyCurrency == NULL)
		ICMTHROW(ERR_INVALID_DATA,"NULL CURRENCY for Super Senior Leverage Pricing!");

	int	iMat;
	iMat = 0;

	its_Vect_of_Maturities.clear();
	its_VectofProducts.clear();
	its_Vect_of_MaturitiesInYF.clear();
	
	its_Vect_of_Used_Losses.clear();
	its_Vect_of_Used_MaturitiesInYF.clear();
	
	while (Current_Maturity > its_ValDate)
	{
		MyCloneProduct	=	(ICM_Mez*)(MyProduct->Clone());

		// MyCloneProduct->SetMaturity(Current_Maturity);	
		MyCloneProduct->CptCashFlowDatesCredit(Current_Maturity);

		// ---------------
		// VIEW
//		MyCloneProduct->View("",stream);
		// ---------------
			
		its_Vect_of_Maturities.push_back(Current_Maturity);
		its_VectofProducts.push_back(MyCloneProduct);
		its_Vect_of_MaturitiesInYF.push_back(MATHTIME(Current_Maturity - its_ValDate));

		// somehow ad-hoc, because AdjDateToCdsCalendard gives the next CDS Roll date, even if the
		// current Date is a CDS Roll Date
		Current_Maturity	= AddPeriod(Current_Maturity - 1, Matu, MyCurrency->GetCcyName(), true, qCredit_Adjust20);
		ICMLOG("Current Maturity: " << Current_Maturity);

		iMat++;
	}
	
//	fclose(stream);

	// -------------------------------------------------------------------------------------------------
	// MATURITIES
	// -------------------------------------------------------------------------------------------------
	
	ARM_Date	MyDate;

	its_NbMat	=	iMat;

	current_loss	=	0.0;

	its_NbLosses	=	(int) (its_loss_level / its_loss_step) + 1; 

	// just one loss level
	if (its_cutoff_maturity)
		its_NbMat	=	FMIN(its_NbMat, its_cutoff_maturity);

	if (its_cutoff_loss)
		its_NbLosses	=	FMIN(its_NbLosses, its_cutoff_loss);

	for (i=0; i<its_NbLosses; i++)
		its_Vect_of_Used_Losses.push_back(its_loss_step * i);

	for (i=0; i<its_NbMat; i++)
		its_Vect_of_Used_MaturitiesInYF.push_back(its_Vect_of_MaturitiesInYF[i]);
	//	TheLosses	=	new double[its_NbLosses];

	// --------------------------------------
	// MAYBE I WILL HAVE TO KEEP RATHER THAN
	
	if (its_Matrix_Triggers)
		delete its_Matrix_Triggers;
	its_Matrix_Triggers	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat);
	
	// --------------------------------------
/*
	bool	TestCalibrationFlag;

	TestCalibrationFlag	=	false;

	if (TestCalibrationFlag)
		TestCalibrationWithGuess();
*/
	bool	BrentFlag;

	BrentFlag = (its_Pricing_Method == SSL_PM_BRENT) ? true:false;


	// -------------------------------------------------------------------------
	switch (its_Computation_Type)
	{
	case SSL_CC_MULTIPLES:
		break;

	case SSL_CC_TRIGGERS:
			ComputeTriggersFromMultiples();
			return 0.0;

		break;
	}

	// -------------------------------------------------------------------------
	double	tmpvalue;

	if (its_min_down_multiple > its_max_down_multiple)
	{
		tmpvalue	=	its_min_down_multiple;
		its_min_down_multiple	=	its_max_down_multiple;
		its_max_down_multiple	=	tmpvalue;
	}
	
	if (its_min_high_multiple > its_max_high_multiple)
	{
		tmpvalue	=	its_min_high_multiple;
		its_min_high_multiple	=	its_max_high_multiple;
		its_max_high_multiple	=	tmpvalue;
	}
	// -------------------------------------------------------------------------

	if (its_Matrix_Multiples)
		delete its_Matrix_Multiples;

	its_Matrix_Multiples	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat);

	if (its_Min_Matrix_Multiples)
		delete its_Min_Matrix_Multiples;

	its_Min_Matrix_Multiples	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, its_min_down_multiple);

	if (its_Max_Matrix_Multiples)
		delete its_Max_Matrix_Multiples;

	its_Max_Matrix_Multiples	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, its_max_high_multiple);

//	fprintf(its_fOut, "Nb Maturities:\t%u\t\tNb Losses:\t%u\n", its_NbMat, its_NbLosses);

	if (its_Matrix_Multiples_Flags)
		delete its_Matrix_Multiples_Flags;

	its_Matrix_Multiples_Flags	=	new ICM_QMatrix<bool>(its_NbLosses, its_NbMat, false);

	if (its_Matrix_Multiples_Products)
		delete its_Matrix_Multiples_Products;

	its_Matrix_Multiples_Products	=	new ICM_QMatrix<ICM_Mez*>(its_NbLosses, its_NbMat, NULL);

	if (its_Min_Matrix_Control_Multiples)
		delete its_Min_Matrix_Control_Multiples;

	its_Min_Matrix_Control_Multiples	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, its_min_down_multiple);

	if (its_Max_Matrix_Control_Multiples)
		delete its_Max_Matrix_Control_Multiples;

	its_Max_Matrix_Control_Multiples	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, its_max_high_multiple);

	if (its_Min_Matrix_Control_Values)
		delete its_Min_Matrix_Control_Values;

	its_Min_Matrix_Control_Values	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, 0.0);

	if (its_Max_Matrix_Control_Values)
		delete its_Max_Matrix_Control_Values;

	its_Max_Matrix_Control_Values	=	new ICM_QMatrix<double>(its_NbLosses, its_NbMat, 0.0);

	// -------------------------------------------------------------------------
	// --------------------------------------

	bool IsDecreasingMaturity;
	bool IsIncresingLoss;

	IsDecreasingMaturity	=	false;
	IsIncresingLoss			=	true;

	// -------------------------------------------------------------------------
	// THIS PROCESS MUST ENABLE US TO ADJUST
	// MIN and MAX more precisely

//	ComputeMinAndMaxMultipleForARange(its_NbMat, its_NbLosses, new_min_multiple, new_max_multiple);

	// -------------------------------------------------------------------------
	// ALGORITHM
	// -------------------------------------------------------------------------
	// RAZ
//	its_Max_Matrix_Multiples->SetValue(
	// FIRST INITIALIZATION
//	UpdateMaxMultiples(0, its_NbMat-1, 0, its_NbLosses-1, its_max_high_multiple);
//	UpdateMinMultiples(0, its_NbMat-1, 0, its_NbLosses-1, its_min_down_multiple);
	// and then just specific
	if (its_NbMat > 1)
	{
		UpdateMinMultiples(its_NbMat-1, its_NbMat-1, 0, its_NbLosses-1, its_min_high_multiple);
		UpdateMaxMultiples(0, 0, 0, its_NbLosses-1, its_max_down_multiple);
	}

	// -------------------------------------------------------------------------
	// FIRST POINT
	// -------------------------------------------------------------------------
	// first point (0, its_NbMat-1) --> MAX
	iloss	=	0;
	iMat	=	its_NbMat-1;

	PrepareProductForPricing(iloss, iMat);

	// ---------------
	// TIMER BEGIN
	chrono.Start();
	// ---------------

//	its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//	fprintf(its_fOut, "Maturity:\t%u\t\tLoss:\t%u\n", iMat, iloss);
	ICMLOG("Maturity " << iMat << " - Loss " << iloss);

	// ---------------
	// PRICING CORE
	min_multiple	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
	max_multiple	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

	if (BrentFlag)
		GetOneMultiple_Brent(min_multiple, max_multiple, nbiter, multiple, f_multiple);
	else
		GetOneMultiple(min_multiple, max_multiple, nbiter, multiple);
	
	// ---------------

	(*its_Matrix_Multiples)(iloss, iMat)	=	multiple;
	(*its_Matrix_Triggers)(iloss, iMat)		=	f_multiple;
	
	// UPDATE
	UpdateMaxMultiples(1, its_NbMat-1, 0, its_NbLosses-1, multiple);

	// ---------------
	// TIMER END
	chrono.Stop();
	// ---------------
//	fprintf(its_fOut, "Computation done in %.0lf seconds\n", chrono.GetDurationInSeconds());
	ICMLOG("Computation done in %.0lf seconds" << chrono.GetDurationInSeconds() << ".");

//	fclose(its_fOut);

	// ------------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------	

	// Only one point
	if ((its_NbLosses == 1) && (its_NbMat == 1))
		return multiple;

	// -------------------------------------------------------------------------
	// SECOND POINT
	// -------------------------------------------------------------------------
	// (NbLosses-1, 0) --> MIN
	iloss	=	its_NbLosses-1;
	iMat	=	0;

	// ------------------------------------------------------------------------------------------------
	PrepareProductForPricing(iloss, iMat);

	// ---------------
	// TIMER BEGIN
	chrono.Start();
	// ---------------

//	its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//	fprintf(its_fOut, "Maturity:\t%u\t\tLoss:\t%u\n", iMat, iloss);
	ICMLOG("Maturity " << iMat << " - Loss " << iloss);

	// ---------------
	// PRICING CORE
	min_multiple	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
	max_multiple	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

	if (BrentFlag)
		GetOneMultiple_Brent(min_multiple, max_multiple, nbiter, multiple, f_multiple);
	else
		GetOneMultiple(min_multiple, max_multiple, nbiter, multiple);
	// ---------------

	(*its_Matrix_Multiples)(iloss, iMat)	=	multiple;
	(*its_Matrix_Triggers)(iloss, iMat)		=	f_multiple;
	
	// UPDATE
	UpdateMinMultiples(0, its_NbMat-2, 0, its_NbLosses-1, multiple);

	// ---------------
	// TIMER END
	chrono.Stop();
	// ---------------
//	fprintf(its_fOut, "Computation done in %.0lf seconds\n", chrono.GetDurationInSeconds());
	ICMLOG("Computation done in %.0lf seconds" << chrono.GetDurationInSeconds() << ".");

//	fclose(its_fOut);

	// ------------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------
			
	// LAST MATURITY for all remaining LOSSES
	iMat	=	its_NbMat-1;
	for (iloss=its_NbLosses-1; iloss>0; iloss--)
	{
		// ------------------------------------------------------------------------------------------------
		PrepareProductForPricing(iloss, iMat);

		// ---------------
		// TIMER BEGIN
		chrono.Start();
		// ---------------

//		its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//		fprintf(its_fOut, "Maturity:\t%u\t\tLoss:\t%u\n", iMat, iloss);
		ICMLOG("Maturity " << iMat << " - Loss " << iloss);

		// ---------------
		// PRICING CORE
		min_multiple	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
		max_multiple	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

		if (BrentFlag)
			GetOneMultiple_Brent(min_multiple, max_multiple, nbiter, multiple, f_multiple);
		else
			GetOneMultiple(min_multiple, max_multiple, nbiter, multiple);
		// ---------------

		(*its_Matrix_Multiples)(iloss, iMat)	=	multiple;
		(*its_Matrix_Triggers)(iloss, iMat)		=	f_multiple;
		
		// UPDATE
		UpdateMaxMultiples(1, iMat, iloss, iloss, multiple);
		UpdateMinMultiples(iMat, iMat, 0, iloss-1, multiple);

		// ---------------
		// TIMER END
		chrono.Stop();
		// ---------------
//		fprintf(its_fOut, "Computation done in %.0lf seconds\n", chrono.GetDurationInSeconds());
		ICMLOG("Computation done in %.0lf seconds" << chrono.GetDurationInSeconds() << ".");

//		fclose(its_fOut);
		// ------------------------------------------------------------------------------------------------
	}

	// ------------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------
	// ALL REMAINING MATURITIES EXCEPT THE FIRST ONE
	// ------------------------------------------------------------------------------------------------
	for (iMat=its_NbMat-2; iMat>0; iMat--)
	{
		for (iloss=its_NbLosses-1; iloss>=0; iloss--)
		{
			// ------------------------------------------------------------------------------------------------
			PrepareProductForPricing(iloss, iMat);

			// ---------------
			// TIMER BEGIN
			chrono.Start();
			// ---------------

//			its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//			fprintf(its_fOut, "Maturity:\t%u\t\tLoss:\t%u\n", iMat, iloss);
			ICMLOG("Maturity " << iMat << " - Loss " << iloss);

			// ---------------
			// PRICING CORE
			min_multiple	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
			max_multiple	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

			if (BrentFlag)
				GetOneMultiple_Brent(min_multiple, max_multiple, nbiter, multiple, f_multiple);
			else
				GetOneMultiple(min_multiple, max_multiple, nbiter, multiple);
			// ---------------

			(*its_Matrix_Multiples)(iloss, iMat)	=	multiple;
			(*its_Matrix_Triggers)(iloss, iMat)		=	f_multiple;
			
			// UPDATE
			UpdateMaxMultiples(0, iMat-1, iloss, iloss, multiple);
			UpdateMinMultiples(iMat, iMat, 1, iloss-1, multiple);

			// ---------------
			// TIMER END
			chrono.Stop();
			// ---------------
//			fprintf(its_fOut, "Computation done in %.0lf seconds\n", chrono.GetDurationInSeconds());
			ICMLOG("Computation done in %.0lf seconds" << chrono.GetDurationInSeconds() << ".");

//			fclose(its_fOut);
			// ------------------------------------------------------------------------------------------------
		}
	}


	// ------------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------
	// FIRST MATURIY, ALL LOSSES, except the LAST ONE
	// ------------------------------------------------------------------------------------------------
	iMat=0;

	for (iloss=its_NbLosses-2; iloss>=0; iloss--)
	{
		// ------------------------------------------------------------------------------------------------
		PrepareProductForPricing(iloss, iMat);

		// ---------------
		// TIMER BEGIN
		chrono.Start();
		// ---------------

//		its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//		fprintf(its_fOut, "Maturity:\t%u\t\tLoss:\t%u\n", iMat, iloss);
		ICMLOG("Maturity " << iMat << " - Loss " << iloss);

		// ---------------
		// PRICING CORE
		min_multiple	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
		max_multiple	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

		if (BrentFlag)
			GetOneMultiple_Brent(min_multiple, max_multiple, nbiter, multiple, f_multiple);
		else
			GetOneMultiple(min_multiple, max_multiple, nbiter, multiple);
		// ---------------

		(*its_Matrix_Multiples)(iloss, iMat)	=	multiple;
		(*its_Matrix_Triggers)(iloss, iMat)		=	f_multiple;
		
		// UPDATE
		UpdateMinMultiples(iMat, iMat, 0, iloss-1, multiple);

		// ---------------
		// TIMER END
		chrono.Stop();
		// ---------------
//		fprintf(its_fOut, "Computation done in %.0lf seconds\n", chrono.GetDurationInSeconds());
		ICMLOG("Computation done in %.0lf seconds" << chrono.GetDurationInSeconds() << ".");

//		fclose(its_fOut);
		// ------------------------------------------------------------------------------------------------	
	}

//	its_fOut = fopen("c:\\test\\SuperSeniorLeverage.txt", "a+");

//	fprintf(its_fOut, "----------------- Matrix of Multiples -----------------\n");
/*
	for (iloss=0; iloss<its_NbLosses; iloss++)
	{
		for (iMat=0; iMat<its_NbMat; iMat++)
			fprintf(its_fOut, "\t%.5lf\t", (*its_Matrix_Multiples)(iloss, iMat));
		fprintf(its_fOut, "\n");
	}

	fclose(its_fOut);
*/	
	return	multiple;

//return 0.0;	
}


void	ICM_SuperSeniorLeverage::ComputeMultipleRanges(bool IsDecreasingMaturity, bool IsIncreasingLoss, int iloss, int iMat, int its_NbLosses, int its_NbMat, double& min_multiple, double& max_multiple)
{
	// ONLY INCREASING LOSS HAVE BEEN IMPLEMENTED
	if (!IsIncreasingLoss)
		ICMTHROW(ERR_INVALID_DATA,"ONLY INCREASING LOSS HAVE BEEN IMPLEMENTED!");

	if (IsDecreasingMaturity)
	{
		if (iloss > 0)
		{
			if (iMat > 0)
			{
				// iloss > 0 && iMat > 0
				min_multiple	=	(*its_Matrix_Multiples)(iloss, iMat - 1);
			}
			else
			{
				// iloss > 0 && iMat = 0
				min_multiple	=	its_min_down_multiple;
			}
			max_multiple	=	(*its_Matrix_Multiples)(iloss - 1, iMat);
		}
		else
		{				
			if (iMat > 0)
				// iloss = 0 && iMat > 0
				min_multiple	=	(*its_Matrix_Multiples)(iloss, iMat - 1);
			else
				// iloss = 0 && iMat = 0
				min_multiple	=	its_min_down_multiple;
			max_multiple	=	its_max_down_multiple * (its_Vect_of_MaturitiesInYF[iMat]);
		}
	}
	else
	{
		// STARTING FROM THE LAST MATURITY WITH INCREASING LOSSES
		if (iloss > 0)
		{
			if (iMat < its_NbMat - 1)
			{
				// iloss > 0 && iMat < its_NbMat - 1
				max_multiple	=	FMIN((*its_Matrix_Multiples)(iloss, iMat + 1), (*its_Matrix_Multiples)(iloss - 1, iMat));
			}
			else
			{
				// iloss > 0 && iMat = its_NbMat - 1
				max_multiple	=	(*its_Matrix_Multiples)(iloss - 1, iMat);;
			}
			min_multiple	=	its_min_down_multiple;			
		}
		else
		{				
			if (iMat == its_NbMat - 1)
				// iloss = 0 && iMat = its_NbMat - 1
				max_multiple	=	its_max_down_multiple;
			else
				// iloss = 0 && iMat < its_NbMat - 1
				max_multiple	=	(*its_Matrix_Multiples)(iloss, iMat + 1);

			// should be enhanced
			min_multiple	=	its_min_down_multiple;

		}
	}

}



void	ICM_SuperSeniorLeverage::GetOneMultiple(double min_multiple, double max_multiple, int& nbiter, double& multiple)
{
	double	delta_value;
	double	prev_delta_value;
	double	prev_multiple;

	nbiter	=	0;

	delta_value	=	0.0;
	prev_delta_value	=	-1.0;
	prev_multiple	=	0.0;

	double	init_min_multiple;
	double	init_max_multiple;

	init_min_multiple	=	min_multiple;
	init_max_multiple	=	max_multiple;

	while ((fabs((fabs(delta_value) - its_trigger_pct)) > its_trigger_tol) && (nbiter < NB_ITER_MAX))
	{
		if (fabs(prev_delta_value - delta_value) < 1e-4)
		{
			// bad initial interval, which bound?
			if (fabs(prev_multiple - init_min_multiple) < 1e-3)
				min_multiple	=	init_min_multiple / 2.0;
			else
				max_multiple	=	init_max_multiple * 1.5;
		}
		multiple = (min_multiple + max_multiple) / 2.0;

		prev_delta_value	=	delta_value;
		GetOnePrice(multiple, delta_value);

		// Dichotomy
		if (fabs(delta_value) < its_trigger_pct)
			min_multiple	=	multiple;
		else
			max_multiple	=	multiple;

//		fprintf(its_fOut, "Iter:\t%u\t\tDelta MtM:\t%.6lf\t\tMultiple:\t%.6lf\n", nbiter, delta_value, multiple);
		ICMLOG("Iter: " << nbiter << " - Delta Value: " << delta_value << " - multiple: " << multiple);

		nbiter++;
		prev_multiple		=	multiple;
	}
}


void	ICM_SuperSeniorLeverage::PrepareProductForPricing(int iloss, int iMat)
{
	double	SubAmount;
	double	MezzAmount;
	double	current_loss;

	ARM_Date	Current_Maturity;

	ICM_Mez*	MyCloneProduct;
	
	// --------------------------------------

	bool IsDecreasingMaturity;
	bool IsIncresingLoss;

	IsDecreasingMaturity	=	false;
	IsIncresingLoss			=	true;

	// --------------------------------------

	MyCloneProduct	=	its_VectofProducts[iMat];
	current_loss	=	iloss * its_loss_step;
	
	Current_Maturity	=	its_Vect_of_Maturities[iMat];

	SubAmount		=	FMAX(its_Notional * (its_Init_Detach - current_loss), 0.0);
	MezzAmount		=	its_Notional * (FMAX(its_Init_Attach + its_Init_Detach - current_loss, 0.0) - FMAX(its_Init_Detach - current_loss, 0.0));

	MyCloneProduct->SetSubAmount(SubAmount);
	MyCloneProduct->SetMezzAmount(MezzAmount);

	MyCloneProduct->CptCashFlowDatesCredit(Current_Maturity);

	its_Pricer->SetSecurity(MyCloneProduct);

}


void	ICM_SuperSeniorLeverage::UpdateMaxMultiples(int from_Mat, int to_Mat, int from_Loss, int to_Loss, double max_multiple)
{
	int	iMat;
	int	iloss;
	double	value;

	for (iMat=from_Mat; iMat<=to_Mat; iMat++)
		for (iloss=from_Loss; iloss<=to_Loss; iloss++)
		{
			value	=	(*its_Max_Matrix_Multiples)(iloss, iMat);

			(*its_Max_Matrix_Multiples)(iloss, iMat)	=	FMIN(value, max_multiple);
		}

}


void	ICM_SuperSeniorLeverage::UpdateMinMultiples(int from_Mat, int to_Mat, int from_Loss, int to_Loss, double min_multiple)
{
	int	iMat;
	int	iloss;
	double	value;

	for (iMat=from_Mat; iMat<=to_Mat; iMat++)
		for (iloss=from_Loss; iloss<=to_Loss; iloss++)
		{
			value	=	(*its_Min_Matrix_Multiples)(iloss, iMat);

			(*its_Min_Matrix_Multiples)(iloss, iMat)	=	FMAX(value, min_multiple);
		}

}



/**  
JLA old test function
void ICM_SuperSeniorLeverage::TestCalibrationWithGuess()
{
	double	multiple;
	const ICM_DefaultCurve*		Current_DefaultCurve = NULL;
	ARM_Vector*	CurrentRates;
	
	FILE*	fOut;

	// ----------------------------------------------------------------------------
	if ((fOut = fopen("c:\\test\\CalibrationTest.txt", "a+")) == NULL) return;
	fprintf(fOut, " ----------------- Calibration Tests ----------------- \n");

	multiple	=	1.;
	fprintf(fOut, "Multiple:\t%f\n", multiple);

	for (int i=0; i<its_NbCredits; i++)
	{
		Current_DefaultCurve	=	its_MultiCurvesModel->GetDefaultCurve(i);
	
		CurrentRates	=	its_Initial_Spreads->RowAsVector(i);
		(*CurrentRates)	*=	multiple;

		Current_DefaultCurve->SetRates(CurrentRates);
	}

	PreciseChrono chrono;

	fclose(fOut);
} 
**/ 



void	ICM_SuperSeniorLeverage::GetOneMultiple_Brent(double min_multiple, double max_multiple, int& nbiter, double& multiple, double& f_multiple)
{
	double a,b,c,fa,fb,fc;
	double Newb;

	a = min_multiple;
	c = max_multiple;

	double	fx1;
	double	fx2;

	GetOnePrice(a, fx1);
	GetOnePrice(c, fx2);

	double	Target	=	its_trigger_pct;
	double	tol		=	its_trigger_tol;

	fa = fabs(fx1) - Target;
	fc = fabs(fx2) - Target;

	if (fabs(fa) < tol)
	{
		multiple = a; 
		f_multiple	=	fa;
		return;
	}
	else
  	{
		if (fabs(fc) < tol) 
		{
			multiple = c; 
			f_multiple	=	fc;
			return;
		}
  		else 
			b = (a * fc - c * fa) / (fc - fa);
		if (fabs(b - a) < 1e-10 || fabs( b - c)< 1e-10)
			b = (a + c) * .5;
	}

	for (nbiter = 1; nbiter <= NB_ITER_MAX; nbiter++)
	{
		GetOnePrice(b, fb);
		fb	= fabs(fb) - Target;
		
		if (fabs(fb) >= tol)
		{ 
			Newb = ZBrent(a, fa, b, fb, c, fc);                                             
			if (Newb == b) 
	    	{ 
				Newb = 0.5 * (a + c);
			}
			else if (((Newb < a) && (Newb < c)) || ((Newb > a) && (Newb > c)))
			{
					Newb = 0.5 * (a + c);
					if (Newb == b) throw;
			}
			if (fa * fb > 0.0) 
			{ 
				a = b; 
				fa = fb;
			}
			else 
			{ 
				c = b; 
				fc = fb;
			}
			b = Newb;
		}
		else 
	  	{ 
			multiple = b; 
			f_multiple	=	fb;
			return;
		}
	}
	
}

double ICM_SuperSeniorLeverage :: ZBrent(double a, double fa, double b, double fb, double c, double fc) const
{
	double z;

	if ((fa == 0.0) || (fb == 0.0) || (fc == 0.0))
	{
		return( 0.0 );
	}
	else if ((fa == fc) || (fb == fc) || (fb == fa))
	{
		return( b );
	}
	else
	{
		z = (fb / fa) * ((fa / fc) * ((fb / fc) - (fa / fc)) * (c - b) - (1.0 - (fb / fc)) * (b - a));
		z /= ((fa / fc) - 1.0) * ((fb / fc) - 1.0) * ((fb / fa) - 1.0);
		z += b;
		if ((fabs(z-a) < 1e-12) || (fabs(z-b) < 1e-12) )
			z = (a + b) * 0.5;
		return(z);
	}
}

void	ICM_SuperSeniorLeverage::GetOnePrice(double multiple, double& result, bool Imposed_NPV_Flag)
{
	int	i;
	const ICM_DefaultCurve*	Current_DefaultCurve = NULL;
	ARM_Vector*			CurrentRates	=	NULL;

	double	delta_value;
	double	current_value;

	// Spreads impact
	//JLA added:
	// ICM_DefaultCurve** newDefCurves = new ICM_DefaultCurve*[its_NbCredits]; 
	std::vector<const ICM_DefaultCurve*> newDefCurves (its_NbCredits); 
	for (i=0; i<its_NbCredits; i++)
	{
		/** Current_DefaultCurve	=	its_MultiCurvesModel->GetDefaultCurve(i);
		
		CurrentRates	=	its_Initial_Spreads->RowAsVector(i);
		(*CurrentRates)	*=	multiple;

		Current_DefaultCurve->SetRates(CurrentRates);

		Current_DefaultCurve->Calibrate_Stress_Test_Guess_Brent();
		**/ 
		ICM_DefaultCurve* item = dyn_clone(its_MultiCurvesModel->GetDefaultCurve(i)); 
		
		CurrentRates	=	its_Initial_Spreads->RowAsVector(i);
		(*CurrentRates)	*=	multiple;

		item->SetRates(CurrentRates);

		// newDefCurves[i]->Calibrate_Stress_Test_Guess_Brent();
		item->Calibrate();
		newDefCurves[i]=item; 
	}
	its_MultiCurvesModel->SetDefaultCurves(newDefCurves); 
	// Price
	its_Pricer->ResetRootPricer();
	((ICM_Pricer_Distrib_Smile*) its_Pricer)->clearDistribLoss();

	ENUM_CMPMETH	ThePricingMode;

	if (Imposed_NPV_Flag)
		ThePricingMode	=	qCMPPRICE;
	else
		ThePricingMode	=	its_PricingMode;

	current_value	=	((ICM_Pricer_Distrib_Smile*) its_Pricer)->ComputePrice( ThePricingMode);

	delta_value	=	current_value / its_Tranche_Size * its_leverage;

	result	=	delta_value;
}


void	ICM_SuperSeniorLeverage::UpdateAllPrices(int curr_iloss, int curr_iMat, double multiple)
{
	// from the Current Loss Distribution
	// I want to impact all the product
	
	int	iloss, iMat;
	double	current_min, current_max;

	ICM_Mez*	TheProduct;
	
	double	delta_value;
	double	current_value;

	for (iMat=curr_iMat; iMat>=0; iMat--)
	{
		for (iloss=its_NbLosses-1; iloss>=0; iloss--)
		{
			if ((iMat != curr_iMat) || (iloss <= curr_iloss))
			{
				if (!(*its_Matrix_Multiples_Flags)(iloss, iMat))
				{
					TheProduct	=	(*its_Matrix_Multiples_Products)(iloss, iMat);

					// I have to set the current Product to the Pricer
					// price it all
					its_Pricer->SetSecurity(TheProduct);
					its_Pricer->ResetRootPricer();

					current_value	=	((ICM_Pricer_Distrib_Smile*) its_Pricer)->ComputePrice( its_PricingMode);
					delta_value		=	current_value / its_Tranche_Size * its_leverage;

					if (fabs(delta_value) < its_trigger_pct)
					{
						// is it a new min for this product?
						current_min	=	(*its_Min_Matrix_Multiples)(iloss, iMat);
						if (multiple > current_min)
						{
							(*its_Min_Matrix_Control_Multiples)(iloss, iMat)	=	multiple;
							(*its_Min_Matrix_Control_Values)(iloss, iMat)		=	delta_value;
						}

					}
					else
						if (fabs(delta_value) > its_trigger_pct)
						{
							// is it a new max for this product?
							current_max	=	(*its_Max_Matrix_Multiples)(iloss, iMat);
							if (multiple < current_max)
							{
								(*its_Max_Matrix_Control_Multiples)(iloss, iMat)	=	multiple;
								(*its_Max_Matrix_Control_Values)(iloss, iMat)		=	delta_value;
							}
						}
						else
						{
							(*its_Matrix_Multiples_Flags)(iloss, iMat)	=	true;
							(*its_Matrix_Multiples)(iloss, iMat)		=	multiple;
						}
				}
			}
		}
	}
}




void	ICM_SuperSeniorLeverage::ComputeTriggersFromMultiples()
{
	// SCAN FLAGS MATRIX
	int	iMat, iloss;
	double	TheFlag;

	// should i check the dimensions
	int	NbRows;
	int NbCols;
/*
	NbRows	=	its_Matrix_Triggers->Getnbrows();
	NbCols	=	its_Matrix_Triggers->Getnbcols();

	if ((NbRows < its_NbLosses) || (NbCols < its_NbMat))
		// ERROR
		ICMTHROW(ERR_INVALID_DATA,"BAD size for Triggers Matrix!");
*/
	NbRows	=	its_Matrix_MultiplesInput->Getnbrows();
	NbCols	=	its_Matrix_MultiplesInput->Getnbcols();

	if ((NbRows < its_NbLosses) || (NbCols < its_NbMat))
		// ERROR
		ICMTHROW(ERR_INVALID_DATA,"BAD size for Multiples Input Matrix!");

	double	Current_MtM;
	double	Current_Multiple;
	double	Implied_Leverage;

	for (iMat=0; iMat<its_NbMat; iMat++)
	{
		for (iloss=0; iloss<its_NbLosses; iloss++)
		{
			TheFlag	=	(*its_Matrix_Flags)(iloss, iMat);

			if (TheFlag == 1.0)
			{
				// Get the MULTIPLE
				Current_Multiple	=	(*its_Matrix_MultiplesInput)(iloss, iMat);

				// Modify Spreads
				// DO THE PRICE
				GetOnePrice(Current_Multiple, Current_MtM); //, true);		// IMPOSE NPV
				
				// Get the implied Trigger level
				Implied_Leverage	=	fabs(Current_MtM); // its_leverage * Current_MtM / its_Tranche_Size;
				
				// FILL THE MATRIX
				(*its_Matrix_Triggers)(iloss, iMat)	=	Implied_Leverage;
			}
		}
	}
}