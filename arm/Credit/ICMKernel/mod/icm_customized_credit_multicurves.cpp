#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CUSTOMIZED_CREDIT_MULTICURVES.CPP
	PROJECT:	MOD
	
	DESCRIPTION:	this class provides a basic Multi Curves Market Data Container


   -----------------------------------------------------------------
   
	ICM CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICMKernel\mod\icm_customized_credit_multicurves.h"
#include "ICMKernel\glob\icm_correlation.h"


void ICM_Customized_Credit_MultiCurves::Init(void)
{
	ICM_ModelMultiCurves::Init();	// is it necessary?

	SetName(ICM_CUSTOMIZED_CREDIT_MULTICURVES);

	its_Matrix_Data_Description	=	NULL;
	its_Matrix_Market_Parameters	=	NULL;
	its_Matrix_CDO_Square_Data	=	NULL;
	its_Matrix_CDO_Square_Parameters	=	NULL;

	its_CreditDataSpreads	=	NULL;

	//its_CreditsLabelsAsChar	=	NULL;
	//its_Maturities_AsChar	=	NULL;

	its_CDO_Square_Nb_Underlyings	=	0;

	its_ValDateAsDouble	=	0.0;

	its_Categories.clear();
	its_Currencies.clear();
	its_Accrueds.clear();
	its_Recoveries.clear();
	its_Notionals.clear();
	its_Losses.clear();
	its_DefaultDates.clear();
	its_AmortizationDates.clear();

	its_RollDateFlag		=	true;
}



void ICM_Customized_Credit_MultiCurves::Set(
							ARM_ZeroCurve*	ZC_Curve,
							ICM_Matrix<ARM_Vector>* Data_Description_matrix,
							ICM_Matrix<ARM_Vector>* CDO_Square_Data_matrix,
							ICM_Matrix<ARM_Vector>* CDO_Square_Parameters_matrix,
							ICM_Matrix<ARM_Vector>* Market_Parameters,
							ICM_QMatrix<double>*	spreads,
							const vector<string>& Labels,
							//char**	Labels_As_Char,
							const vector<string>& Maturities,
							//char**	Maturities_As_Char,
							ICM_Correlation* Correlation
							)
{
	SetName(ICM_CUSTOMIZED_CREDIT_MULTICURVES);

	// ---------------------------------------------------------------------------	
	// ---------------------------------------------------------------------------

	// ---------------------------------------------------------------------------
	if (ZC_Curve)
		SetZeroCurve(ZC_Curve);
	// ---------------------------------------------------------------------------

	// ---------------------------------------------------------------------------
	if (its_Matrix_Data_Description)
		delete its_Matrix_Data_Description;

	if (Data_Description_matrix)
		its_Matrix_Data_Description = (ICM_Matrix<ARM_Vector>*) Data_Description_matrix->Clone();
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	if (its_Matrix_Market_Parameters)
		delete its_Matrix_Market_Parameters;

	if (Market_Parameters)
		its_Matrix_Market_Parameters = (ICM_Matrix<ARM_Vector>*) Market_Parameters->Clone();
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	if (its_Matrix_CDO_Square_Data)
		delete its_Matrix_CDO_Square_Data;

	if (CDO_Square_Data_matrix)
		its_Matrix_CDO_Square_Data = (ICM_Matrix<ARM_Vector>*) CDO_Square_Data_matrix->Clone();
	// ---------------------------------------------------------------------------


	// ---------------------------------------------------------------------------
	if (its_Matrix_CDO_Square_Parameters)
		delete its_Matrix_CDO_Square_Parameters;

	if (CDO_Square_Parameters_matrix)
		its_Matrix_CDO_Square_Parameters = (ICM_Matrix<ARM_Vector>*) CDO_Square_Parameters_matrix->Clone();
	// ---------------------------------------------------------------------------

	// ---------------------------------------------------------------------------
	// ---------------------------------------------------------------------------
	if (its_CreditDataSpreads)
		delete its_CreditDataSpreads;

	if (spreads)
		its_CreditDataSpreads = (ICM_QMatrix<double>*) spreads->Clone();
	// ---------------------------------------------------------------------------

	its_CreditsLabels = Labels;
	its_Maturities = Maturities;

	if (Correlation)
		SetCorrelation((ICM_Correlation*) Correlation); // ?? ->Clone()

}


void ICM_Customized_Credit_MultiCurves::BitwiseCopy(const ARM_Object* src)
{
    ICM_Customized_Credit_MultiCurves* gen = (ICM_Customized_Credit_MultiCurves*) src;

	if (gen->its_Matrix_Data_Description)
	{
		if (its_Matrix_Data_Description)
			delete its_Matrix_Data_Description;
		its_Matrix_Data_Description = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_Data_Description->Clone();
	}

	if (gen->its_Matrix_Market_Parameters)
	{
		if (its_Matrix_Market_Parameters)
			delete its_Matrix_Market_Parameters;
		its_Matrix_Market_Parameters = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_Market_Parameters->Clone();
	}

	if (gen->its_Matrix_CDO_Square_Data)
	{
		if (its_Matrix_CDO_Square_Data)
			delete its_Matrix_CDO_Square_Data;
		its_Matrix_CDO_Square_Data = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_CDO_Square_Data->Clone();
	}

	if (gen->its_Matrix_CDO_Square_Parameters)
	{
		if (its_Matrix_CDO_Square_Parameters)
			delete its_Matrix_CDO_Square_Parameters;
		its_Matrix_CDO_Square_Parameters = (ICM_Matrix<ARM_Vector>*) gen->its_Matrix_CDO_Square_Parameters->Clone();
	}

	if (gen->its_CreditDataSpreads)
	{
		if (its_CreditDataSpreads)
			delete its_CreditDataSpreads;
		its_CreditDataSpreads = (ICM_QMatrix<double>*) gen->its_CreditDataSpreads->Clone();
	}

	/* Correlation
	if (gen->itsCorrelation)
		itsCorrelation = (ICM_Correlation*) gen->itsCorrelation->Clone();*/
	///  ??? SetCorrelation(src->GetCorrelation()->Clone());

	// DATA
	its_CDO_Square_Nb_Underlyings			=	gen->its_CDO_Square_Nb_Underlyings;

	its_ValDateAsDouble		=	gen->its_ValDateAsDouble;
	its_RollDateFlag		=	gen->its_RollDateFlag;
	
	its_BumpSpread			=	gen->its_BumpSpread;
	its_BumpSpread_Type		=	gen->its_BumpSpread_Type;
	its_BumpRecovery		=	gen->its_BumpRecovery;
	its_BumpCorrelation		=	gen->its_BumpCorrelation;

	// CHARACTERS
	its_CreditsLabels =	gen->its_CreditsLabels;
	its_Maturities	=	gen->its_Maturities;

	its_Categories		=	gen->its_Categories;
	its_Currencies		=	gen->its_Currencies;
	its_Accrueds		=	gen->its_Accrueds;
	its_Recoveries		=	gen->its_Recoveries;
	its_Notionals		=	gen->its_Notionals;
	its_Losses			=	gen->its_Losses;
	its_DefaultDates	=	gen->its_DefaultDates;
	its_AmortizationDates		=	gen->its_AmortizationDates;


}


void ICM_Customized_Credit_MultiCurves::Copy(const ARM_Object* src)
{
	ICM_ModelMultiCurves::Copy(src);
	BitwiseCopy(src);
}


ARM_Object* ICM_Customized_Credit_MultiCurves::Clone(void)
{
     ICM_Customized_Credit_MultiCurves* theClone = new ICM_Customized_Credit_MultiCurves();

     theClone->Copy(this);
 
     return(theClone);
}

// *************************************************************
// View Matrix 
// *************************************************************

void ICM_Customized_Credit_MultiCurves::View(char* id, FILE* ficOut)
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

	int	i;

	// -----------------------------------------------------------------------------
	// Nb Credits
	// -----------------------------------------------------------------------------
	fprintf(fOut, "Number of Credits:\t%u\nNumber of CDS Maturities:\t%u\n\n", its_CreditsLabels.size(), its_Maturities.size());

	// ------------------------------------------------------------------------------
	// LABELS
	// ------------------------------------------------------------------------------
	fprintf(fOut, "Labels:\n");
	for (i=0; i<its_CreditsLabels.size(); i++)
		fprintf(fOut,"\t%s\n", its_CreditsLabels[i].c_str());
		
	fprintf(fOut, "Maturities:\n");
	for (i=0; i<its_Maturities.size(); i++)
		fprintf(fOut,"\t%s\t\t", its_Maturities[i].c_str());
	fprintf(fOut,"\n\n");


	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Data Description ----------------- \n");

	if (Get_Matrix_Data_Description())
		Get_Matrix_Data_Description()->View(id, fOut);

	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Data Spreads ----------------- \n");

	if (Get_CreditDataSpreads())
		Get_CreditDataSpreads()->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- Market Parameters ----------------- \n");

	if (Get_Matrix_Market_Parameters())
		Get_Matrix_Market_Parameters()->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- CDO Square Data ----------------- \n");

	if (Get_Matrix_CDO_Square_Data())
		Get_Matrix_CDO_Square_Data()->View(id, fOut);

	fprintf(fOut, "\t\t\t ----------------- CDO Square Parameters ----------------- \n");

	if (Get_Matrix_CDO_Square_Parameters())
		Get_Matrix_CDO_Square_Parameters()->View(id, fOut);

	fprintf(fOut,"\t\t\t ----------------- Correlation Matrix ----------------- \n\n");

	ICM_Correlation* matrix = GetCorrelation();
	matrix->View(id,fOut);
	
	if (ficOut == NULL)
	{
		fclose(fOut);
	}
}


// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------
//	GET DEFAULT LEG INFORMATION
// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------

void	ICM_Customized_Credit_MultiCurves::Set_CreditsLabels(const vector<string> & vLabels)
{
	its_CreditsLabels = vLabels;
}

void	ICM_Customized_Credit_MultiCurves::Set_CreditDataMaturities(const vector<string> & vMat)
{
	its_Maturities = vMat;
}
//--------------------------------------------
//	CREDIT DATA PARAMETERS
//--------------------------------------------

void	ICM_Customized_Credit_MultiCurves::SetCreditDataParameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	int	TheValue;

	ARM_Vector* THE_VECT=NULL;

	// ------------------------------------------------------
	// FIRST
	// ------------------------------------------------------
	// DOUBLE TYPE: ValDate from Excel
	
	THE_VECT	= parameters->GetColVect("EXCEL_VALDATE");
	if (!THE_VECT)
		its_ValDateAsDouble	=	0.0;
	else
		its_ValDateAsDouble = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// ROLL_DATE:	BOOLEAN

	THE_VECT	= parameters->GetColVect("ROLL_DATE");
	if (!THE_VECT)
		its_RollDateFlag	=	true;
	else
		its_RollDateFlag = (THE_VECT->Elt(0) ? false : true);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQUARE_NB_UNDERLYINGS:	INT

	THE_VECT	= parameters->GetColVect("CDO_SQUARE_NB_UNDERLYINGS");
	if (!THE_VECT)
		its_CDO_Square_Nb_Underlyings	=	0;
	else
		its_CDO_Square_Nb_Underlyings = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_SPREAD");
	if (!THE_VECT)
		its_BumpSpread	=	10.0;
	else
		its_BumpSpread = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD_TYPE:	int

	THE_VECT	= parameters->GetColVect("BUMP_SPREAD_TYPE");
	if (!THE_VECT)
		its_BumpSpread_Type	=	BT_ADD;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < BT_ADD) || (TheValue > BT_MULT)) 
			its_BumpSpread_Type	=	BT_ADD;
		else	
			its_BumpSpread_Type = (BumpType) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_RECOVERY:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_RECOVERY");
	if (!THE_VECT)
		its_BumpRecovery	=	0.1;
	else
		its_BumpRecovery = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_CORRELATION:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_CORRELATION");
	if (!THE_VECT)
		its_BumpCorrelation	=	0.1;
	else
		its_BumpCorrelation = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

}


//--------------------------------------------
//	DATA DESCRIPTION
//--------------------------------------------

void	ICM_Customized_Credit_MultiCurves::SetCreditDataDescription(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size, tmp_size;

	// ------------------------------------------------------
	// CATEGORY:	VECTOR OF INT --> TO MATCH WITH ENUM TYPE

	THE_VECT	= parameters->GetColVect("CATEGORY");
	if (!THE_VECT)
		size	=	0;
	else
		size =	THE_VECT->GetSize();	

	its_Categories.resize(size);

	for (i=0;i<size;i++)
		its_Categories[i]	=	(CreditCategory) (int) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CURRENCY:	VECTOR OF INT --> TO MATCH WITH ENUM TYPE

	THE_VECT	= parameters->GetColVect("CURRENCY");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;

	its_Currencies.resize(size);

	for (i=0;i<size;i++)
		its_Currencies[i]	=	(CurrencyName) (int) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// ACCRUED:	VECTOR OF INT

	THE_VECT	= parameters->GetColVect("ACCRUED");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_Accrueds.resize(size);

	for (i=0;i<size;i++)
		its_Accrueds[i]	=	(int) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// RECOVERY:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("RECOVERY");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_Recoveries.resize(size);

	for (i=0;i<size;i++)
		its_Recoveries[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NOTIONAL:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("NOTIONAL");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
		
	its_Notionals.resize(size);

	for (i=0;i<size;i++)
		its_Notionals[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_Losses.resize(size);

	for (i=0;i<size;i++)
		its_Losses[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DEFAULT_DATE:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("DEFAULT_DATE");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_DefaultDates.resize(size);

	for (i=0;i<size;i++)
		its_DefaultDates[i]	=	(RelativeDate) ((*THE_VECT)[i] - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// AMORTIZATION_DATE:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("AMORTIZATION_DATE");
	if (!THE_VECT)
		tmp_size	=	0;
	else
		tmp_size	=	THE_VECT->GetSize();
	
	if (size !=	tmp_size)
		tmp_size	=	0;
	
	its_AmortizationDates.resize(size);

	for (i=0;i<size;i++)
		its_AmortizationDates[i]	=	(RelativeDate) ((*THE_VECT)[i] - its_ValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------
}


const vector<string>& ICM_Customized_Credit_MultiCurves::Get_CreditsLabels() const 
{
	return	its_CreditsLabels;
}


const vector<string>&	ICM_Customized_Credit_MultiCurves::Get_CreditDataMaturities() const
{
	return its_Maturities;
}

