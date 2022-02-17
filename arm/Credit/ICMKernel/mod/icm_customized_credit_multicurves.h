/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CUSTOMIZED_CREDIT_MULTICURVES.H
	PROJECT:	MOD
	
	DESCRIPTION:	this class provides a basic Multi Curves Market Data Container


   -----------------------------------------------------------------
   
	ICM CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _ICM_CUSTOMIZED_CREDIT_MULTICURVES_H_
# define _ICM_CUSTOMIZED_CREDIT_MULTICURVES_H_

# define	_SIZE_CREDIT_LABEL_		60

#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\util\icm_matrix.h"

#include "ICMKernel\cair\types.h"
#include "ICMKernel\util\icm_utils.h"


class ICM_Customized_Credit_MultiCurves : public ICM_ModelMultiCurves        
{        

	private: 

		ICM_Matrix<ARM_Vector>* its_Matrix_Data_Description;		
		
		ICM_Matrix<ARM_Vector>* its_Matrix_CDO_Square_Data;
		ICM_Matrix<ARM_Vector>* its_Matrix_CDO_Square_Parameters;

		ICM_Matrix<ARM_Vector>* its_Matrix_Market_Parameters;
		
		ICM_QMatrix<double>*	its_CreditDataSpreads;			// Credit Data Spreads
		
		// --------------------------------------------------------------------------
		// Credit Data from Interface
		// --------------------------------------------------------------------------
		
		//int		its_NbCredits;
		//int		its_NbMaturities;
		vector<string> its_CreditsLabels;
		vector<string> its_Maturities;
		//char**	its_CreditsLabelsAsChar;	// labels
		//char**	its_Maturities_AsChar;

		// CDO Square
		int		its_CDO_Square_Nb_Underlyings;

		// useful
		double	its_ValDateAsDouble;

	public:

		// --------------------------------------------------------------------------------
		// CONSTRUCTOR
		// --------------------------------------------------------------------------------
		
		ICM_Customized_Credit_MultiCurves()
		{
			Init();
		}

		ICM_Customized_Credit_MultiCurves(
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
							ICM_Correlation* Correlation = NULL)
		{
			Init();
			Set(ZC_Curve, Data_Description_matrix, CDO_Square_Data_matrix, CDO_Square_Parameters_matrix, Market_Parameters, spreads, Labels, 
							 Maturities, Correlation);
		}

	
		void	ICM_Customized_Credit_MultiCurves::Set(
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
							ICM_Correlation* Correlation = NULL
							);


		ICM_Customized_Credit_MultiCurves(ICM_ModelMultiCurves* data)
		{
			Init();
			Set(data);
		}

		void	ICM_Customized_Credit_MultiCurves::Set(ICM_ModelMultiCurves* data)
		{};

		// --------------------------------------------------------------------------------
		// DESTRUCTOR
		// --------------------------------------------------------------------------------

		virtual ~ICM_Customized_Credit_MultiCurves()
		{
			if (its_Matrix_Data_Description)
				delete its_Matrix_Data_Description;
			its_Matrix_Data_Description = NULL;

			if (its_Matrix_Market_Parameters)
				delete its_Matrix_Market_Parameters;
			its_Matrix_Market_Parameters = NULL;

			if (its_Matrix_CDO_Square_Data)
				delete its_Matrix_CDO_Square_Data;
			its_Matrix_CDO_Square_Data = NULL;

			if (its_Matrix_CDO_Square_Parameters)
				delete its_Matrix_CDO_Square_Parameters;
			its_Matrix_CDO_Square_Parameters = NULL;


			if (its_CreditDataSpreads)
				delete its_CreditDataSpreads;
			its_CreditDataSpreads	=	NULL;

			/*if (its_CreditsLabelsAsChar)
				FreePointerTabChar(its_CreditsLabelsAsChar,	its_NbCredits);

			if (its_Maturities_AsChar)
				FreePointerTabChar(its_Maturities_AsChar,	its_NbMaturities);	*/	

			its_Categories.clear();
			its_Currencies.clear();
			its_Accrueds.clear();
			its_Recoveries.clear();
			its_Notionals.clear();
			its_Losses.clear();
			its_DefaultDates.clear();
			its_AmortizationDates.clear();

		}

		// --------------------------------------------------------------------------------
		// UTILS
		// --------------------------------------------------------------------------------

		void BitwiseCopy(const ARM_Object* src);
		void Copy(const ARM_Object* src);
		ARM_Object* Clone(void);
		void View(char* id, FILE* ficOut);
		void Init();
				

	// Methodes Set and Get ----------------------------

	void Set_Matrix_Data_Description(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_Data_Description)
			delete its_Matrix_Data_Description;
		its_Matrix_Data_Description = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_Data_Description(void)	{return its_Matrix_Data_Description;}


	void Set_Matrix_Market_Parameters(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_Market_Parameters)
			delete its_Matrix_Market_Parameters;
		its_Matrix_Market_Parameters = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_Market_Parameters(void)	{return its_Matrix_Market_Parameters;}


	void Set_Matrix_CDO_Square_Data(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_CDO_Square_Data)
			delete its_Matrix_CDO_Square_Data;
		its_Matrix_CDO_Square_Data = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_CDO_Square_Data(void)	{return its_Matrix_CDO_Square_Data;}


	void Set_Matrix_CDO_Square_Parameters(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_CDO_Square_Parameters)
			delete its_Matrix_CDO_Square_Parameters;
		its_Matrix_CDO_Square_Parameters = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_CDO_Square_Parameters(void)	{return its_Matrix_CDO_Square_Parameters;}

	// ----------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------
	// DATA
	// ----------------------------------------------------------------------------------------------

	public:

		void	SetCreditDataParameters(ICM_Matrix<ARM_Vector>* parameters);
		void	SetCreditDataDescription(ICM_Matrix<ARM_Vector>* parameters);

		void	FillCreditDataParameters()
		{
			if (its_Matrix_Market_Parameters)
				SetCreditDataParameters(its_Matrix_Market_Parameters);
		}

		void	FillCreditDataDescription()
		{
			if (its_Matrix_Data_Description)
				SetCreditDataDescription(its_Matrix_Data_Description);
		}

		void	SetCreditCDO_Square_Data(ICM_Matrix<ARM_Vector>* parameters);
		void	SetCreditCDO_Square_Parameters(ICM_Matrix<ARM_Vector>* parameters);

		void	FillCreditCDO_Square_Data()
		{
			if (its_Matrix_CDO_Square_Data)
				SetCreditCDO_Square_Data(its_Matrix_CDO_Square_Data);
		}

		void	FillCreditCDO_Square_Parameters()
		{
			if (its_Matrix_CDO_Square_Parameters)
				SetCreditCDO_Square_Parameters(its_Matrix_CDO_Square_Parameters);
		}


		// -----------------------------------------------------------
	// -----------------------------------------------------------
	// CLASS USE
	// -----------------------------------------------------------
	// -----------------------------------------------------------




	public:

		int	Get_NbCredits()	const {return its_CreditsLabels.size();}

		void Set_CreditDataSpreads(ICM_QMatrix<double>* value) 
		{ 
			if (its_CreditDataSpreads)
				delete its_CreditDataSpreads;
			its_CreditDataSpreads = value; 
		}

		ICM_QMatrix<double>* Get_CreditDataSpreads(void) { return its_CreditDataSpreads;}

		void	Get_Categories(vector<CreditCategory>& data)	{data = its_Categories;}
		void	Set_Categories(const vector<CreditCategory>& data)	{its_Categories = data;}

		void	Get_Currencies(vector<CurrencyName>& data)	{data = its_Currencies;}
		void	Set_Currencies(const vector<CurrencyName>& data)	{its_Currencies = data;}

		void	Get_Accrueds(IntVector& data)	{data = its_Accrueds;}
		void	Set_Accrueds(const IntVector& data)	{its_Accrueds = data;}

		void	Get_Recoveries(DoubleVector& data)	{data = its_Recoveries;}
		void	Set_Recoveries(const DoubleVector& data)	{its_Recoveries = data;}

		void	Get_Notionals(DoubleVector& data)	{data = its_Notionals;}
		void	Set_Notionals(const DoubleVector& data)	{its_Notionals= data;}
		
		void	Get_Losses(DoubleVector& data)	{data = its_Losses;}
		void	Set_Losses(const DoubleVector& data)	{its_Losses= data;}

		void	Get_DefaultDates(DateVector& data)	{data = its_DefaultDates;}
		void	Set_DefaultDates(const DateVector& data)	{its_DefaultDates = data;}

		void	Get_AmortizationDates(DateVector& data)	{data = its_AmortizationDates;}
		void	Set_AmortizationDates(const DateVector& data)	{its_AmortizationDates = data;}

		int Get_NbMaturities() const	{return its_Maturities.size();}
		
		void	Get_RollDateFlag(bool& data)	{data = its_RollDateFlag;}
		void	Get_CDO_Square_Nb_Underlyings(int& data)	{data = its_CDO_Square_Nb_Underlyings;}


		const vector<string>&	Get_CreditsLabels() const ;
		const vector<string>&	Get_CreditDataMaturities() const ;

		
		void	Set_CreditsLabels(const vector<string>& vLabels);
		void	Set_CreditDataMaturities(const vector<string>& vMat);

		// -----------------------------------------------------------------------------

		void	Get_BumpSpread(double& data)	{data = its_BumpSpread;}
		void	Set_BumpSpread(double data)	{its_BumpSpread	=	data;}

		void	Get_BumpType(BumpType& data)	{data = its_BumpSpread_Type;}
		void	Set_BumpType(BumpType data)	{its_BumpSpread_Type	=	data;}

		void	Get_BumpRecovery(double& data)	{data = its_BumpRecovery;}
		void	Set_BumpRecovery(double data)	{its_BumpRecovery	=	data;}

		void	Get_BumpCorrelation(double& data)	{data = its_BumpCorrelation;}
		void	Set_BumpCorrelation(double data)	{its_BumpCorrelation	=	data;}

	private:

		bool			its_RollDateFlag;			// use Roll Dates or Not

		// Credit Data Description
		vector<CreditCategory>			its_Categories;			// Sectorial or Geographical Category
		vector<CurrencyName>			its_Currencies;			// Currencies

		IntVector			its_Accrueds;			// Accrued Flags
		DoubleVector		its_Recoveries;			// Recoveries
		DoubleVector		its_Notionals;			// Notionals
		DoubleVector		its_Losses;				// Losses
		DateVector			its_DefaultDates;		// Imposed Default Dates
		DateVector			its_AmortizationDates;		// Imposed Amortization Dates

		// surely not the better place to be
		double			its_BumpSpread;			// Bump Spread Value
		BumpType		its_BumpSpread_Type;	// Enum Type Add or Mult
		double			its_BumpRecovery;			// Bump Recovrery Value (in %: 0.10 for 10ù)
		double			its_BumpCorrelation;		// Bump Correlation Value (in %: 0.10 for 10ù)

};

#endif
