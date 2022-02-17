
//////////////////////////////////////////////////////////////////////
// ICM_GEN.h: interface for the ICM_Customized_CDO class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ICM_CUSTOMIZED_CDO_H_)
#define _ICM_CUSTOMIZED_CDO_H_


#include "ICMKernel\inst\icm_ftd.h"
#include "ICMKernel\util\icm_matrix.h"

#include "ICMKernel\cair\types.h"

/*********************************************************************************/
/*! \class  ICM_Customized_CDO icm_gen.h "icm_gen.h"
 *  \author	Laurent Jacquel 
 *	\version 1.0
 *	\date   December 2005
 *	\brief  Pricing <B> Customized CDO Cash Flow </B> */
/***********************************************************************************/


class ICM_Customized_CDO : public ICM_Ftd  
{
private: 

	ICM_Matrix<ARM_Vector>* its_Matrix_PL_Data;
	ICM_Matrix<ARM_Vector>* its_Matrix_DL_Data;
	ICM_Matrix<ARM_Vector>* its_Matrix_Parameters_Data;
/*
	ICM_Parameters*		its_Default_Parameters,
	ICM_Parameters*		its_Premium_Parameters,
	ICM_Parameters*		its_Data_Parameters,
*/
public:

	void Init();


	ICM_Customized_CDO()
	{
		Init();
	}	


	ICM_Customized_CDO(
							ICM_Parameters*	Default_Parameters,
							ICM_Parameters* Premium_Parameters,
							ICM_Parameters* Data_Parameters,
							const int& NbIssuers =-1,
							char**	IssuersLabels = NULL,
							double*	IssuersNotionals = NULL,
							const std::string& ccy = "EUR", 
							ARM_Date* EffectiveDate = NULL
							);

	ICM_Customized_CDO(
							ICM_Matrix<ARM_Vector>* PL_matrix,
							ICM_Matrix<ARM_Vector>* DL_matrix,
							ICM_Matrix<ARM_Vector>* Pricing_matrix,
							const int& NbIssuers =-1,
							char**	IssuersLabels = NULL,
							double*	IssuersNotionals = NULL,
							const std::string& ccy = "EUR", 
							ARM_Date* EffectiveDate = NULL
							);

	void Set(
				ICM_Matrix<ARM_Vector>* PL_matrix,
				ICM_Matrix<ARM_Vector>* DL_matrix,
				ICM_Matrix<ARM_Vector>* Pricing_matrix,
				const int& NbIssuers =-1,
				char**	IssuersLabels = NULL,
				double*	IssuersNotionals = NULL,
				const std::string& ccy = "EUR", 
				ARM_Date* EffectiveDate = NULL);

	virtual ~ICM_Customized_CDO()
	{

		if (its_Matrix_PL_Data)
			delete its_Matrix_PL_Data;
		its_Matrix_PL_Data = NULL;

		if (its_Matrix_DL_Data)
			delete its_Matrix_DL_Data;
		its_Matrix_DL_Data = NULL;

		if (its_Matrix_Parameters_Data)
			delete its_Matrix_Parameters_Data;
		its_Matrix_Parameters_Data = NULL;

/*
		if (its_Default_Parameters)
			delete its_Default_Parameters;
		its_Default_Parameters = NULL;

		if (its_Premium_Parameters)
			delete its_Premium_Parameters;
		its_Premium_Parameters = NULL;

		if (its_Data_Parameters)
			delete its_Data_Parameters;
		its_Data_Parameters = NULL;
*/		

		// Premium Leg Data
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
	}


	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	void View(char* id, FILE* ficOut);

	// Methodes Set and Get ----------------------------

	void Set_Matrix_PL_Data(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_PL_Data)
			delete its_Matrix_PL_Data;
		its_Matrix_PL_Data = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_PL_Data(void)	{return its_Matrix_PL_Data;}

	void Set_Matrix_DL_Data(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_DL_Data)
			delete its_Matrix_DL_Data;
		its_Matrix_DL_Data = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_DL_Data(void)	{return its_Matrix_DL_Data;}

	void Set_Matrix_Parameters_Data(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (its_Matrix_Parameters_Data)
			delete its_Matrix_Parameters_Data;
		its_Matrix_Parameters_Data = matrix;
	}

	ICM_Matrix<ARM_Vector>* Get_Matrix_Parameters_Data(void)	{return its_Matrix_Parameters_Data;}

	// ----------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------
	// DATA
	// ----------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------

	public:

		void	Get_PL_NbFlows(int& data)	{data = its_PL_NbFlows;}
		void	Set_PL_NbFlows(const int& data)	{its_PL_NbFlows = data;}

		void	Get_ValDate(ARM_Date& data)	{data = its_ValDate;}
		void	Get_ValDateAsDouble(double& data)	{data = its_ValDateAsDouble;}

	private:

		ARM_Date		its_ValDate;		// AsOfDate
		double			its_ValDateAsDouble;

		int				its_PL_NbFlows;

		// -----------------------------------------------------------------
		// Default Leg Data
		// -----------------------------------------------------------------

	public:

		void	Get_DL_CreditWindowLow(RelativeDate& data)	{data = its_DL_CreditWindowLow;}
		void	Set_DL_CreditWindowLow(const RelativeDate& data)	{its_DL_CreditWindowLow = data;}

		void	Get_DL_CreditWindowUp(RelativeDate& data)	{data = its_DL_CreditWindowUp;}
		void	Set_DL_CreditWindowUp(const RelativeDate& data)	{its_DL_CreditWindowUp = data;}

		void	Get_DL_PaymentDate(RelativeDate& data)	{data = its_DL_PaymentDate;}
		void	Set_DL_PaymentDate(const RelativeDate& data)	{its_DL_PaymentDate = data;}

		void	Get_DL_LossMin(double& data)	{data = its_DL_LossMin;}
		void	Set_DL_LossMin(const double& data)	{its_DL_LossMin= data;}
		
		void	Get_DL_LossMax(double& data)	{data = its_DL_LossMax;}
		void	Set_DL_LossMax(const double& data)	{its_DL_LossMax= data;}

		void	Get_DL_PaymentLag(int& data)	{data = its_DL_PaymentLag;}
		void	Set_DL_PaymentLag(const int& data)	{its_DL_PaymentLag = data;}

		void	Get_DL_PaymentType(CreditEventPayment& data)	{data = its_DL_PaymentType;}

	private:
		

		RelativeDate	its_DL_CreditWindowLow;			// Credit Observation Window Low
		RelativeDate	its_DL_CreditWindowUp;			// Credit Observation Window High
		double			its_DL_LossMin;					// Credit Loss Min
		double			its_DL_LossMax;					// Credit Loss Max
		int				its_DL_NbDefMin;					// Credit Nb Def Min
		int				its_DL_NbDefMax;					// Credit Nb Def Max
		CreditEventPayment	its_DL_PaymentType;				// Credit Payment (Default, Maturity, etc.)
		RelativeDate	its_DL_PaymentDate;					// if required the Date
		int				its_DL_PaymentLag;					// the lag
		// -----------------------------------------------------------------


		// -----------------------------------------------------------------
		// Premium Leg Data
		// -----------------------------------------------------------------
		
	public:

		void	Get_PL_CreditWindowLows(DateVector& data)	{data = its_PL_CreditWindowLows;}
		void	Set_PL_CreditWindowLows(const DateVector& data)	{its_PL_CreditWindowLows = data;}

		void	Get_PL_CreditWindowUps(DateVector& data)	{data = its_PL_CreditWindowUps;}
		void	Set_PL_CreditWindowUps(const DateVector& data)	{its_PL_CreditWindowUps = data;}

		void	Get_PL_StartDates(DateVector& data)	{data = its_PL_StartDates;}
		void	Set_PL_StartDates(const DateVector& data)	{its_PL_StartDates = data;}

		void	Get_PL_EndDates(DateVector& data)	{data = its_PL_EndDates;}
		void	Set_PL_EndDates(const DateVector& data)	{its_PL_EndDates = data;}

		void	Get_PL_PaymentDates(DateVector& data)	{data = its_PL_PaymentDates;}
		void	Set_PL_PaymentDates(const DateVector& data)	{its_PL_PaymentDates = data;}

		void	Get_PL_LossMins(DoubleVector& data)	{data = its_PL_LossMins;}
		void	Set_PL_LossMins(const DoubleVector& data)	{its_PL_LossMins= data;}
		
		void	Get_PL_LossMaxs(DoubleVector& data)	{data = its_PL_LossMaxs;}
		void	Set_PL_LossMaxs(const DoubleVector& data)	{its_PL_LossMaxs= data;}

		void	Get_PL_Ratios(DoubleVector& data)	{data = its_PL_Ratios;}
		void	Set_PL_Ratios(const DoubleVector& data)	{its_PL_Ratios = data;}

		void	Get_PL_Notios(DoubleVector& data)	{data = its_PL_Notios;}
		void	Set_PL_Notios(const DoubleVector& data)	{its_PL_Notios = data;}

		void	Get_PL_Spreads(DoubleVector& data)	{data = its_PL_Spreads;}
		void	Set_PL_Spreads(const DoubleVector& data)	{its_PL_Spreads = data;}
		
		void	Get_PL_CreditFlags(vector<CreditPremiumLegType>& data)	{data = its_PL_CreditFlags;}
		void	Set_PL_CreditFlags(const vector<CreditPremiumLegType>& data)	{its_PL_CreditFlags = data;}

		void	Get_PL_CreditSpreadCaps(DoubleVector& data)	{data = its_PL_CreditSpreadCaps;}
		void	Set_PL_CreditSpreadCaps(const DoubleVector& data)	{its_PL_CreditSpreadCaps = data;}

		void	Get_PL_Redemptions(DoubleVector& data)	{data = its_PL_Redemptions;}
		void	Set_PL_Redemptions(const DoubleVector& data)	{its_PL_Redemptions = data;}

		void	Get_PL_PaymentType(CreditEventPayment& data)	{data = its_PL_PaymentType;}

	private:

		// first Vectors
		DateVector		its_PL_CreditWindowLows;			// Credit Observation Window Low
		DateVector		its_PL_CreditWindowUps;			// Credit Observation Window Low
		DateVector		its_PL_StartDates;					// Start Dates for Period
		DateVector		its_PL_EndDates;					// End Dates for Period
		DateVector		its_PL_PaymentDates;				// Payment Dates for Period

		DoubleVector	its_PL_LossMins;
		DoubleVector	its_PL_LossMaxs;
		DoubleVector	its_PL_NbDefMins;
		DoubleVector	its_PL_NbDefMaxs;
		DoubleVector	its_PL_Ratios;
		DoubleVector	its_PL_Notios;
		DoubleVector	its_PL_Spreads;

		vector<CreditPremiumLegType>		its_PL_CreditFlags;		// Guaranteed, Fix Nominal, Outstanding Nominal
		DoubleVector	its_PL_CreditSpreadCaps;
		DoubleVector	its_PL_Redemptions;
		
		CreditEventPayment		its_PL_PaymentType;			// only in case of Prorata


		// -----------------------------------------------------------------
		// PRODUCT PARAMETERS
		// -----------------------------------------------------------------

	public:
	
		void	Get_UpFront_RunningSpread(double&	data)	{data	=	its_UpFront_RunningSpread;}
		void	Get_UpFront_Premium(double&	data)	{data	=	its_UpFront_Premium;}
		
		void	Get_PricingLegsType(CreditBasketNPV&	data)	{data	=	its_PricingLegsType;}

		void	Get_DefNPVFlag(double&	data)	{data	=	its_DefNPVFlag;}
		void	Get_PremNPVFlag(double&	data)	{data	=	its_PremNPVFlag;}

		void	Get_CreditObservationType(CreditObservation&	data)	{data	=	its_CreditObservationType;}
		void	Get_CreditAccruedPayment(CreditAccruedPayment&	data)	{data	=	its_CreditPremiumLegAccrued;}
		void	Get_ATMDataFlag(CreditPremiumLegATMData&	data)	{data	=	its_ATMDataFlag;}
		void	Get_NPV_Type(CreditNPV&	data)	{data	=	its_NPV_Type;}
		void	Get_CDO_Type(CreditProductType&	data)	{data	=	its_CDO_Type;}

		void	Get_TimeStep_ProrataF(double&	data)	{data	=	its_Time_Step_Prorata;}
		void	Set_TimeStep_ProrataF(const double&	data)	{its_Time_Step_Prorata	=	data;}

		bool	IsCDOSquare() {return (its_CDO_Type == CPT_CDO_SQUARE);}
		bool	IsCDOStandard() {return (its_CDO_Type == CPT_STD_CDO);}

	private:

		double		its_UpFront_RunningSpread;
		double		its_UpFront_Premium;

		// in order to know whether it is Def-Prem vs. Prem-Def or Prem or Def alone
		CreditBasketNPV			its_PricingLegsType;
		double					its_DefNPVFlag;
		double					its_PremNPVFlag;

		CreditObservation		its_CreditObservationType;	// Do observe Losses or Nb Defaults?
		CreditAccruedPayment	its_CreditPremiumLegAccrued;
		CreditPremiumLegATMData its_ATMDataFlag;

		CreditNPV				its_NPV_Type;


		double		its_Time_Step_Prorata;

		CreditProductType		its_CDO_Type;

	// ----------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------
	// METHODS
	// ----------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------

	public:
		// Set Cash Flows Matrix
		void	SetCreditProductDefaultMatrix(ICM_Matrix<ARM_Vector>* parameters);
		void	SetCreditProductPremiumMatrix(ICM_Matrix<ARM_Vector>* parameters);
		void	SetCreditProductPricingParametersMatrix(ICM_Matrix<ARM_Vector>* parameters);

		// with its own parameters
		void	FillCreditProductDefaultMatrix()
		{
			if (its_Matrix_PL_Data)
				SetCreditProductDefaultMatrix(its_Matrix_DL_Data);
		}

		void	FillCreditProductPremiumMatrix()
		{
			if (its_Matrix_DL_Data)
				SetCreditProductPremiumMatrix(its_Matrix_PL_Data);
		}

		void	FillCreditProductPricingParametersMatrix()
		{
			if (its_Matrix_Parameters_Data)
				SetCreditProductPricingParametersMatrix(its_Matrix_Parameters_Data);
		}


};

#endif
