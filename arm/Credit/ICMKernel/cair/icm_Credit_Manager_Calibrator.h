/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CREDIT_MANAGER_CALIBRATOR.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides a basic Credit Index Calibrator


   -----------------------------------------------------------------
   
	ICM CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _ICM_CREDIT_MANAGER_CALIBRATOR_H_
# define _ICM_CREDIT_MANAGER_CALIBRATOR_H_

# include "ICMKernel\cair\types.h"
//#include "ICMKernel\cair\enums.h"
# include "ICMKernel\cair\credit_manager.h"
# include "ICMKernel\inst\icm_mez.h"

# include "ICMKernel\util\icm_brentsolver.h"

class	CreditManager;

class ICM_Credit_Manager_Calibrator
{
	// -----------------------------------------------------------
	// CONSTRUCTORS AND DESTRUCTORS
	// -----------------------------------------------------------

	public:

		// constructors
		ICM_Credit_Manager_Calibrator() {Init();}
		ICM_Credit_Manager_Calibrator(const ICM_Credit_Manager_Calibrator& Data) {Init(); Copy(Data);}
		
		// destructor
		~ICM_Credit_Manager_Calibrator()	{Reset();}

	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		virtual ICM_Credit_Manager_Calibrator*	Clone() const;
		virtual void	Destroy() {delete this;};
		virtual void	Copy(const ICM_Credit_Manager_Calibrator& Data);
		virtual void	Reset();

	private:

		void Init()
		{

			its_Attach_Points.clear();
			its_Detach_Points.clear();
			its_Market_Correlations.clear();
			its_Market_UpFrontPremiums.clear();
			its_Market_Spreads.clear();

			its_Credit_Tranches.clear();
			its_Data.clear();

			its_Correlation_Calibration_Type	=	CMCC_FROM_LHP_COMPOUND_TO_LHP_BASE;

			its_Market_Model	=	CMM_LHP_GAUSSIAN;
			its_Pricing_Model	=	CMT_ANALYTIC_RECURSIVE_1F;

			TheX	=	NULL;
			its_Steps.clear();

			its_ModelParameters.clear();
			its_CorrelationParametersToCalibrate.clear();
			its_ProductParameters.clear();

			its_Credit_Manager	=	NULL;

			its_Brent_Solver	=	NULL;

			its_Implied_Correlations.clear();
			its_Implied_Spreads.clear();
			its_Implied_UpFront_Premiums.clear();

			its_Outputs	=	NULL;
		}

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// CLASS USE
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		void	Calibrate();

		void	Compute_Market_Implied_Spreads();

		void	Compute_Implied_Correlation(int Num, double& Result);

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// INTERFACE DATA
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	public:

		void	Get_Attach_Points(DoubleVector& data)	{data = its_Attach_Points;}
		void	Set_Attach_Points(const DoubleVector& data)	{its_Attach_Points = data;}

		void	Get_Detach_Points(DoubleVector& data)	{data = its_Detach_Points;}
		void	Set_Detach_Points(const DoubleVector& data)	{its_Detach_Points = data;}

		void	Get_Market_Correlations(DoubleVector& data)	{data = its_Market_Correlations;}
		void	Set_Market_Correlations(const DoubleVector& data)	{its_Market_Correlations = data;}

		void	Get_Market_UpFrontPremiums(DoubleVector& data)	{data = its_Market_UpFrontPremiums;}
		void	Set_Market_UpFrontPremiums(const DoubleVector& data)	{its_Market_UpFrontPremiums = data;}

		void	Get_Market_Spreads(DoubleVector& data)	{data = its_Market_Spreads;}
		void	Set_Market_Spreads(const DoubleVector& data)	{its_Market_Spreads = data;}


		void	Get_Market_Model_Type(CreditManager_Model&	data) const	{data = its_Market_Model;}
		void	Set_Market_Model_Type(CreditManager_Model	data) {its_Market_Model = data;}

		void	Get_Pricing_Model_Type(CreditModelType&	data) const	{data = its_Pricing_Model;}
		void	Set_Pricing_Model_Type(CreditModelType	data) {its_Pricing_Model = data;}

		void	Get_Calibration_Type(CreditManager_Correlation_Calibration&	data) const	{data = its_Correlation_Calibration_Type;}
		void	Set_Calibration_Type(CreditManager_Correlation_Calibration	data) {its_Correlation_Calibration_Type = data;}

		void	Get_Input_Type(CreditManager_Correlation_Calibration_Input&	data) const	{data = its_Input_Type;}
		void	Set_Input_Type(CreditManager_Correlation_Calibration_Input	data) {its_Input_Type = data;}

		bool	IsBootStrapMethod() {return its_Bootstrap_Flag;}	// BRENT
		bool	IsGlobalFitMethod() {return !its_Bootstrap_Flag;}	// BFGS

		CreditManager* Get_CreditManager() {return its_Credit_Manager;}
		void	Set_CreditManager(CreditManager* data)	{its_Credit_Manager = data;}

		void	Get_Credit_Tranches(vector<ICM_Mez*>& data)	{data = its_Credit_Tranches;}
		void	Set_Credit_Tranches(const vector<ICM_Mez*>& data)	{its_Credit_Tranches = data;}


	private:

		DoubleVector	its_Attach_Points;
		DoubleVector	its_Detach_Points;
		
		DoubleVector	its_Market_Correlations;

		vector<ICM_Mez*>	its_Credit_Tranches;

		DoubleVector	its_Market_UpFrontPremiums;
		DoubleVector	its_Market_Spreads;

		// -------------------------------------------------------

		CreditManager_Model		its_Market_Model;
		
		CreditManager_Correlation_Calibration	its_Correlation_Calibration_Type;

		CreditModelType			its_Pricing_Model;
		CreditCopulaType		its_Pricing_Copula;

		CreditManager*			its_Credit_Manager;

		DoubleVector			its_Data;

		CreditManager_Correlation_Calibration_Input	its_Input_Type;
	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// INTERNAL DATA
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	private:

		void	ResizeAndAllocateData();
		void	DetermineModel();
		
		void	GetPricesFormMarketData();
		void	GetMarketParameters(DoubleVector& MarketParameters);

		// BFGS
		void	PrepareSolver();

		void	StoreBFGSResults(bool res, int nbiter, double val);

		void	RetrieveData();

		// BRENT
		void	CheckCalibrationSetForBootStrap();

		void	PrepareBrentSolver();
		void	LaunchBrentSolver(DoubleVector& Outputs);
		void	StoreBrentResults(DoubleVector& OutPuts);


	private:

		CreditManager_Model		its_Model_Choice;		// just for me to know what it is
		DoubleVector			its_Implied_UpFront_Premiums;
		DoubleVector			its_Implied_Spreads;
		DoubleVector			its_Implied_Correlations;

		DoubleVector			its_ModelParameters;
		DoubleVector			its_CorrelationParametersToCalibrate;
		DoubleVector			its_ProductParameters;

		bool	its_Bootstrap_Flag;	// (cf. Brent) otherwise BFGS

		BrentSolver*	its_Brent_Solver;
	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// SOLVER DATA
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	protected:
		
		// USED FOR BFGS algorithm
		double	ThePricingFunction(double* X);
		void	ThePricingGradientFunction(double* X, double* Gradient);

		void	ExtractRoots(double* X);

	public:

		double	GetTheImpliedSpread(double* X, int Num);
		double	GetTheNPV(double* TheX, int Num);

	protected:

		int				its_Dim_Data;
		int				its_Dim_Model;		// number of parameters to calibrate
		
		DoubleVector	its_Steps;	// to move from one point to one another
		double*			TheX;

		// OUTPUTS
		string			its_Calibration_Info;		// Least Squares final value
		int				its_Calibration_NbIter;		// NbIter
		double			its_Calibration_Val;		// Minimization Value
		bool			its_Calibration_Result;		// SUCCESS OR FAILURE?

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// MODEL DATA
	// -----------------------------------------------------------
	// -----------------------------------------------------------

		// CMM_LHP_GAUSSIAN model parameters
	
	public:

		inline void	Set_LHP_Maturity(const double& data) {its_LHP_Maturity = data;}
		inline double Get_LHP_Maturity(void) {return its_LHP_Maturity;}

		inline void	Set_LHP_Spread(const double& data) {its_LHP_Spread = data;}
		inline double Get_LHP_Spread(void) {return its_LHP_Spread;}

		inline void	Set_LHP_Recovery(const double& data) {its_LHP_Recovery = data;}
		inline double Get_LHP_Recovery(void) {return its_LHP_Recovery;}
	
	private:


		double	its_LHP_Maturity;
		double	its_LHP_Spread;
		double	its_LHP_Recovery;


		// CMM_HT_JPM_GAUSS model parameters
	
	public:

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// OUTPUTS
	// -----------------------------------------------------------
	// -----------------------------------------------------------

	private:

		ICM_QMatrix<double>*	its_Outputs;

	public:

		void	Get_Correlation_Calibrator_Outputs(ICM_QMatrix<double>*& Outputs);

};

# endif