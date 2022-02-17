#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\inst\icm_collateral.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_CREDIT_MANAGER_CALIBRATOR.CXX
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides a basic Credit Index Calibrator


   -----------------------------------------------------------------
   
	ICM CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# include "ICM_Credit_Manager_Calibrator.h"

// # include <stdlib.h>

# include "ICMKernel/glob/icm_maths.h"
# include "ICMKernel/util/icm_BFGS.h"


ICM_Credit_Manager_Calibrator* ICM_Credit_Manager_Calibrator :: Clone() const
{
	return NULL;
}


void ICM_Credit_Manager_Calibrator :: Copy(const ICM_Credit_Manager_Calibrator& Data)
{
	return;
}


void ICM_Credit_Manager_Calibrator :: Reset()
{
	its_Steps.clear();

	if (TheX)
		delete TheX;
	TheX	=	NULL;

	// is it really what I want?
//	if (its_Credit_Manager)
//		delete its_Credit_Manager;
	its_Credit_Manager	=	NULL;

	if (its_Outputs)
		delete	its_Outputs;
	its_Outputs	=	NULL;
}



//  -----------------------------------------------------------------
//  -----------------------------------------------------------------
//	CALIBRATION
//  -----------------------------------------------------------------
//  -----------------------------------------------------------------

void ICM_Credit_Manager_Calibrator :: Calibrate()
{
	// ----------------------------------------------
	// SOLVER
	// ----------------------------------------------
	bool	res;
	double	val;
	int		nbiter;

	// ----------------------------------------------
	
	// -------------------------
	// STEP 0: Set the models
	// -------------------------

	// -------------------------
	// STEP 1: Set the dimensionalities
	// -------------------------
	its_Dim_Data	=	its_Attach_Points.size();

	if (its_Dim_Data <= 0)
		ICMTHROW(ERR_INVALID_ARGUMENT,"DATA dimension is Null. No points to calibrate!");

	ResizeAndAllocateData();

	// -----------------------
	// STEP 2: GET MARKET PRICES
	// -----------------------
	GetPricesFormMarketData();

	if (IsGlobalFitMethod())
	{
		// -----------------------
		// STEP 3: PREPARE SOLVER - Model dependant
		// -----------------------
		PrepareSolver();

		// -----------------------
		// STEP 4: LAUNCH THE ROUTINE
		// -----------------------
		res = BFGS(	ff1::mem_call(&ICM_Credit_Manager_Calibrator::ThePricingFunction, (*this)),
							ff1::mem_call(&ICM_Credit_Manager_Calibrator::ThePricingGradientFunction, (*this))
						).Minimize(TheX, its_Dim_Model, 0.000001, &nbiter, &val );	

		// -----------------------
		// STEP 5: GET RESULTS: and DISPLAY
		// -----------------------
		StoreBFGSResults(res, nbiter, val);

		// -----------------------
		// STEP 6: RETRIEVE DATA
		// -----------------------
		RetrieveData();
	}
	else
	{
		DoubleVector	OutPuts;

		// BRENT SOLVER		
		PrepareBrentSolver();

		LaunchBrentSolver(OutPuts);

		StoreBrentResults(OutPuts);
	}
}



void ICM_Credit_Manager_Calibrator :: ResizeAndAllocateData()
{
	// its_dim has been validated
	its_Implied_Spreads.resize(its_Dim_Data);
	its_Implied_Correlations.resize(its_Dim_Data);
	its_Implied_UpFront_Premiums.resize(its_Dim_Data);

	// Determine Model
	// According to the Pricing Model
	// computes its_Dim_Model
	DetermineModel();

	if (its_Dim_Model)
	{
		its_Steps.resize(its_Dim_Model);

		if (TheX)
			delete	TheX;
		TheX	=	new double [its_Dim_Model];

		// the real Data they may be initialized from interface
		its_Data.resize(its_Dim_Model);
	}

	// vector of Mezz tranches
	its_Credit_Tranches.resize(its_Dim_Data);

}


void ICM_Credit_Manager_Calibrator :: 	DetermineModel()
{
	switch (its_Correlation_Calibration_Type)
	{
	case CMCC_FROM_LHP_COMPOUND_TO_LHP_BASE:

		its_Market_Model	=	CMM_LHP_GAUSSIAN;
		its_Pricing_Model	=	CMT_ANALYTIC_LHP;
		its_Model_Choice	=	CMM_LHP_GAUSSIAN;

		its_ModelParameters.clear();
		its_ModelParameters.push_back(its_LHP_Maturity);		// or Mty?
		its_ModelParameters.push_back(its_LHP_Spread);	
		its_ModelParameters.push_back(its_LHP_Recovery);

		its_CorrelationParametersToCalibrate.clear();
		its_CorrelationParametersToCalibrate.resize(2);

		its_Bootstrap_Flag	=	true;
		
		its_Dim_Model	=	1;

		break;

	case CMCC_FROM_LHP_JPM_COMPOUND_TO_LHP_JPM_BASE:

		its_Market_Model	=	CMM_LHP_JPM_GAUSSIAN;
		its_Pricing_Model	=	CMT_ANALYTIC_LHP_JPM;
		its_Model_Choice	=	CMM_LHP_JPM_GAUSSIAN;

		its_ModelParameters.clear();
		its_ModelParameters.push_back(its_LHP_Maturity);		// or Mty?
		its_ModelParameters.push_back(its_LHP_Spread);	
		its_ModelParameters.push_back(its_LHP_Recovery);

		its_CorrelationParametersToCalibrate.clear();
		its_CorrelationParametersToCalibrate.resize(2);

		its_Bootstrap_Flag	=	true;
		
		its_Dim_Model	=	1;

		break;

	case CMCC_FROM_HT_JPM_COMPOUND_TO_HT_JPM_BASE:

		its_Market_Model	=	CMM_ANDERSEN_GAUSSIAN;
		its_Pricing_Model	=	CMT_ANALYTIC_RECURSIVE_1F;
		its_Model_Choice	=	CMM_ANDERSEN_GAUSSIAN;

		its_ModelParameters.clear();

		its_CorrelationParametersToCalibrate.clear();
		its_CorrelationParametersToCalibrate.resize(2);

		its_Bootstrap_Flag	=	true;
		
		its_Dim_Model	=	1;

		break;

	case CMCC_FROM_LHP_JPM_COMPOUND_TO_LHP_JPM_NIG:

		its_Market_Model	=	CMM_LHP_JPM_GAUSSIAN;
		its_Pricing_Model	=	CMT_ANALYTIC_LHP_JPM;
		its_Model_Choice	=	CMM_LHP_JPM_NIG;

		its_ModelParameters.clear();
		its_ModelParameters.push_back(its_LHP_Maturity);		// or Mty?
		its_ModelParameters.push_back(its_LHP_Spread);	
		its_ModelParameters.push_back(its_LHP_Recovery);

		its_CorrelationParametersToCalibrate.clear();
		its_CorrelationParametersToCalibrate.resize(4);

		its_Bootstrap_Flag	=	false;

		its_Dim_Model	=	3;

		break;

	}
/*
	switch (its_Pricing_Model)
	{
	case CMT_ANALYTIC_RECURSIVE_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"CMT_ANALYTIC_RECURSIVE_1F not done!");
		
		switch (its_Pricing_Copula)
		{
		case CCT_GAUSSIAN:
			break;

		case CCT_STUDENT:
			break;

		case CCT_NIG:
			break;
		}

		break;

	case CMT_ANALYTIC_LHP:

		ICMTHROW(ERR_INVALID_ARGUMENT,"CMT_ANALYTIC_LHP not done!");
			
		switch (its_Pricing_Copula)
		{
		case CCT_GAUSSIAN:
			its_Model_Choice	=	CMM_LHP_GAUSSIAN;
			break;

		case CCT_STUDENT:
			break;

		case CCT_NIG:
			break;
		}
			
		break;

	case CMT_ANALYTIC_LHP_JPM:

		its_Model_Choice	=	CMM_HT_JPM_GAUSS;
		ICMTHROW(ERR_INVALID_ARGUMENT,"CMT_ANALYTIC_LHP not done!");
			
		break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:

		// Copula is Gaussian
		if (its_Pricing_Copula != CCT_GAUSSIAN)
			ICMTHROW(ERR_INVALID_ARGUMENT,"GAUSSIAN Copula expected for CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F!");
		// two correlation levels + One Bernoulli proba
		
		its_Model_Choice	=	CMM_RSCB;
		its_Dim_Model	=	3;
		
		break;
	
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"Unknown Model!");
	}
*/
	its_ProductParameters.clear();
	its_ProductParameters.resize(2);

}
	

void ICM_Credit_Manager_Calibrator :: GetPricesFormMarketData()
{
	int	i;
	double	Current_Implied_Spread;
	double	Current_Implied_Correlation;
	double	UpFrontPremium;
	double	RPV01;

	ICM_Mez			*Current_Product;
	DoubleVector	MarketParameters;

	if (its_Credit_Manager == NULL)
		ICMTHROW(ERR_INVALID_ARGUMENT,"NULL Credit Manager!");

	GetMarketParameters(MarketParameters);

	switch (its_Input_Type)
	{
	case CMCCI_COMPOUND_CORRELATION:

		its_CorrelationParametersToCalibrate.resize(2);

		// the first parameter is the type of Correlation
		its_CorrelationParametersToCalibrate[0]	=	(double) CMCT_COMPOUND_CORRELATION;

		// loop all Attachment and Detachment points
		for (i=0; i<its_Dim_Data; i++)
		{
			Current_Product	=	its_Credit_Tranches[i];
		
			UpFrontPremium	=	its_Market_UpFrontPremiums[i];

			if (UpFrontPremium > 0.0)
			{
				its_ProductParameters[0]	=	1.0;			// up-front premium FLAG
				its_ProductParameters[1]	=	UpFrontPremium;	// up-front premium VALUE
			}
			else
			{
				its_ProductParameters[0]	=	0.0;
				its_ProductParameters[1]	=	0.0;
			}

			// the implied (compound) correlation
			its_CorrelationParametersToCalibrate[1]	=	its_Market_Correlations[i];

			if (Current_Product == NULL)
				ICMTHROW(ERR_INVALID_ARGUMENT,"NULL Product!");

			its_Credit_Manager->Price_CreditProduct(Current_Product, its_Market_Model, MarketParameters, its_CorrelationParametersToCalibrate, its_ProductParameters);

			its_Credit_Manager->Get_ATM_Spread(Current_Implied_Spread);
			its_Implied_Spreads[i]	=	Current_Implied_Spread;
			its_Implied_Correlations[i]	=	its_Market_Correlations[i];
		}

		break;

	case CMCCI_SPREADS:

		// ---------------------------------------------
		// Market Implied Correlations are not relevant
		// ---------------------------------------------
		
		its_CorrelationParametersToCalibrate.resize(2);

		// the first parameter is the type of Correlation
		its_CorrelationParametersToCalibrate[0]	=	(double) CMCT_COMPOUND_CORRELATION;

		// loop all Attachment and Detachment points
		for (i=0; i<its_Dim_Data; i++)
		{
			Current_Product	=	its_Credit_Tranches[i];

			Current_Implied_Spread	=	its_Market_Spreads[i];
			its_Implied_Spreads[i]	=	Current_Implied_Spread;
			Current_Product->SetCreditSpread(Current_Implied_Spread / 10000.0);	// take care of the form : comes in bps

			UpFrontPremium			=	its_Market_UpFrontPremiums[i];

			if (UpFrontPremium > 0.0)
			{
				its_ProductParameters[0]	=	1.0;			// up-front premium FLAG
				its_ProductParameters[1]	=	UpFrontPremium;	// up-front premium VALUE
			}
			else
			{
				its_ProductParameters[0]	=	0.0;
				its_ProductParameters[1]	=	0.0;
			}

			Compute_Implied_Correlation(i, Current_Implied_Correlation);

			its_Implied_Correlations[i]	=	Current_Implied_Correlation;
			
			// the implied (compound) correlation (what we are looking for)
			its_CorrelationParametersToCalibrate[1]	=	0.0;

			// if I want to have the fair spread with no up-Front Premium
			if (UpFrontPremium > 0.0)
			{
				its_Credit_Manager->Get_ATM_PremiumLegWithoutNotio(RPV01);
				its_Implied_Spreads[i]	=	its_Market_Spreads[i] + UpFrontPremium / RPV01;
			}
			
		}

		break;

	case CMCCI_BASE_CORRELATION:
		break;
	}
}


void ICM_Credit_Manager_Calibrator :: GetMarketParameters(DoubleVector& MarketParameters)
{
	ARM_Date maturity;
//	double Mty;

//	ARM_Date Asof = GetAsOfDate();

	switch (its_Market_Model)
	{
	case CMM_LHP_GAUSSIAN:
	case CMM_LHP_JPM_GAUSSIAN:
	case CMM_LHP_JPM_NIG:
		MarketParameters.clear();
		
		// Maturity
//		maturity	=	its_Credit_Tranches[0]->GetMaturity();
//		Mty			=	CountYears(KACTUAL_360, Asof, maturity);

		MarketParameters.push_back(its_LHP_Maturity);		// or Mty?
		MarketParameters.push_back(its_LHP_Spread);	
		MarketParameters.push_back(its_LHP_Recovery);

		break;

	case CMM_HT_JPM_GAUSSIAN:
		MarketParameters.clear();
		break;

	case CMM_RSCB:
		ICMTHROW(ERR_INVALID_ARGUMENT,"Recursive Stochastic Correlation Bernoulli is not a Market Model yet!");
		break;

	case CMM_ANDERSEN_GAUSSIAN:
		MarketParameters.clear();
		break;
	}
}


void ICM_Credit_Manager_Calibrator :: PrepareSolver()
{
	// minimization Step
	double	TheStep	=	1e-6;

	switch (its_Model_Choice)
	{
		case CMM_LHP_GAUSSIAN:
		case CMM_LHP_JPM_GAUSSIAN:
		case CMM_HT_JPM_GAUSSIAN:
		case CMM_ANDERSEN_GAUSSIAN:
			
			its_Steps[0]	=	TheStep;
			TheX[0]	=	0.5;
			break;

		case CMM_LHP_JPM_NIG:

			// first is ALPHA 0 <= |BETA| < ALPHA
			// second is BETA in ]-inf ; +inf[
			// third is RHO in [0;1]

			its_CorrelationParametersToCalibrate[0]	=	(double) CMCT_NIG_CORRELATION;

			// ALPHA
			its_Steps[0]	=	TheStep;
			TheX[0]			=	0.5;		//	alpha(0) =	0.0

			// BETA
			its_Steps[1]	=	TheStep;
			TheX[1]			=	0.0;		//	beta(0)	=	0.0

			// RHO
			its_Steps[2]	=	TheStep;
			TheX[2]			=	(acos(0.5) + 1.0) / 2.0;	// start with 50% correlation 1.04719755119660.5;		//	(acos(0.5) + 1.0) / 2.0;

			break;
			
		case CMM_RSCB:

			// first is RHO	in [0;1]
			// second is BETA in [0;1]
			// third is q = Prob(Bernoulli = 1) in [0;1]

			its_Steps[0]	=	TheStep;
			TheX[0]			=	0.5;		//	(acos(0.5) + 1.0) / 2.0;

			its_Steps[1]	=	TheStep;
			TheX[1]			=	0.5;		//	(acos(0.5) + 1.0) / 2.0;

			its_Steps[2]	=	TheStep;
			TheX[2]			=	0.5;		//	(acos(0.5) + 1.0) / 2.0;
		break;

	}
	
}


void ICM_Credit_Manager_Calibrator :: StoreBFGSResults(bool res, int nbiter, double val)
{
	its_Calibration_Info.empty();

	if (res)
		ICMLOG("BFGS:: SUCCEEDED")
	else
		ICMLOG("BFGS:: FAILED")

	ICMLOG("BFGS::Results obtained after " << nbiter << " iterations and Val " << val)

	// store OUTPUTS

	its_Calibration_NbIter	=	nbiter;
	its_Calibration_Val		=	val;
	its_Calibration_Result	=	res;

	ostringstream	temp;
	if (res)
		temp	<< "Calibration SUCCEEDED!";
	else
		temp	<< "Calibration FAILED!";
	temp	<< " - Results obtained after " << nbiter << " iterations with a Minimization value of " << val << ".";
	
	its_Calibration_Info	=	temp.str();

}


void ICM_Credit_Manager_Calibrator :: RetrieveData()
{
	switch (its_Model_Choice)
	{
		case CMM_LHP_GAUSSIAN:
		case CMM_LHP_JPM_GAUSSIAN:
		case CMM_HT_JPM_GAUSSIAN:
			its_Data[0]	=	TheX[0];
			break;

			ICMTHROW(ERR_INVALID_ARGUMENT,"CMT_ANALYTIC_RECURSIVE_1F not done!");
			break;

		case CMM_LHP_JPM_NIG:
			// Beta
			its_Data[1]	=	TheX[1];

			// Alpha
			its_Data[0]	=	sqrt(exp(TheX[0]) + TheX[1] * TheX[1]);

			// Rho
			its_Data[2]	=	(acos(TheX[2]) + 1.0) / 2.0;

			break;
		case CMM_RSCB:

			// first is RHO	in [0;1]
			// second is BETA in [0;1]
			// third is q = Prob(Bernoulli = 1) in [0;1]

			its_Data[0]	=	cos(2*TheX[0] - 1.0);

			its_Data[1]	=	cos(2*TheX[1] - 1.0);

			its_Data[2]	=	cos(2*TheX[2] - 1.0);
		
		break;

	}
}



double	ICM_Credit_Manager_Calibrator::ThePricingFunction(double* X)
{
	int	i;
	double	TheImpliedSpread, TheTargetSpread;
	double	Sum	=	0.0;

	for (i=0; i<its_Dim_Data; i++)
	{
		TheImpliedSpread	=	GetTheImpliedSpread(X, i);
		TheTargetSpread		=	its_Implied_Spreads[i];

		Sum	+=	(TheImpliedSpread - TheTargetSpread) * (TheImpliedSpread - TheTargetSpread);
	}

	return	Sum;
}


void	ICM_Credit_Manager_Calibrator::ThePricingGradientFunction(double* X, double* Gradient)
{
	int	i;
	double	ResultUp, ResultDown, TheStep;
	
	for (i=0; i<its_Dim_Model; i++)
	{
		TheStep	=	its_Steps[i];
		// function evaluation
		X[i] += TheStep;		
		ResultUp	=	ThePricingFunction(X);

		// second function evluation
		X[i] -= 2.0 * TheStep;
		ResultDown	=	ThePricingFunction(X);

		// restore current point
		X[i] += TheStep;

		// set the gradient
		Gradient[i] = (ResultUp - ResultDown) / (2.0 * TheStep);
	}
}


void	ICM_Credit_Manager_Calibrator::ExtractRoots(double* X)
{
	double	tmpValue;

	switch (its_Model_Choice)
	{
		case CMM_LHP_GAUSSIAN:
		case CMM_LHP_JPM_GAUSSIAN:
		case CMM_HT_JPM_GAUSSIAN:
		case CMM_ANDERSEN_GAUSSIAN:
			its_CorrelationParametersToCalibrate[1]	=	X[0];	// nothing to do
			break;

		case CMM_LHP_JPM_NIG:

			// ALPHA
			tmpValue	=	X[0]*X[0] - X[1]*X[1];
			if (X[0]*X[0] - X[1]*X[1] > 0.0)
				its_CorrelationParametersToCalibrate[1]	=	log(tmpValue);
			else
				its_CorrelationParametersToCalibrate[1]	=	-10.0;

			// BETA
			its_CorrelationParametersToCalibrate[2]	=	X[1];
			// RHO
			its_CorrelationParametersToCalibrate[3]	=	X[2];

			break;

		case CMM_RSCB:

			// first is RHO	in [0;1]
			// second is BETA in [0;1]
			// third is q = Prob(Bernoulli = 1) in [0;1]

			its_CorrelationParametersToCalibrate[1]	=	cos(2*X[0] - 1.0);

			its_CorrelationParametersToCalibrate[2]	=	cos(2*X[1] - 1.0);

			its_CorrelationParametersToCalibrate[3]	=	cos(2*X[2] - 1.0);
		
		break;

	}
}



double	ICM_Credit_Manager_Calibrator::GetTheImpliedSpread(double* TheX, int Num)
{
	double	Current_Implied_Spread;
	ICM_Mez*	Current_Product;

	// the tests have already been carried out
	Current_Product	=	its_Credit_Tranches[Num];

//	ModelParameters[1]	=	TheX[0];
	ExtractRoots(TheX);	// from real values (cf. BFGS algorithm) to model values

	its_Credit_Manager->Price_CreditProduct(Current_Product, its_Model_Choice, its_ModelParameters, its_CorrelationParametersToCalibrate, its_ProductParameters);

	its_Credit_Manager->Get_ATM_Spread(Current_Implied_Spread);

	return	Current_Implied_Spread;
}


double	ICM_Credit_Manager_Calibrator::GetTheNPV(double* TheX, int Num)
{
	double		Current_NPV;
	ICM_Mez*	Current_Product;

	// the tests have already been carried out
	Current_Product	=	its_Credit_Tranches[Num];

//	ModelParameters[1]	=	TheX[0];
	ExtractRoots(TheX);	// from real values (cf. BFGS algorithm) to model values

	its_Credit_Manager->Price_CreditProduct(Current_Product, its_Model_Choice, its_ModelParameters, its_CorrelationParametersToCalibrate, its_ProductParameters);

	its_Credit_Manager->Get_NPV(Current_NPV);

	return	Current_NPV;
}



void	ICM_Credit_Manager_Calibrator::PrepareBrentSolver()
{
	if (its_Brent_Solver)
		delete	its_Brent_Solver;

	its_Brent_Solver	=	new BrentSolver();

//	int	NbMaxIter;
	// SET NB MAX ITER
//	its_Brent_Solver->GetNbMaxIter(NbMaxIter);


}

void	ICM_Credit_Manager_Calibrator::CheckCalibrationSetForBootStrap()
{
	int	i;
	int	NbTranches;

	double	Last_Attach;
	double	Last_Detach;
	double	Current_Attach;
	double	Current_Detach;

	Last_Attach	=	its_Attach_Points[0];
	Last_Detach	=	its_Detach_Points[0];

	if (Last_Attach != 0.0)
		ICMTHROW(ERR_INVALID_ARGUMENT,"FIRST ATTACHMENT point must be 0.0!");

	if ((Last_Detach >= 1.0) || (Last_Detach <= 0.0))
		ICMTHROW(ERR_INVALID_ARGUMENT,"FIRST DETACHMENT point must lie in ]0;1[!");

	NbTranches	=	its_Market_Correlations.size();

	for (i=1; i<NbTranches; i++)
	{
		Current_Attach	=	its_Attach_Points[i];
		Current_Detach	=	its_Detach_Points[i];

		if (Current_Attach != Last_Detach)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ATTACHMENT POINT must be recovered with last DETACHMENT POINT!");

		if (Current_Attach >= Current_Detach)
			ICMTHROW(ERR_INVALID_ARGUMENT,"CHECK that ATTACHMENT POINT < DETACHMENT POINT! for Attachment input " << Current_Attach);

		if (Current_Detach > 1.0)
			ICMTHROW(ERR_INVALID_ARGUMENT,"DETACHMENT POINT must be <= 1.0! Problem with input row number " << i+1);

	}

}


// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

extern "C" unsigned short NPVFunction(void* Param, double X, double& Result)
{
	int	*Current_Num;
	ICM_Credit_Manager_Calibrator	*TheCorrelationCalibrator;

	TheCorrelationCalibrator	=	(ICM_Credit_Manager_Calibrator*)((*((vector<void*>*)Param))[0]);
	Current_Num					=	(int*) ((*((vector<void*>*)Param))[1]);

	double*	tmpX;
	tmpX	=	new double[1];
	
	tmpX[0]	=	X;

	// BASE CORRELATION FITTING, match NPV
	Result = TheCorrelationCalibrator->GetTheNPV(tmpX, *Current_Num);

	delete [] tmpX;

	return	RetOk;	// 1
}


void	ICM_Credit_Manager_Calibrator::LaunchBrentSolver(DoubleVector& Outputs)
{
	int		i, j, NbTranches;
	int		nb_iter;

	bool	CalibSuccessFlag;

	double	Current_Correl;
	double	Last_Attach;
	double	Last_Detach;
	double	Current_Attach;
	double	Current_Detach;

	NbTranches	=	its_Market_Correlations.size();
	
	Outputs.clear();
	Outputs.resize(NbTranches);

	// The first Equity is the same...
	Current_Correl	=	its_Implied_Correlations[0];	// must have been computed (if Spreads are MARKET inputs)!
	Outputs[0]		=	Current_Correl;

 	double	SolverBrentMin;
	double	SolverBrentMax;
	double	SolverBrentTol	=	1e-4;

	int	NbMaxIter;
	its_Brent_Solver->GetNbMaxIter(NbMaxIter);

	// MUST BE SORTED!
	Last_Attach	=	its_Attach_Points[0];
	Last_Detach	=	its_Detach_Points[0];

	// for parameters
	vector<void*>	TheVector;
	
	TheVector.resize(2);

	i	=	0;

	TheVector[0]	=	this;
	TheVector[1]	=	&i;

	double	KeepSubAmount;
	double	KeepMezzAmount;
	double	Notional;
	double	Current_NPV;
	double	Current_ImpliedSpread;

	ICM_Mez*	Current_Product;

	its_Brent_Solver->SetNbMaxIter(10);

	for (i=1; i<NbTranches; i++)
	{
//		TheVector[1]	=	&i;

		CalibSuccessFlag	=	true;

		Current_Attach	=	its_Attach_Points[i];
		Current_Detach	=	its_Detach_Points[i];

 		SolverBrentMin	=	Current_Correl;				// can not decrease
		SolverBrentMax	=	FMIN(2 * Current_Correl * Current_Detach / Last_Detach, 0.999);	// may be pushed higher for steep curves, otherwise 

		its_Brent_Solver->SetTarget(its_Implied_Spreads[i]);

		// ------------
		// FIRST PART, evaluates [0%; Current_Attach] at [Current_Attach; Current_Detach] fair spread
		// ------------
		Current_Product	=	its_Credit_Tranches[i];
		
		KeepSubAmount	=	Current_Product->GetSubAmount(Current_Product->GetFeeLeg()->GetStartDate());
		KeepMezzAmount	=	Current_Product->GetMezzAmount(Current_Product->GetFeeLeg()->GetStartDate());
		Notional		=	Current_Product->GetCollateral()->SumNotionals(Current_Product->GetStartDateNA());

		Current_ImpliedSpread	=	its_Implied_Spreads[i];
		Current_Product->SetCreditSpread(Current_ImpliedSpread / 10000.0);	// take care of the form : comes in bps

		Current_Product->SetSubAmount(0.0);
		Current_Product->SetMezzAmount(Current_Attach * Notional);
		
		
		// the implied (compound) correlation
		its_CorrelationParametersToCalibrate[1]	=	Outputs[i-1];

		its_Credit_Manager->Price_CreditProduct(Current_Product, its_Market_Model, its_ModelParameters, its_CorrelationParametersToCalibrate, its_ProductParameters);

		its_Credit_Manager->Get_NPV(Current_NPV);


		// create product [0% ; Current_Detach]
		Current_Product->SetMezzAmount(Current_Detach * Notional);

		its_Brent_Solver->SetTarget(Current_NPV);

		try
		{
			// SOLVER
			its_Brent_Solver->Solve(SolverBrentMin, SolverBrentMax, SolverBrentTol, nb_iter, NPVFunction, &TheVector, Current_Correl);
			
			if ((nb_iter > NbMaxIter) || (Current_Correl <= 0.0))
			{
				// maybe one more try with a bigger range!

				SolverBrentMax	=	FMIN(SolverBrentMax * 2.0, 0.999);

				its_Brent_Solver->Solve(SolverBrentMin, SolverBrentMax, SolverBrentTol, nb_iter, NPVFunction, &TheVector, Current_Correl);

				if ((nb_iter > NbMaxIter) || (Current_Correl <= 0.0))
					CalibSuccessFlag	=	false;

				// do not forget!
				nb_iter += NbMaxIter;
			}
		}

		catch (...)
		{
			CalibSuccessFlag	=	false;
		}


		// FINALLY
		if (!CalibSuccessFlag)
		{
			for (j=i; j<NbTranches; j++)
				Outputs[j]	=	-1.0;

			ICMLOG("CORRELATION Calibrate failed for row input number " << i+1);

			return;
		}

		Current_Product->SetSubAmount(KeepSubAmount);

		Outputs[i]	=	Current_Correl;

		// NEXT POINT
		Last_Attach	=	Current_Attach;
		Last_Detach	=	Current_Detach;

	}
}


void	ICM_Credit_Manager_Calibrator::StoreBrentResults(DoubleVector& OutPuts)
{
	its_Outputs	=	new ICM_QMatrix<double>(its_Dim_Data, 1);

	for (int i=0; i<its_Dim_Data; i++)
	{
		its_Outputs->SetValue(i, 0, OutPuts[i]);
	}

}



// -----------------------------------------------------------
// -----------------------------------------------------------
// OUTPUTS
// -----------------------------------------------------------
// -----------------------------------------------------------

void	ICM_Credit_Manager_Calibrator::Get_Correlation_Calibrator_Outputs(ICM_QMatrix<double>*& Outputs)
{
	int	i;

	// allocations
	// number of rows: Nb Tranches
	// number of columns: 3 (correl, up-front, spread // model dependance)
	int	NbRows;
	int	NbCols;

	NbRows	=	its_Dim_Data;
	NbCols	=	its_Dim_Model;

	if (Outputs)
		delete	Outputs;
	Outputs	=	new ICM_QMatrix<double>(NbRows, 3 + NbCols);

	for (i=0; i<NbRows; i++)
	{
		// IMPLIED CORREL
		Outputs->SetValue(i, 0, its_Implied_Correlations[i]);

		// IMPLIED UP-FRONT PREMIUM
		Outputs->SetValue(i, 1, its_Implied_UpFront_Premiums[i]);

		// IMPLIED SPREAD
		Outputs->SetValue(i, 2, its_Implied_Spreads[i]);
	}

	// MODEL DEPENDANT
	switch (its_Correlation_Calibration_Type)
	{
	case CMCC_FROM_LHP_COMPOUND_TO_LHP_BASE:
	case CMCC_FROM_LHP_JPM_COMPOUND_TO_LHP_JPM_BASE:
	case CMCC_FROM_HT_JPM_COMPOUND_TO_HT_JPM_BASE:

		for (i=0; i<NbRows; i++)
		{
			// BASE CORREL
			Outputs->SetValue(i, 3, (*its_Outputs)(i, 0));
		}

		break;
	}
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// COMPUTE IMPLIED CORRELATION
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void	ICM_Credit_Manager_Calibrator::Compute_Implied_Correlation(int Num, double& Result)
{
	Result	=	1.0;

	// --------------------------------------------------------
	// STORE
	DoubleVector	Kept_CorrelationParametersToCalibrate;

	Kept_CorrelationParametersToCalibrate	=	its_CorrelationParametersToCalibrate;
	// --------------------------------------------------------

	if (its_Brent_Solver)
		delete	its_Brent_Solver;

	its_Brent_Solver	=	new BrentSolver();

	int	NbMaxIter;
	its_Brent_Solver->GetNbMaxIter(NbMaxIter);

	// Get 0.0 NPV
	its_Brent_Solver->SetTarget(0.0);

	int	nb_iter;
	double	SolverBrentMin;
	double	SolverBrentMax;
	double	SolverBrentTol;

	SolverBrentMin	=	0.0;
	SolverBrentMax	=	0.99999;
	SolverBrentTol	=	1e-6;

	// for parameters
	vector<void*>	TheVector;
	
	TheVector.resize(2);

	TheVector[0]	=	this;
	TheVector[1]	=	&Num;

	its_CorrelationParametersToCalibrate.clear();
	its_CorrelationParametersToCalibrate.resize(2);
	
	its_CorrelationParametersToCalibrate[0]	=	(double) CMCT_COMPOUND_CORRELATION;
	its_CorrelationParametersToCalibrate[1]	=	0.25;

	// SOLVER
	// with Zero Bracket?
	its_Brent_Solver->Solve(SolverBrentMin, SolverBrentMax, SolverBrentTol, nb_iter, NPVFunction, &TheVector, Result, true, 2, 10);
	
	if ((nb_iter > NbMaxIter) || (Result < 0.0) || (Result > 1.0))
		ICMLOG("Compute Implied Correlation failed!");

	// --------------------------------------------------------
	// RESTORE

	its_CorrelationParametersToCalibrate	=	Kept_CorrelationParametersToCalibrate;
	// --------------------------------------------------------
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------