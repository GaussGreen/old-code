/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file numerical.cpp
 *
 *  \brief model fitter that makes the model fit a whole smile
 *	\author  A Schauly
 *	\version 1.0
 *	\date August 2005
 */


/// gpbase
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/datestrip.h"
#include "gpbase/curve.h"
#include "gpbase/surface.h"
#include "gpbase/warningkeeper.h"

/// gpcalib
#include "gpcalib/numerical.h"
#include "gpcalib/vanillasecuritydensity.h"
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"


/// Will be deleted on the long run
#include "gpcalib/densityfunctors.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"

/// kernel
#include <inst/portfolio.h>

/// gpmodels
#include "gpmodels/MarkovFunctional.h"
#include "gpmodels/ModelParamsMF.h"

/// gpnummethods
#include "gpnummethods/pdemethod.h"

///gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/nagsolver.h"


/// code for variance squeeze exception ...
#define VAR_SQUEEZE_EXCEPTION 12


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Constructor (default)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NumericalModelFitter::ARM_NumericalModelFitter()
:	ARM_ModelFitter (),  
	itsCalibSecDensities( 0 ),
	itsCalibSchedule	( NULL ),
	itsCalibStates		( NULL ),
	itsFunction			( ),
	itsFunctionWithNumDerivative( NULL ),
	itsSolver			( NULL )
{
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Constructor (contextual)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NumericalModelFitter::ARM_NumericalModelFitter(	ARM_PricingModel* model,
													ARM_DateStrip* calibSchedule, 
													const ARM_VanillaSecDensityPtrVector& calibSecDensities,
													const ARM_StdPortfolioPtr portfolio,

													/// -> Bootstrap specific
													const size_t	max_iter,
													bool			printLevel,
													double			fxTolerance,
													double			xTolerance,
													const size_t	dicho_max_iter,
													double			dichoXTolerance,
													ARM_SolverType	type,
													size_t			factorNb,
													size_t			NbCalibIter,
													ARM_MktTargetType  targetType)
													
:	ARM_ModelFitter (	model, 
						portfolio, 
						ARM_ModelParamVector(0), 
						ARM_ModelFitterPtr(NULL), 
						ARM_ModelFitterPtr(NULL),
						max_iter,
						printLevel, 
						factorNb, 
						NbCalibIter,
						targetType   ), 
	
	itsCalibSchedule		( CreateClonedPtr(calibSchedule) ),
	itsCalibSecDensities	( calibSecDensities ),
	itsCalibStates			( NULL ),
	itsFunction				(  ),
	itsFunctionWithNumDerivative( NULL ),
	itsSolver				( NULL )
{
	/// a little check...
	if(itsCalibSecDensities.size() == 0 && itsCalibSchedule->GetResetDates()->size() > 0)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : the nb of securities is required to match calibration schedule size..." );
	}

	if ( itsCalibSchedule->GetResetDates()->size() != itsCalibSecDensities.size() )
	{
		if(dynamic_cast<ARM_VanillaSecurityDensitySpread*>(&*itsCalibSecDensities[0]) != NULL && itsCalibSchedule->GetResetDates()->size() > itsCalibSecDensities.size())
		{
			// pas grave pour les spread
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : the nb of securities is required to match calibration schedule size..." );
	}

	/// Associate calibration schedule
	for (size_t i(0); i<itsCalibSecDensities.size(); i++)
	{
		itsCalibSecDensities[i]->setZeroCurve( GetPricingModel()->GetZeroCurve() ); /// Ok, c'est un peu du bricolage
		itsCalibSecDensities[i]->AssociateCalibSchedule(itsCalibSchedule);
	}

	/// Bootstrap specific
	if (!GetPortfolio().IsNull())
	{		
		StoreDetailsInFileIfRequired();
		Validate();

		/// eps for numerical derivative
		double eps = 1.0e-8;
		itsFunctionWithNumDerivative = new UnaryFuncWithNumDerivative<double>( itsFunction, eps);
			
		const double target   = 0.0;	

		switch( type )
		{
		case ARM_ModelFitterSolverType::NewtonRaphson:
			itsSolver   =  new T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> >(
				*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		case  ARM_ModelFitterSolverType::NewtonRaphsonWithRetrial:
			itsSolver   = new T_NewtonRaphsonSolverRetrial< UnaryFuncWithNumDerivative<double> >(
				*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		case ARM_ModelFitterSolverType::SmoothNewthonRhaphson:
			itsSolver   = new T_SmoothNewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> >(
				*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		case ARM_ModelFitterSolverType::NewthonRhaphsonNoThrow:
			itsSolver   = new T_NewtonRaphsonSolverNoThrow< UnaryFuncWithNumDerivative<double> >(
				*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		case ARM_ModelFitterSolverType::NewtonRaphsonWithDichotomy:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : NewtonRaphsonWithDichotomy not implemented..." );
			break;
		case ARM_ModelFitterSolverType::Dichotomy:
			itsSolver   = new T_DichotomySolver< UnaryFuncWithNumDerivative<double> >(
					*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		case ARM_ModelFitterSolverType::NagSolver:
			itsSolver   = new T_NagSolver< UnaryFuncWithNumDerivative<double> >(
					*itsFunctionWithNumDerivative,target, fxTolerance,xTolerance, max_iter, printLevel );
			break;
		default:
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,	"ARM_SolverType::unknown type!" );
		}

		Initialise();
	}

}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Constructor
///	Returns: 
///	Action : copy-const
////////////////////////////////////////////////////
ARM_NumericalModelFitter::ARM_NumericalModelFitter( const ARM_NumericalModelFitter& rhs )
:	ARM_ModelFitter(rhs),
	itsCalibSchedule (rhs.itsCalibSchedule), /// ASSOC
	itsCalibSecDensities(rhs.itsCalibSecDensities), /// ASSOC
	itsCalibStates(CreateClonedPtr(&*rhs.itsCalibStates)),
	itsFunction(rhs.itsFunction),
	/// itsFunctionWithNumDerivative(new UnaryFuncWithNumDerivative<double> (*rhs.itsFunctionWithNumDerivative)),
	itsSolver (rhs.itsSolver ? rhs.itsSolver->Clone() : NULL)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NumericalModelFitter::~ARM_NumericalModelFitter()
{
	if (itsSolver)
		delete itsSolver;

	if (itsFunctionWithNumDerivative)
		delete itsFunctionWithNumDerivative;
}



////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_NumericalModelFitter& ARM_NumericalModelFitter::operator=(const ARM_NumericalModelFitter& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_NumericalModelFitter();
		new (this) ARM_NumericalModelFitter (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: toString
///	Returns: string
///	Action : 
////////////////////////////////////////////////////
string ARM_NumericalModelFitter::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << indent <<  "\n";
    os << indent << "ARM_NumericalModelFitter\n";
    os << indent << "---------------------------\n";
	os << indent << CC_NS(std,endl);

	for (size_t i(0); i<itsCalibSecDensities.size(); i++)
	{
		os << indent << "*********************************************\n";
		os << indent << "----> CALIB SECURITY #"<< i + 1 << " <----\n";
		os << indent << "*********************************************\n";
		os << itsCalibSecDensities[i]->toString(indent + "    ", nextIndent) << CC_NS(std,endl);
	}

	if ( !GetPortfolio().IsNull() )
	{
		os << indent << "BOOTSTRAP PORTFOLIO : " << CC_NS(std,endl);
		os << GetPortfolio()->toString() << CC_NS(std,endl);
	}
	else
	{
		os << indent << "BOOTSTRAP PORTFOLIO : NONE" << CC_NS(std,endl);
	}

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: ModelFitterName
///	Returns: string
///	Action : try to guess... Tu peux le faire
////////////////////////////////////////////////////
string ARM_NumericalModelFitter::ModelFitterName() const
{
	return "ARM_NumericalModelFitter";
}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Calibrate
///	Returns: void 
///	Action : Does the calibration...
////////////////////////////////////////////////////
void ARM_NumericalModelFitter::Calibrate()
{
	/// -------------------------------------------------------
	/// There's a PF --> bootstrap on numerical PDE calibration
	/// -------------------------------------------------------
	if ( !GetPortfolio().IsNull() )
	{
		BootstrapCalibration();
		return;
	}

	/// --------------------------------------------------
	/// No PF --> no boostrap, numerical calibration only
	/// --------------------------------------------------
	
	/// Gets the pricingmodel to calibrate
	ARM_PricingModel* model = GetPricingModel();

	/// Sets the calibration stuff
	model->setCalibrationStatus(true);
	model->setNumericalModelFitter(this);

	/// Inits the model
	if( itsCalibStates == ARM_PricingStatesPtr(NULL) )
		itsCalibStates = model->Init("EUR", ARM_TimeInfoPtrVector(0));

	/// resiezs unneeded itsCalibStates stuff 
	itsCalibStates->resizePayoffStatesVector(0);
	itsCalibStates->SetOtherPayoffsFlag(0);

	///
	double asOf = model->GetAsOfDate().GetJulian();
	ARM_GP_Vector* resetDates = itsCalibSchedule->GetResetDates();
	
	// start index
	int begin = resetDates->size() - 1 ;

	/// Backward induction
	for(int i = begin; i>=0 ; i-- )
	{	
			model->Induct( itsCalibStates, (*resetDates)[i] - asOf);
	}

	/// The model is no longer calibrating
	model->setCalibrationStatus( false );
	model->setNumericalModelFitter(NULL);

}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: CalibrateWithBootstrap
///	Returns: void 
///	Action : Does the calibration when a porfolio
///			 is provided (bootstrap on global var)
////////////////////////////////////////////////////
void ARM_NumericalModelFitter::BootstrapCalibration()
{	
	/// sanity check
	if ( GetPortfolio().IsNull() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : CalibrateWithBootstrap cannot be performed : no portfolio provided" );
	
	/// Gets the pricingmodel to calibrate
	ARM_PricingModel* model = GetPricingModel();

	/// oups, it only works for MarkovFunctional ; apologies...
	ARM_MarkovFunctional* mfModel = dynamic_cast< ARM_MarkovFunctional* > (model);
	if (!mfModel)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : CalibrateWithBootstrap requires a MarkovFunctional model" );

	ARM_ModelParamsMF* mfParams     = (ARM_ModelParamsMF*)mfModel->GetModelParams();
	ARM_CurveModelParam* curveParam	= dynamic_cast<ARM_CurveModelParam*>(&mfParams->GetModelParam(ARM_ModelParamType::Volatility));

	if (!curveParam)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : CalibrateWithBootstrap : model param is not a curve" );
	
	/// Sets the calibration stuff
	model->setCalibrationStatus(true);
	model->setNumericalModelFitter(this);

	/// Inits the model
	if( itsCalibStates == ARM_PricingStatesPtr(NULL) )
		itsCalibStates = model->Init("EUR", ARM_TimeInfoPtrVector(0));

	/// resiezs unneeded itsCalibStates stuff 
	itsCalibStates->resizePayoffStatesVector(0);
	itsCalibStates->SetOtherPayoffsFlag(0);

	///
	double asOf = model->GetAsOfDate().GetJulian();
	ARM_GP_Vector* resetDates = itsCalibSchedule->GetResetDates();

	/// Get pde scheme
	ARM_NumMethodPtr numMethod = model->GetNumMethod();
	if( numMethod == ARM_NumMethodPtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set!");

	/// check that num method is a PDE
	ARM_PDEMethod* pdeMethod = dynamic_cast<ARM_PDEMethod*>(&*numMethod);
	if (!pdeMethod)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CalibrateWithBootstrap : calibration requires a PDE method");

	ARM_PDENumericalScheme* pdeScheme = pdeMethod->GetPDENumericalScheme();

	/// Initialise model param
	Initialise();

	/// Calibration securitities in terms of resetDates indexes
	ARM_IntVectorPtr portfolioIdx = ComputePortfolioSecuritiesIdx(); 
	size_t sizePf = GetPortfolio()->size();

	/// get time steps from num method
	ARM_GP_VectorPtr timeSteps = ARM_GP_VectorPtr((ARM_GP_Vector*)numMethod->GetTimeSteps()->Clone());
	
	/// please, no instrument on last reset date !
	if ((*portfolioIdx)[sizePf-1] == resetDates->size() - 1) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BootstrapCalibration : cannot boostrap calib on last reset date");
	
	/// Backward induction down to last Pf reset date
	for(int i = resetDates->size() - 1 ; i>(*portfolioIdx)[sizePf-1] ; i-- )
		model->Induct( itsCalibStates, (*resetDates)[i] - asOf);
	
	/// initial guess & bounds
	double lambda = mfParams->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
	double T   = (*resetDates)[(*portfolioIdx)[sizePf-1]] - asOf ;
	double var = mfParams->StateLocalVariance(0, T, T);
	double initialGuess;
	
	if (fabs(lambda)<1e-6)
	{
		if (var>=0)
			initialGuess = sqrt(var / (T/K_YEAR_LEN) );
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BootstrapCalibration : unresolved variance problem");
	}
	else
	{
		initialGuess = var * 2.0 * lambda;
		initialGuess /= 1.- exp( -2. * lambda * T / K_YEAR_LEN ) ;
		if (initialGuess>=0)
			initialGuess  = sqrt(initialGuess);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BootstrapCalibration : unresolved variance problem");
	}


	double lowerbound   = curveParam->GetLowerBound()->Elt(0);
	double upperbound   = curveParam->GetUpperBound()->Elt(0);

	/// Bootstrap on Numerical calibration
	for(i = sizePf-1 ; i>=0 ; i-- )
	{
		/// define indexes in terms of resetDates
		int nextIndex  = (i == sizePf-1) ? (*portfolioIdx)[sizePf-1] + 1 : (*portfolioIdx)[i+1];
		int curIndex   = (*portfolioIdx)[i];

		if ( fabs( GetPortfolio()->GetWeights()->Elt(i) ) > K_DOUBLE_TOL )
		{
			/// Induct down to first following reset date
			for (int j = nextIndex-1; j>=curIndex+1; j--)
			{
				double startTime = (*resetDates)[j] - asOf;
				double endTime   = (*resetDates)[j+1] - asOf;
				pdeScheme->BuildTransitionMatrixes(*model, timeSteps, itsCalibStates, &startTime, &endTime);
				model->Induct( itsCalibStates, (*resetDates)[j] - asOf);
			}

			/// Build vanilla arg
			ARM_VanillaArg* arg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(GetPortfolio()->GetAsset(i), asOf );
			arg->SetIndex(i);
			arg->SetMktPrice( GetPortfolio()->GetMktPrices()->Elt(i) );
			arg->SetExpiry((*resetDates)[curIndex] - asOf);
			arg->SetTenor(0.0);
			arg->SetCurveName( model->GetModelNamePerFactor( GetFactorNb() ) );

			/// define function to solve
			itsFunction.Init(model, mfParams, arg, pdeScheme, timeSteps, itsCalibStates, resetDates, curIndex);
									
			double vol;
			
			try 
			{	
				/// set initial guess
				itsSolver->setInitialGuess(initialGuess, lowerbound, upperbound);
							
				/// solve
				vol = itsSolver->Solve();

				/// update initial guess
				initialGuess = vol;

				/// don't forget to RAZ model (last target computed has been a derivative)
				itsFunction (vol);
			}
			catch (int code)
			{
				if (code == VAR_SQUEEZE_EXCEPTION)
				{
					itsFunction.setVarianceSqueeze(lowerbound);

					AddWarningMessage( "\n" );
					AddWarningMessage( "        ****************************************************************************************\n" );
					AddWarningMessage( "        *** Bootstrap solution is lower than lowerbound, set the solution to the lower bound ***\n" );
					AddWarningMessage( "        ****************************************************************************************\n" );
					AddWarningMessage( "\n" );
				}
				else
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BootstrapCalibration : unresolved exception");
				}
			}
			
			/// free memory
			delete arg;
		}
		else
		{	/// calibrate functionals without boostrap
			for (int j = nextIndex-1; j>=curIndex; j--)
			{
				double startTime = (*resetDates)[j] - asOf;
				double endTime   = (*resetDates)[j+1] - asOf;
				pdeScheme->BuildTransitionMatrixes(*model, timeSteps, itsCalibStates, &startTime, &endTime);
				model->Induct( itsCalibStates, (*resetDates)[j] - asOf);
			}
		}
	}

	/// Backward induction down to first Pf reset date
	for(i = (*portfolioIdx)[0] - 1 ; i>=0 ; i-- )
	{
		double startTime = (*resetDates)[i] - asOf;
		double endTime   = (*resetDates)[i+1] - asOf;
		pdeScheme->BuildTransitionMatrixes(*model, timeSteps, itsCalibStates, &startTime, &endTime);
		model->Induct( itsCalibStates, (*resetDates)[i] - asOf);
	}
		

	/// The model is no longer calibrating
	model->setCalibrationStatus( false );
	model->setNumericalModelFitter(NULL);

}

////////////////////////////////////////////////////
///	Class  : ARM_Bootstrap1D 
///	Routine: Validate
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_NumericalModelFitter::Validate()
{
	if ( ! GetPortfolio().IsNull() )
	{
		if(!GetPortfolio()->IsSameAssetsName())
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				" It is not authorised to calibrate a portofolio with differents  product Name\
				Classe. In case of need, please, use previous calibration principle or ask R&D team to help you!" );

		if(!GetPortfolio()->IsGrowingByExpiry())
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				" Only portfolio increasing by maturity is avalaible, please advise!" );
    
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: Initialise
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////
void ARM_NumericalModelFitter::Initialise() const
{
	/// change the calibParam to the right type!

	/// we are necesarily calibrating a vol...
	ARM_ModelParam* param = &GetPricingModel()->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility) ;
	
	if ( ARM_CurveModelParam* curveParam = dynamic_cast<ARM_CurveModelParam*>(param) )
	{
		double asOf = GetPricingModel()->GetAsOfDate().GetJulian();
		
		if (!itsCalibSchedule.IsNull())
		{
			ARM_GP_Vector* resetDates     = itsCalibSchedule->GetResetDates();
			ARM_GP_Vector asOfs	(resetDates->size(), asOf);
			ARM_GP_Vector breakPointTimes = *resetDates - asOfs;
			curveParam->UpdateValues(&breakPointTimes);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: ComputePortfolioSecuritiesIdx
///	Returns: ARM_IntVectorPtr 
///	Action : Calibration securitities in terms of 
///			 resetDates indexes
////////////////////////////////////////////////////
ARM_IntVectorPtr ARM_NumericalModelFitter::ComputePortfolioSecuritiesIdx() const 
{
	ARM_GP_Vector* resetDates = itsCalibSchedule->GetResetDates();
	size_t sizeSched  = resetDates->size();
	size_t sizePf	  = GetPortfolio()->size();
	ARM_IntVector* indexes = new ARM_IntVector(sizePf);
	double asOf = GetPricingModel()->GetAsOfDate().GetJulian() ;

	for (size_t i(0); i<sizePf; i++)
	{
		int idx = -1;
		double secResetDate = GetPortfolio()->GetAsset(i)->GetResetDates()->Elt(0);

		/// optimisons un petit peu
		int begin = (i==0) ? 0 : (*indexes)[i-1] ;

		for (int j(0); j<sizeSched; j++)
		{
			if ( fabs( (*resetDates)[j] - secResetDate ) < 0.5 )
			{
				idx = j;
				break;
			}
		}
		
		if (idx == -1)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_NumericalModelFitter : security reset date doesn't match any of schedule reset dates." );
		else
			(*indexes)[i] = idx;
	}

	return ARM_IntVectorPtr( indexes ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: UpdateStates
///	Returns: void 
///	Action : Does part of the calibration...
////////////////////////////////////////////////////
void ARM_NumericalModelFitter::UpdateStates( const ARM_PricingStatesPtr& states, ARM_GP_VectorPtr& probas, size_t toResetIdx ) const
{		
	itsCalibSecDensities[toResetIdx]->UpdateStates( states, probas, toResetIdx);
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: getCalibSecPayDatesRelIndexes
///	Returns: ARM_IntVector
///	Action :
////////////////////////////////////////////////////
const ARM_IntVector& ARM_NumericalModelFitter::getCalibSecPayDatesRelIndexes (size_t resetIndex) const
{		
	return itsCalibSecDensities[resetIndex]->getPayDatesRelIndexes();
}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: getCalibSecPayDatesAbsIndexes
///	Returns: ARM_IntVector
///	Action :
////////////////////////////////////////////////////
const ARM_IntVector& ARM_NumericalModelFitter::getCalibSecPayDatesAbsIndexes (size_t resetIndex) const
{		
	return itsCalibSecDensities[resetIndex]->getPayDatesAbsIndexes();
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalModelFitter
///	Routine: getCalibSecInterestTerms
///	Returns: ARM_GP_Vector
///	Action :
////////////////////////////////////////////////////
const ARM_GP_Vector& ARM_NumericalModelFitter::getCalibSecInterestTerms (size_t resetIndex) const
{		
	return itsCalibSecDensities[resetIndex]->getInterestTerms();
}


/// =================================
/// ***** BootstrapPriceFunc ********
/// =================================

///----------------------------------------------------------------
void
BootstrapPriceFunc::Init (	ARM_PricingModel*			model,
							ARM_ModelParamsMF*			modelParams,
							ARM_VanillaArg*				argVanilla,
							ARM_PDENumericalScheme*		pdeScheme,
							const ARM_GP_VectorPtr&		timeSteps,
							const ARM_PricingStatesPtr&	pricingStates,
							ARM_GP_Vector*				resetDates,
							int							curIndex)

{
	itsModel			= model;
	itsModelParams		= modelParams;
	itsArgVanilla		= argVanilla;
	itsPdeScheme		= pdeScheme;
	itsTimeSteps		= timeSteps;
	itsPricingStates	= pricingStates,
	itsResetDates		= resetDates;
	itsCurIndex			= curIndex;
}

///----------------------------------------------------------------
double BootstrapPriceFunc::operator () (double sigma) const 
{
	/// get as of date
	double asOf = itsModel->GetAsOfDate().GetJulian();
	
	
	/// move vol keeping global var constant
	bool varSqueeze = false;
	itsModelParams->SetVolUpToT1AndFreezeGlobVarUpToT2 ((*itsResetDates)[itsCurIndex] - asOf, 
														(*itsResetDates)[itsCurIndex + 1] - asOf, 
														sigma,
														varSqueeze);

	if( varSqueeze )
	{
		throw VAR_SQUEEZE_EXCEPTION;
	}

	/// startTime and endTime to update impacted PDE transition matrixes
	double startTime = (*itsResetDates)[itsCurIndex]     - asOf;
	double endTime   = (*itsResetDates)[itsCurIndex + 1] - asOf;
	
	/// update transition matrixes in pde scheme
	itsPdeScheme->BuildTransitionMatrixes(*itsModel, itsTimeSteps, itsPricingStates, &startTime, &endTime);

	/// calibrate functionals
	itsModel->Induct( itsPricingStates, (*itsResetDates)[itsCurIndex] - asOf);
	
	/// return target function
	return itsArgVanilla->Price(itsModel) - itsArgVanilla->GetMktPrice();
}


///----------------------------------------------------------------
void BootstrapPriceFunc::setVarianceSqueeze (double lowerBound) const 
{
	/// get as of date
	double asOf = itsModel->GetAsOfDate().GetJulian();
	
	
	/// move vol keeping global var constant
	bool varSqueeze = false;
	itsModelParams->SetVolFromT1toT2AndFreezeGlobVarUpToT2(	(*itsResetDates)[itsCurIndex] - asOf, 
															(*itsResetDates)[itsCurIndex + 1] - asOf, 
															lowerBound,
															varSqueeze);

	if( varSqueeze )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Unresolved variance squeeze");
	}

	/// startTime and endTime to update impacted PDE transition matrixes
	double startTime = (*itsResetDates)[itsCurIndex]     - asOf;
	double endTime   = (*itsResetDates)[itsCurIndex + 1] - asOf;
	
	/// update transition matrixes in pde scheme
	itsPdeScheme->BuildTransitionMatrixes(*itsModel, itsTimeSteps, itsPricingStates, &startTime, &endTime);

	/// calibrate functionals
	itsModel->Induct( itsPricingStates, (*itsResetDates)[itsCurIndex] - asOf);
	
}





CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
