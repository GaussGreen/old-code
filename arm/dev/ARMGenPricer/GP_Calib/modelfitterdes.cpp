/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *  \file modelfitterdes.cpp
 *  \brief file for the definition of templated curves
 *	\author  Amine TRIKI & Arnaud Schauly
 *	\version 1.0
 *	\date October 2004
 */

#include "gpcalib/modelfitterdes.h"
#include "gpcalib/argconvdefault.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include <inst/portfolio.h>
#include "gpbase/gpmatrix.h"

/// gpcalib
#include "gpcalib/bootstrap1d.h"
#include "gpcalib/bootstrapnd.h"
#include "gpcalib/optimise.h"
#include "gpcalib/optimise1d.h"
#include "gpcalib/optimisewithbrent.h"
#include "gpcalib/numerical.h"
#include "gpcalib/vanillaarg.h"

#include "gpmodels/modelfitterhw2f.h"

#include <string>
CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: Constructor
///	Returns: 
///	Action : Optimizer Model Fitter building
////////////////////////////////////////////////////

ARM_ModelFitterDes::ARM_ModelFitterDes(
	ARM_OptimizerType algoType,
	size_t max_iter, 
	double Precision, 
	double stepMax, 
	bool localSearch, 
	bool printLevel)
:
	ARM_RootObject(),
	itsMax_iter(max_iter),
	itsFxTolerance (Precision),
	itsStepMax(stepMax),
	itsLocalSearch(localSearch),
	itsSolverType( ARM_ModelFitterSolverType::NoSolverType ),
	itsOptimizerType ( algoType ),
	itsGetDetails(printLevel)
{
};

	////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: Constructor
///	Returns: 
///	Action : Solver Model Fitter building
////////////////////////////////////////////////////

ARM_ModelFitterDes::ARM_ModelFitterDes(
	ARM_SolverType algoType,	
	size_t max_iter,			
	double fxTol	,			
	double xTol,			
	double gradTol,
	size_t dichoMaxIter,
	double dichoXTol,
	bool printLevel	)	
:
	ARM_RootObject(),
	itsSolverType( algoType ),
	itsMax_iter(max_iter),
	itsFxTolerance (fxTol),
	itsXTolerance(xTol),
	itsGradTolerance(gradTol),
	itsGetDetails(printLevel),
	itsOptimizerType (ARM_ModelFitterOptimizerType::NoOptimizerType ),
	itsDichoMaxIter(dichoMaxIter),
	itsDichoXTol(dichoXTol)

{
};

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////

ARM_Object* ARM_ModelFitterDes::Clone() const
{
	return new ARM_ModelFitterDes( *this ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : copies the additional member variables defined
///	in this object
////////////////////////////////////////////////////

void ARM_ModelFitterDes::CopyNoCleanUp(const ARM_ModelFitterDes& rhs )
{
	itsMax_iter			= rhs.itsMax_iter;
	itsGetDetails		= rhs.itsGetDetails;
	itsFxTolerance		= rhs.itsFxTolerance;
	itsXTolerance       = rhs.itsXTolerance;
	itsGradTolerance    = rhs.itsGradTolerance;
	itsDichoMaxIter		= rhs.itsDichoMaxIter;
	itsDichoXTol		= rhs.itsDichoXTol;
	itsStepMax			= rhs.itsStepMax;
	itsLocalSearch		= rhs.itsLocalSearch;
	itsSolverType		= rhs.itsSolverType;
	itsOptimizerType	= rhs.itsOptimizerType;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: CreateModelFitter
///	Returns: Model Fitter
///	Action : 
////////////////////////////////////////////////////
ARM_ModelFitterPtr ARM_ModelFitterDes::CreateModelFitter(
	ARM_PricingModel* model,
	const ARM_StdPortfolioPtr portfolio, 
	const ARM_ModelParamVector&  calibParams,
	ARM_MethodType methodType,
	ARM_ModelFitterPtr&  linkedmodelfitter,
	ARM_ModelFitterPtr&  previousmodelfitter,
	ARM_MktTargetType targetType,
	size_t FactorNb,
	size_t nbIter,
	ARM_DateStrip* numSchedule,
	const ARM_VanillaSecDensityPtrVector& numSecDensities) const
{
	ARM_ModelFitterPtr result;

	switch(methodType)
	{
	case ARM_CalibMethodType::Bootstrap1D:
		result = ARM_ModelFitterPtr( new ARM_Bootstrap1D( model,
			portfolio, 
			calibParams, 
			linkedmodelfitter, 
			previousmodelfitter, 
			itsMax_iter, 
			itsGetDetails,
			itsFxTolerance,
			itsXTolerance,
			itsDichoMaxIter,
			itsDichoXTol,
			itsSolverType,
			FactorNb, 
			nbIter,
			targetType ) );
		break;
		
	case ARM_CalibMethodType::Optimize:
		if(itsOptimizerType == ARM_ModelFitterOptimizerType::OptimiseWithBrent)
			result = ARM_ModelFitterPtr( new ARM_OptimiseWithBrent(model,
				portfolio, 
				calibParams,
				linkedmodelfitter, 
				previousmodelfitter, 
				itsMax_iter, 
				itsGetDetails, 
				itsFxTolerance, 
				itsStepMax, 
				FactorNb,
				nbIter,
				targetType) );
		else
			result = ARM_ModelFitterPtr( new ARM_Optimise(model,
				portfolio, 
				calibParams,
				linkedmodelfitter,
				previousmodelfitter,
				itsMax_iter, 
				itsGetDetails, 
				itsFxTolerance, 
				itsStepMax,
				itsLocalSearch, 
				FactorNb, 
				nbIter, 
				itsOptimizerType,
				targetType ) );
		break;
		
	case ARM_CalibMethodType::Optimize1D:
		result = ARM_ModelFitterPtr( new ARM_Optimise1D(model,
			portfolio,
			calibParams, 
			linkedmodelfitter, 
			previousmodelfitter, 
			itsMax_iter, 
			itsGetDetails, 
			FactorNb,
			ARM_OptimiseBase::DefaultRelativeTolerance,
			ARM_OptimiseBase::DefaultAbsoluteTolerance, 
			nbIter,
			targetType) );
		break;
		
	case ARM_CalibMethodType::BootstrapND:
		result = ARM_ModelFitterPtr( new ARM_BootstrapND(model,
			portfolio, 
			calibParams, 
			linkedmodelfitter,
			previousmodelfitter, 
			itsMax_iter, 
			itsGetDetails, 
			FactorNb,
			nbIter,
			targetType ) );
		break;

	case ARM_CalibMethodType::Numerical:
		result = ARM_ModelFitterPtr( new ARM_NumericalModelFitter(
													model,
													numSchedule,
													numSecDensities,
													portfolio,
													itsMax_iter, 
													itsGetDetails,
													itsFxTolerance,
													itsXTolerance,
													itsDichoMaxIter,
													itsDichoXTol,
													itsSolverType,
													FactorNb, 
													nbIter,
													targetType ) );
		break;
		
	case ARM_CalibMethodType::HW2FOnly:
		result = ARM_ModelFitterPtr( new ARM_ModelFitterHW2F(model,
			portfolio, 
			calibParams,
			itsMax_iter,
			linkedmodelfitter) );
		break;
		

	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"Calib Method with no calib method type! Cannot create model fitter" );
	}
	
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterDes
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ModelFitterDes::toString( const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "Model Fitter \t:\n";
	if(itsSolverType != ARM_ModelFitterSolverType::NoSolverType)
	{
		os << indent << "Solver Type	= " << ARM_ArgConvReverse_SolverTypeMethod.GetString(itsSolverType)<< "\n";
	}
	else if(itsOptimizerType != ARM_ModelFitterOptimizerType::NoOptimizerType)
	{
		os << indent << "Optimise Type	= " << ARM_ArgConvReverse_OptimizerTypeMethod.GetString(itsOptimizerType)<< "\n";
		os << indent << "Step Max	= " << itsStepMax << "\n";
		os << indent << "LocalSearch?	= " << string(itsLocalSearch? "on" : "off") << "\n";
	}
	os << indent << "Max iteration	= " << itsMax_iter << "\n";
	os << indent << "FxTolerance	= " << itsFxTolerance << "\n";
	if (itsSolverType == ARM_ModelFitterSolverType::NewtonRaphsonWithDichotomy){
		os << indent << "DichoMaxIter	= " << itsDichoMaxIter << "\n";
		os << indent << "DichoXTol		= " << itsDichoXTol << "\n";
	}
	os << indent << "Print Level?	= " << string(itsGetDetails? "on" : "off") << "\n";
	return os.str();
}

CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
