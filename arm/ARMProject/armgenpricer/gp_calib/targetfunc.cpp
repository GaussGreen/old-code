/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file targetfunc.cpp
 *
 *  \brief file for the various objective or target funciton
 *		in the generic calibration
 *	\author  E.M Ezzine E.Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/targetfunc.h"

/// gpCalib
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/optimisebase.h"
#include "gpcalib/vanillaportfolio.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/gpvector.h"
#include "gpbase/surface.h"

/// kernel
//#include <inst/portfolio.h>

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"


/// STL
#include <algorithm>
#include <functional>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///           FunctionToSolve
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: Constructor
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
FunctionToSolve::FunctionToSolve(ARM_PricingModel* model,
    ARM_ModelFitter* modelFitter,
    ARM_ModelParamType::ParamNb paramType,
    CloneableDbleBinaryFunctor* Func )
:	
    CC_NS( ARM_GP, UnaryFunc)<double,double>(),
    itsPricingModel(model), 
    itsModelFitter(modelFitter),
    itsParamType( paramType),
    itsBinaryFunc(Func),
	itsArgVanilla(NULL)
{}


    
////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: Copy Constructor and assignment operator
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
FunctionToSolve::FunctionToSolve( const FunctionToSolve& rhs )
:
    itsPricingModel(rhs.itsPricingModel ),
    itsModelFitter(  rhs.itsModelFitter ),
    itsParamType(   rhs.itsParamType    ),
    itsBinaryFunc(  rhs.itsBinaryFunc? rhs.itsBinaryFunc->Clone() : NULL	), 
    itsArgVanilla(  rhs.itsArgVanilla   )
{}


FunctionToSolve& FunctionToSolve::operator=( const FunctionToSolve& rhs )
{
    if( this != &rhs )
    {
        itsPricingModel = rhs.itsPricingModel;
        itsModelFitter  = rhs.itsModelFitter;
        itsParamType    = rhs.itsParamType;
	    itsBinaryFunc	= rhs.itsBinaryFunc? rhs.itsBinaryFunc->Clone() : NULL;
        itsBinaryFunc   = rhs.itsBinaryFunc; 
        itsArgVanilla   = rhs.itsArgVanilla;
    }
    return *this;
}

////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
FunctionToSolve::~FunctionToSolve()
{
    delete itsBinaryFunc;
}

////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: InitArgument
///	Returns: 
///	Action : init the various arguments!
////////////////////////////////////////////////////
void FunctionToSolve::InitArgument(ARM_Security* Security, 
								   double price,
								   int index,
								   double expiry,
								   double tenor)
{
    double asOfDate = itsPricingModel->GetAsOfDate().GetJulian();
    itsArgVanilla   = ARM_VanillaArgPtr(ARM_ConverterFromKernel::ConvertSecuritytoArgObject( Security, asOfDate ));
    itsArgVanilla->SetIndex(index);
    itsArgVanilla->SetMktPrice(price);
	
	/// warning, the expiry and the tenor is the expiry of the model parameter to calibrate
	/// it can be different from the expiry of the calibration product
	/// I REPEAT IT CAN BE DIFFERENT FROM THE EXPIRY OF THE CALIBRATION PRODUCT
    itsArgVanilla->SetExpiry(expiry);
    itsArgVanilla->SetTenor(tenor);
	
	itsArgVanilla->SetCurveName( itsPricingModel->GetModelNamePerFactor( itsModelFitter->GetFactorNb() ) );
}


////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: operator()( double x)
///	Returns: 
///	Action : compute the value of the target function!
////////////////////////////////////////////////////

double FunctionToSolve::operator()( double x ) const
{
	/// set to the volatility curve on the index itsArgVanilla->itsIndex the value x!
    itsPricingModel->GetModelParams()->SetModelParamValue( itsParamType, itsArgVanilla->GetIndex(), x, itsArgVanilla->GetExpiry(), itsArgVanilla->GetTenor(), itsModelFitter->GetFactorNb() );
	itsPricingModel->AdviseCurrentCalib(*itsModelFitter);

	if(itsModelFitter->GetLinkedModelFitter() != ARM_ModelFitterPtr(NULL))
		itsModelFitter->GetLinkedModelFitter()->DoCalibrationProcess();
	
	/// compute the binary func with the first argument being the model price and the second one being the market price

	double price = 0.0;
	if(itsModelFitter->GetTargetType() == ARM_CalibrationTarget::PriceTarget)
		price = (*itsBinaryFunc)( itsArgVanilla->Price(itsPricingModel), itsArgVanilla->GetMktPrice() );
	else if(itsModelFitter->GetTargetType() == ARM_CalibrationTarget::ImpliedVolatilityTarget)
		price = (*itsBinaryFunc)( itsArgVanilla->ImpliedVol(itsPricingModel), itsArgVanilla->GetMktPrice() );
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : ARM_CalibrationTarget is not valid" ); 
	return price;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///           FunctionToSolveWithDerivative
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : FunctionToSolveWithDerivative
///	Routine: Constructor
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
FunctionToSolveWithDerivative::FunctionToSolveWithDerivative(ARM_PricingModel* model,
    ARM_ModelFitter* modelFitter,
    ARM_ModelParamType::ParamNb paramType,
    CloneableDbleBinaryFunctor* Func )
:	
	UnaryFuncWithDerivative<double,double>(),
	itsFunctionToSolve(model, modelFitter, paramType, Func )
{
    itsDerivative = new ARM_NumDerivativeDbleToDbleFunctor(*this);
}


    
////////////////////////////////////////////////////
///	Class  : FunctionToSolveWithDerivative
///	Routine: Copy Constructor and assignment operator
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
FunctionToSolveWithDerivative::FunctionToSolveWithDerivative( const FunctionToSolveWithDerivative& rhs )
:
    itsFunctionToSolve(    rhs.itsFunctionToSolve ),
    itsDerivative(  new ARM_NumDerivativeDbleToDbleFunctor( *rhs.itsDerivative ) )
{}


FunctionToSolveWithDerivative& FunctionToSolveWithDerivative::operator=( const FunctionToSolveWithDerivative& rhs )
{
    if( this != &rhs )
    {
        UnaryFuncWithDerivative<double,double>::operator=(rhs );
        itsDerivative   = new ARM_NumDerivativeDbleToDbleFunctor( *rhs.itsDerivative );
    }
    return *this;
}

////////////////////////////////////////////////////
///	Class  : FunctionToSolveWithDerivative
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
FunctionToSolveWithDerivative::~FunctionToSolveWithDerivative()
{
    delete itsDerivative;
}

////////////////////////////////////////////////////
///	Class  : FunctionToSolveWithDerivative
///	Routine: InitArgument
///	Returns: 
///	Action : init the various arguments!
////////////////////////////////////////////////////
void FunctionToSolveWithDerivative::InitArgument( ARM_Security* Security, 
	double price,
	int index,
	double expiry,
	double tenor )
{
    itsFunctionToSolve.InitArgument(Security, price, index, expiry, tenor);
}


////////////////////////////////////////////////////
///	Class  : FunctionToSolveWithDerivative
///	Routine: operator()( double x)
///	Returns: 
///	Action : compute the value of the target function!
////////////////////////////////////////////////////

double FunctionToSolveWithDerivative::operator()( double x ) const
{	
	/// compute the binary func with the first argument being the model price and the second one being the market price
	return itsFunctionToSolve(x);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///           MultiDimFunc            
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : MultiDimFunc 
///	Routine: Constructor
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
MultiDimFunc ::MultiDimFunc (
	ARM_PricingModel* model, 
	ARM_OptimiseBase* optimise )
:  
    itsPricingModel(model),
	itsOptimise(optimise),
	itsPortfolio(NULL)
{
    double asOfDate = model->GetAsOfDate().GetJulian();
	string modelName = model->GetModelNamePerFactor( optimise->GetFactorNb() );
	//itsPortfolio = ARM_ConverterFromKernel::ConvertPortoflio( &*itsOptimise->GetPortfolio(), asOfDate, modelName );
}


////////////////////////////////////////////////////
///	Class  : MultiDimFunc 
///	Routine: copy constructor
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
MultiDimFunc::MultiDimFunc(const MultiDimFunc & rhs )
:
    itsPricingModel(rhs.itsPricingModel ),
	itsOptimise( rhs.itsOptimise		),
	itsPortfolio( rhs.itsPortfolio		)
{}


////////////////////////////////////////////////////
///	Class  : MultiDimFunc 
///	Routine: Assignment operator
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
MultiDimFunc& MultiDimFunc ::operator=(const MultiDimFunc & rhs )
{
    if( this != &rhs )
    {
        itsPricingModel = rhs.itsPricingModel;
		itsOptimise		= rhs.itsOptimise;
		delete itsPortfolio;
		itsPortfolio	= static_cast<ARM_VanillaPortfolio*>( rhs.itsPortfolio->Clone() );
    }
    return *this;
}

////////////////////////////////////////////////////
///	Class  : MultiDimFunc 
///	Routine: Destructor
///	Returns: 
///	Action : Destroy the object!
////////////////////////////////////////////////////
MultiDimFunc :: ~MultiDimFunc ()
{ 
	delete itsPortfolio;
} 

////////////////////////////////////////////////////
///	Class  : MultiDimFunc 
///	Routine: operator()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double MultiDimFunc :: operator()(const std::vector<double>& x ) const
{ ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : unimplemented function" ); }

void MultiDimFunc :: operator()(const std::vector<double>& x, std::vector<double>& fx) const
{ ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : unimplemented function" ); }

void MultiDimFunc :: operator()(const std::vector<double>& x, double* fxVals, double* fxDeriv) const
{ ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : unimplemented function" ); }



////////////////////////////////////////////////////
///	Class  : MultiDimWeightedSquareFunc 
///	Routine: functor
///	Returns: double
///	Action : Calculate the quadratic sum
////////////////////////////////////////////////////
void  MultiDimWeightedSquareFunc::operator()( const std::vector<double>& x, std::vector<double>& result) const
{
    int secIdx, nbSec = GetPortfolio()->size();
	itsOptimise->SplitVariables(ARM_GP_Vector(x));
	
	if(itsOptimise->GetLinkedModelFitter() != ARM_ModelFitterPtr(NULL))
		itsOptimise->GetLinkedModelFitter()->DoCalibrationProcess();

    double fx=0.0, price, wSum=0.0;
    for( secIdx=0; secIdx<nbSec; ++secIdx)
    { 
		itsPricingModel->AdviseCurrentCalibSecIndex(secIdx,*itsOptimise);

		if(itsOptimise->GetTargetType() == ARM_CalibrationTarget::PriceTarget)
			price = GetPortfolio()->GetAsset(secIdx)->Price(itsPricingModel);
        else if(itsOptimise->GetTargetType() == ARM_CalibrationTarget::ImpliedVolatilityTarget)
            price = GetPortfolio()->GetAsset(secIdx)->ImpliedVol(itsPricingModel);
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : ARM_CalibrationTarget is not valid" ); 

		
		price-= GetPortfolio()->GetMktPrice(secIdx);
		price /= GetPortfolio()->GetPrecisiont(secIdx)/0.001;
		double w = GetPortfolio()->GetWeight(secIdx);
		wSum += w;
		fx= w*price*price;
		result[secIdx] = fx;
    }
    
	/// renormalize!
	wSum = fabs(wSum) < K_FRM_TOL ? K_FRM_TOL : wSum;
    for( secIdx=0; secIdx<result.size(); ++secIdx)
		result[secIdx] /= wSum;
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
///           WeightedSquareFunc            
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : WeightedSquareFunc 
///	Routine: functor
///	Returns: double
///	Action : Caloculte the quadratic sum
////////////////////////////////////////////////////
double WeightedSquareFunc ::operator () ( const std::vector<double>& x ) const
{
    int nbSec = GetPortfolio()->size();
	std::vector<double> result(nbSec ,0.0);
	MultiDimWeightedSquareFunc::operator ()( x, result );

	double resultDble = 0;
	for( size_t i=0; i<result.size(); ++i )
		resultDble += result[i];
	return resultDble;
}




////////////////////////////////////////////////////
////////////////////////////////////////////////////
///           MultiDimWithDerivativeFunc            
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : MultiDimWithDerivativeFunc 
///	Routine: functor
///	Returns: double
///	Action : Caloculte the quadratic sum
////////////////////////////////////////////////////
void MultiDimWithDerivativeFunc ::operator () ( const std::vector<double>& x, 
											   double* f, 
											   double* fjac ) const
{
    int secIdx, nbSec = GetPortfolio()->size();
    int nbVar = x.size();
	itsOptimise->SplitVariables(ARM_GP_Vector(x));
    
	if(itsOptimise->GetLinkedModelFitter() != ARM_ModelFitterPtr(NULL))
		itsOptimise->GetLinkedModelFitter()->DoCalibrationProcess();

	ARM_ModelParamVector modelParamVector = itsOptimise->GetCalibParams();
	size_t size = modelParamVector.size();

	double wSum=0.0;
	for( secIdx=0; secIdx<nbSec; ++secIdx)
		wSum += GetPortfolio()->GetWeight(secIdx);
    double price; 
    size_t jacobianIdx=0;
    for( secIdx=0; secIdx<nbSec; ++secIdx)
    { 
		/// Compute price and part in fct to optimise
		itsPricingModel->AdviseCurrentCalibSecIndex(secIdx,*itsOptimise);
        if(itsOptimise->GetTargetType() == ARM_CalibrationTarget::PriceTarget)
            price = GetPortfolio()->GetAsset(secIdx)->Price(itsPricingModel);
        else if(itsOptimise->GetTargetType() == ARM_CalibrationTarget::ImpliedVolatilityTarget)
            price = GetPortfolio()->GetAsset(secIdx)->ImpliedVol(itsPricingModel);
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : ARM_CalibrationTarget is not valid" ); 

		price-= GetPortfolio()->GetMktPrice(secIdx);
		double precision = GetPortfolio()->GetPrecisiont(secIdx);
		double w = sqrt(GetPortfolio()->GetWeight(secIdx)/wSum)/precision;
		
		price*=w;
        f[secIdx] = price;
		
		/// Derivative
		for( size_t paramIdx=0; paramIdx<size; ++paramIdx)
		{
			if (ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(modelParamVector[paramIdx]))
			{
				size_t paramSize = modelParamVector[paramIdx]->size();
				for( size_t varIdx(0); varIdx<paramSize; ++varIdx)
				{
					fjac[jacobianIdx] = GetPortfolio()->GetAsset(secIdx)->Derivative( itsPricingModel,
						*modelParamVector[paramIdx], varIdx,0,0, itsOptimise->GetFactorNb(),itsOptimise->GetTargetType() )*w;

				++jacobianIdx;
				}
			}
			else if (ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(modelParamVector[paramIdx]))
			{
				size_t nbrow = modelParamVector[paramIdx]->rows();
				size_t nbcol = modelParamVector[paramIdx]->cols();
				for( size_t xIndex(0); xIndex<nbrow; ++xIndex)
				{
					for( size_t yIndex(0); yIndex<nbcol; ++yIndex)
					{
						fjac[jacobianIdx] = GetPortfolio()->GetAsset(secIdx)->Derivative( itsPricingModel,
							*modelParamVector[paramIdx], xIndex,yIndex,0, itsOptimise->GetFactorNb(),itsOptimise->GetTargetType() )*w;
					}

				++jacobianIdx;
				}
			}
			else if( ARM_SurfaceListModelParam* surfaceListCalibParam = dynamic_cast<ARM_SurfaceListModelParam*>(modelParamVector[paramIdx]))
				throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "MultiDimWithDerivativeFunc : Unimplemented method for ARM_SurfaceListModelParam !" );
		}


    }
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

