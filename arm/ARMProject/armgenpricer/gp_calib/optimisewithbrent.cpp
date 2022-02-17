/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file optimisewithbrent.cpp
 *
 *  \brief optimisation with derivation
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/optimisewithbrent.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

///gpnumlib
#include "gpnumlib/optimizer.h"

/// gpcalib
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"

/// kernel
//#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_OptimiseWithBrent::ARM_OptimiseWithBrent(const ARM_OptimiseWithBrent& rhs)
:   ARM_OptimiseBase(rhs)
{}

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_OptimiseWithBrent:: ~ARM_OptimiseWithBrent()
{}

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithDeriv
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_OptimiseWithBrent& ARM_OptimiseWithBrent::operator = (const ARM_OptimiseWithBrent& rhs)
{
    if( this != &rhs )
		ARM_OptimiseBase::operator=( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_OptimiseWithBrent::Clone() const
{
	return new ARM_OptimiseWithBrent(*this); 
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_OptimiseWithBrent::ARM_OptimiseWithBrent( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    const size_t max_iter,
	bool getDetails,
	double tolerance,
	double stepMax,
	size_t factorNb,
	size_t NbCalibIter,
    ARM_MktTargetType  targetType)
: 
    ARM_OptimiseBase(model,
        portfolio,
        modelParam,
        linkedModelFitter,
        previousModelFitter, 
        max_iter,
        getDetails,
        tolerance,
        stepMax,
		ARM_OptimiseBase::LocalSearch,
		factorNb,NbCalibIter,
        targetType)
{
	itsPrecision=tolerance;
	SetFunctionNoClone( new WeightedSquareFunc(GetPricingModel(), this ) );
}



////////////////////////////////////////////////////
///	Class  : WeightedSquareCalculate
///	Routine: Call Brent Optimizer
///	Returns: 
///	Action : optimisation routine
////////////////////////////////////////////////////

WeightedSquareCalculate::WeightedSquareCalculate(const MultiDimFunc* Func)
:	itsFunction ( Func->Clone() ) 
{}

double WeightedSquareCalculate::operator () ( double x ) const
{ return (*itsFunction)( std::vector<double>(1,x) );}


WeightedSquareCalculate::WeightedSquareCalculate(const WeightedSquareCalculate& rhs )
:	itsFunction( rhs.itsFunction? rhs.itsFunction->Clone() : NULL )
{
}

WeightedSquareCalculate& WeightedSquareCalculate::operator=( const WeightedSquareCalculate& rhs )
{
	if( this != & rhs )
	{
		delete itsFunction;
		itsFunction = rhs.itsFunction? rhs.itsFunction->Clone() : NULL;
	}
	return *this;
}

WeightedSquareCalculate::~WeightedSquareCalculate()
{
	delete itsFunction;
}





////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: Call Brent Optimizer
///	Returns: 
///	Action : optimisation routine
////////////////////////////////////////////////////
void ARM_OptimiseWithBrent::Optimizer1DWithBrent(double Var,double boundLower,double boundUpper,double Max_Iter,double Precision)
{
	
	/// Does the optimization using Brent
	WeightedSquareCalculate WSC(GetFunction());
	T_BrentOptimizer<WeightedSquareCalculate> Optimizer( WSC, boundLower, 0.5*(boundLower+boundUpper), boundUpper, Precision, Max_Iter );
	Var = Optimizer.Optimize();

    /// Update all parameters and Call
	SplitVariables(ARM_GP_Vector(1,Var));
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithDeriv
///	Routine: Calibrate
///	Returns: 
///	Action :  call nag optimiser and calibrate
////////////////////////////////////////////////////
void ARM_OptimiseWithBrent::Calibrate()
{
	/// get the all parmeters of optimisation multidimen using Nag soft
    ARM_VectorVector variablesvector = MergeVariables();
    int size = variablesvector[0]->size();   
    Optimizer1DWithBrent((*variablesvector[0])[0],(*variablesvector[1])[0],(*variablesvector[2])[0],GetMax_Iter(),itsPrecision);    

    for(int i = 0; i < variablesvector.size(); ++i)
    {
		delete variablesvector[i];
		variablesvector[i]=NULL;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_OptimiseWithBrent
///	Routine: ModelFitterName
///	Returns: string
///	Action : returns the model fitter name
////////////////////////////////////////////////////
string ARM_OptimiseWithBrent::ModelFitterName() const
{
	return "Brent Optimizer";
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

