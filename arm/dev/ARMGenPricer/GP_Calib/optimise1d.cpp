/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file optimise1d.cpp
 *
 *  \brief optimiser1d model fitter
 *	\author  E.M Ezzine
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/optimise1d.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparam.h"

/// gpcalib
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"

/// kernel
#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_Optimise1D::ARM_Optimise1D(const ARM_Optimise1D& rhs)
:   ARM_OptimiseBase(rhs),
	itsRelativeTolerance(rhs.itsRelativeTolerance),
	itsAbsoluteTolerance(rhs.itsAbsoluteTolerance)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Optimise1D:: ~ARM_Optimise1D()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Optimise1D& ARM_Optimise1D::operator = (const ARM_Optimise1D& rhs)
{
    if( this != &rhs )
	{
		ARM_OptimiseBase::operator=( rhs );
		itsRelativeTolerance= rhs.itsRelativeTolerance;
		itsAbsoluteTolerance= rhs.itsAbsoluteTolerance;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_Optimise1D::Clone() const
{
	return new ARM_Optimise1D(*this); 
}



////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_Optimise1D::ARM_Optimise1D( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    const size_t max_iter,
	bool getDetails,
	size_t factorNb,
	double relativeTolerance,
	double absoluteTolerance,
	size_t NbCalibIter,
    ARM_MktTargetType  targetType)
:   
    ARM_OptimiseBase(model,portfolio,modelParam,linkedModelFitter,previousModelFitter,max_iter,
		ARM_OptimiseBase::DefaultGetDetails, 
		ARM_OptimiseBase::DefaultTolerance, 
		ARM_OptimiseBase::DefaultStepMax,
		ARM_OptimiseBase::LocalSearch,
		factorNb,
		NbCalibIter,
        targetType),
	itsRelativeTolerance(relativeTolerance),
	itsAbsoluteTolerance(absoluteTolerance)
{
    Validate();  
    SetFunctionNoClone( new WeightedSquareFunc(GetPricingModel(), this ) );
};



////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: Validate
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_Optimise1D::Validate()
{
    if(GetCalibParams().size() != 1)
    {
         CC_Ostringstream os;
            os  << ARM_USERNAME  << ": ARM_Optimise1D has to take a one calib param ";
            throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,  os.str() );
    }
    if(GetCalibParams()[0]->size() != 1)
    {
         CC_Ostringstream os;
            os  << ARM_USERNAME  << ": ARM_Optimise1D has to take a monodimentionnel caliparam";
            throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,  os.str() );
    }
}


////////////////////////////////////////////////////
///	Routine: WeightedSquaredND
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL WeightedSquared1D(double x, double* fx,Nag_Comm *comm )
{
    ARM_Optimise1D* Optimise = (ARM_Optimise1D*)(comm->p);
	ARM_GP_Vector var(1,x);
    *fx = (*Optimise->GetFunction())( var);
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: OptimizerND
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_Optimise1D::Optimizer1D(double Var,
	 double boundLower,
	 double boundUpper)
{
    double err;
    Nag_Comm comm;
    comm.p = this;
    static NagError fail;
	fail.print = FALSE;

    /// Call to NAG 1D optimizer (e04abc)
    nag_opt_one_var_no_deriv(&WeightedSquared1D,itsRelativeTolerance,itsAbsoluteTolerance,
		&boundLower, &boundUpper,GetMax_Iter(),&Var,&err,&comm,&fail);

	/// process nag warning
	ProcessNagWarning(fail);

    /// Update all parameters and Call
    SplitVariables(ARM_GP_Vector(1,Var));
}



////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: Calibrate
///	Returns: 
///	Action :  call nag optimiser and calibrate
////////////////////////////////////////////////////
void ARM_Optimise1D::Calibrate()
{
    /// get the all parmeters of optimisation multidimen using Nag soft
    ARM_VectorVector variablesvector = MergeVariables();
    int size = variablesvector[0]->size();   
    Optimizer1D((*variablesvector[0])[0],(*variablesvector[1])[0],(*variablesvector[2])[0]);    

    for(int i = 0; i < variablesvector.size(); ++i)
    {
		delete variablesvector[i];
		variablesvector[i]=NULL;
	}

}



////////////////////////////////////////////////////
///	Class  : ARM_Optimise1D
///	Routine: ModelFitterName
///	Returns: string
///	Action : returns the model fitter name
////////////////////////////////////////////////////
string ARM_Optimise1D::ModelFitterName() const
{
	return "Optimizer1D Model Fitter";
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

