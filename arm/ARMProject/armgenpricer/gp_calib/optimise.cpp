/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file optimise.cpp
 *
 *  \brief model fitter that optimizes
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/optimise.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/targetfunc.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/nagfunction.h"

/// kernel
//#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_Optimise::ARM_Optimise(const ARM_Optimise& rhs)
:	ARM_OptimiseBase(rhs), itsAlgorithm(rhs.itsAlgorithm)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Optimise:: ~ARM_Optimise()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_Optimise& ARM_Optimise::operator = (const ARM_Optimise& rhs)
{
    if( this != &rhs )
	{
		ARM_OptimiseBase::operator=( rhs );
		itsAlgorithm = rhs.itsAlgorithm;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_Optimise::Clone() const
{
	return new ARM_Optimise(*this); 
}

////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_Optimise::ARM_Optimise( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    const size_t max_iter,
	bool getDetails,
	double tolerance,
	double stepMax,
	bool localSearch,
	size_t factorNb,
	size_t NbCalibIter,
	int algorithm,
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
		localSearch,
		factorNb,
		NbCalibIter,
		targetType ),
	itsAlgorithm((Algorithm)algorithm)
{
	SetTargetFunction();
};

////////////////////////////////////////////////////
///	Routine: WeightedSquaredND
///	Returns: void
///	Action : global function used to call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL WeightedSquaredND(
	Integer n,
	double* x, 
    double* fx,
	double* g,
    Nag_Comm *comm )
{
    ARM_Optimise* Optimise = (ARM_Optimise*)(comm->p);
	std::vector<double> var(n, *x);
	*fx = (*Optimise->GetFunction())(var);
}

////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: OptimizerND
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_Optimise::OptimizerND(
	ARM_GP_Vector* Var,
	ARM_GP_Vector* boundLower,
	ARM_GP_Vector* boundUpper)
{
    double err;
    Nag_BoundType bound;
	bound = Nag_Bounds;
	int NbVar   = Var->size();  
    double* x	= new double [NbVar];
	double* g	= new double [NbVar];
	double* bl	= new double [NbVar];
	double* bu	= new double [NbVar];

	size_t i;
	for (i = 0; i<NbVar; ++i)
	{
		bl[i] = (*boundLower)[i];
		bu[i] = (*boundUpper)[i];
		x[i] = (*Var)[i];
	}

	/// initialise NAG (e04xxc)
    Nag_E04_Opt options;	
	nag_opt_init(&options);

    Nag_Comm comm;
    comm.p = this;
    static NagError fail,fail2;
    fail.print = FALSE;
	
	/// should we get some details!
	StoreDetailsInFileIfRequired(options,fail);

    options.max_iter	= GetMax_Iter();
    options.optim_tol	= GetTolerance();
    options.step_max	= GetStepMax();
    options.local_search= GetLocalSearch();

	/// optimisation with boundary using no derivatives (e04jbc)
    nag_opt_bounds_no_deriv(NbVar,&WeightedSquaredND,bound,bl,bu,x, &err, g,&options, &comm, &fail);

	ProcessNagWarning(fail);

	/// free the option (e04xzc)
    nag_opt_free(&options,"all",&fail2);

	/// Update all parameters
	ARM_GP_Vector result(NbVar,*x);
    SplitVariables(result);

	/// free memory
    delete g;
    delete bl;
    delete bu;
    delete x;
}

////////////////////////////////////////////////////
///	Routine: ObjectiveFunction
///	Returns: void
///	Action : global function used for the call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL WeightedSquaredNDWithDeriv(
	Integer m,		/// number of subfunctions
	Integer n,		/// number of parameters
	double* x,		/// input
	double* f,		/// output (f(x))
	double* fjac,	/// output  (Df(x,i))
	Integer tdfjac,
	Nag_Comm *comm)	/// context 
{
    ARM_Optimise* Optimise	= (ARM_Optimise*)(comm->p);
    ARM_GP_Vector var(n, x);

	(*Optimise->GetFunction())(var.GetValues(),f,fjac);
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: OptimizerNDWithDeriv
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_Optimise::OptimizerNDWithDeriv(
	ARM_GP_Vector* Var,
	ARM_GP_Vector* boundLower,
	ARM_GP_Vector* boundUpper)
{
	Nag_Comm comm;
	int m = 1;//GetPortfolio()->size();
	int n = Var->size();
	int tdfjac	= n;

	int nclin			= 0;		/// no linear constraints
	int ncnlin			= 0;		/// no non linear constraints
	double *a			= NULL;		/// array of constraints
	int tda				= 0;		///  no linear constraints
	double objf;					/// value of the objective function at the last loop! or value of the sum of residuals 

	int i;
    double* x	= new double [n];			///		input  initial valiues , ouput , result of the optimization
    double* y	= new double [m];			///		input : market prices (y=f(x)
    double* f	= new double [m];			///		output f function values at optimal state 
    double* fjac= new double [m*n];	///		output  Jacobian values	at optimal state	 

	/// boundary for the functions
	double* bl	= new double [n];			///		input :lower bounds of the parameters
	double* bu	= new double [n];			///		input :upper bounds of the parameters
	
	for (i = 0; i<n; ++i)
	{
		bl[i] = (*boundLower)[i];
		bu[i] = (*boundUpper)[i];
		x[i]  = (*Var)[i];
	}

	for (i=0; i<m; ++i)  
		y[i] = 0.0;    /// the target is already taken into account

	/// initialise NAG (e04xxc)
    Nag_E04_Opt options;	
	nag_opt_init(&options);
   
    /// initialise NAG (e04xxc)    
	static NagError fail,fail2;
	INIT_FAIL(fail);
    fail.print = FALSE;
	fail.handler=&NagErrorhandler;
	
	/// should we get some details!
	StoreDetailsInFileIfRequired(options,fail);
	 
	options.max_iter	= GetMax_Iter();
    options.optim_tol	= GetTolerance();
    options.step_max	= GetStepMax();
    options.local_search= GetLocalSearch();

	comm.p = this;
	switch (itsAlgorithm) 
	{
	case NAG_OPT_NLIN_LSQ :
		{
			/// = e04unc  =function to do a non linear least square with boundaries
			nag_opt_nlin_lsq (m, n, nclin, ncnlin, a, tda, bl, bu,y, &WeightedSquaredNDWithDeriv, &NagConfunThrowExceptionIfUsed, x, &objf, f, fjac,tdfjac, &options, &comm, &fail);
			
			break;
		}
	case NAG_OPT_LSQ_DERIV :
		
		{
			/// = e04gbc  =function to do a non linear least square 
			nag_opt_lsq_deriv (m, n, &WeightedSquaredNDWithDeriv,  x, &objf, f, fjac,tdfjac, &options, &comm, &fail);
			
			break;
		}
	case NAG_OPT_LSQ_CHECK_DERIV :
		
		{
			/// = e04yac  = check derivatives
			nag_opt_lsq_check_deriv (m, n, &WeightedSquaredNDWithDeriv,  x,  f, fjac,tdfjac, &comm, &fail);
			
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " OptimizeWithDerivatives : algorithm  bad input :");
			break;
		}
	}

	ProcessNagWarning(fail);

	///free the option (e04xzc)
    nag_opt_free(&options,"all",&fail2);

	/// Update all parameters
	ARM_GP_Vector result(n,x);
    SplitVariables(result);

	/// free memory
    delete bl;
    delete bu;
    delete fjac;
    delete f;
	delete y;
	delete x;
}



////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: Calibrate
///	Returns: 
///	Action :  call nag optimiser and calibrate
////////////////////////////////////////////////////
void ARM_Optimise::Calibrate()
{
    /// get the all parmeters of optimisation multidimen using Nag soft
    ARM_VectorVector variablesvector = MergeVariables();
	switch( itsAlgorithm )
	{
		/// optimizer with no derivatives
		case NAG_OPT_BOUNDS_NO_DERIV:
			OptimizerND(variablesvector[0],variablesvector[1],variablesvector[2]);
			break;

		/// optimizer with derivatives
		case NAG_OPT_NLIN_LSQ: 
		case NAG_OPT_LSQ_DERIV: 
		case NAG_OPT_LSQ_CHECK_DERIV:
			OptimizerNDWithDeriv(variablesvector[0],variablesvector[1],variablesvector[2]);
			break;

		default:
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " unknonw Algorithm!" );
			break;

	}

    for(size_t i=0; i<variablesvector.size(); ++i)
    {
		delete variablesvector[i];
		variablesvector[i]=NULL;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: ModelFitterName
///	Returns: string
///	Action : returns the model fitter name
////////////////////////////////////////////////////
string ARM_Optimise::ModelFitterName() const
{
	return "Optimizer Model Fitter";
}


////////////////////////////////////////////////////
///	Class  : ARM_Optimise
///	Routine: SetTargetFunction
///	Returns: MultiDimFunc*
///	Action : returns the corresponding multi dimension function
////////////////////////////////////////////////////

void ARM_Optimise::SetTargetFunction()
{
	switch( itsAlgorithm )
	{
		/// optimizer with no derivatives
		case NAG_OPT_BOUNDS_NO_DERIV:
			SetFunctionNoClone( new WeightedSquareFunc(GetPricingModel(), this ) );
			break;

		/// optimizer with derivatives
		case NAG_OPT_NLIN_LSQ: 
		case NAG_OPT_LSQ_DERIV: 
		case NAG_OPT_LSQ_CHECK_DERIV:
			SetFunctionNoClone( new MultiDimWithDerivativeFunc(GetPricingModel(), this ) );
			break;

		default:
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " unknonw Algorithm!" );
			break;
	}
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

