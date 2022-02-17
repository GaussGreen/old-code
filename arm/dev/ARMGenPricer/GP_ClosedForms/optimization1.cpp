/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file Optimization.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpbase/numericconstant.h"

#include "gpclosedforms/gamma.h"
#include "gpclosedforms/optimization1.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////
///
///    Optimization with derivatves and NAG function
///
/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_OptimiseBase
///	Routine: StoreDetailsInFile
///	Returns: void
///	Action : store the details in the file
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Routine: ObjectiveFunction
///	Returns: void
///	Action : global function used for the call Nag optimizer 
////////////////////////////////////////////////////
void NAG_CALL ObjectiveFunction(
	Integer m, Integer n, 
	double x[],		/// input
	double f[],		/// output (f(x))
	double fjac[],  /// output  (Df(x,i))
	Integer tdfjac, Nag_Comm *comm )
{
	int i,j;
	CommSet* cm = (CommSet*) comm->p;
	ARM_GP_Vector weigths(cm->weigthlist);
	ARM_GP_Vector prices(cm->pricelist);
	bool addObj_flag=cm->additional_Obj_flag;
	bool rem_price_flag=cm->remove_price_flag;
	Optimization_ObjectiveFuntion* obj=(Optimization_ObjectiveFuntion*) (cm->ptr);
	(*obj)( m,n,x,f,fjac,tdfjac,comm);

	ARM_GP_Vector* deriv0=new(ARM_GP_Vector)(m*n);
	for(i=0;i<m;i++) for(j=0;j<n;j++) (*deriv0)[i*n+j]=fjac[i*n+j];
	delete deriv0;
	deriv0 =  NULL;
	///

	if(rem_price_flag == TRUE)
	{
		for(i=0;i<m;i++)				/// to remove the objective price if the NAG optimization routine does not do it 
			f[i]-=prices[i];
	}
	
	for(i=0;i<m;i++)					///
	{									///	
		f[i]*=weigths[i];			///		we take into account the weights
		for(j=0;j<n;j++)				///
			fjac[i*n+j]*=weigths[i];	///
	}
	if (addObj_flag==TRUE)				/// if necessairy, we add to the objective  the part associated with the Link
	{
		ARM_GP_Vector param_weigths(cm->paramweightlist);
		ARM_GP_Vector param_previousvalues(cm->previousparamlist);
		int nbp=param_weigths.size();
		if(param_previousvalues.size() != nbp)
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"ObjectiveFunction:  size(param_weigths) # size(param_previousvalues)!" );

		double add_obj=0.0;
		double xx;
		for (i=0;i<nbp;i++)
		{
			xx=(param_previousvalues)[i]-x[i];
			add_obj+=xx*xx*param_weigths[i];
		}
		add_obj /=(double)m;
	}
}

////////////////////////////////////////////////////
///	Routine: confunThrowExceptionIfUsed
///	Returns: void
///	Action : global function used for the call Nag optimizer 
////////////////////////////////////////////////////

void NAG_CALL confunThrowExceptionIfUsed(Integer n, Integer m, Integer needc[], double x[],
	double conf[], double cjac[], Nag_Comm* comm)
{
	throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +": you should not come here!" );
}


////////////////////////////////////////////////////
///	Routine: confunConstraint
///	Returns: void
///	Action : global function used for the call Nag optimizer 
////////////////////////////////////////////////////

void NAG_CALL confunConstraint(Integer n, Integer ncnl, Integer needc[], double x[],
	double conf[], double cjac[], Nag_Comm* comm)
{
	CommSet* cm = (CommSet*) comm->p;
	Optimization_ConstraintFuntion* obj=(Optimization_ConstraintFuntion*) (cm->ptr1);
	(*obj)( ncnl,n,x,conf,cjac,comm);
}


////////////////////////////////////////////////////////////////////
/// \function OptimizeWithDerivatives
/// \brief
/// 
////////////////////////////////////////////////////////////////////

void NAG_CALL errorhandler( char *strng, int code,  char *name)
{
//	if ((code != NE_NOERROR)&&(code != NW_KT_CONDITIONS)&&(code != NW_NOT_CONVERGED))
	{
		string ss1("Error or warning from ");
		string ss2(name);
		string ss3(strng);
		string ss=ss1+ss2+ss3;
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +ss.c_str());
	}
}





Optimization_Result_Set* OptimizeWithDerivatives(ARM_GP_Vector* strikelist,
								 ARM_GP_Vector* pricelist,
								 ARM_GP_Vector* weightlist,
								 Optimization_ObjectiveFuntion* func,
								 ARM_GP_Vector* initialparamlist,
								 ARM_GP_Vector* lowerboundaryparamlist,
								 ARM_GP_Vector* upperboundaryparamlist,
								 int algorithm,
								 bool trace_flag,
								 string trace_file,
								 double tolerance,
								 int max_iter
								 )
{
	
	Nag_Comm comm;
	int m=strikelist->size();			///		input : number of subfunctions = nb of market points

	int n=initialparamlist->size();		///		input : number of parameters
	int tdfjac = n;						///

	int ncnlin = 0;		///
    int nclin  = 0;		///  no linear constraints
	double *a = NULL;	///
	int tda    = 0;		/// 
   
	double objf;		/// value of the objective function at the last loop! or value of the sum of residuals 

	int i;
	double* x   = new double[n];		///		input  initial valiues , ouput , rsult of the optimization
	CC_NS(std,auto_ptr)<double> holdx(x);
    double* y	= new double[m];		///		input : market prices
	CC_NS(std,auto_ptr)<double> holdy(y);
    double* f	= new double[m];		///		output f function values at optimal state 
	CC_NS(std,auto_ptr)<double> holdf(f);
    double* fjac= new double[m*n];		///		output  Jacobian values	at optimal state
	CC_NS(std,auto_ptr)<double> holdfjac(fjac);
    
	
	/// boundary for the functions
	double* bl = new double[n];			///		input :lower bounds of the parameters
	CC_NS(std,auto_ptr)<double> holdbl(bl);
	double* bu = new double[n];			///		input :upper bounds of the parameters
	CC_NS(std,auto_ptr)<double> holdbu(bu);
	for (i=0; i<n; ++i)  
	{
		bl[i] = (*lowerboundaryparamlist)[i];
		bu[i] = (*upperboundaryparamlist)[i];
		x[i] = (*initialparamlist)[i];
    }
	for (i=0; i<m; ++i)  
		y[i] = (*pricelist)[i] * (*weightlist)[i];		/// taking the weights into account for the constants y

	/// initialise NAG Options
    Nag_E04_Opt options;
	nag_opt_init(&options);
	options.max_iter    = max_iter;
	options.optim_tol	= tolerance;

       /// initialise NAG (e04xxc)     
    static NagError fail;
	INIT_FAIL(fail);
    fail.print			= FALSE;
	fail.handler=&errorhandler;

	options.list		= FALSE;
	options.print_level = Nag_NoPrint;
    options.output_level= Nag_NoOutput;
	options.minor_print_level= Nag_NoPrint;
	options.print_deriv=Nag_D_NoPrint;
	options.verify_grad = Nag_NoCheck;
	
	
	strcpy( options.outfile, trace_file.c_str() );
	
		/// Full Print in File of NAG for error tracking
	if(trace_flag==TRUE)
	{
		fail.print			= TRUE;
		options.list		= TRUE;
		options.print_level = Nag_Soln_Iter_Full;
		options.output_level= Nag_MPS_List;
		options.print_deriv=Nag_D_Full;
		options.verify_grad = Nag_SimpleCheck;
	}
	try
	{
		switch (algorithm) 
		{
			/// = e04unc  =function to do a non linear least square with boundaries
		case Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ :
			{
				CommSet cm(func,*weightlist,*pricelist,FALSE);
				comm.p = &cm;
				nag_opt_nlin_lsq (m, n, nclin, ncnlin, a, tda, bl, bu,y, &ObjectiveFunction, &confunThrowExceptionIfUsed, x, &objf, f, fjac,tdfjac, &options, &comm, &fail);
				
				break;
			}
			/// = e04gbc  =function to do a non linear least square 
		case Optimization_ObjectiveFuntion::NAG_OPT_LSQ_DERIV :
			
			{
				CommSet cm(func,*weightlist,*pricelist,TRUE);
				comm.p = &cm;
				nag_opt_lsq_deriv (m, n, &ObjectiveFunction,  x, &objf, f, fjac,tdfjac, &options, &comm, &fail);
				
				break;
			}
			/// = e04yac  = check derivatives
		case Optimization_ObjectiveFuntion::NAG_OPT_LSQ_CHECK_DERIV :
			{
				CommSet cm(func,*weightlist,*pricelist,TRUE);
				comm.p = &cm;
				nag_opt_lsq_check_deriv (m, n, &ObjectiveFunction,  x,  f, fjac,tdfjac, &comm, &fail);
				
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OptimizeWithDerivatives : algorithm  bad input :");
				break;
			}
		}
	}
	catch(Exception&)
	{
	
		if(	fail.code != NE_NOERROR				&& 
			fail.code != NW_NOT_CONVERGED 		&& 
			fail.code != NW_LIN_NOT_FEASIBLE	&&
			fail.code != NW_NONLIN_NOT_FEASIBLE && 
			fail.code != NW_TOO_MANY_ITER		&& 
			fail.code != NW_KT_CONDITIONS		&& 
			fail.code != NW_OVERFLOW_WARN)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OptimizeWithDerivatives : error Nag");
		}
		Optimization_Result_Set* result=new Optimization_Result_Set(n);
		result->OptimalObjective=objf;
		for(i=0;i<n;i++) (*(result->OptimalParamSet))[i]=x[i];
		return (result);

	}
	
	/// free the option (e04xzc)
    nag_opt_free(&options,"all",&fail);

	Optimization_Result_Set* result=new Optimization_Result_Set(n);
	result->OptimalObjective=objf;
	for(i=0;i<n;i++) (*(result->OptimalParamSet))[i]=x[i];
	return (result);
}


Optimization_Result_Set* OptimizeWithDerivatives_Constraint(ARM_GP_Vector* strikelist,
								 ARM_GP_Vector* pricelist,
								 ARM_GP_Vector* weightlist,
								 Optimization_ObjectiveFuntion* func,
								 Optimization_ConstraintFuntion* confunc,
								 ARM_GP_Vector* initialparamlist,
								 ARM_GP_Vector* lowerboundaryparamlist,
								 ARM_GP_Vector* upperboundaryparamlist,
								 ARM_GP_Vector* constraintlowerboundaryparamlist,
								 ARM_GP_Vector* constraintupperboundaryparamlist,
								 int algorithm,
								 bool trace_flag,
								 string trace_file,
								 double tolerance,
								 int max_iter,
								 double progression_step
								 )
{
	
	Nag_Comm comm;
	int m=strikelist->size();			///		input : number of subfunctions = nb of market points

	int n=initialparamlist->size();		///		input : number of parameters

	int nl=constraintlowerboundaryparamlist->size();		///		input : number of nonlinear constraints

	int tdfjac = n;						///

	int ncnlin = nl;	///  nl non linear constraints
    int nclin  = 0;		///  no linear constraints
	double *a = NULL;	///
	int tda    = 0;		/// 
   
	double objf;		/// value of the objective function at the last loop! or value of the sum of residuals 

	int i;
	double* x   = new double[n];		///		input  initial valiues , ouput , rsult of the optimization
	CC_NS(std,auto_ptr)<double> holdx(x);
    double* y	= new double[m];		///		input : market prices
	CC_NS(std,auto_ptr)<double> holdy(y);
    double* f	= new double[m];		///		output f function values at optimal state 
	CC_NS(std,auto_ptr)<double> holdf(f);
    double* fjac= new double[m*n];		///		output  Jacobian values	at optimal state
	CC_NS(std,auto_ptr)<double> holdfjac(fjac);
    
	
	/// boundary for the functions
	double* bl = new double[n+nl];			///		input :lower bounds of the parameters
	CC_NS(std,auto_ptr)<double> holdbl(bl);
	double* bu = new double[n+nl];			///		input :upper bounds of the parameters
	CC_NS(std,auto_ptr)<double> holdbu(bu);
	for (i=0; i<n; ++i)  
	{
		bl[i] = (*lowerboundaryparamlist)[i];
		bu[i] = (*upperboundaryparamlist)[i];
		x[i] = (*initialparamlist)[i];
    }
	for (i=n; i<n+nl; ++i)  
	{
		bl[i] = (*constraintlowerboundaryparamlist)[i-n];///		input :lower bounds of the non linear constraints
		bu[i] = (*constraintupperboundaryparamlist)[i-n];///		input :upper bounds of the non linear constraints
    }
	for (i=0; i<m; ++i)  
		y[i] = (*pricelist)[i] * (*weightlist)[i];		/// taking the weights into account for the constants y

	/// initialise NAG Options
    Nag_E04_Opt options;
	nag_opt_init(&options);
	options.max_iter    = max_iter;
	options.optim_tol	= tolerance;
	options.step_limit  =progression_step;

       /// initialise NAG (e04xxc)     
    static NagError fail;
	INIT_FAIL(fail);
    fail.print			= FALSE;
	fail.handler=&errorhandler;

	options.list		= FALSE;
	options.print_level = Nag_NoPrint;
    options.output_level= Nag_NoOutput;
	options.minor_print_level= Nag_NoPrint;
	options.print_deriv=Nag_D_NoPrint;
	options.verify_grad = Nag_NoCheck;
	
	
	strcpy( options.outfile, trace_file.c_str() );
	
		/// Full Print in File of NAG for error tracking
	if(trace_flag==TRUE)
	{
		fail.print			= TRUE;
		options.list		= TRUE;
		options.print_level = Nag_Soln_Iter_Full;
		options.output_level= Nag_MPS_List;
		options.print_deriv = Nag_D_Full;
		options.verify_grad = Nag_CheckObjCon;
	}
	try
	{
		switch (algorithm) 
		{
			/// = e04unc  =function to do a non linear least square with constraints
		case Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ_CONSTRAINT :
			
			{
				
				CommSet cm(func,confunc,*weightlist,*pricelist,FALSE);
				comm.p = &cm;
				nag_opt_nlin_lsq (m, n, nclin, ncnlin, a, tda, bl, bu,y, &ObjectiveFunction, &confunConstraint, x, &objf, f, fjac,tdfjac, &options, &comm, &fail);
				
				break;
			}

		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OptimizeWithDerivatives_Constraint : algorithm  bad input :");
				break;
			}
		}
	}
	catch(Exception&)
	{
	
		if(	fail.code != NE_NOERROR				&& 
			fail.code != NW_NOT_CONVERGED 		&& 
			fail.code != NW_LIN_NOT_FEASIBLE	&&
			fail.code != NW_NONLIN_NOT_FEASIBLE && 
			fail.code != NW_TOO_MANY_ITER		&& 
			fail.code != NW_KT_CONDITIONS		&& 
			fail.code != NW_OVERFLOW_WARN)
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"OptimizeWithDerivatives_Constraint : error Nag");
		}
		Optimization_Result_Set* result=new Optimization_Result_Set(n);
		result->OptimalObjective=objf;
		for(i=0;i<n;i++) (*(result->OptimalParamSet))[i]=x[i];
		return (result);

	}
	
	/// free the option (e04xzc)
    nag_opt_free(&options,"all",&fail);

	Optimization_Result_Set* result=new Optimization_Result_Set(n);
	result->OptimalObjective=objf;
	for(i=0;i<n;i++) (*(result->OptimalParamSet))[i]=x[i];
	return (result);
}


CC_END_NAMESPACE()


 
#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
