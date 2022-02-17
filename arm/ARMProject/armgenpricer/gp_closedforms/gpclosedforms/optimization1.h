/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  input extender for the closed form framework 
 *
 *	\file input_extender.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_OPTIMIZATION_H
#define _GP_CF_OPTIMIZATION_H

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

CC_BEGIN_NAMESPACE(ARM)


#define ARM_CF_DEFAULT_TOLERANCE 1.0e-9
#define ARM_CF_DEFAULT_MAXIT 2000
////////////////////////////////////////////////////////////////////
/// \class Optimization_Result_Set
/// \brief
/// 
////////////////////////////////////////////////////////////////////

class Optimization_Result_Set
{
	public:
	double OptimalObjective;
	std::vector<double>*	OptimalParamSet;
	Optimization_Result_Set(int n)
	{
		OptimalParamSet= new std::vector<double>(n);
	}
	~Optimization_Result_Set()
	{
		delete OptimalParamSet;
		OptimalParamSet=NULL;
	}

};

class CommSet
{
public:
	void* ptr;			/// used to transmit objective functions
	void* ptr1;			/// used to transmit constraint functions
	std::vector<double> weigthlist;
	std::vector<double> pricelist;
	bool additional_Obj_flag;
	bool remove_price_flag;
	std::vector<double> previousparamlist;
	std::vector<double> paramweightlist;
	CommSet(void* ptr0,
		const std::vector<double>& weigthlist0,
		const std::vector<double>& pricelist0,
		bool removepriceflag)

	:	ptr(ptr0),
		weigthlist(weigthlist0),
		pricelist(pricelist0),
		remove_price_flag(removepriceflag),
		previousparamlist(NULL),
		paramweightlist(NULL),
		additional_Obj_flag(false)
	{}

	CommSet(void* ptr0,void* ptr01,
		const std::vector<double>& weigthlist0,
		const std::vector<double>& pricelist0,
		bool removepriceflag)

	:	ptr(ptr0),ptr1(ptr01),
		weigthlist(weigthlist0),
		pricelist(pricelist0),
		remove_price_flag(removepriceflag),
		previousparamlist(NULL),
		paramweightlist(NULL),
		additional_Obj_flag(false)
	{}


	CommSet(void* ptr0,
		const std::vector<double>& weigthlist0,
		const std::vector<double>& pricelist0,
		bool removepriceflag,
		const std::vector<double>& previousparamlist0,
		const std::vector<double>& paramweightlist0)

	:	ptr(ptr0),
		weigthlist(weigthlist0),
		pricelist(pricelist0),
		previousparamlist(previousparamlist0),
		paramweightlist(paramweightlist0),
		additional_Obj_flag(true)
	{}
	~ CommSet()
	{
	}

};

class Optimization_ObjectiveFuntion 
{
public:
	virtual void NAG_CALL operator() (Integer m, Integer n, 
		double x[], /// input
		double f[],	/// output (f(x))
		double fjac[],  /// output  (Df(x,i))
		Integer tdfjac, Nag_Comm *comm)=0;
	enum Algorithm
	{
		NAG_OPT_NLIN_LSQ,
		NAG_OPT_LSQ_DERIV,
		NAG_OPT_NLIN_LSQ_CONSTRAINT,
		NAG_OPT_LSQ_CHECK_DERIV
	};
	
};

class Optimization_ConstraintFuntion 
{
public:
	virtual void NAG_CALL operator() (Integer m, Integer n, 
		double x[], /// input
		double conf[],	/// output (f(x))
		double cjac[],  /// output  (Df(x,i))
		 Nag_Comm *comm)=0;
	enum Algorithm
	{
		NAG_OPT_NLIN_LSQ,
		NAG_OPT_LSQ_DERIV,
		NAG_OPT_NLIN_LSQ_CONSTRAINT,
		NAG_OPT_LSQ_CHECK_DERIV
	};
	
};


Optimization_Result_Set* OptimizeWithDerivatives(std::vector<double>* strikelist,
								 std::vector<double>* pricelist,
								 std::vector<double>* weightlist,
								 Optimization_ObjectiveFuntion* func,
								 std::vector<double>* initialparamlist,
								 std::vector<double>* lowerboundaryparamlist,
								 std::vector<double>* upperboundaryparamlist,
								 int OptimizationAlgorithm,
								 bool trace_flag,
								 string trace_file,
								 double tolerance = ARM_CF_DEFAULT_TOLERANCE,
								 int max_iter = ARM_CF_DEFAULT_MAXIT
								 );

Optimization_Result_Set* OptimizeWithDerivatives_Constraint(std::vector<double>* strikelist,
								 std::vector<double>* pricelist,
								 std::vector<double>* weightlist,
								 Optimization_ObjectiveFuntion* func,
								 Optimization_ConstraintFuntion* confunc,
								 std::vector<double>* initialparamlist,
								 std::vector<double>* lowerboundaryparamlist,
								 std::vector<double>* upperboundaryparamlist,
								 std::vector<double>* constraintlowerboundaryparamlist,
								 std::vector<double>* constraintupperboundaryparamlist,
								 int algorithm,
								 bool trace_flag,
								 string trace_file,
								 double tolerance,
								 int max_iter,
								 double step
								 );
	
CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

