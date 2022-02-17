/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file targetfunc.h
 *
 *  \brief file for the objective or target function
 *		in the calibration problem
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_TARGETFUNC_H
#define _INGPCALIB_TARGETFUNC_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/functor.h"
#include "gpbase/port.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"

#include "typedef.h"
#include <functional>

/// forward declaration
/// kernel object in the global namespace!
class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;
struct ARM_VanillaArg;
struct ARM_VanillaPortfolio;
class FunctionToSolve;
class ARM_OptimiseBase;


///////////////////////////////////////////////////////////////
/// \class FunctionToSolve
/// \brief
///  Function to solve using a solver monodim.
///////////////////////////////////////////////////////////////

class FunctionToSolve : public CC_NS( ARM_GP, UnaryFunc)<double,double>
{
private:
	typedef ARM_ModelParamType::ParamNb ParamType;
    ARM_PricingModel*			itsPricingModel;
    ARM_ModelFitter*			itsModelFitter;
    ParamType					itsParamType;
    CloneableDbleBinaryFunctor*	itsBinaryFunc;
    ARM_VanillaArgPtr			itsArgVanilla;
	
public:
    FunctionToSolve(ARM_PricingModel* model, 
        ARM_ModelFitter* modelFitter		= NULL,
		ParamType	paramType				= ARM_ModelParamType::Volatility, 
        CloneableDbleBinaryFunctor* Func	= new ARM_BinFuncMinus );
    FunctionToSolve( const FunctionToSolve& rhs );
    FunctionToSolve& operator=( const FunctionToSolve& rhs );
	
	/// this is not intented to be a derived class hence the destructor is not private
	/// I repeat THE DESTRUCTOR IS NOT PRIVATE
    ~FunctionToSolve();

    void InitArgument(ARM_Security* Security, 
        double price,
        int index		= 0, 
        double expiry	= 0.0,
		double tenor	= 0.0 );

	ARM_ModelFitter* GetModelFitter() const {return itsModelFitter;};

    inline ARM_VanillaArgPtr& GetVanillaProduct()			{ return itsArgVanilla;			}
    void SetVanillaProduct(ARM_VanillaArgPtr argVanilla)	{ itsArgVanilla = argVanilla;	}
	virtual double operator () ( double x ) const ;

};

///////////////////////////////////////////////////////////////
/// \class FunctionToSolveWithDerivative
/// \brief
///  Function to solve using a solver monodim using derivative.
///////////////////////////////////////////////////////////////

class FunctionToSolveWithDerivative : public UnaryFuncWithDerivative<double,double>
{
private:
	typedef ARM_ModelParamType::ParamNb ParamType;
	FunctionToSolve itsFunctionToSolve;
	ARM_DbleToDbleFunctor* itsDerivative;
	
public:
	FunctionToSolveWithDerivative(ARM_PricingModel* model,
		ARM_ModelFitter* modelFitter		= NULL,
		ParamType	paramType				= ARM_ModelParamType::Volatility, 
        CloneableDbleBinaryFunctor* Func	= new ARM_BinFuncMinus );

    FunctionToSolveWithDerivative( const FunctionToSolveWithDerivative& rhs );
    FunctionToSolveWithDerivative& operator=( const FunctionToSolveWithDerivative& rhs );
    virtual ~FunctionToSolveWithDerivative();

    void InitArgument(ARM_Security* Security, 
        double price,
        int index		= 0, 
        double expiry	= 0.0,
		double tenor	= 0.0 );

	ARM_ModelFitter* GetModelFitter() const {return itsFunctionToSolve.GetModelFitter();};

    inline ARM_VanillaArgPtr& GetVanillaProduct()			{ return itsFunctionToSolve.GetVanillaProduct();	};
    void SetVanillaProduct(ARM_VanillaArgPtr argVanilla)	{ itsFunctionToSolve.SetVanillaProduct(argVanilla);	};
   
	virtual double operator () ( double x ) const ;
    inline virtual ARM_DbleToDbleFunctor* Derivative() const{ return itsDerivative;			}
};


///////////////////////////////////////////////////////////////
/// \class MultiDimFunc
/// \brief
///  base class for all the optimization objective function
///////////////////////////////////////////////////////////////

class MultiDimFunc
{
protected:    
    ARM_PricingModel*	 itsPricingModel;
    ARM_OptimiseBase*    itsOptimise;
	ARM_VanillaPortfolio* itsPortfolio;

public:
    virtual ~MultiDimFunc();
    MultiDimFunc(ARM_PricingModel* model, ARM_OptimiseBase* optimise = NULL );
    MultiDimFunc( const MultiDimFunc& rhs );
    MultiDimFunc& operator=(const MultiDimFunc& rhs );
	ARM_VanillaPortfolio* GetPortfolio() const {return itsPortfolio; }
	
	virtual double operator()(const ARM_GP_Vector& x ) const;
	virtual void operator()(const ARM_GP_Vector& x, ARM_GP_Vector* fx) const;
	virtual void operator()(const ARM_GP_Vector& x, double* f, double* fjac) const;
	virtual MultiDimFunc* Clone() const = 0;
};


///////////////////////////////////////////////////////////////
/// \class MultiDimWeightedSquareFunc
/// \brief
///  Computes the vector of Weight[i]*(MktPrice[i]-ModelPrice[i])^2
///////////////////////////////////////////////////////////////

class MultiDimWeightedSquareFunc: public MultiDimFunc
{
public:
    virtual ~MultiDimWeightedSquareFunc() {};
    MultiDimWeightedSquareFunc(ARM_PricingModel* model, ARM_OptimiseBase* optimise = NULL )
	:	MultiDimFunc( model, optimise ) {}
    MultiDimWeightedSquareFunc( const MultiDimWeightedSquareFunc& rhs ): MultiDimFunc( rhs ) {}
    MultiDimWeightedSquareFunc& operator=(const MultiDimWeightedSquareFunc& rhs )
	{
		if( this != &rhs )
			MultiDimFunc::operator =(rhs);
		return *this;
	}

	virtual void operator()(const ARM_GP_Vector& x, ARM_GP_Vector* fx) const;
	virtual MultiDimFunc* Clone() const{ return new MultiDimWeightedSquareFunc(*this); }
};


///////////////////////////////////////////////////////////////
/// \class WeightedSquareFunc
/// \brief
///  Computes the Sum of Weight[i]*(MktPrice[i]-ModelPrice[i])^2
///////////////////////////////////////////////////////////////
class WeightedSquareFunc : public MultiDimWeightedSquareFunc
{
public:
    virtual ~WeightedSquareFunc() {};
    WeightedSquareFunc(ARM_PricingModel* model,ARM_OptimiseBase* optimiseBase = NULL )
	: MultiDimWeightedSquareFunc( model, optimiseBase ) {}
    WeightedSquareFunc( const WeightedSquareFunc& rhs ): MultiDimWeightedSquareFunc( rhs ) {}
    WeightedSquareFunc& operator=(const WeightedSquareFunc& rhs )
	{
		if( this != &rhs )
			MultiDimWeightedSquareFunc::operator =(rhs);
		return *this;
	}

	virtual double operator()(const ARM_GP_Vector& x ) const;
	virtual MultiDimFunc* Clone() const{ return new WeightedSquareFunc(*this); }
};


///////////////////////////////////////////////////////////////
/// \class MultiDimWeightedSquareFunc
/// \brief
///  Computes the vector of Weight[i]*MktPrice[i] and derivatives
///////////////////////////////////////////////////////////////

class MultiDimWithDerivativeFunc: public MultiDimFunc
{
public:
    virtual ~MultiDimWithDerivativeFunc() {};
    MultiDimWithDerivativeFunc(ARM_PricingModel* model, ARM_OptimiseBase* optimise = NULL )
	:	MultiDimFunc( model, optimise ) {}
    MultiDimWithDerivativeFunc( const MultiDimWeightedSquareFunc& rhs ): MultiDimFunc( rhs ) {}
    MultiDimWithDerivativeFunc& operator=(const MultiDimWithDerivativeFunc& rhs )
	{
		if( this != &rhs )
			MultiDimFunc::operator =(rhs);
		return *this;
	}

	virtual void operator()(const ARM_GP_Vector& x, double* f,double* fjac) const;
	virtual MultiDimFunc* Clone() const{ return new MultiDimWithDerivativeFunc(*this); }
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
