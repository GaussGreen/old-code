/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file numfunction.h
 *
 *  \brief template function to compute easily the derivative
 *		of a one dimensional function
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date Ocotber 2004
 */


#ifndef _INGPNUMLIB_NUMFUNCTION_H
#define _INGPNUMLIB_NUMFUNCTION_H

#include "gpbase/port.h"
#include "gpbase/functor.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////
/// TEMPLATE UnaryFuncWithDerivative
/// an unary function that has a method Derivative that 
///		returns a pointor to another unary function
///		that represents its derivative
//////////////////////////////////////////////////////

template<typename T,typename U=T>
	struct UnaryFuncWithDerivative : CC_NS( ARM_GP, UnaryFunc)<T,U>
{
	virtual CC_NS( ARM_GP, UnaryFunc)<T,U>* Derivative() const = 0;
	virtual ~UnaryFuncWithDerivative(){};
};


//////////////////////////////////////////////////////
/// template class to compute automatically the numerical derivative
/// from a function!
//////////////////////////////////////////////////////
template<typename T,typename U=T>
	struct NumDerivativeFunc : public CC_NS( ARM_GP, UnaryFunc)<T,U>
{
public:
    NumDerivativeFunc( const CC_NS( ARM_GP, UnaryFunc)<T,U>& f,double epsilon = 0.00001 )
    :	itsFunction(f),  itsEpsilon(epsilon)
	{}; 

    /// copy constructor
    NumDerivativeFunc( const NumDerivativeFunc<T,U>& rhs )
    :	CC_NS( ARM_GP, UnaryFunc)<T,U>( rhs ), itsFunction( rhs.itsFunction ), itsEpsilon( rhs.itsEpsilon )
    {}

	virtual U operator () ( T _X ) const
    {
        return (itsFunction(_X+itsEpsilon)-itsFunction(_X-itsEpsilon))/(2.0*itsEpsilon);
    }

	virtual ~NumDerivativeFunc(){};

private:    
	const CC_NS( ARM_GP, UnaryFunc)<T,U>& itsFunction;
	double itsEpsilon;
    /// assignment operator (unavailable)
    NumDerivativeFunc<T,U>& operator=( const NumDerivativeFunc<T,U>& rhs );
};


//////////////////////////////////////////////////////
/// template function to create an unary function with a numeric derivative
/// so to use it, give as an input any object that derives already from CC_NS( ARM_GP, UnaryFunc)<T,U>
/// useful for the templated solvers!
//////////////////////////////////////////////////////
template<typename T,typename U=T>
	struct UnaryFuncWithNumDerivative : CC_NS( ARM, UnaryFuncWithDerivative)<T,U>
{
	UnaryFuncWithNumDerivative( const CC_NS( ARM_GP, UnaryFunc)<T,U>& f, double epsilon = 0.00001 )
    :	CC_NS(ARM,UnaryFuncWithDerivative)<T,U>(), itsFunction(f),  itsDerivative( new NumDerivativeFunc<T,U>(f,epsilon) )
	{}

	virtual ~UnaryFuncWithNumDerivative()
	{	delete itsDerivative; }

	virtual U operator () ( T _X ) const
	{	return itsFunction(_X); }


	virtual CC_NS( ARM_GP, UnaryFunc)<T,U>* Derivative() const
	{	return itsDerivative; }

private:
	// It doesn't work we dont know why !!!!!
	UnaryFuncWithNumDerivative( const UnaryFuncWithNumDerivative<T,U>& rhs )
	:	CC_NS(ARM,UnaryFuncWithDerivative)<T,U>(rhs), 
		itsFunction(rhs.itsFunction),  
		itsDerivative( NULL )
	{
		itsDerivative = new NumDerivativeFunc<T,U>(rhs.itsDerivative);
	}

	UnaryFuncWithNumDerivative<T,U>& operator=( const UnaryFuncWithNumDerivative<T,U>& rhs )
	{
		if( this != &rhs )
		{
		///	CC_NS( ARM, UnaryFuncWithDerivative)<T,U>::operator(rhs);
			itsFunction = rhs.itsFunction;
			itsDerivative = new NumDerivativeFunc<T,U>(rhs.itsDerivative);
		}
		return *this;
	}

private:
	const CC_NS( ARM_GP, UnaryFunc)<T,U>& itsFunction;
	CC_NS( ARM_GP, UnaryFunc)<T,U>* itsDerivative;
};

//////////////////////////////////////////////////////
/// functor to return the difference between two values
//////////////////////////////////////////////////////
struct CloneableDbleBinaryFunctor : public DbleBinaryFunctor
{
	virtual CloneableDbleBinaryFunctor* Clone() const = 0;
};

struct ARM_BinFuncMinus : public CloneableDbleBinaryFunctor
{
    double operator()( double x, double y ) const { return x-y; }
	virtual CloneableDbleBinaryFunctor* Clone() const { return new ARM_BinFuncMinus(*this); }
};

//////////////////////////////////////////////////////
/// struct that implements a linear combination of 
/// ARM_VectorPtrDoubleFunc
//////////////////////////////////////////////////////

struct ARM_VectorPtrDbleFuncLinearCombination : public ARM_VectorPtrDbleFunc
{
public: 
	double operator()( ARM_VectorPtr x ) const;
	ARM_VectorPtrDbleFuncLinearCombination( const ARM_VectorPtr& Coefficients, 
		const ARM_VectorPtrDbleFuncPtrVector& FuncVect );
private: 
	ARM_VectorPtr itsCoefficients;
	ARM_VectorPtrDbleFuncPtrVector itsFuncVect;
};

//////////////////////////////////////////////////////
/// struct that implements a multivariable monomial
//////////////////////////////////////////////////////

struct ARM_Monomial : public ARM_VectorPtrDbleFunc
{
public: 
	double operator()( ARM_VectorPtr x ) const;
	ARM_Monomial( const ARM_IntVectorPtr& RiseToPower );
private: 
	ARM_IntVectorPtr itsRiseToPower;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/




