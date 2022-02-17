/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file nummultidimfunction.h
 *
 *  \brief template function to compute easily the gradient
 *		of multi dimensional function
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date Ocotber 2004
 */


#ifndef _INGPNUMLIB_NUMMULTIDIMFUNCTION_H
#define _INGPNUMLIB_NUMMULTIDIMFUNCTION_H

#include "gpbase/port.h"
#include "gpbase/functor.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gplinalgtypedef.h"

CC_BEGIN_NAMESPACE( ARM )

const double DefaultEpsilon = 1.0e-3;

/// Gradient function from a multi dimensional function
/// a multidimensional function should have the following function
/// size_t resultSize() const; 
/// double operator()( const std::vector<double>& x ) const;

template <typename T>
	struct GradientFunction : public CC_NS( ARM_GP, UnaryFunc)<std::vector<double>,ARM_GP_Matrix> 
{
	GradientFunction( const T& function, double epsilon = DefaultEpsilon )
	:	itsFunction(function), itsEpsilon(epsilon)
	{};

	virtual ~GradientFunction(){};
	virtual ARM_GP_Matrix operator()( const std::vector<double>& x ) const;

private:
	const T& itsFunction;
	double itsEpsilon;
};

template <typename T>
	ARM_GP_Matrix GradientFunction<T>::operator()( const std::vector<double>& x ) const
{
	std::vector<double> tmpX(x), fx_plus, fx_minus;
	ARM_GP_Matrix result(x.size(),itsFunction.resultSize());
	double bump;

	for( size_t i=0; i<x.size(); ++i )
	{
		bump = itsEpsilon*CC_Max(fabs(x[i]),1.0);
		tmpX[i] += bump;
		fx_plus  = itsFunction(tmpX);
		tmpX[i] -= 2*bump;
		fx_minus = itsFunction(tmpX);
		tmpX[i] += bump;

#ifdef __GP_STRICT_VALIDATION
		if( fx_plus.size() !=  itsFunction.resultSize() )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "fx_plus.size() !=  itsFunction.resultSize()" );
		if( fx_minus.size() !=  itsFunction.resultSize() )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "fx_minus.size() !=  itsFunction.resultSize()" );
#endif

		for( size_t j=0; j<itsFunction.resultSize(); ++j )
			result(i,j) = (fx_plus[j]-fx_minus[j])/(2.0*bump);
	}
	return result;
}





CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/




