/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gpmatrixlinalg.h
 *
 *  \brief gp matrix linar algebra
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPINFRA_GPMATRIXLINALG_H
#define _INGPINFRA_GPMATRIXLINALG_H

#include "port.h"
#include "env.h"

#include "gplinalgtypedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// linear solving of a system a X = y with chek of the determinant
void LinSolve( ARM_GP_Matrix* a, std::vector<double>* y );
/// Same thing, but y is a vector in form of a matrix.
void LinSolve( ARM_GP_Matrix* a, ARM_GP_Matrix* y );

void SingularValuesDecomposition(ARM_GP_Matrix& a,std::vector<double>& sv,ARM_GP_Matrix& v);
void SingularValuesRecomposition(ARM_GP_Matrix& u, std::vector<double>& w, ARM_GP_Matrix& v, const std::vector<double>& b, std::vector<double> *x);

ARM_GP_Matrix*	ACPTransformation(ARM_GP_Matrix* matrix,	ARM_GP_Vector& eigenvalues, int nbFactors = -1);
double ACPTransformationWithRescalling(ARM_GP_Matrix* matrix, ARM_GP_Vector& eigenvalues, ARM_GP_Matrix& ACPMatrix);
ARM_GP_Matrix*	JacobiTransformation(ARM_GP_Matrix& matrix, ARM_GP_Vector& eigenvalues, int& nrot);
ARM_GP_Matrix*	SortedWEigenValues( ARM_GP_Matrix* matrix,	ARM_GP_Vector& eigenvalues);
std::vector<double>*  LeastSquareRegression( const ARM_GP_Matrix& X , const std::vector<double>& Y );
std::vector<double>*  LeastSquareRegressionSVD(const ARM_GP_Matrix& X , const std::vector<double>& Y );

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

