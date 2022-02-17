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
void LinSolve( ARM_GP_Matrix* a, ARM_GP_Vector* y );
/// Same thing, but y is a vector in form of a matrix.
void LinSolve( ARM_GP_Matrix* a, ARM_GP_Matrix* y );

void SingularValuesDecomposition(ARM_GP_Matrix& a,ARM_GP_Vector& sv,ARM_GP_Matrix& v);
void SingularValuesRecomposition(ARM_GP_Matrix& u, ARM_GP_Vector& w, ARM_GP_Matrix& v, const ARM_GP_Vector& b, ARM_GP_Vector* x);

ARM_GP_Matrix*	ACPTransformation(ARM_GP_Matrix* matrix,	ARM_GP_Vector& eigenvalues, int nbFactors = -1);
double ACPTransformationWithRescalling(ARM_GP_Matrix* matrix, ARM_GP_Vector& eigenvalues, ARM_GP_Matrix& ACPMatrix);
ARM_GP_Matrix*	JacobiTransformation(ARM_GP_Matrix& matrix, ARM_GP_Vector& eigenvalues, int& nrot);
ARM_GP_Matrix*	SortedWEigenValues( ARM_GP_Matrix* matrix,	ARM_GP_Vector& eigenvalues);
ARM_GP_Vector*  LeastSquareRegression( const ARM_GP_Matrix& X , const ARM_GP_Vector& Y );
ARM_GP_Vector*  LeastSquareRegressionSVD(const ARM_GP_Matrix& X , const ARM_GP_Vector& Y );

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

