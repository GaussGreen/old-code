/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: triangularmatrix.h,v $
 * Revision 1.1  2003/10/13 07:51:33  jmprie
 * Initial revision
 *
 *
 */

/*! \file triangularmatrix.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPBASE_TRIANGULARMATRIX_H
#define _INGPBASE_TRIANGULARMATRIX_H

#include "env.h"
#include "port.h"
#include "gpmatrix.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_GP_T_TriangularMatrix
/// \brief
/// ARM_GP_T_TriangularMatrix class defines the
/// triangular mapping between correlated
/// state variables and uncorrelated diffused
/// variables :
/// Xk = a(k,1).Z1 + a(k,2).Z2 +...+ a(k,k).Zk
/// a lot of code is put in the header for FAST ACCESS!
//////////////////////////////////////////////

template <typename T> class ARM_GP_T_TriangularMatrix : public ARM_GP_T_Matrix<T>
{
public:
	/// constructor
	/// explicit to avoid stupid conversion!
    inline explicit ARM_GP_T_TriangularMatrix(size_t rowSize=1, double value=0.0)
	:	ARM_GP_T_Matrix<T>(rowSize, rowSize, value) { CheckSquaredMatrix(); }

	inline explicit ARM_GP_T_TriangularMatrix( ARM_GP_T_Matrix<T> const & matrix )
	:	ARM_GP_T_Matrix<T>(matrix) { CheckSquaredMatrix(); }

	/// copy, assignment and destructor
    inline ARM_GP_T_TriangularMatrix(const ARM_GP_T_TriangularMatrix<T>& rhs)
	:	 ARM_GP_T_Matrix<T>( rhs ){}

	inline ARM_GP_T_TriangularMatrix& operator = (const ARM_GP_T_TriangularMatrix& rhs)
	{
		if(this != &rhs)
			ARM_GP_T_Matrix<T>::operator=(rhs);
		return *this;
	}
    
	virtual ~ARM_GP_T_TriangularMatrix(){};

	inline ARM_GP_T_TriangularMatrix& Inverse();
	inline ARM_GP_T_TriangularMatrix<T>& CholeskyDecompose();

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_GP_T_TriangularMatrix<T>(*this);}
};




////////////////////////////////////////////////////
///	Cholesky decomposition on triangular matrix A (since Cholesky is only on
///	symetric symmetric positive definite matrix with real entries
///	Basically decomposes A = L LT 
///	where L is a lower triangular matrix with positive diagonal entries, and LT denotes the transpose of L. 
///	the result is stored in the initial matrix.
////////////////////////////////////////////////////

template <typename T> ARM_GP_T_TriangularMatrix<T>& ARM_GP_T_TriangularMatrix<T>::CholeskyDecompose( )
{
	double sum;

    for(size_t j=0;j<GetRowsNb();++j)
    {
        for(size_t i=0;i<=j;++i)
        {
            sum=0.0;
            for(size_t k=0;k<i;++k)
                sum += Elt(i,k) * Elt(j,k);

            if(i<j)
            {
                if( - K_NEW_DOUBLE_TOL <= Elt(i,i) &&  Elt(i,i) <= K_NEW_DOUBLE_TOL )
					Elt(i,i) = K_NEW_DOUBLE_TOL;

                  // throw Exception(__LINE__, __FILE__, ERR_TOL_PB,
                    //"Negative or null eigen value in variance/covariance matrix");
                Elt(j,i) = (Elt(j,i)-sum)/Elt(i,i);
            }
            else
            {
                sum = Elt(j,j)-sum;
                if(sum < 0.0 )
                {
                    if(sum < - K_NEW_DOUBLE_TOL )
  						//throw Exception(__LINE__, __FILE__, ERR_TOL_PB,
							//"Negative or null eigen value in variance/covariance matrix");
					sum = K_NEW_DOUBLE_TOL;
                }
                Elt(j,j) = sqrt(sum);
            } /// j=i
        } /// loop i<=j
    } /// loop j < RowSize

	return *this;
}



////////////////////////////////////////////////////
/// Inversion of a matrix
////////////////////////////////////////////////////
template <typename T> ARM_GP_T_TriangularMatrix<T>& ARM_GP_T_TriangularMatrix<T>::Inverse( )
{
    int i,l,k,size=GetRowsNb();
    double invDiag,sum;
    for(i=0;i<size;++i)
    {
        invDiag=Elt(i,i);
        if(invDiag < K_DOUBLE_TOL)
           throw Exception(__LINE__, __FILE__, ERR_MATRIX_INVER_PB, "Can't invert this matrix !");

        invDiag	= 1.0/invDiag;
        Elt(i,i)= invDiag;
        for(l=0;l<i;++l)
        {
            sum=0.0;
            for(k=l;k<i;++k)
                sum += Elt(i,k) * Elt(k,l);
            Elt(i,l) = - invDiag*sum;
        }
    }

	return *this;
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

