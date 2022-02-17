/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file regression.cpp
 *
 *  \brief 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#include "gpnumlib/regression.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpVector.h"
#include "gpbase/utilityport.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_LSRegression
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////

ARM_LSRegression::ARM_LSRegression(const ARM_GP_VectorPtr& y, const ARM_GP_MatrixPtr& x)
{
	itsP = x->GetColsNb();
	std::vector<double> * regressedCoeffs = LeastSquareRegression( *x, y->GetValues());
	itsRegressedCoeffs = ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>(*regressedCoeffs));
}

////////////////////////////////////////////////////
///	Class  : ARM_LSRegression
///	Routine: ComputeValues
///	Returns: 
///	Action : Compute values from the regression 
/// functions
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_LSRegression::ComputeValues(const ARM_GP_MatrixPtr& x) const
{
	if (x->GetColsNb() != itsP)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": in the interpolation x don't have the right number of factors." );
	}

	ARM_GP_T_Vector<double>* y = (*x)*(*itsRegressedCoeffs);

	ARM_GP_VectorPtr py( y );

	return py;
}

////////////////////////////////////////////////////
///	Class  : ARM_LOESSRegression
///	Routine: constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////

ARM_LOESSRegression::ARM_LOESSRegression(const ARM_GP_VectorPtr& y, const ARM_GP_MatrixPtr& x, double span)
{
	itsN = x->GetRowsNb();
	itsP = x->GetColsNb();

	double* factors = new double[itsN*itsP];
	double* response = new double[itsN];

	itsMinFactors.resize(itsP);
	itsMaxFactors.resize(itsP);

	if (y->size() != itsN)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": in the regression the response and the factor don't have the same size." );
	}

	int i,j;

	for (i = 0; i < itsP; ++i)
	{
		for (j = 0; j < itsN; ++j)
		{
			factors[i*itsN+j] = (*x)(j,i);

			if (j == 0)
			{
				itsMinFactors[i] = (*x)(j,i);
				itsMaxFactors[i] = (*x)(j,i);
			}
			else
			{
				itsMinFactors[i] = CC_Min((*x)(j,i),itsMinFactors[i]);
				itsMaxFactors[i] = CC_Max((*x)(j,i),itsMaxFactors[i]);
			}
		}
	}

	for (i = 0; i < itsN; ++i)
		response[i] = (*y)[i];

	itsLoess = new loess_struct;

	loess_setup(factors, response, itsN, itsP, itsLoess);
    itsLoess->model.span = span;
    loess(itsLoess);

	delete factors;
	delete response;
}

////////////////////////////////////////////////////
///	Class  : ARM_LOESSRegression
///	Routine: destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////

ARM_LOESSRegression::~ARM_LOESSRegression()
{
	loess_free_mem(itsLoess);
	delete itsLoess;
}

////////////////////////////////////////////////////
///	Class  : ARM_LOESSRegression
///	Routine: ComputeValues
///	Returns: 
///	Action : Compute values from the regression 
/// functions
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_LOESSRegression::ComputeValues(const ARM_GP_MatrixPtr& x) const
{
	int m = x->GetRowsNb();

	if (x->GetColsNb() != itsP)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": in the interpolation x don't have the right number of factors." );
	}

	long se_fit = false;

	double* newFactors = new double[m*itsP];

	int i,j;

	for (i = 0; i < itsP; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			if ( ((*x)(j,i) >= itsMinFactors[i]) && ((*x)(j,i) <= itsMaxFactors[i]) )
			{
				newFactors[i*m+j] = (*x)(j,i);
			}
			else if ((*x)(j,i) < itsMinFactors[i])
			{
				newFactors[i*m+j] = itsMinFactors[i];
			}
			else if ((*x)(j,i) > itsMaxFactors[i])
			{
				newFactors[i*m+j] = itsMaxFactors[i];
			}
		}
	}

	struct pred_struct *loessPred = new pred_struct;

	predict(newFactors, m, itsLoess, loessPred, se_fit);

	ARM_GP_VectorPtr ret(new ARM_GP_T_Vector<double>(m));

	for (i = 0; i < m; ++i)
	{
		ret->Elt(i) = loessPred->fit[i];
	}

	delete newFactors;
	pred_free_mem(loessPred);
	delete loessPred;

	return ret;
}

CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/