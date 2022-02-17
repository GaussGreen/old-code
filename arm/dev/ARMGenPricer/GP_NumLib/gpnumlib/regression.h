/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file regression.h
 *
 *  \brief 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#ifndef REGRESSION_H
#define REGRESSION_H

#include "gpbase\typedef.h"
#include "gpnumlib\loess.h"

CC_BEGIN_NAMESPACE( ARM )

struct loess_struct;

// This abstract class 
class ARM_Regression
{
public:
	typedef enum {LS, LOESS} RegressionMode;
	
	virtual ~ARM_Regression() {};

	virtual ARM_GP_VectorPtr ComputeValues(const ARM_GP_MatrixPtr& x) const = 0;
	virtual int GetFactorNb() = 0;
};


class ARM_LSRegression : public ARM_Regression
{
public:
	ARM_LSRegression(const ARM_GP_VectorPtr& y, const ARM_GP_MatrixPtr& x);

	virtual ARM_GP_VectorPtr ComputeValues(const ARM_GP_MatrixPtr& x) const;
	virtual int GetFactorNb() { return itsP; }
	ARM_GP_VectorPtr GetRegressedCoeffs() const { return itsRegressedCoeffs; }

private:
	ARM_GP_VectorPtr itsRegressedCoeffs;
	int itsP;
};

class ARM_LOESSRegression : public ARM_Regression
{
public:
	ARM_LOESSRegression(const ARM_GP_VectorPtr& y, const ARM_GP_MatrixPtr& x, double span);
	virtual ~ARM_LOESSRegression();

	virtual ARM_GP_VectorPtr ComputeValues(const ARM_GP_MatrixPtr& x) const;
	virtual int GetFactorNb() { return itsP; }

private:
	struct loess_struct* itsLoess;
	int itsN;
	int itsP;
	vector<double> itsMinFactors;
	vector<double> itsMaxFactors;
};

CC_END_NAMESPACE()

#endif