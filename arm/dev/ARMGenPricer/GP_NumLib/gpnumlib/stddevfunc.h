/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file stddevfunc.h
 *
 *  \brief General file for the computation of std deviation
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_STDDEVFUNC_H
#define _INGPNUMLIB_STDDEVFUNC_H

#include "gpbase/port.h"
#include "random.h"
#include "typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_StdDevChecker
{
	static void TestVarianceIsNull( double& Variance );
};


// Those 3 class define the calculation of the 
// expectation and std dev

class ARM_MomentFunc
{
protected:
	ARM_GP_Vector itsMeanPerBucket;
	ARM_GP_Vector itsStdDevPerBucket;
	ARM_GP_Vector itsSizePerBucket;
	double itsMean;
	double itsStdDev;
	ARM_GP_Vector itsRunningAverage;
	double itsEmpiricalStdDev;
	int itsTotalSize;
	int itsValuesIndex;
	int itsRunningAverageStep;
	
public:
	// Process method
	virtual void init(int totalSize, int runningAverageStep) = 0;
	virtual void addBucket(const ARM_GP_Vector& vec) = 0;
	virtual void finalize() = 0;
	virtual double ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const = 0;

	virtual ARM_MomentFunc* Clone() const = 0;

	// Accessors
	const ARM_GP_VectorPtr GetMeanPerBucket() const { return ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(itsMeanPerBucket.Clone())); };
	const ARM_GP_VectorPtr GetStdDevPerBucket() const { return ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(itsStdDevPerBucket.Clone())); };
	const ARM_GP_VectorPtr GetRunningAverage() const { return ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(itsRunningAverage.Clone())); };
	double GetMean() const { return itsMean; };
	double GetStdDev() const { return itsStdDev; };
	double GetEmpiricalStdDev() const const { return itsEmpiricalStdDev; };
};

class ARM_StdMomentFunc : public ARM_MomentFunc
{
protected:
	void finalizeRunningAverage();

public:
	virtual void init(int totalSize, int runningAverageStep);
	virtual void addBucket(const ARM_GP_Vector& vec);
	virtual void finalize();
	virtual double ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const;
	virtual ARM_MomentFunc* Clone() const { return new ARM_StdMomentFunc(*this); }
};

class ARM_AntitheticMomentFunc : public ARM_StdMomentFunc
{
private:
	double itsMeanOdd;
	double itsMeanEven;
	double itsStdDevOdd;
	double itsStdDevEven;
	double itsCov;
	ARM_GP_Vector itsMeanOddPerBucket;
	ARM_GP_Vector itsMeanEvenPerBucket;
	ARM_GP_Vector itsStdDevOddPerBucket;
	ARM_GP_Vector itsStdDevEvenPerBucket;
	ARM_GP_Vector itsCovPerBucket;
public:
	virtual void init(int totalSize, int runningAverageStep);
	virtual void addBucket(const ARM_GP_Vector& vec);
	virtual void finalize();
	virtual double ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const;
	virtual ARM_MomentFunc* Clone() const { return new ARM_AntitheticMomentFunc(*this); }
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
