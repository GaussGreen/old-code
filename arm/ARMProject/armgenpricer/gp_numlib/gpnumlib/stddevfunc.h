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
	std::vector<double> itsMeanPerBucket;
	std::vector<double> itsStdDevPerBucket;
	std::vector<double> itsSizePerBucket;
	double itsMean;
	double itsStdDev;
	std::vector<double> itsRunningAverage;
	double itsEmpiricalStdDev;
	int itsTotalSize;
	int itsValuesIndex;
	int itsRunningAverageStep;
	
public:
	// Process method
	virtual void init(int totalSize, int runningAverageStep) = 0;
	virtual void addBucket(const std::vector<double>& vec) = 0;
	virtual void finalize() = 0;
	virtual double ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const = 0;

	virtual ARM_MomentFunc* Clone() const = 0;

	// Accessors
	const ARM_GP_VectorPtr GetMeanPerBucket() const { return ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>(itsMeanPerBucket.begin(),itsMeanPerBucket.end())); };
	const ARM_GP_VectorPtr GetStdDevPerBucket() const { return ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>((itsStdDevPerBucket))); };
	const ARM_GP_VectorPtr GetRunningAverage() const { return ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>((itsRunningAverage))); };
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
	virtual void addBucket(const std::vector<double>& vec);
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
	std::vector<double> itsMeanOddPerBucket;
	std::vector<double> itsMeanEvenPerBucket;
	std::vector<double> itsStdDevOddPerBucket;
	std::vector<double> itsStdDevEvenPerBucket;
	std::vector<double> itsCovPerBucket;
public:
	virtual void init(int totalSize, int runningAverageStep);
	virtual void addBucket(const std::vector<double>& vec);
	virtual void finalize();
	virtual double ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const;
	virtual ARM_MomentFunc* Clone() const { return new ARM_AntitheticMomentFunc(*this); }
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
