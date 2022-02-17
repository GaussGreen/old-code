/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file stddevfunc.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#include "gpnumlib/stddevfunc.h"

#include <cmath>
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_StdStdDev
///	Routine: TestVarianceIsNull
///	Returns: 
///	Action : checks that the variance is null up to the
///				precision of K_NEW_DOUBLE_TOL
///				in strict validation mode, checks that 
///				the variance is really not negative!
////////////////////////////////////////////////////
void ARM_StdDevChecker::TestVarianceIsNull( double& variance ) 
{
	if( fabs(variance) < K_NEW_DOUBLE_TOL*K_NEW_DOUBLE_TOL )
		variance = 0.0;
	if( variance < 0 )
		variance = 0;
	
#if defined(__GP_STRICT_VALIDATION)
	if( variance < -K_NEW_DOUBLE_TOL )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": a variance cannot be negative!" );
#endif
}


////////////////////////////////////////////////////
///	Class  : ARM_StdMomentFunc
///	Routine: operator
///	Returns: 
///	Action : initialize the calculation
////////////////////////////////////////////////////
void ARM_StdMomentFunc::init(int totalSize, int runningAverageStep)
{
	itsTotalSize = totalSize;
	itsMeanPerBucket.resize(0);
	itsStdDevPerBucket.resize(0);
	itsSizePerBucket.resize(0);
	itsMean = 0;
	itsStdDev = 0;
	itsRunningAverage.resize(0);
	itsEmpiricalStdDev = 0;
	itsTotalSize = totalSize;
	itsValuesIndex = 0;
	itsRunningAverageStep = runningAverageStep;
}


////////////////////////////////////////////////////
///	Class  : ARM_StdMomentFunc
///	Routine: addBucket
///	Returns: 
///	Action : add a bucket
////////////////////////////////////////////////////
void ARM_StdMomentFunc::addBucket(const ARM_GP_Vector& vec)
{
	itsMeanPerBucket.resize(itsMeanPerBucket.size()+1);
	itsStdDevPerBucket.resize(itsStdDevPerBucket.size()+1);
	itsSizePerBucket.resize(itsSizePerBucket.size()+1);

	itsSizePerBucket[itsSizePerBucket.size()-1] = vec.size();

	double invN	= 1./itsTotalSize;
	
	for( size_t i=0;i<vec.size(); ++i )
	{
		itsMeanPerBucket[itsMeanPerBucket.size()-1]		+= vec[i]*invN;
		itsStdDevPerBucket[itsStdDevPerBucket.size()-1]	+= vec[i]*vec[i]*invN;
		if ((itsValuesIndex+1)%itsRunningAverageStep == 0)
			itsRunningAverage.push_back((itsMean + itsMeanPerBucket[itsMeanPerBucket.size()-1])/(itsValuesIndex+1)/invN);
		itsValuesIndex++;
	}

	itsMean += itsMeanPerBucket[itsMeanPerBucket.size()-1];
	itsStdDev += itsStdDevPerBucket[itsStdDevPerBucket.size()-1];
}
	
////////////////////////////////////////////////////
///	Class  : ARM_StdMomentFunc
///	Routine: finalize
///	Returns: 
///	Action : finalize the calculation
////////////////////////////////////////////////////
void ARM_StdMomentFunc::finalize()
{
	itsStdDev	-= itsMean*itsMean;
	itsStdDev = itsStdDev/(itsTotalSize-1);
	ARM_StdDevChecker::TestVarianceIsNull(itsStdDev);
	itsStdDev = sqrt(itsStdDev);

	double invN	= 1./itsTotalSize;

	for (size_t i = 0; i < itsSizePerBucket.size(); ++i)
	{
		itsMeanPerBucket[i] /= (itsSizePerBucket[i]*invN);
		itsStdDevPerBucket[i] /= (itsSizePerBucket[i]*invN);

		itsStdDevPerBucket[i]	-= itsMeanPerBucket[i]*itsMeanPerBucket[i];

		itsStdDevPerBucket[i] = itsStdDevPerBucket[i]/itsSizePerBucket[i];
		ARM_StdDevChecker::TestVarianceIsNull(itsStdDevPerBucket[i]);
		itsStdDevPerBucket[i]=sqrt(itsStdDevPerBucket[i]);
	}

	finalizeRunningAverage();

}

////////////////////////////////////////////////////
///	Class  : ARM_StdMomentFunc
///	Routine: ComputeCov
///	Returns: 
///	Action : computes the stdDeviation (std method)
////////////////////////////////////////////////////
double ARM_StdMomentFunc::ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( vec1->size() != vec2->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vec1->size() != vec2->size()!" );
#endif

	double n = vec1->size();
	double mean1	= 0.0,
		   mean2	= 0.0,
		   cov		= 0.0;

	double	invN	= 1./n,
			invN2	= 1./(n-1);

	for( size_t i=0;i<n; ++i )
	{
		mean1 += (*vec1)[i]*invN;
		mean2 += (*vec2)[i]*invN;
		cov   += (*vec1)[i]*(*vec2)[i]*invN;
	}

	cov   -= mean1*mean2;
	cov   *= invN2;
	return cov;
}

////////////////////////////////////////////////////
///	Class  : ARM_StdMomentFunc
///	Routine: finalizeRunningAverage
///	Returns: 
///	Action : finalize the running average calculation
////////////////////////////////////////////////////
void ARM_StdMomentFunc::finalizeRunningAverage()
{
	size_t i;

	int nbEmpiricalStdDevStep = itsRunningAverage.size()/2;
	int beginEmpiricalStdDev = itsRunningAverage.size() - nbEmpiricalStdDevStep;

	itsEmpiricalStdDev = 0.0;

	for(i=beginEmpiricalStdDev; i<itsRunningAverage.size(); ++i )
	{
		itsEmpiricalStdDev += (itsRunningAverage[i]-itsMean)*(itsRunningAverage[i]-itsMean);
	}

	itsEmpiricalStdDev /= nbEmpiricalStdDevStep;
	itsEmpiricalStdDev = sqrt(itsEmpiricalStdDev);
}

////////////////////////////////////////////////////
///	Class  : ARM_AntitheticMomentFunc
///	Routine: init
///	Returns: 
///	Action : initialize the calculation
////////////////////////////////////////////////////
void ARM_AntitheticMomentFunc::init(int totalSize, int runningAverageStep)
{
	ARM_StdMomentFunc::init(totalSize, runningAverageStep);

	itsTotalSize >>= 1;

	itsMeanOdd = 0;
	itsMeanEven = 0;
	itsStdDevOdd = 0;
	itsStdDevEven = 0;
	itsCov = 0;

	itsMeanOddPerBucket.resize(0);
	itsMeanEvenPerBucket.resize(0);
	itsStdDevOddPerBucket.resize(0);
	itsStdDevEvenPerBucket.resize(0);
	itsCovPerBucket.resize(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_AntitheticMomentFunc
///	Routine: addBucket
///	Returns: 
///	Action : add a bucket
///
///	\latexonly
///	An unbiased estimator of $I$ with $N$ trials is defined by: 
///	\[
///	\theta _{N}=\frac{1}{2N}\sum_{i=1}^{N}\left( \psi (x_{i})+\psi
///	(x_{i}^{\prime })\right) 
///	\]%
///	with $x_{i}$ $i.i.d$ to $\mu $. \newline
///	\newline
///	- Variance of the estimator is given as: 
///	\[
///	\sigma _{N}^{2}=\frac{1}{2N}\left( Var[\psi (X)]+Cov[\psi (X),\psi
///	(X^{\prime })]\right) 
///	\]
///	\endlatexonly
////////////////////////////////////////////////////
void ARM_AntitheticMomentFunc::addBucket(const ARM_GP_Vector& vec)
{
	itsMeanPerBucket.resize(itsMeanPerBucket.size()+1);
	itsMeanOddPerBucket.resize(itsMeanOddPerBucket.size()+1);
	itsMeanEvenPerBucket.resize(itsMeanEvenPerBucket.size()+1);
	itsStdDevPerBucket.resize(itsStdDevPerBucket.size()+1);
	itsStdDevOddPerBucket.resize(itsStdDevOddPerBucket.size()+1);
	itsStdDevEvenPerBucket.resize(itsStdDevEvenPerBucket.size()+1);
	itsSizePerBucket.resize(itsSizePerBucket.size()+1);
	itsCovPerBucket.resize(itsCovPerBucket.size()+1);

	itsSizePerBucket[itsSizePerBucket.size()-1] = vec.size() >> 1;

	double invN	= 1./itsTotalSize;
	
	for( size_t i=0;i<vec.size()/2; ++i )
	{
		itsMeanOddPerBucket[itsMeanOddPerBucket.size()-1]		+= vec[i*2]*invN;
		itsStdDevOddPerBucket[itsStdDevOddPerBucket.size()-1]	+= vec[i*2]*vec[i*2]*invN;

		if ((itsValuesIndex+1)%itsRunningAverageStep == 0)
			itsRunningAverage.push_back((itsMeanOdd \
										+ itsMeanEven \
										+ itsMeanOddPerBucket[itsMeanOddPerBucket.size()-1] \
										+ itsMeanEvenPerBucket[itsMeanEvenPerBucket.size()-1])/(itsValuesIndex+1)/invN);

		itsValuesIndex++;

		itsMeanEvenPerBucket[itsMeanEvenPerBucket.size()-1]		+= vec[i*2+1]*invN;
		itsStdDevEvenPerBucket[itsStdDevEvenPerBucket.size()-1]	+= vec[i*2+1]*vec[i*2+1]*invN;
		itsCovPerBucket[itsStdDevEvenPerBucket.size()-1]		+= vec[i*2+1]*vec[i*2]*invN;
		
		if ((itsValuesIndex+1)%itsRunningAverageStep == 0)
			itsRunningAverage.push_back((itsMeanOdd \
										+ itsMeanEven \
										+ itsMeanOddPerBucket[itsMeanOddPerBucket.size()-1] \
										+ itsMeanEvenPerBucket[itsMeanEvenPerBucket.size()-1])/(itsValuesIndex+1)/invN);
		itsValuesIndex++;
	}

	itsMeanOdd += itsMeanOddPerBucket[itsMeanOddPerBucket.size()-1];
	itsMeanEven += itsMeanEvenPerBucket[itsMeanEvenPerBucket.size()-1];
	itsStdDevOdd += itsStdDevOddPerBucket[itsStdDevOddPerBucket.size()-1];
	itsStdDevEven += itsStdDevEvenPerBucket[itsStdDevEvenPerBucket.size()-1];
	itsCov += itsCovPerBucket[itsCovPerBucket.size()-1];
}

////////////////////////////////////////////////////
///	Class  : ARM_AntitheticMomentFunc
///	Routine: finalize
///	Returns: 
///	Action : finalize the calculation
////////////////////////////////////////////////////	
void ARM_AntitheticMomentFunc::finalize()
{
	int n = itsTotalSize;
	
	itsStdDevOdd -= itsMeanOdd*itsMeanOdd;
	itsStdDevEven-= itsMeanEven*itsMeanEven;
	itsCov		  -= itsMeanOdd*itsMeanEven;

	itsMean = (itsMeanOdd + itsMeanEven)/2;

	itsStdDev = (itsStdDevOdd+itsStdDevEven+2.*itsCov)/(4.*itsTotalSize);
	ARM_StdDevChecker::TestVarianceIsNull(itsStdDev);

	itsStdDev =  sqrt(itsStdDev);

	double invN	= 1./itsTotalSize;

	for (size_t i = 0; i < itsSizePerBucket.size(); ++i)
	{
		itsMeanOddPerBucket[i]		/= (itsSizePerBucket[i]*invN);
		itsMeanEvenPerBucket[i]		/= (itsSizePerBucket[i]*invN);
		itsStdDevOddPerBucket[i]	/= (itsSizePerBucket[i]*invN);
		itsStdDevEvenPerBucket[i]	/= (itsSizePerBucket[i]*invN);
		itsCovPerBucket[i]			/= (itsSizePerBucket[i]*invN);

		itsStdDevOddPerBucket[i]	-= itsMeanOddPerBucket[i]*itsMeanOddPerBucket[i];
		itsStdDevEvenPerBucket[i]	-= itsMeanEvenPerBucket[i]*itsMeanEvenPerBucket[i];
		itsCovPerBucket[i]			-= itsMeanOddPerBucket[i]*itsMeanEvenPerBucket[i];

		itsMeanPerBucket[i] = (itsMeanOddPerBucket[i] + itsMeanEvenPerBucket[i])/2;

		itsStdDevPerBucket[i] = (itsStdDevOddPerBucket[i]+itsStdDevEvenPerBucket[i]+2.*itsCovPerBucket[i])/(4.*itsSizePerBucket[i]);
		ARM_StdDevChecker::TestVarianceIsNull(itsStdDevPerBucket[i]);

		itsStdDevPerBucket[i] =  sqrt(itsStdDevPerBucket[i]);
	}

	finalizeRunningAverage();
}

////////////////////////////////////////////////////
///	Class  : ARM_AntitheticStdDev
///	Routine: ComputeCov
///	Returns: 
///	Action : computes the covariance of antithetic variates
///	
/// the formula is
///		cov( Sum(Xi)/2n + Sum(Xbi)/2n, Sum(Yi)/2n + Sum(Ybi)/2n )
///   = 1/4n^2 * (Cov(X,Y)+Cov(Xb,Yb)+Cov(Xb,Y)+Cov(X,Yb) )
///   = 1/4n^2 * ( E[X-Xmean)*(Y-Ymean)] + E[Xb-Xbmean)*(Yb-Ybmean)] + E[Xb-Xbmean)*(Y-Ymean)] + E[X-Xmean)*(Yb-Ybmean)] )
////////////////////////////////////////////////////
double ARM_AntitheticMomentFunc::ComputeCov( ARM_VectorPtr vec1, ARM_VectorPtr vec2 ) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( vec1->size() != vec2->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vec1->size() != vec2->size()!" );
#endif

	/// divides by two
	double n = vec1->size() >> 1;
	double	mean1Odd	= 0.0,
			mean2Odd	= 0.0,
			mean1Even	= 0.0,
			mean2Even	= 0.0,
			covOddOdd	= 0.0,
			covOddEven	= 0.0,
			covEvenOdd	= 0.0,
			covEvenEven	= 0.0,
			invN		= 1./n;

	for( size_t i=0;i<n; ++i )
	{
		mean1Odd   += (*vec1)[i*2]*invN;
		mean1Even  += (*vec1)[i*2+1]*invN;
		mean2Odd   += (*vec2)[i*2]*invN;
		mean2Even  += (*vec2)[i*2+1]*invN;
		covOddOdd  += (*vec1)[i*2]*(*vec2)[i*2]*invN;
		covOddEven += (*vec1)[i*2]*(*vec2)[i*2+1]*invN;
		covEvenOdd += (*vec1)[i*2+1]*(*vec2)[i*2]*invN;
		covEvenEven+= (*vec1)[i*2+1]*(*vec2)[i*2+1]*invN;
	}

	covOddOdd  -= mean1Odd  * mean2Odd;
	covOddEven -= mean1Odd  * mean2Even;
	covEvenOdd -= mean1Even * mean2Odd;
	covEvenEven-= mean1Even * mean2Even;

	double cov = (covOddOdd+covOddEven+covEvenOdd+covEvenEven)/(4.*n);
	return cov;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

