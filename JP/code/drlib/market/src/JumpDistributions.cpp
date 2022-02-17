#include "edginc/config.hpp"
#include "edginc/JumpDistributions.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE



IObject* ExponentialJumpDistribution::defaultConstructor()
{
	return new ExponentialJumpDistribution();
}

CClassConstSP const ExponentialJumpDistribution::TYPE =
CClass::registerClassLoadMethod("ExponentialJumpDistribution",
								typeid(ExponentialJumpDistribution),
								ExponentialJumpDistribution::load);

void ExponentialJumpDistribution::load(CClassSP& clazz)
{
	clazz->setPublic();
	REGISTER(ExponentialJumpDistribution, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(IJumpDistribution);
	EMPTY_SHELL_METHOD(defaultConstructor);
}

double ExponentialJumpDistribution::meanAtTime(double u) const 
{
	int index = lower_bound(pieceT->begin(),pieceT->end(),u) - pieceT->begin();
	if (index==pieceT->size()) --index;
	return mean[index];
}

double ExponentialJumpDistribution::simulateOneJump (double jumpTime, 
													 Sequence *rng, 
													 long& iPath) const
{
	rng->populateVector(iPath);
	iPath++;
	return -meanAtTime(jumpTime) * log( *(rng->getVector()) );
}

double ExponentialJumpDistribution::laplace(int index, double u, double x) const
{
	int _index = lower_bound(pieceT->begin(),pieceT->end(),u) - pieceT->begin();
	if (_index==pieceT->size()) --_index;
	if (_index!=index)
	{
		index=index;  // for debugging purposes.. to remove
	}
//	return 1.0/(1.0 + meanAtTime(u)*x);
	double value = 1.0 / (1.0 + mean[index]*x);
	return value;

}

void ExponentialJumpDistribution::multiplyJump(int index, double _mult)
{
	mean[index] *= _mult;
}

void ExponentialJumpDistribution::setPiecewiseTimes(DoubleArraySP &_pieceT) 
{ 
	pieceT = _pieceT; 
	if (pieceT->size()!=mean.size())
	{
		double m;
		if (mean.size()==0) m=1.0;
		else m=mean[0];

		mean.resize(_pieceT->size());
		fill(mean.begin(),mean.end(),m);
	}
}


NonParametricDistribution::NonParametricDistribution(DoubleArraySP &_pieceT, DoubleArrayArray &_means, DoubleArrayArray &_weights)
:CObject(TYPE)
{
	initialize(_pieceT,_means,_weights);
}

void NonParametricDistribution::initialize(DoubleArraySP &_pieceT, DoubleArrayArray &_means, DoubleArrayArray &_weights)
{
	pieceT = _pieceT;
	means = _means;
	weights = _weights;
	{
		QLIB_VERIFY(weights.size()>0, "empty data");
		for (int i=0; i<weights.size(); ++i)
		{
			double sumWeight = 0.0;
			for (int j=0; j<weights[i].size(); ++j)
			{
				sumWeight += weights[i][j];
				QLIB_VERIFY(weights[i][j]>=0, "negative weight!");
			}
			QLIB_VERIFY(sumWeight>=0, "non-positive sum of weight!");
			for (int j=0; j<weights[i].size(); ++j)
				weights[i][j] /= sumWeight;
		}
	}
}

IObject* NonParametricDistribution::defaultConstructor()
{
	return new NonParametricDistribution();
}

CClassConstSP const NonParametricDistribution::TYPE =
CClass::registerClassLoadMethod("NonParametricDistribution",
								typeid(NonParametricDistribution),
								NonParametricDistribution::load);

void NonParametricDistribution::load(CClassSP& clazz)
{
	clazz->setPublic();
	REGISTER(NonParametricDistribution, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(IJumpDistribution);
	EMPTY_SHELL_METHOD(defaultConstructor);
}

int NonParametricDistribution::timeInterval(double u) const 
{
	int index = lower_bound(pieceT->begin(),pieceT->end(),u) - pieceT->begin();
	if (index==pieceT->size()) --index;
	return index;
}

double NonParametricDistribution::simulateOneJump (double jumpTime, 
													 Sequence *rng, 
													 long& iPath) const
{
	int index = timeInterval(jumpTime);
	rng->populateVector(iPath);
	iPath++;
	double unif = *(rng->getVector());
	// first define which exponential random variable are we going to use
	double sumWeight = weights[index][0];
	int pos = 0;
	while (unif>sumWeight)
	{
		++pos;
		sumWeight += weights[index][pos];
	}
	rng->populateVector(iPath);
	iPath++;
	unif = *(rng->getVector());
	return -means[index][pos] * log(unif);
}

double NonParametricDistribution::laplace(int index, double u, double x) const
{
	double value = 0;
	for (int i=0; i<weights[index].size(); ++i)
		value += weights[index][i] / (1.0 + means[index][i]*x);
	return value;
}

void NonParametricDistribution::multiplyJump(int index, double _mult)
{
	for (int pos=0; pos<weights[index].size(); ++pos)
		means[index][pos] *= _mult;
}

void NonParametricDistribution::setPiecewiseTimes(DoubleArraySP &_pieceT) 
{ 
	pieceT = _pieceT; 
	if (pieceT->size()!=means.size())
	{
		DoubleArray m(weights.size());
		if (means.size()==0) fill(m.begin(),m.end(),1.0);
		else m=means[0];

		means.resize(_pieceT->size());
		fill(means.begin(),means.end(),m);
	}
}



DRLIB_END_NAMESPACE
