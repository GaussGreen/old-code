#ifndef __JUMPDISTRIBUTIONS__HPP
#define __JUMPDISTRIBUTIONS__HPP

#include "edginc/IAffineProcesses.hpp"
#include "edginc/Maths.hpp"
#include "edginc/UniformRandomSequence.h"
#include <numeric>


DRLIB_BEGIN_NAMESPACE


class MARKET_DLL ExponentialJumpDistribution: public CObject,
											  public virtual IJumpDistribution
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;
	ExponentialJumpDistribution():CObject(TYPE){}
	ExponentialJumpDistribution(double _mean):
			pieceT(DoubleArraySP(new DoubleArray(1,0.0))),
			mean(1,_mean),
			CObject(TYPE){}
	ExponentialJumpDistribution(DoubleArraySP &_pieceT, DoubleArray &_mean):
			mean(_mean),
			pieceT(_pieceT),
			CObject(TYPE){}
	
	virtual void setPiecewiseTimes(DoubleArraySP &_pieceT);
	virtual DoubleArraySP getPiecewiseTimes() const { return pieceT; }  // nothing to do as there are no time dependency in this example
	virtual void   storeParameters(int index) { meanCopy = mean[index]; }
	virtual void   restoreParameters(int index) { mean[index] = meanCopy; };
	virtual double simulateOneJump (double jumpTime, 
									Sequence *rng, 
									long& iPath) const;
	virtual double laplace(int index, double u, double x) const; // E(exp(-xJ_u))
	virtual void   multiplyJump(int index, double _mult);

private:
	double meanAtTime(double u) const;
	DoubleArraySP pieceT;
	DoubleArray mean;
	double meanCopy;
	static void load(CClassSP& clazz);
	static IObject* defaultConstructor();
};

DECLARE(ExponentialJumpDistribution);

class MARKET_DLL NonParametricDistribution: public CObject,
	public virtual IJumpDistribution
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;
	NonParametricDistribution():CObject(TYPE){}
	NonParametricDistribution(DoubleArraySP &_pieceT, DoubleArrayArray &_means, DoubleArrayArray &_weights);

	void initialize(DoubleArraySP &_pieceT, DoubleArrayArray &_means, DoubleArrayArray &_weights);

	virtual void   setPiecewiseTimes(DoubleArraySP &_pieceT);
	virtual DoubleArraySP getPiecewiseTimes() const { return pieceT; }  // nothing to do as there are no time dependency in this example
	virtual void   storeParameters(int index) { meansCopy = means[index]; }
	virtual void   restoreParameters(int index) { means[index] = meansCopy; };
	void    changeWeights(int index, DoubleArray newWeights) { weights[index] = newWeights; }
	DoubleArrayArray getWeights() { return weights; }
	DoubleArrayArray getMeans() { return means; }

	virtual double simulateOneJump (double jumpTime, 
									Sequence *rng, 
									long& iPath) const;
	virtual double laplace(int index, double u, double x) const; // E(exp(-xJ_u))
	virtual void   multiplyJump(int index, double _mult);

private:
	int timeInterval(double u) const;
	DoubleArraySP pieceT;
	DoubleArrayArray means, weights;
	DoubleArray meansCopy;
	static void load(CClassSP& clazz);
	static IObject* defaultConstructor();
};

DECLARE(NonParametricDistribution);




DRLIB_END_NAMESPACE
#endif // __JUMPDISTRIBUTIONS__HPP
