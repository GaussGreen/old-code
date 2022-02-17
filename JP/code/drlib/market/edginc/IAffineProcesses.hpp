#ifndef AFFINE_PROCESSES_HPP
#define AFFINE_PROCESSES_HPP

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Sequence.h"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


 

class MARKET_DLL IAffineProcess: public virtual IObject
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;

	virtual void   setRungeKuttaDiscretization(double _h) = 0;
	virtual double getRungeKuttaDiscretization() const = 0;
	virtual void   setPiecewiseTimes(DoubleArraySP &_pieceT) = 0;
	virtual DoubleArraySP getPiecewiseTimes() const = 0;

	/** When implementing a concrete Affine Process, 
	give the Ricatti equation below*/
	virtual void    ricattiVF(int index, double u, double beta, double &fa, double &fb) const = 0; // u belong to pieceT[i], pieceT[i+1]
};

DECLARE(IAffineProcess);

double minusLogSurvProbability(IAffineProcess *ricattiEq,
						 double bigT, 
						 double t, 
						 double lambdat);

void ricatti(IAffineProcess *ricattiEq,
			 int index, 
			 double bigT, 
			 double t,
			 double &alpha, 
			 double &beta,
			 double h);


class MARKET_DLL IJumpDistribution: public virtual IObject
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;

	virtual void   setPiecewiseTimes(DoubleArraySP &_pieceT) = 0;
	virtual DoubleArraySP getPiecewiseTimes() const = 0;
	virtual void   storeParameters(int index) = 0;
	virtual void   restoreParameters(int index) = 0;

	virtual double simulateOneJump (double jumpTime, 
									Sequence *rng, 
									long& iPath) const = 0;
	virtual double laplace(int index, double u, double x) const = 0; // E(exp(-xJ_u))
	virtual void   multiplyJump(int index, double _mult) = 0;
};

DECLARE(IJumpDistribution);


DRLIB_END_NAMESPACE

#endif
