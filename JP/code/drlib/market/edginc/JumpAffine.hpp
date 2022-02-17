#ifndef Jump_AFFINE_HPP
#define Jump_AFFINE_HPP

#include "edginc/IAffineProcesses.hpp"
#include "edginc/DoubleMatrix.hpp"	
#include "edginc/Maths.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/UniformRandomSequence.h"
#include "edginc/RootFinder.hpp"
#include <numeric>


DRLIB_BEGIN_NAMESPACE


/** some help to simulate Poisson processes
 if bool is true, simulate conditional on at least one jump in all the Poisson processes up to time T */

class MARKET_DLL PoissonDiffuseHelp
{
public:
	void Initialize(DoubleArray &_avFreq);

	vector<int> nbJumps;
	vector< DoubleArraySP > jumpTimes;

	int GetPoissonNbJumps(double T, 
						  bool cond,
						  double unif,  // unif is the uniform rv which will give the total nb of jumps
						  Sequence *rng, // random number then used to select which Poisson process jumped
						  long& iPath,
						  double t=0.0); // if cond false, nb Jumps conditionned on at least one jump...
	void simulateJumpTimeCondNbJumps(double T, int index, Sequence *rng, long& iPath, double t=0.0);
	double getSumFreq() { return sumFreq; }
	double getFirstJump() { return firstJump; }
private:
	DoubleArray avFreq;
	double sumFreq, firstJump;
	int simulatePoisson(double freq, bool cond, double unif);
	int simulateBinomial(double p, int N, double unif);
};



// defines the jump affine model f(t)*lambda_t
// where d lambda = -k lambda dt + J_tdN_t,
// where J has a given Laplace transform
// and N has a given frequency
// still virtual as we leave free the definition of J
// all the parameters models can very easily be time dependent 
// (for simplicity, at least for now, only J and f are given to be time dependent)
// the discontinuity of f' and J are assumed to be within pieceT

class MARKET_DLL JumpAffine: public CObject,
							 public virtual IAffineProcess
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;

	JumpAffine();
	JumpAffine(double _lambda0, double _kappa, double _freq, IJumpDistributionSP _jump);
	virtual ~JumpAffine(){};

	virtual void setParameters(double _lambda0, double _kappa, double _freq, IJumpDistributionSP _jump);
	void setFs(LinearInterpolantSP &_f);

	virtual void   setRungeKuttaDiscretization(double _h) { h=_h; }
	virtual double getRungeKuttaDiscretization() const { return h; };
	virtual void   setPiecewiseTimes(DoubleArraySP &_pieceT) { pieceT = _pieceT; }
	virtual DoubleArraySP getPiecewiseTimes() const { return pieceT; }

	void setSimulatedJumpTimes(DoubleArraySP &_jumpTimes, int _nbJumps) { jumpTimes = _jumpTimes; nbJumps = _nbJumps;}
	void setNbSimulatedJumpPoisson(int _nbJumps) { nbJumps = _nbJumps; }
	void simulateJumpSizes(UniformRandomSequence *rng,  // conditional on nbJumps
						   long& iPath);

	/**  replace J_t by alpha(t)J so that it calibrate to spt^ratio 
	This is generic, i.e. will work for any J_t (when it can!) */
	void calibrate(DoubleArray &t, double *spt, double ratio, const double MAXMULTJUMP = 200);
	
	double minusLogSurvProba(double bigT, 
		double futureTime = 0.0, 
		double lambdaFutureTime = 0.0); 

	double readLambdaT(double T, 
		 			   double futureTime = 0, 
					   double lambdaFutureTime = 0);
	double readMinusLogSurvProba(bool weKnowJumpSizes,
								 double bigT, 
								 double futureTime = 0, 
								 double lambdaFutureTime = 0);
	
	double read_F_LambdaT(double T, 
						double futureTime = 0, 
						double lambdaFutureTime = 0);
	double read_F_MinusLogSurvProba(bool weKnowJumpSizes,
									double T, 
									double futureTime = 0, 
									double lambdaFutureTime = 0);

	virtual void ricattiVF(int index, double u, double beta, double &fa, double &fb) const; // u belong to pieceT[i], pieceT[i+1]

protected:
	// model parameters
	double kappa, freq, lambda0;
	LinearInterpolantSP f;   
	mutable int guessPosf;

	DoubleArraySP pieceT;
	double h;

	IJumpDistributionSP jump;

	// for mc
	DoubleArray jumpSizes;
	int nbJumps;
	DoubleArraySP jumpTimes;

	double integrateExp(double t);
	double integrateIntegrateExp(double t);
	double integrateExpF(double s, double t);

	double laplace(double u, double x) const;


	static void load(CClassSP& clazz);
	static IObject* defaultConstructor();

};

DECLARE(JumpAffine);


DRLIB_END_NAMESPACE
#endif
