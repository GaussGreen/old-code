//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : sCIDaffineProcesses.hpp
//
//   Description : CIR, CIR2, JumpProcesses
//
//   Date        : Nov 2005
//
//
//----------------------------------------------------------------------------

#ifndef CIR_AFFINE_HPP
#define CIR_AFFINE_HPP

#include "edginc/IAffineProcesses.hpp"
#include "edginc/DoubleMatrix.hpp"	
#include "edginc/Maths.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/GaussianRandomSequence.h"
#include "edginc/RootFinder.hpp"


DRLIB_BEGIN_NAMESPACE


// MOST LIKELY TO BE MOVED BE PUT IN UTIL 
class MARKET_DLL LinearIntegrator
{
public:
	static void createArrayFromInterpolator(const DoubleArray &t,const LinearInterpolantSP &f,DoubleArray &sol);
	static LinearInterpolantSP createInterpolatorFromArray(const DoubleArray &newXaxis,const DoubleArray &x,const DoubleArray &y);
	static LinearInterpolantSP createNewInterpolator( const DoubleArray &newXaxis, const LinearInterpolantSP &f);
	
	static void integrateYZ(
		const DoubleArray &t, 
		const DoubleArray &y, 
		const LinearInterpolantSP &z, 
		double *integral, // (0) integral of y between t0 and t1.. assume t1 increasing and >= t0
		double t0,
		const DoubleArray &t1,
		bool addToPreviousValues = false); 
	static void integrateYZ(
		const DoubleArray &t, 
		const DoubleArray &y, 
		const DoubleArray &z, 
		double *integral, // (0) integral of y between t0 and t1.. assume t1 increasing and >= t0
		double t0,
		const DoubleArray &t1,
		bool addToPreviousValues = false); 
	static void integrateY(
		const DoubleArray &t, 
		const DoubleArray &y, 
		double *integral, // (0) integral of y between t0 and t1.. assume t1 increasing
		double t0,
		const DoubleArray &t1,
		bool addToPreviousValues = false); // assume t1 increasing and greater than t0

private:
	static double integrateHelp(double s, double u, double v, double t, double fs, double ft) 
		// integrate f between u and v, with f linear, and given at time s and t.
	{
		return (v-u)*(fs+(ft-fs)/(t-s)*(0.5*(u+v)-s)); 
	};
	static double integrateHelp(double s, double u, double v, double t, double fs, double gs, double ft, double gt) 
		// integrate f*g between u and v, with f and g linear, and given at time s and t.
	{
		double df = (ft-fs)/(t-s);
		double dg = (gt-gs)/(t-s);
		return 0.166666666666666666666*(v-u)*(6*fs*gs+3*(fs*dg+gs*df)*(u+v-2*s)+2*dg*df*(u*(u+v)+v*v+3*s*(s-u-v)) );
	};
	static inline double integrateHelpSimple(double s, double t, double fs, double ft) 
		// integrate f between s and t, with f linear, and given at time s and t.
	{
		return 0.5*(t-s)*(fs+ft);
	};
	static inline double integrateHelpSimple(double s, double t, double fs, double gs, double ft, double gt) 
		// integrate f*g between s and t, with f and g linear, and given at time s and t.
	{
		return 0.166666666666666666666*(t-s)*(fs*(2*gs+gt)+ft*(gs+2*gt));
	};
};


// These affine processes are  f(t)g(t)lambda_t, where lambda is a classical CIR
// g should be defined on a finer grid as f

class MARKET_DLL CIRAffine: public CObject,
							public virtual IAffineProcess
{
public:
	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;

	CIRAffine();
	CIRAffine(double _lambda0, double _sigma, double _kappa, double _theta);
	virtual ~CIRAffine(){};
 
	void setParameters(double _lambda0, double _sigma, double _kappa, double _theta);

	void setFs(LinearInterpolantSP &_f);

	void setGs(LinearInterpolantSP &_g, size_t pos);
	void setGs(vector<LinearInterpolantSP> &_g);
	void setWhichg(size_t _whichg, bool test = true);

	virtual void   setRungeKuttaDiscretization(double _h) { h=_h; }
	virtual double getRungeKuttaDiscretization() const { return h; };
	virtual void   setPiecewiseTimes(DoubleArraySP &_pieceT) { pieceT = _pieceT; }
	virtual DoubleArraySP getPiecewiseTimes() const { return pieceT; }


	virtual void ricattiVF(int index, double u, double beta, double &fa, double &fb) const; // u belong to pieceT[i], pieceT[i+1]
	double minusLogSurvProba(double bigT, 
							 double futureTime = 0.0, 
							 double lambdaFutureTime = 0.0); 

// find g piecewise linear (discontinuity of df at t), such that E(exp(-\int g_u CIR_u du)) = spt^ratio
// for now assume theta = 1, lambda0 = 1, and makes a rough approximation
	void calibrate(DoubleArray &t, double *spt, double ratio, DoubleArray *newTimes = 0);

	void setDiscretizationTimes(DoubleArray &_discTimes);
	void setDiscretizationTimes(double lastTime, double delta);
	virtual void diffuse(Sequence *rng,
						 long& iPath,
						 double futureTime = 0.0,
						 double lambdaFutureTime = 0.0);
	
	void readBaseProcess(DoubleArray &_baseProcess) { _baseProcess = baseProcess; }
	void readLambdaT(DoubleArray &lambdat) { lambdat = processes[whichg]; }
	void readMinusLogSurvProba(double futureTime, DoubleArray &bigT, double *minusLogsp, bool add);
	void read_F_LambdaT(DoubleArray &lambdat);
	void read_F_MinusLogSurvProba(double futureTime, const DoubleArray &bigT, double *minusLogsp, bool add);

private:
	static void load(CClassSP& clazz);
	static IObject* defaultConstructor();
	// model parameters
	double sigma, kappa, theta, lambda0;
	vector<LinearInterpolantSP> g;    
	LinearInterpolantSP f;

	DoubleArraySP pieceT;
	double h;

	mutable int guessPosg, guessPosf;
	size_t whichg;


	// for mc
	DoubleArray discTimes;
	DoubleArray baseProcess;
	DoubleArrayArray processes;
	DoubleArray expDecay, sqrtInc;
	double thetaStrato; //  theta - sigma^2/(4kappa);
};

DECLARE(CIRAffine);


DRLIB_END_NAMESPACE
#endif
