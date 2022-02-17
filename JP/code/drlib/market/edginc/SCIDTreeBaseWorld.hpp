#ifndef ___SCIDTREEBASEWORLD_HPP
#define ___SCIDTREEBASEWORLD_HPP

#include "edginc/CIRAffine.hpp"
#include "edginc/JumpAffine.hpp"
#include "edginc/JumpDistributions.hpp"

#include "edginc/UniformRandomSequence.h"
#include "edginc/GaussianRandomSequence.h"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE


class MARKET_DLL baseWorldSpread // to store spread of one given world in MC simulation
{
public:
	void initialize(int _nbNames, DoubleArray _discTimesSpread) 
	{ 
		discTimesSpread = _discTimesSpread;
		nbNames = _nbNames;
		nbTimes = discTimesSpread.size();
		intensity.resize(nbNames,DoubleArray(nbTimes,0.0)); 
		intensityIntegral.resize(nbNames,DoubleArray(nbTimes,0.0)); 
	};
	DoubleArray discTimesSpread;
	IntArray    discTimesLossIndices;
	DoubleArrayArray intensity, intensityIntegral;  // size nbNames * nbTimes
	int nbNames, nbTimes;
};


class MARKET_DLL TreeBaseWorld: public virtual VirtualDestructorBase
{
public:
	TreeBaseWorld(){};

	void InitializeRNG(long seed);

	void setDiscretizationTimesCF(DoubleArray &_times);
	void setDiscretizationTimesCF(double lastTime, double delta);
	void setDiscretizationTimesIDIO(DoubleArray &_times);
	void setDiscretizationTimesIDIO(double lastTime, double delta);

	void setWorlds(vector<LinearInterpolantSP> _f);
	void setUnderlyingProcesses(const int nbNames, 
								const double volDiff, 
								const double decayDiff, 
								const DoubleArray &jumpDecay,
								const DoubleArray &jumpImpact, 
								const DoubleArray &jumpFreq,
								const DoubleArray &jumpRatios,
								const DoubleArrayArrayArray &jumpMeans,
								const DoubleArrayArrayArray &jumpWeights);
	void CalibrateSingleNames ( DateTime today,
								DateTimeArray &singleNameDates,
								DoubleMatrix &survProbas,  // size (nbNames,singleNameDates.size())
								double idioRatio,
								double cfRatio,
								DoubleArray &jumpRatios);

	bool hasCommonFactor() { return hasCF; };
	bool hasIdioVol() { return hasIdVol; };
	int  getNbMarkets() { return nbMarkets; };
	double getSumFreq() { return poissonDiffuseFast.getSumFreq(); }
	size_t nbJumpFactor() { return 0; }; // TO BE MODIFIED WHEN WE PUT JUMPS

	long getIpathFast() { return IpathFast; };
	void setIpathFast(long _iPath ) { IpathFast = _iPath; };
	long getIpathFull() { return IpathFull; };
	void setIpathFull(long _iPath ) { IpathFull = _iPath; };
	double getRandomNb(bool fast, bool uniform);

	void SimulateCFFast(double futureTime=0.0,   // simulate the common factor CIR and the jump times
				    double lambdaFutureTime=0.0);
	double SimulateJumpsFast (double lastTime,
							bool conditional, 
							double unif = -1.0,
							double futureTime = 0.0);

	void SimulateJumpsFull (double lastTime,
							bool conditional);
	void SimulateSpreadFull(baseWorldSpread &bws);


/* retrieve unconditional survival probabilities **/
	void minusLogSurvProba( int world,  // x axis of world assumed to be the same as cfF and idioF x axis
							int name,
							const DoubleArray &bigT,
							double *minusLogSurvProb,  // (O) must be of size bigT
							bool add = false,
							double futureTime = 0.0,
							double lambdaFutureTimeIdio = 0.0,
							double lambdaFutureTimeCf = 0.0,
							double *lambdaFutureTimeJump = 0);
	


	void minusLogSurvProbaIdio( int world,
								int name,
								const DoubleArray &bigT,
								double *minusLogSurvProb,  // (O) must be of size bigT
								bool add = false,
								double futureTime = 0.0,
								double lambdaFutureTimeIdio = 0.0);

	/* retrieve stochastic survival probability for common factor **/
	void minusLogCondSurvProbaCF(int world,
								int name,
								const DoubleArray &bigT,
								double *minusLogSurvProb,  // (O) must be of size bigT
								bool add = true,
								double futureTime = 0.0);

	void minusLogCondSurvProbaJump( int world,
									int name,
									const DoubleArray &bigT,
									double *minusLogSurvProb,  // (O) must be of size bigT
									bool fastMC,
									bool add = true,
									double futureTime = 0.0,
									double lambdaFutureTimeJumpMarket = 0.0);


	vector<LinearInterpolantSP> f;
	vector<double> fweights;
private:

	UniformRandomSequence rngUniformFast, rngUniformFull;
	GaussianRandomSequence rngNormalFast, rngNormalFull;
	long IpathFast, IpathFull;
	PoissonDiffuseHelp poissonDiffuseFast, poissonDiffuseFull;


	bool hasCF, hasIdVol; 
	CIRAffineSP          cfCIR, cfCIRFull;     
	vector<CIRAffineSP>  idioCIR, idioCIRFull;
	vector< vector< JumpAffineSP > > jumpMarket, jumpMarketFull;
	int nbMarkets;
};

typedef smartPtr<TreeBaseWorld> TreeBaseWorldSP;

DRLIB_END_NAMESPACE

#endif
